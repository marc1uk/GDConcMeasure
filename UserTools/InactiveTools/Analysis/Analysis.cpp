#include "Analysis.h"

Analysis::Analysis():Tool(){}


bool Analysis::Initialise(std::string configfile, DataModel &data)
{

	if(configfile!="")
	{
		m_configfile = configfile;
		m_variables.Initialise(configfile);
	}
	//m_variables.Print();

	m_data = &data;
	analysis = Type::undefined;

	//get the order of the fitting polynomial for baseline removal
	m_variables.Get("order", polyOrder);

	return true;
}


bool Analysis::Execute()
{	
	//analyse means data was loaded in datamodel
	//and trace can be averaged
	switch (m_data->mode)
	{
		case state::init:
			iA = iB = -1;	//reset region definition
			pureTrace.clear();
			//darkTrace = AverageTrace(0);
			break;
		case state::calibration:	//taking dark current
			analysis = Type::calibrate;
			darkTrace = AverageTrace(0);
			break;
		case state::measurement:	//taking dark current
			analysis = Type::measure;
			darkTrace = AverageTrace(0);
			break;
		case state::take_spectrum:	//average data
			avgTrace = AverageTrace(1);
			break;
		case state::analyse:	//led off again, we can analyse
			if (analysis == Type::calibrate)
			{
				if (!m_data->calibrationComplete || m_data->calibrationName != m_data->concentrationName)
					CalibrationTrace(m_data->gdconc, m_data->gd_err);
				else
					CalibrationComplete();	//complete missing holes in calibration-concentration
			}
			else if (analysis == Type::measure)
			{
				pureTrace = PureTrace();
				MeasurementTrace();
			}
			else
				std::cout << "Skip" << std::endl;

			DefineRegion(pureTrace);
			break;
		case state::calibration_done:
			LinearFit();		//evaluate linear fit for concentration calibration
			break;
		case state::measurement_done:
			break;
		case state::finalise:
			Finalise();
			break;
		default:
			break;
	}

	return true;
}

bool Analysis::Finalise()
{
	std::cout << "AA finalise " << m_data->SizeGdTree() << std::endl;

	darkTrace.clear();
	pureTrace.clear();
	return true;
}

void Analysis::DefineRegion(const std::vector<double> &trace)
{
	//if there is no pure trace or if region already defined, skip
	if (!trace.size() || reg_posx < 0)
		return;

	int size = trace.size() / 2;
	double posx = 0, perr = 0, weig = 0;
	for (int i = 1; i < size; ++i)
	{
		weig += pow(trace[i], 2);
		posx += pow(trace[i], 2) * m_data->wavelength[i];
		perr += pow(trace[i] * m_data->wavelength[i], 2);
	}

	if (weig > 0)
	{
		posx /= weig;
		perr /= weig;
	}

	perr -= posx*posx;
	std::cout << "got " << perr << " +/- " << sqrt(perr) << std::endl;

	reg_posx = posx;
	reg_perr = std::sqrt(perr);

	/*
	iA = static_cast<int>(floor(posx - 3*sqrt(perr)));
	iB = static_cast<int>(ceil (posx + 3*sqrt(perr)));

	if (iA < 0)
		iA = 0;
	if (iB > size)
		iB = size;
	*/

	std::cout << "Region is 3-sigma defined in: "
		  << reg_posx - 3 * reg_perr << " -> " << reg_posx
		  << " <- " << reg_pos + 3 * reg_perr << std::endl;
}

//obtain puretrace from existing calibration
std::vector<double> Analysis::PureTrace()
{
	GdTree *calib = m_data->GetGdTree(m_data->calibrationName);
	if (!calib)
	{
		std::cout << "Calibration not found/set" << std::endl;
		return std::vector<double>();
	}

	//std::cout << "ENTRIES " << calib->GetEntries() << std::endl;
	std::vector<double> trace;
	for (int i = 0; i < calib->GetEntries(); ++i)
	{
		calib->GetEntry(i);
		if (calib->GdConc == 0)
		{
			trace = calib->Trace;
			trace.insert(trace.end(), calib->T_Err.begin(), calib->T_Err.end());
			break;
		}
	}

	return trace;
}

std::vector<double> Analysis::AverageTrace(bool darkRemove)
{
	int size = m_data->traceCollect.front().size();
	int n = m_data->traceCollect.size();
	std::cout << "averagin " << n << " traces of size " << size << std::endl;
	std::vector<double> avgTrace(2 * size);
	if (size)
	{
		std::vector<std::vector<double> >::iterator it;
		for (it = m_data->traceCollect.begin(); it != m_data->traceCollect.end(); ++it)
		{
			for (int i = 0, e = size; i < size; ++i, ++e)	//loop "size" entries
			{
				avgTrace[i] += (*it)[i] / n;
				avgTrace[e] += (*it)[i] * (*it)[i] / n / (n - 1);
			}
		}

		darkRemove = (darkRemove && darkTrace.size());

		for (int i = 0, e = size; i < size; ++i, ++e)
		{
			avgTrace[e] -= pow(avgTrace[i], 2) / (n - 1);

			if (darkRemove)	//remove background, aka dark current
			{
				avgTrace[i] -= darkTrace[i];
				avgTrace[e] += darkTrace[e];
			}

			if (avgTrace[i] < 0)
				avgTrace[i] = 0.0;
		}
	}
	else
		std::cout << "Houston, we have a problem\n";	//-> no data from spectrometer

	return avgTrace;
}

std::vector<double> Analysis::AbsorbTrace(const std::vector<double> &avgTrace)
{
	if (!pureTrace.size())				//return empty vector;
		return std::vector<double>();

	int size = avgTrace.size() / 2;
	std::vector<double> absTrace(2 * size);
	std::vector<double> xAxis(size), toFit(size);

	//intervals for background fitting
	//adding this sigma values to the map, the value is updated by BinarySearch
	std::map<int, int> iX;
	iX[-4];
	iX[-2];
	iX[+2];
	iX[+4];

	iX = BinarySearch(iX, reg_posx, reg_perr);

	if (pureTrace.size())
	{
		for (int i = 0, e = size; i < size; ++i, ++e)
		{
			//it means it is zero
			//if (std::abs(pureTrace[i]) < std::sqrt(pureTrace[e]))
			//{
			//	absTrace[e] = 0.0;
			//	absTrace[i] = 0.0;
			//}
			//else

			absTrace[e] = pureTrace[e] / pow(pureTrace[i]  * log(10), 2) +
				      avgTrace[e]  / pow(avgTrace[i]   * log(10), 2);

			absTrace[i] = log10(pureTrace[i] / avgTrace[i]); 

			if (!std::isfinite(absTrace[i]))
				absTrace[i] = 0.0;

			if ((i >= iX[-4] && i <= iX[-2]) || (i >= iX[2] && i <= iX[4]))
			{
				xAxis.push_back(m_data->wavelength[i] - reg_posx);
				toFit.push_back(absTrace[i]);
			}
		}

		PolyFit pf(xAxis, toFit, polyOrder);
		std::vector<double> beta = pf.Solve();

		//correct absorbance trace removing baseline
		for (int i = 0; i < size; ++i)
			absTrace[i] -= pf.ff(m_data->wavelength[i] - reg_posx);
	}

	return absTrace;
}

//compute absorbance and find the two gd peaks
std::vector<double> Analysis::Absorbance(const std::vector<double> &avgTrace, int &i1, int &i2)
{
	i1 = -1;
	i2 = -1;

	if (!pureTrace.size())				//return empty vector;
		return std::vector<double>();

	std::vector<int> iPeak, iDeep;
	std::vector<double> absTrace = AbsorbTrace(avgTrace);

	FindPeakDeep(absTrace, iPeak, iDeep);

	if (iPeak.size() > 0)
		i1 = iPeak.at(0);
	if (iPeak.size() > 1)
		i2 = iPeak.at(1);

	//std::cout << "i1 " << i1 << "\t i2 " << i2 << std::endl;

	return absTrace;
}

//////////////////////////////
/* fill calibration tree, with given gdconc e gd_err
 * if gdconc == 0, the pure trace is stored 
 */
void Analysis::CalibrationTrace(double gdconc, double gd_err)
{
	//std::vector<double> avgTrace = AverageTrace(1);
	if (!avgTrace.size())	//empty? skip, it is an error
	{
		std::cout << "No data from spectrometer!!!! DON'T PANIC!!!!" << std::endl;
		return;
	}

	int size = avgTrace.size() / 2;
	GdTree *tree = m_data->GetGdTree(m_data->calibrationName);

	tree->GdConc = gdconc;
	tree->Gd_Err = gd_err;
	tree->Gd_ErrStat = 0;
	tree->Gd_ErrSyst = 0;

	tree->Wavelength = m_data->wavelength;
	tree->Trace = std::vector<double>(avgTrace.begin(), avgTrace.begin() + size);
	tree->T_Err = std::vector<double>(avgTrace.begin() + size, avgTrace.end());

	std::vector<double> absTrace;
	if (gdconc == 0)		//if not empty and gdconc = 0, we are measuring pureTrace!
	{
		pureTrace = avgTrace;	//store this trace
		FillAbsorbance(tree, absTrace, -1, -1, 1);
	}
	else
	{
		int i1, i2;
		std::vector<double> absTrace = Absorbance(avgTrace, i1, i2);
		FillAbsorbance(tree, absTrace, i1, i2, 1);
	}

	tree->Fill();

	/* if there is a puretrace then continue doing calibration and filling the tree
	 * if not, fill the tree partially and complete it later, once the pure water trace
	 * has been taken (hopefully before the end of the calibration)
	 * a few entries are filled with empty vectors or negative numbers which
	 * denotes a incomplete calibration
	 */
	//if (pureTrace.size())
	//{
	//	std::vector<double> absTrace = Absorbance(avgTrace, i1, i2);
	//	tree->Absorbance = std::vector<double>(absTrace.begin(), absTrace.begin() + size);
	//	tree->Absorb_Err = std::vector<double>(absTrace.begin() + size, absTrace.end());
	//}
	//else	//pure trace not taken yet
	//{
	//	//mustComplete = true;			//at the end of calibration, complete it
	//	tree->Absorbance = std::vector<double>(size);
	//	tree->Absorb_Err = std::vector<double>(size);
	//}
}

/////////////////
/* complete calibration
 * run over the tree of calibration measurement
 * and extract abosrbance of entries without one.
 * this can happen when the pure water sample is not measured as first one.
 * Since trees are write-once-only object, a new tree must be recreated
 *
 */
void Analysis::CalibrationComplete()
{
	if (!pureTrace.size())	//this is bad
	{
		std::cout << "Calibration incomplete, due to missing pure water trace" << std::endl;
		return;
	}

	std::cout << "CC completion" << std::endl;

	GdTree *oldTree = m_data->GetGdTree(m_data->calibrationName);	//should be current calibration tree
	GdTree *newTree = new GdTree(m_data->calibrationName);		//new tree for copy/edit

	for (int i = 0; i < oldTree->GetEntries(); ++i)
	{
		oldTree->GetEntry(i);

		newTree->GdConc	     = oldTree->GdConc;
		newTree->Gd_Err	     = oldTree->Gd_Err;
		newTree->Gd_ErrStat = oldTree->Gd_ErrStat;
		newTree->Gd_ErrSyst = oldTree->Gd_ErrSyst;

		newTree->Wavelength = oldTree->Wavelength;
		newTree->Trace	= oldTree->Trace;
		newTree->T_Err	= oldTree->T_Err;

		newTree->Epoch	= oldTree->Epoch;
		newTree->Year	= oldTree->Year;
		newTree->Month	= oldTree->Month;
		newTree->Day	= oldTree->Day;
		newTree->Hour	= oldTree->Hour;
		newTree->Minute	= oldTree->Minute;
		newTree->Second	= oldTree->Second;

		//need both peak to complete calibration
		if (oldTree->Wavelength_1 < 0 || oldTree->Wavelength_2 < 0)
		{
			std::cout << "CC redoing absorbance" << std::endl;
			int size = oldTree->Trace.size();
			std::vector<double> avgTrace = oldTree->Trace;
			avgTrace.insert(avgTrace.end(), oldTree->T_Err.begin(), oldTree->T_Err.end());

			int i1 = -1.0, i2 = -1.0;
			std::vector<double> absTrace = Absorbance(avgTrace, i1, i2);

			FillAbsorbance(newTree, absTrace, i1, i2, 0);	//no timestamp
		}
		else
		{
			std::cout << "CC just copying" << std::endl;
			newTree->Absorbance = oldTree->Absorbance;
			newTree->Absorb_Err = oldTree->Absorb_Err;

			newTree->Wavelength_1 = oldTree->Wavelength_1;
			newTree->Absorbance_1 = oldTree->Absorbance_1;
			newTree->Absorb_Err_1 = oldTree->Absorb_Err_1;

			newTree->Wavelength_2 = oldTree->Wavelength_2;
			newTree->Absorbance_2 = oldTree->Absorbance_2;
			newTree->Absorb_Err_2 = oldTree->Absorb_Err_2;

			newTree->Absorb_Diff = oldTree->Absorb_Diff;
			newTree->AbsDiff_Err = oldTree->AbsDiff_Err;
		}

		newTree->Fill();
	}

	m_data->DeleteGdTree(m_data->calibrationName);	//oldtree
	m_data->AddGdTree(m_data->calibrationName, newTree);	//replace with new tree
}


void Analysis::MeasurementTrace()
{
	//std::vector<double> avgTrace = AverageTrace(1);
	if (!avgTrace.size())	//empty? skip, it is an error
	{
		std::cout << "No data from spectrometer!!!!" << std::endl;
		return;
	}

	int size = avgTrace.size() / 2;
	GdTree *tree = m_data->GetGdTree(m_data->measurementName);

	int i1 = -1, i2 = -1;
	std::vector<double> absTrace = Absorbance(avgTrace, i1, i2);

	//store traces
	tree->Wavelength = m_data->wavelength;
	tree->Trace = std::vector<double>(avgTrace.begin(), avgTrace.begin() + size);
	tree->T_Err = std::vector<double>(avgTrace.begin() + size, avgTrace.end());

	//store peaks
	FillAbsorbance(tree, absTrace, i1, i2, 1);

	//evaluate concetration + errors
	tree->GdConc = m_data->concentrationFunction->Eval(tree->Absorb_Diff);
	tree->Gd_ErrStat = m_data->concentrationFuncStat->Eval(tree->Absorb_Diff, sqrt(tree->AbsDiff_Err));
	tree->Gd_ErrSyst = m_data->concentrationFuncSyst->Eval(tree->Absorb_Diff);	//doesn't depend on abs err
	tree->Gd_Err = tree->Gd_ErrStat + tree->Gd_ErrSyst;
	
	tree->Fill();
}

void Analysis::FillAbsorbance(GdTree *tree, std::vector<double> &abst,
			      int i1, int i2, bool stamp)
{
	int size = abst.size() / 2;
	if (size)
	{
		tree->Absorbance = std::vector<double>(abst.begin(), abst.begin() + size);
		tree->Absorb_Err = std::vector<double>(abst.begin() + size, abst.end());
	}
	else
	{
		tree->Absorbance = std::vector<double>();
		tree->Absorb_Err = std::vector<double>();
	}

	if (i1 >= 0 && i1 < size)	//found one peak
	{
		tree->Wavelength_1 = m_data->wavelength.at(i1);
		tree->Absorbance_1 = abst.at(i1);
		tree->Absorb_Err_1 = abst.at(i1 + size);

		tree->Absorb_Diff = tree->Absorbance_1;
		tree->AbsDiff_Err = tree->Absorb_Err_1;
	}
	else
	{
		std::cout << "Missing main peak" << std::endl;
		tree->Wavelength_1 = -1.0;
		tree->Absorbance_1 = -1.0;
		tree->Absorb_Err_1 = -1.0;

		tree->Absorb_Diff = 0.0;
		tree->AbsDiff_Err = 0.0;
	}

	if (i2 >= 0 && i2 < size)	//found second peak
	{
		tree->Wavelength_2 = m_data->wavelength.at(i2);
		tree->Absorbance_2 = abst.at(i2);
		tree->Absorb_Err_2 = abst.at(i2 + size);

		tree->Absorb_Diff -= tree->Absorbance_2;
		tree->AbsDiff_Err += tree->Absorb_Err_2;
	}
	else
	{
		std::cout << "Missing secondary peak" << std::endl;
		tree->Wavelength_2 = -1.0;
		tree->Absorbance_2 = -1.0;
		tree->Absorb_Err_2 = -1.0;

		tree->Absorb_Diff = 0.0;
		tree->AbsDiff_Err = 0.0;
	}

	if (stamp)
	{
		int year, month, day, hour, minute, second;
		tree->Epoch = TimeStamp(year, month, day, hour, minute, second);

		tree->Year	= year;
		tree->Month	= month;
		tree->Day	= day;
		tree->Hour	= hour;
		tree->Minute	= minute;
		tree->Second	= second;
	}
}

void Analysis::LinearFit()
{
	/* calibration tree will contain n entries,
	 * where n-1 is the number of concentration probed
	 * and 1 is the pure water reference
	 */
	GdTree *tree = m_data->GetGdTree(m_data->concentrationName);
	if (!tree)
		return;

	int n = tree->GetEntries();

	/* creating a graph with calibration points
	 * to be fitted with a+b*x
	 */
	TGraphErrors *gd = new TGraphErrors(n-1);
	int i = 0, j = 0;
	double absMax = -1e6, absMin = 1e6;
	while (i < n)
	{
		++i;
		tree->GetEntry(i);

		if (tree->Absorb_Diff > absMax)
			absMax = tree->Absorb_Diff;

		if (tree->Absorb_Diff < absMin)
			absMin = tree->Absorb_Diff;

		if (tree->GdConc > 0 && tree->Wavelength_1 > 0 && tree->Wavelength_2 > 0)
		{
			gd->SetPoint(j, tree->GdConc, tree->Absorb_Diff);
			gd->SetPointError(j, sqrt(tree->Gd_Err), sqrt(tree->AbsDiff_Err));
			++j;
		}
		else
			std::cout << "Not fitting point at concentration " << tree->GdConc << std::endl;
	}

	//extending range
	--absMin;
	++absMax;
	std::cout << "Fitting " << j << " points in the range [ " << absMin << " : " << absMax << " ]\n";

	//keeps point sorted in X diretion
	gd->Sort();

	/* creat function a + b*x for fit the graph
	 */
	TF1 *line = new TF1("dep", "1 ++ x", absMin, absMax);
	m_data->concentrationFit = gd->Fit(line, "SRNQ");
	int stat = m_data->concentrationFit;

	if (stat >= 0)	//successful fit
	{
		//print fit info
		m_data->concentrationFit->Print("V");

		double a_val = m_data->concentrationFit->Parameter(0);	//[0]
		double b_val = m_data->concentrationFit->Parameter(1);	//[1]
		double a_err = m_data->concentrationFit->ParError(0);	//[2]
		double b_err = m_data->concentrationFit->ParError(1);	//[3]

		double abcov = m_data->concentrationFit->GetCovarianceMatrix()(0,1);

		std::cout << "Covariance a, b " << abcov << std::endl;

		/* invert function line to have
		 * concentration as a function of absorbance x = (y-a)/b
		 */

		if (m_data->concentrationFunction && m_data->concentrationFuncStat && m_data->concentrationFuncSyst)
		{
			m_data->concentrationFunction->SetRange(absMin, absMax);
			m_data->concentrationFunction->SetParameter(0, a_val);
			m_data->concentrationFunction->SetParameter(1, b_val);

			m_data->concentrationFuncStat->SetRange(absMin, absMax);
			m_data->concentrationFuncStat->SetParameter(0, a_val);
			m_data->concentrationFuncStat->SetParameter(1, b_val);

			m_data->concentrationFuncSyst->SetRange(absMin, absMax);
			m_data->concentrationFuncSyst->SetParameter(0, a_val);
			m_data->concentrationFuncSyst->SetParameter(1, b_val);
			m_data->concentrationFuncSyst->SetParameter(2, a_err);
			m_data->concentrationFuncSyst->SetParameter(3, b_err);
			m_data->concentrationFuncSyst->SetParameter(4, abcov);
		}
		else
			std::cerr << "Caution, function errors not created\n";
	}
	else
		std::cout << "Fit not successful\n";

	delete gd;
	delete line;
}

/* this routine can be improved
 * it is not robust and stable
 * but it works quite ok for now
 */
void Analysis::FindPeakDeep(const std::vector<double> &trace, std::vector<int> &iPeak, std::vector<int> &iDeep)
{
	/* Find maximum (and min) value and position
	 * this is needed to find the peaks of a trace
	 */

	//std::vector<double> vx = m_data->wavelength;

	/*
	int size = pureTrace.size() / 2;
	double posx = 0, perr = 0, weig = 0;
	for (int i = 1; i < size; ++i)
	{
		weig += pow(pureTrace[i], 2);
		posx += pow(pureTrace[i], 2) * i;
		perr += pow(pureTrace[i] * i, 2);
	}

	if (weig > 0)
	{
		posx /= weig;
		perr /= weig;
	}

	perr -= posx*posx;
	std::cout << "got " << perr << " +/- " << sqrt(perr) << std::endl;

	iA = static_cast<int>(floor(posx - 3*sqrt(perr)));
	iB = static_cast<int>(ceil (posx + 3*sqrt(perr)));

	std::cout << "region is defined in " << iA << " -> " << posx << " <- " << iB << std::endl;


	if (iA < 0 || iA > size)
		iA = 0;
	if (iB < 0 || iB > size)
		iB = size-1;
	*/

	if (reg_psx < 0)
		return;

	int size = trace.size() / 2;
	std::map<int, int> iX;
	iX[-3];
	iX[+3];

	iX = BinarySearch(iX, reg_posx, reg_perr);

	int iA = iX[0];
	int iB = iX[1];

	//searching max and min between iA and iB
	double fMax = trace[iA], fMin = trace[iA];
	int iMax = iA, iMin = iA;

	for (int i = iA+1; i < iB+1; ++i)
	{
		if (fMax < 0 || fMax < trace[i])
		{
			fMax = trace[i];
			iMax = i;
		}

		if (fMin < 0 || fMin > trace[i])
		{
			fMin = trace[i];
			iMin = i;
		}
	}

	/* a more robust defintion of threshold must be found
	 */
	double thr = 0.10 * (fMax - fMin);

	iPeak.clear();
	iDeep.clear();
	iPeak.push_back(iMax);

	//this depends on the LED, be careful
	//int iA = 0, iB = size;

	/* start peak searching from highest peak
	 * going left first, and then right
	 */ 
	for (int Dir = -1; Dir < 2; Dir += 2)
	{
		double PoV = false;		//F looking for Deep, T looking for P
		int iD = iMax, iS = iMax;
		double fS = fMax;
		/* if the last poinst found was a peak,
		 * the the algorithm searches for a deep
		 * and viceversa
		 */
		while (iD >= iA && iD < iB)
		{
			int Sign = 2*PoV - 1;	//-1 looking for Deep, +1 looking for Peak
			if ( Sign*(trace[iD] - fS) > thr)
			{
				iS = iD;
				fS = trace[iD];
			}

			if (-Sign*(trace[iD] - fS) > thr)
			{
				double fZ = -Sign*fMax;
				int    iZ = -1;
				for (int j = iD; Dir*(iS-j) < 1 && j > -1 && j < size; j -= Dir)
				{
					if (Sign*(trace[j] - fZ) > 0)
					{
						fZ = trace[j];
						iZ = j;
					}
				}
				if (iZ > -1)
				{
					if (PoV)
						iPeak.push_back(iZ);
					else
						iDeep.push_back(iZ);

					fS = fZ;
					iD = iS = iZ;
					PoV = !PoV;
				}
			}
			iD += Dir;
		}
	}

	std::sort(iPeak.begin(), iPeak.end(), Sort(trace,  1));
	std::sort(iDeep.begin(), iDeep.end(), Sort(trace, -1));
}

ULong64_t Analysis::TimeStamp(int &Y, int &M, int &D, int &h, int &m, int &s)
{
	boost::posix_time::ptime current(boost::posix_time::second_clock::local_time());

	//format is YYYYMMDDTHHMMSS (there shouldn't be fractional seconds
	//	    012345678901234
	std::string tc = boost::posix_time::to_iso_string(current);
	//std::string date = tc.substr(0, tc.find_first_of('T'));
	//std::string time = tc.substr(tc.find_first_of('T')+1);

	Y = std::strtol(tc.substr(0,  4).c_str(), NULL, 10);
	M = std::strtol(tc.substr(4,  2).c_str(), NULL, 10);
	D = std::strtol(tc.substr(6,  2).c_str(), NULL, 10);
	h = std::strtol(tc.substr(9,  2).c_str(), NULL, 10);
	m = std::strtol(tc.substr(11, 2).c_str(), NULL, 10);
	s = std::strtol(tc.substr(13, 2).c_str(), NULL, 10);

	boost::posix_time::ptime epoch(boost::gregorian::date(1970,1,1));
	return (current - epoch).total_seconds();
}

void Analysis::BinarySearch(std::vector<int, int> &nSigma, double x, double e)
{
	std::map<int, int>::iterator is = nSigma.begin();
	for ( ; is != nSigma.end(); ++is)
	{
		int jL = 0, jH = m_data->wavelength.size() - 1, iN = 0;

		while (std::abs(jL - jH) <= 1)
		{
			if (is->first <= 0)
				iN = int((jL + jH) / 2);
			else
				iN = int((jL + jH + 1) / 2);

			if (m_data->wavelength[jL] < reg_posx + is->first * reg_perr &&
			    reg_posx + is->first * reg_perr < m_data->wavelength[iN])
				jH = iN;
			else if (m_data->wavelength[iN] < reg_posx + is->first * reg_perr &&
				 reg_posx + is->first * reg_perr < m_data->wavelength[jH])
				jL = iN;
		}

		is->second = iN;
	}
}

