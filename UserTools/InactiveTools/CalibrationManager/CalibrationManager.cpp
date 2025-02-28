#include "CalibrationManager.h"

CalibrationManager::CalibrationManager():Tool(){}


bool CalibrationManager::Initialise(std::string configfile, DataModel &data)
{
	if (configfile != "")
	{
		m_configfile = configfile;
		m_variables.Initialise(configfile);
	}
	//m_variables.Print();

	m_data = &data;


	double skip;

	m_variables.Get("verbose", verbose);

	m_variables.Get("calibration",	calibFile);	//file with list of calibration/measurement
	m_variables.Get("output",	m_data->calibrationFile);	//file in which calibration is saved
	m_variables.Get("base_name",	m_data->calibrationBase);	//base_name for calibation
	m_variables.Get("concfunction",	concFuncName);	//name of concentration function
	m_variables.Get("statfunction",	statFuncName);	//name of uncertainity function
	m_variables.Get("systfunction",	systFuncName);	//name of uncertainity function
	m_variables.Get("fit_function",	m_data->fitFuncName);	//name of uncertainity function
	m_variables.Get("skip_calibration",	skip);	//force to use this calibration
	std::cout << " skip :" << skip << "." << std::endl;

	if (skip)
		skipCalib = true;
	else
		skipCalib = false;

	if (skip)
		std::cout << "Calibration will be skipped" << std::endl;
	else
		std::cout << "Calibration will not be skipped" << std::endl;

	m_data->concentrationFunction = 0;
	m_data->concentrationFuncStat = 0;
	m_data->concentrationFuncSyst = 0;
	m_data->calibrationTime = "";

	m_data->gdconc = 0;
	m_data->gd_err = 0;

	m_data->isCalibrationTool = true;

	return true;
}


bool CalibrationManager::Execute()
{
	switch (m_data->mode)
	{
		case state::init:	//about to make measurement, check LED mapping
			Initialise(m_configfile, *m_data);
			calibList = LoadCalibration(allList);

			std::cout << "Number of calibrations out of date or not existing " << calibList.size() << std::endl;

			if (calibList.size()) //calibration out of date!
			{
				m_data->isCalibrated = false;
				m_data->calibrationDone = false;
			}
			else
				m_data->isCalibrated = true;

			ic = calibList.begin();	//global iterator
			break;
		case state::calibration:		//turn on led set
			if (ic != calibList.end())
			{
				Calibrate();
				++ic;
			}

			if (ic == calibList.end())
			{
				m_data->calibrationDone = true;
				std::cout << "Calibration completed\n";
			}
			break;
		case state::calibration_done:
			m_data->isCalibrated = true;
			break;
		case state::finalise:
			Finalise();
			break;
		default:
			break;
	}

	return true;
}

bool CalibrationManager::Finalise()
{
	//clean
	for (ic = allList.begin(); ic != allList.end(); ++ic)
		m_data->DeleteGdTree(m_data->calibrationBase + "_" + *ic);

	calibList.clear();
	allList.clear();

	std::cout << "CM finalise " << m_data->SizeGdTree() << std::endl;

	delete m_data->concentrationFunction;
	delete m_data->concentrationFuncStat;
	delete m_data->concentrationFuncSyst;
	m_data->concentrationFunction = 0;
	m_data->concentrationFuncStat = 0;
	m_data->concentrationFuncSyst = 0;

	timeUpdate.clear();

	return true;
}

//get list of calibration needed
//eg.
//*led0		time1000		#measurement of led0, to be calibrated every 1000 days.
//					#The star indicates this is the concentration measurement
//led0, led1, led2, time54321		#measurement of led0_led1_led2, to be calibrated every 54321 days
//led1 led3, 	time99999		#measurement of led1_led3, to be calibrated every 99999 days
//
std::vector<std::string> CalibrationManager::CalibrationList()
{
	std::vector<std::string> vList;

	std::ifstream cf(calibFile.c_str());
	std::string line;
	while (std::getline(cf, line))
	{
		if (line.find_first_of('#') != std::string::npos)
			line.erase(line.find_first_of('#'));

		if (line.empty())
			continue;

		bool ct;
		size_t ps = line.find_first_of('*');
		if (ps != std::string::npos)
		{
			ct = true;
			line.erase(ps, 1);
		}
		else
			ct = false;
		
		//replace commas into spaces
		std::replace(line.begin(), line.end(), ',', ' ');
		std::replace(line.begin(), line.end(), ';', ' ');

		std::stringstream ssl(line);
		std::string name, word;
		int tu;

		while (ssl)
		{
			if (!std::isalnum(ssl.peek()) && ssl.peek() != '_')
				ssl.ignore();
			else if (ssl >> word)
			{
				if (word.find("time") != std::string::npos)
					tu = std::strtol(word.substr(word.find("time")+4).c_str(), NULL, 10);
				else
					name += word + "_";
			}
		}

		name.erase(name.size()-1);
		vList.push_back(name);

		timeUpdate[name] = tu;

		if (ct)
			concTreeName = name;
	}

	return vList;
}

/* this routine will load all the calibration listed in calibFile
 * and the trees are saved in the data model
 * returns a list of calibrations to do and a list of all calibration, for memory mgmt
 */
std::vector<std::string> CalibrationManager::LoadCalibration(std::vector<std::string> &cList)
{
	//get list of calibrations
	std::vector<std::string> uList;
	cList = CalibrationList();

	std::cout << "loaded " << cList.size() << std::endl;
	//open read only calibration file
	TFile f(m_data->calibrationFile.c_str(), "OPEN");
	if (f.IsZombie())
	{
		std::cout << "No calibration file found at " << m_data->calibrationFile << std::endl;
		if (!skipCalib)
			uList = cList;

		std::cout << cList.size() << "\t" << uList.size() << std::endl;
		for (ic = uList.begin(); ic != uList.end(); )
		{
			std::string cal;
			std::cout << "Do you want to calibrate " << *ic << "? [y/n]" << std::endl;
			std::getline(std::cin, cal);

			while (cal[0] != 'y' && cal[0] != 'n')
			{
				std::cout << "Reply 'yes' or 'no'" << std::endl;
				std::getline(std::cin, cal);
			}

			if (cal[0] == 'y')
			{
				Create(*ic);
				++ic;
			}
			else if (cal[0] == 'n')
			{
				ic = uList.erase(ic);
				std::cout << "Not calibrating" << std::endl;
			}
		}
	}
	else
	{
		std::cout << "Calibration file exists!\n";
		//list of calibration to update/create

		//find latest calibration of each type
		TList *lk = f.GetListOfKeys();
		lk->Sort(kSortDescending);

		//loop on calibration directories
		//they are sorted from last to first
		//directory names are "calibname_timestamp"
		//they contain a GdTree named "calib" and possily functions
		std::string type;
		for (int l = 0; l < lk->GetEntries(); ++l)
		{
			//skip forward because this type is already done
			std::string dir = lk->At(l)->GetName();
			if (!type.empty() && dir.find(type) != std::string::npos)
				continue;

			//loop on list of calibration types, which are "calibname"
			for (ic = cList.begin(); ic != cList.end(); ++ic)
			{
				//name match to find it
				if (dir.find(*ic) != std::string::npos)
				{
					type = *ic;	//store type, needed to fast forward loop


					if (IsUpdate(dir, timeUpdate[type]) || skipCalib)	//calibration is ok
					{
						std::cout << "Calibration of " << type << " up to date!\n";
						Load(f, dir, type);
					}
					else
					{
						std::string cal;
						std::cout << "Calibration of " << type << " out of date!\n";
						std::cout << "Do you want to calibrate it? [y/n]" << std::endl;
						std::getline(std::cin, cal);

						while (cal[0] != 'y' && cal[0] != 'n')
						{
							std::cout << "Reply 'yes' or 'no'" << std::endl;
							std::getline(std::cin, cal);
						}

						if (cal[0] == 'y')
						{
							//Create(type);
							Load(f, dir, type);
							uList.push_back(type);
						}
						else if (cal[0] == 'n' || skipCalib)
						{
							std::cout << "Forcing old calibration" << std::endl;
							Load(f, dir, type);
						}
					}

					//found it, so break
					break;
				}
			}
		}

		f.Close();

		for (ic = cList.begin(); ic != cList.end(); ++ic)
		{
			bool notFound = true;
			std::vector<std::string>::iterator iu;
			for (iu = uList.begin(); iu != uList.end(); ++iu)
				if (*ic == *iu)
				{
					notFound = false;
					break;
				}

			if (notFound && !skipCalib)
			{
				std::string cal;
				std::cout << "Calibration of " << *ic << " not found\n";
				std::cout << "Do you want to calibrate it? (y/n)" << std::endl;
				std::getline(std::cin, cal);

				while (cal[0] != 'y' && cal[0] != 'n')
				{
					std::cout << "Reply 'yes' or 'no'" << std::endl;
					std::getline(std::cin, cal);
				}

				if (cal[0] == 'y')
				{
					Create(type);
					uList.push_back(type);
				}
				else if (cal[0] == 'n')
					std::cout << "not calibrating" << std::endl;
			}
		}
	}

	return uList;
}

void CalibrationManager::Load(TFile &f, std::string dir, std::string type)
{
	std::string name = m_data->calibrationBase + "_" + type;
	std::string full = dir + "/" + m_data->calibrationBase;
	//std::string full = dir + "/" + name;

	TTree *gdt = static_cast<TTree*>(f.Get(full.c_str()));
	GdTree *calib = new GdTree(gdt);

	m_data->AddGdTree(name, calib);

	if (type == concTreeName)	//happens once
	{
		m_data->concentrationName = name;

		std::string f1name = dir + "/" + concFuncName;
		std::string f2name = dir + "/" + statFuncName;
		std::string f3name = dir + "/" + systFuncName;

		//clone functions (no SetDirectory)
		m_data->concentrationFunction = static_cast<TF1*>(f.Get(f1name.c_str())->Clone());
		m_data->concentrationFuncStat = static_cast<TF2*>(f.Get(f2name.c_str())->Clone());
		m_data->concentrationFuncSyst = static_cast<TF2*>(f.Get(f3name.c_str())->Clone());
	}
}

void CalibrationManager::Create(std::string type)
{
	std::string name = m_data->calibrationBase + "_" + type;
	GdTree *calib = new GdTree(name.c_str());
	m_data->AddGdTree(name.c_str(), calib);

	if (type == concTreeName)	//happens once
	{
		m_data->concentrationName = name;

		//create functions for absorbance, stat and syst error
		m_data->concentrationFunction = new TF1(concFuncName.c_str(),
				"(x - [0])/[1]", 0, 1);
		m_data->concentrationFuncStat = new TF2(statFuncName.c_str(),
				"((x - [0])/[1])^2 * ( (y/(x - [0]))^2 )", 0, 1);
		m_data->concentrationFuncSyst = new TF1(systFuncName.c_str(),
				"((x - [0])/[1])^2 * ( ([2]/(x - [0]))^2 + ([3]/[1])^2 + 2*[4] / (x - [0]) / [1] )", 0, 1);
	}

	std::cout << "time " << m_data->calibrationTime << std::endl;
	if (m_data->calibrationTime.empty())
	{
		//	    012345678901234
		//format is YYYYMMDDTHHMMSS (there shouldn't be fractional seconds
		boost::posix_time::ptime current(boost::posix_time::second_clock::local_time());
		std::string tc = boost::posix_time::to_iso_string(current);
		m_data->calibrationTime = tc.substr(0, tc.find_first_of('T'));
		std::cout << "time " << m_data->calibrationTime << std::endl;
	}
}

//directory names are "calibname_timestamp"
//with timestamp an iso string YYYYMMDDThhmmss
bool CalibrationManager::IsUpdate(std::string name, int timeUpdate)
{
	std::string ts = name.substr(name.find_last_of('_') + 1);
	ts += "T000000";
	std::cout << "calibtraion found " << ts << std::endl;

	//get current time
	boost::posix_time::ptime last(boost::posix_time::from_iso_string(ts));
	boost::posix_time::ptime current(boost::posix_time::second_clock::local_time());
	boost::posix_time::time_duration lapse(current - last);

	std::cout << "lapse " << lapse.hours() << "\t" << timeUpdate << std::endl;
	if (lapse.hours()/24.0 > timeUpdate)
		return false;
	else
		return true;
}

//calibrate is like a measurement
bool CalibrationManager::Calibrate()
{
	m_data->depleteWater = false;
	m_data->circulateWater = false;
	m_data->calibrationName = m_data->calibrationBase + "_" + *ic;
	m_data->turnOnLED = *ic;

	std::string input;
	double gdc, gde;
	std::cout << "Calibrating " << *ic << std::endl;

	std::string cal;
	std::cout << "\nDo you need to change water? [y/n]" << std::endl;
	std::getline(std::cin, cal);

	while (cal[0] != 'y' && cal[0] != 'n')
	{
		std::cout << "Reply 'yes' or 'no'" << std::endl;
		std::getline(std::cin, cal);
	}

	if (cal[0] == 'y')
	{
		m_data->depleteWater = true;
		m_data->circulateWater = false;
		--ic;

		return true;
	}
	else if (cal[0] == 'n')
	{
		if (*ic != concTreeName)
		{
			m_data->depleteWater = false;
			m_data->circulateWater = true;

			std::cout << "Loading water in system." << std::endl;

			m_data->gdconc = 0.0;
			m_data->gd_err = 0.0;
			//m_data->calibrationComplete = true;
		}
		else
		{
			std::cout << "Enter next concentration (-1 to quit, 0.0 for pure water).";
			std::cout << " Last concentration: " << m_data->gdconc << std::endl;
			std::getline(std::cin, input);
			gdc = std::strtod(input.c_str(), NULL);
			m_data->gdconc = gdc;
			if (!(gdc < 0))	//not finished, repeat this step
			{
				std::cout << "Enter concentration error (0.0 is a fine value)\n";
				std::getline(std::cin, input);
				gde = std::strtod(input.c_str(), NULL);

				m_data->gd_err = gde*gde;
				//m_data->calibrationComplete = false;

				--ic;

				m_data->depleteWater = false;
				m_data->circulateWater = true;

				std::cout << "Loading water in system." << std::endl;
			}
			else
			{
				//m_data->calibrationComplete = true;
				m_data->gd_err = 0.0;;

				m_data->depleteWater = false;
				m_data->circulateWater = false;
			}
		}
	}

	return true;
}
