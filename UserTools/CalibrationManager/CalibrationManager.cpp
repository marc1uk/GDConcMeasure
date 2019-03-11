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

	m_variables.Get("verbose", verbose);

	m_variables.Get("calibration",	calibFile);	//file with list of calibration/measurement
	m_variables.Get("output",	m_data->calibrationFile);	//file in which calibration is saved
	m_variables.Get("base_name",	m_data->calibrationBase);	//base_name for calibation
	m_variables.Get("concfunction",	concFuncName);	//name of concentration function
	m_variables.Get("err_function",	err_FuncName);	//name of uncertainity function


	m_data->concentrationFunction = 0;
	m_data->concentrationFunc_Err = 0;

	return true;
}


bool CalibrationManager::Execute()
{
	switch (m_data->mode)
	{
		case state::init:	//about to make measurement, check LED mapping
			Initialise(m_configfile, *m_data);
			calibList = LoadCalibration();

			std::cout << "TODO number of calibrations " << calibList.size() << std::endl;

			if (calibList.size()) //calibration out of date!
			{
				m_data->isCalibrated = false;
				m_data->calibrationDone = false;
			}
			else
				m_data->isCalibrated = true;

			ic = calibList.begin();	//global iterator
			std::cout << "I3\n";
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
	for (ic = calibList.begin(); ic != calibList.end(); ++ic)
		m_data->DeleteGdTree(*ic);

	calibList.clear();

	delete m_data->concentrationFunction;
	delete m_data->concentrationFunc_Err;
	m_data->concentrationFunction = 0;
	m_data->concentrationFunc_Err = 0;

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
 * returns a list of calibrations to do
 */
std::vector<std::string> CalibrationManager::LoadCalibration()
{
	//get list of calibrations
	std::vector<std::string> cList = CalibrationList();
	std::vector<std::string> uList;

	//open read only calibration file
	TFile f(m_data->calibrationFile.c_str(), "OPEN");
	if (f.IsZombie())
	{
		std::cerr << "Calibration file does not exist! A new one will be created\n";
		uList = cList;
		for (ic = uList.begin(); ic != uList.end(); ++ic)
			Create(*ic);
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
			std::cout << "dir found " << dir << std::endl;
			if (!type.empty() && dir.find(type) != std::string::npos)
				continue;

			//loop on list of calibration types, which are "calibname"
			for (ic = cList.begin(); ic != cList.end(); ++ic)
			{
				std::cout << "match " << *ic << std::endl;
				//name match to find it
				if (dir.find(*ic) != std::string::npos)
				{
					type = *ic;	//store type, needed to fast forward loop
					std::cout << "latest " << type << std::endl;

					if (IsUpdate(dir, timeUpdate[type]))	//calibration is ok
					{
						std::cout << "is update" << std::endl;
						Load(f, dir, type);
					}
					else
					{
						std::cout << "not update" << std::endl;
						Create(type);
						uList.push_back(type);
					}

					//found it, so break
					break;
				}
			}
		}

		f.Close();
	}

	return uList;
}

void CalibrationManager::Load(TFile &f, std::string dir, std::string type)
{
	std::string name = m_data->calibrationBase + "_" + type;
	std::string full = dir + "/" + name;

	std::cout << "LOOKING " << full << std::endl;
	TTree *gdt = static_cast<TTree*>(f.Get(full.c_str()));
	std::cout << "LOADED GDT " << gdt->GetEntries() << std::endl;
	GdTree *calib = new GdTree(gdt);

	std::cout << "LOADED calibraiton " << calib->GetTree()->GetEntries() << " w/ " << name << std::endl;
	std::cout << "LOADED entries " << calib->GetEntries() << std::endl;

	m_data->AddGdTree(name, calib);

	if (type == concTreeName)	//happens once
	{
		m_data->concentrationName = name;

		std::string f1name = name + "/" + concFuncName;
		std::string f2name = name + "/" + err_FuncName;

		//clone tree and functions
		m_data->concentrationFunction = static_cast<TF1*>(f.Get(f1name.c_str())->Clone());
		m_data->concentrationFunc_Err = static_cast<TF2*>(f.Get(f2name.c_str())->Clone());
	}
}

void CalibrationManager::Create(std::string type)
{
	std::string name = m_data->calibrationBase + "_" + type;
	GdTree *calib = new GdTree(name.c_str());
	m_data->AddGdTree(name.c_str(), calib);
	std::cout << "creating " << name << std::endl;

	if (type == concTreeName)	//happens once
	{
		m_data->concentrationName = name;

		//create functions for absorbance, value and error
		m_data->concentrationFunction = new TF1(concFuncName.c_str(),
				"(x - [0])/[1]", 0, 1);
		m_data->concentrationFunc_Err = new TF2(err_FuncName.c_str(),
				"((x - [0])/[1])^2 * ( (y^2 + [2]^2)/(x - [0])^2 + ([3]/[1])^2 )", 0, 1);
	}
}

//directory names are "calibname_timestamp"
//with timestamp an iso string YYYYMMDDThhmmss
bool CalibrationManager::IsUpdate(std::string name, int timeUpdate)
{
	std::cout << "update?" << std::endl;
	std::string ts = name.substr(name.find_last_of('_') + 1);
	std::cout << "U0 " << ts << std::endl;

	//get current time
	boost::posix_time::ptime last(boost::posix_time::from_iso_string(ts));
	boost::posix_time::ptime current(boost::posix_time::second_clock::local_time());
	boost::posix_time::time_duration lapse(current - last);

	std::cout << "lapse " << lapse.hours() << "\t" << timeUpdate << std::endl;
	if (lapse.hours()/24.0 > timeUpdate)
	{
		std::cout << "Calibration out of date!\n";
		return false;
	}
	else
		return true;
}

//calibrate is like a measurement
bool CalibrationManager::Calibrate()
{
	std::string acknowledge;	// <--- this is some form of aknoledgment
	double gdc, gde;
	if (*ic != concTreeName)
	{
		std::cout << "Calibrating " << *ic << std::endl;
		std::cout << "Have you put pure water? Enter to continue...\n";
		std::cin.get();
		m_data->gdconc = 0.0;
		m_data->gd_err = 0.0;
		m_data->calibrationComplete = true;
	}
	else
	{
		std::cout << "Calibrating " << *ic << std::endl;
		std::cout << "Enter concentration loaded in cell now (-1 to quit, 0 for pure water) \n";
		std::cin >> gdc;
		if (!(gdc < 0))	//not finished, repeat this step
		{
			std::cout << "Enter concentration error (0 is a fine value)\n";
			std::cin >> gde;
			m_data->gd_err = gde;
			m_data->calibrationComplete = false;

			--ic;
		}
		else
		{
			m_data->calibrationComplete = true;
			m_data->gd_err = 0.0;;
		}

		m_data->gdconc = gdc;
	}

	m_data->turnOnPump = true;	//must change water
	m_data->turnOnLED = *ic;
	m_data->calibrationName = m_data->calibrationBase + "_" + *ic;
}