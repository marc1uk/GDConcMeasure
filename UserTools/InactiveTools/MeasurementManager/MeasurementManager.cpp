#include "MeasurementManager.h"

MeasurementManager::MeasurementManager() : Tool()
{
}

bool MeasurementManager::Initialise(std::string configfile, DataModel &data)
{
	if (configfile != "")
	{
		m_configfile = configfile;
		m_variables.Initialise(configfile);
	}
	//m_variables.Print();

	m_data= &data;

	m_variables.Get("verbose", verbose);

	m_variables.Get("measurement",	measureFile);	//file with list of calibration/measurement
	m_variables.Get("output",	m_data->measurementFile);	//file in which calibration is saved
	m_variables.Get("base_name",	m_data->measurementBase);	//base_name for calibation
	//m_variables.Get("tree_name", treeName);		//name of output tree

	m_data->isMeasurementTool = true;
	//HACK
	//m_data->calibrationBase = m_data->measurementBase;

	return true;
}


bool MeasurementManager::Execute()
{
	switch (m_data->mode)
	{
		case state::init:	//about to make measurement, check LED mapping
			Initialise(m_configfile, *m_data);
			measureList = LoadMeasurement();

			if (measureList.size())
				m_data->measurementDone = false;

			im = measureList.begin();	//global iterator
			//std::cout << "Set up for " << measureList.size() << " measurements.\n";
			break;
		case state::measurement:		//turn on led set
			if (im != measureList.end())
			{
				Measure();
				++im;
			}

			if (im == measureList.end())
			{
				m_data->measurementDone = true;
				std::cout << "Measurements completed\n";
			}
			break;
		case state::measurement_done:
			break;
		case state::finalise:
			Finalise();
			break;
	}
	return true;
}


bool MeasurementManager::Finalise()
{
	for (im = measureList.begin(); im != measureList.end(); ++im)
		m_data->DeleteGdTree(m_data->measurementBase + "_" + *im);

	measureList.clear();

	return true;
}

std::vector<std::string> MeasurementManager::MeasurementList()
{
	std::vector<std::string> vList;

	std::ifstream cf(measureFile.c_str());
	std::string line;
	while (std::getline(cf, line))
	{
		if (line.find_first_of('#') != std::string::npos)
			line.erase(line.find_first_of('#'));

		if (line.empty())
			continue;

		//replace commas into spaces
		std::replace(line.begin(), line.end(), ',', ' ');
		std::replace(line.begin(), line.end(), ';', ' ');

		std::stringstream ssl(line);
		std::string name, word;

		while (ssl)
		{
			if (!std::isalnum(ssl.peek()) && ssl.peek() != '_')
				ssl.ignore();
			else if (ssl >> word)
			{
				if (word.find("time") == std::string::npos)	//no time
					name += word + "_";
			}
		}

		name.erase(name.size()-1);
		vList.push_back(name);
	}

	return vList;
}

std::vector<std::string> MeasurementManager::LoadMeasurement()
{
	//get list of calibrations
	std::vector<std::string> mList = MeasurementList();

	//open read only calibration file
	TFile f(m_data->measurementFile.c_str(), "OPEN");
	if (f.IsZombie())
		std::cerr << "A new file will be created!\n";
	else
		std::cout << "Opened file " << f.GetName() << std::endl;
	
	//m_data->measurementFile = 0;

	//find the measurement tree
	//
	for (im = mList.begin(); im != mList.end(); ++im)
	{
		TTree *gdt = 0;

		std::string name = m_data->measurementBase + "_" + *im;
		std::string full = "d_" + *im + "/" + m_data->measurementBase;
		//if exists, get it
		//std::cout << "LOOKING " << full << std::endl;
		if (f.IsOpen())
			gdt = static_cast<TTree*>(f.Get(full.c_str()));
		if (gdt)
			std::cout << "LOADED GDT " << gdt->GetEntries() << std::endl;
		//else create new
		
		//std::cout << "\t--> " << name << std::endl;
		GdTree *measure = gdt ? new GdTree(gdt) : new GdTree(name);

		m_data->AddGdTree(name.c_str(), measure);
	}

	if (f.IsOpen())
		f.Close();

	return mList;
}

bool MeasurementManager::Measure()
{
	m_data->measurementName = m_data->measurementBase + "_" + *im;
	m_data->calibrationName = m_data->calibrationBase + "_" + *im;
	std::cout << "Measure " << *im << " (" << m_data->measurementName << ") with cal " << m_data->calibrationName << std::endl;

	//HACK
	//std::cout << "Ready?\n";
	//std::cin.get();
	m_data->turnOnLED = *im;
	//m_data->turnOnLED = "";
}
