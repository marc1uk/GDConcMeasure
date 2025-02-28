#include "Writer.h"

Writer::Writer() : Tool()
{
}

bool Writer::Initialise(std::string configfile, DataModel &data)
{
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();

	m_data= &data;

	return true;
}

bool Writer::Execute()
{
	Loop();

	switch (m_data->mode)
	{
		case state::init:
			break;
		case state::analyse:
		case state::calibration_done:
			WriteCalibration();
			break;
		case state::measurement_done:
			WriteMeasurement();
			break;
		default:
			break;
	}

	return true;
}

bool Writer::Finalise()
{
	Loop();
	return true;
}

/* this method saves calibration trees in the calibration output
 * givin as a name hte calibration base element
 * it is arranged in sub-directories with as name "d_" + calibration name + timestamp
 * the sub-directories will contain the calibration and the fitting functions, if present
 * the method creates the subdir, unless they exist already
 */
void Writer::WriteCalibration()
{
	/*
	//	    012345678901234
	//format is YYYYMMDDTHHMMSS (there shouldn't be fractional seconds
	boost::posix_time::ptime current(boost::posix_time::second_clock::local_time());
	std::string tc = boost::posix_time::to_iso_string(current);
	*/

	std::cout << "OPENING " << m_data->calibrationFile << std::endl;
	TFile file(m_data->calibrationFile.c_str(), "UPDATE");

	std::map<std::string, GdTree*>::iterator it = m_data->m_gdtrees.begin();
	for (; it != m_data->m_gdtrees.end(); ++it)
	{
		if (it->first.find(m_data->calibrationBase) != 0)	//string starts with base_name
			continue;

		std::string dirname = it->first;
		dirname.replace(0, m_data->calibrationBase.length(), "d");
		dirname += "_" + m_data->calibrationTime;

		if (!file.GetDirectory(dirname.c_str()))	//no sub-dir with this name
			file.mkdir(dirname.c_str());

		//cd into directory
		file.cd(dirname.c_str());

		//write Tree
		//it->second->Write();
		it->second->Write(m_data->calibrationBase);

		int stat = m_data->concentrationFit;
		//if concentration, write also functions
		if ((it->first == m_data->concentrationName) &&
				stat >= 0)
				//m_data->concentrationFunction &&
				//m_data->concentrationFuncStat &&
				//m_data->concentrationFuncSyst &&
				//stat >= 0)
		{
			m_data->concentrationFunction->Write("", TObject::kWriteDelete);
			m_data->concentrationFuncStat->Write("", TObject::kWriteDelete);
			m_data->concentrationFuncSyst->Write("", TObject::kWriteDelete);
			m_data->concentrationFit->Write(m_data->fitFuncName.c_str(), TObject::kWriteDelete);
		}

		//if (it->first.find(m_data->calibrationBase) == 0)
		//	it->second->Write();
	}

	file.Close();
}

/* this method saves measurement trees in the measurment output
 * giving as a name the measurement base element
 * it is arranged in sub-directories with as name "d_" + measurement name
 * the sub-directories will contain also the calibration used for that measurement
 * the method creates the subdir, unless they exist already
 */
void Writer::WriteMeasurement()
{
	TFile file(m_data->measurementFile.c_str(), "UPDATE");

	std::map<std::string, GdTree*>::iterator it = m_data->m_gdtrees.begin();
	for (; it != m_data->m_gdtrees.end(); ++it)
	{
		std::string dirname = it->first;
		std::string name;
		if (dirname.find(m_data->calibrationBase) == 0)
		{
			dirname.replace(0, m_data->calibrationBase.length(), "d");
			name = m_data->calibrationBase;
		}
		else if (dirname.find(m_data->measurementBase) == 0)
		{
			dirname.replace(0, m_data->measurementBase.length(), "d");
			name = m_data->measurementBase;
		}
		else
			continue;

		//std::string dirname = "d" + it->first;

		//if (dirname.find(m_data->calibrationBase) == 1)
		//	dirname.replace(1, m_data->calibrationBase.length(), m_data->measurementBase);

		std::cout << "directory name " << dirname << std::endl;
		if (!file.GetDirectory(dirname.c_str()))	//no sub-dir with this name
			file.mkdir(dirname.c_str());

		//cd into sub-dir
		file.cd(dirname.c_str());

		//write tree
		//it->second->Write();
		it->second->Write(name);
	}

	file.Close();
}

void Writer::Loop()
{
	std::cout << "WRITER -- size GdTree " << m_data->SizeGdTree() << std::endl;
	std::map<std::string, GdTree*>::iterator it = m_data->m_gdtrees.begin();
	for (; it != m_data->m_gdtrees.end(); ++it)
	{
		std::cout << "  ---   " << it->first << ",\t" << it->second;
		std::cout << std::endl;
		//if (it->second)
		//	std::cout << ":\t" << it->second->GetEntries() << std::endl;
		//else
		//	std::cout << std::endl;
	}
}
