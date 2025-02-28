#include "FakeSpectrometer.h"

FakeSpectrometer::FakeSpectrometer():Tool(){}


bool FakeSpectrometer::Initialise(std::string configfile, DataModel &data){

	if(configfile!="")
		m_variables.Initialise(configfile);
	m_configfile = configfile;
	//m_variables.Print();

	m_data= &data;

	m_variables.Get("noise",	noise);
	m_variables.Get("pedestal",	pedestal);	//file in which calibration is saved
	m_variables.Get("error", 	error);
	//m_variables.Get("base_name",	m_data->measurementBase);	//base_name for calibation

	std::random_device rd;
	mt = std::mt19937(rd());
	ran = std::uniform_real_distribution<>(-noise, noise);
	err = std::uniform_real_distribution<>(-error, error);

	nnn = 500;
	conc = 0.2;

	return true;
}


bool FakeSpectrometer::Execute()
{
	switch (m_data->mode)
	{
		case state::init:	//wake up, connect spectrometer on USB
			Initialise(m_configfile, *m_data);
			GetX();
			break;
		case state::take_dark:
			GetDark();	//measure dark noise on wake up
			break;
		case state::take_spectrum:
			std::cout << "conc " << m_data->gdconc << std::endl;
			if (m_data->gdconc > 0.0)
				GetData(m_data->gdconc+Error());
			else
				GetPure();
			//Wait();
			break;
		case state::measurement:
			//GetData(conc+Error());
			//Wait();
			break;
		case state::finalise:
			Finalise();
			break;
		default:
			break;
	}
	return true;
}


bool FakeSpectrometer::Finalise()
{
	m_data->wavelength.clear();
	m_data->traceCollect.clear();

	return true;
}

void FakeSpectrometer::Wait()
{
	//sleep(5);
}

bool FakeSpectrometer::GetX()
{
	m_data->wavelength.clear();
	for (int i = 0; i < 1000; ++i)
		m_data->wavelength.push_back(X(i));
}

bool FakeSpectrometer::GetDark()
{
	std::cout << "FS get dark" << std::endl;
	m_data->traceCollect.clear(); 
	for (int i = 0; i < nnn; ++i)
	{
		std::vector<double> tmp;
		for (int j = 0; j < 1000; ++j)
			tmp.push_back(Dark(j));

		m_data->traceCollect.push_back(tmp);
	}


	return true;
}

bool FakeSpectrometer::GetPure()
{
	std::cout << "FS get pure" << std::endl;
	m_data->traceCollect.clear(); 
	for (int i = 0; i < nnn; ++i)
	{
		std::vector<double> tmp;
		for (int j = 0; j < 1000; ++j)
			tmp.push_back(Pure(j));

		m_data->traceCollect.push_back(tmp);
	}

	return true;
}

bool FakeSpectrometer::GetData(double cc)
{
	std::cout << "FS get data" << std::endl;
	m_data->traceCollect.clear(); 
	for (int i = 0; i < nnn; ++i)
	{
		std::vector<double> tmp;
		for (int j = 0; j < 1000; ++j)
			tmp.push_back(Data(j, cc));//-noise*10);

		m_data->traceCollect.push_back(tmp);
	}


	return true;
}

double FakeSpectrometer::Dark(int i)
{
	return pedestal + ran(mt);
}

double FakeSpectrometer::Error()
{
	return err(mt);
}

double FakeSpectrometer::Pure(int i)
{
	return PureFunc(i) + Dark(i);
}

double FakeSpectrometer::Data(int i, double cc)
{
	return DataFunc(i, cc) + Dark(i);
}

double FakeSpectrometer::PureFunc(int i)
{
	return 0.1 * std::exp(- pow((X(i) - 400.0)/5.0, 2));
}

double FakeSpectrometer::DataFunc(int i, double cc)
{
	return PureFunc(i) - Absorb(i, cc);
}

double FakeSpectrometer::Absorb(int i, double cc)
{
	return cc * 0.1  * std::exp(- pow((X(i) - 399.0)/0.5, 2)) + 
	       cc * 0.06 * std::exp(- pow((X(i) - 401.0)/0.5, 2));
}

double FakeSpectrometer::X(int i)
{
	return 200.0 + 2.0 * i / 5.0;
}
