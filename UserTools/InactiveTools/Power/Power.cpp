#include "Power.h"

Power::Power() : Tool()
{
}

bool Power::Initialise(std::string configfile, DataModel &data)
{
	if(configfile!="")
		m_variables.Initialise(configfile);
	//m_variables.Print();

	m_data= &data;

	m_variables.Get("wake_pin",powerpin);
	
	std::stringstream command;
	command<<"echo \""<<powerpin<<"\" > /sys/class/gpio/export";
	system(command.str().c_str());

	command.str("");
	command<<"echo \"out\" > /sys/class/gpio/gpio"<<powerpin<<"/direction";
	system(command.str().c_str());

	TurnOff();
	sleep(1);	
	return true;
}

bool Power::Execute()
{
	switch (m_data->mode)
	{
		case state::power_up:
			TurnOn();
			//sleep(4);	
			//sleep(10);
			break;
		case state::power_down:
			TurnOff();
			break;
	}
	return true;
}


bool Power::Finalise()
{
	TurnOff();
	return true;
}

void Power::TurnOn()
{
	//write to GPIO
 	std::stringstream command;
	command<<"echo \"0\" > /sys/class/gpio/gpio"<<powerpin<<"/value";
	system(command.str().c_str());

}

void Power::TurnOff()
{
	//write to GPIO
 	std::stringstream command;
	command<<"echo \"1\" > /sys/class/gpio/gpio"<<powerpin<<"/value";
	system(command.str().c_str());


}
