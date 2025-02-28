#include "Scheduler.h"

Scheduler::Scheduler():Tool(){}

bool Scheduler::Initialise(std::string configfile, DataModel &data)
{
	if(configfile!="")
		m_variables.Initialise(configfile);
	//m_variables.Print();

	m_data = &data;

	m_data->mode = state::idle;

	m_variables.Get("verbose", verbose);

	m_variables.Get("idle",		idle_time);
	m_variables.Get("power_up",	power_up_time);
	m_variables.Get("power_down",	power_down_time);
	m_variables.Get("change_water", change_water_time);
	m_variables.Get("settle_water", settle_water_time);

	if (idle_time < settle_water_time)
		idle_time = settle_water_time;

	last = boost::posix_time::second_clock::local_time();

	stateName[state::idle]			= "idle";
	stateName[state::power_up]              = "power_up";
	stateName[state::power_down]            = "power_down";
	stateName[state::init]                  = "init";
	stateName[state::calibration]           = "calibration";
	stateName[state::calibration_done]      = "calibration_done";
	stateName[state::measurement]           = "measurement";
	stateName[state::measurement_done]      = "measurement_done";
	stateName[state::finalise]              = "finalise";
	stateName[state::take_dark]		= "take_dark";
	stateName[state::take_spectrum]		= "take_spectrum";
	stateName[state::analyse]		= "analyse";
	stateName[state::change_water]		= "change_water";
	stateName[state::settle_water]		= "settle_water";
	stateName[state::manual_on]		= "manual_on";
	stateName[state::manual_off]		= "manual_off";

	nextState = state::idle;
	lastState = state::idle;
	rest_time = 0;

	m_data->isCalibrationTool = false;
	m_data->isMeasurementTool = false;
	
	std::cout << "Size GdTree " << m_data->SizeGdTree() << std::endl;
	m_data->m_gdtrees.clear();
	std::cout << "Size GdTree " << m_data->SizeGdTree() << std::endl;
	//HACK
	//m_data->isCalibrated = true;

	return true;
}

bool Scheduler::Execute()
{
	std::cout << "mode start = " << stateName[m_data->mode] << std::endl;

	switch (m_data->mode)
	{
		case state::idle:
			last = Wait(rest_time);			//ToolDAQ is put to sleep until next measurement
			m_data->mode = state::power_up;		//time to power up the PSU!	in power up, Power tool turns on PSU
			break;
		case state::power_up:				//PSU has been turned on
			last = Wait(power_up_time);		//small wait to let PSU be stable
			m_data->mode = state::init;		//init the other Tools
			break;
		case state::init:				//ToolChain has been configured
			last = Wait();
			if (Calibrate() && Measure())
			{
				if (IsCalibrated())
					m_data->mode = state::measurement;
				else
					m_data->mode = state::calibration;	//if there isn't, this one must be done first!
			}
			else if (Calibrate() && !Measure())
			{
				if (IsCalibrated())
					m_data->mode = state::finalise;
				else
					m_data->mode = state::calibration;	//if there isn't, this one must be done first!
			}
			else if (!Calibrate() && !Measure())	//turn on LED?
				m_data->mode = manual_on;
			else
			{
				std::cout << "No Calibration tool, I have to quit" << std::endl;
				m_data->mode = state::finalise;
			}
			break;
		case state::calibration:
			last = Wait();

			if (DepleteWater())	//override cause water must be removed
			{
				nextState = state::calibration;
				m_data->mode = state::change_water;
				break;
			}

			if (CirculateWater())	//override cause water must be renewd
			{
				if (IsCalibrationDone())
					lastState = state::calibration_done;
				else
					lastState = state::calibration;
				nextState = state::take_dark;
				m_data->mode = state::change_water;
			}
			else 
			{
				if (IsCalibrationDone())
					lastState = state::calibration_done;
				else
					lastState = state::calibration;
				m_data->mode = state::analyse;	//turn LED on
			}
			break;
		case state::calibration_done:
			last = Wait();
			m_data->mode = state::finalise;	//calibration done and saved, finalise
			break;
		case state::measurement:
			last = Wait();
			if (IsMeasurementDone())			//check if calibration is completed
				lastState = state::measurement_done;	//if so, it can be saved to disk
			else
				lastState = state::measurement;	//if not, repeat calibration loop
			m_data->mode = state::take_dark;	//turn LED on
			break;
		case state::measurement_done:
			last = Wait();
			m_data->mode = state::finalise;	//measurement done, turn on pump and change water for next measurement and finalise
			break;
		//////////////////MINI LOOP FOR ANALYSIS	//////////
		//take dark rate and save it
		//take spectrum (LED is on) and average it
		//analyse data collected
		//return to current mode, stored in nextState
		case state::take_dark:	//spectrum taken ->  turn off led, analyse data
			m_data->mode = state::take_spectrum;
			break;
		case state::take_spectrum:	//spectrum taken ->  turn off led, analyse data
			m_data->mode = state::analyse;
			break;
		case state::analyse:	//go back to state (calibration or measurement)
			m_data->mode = lastState;
			break;
		////////////////END OF MINI LOOP
		case state::change_water:			//only used for calibration
			last = Wait(change_water_time);			//wait for pump to complete if needed
			if (DepleteWater())
				m_data->mode = nextState;
			else
				m_data->mode = state::settle_water;
			break;
		case state::settle_water:
			last = Wait(settle_water_time);			//wait for water to settle if needed
			m_data->mode = nextState;
			break;
		case state::finalise:
			last = Wait(change_water_time);			//wait for pump to complete if needed
			m_data->mode = state::power_down;		//turn off PSU
			break;
		case state::power_down:
			last = Wait(power_down_time);			//small wait to let PSU be down
			rest_time = idle_time;
			m_data->mode = state::idle;			//and move to idle
			break;
			//////////
			/////special states for LED 
		case state::manual_on:
			m_data->mode = state::manual_off;
			break;
		case state::manual_off:
			m_data->mode = state::power_down;
			break;
	}

	std::cout << "mode end = " << stateName[m_data->mode] << std::endl;

	return true;
}

bool Scheduler::Finalise()
{
	return true;
}

bool Scheduler::DepleteWater()
{
	return m_data->depleteWater;
}

bool Scheduler::CirculateWater()
{
	return m_data->circulateWater;
}

bool Scheduler::Calibrate()
{
	return m_data->isCalibrationTool;
}

bool Scheduler::IsCalibrated()
{
	return m_data->isCalibrated;
}

bool Scheduler::IsCalibrationDone()
{
	return m_data->calibrationDone;
}

bool Scheduler::Measure()
{
	return m_data->isMeasurementTool;
}

bool Scheduler::IsMeasurementDone()
{
	return m_data->measurementDone;
}

//check if enough time has passed since last status change
//wait remaining time so that in total t seconds have passed since last
//with t = 0, no wait happens and current time is passed
boost::posix_time::ptime Scheduler::Wait(double t)
{
	std::cout << "waiting " << t << std::endl;

	boost::posix_time::ptime current(boost::posix_time::second_clock::local_time());
	boost::posix_time::time_duration period(0, 0, t, 0);
	boost::posix_time::time_duration lapse(period - (current - last));

	if (!lapse.is_negative())			//if positive, wait!
		usleep(lapse.total_microseconds());

	return boost::posix_time::second_clock::local_time();
}
