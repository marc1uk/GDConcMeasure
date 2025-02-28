#ifndef Scheduler_H
#define Scheduler_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/chrono.hpp"

class Scheduler: public Tool
{
	public:
		Scheduler();

		bool Initialise(std::string configfile,DataModel &data);
		bool Execute();
		bool Finalise();

		bool DepleteWater();
		bool CirculateWater();
		bool Calibrate();
		bool IsCalibrated();
		bool IsCalibrationDone();
		bool Measure();
		bool IsMeasurementDone();

		boost::posix_time::ptime Wait(double t = 0.0);

	private:

		boost::posix_time::ptime last;

		int verbose;

		int idle_time;
		int power_up_time;
		int power_down_time;
		int change_water_time;
		int settle_water_time;

		int rest_time;

		state nextState, lastState;

		std::map<state, std::string> stateName;

};


#endif
