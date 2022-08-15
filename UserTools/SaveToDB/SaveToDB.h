#ifndef SaveToDB_H
#define SaveToDB_H

#include <string>
#include <iostream>

#include "Tool.h"

class SaveToDB: public Tool {
	
public:
	
	SaveToDB();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
private:
	int max_webpage_records=5;
	
	// one method per Tool to put that Tool's results into the DB
	bool MatthewAnalysis();
	bool MarcusAnalysis();
	bool BenPower();
	bool Valve();
	bool BenLED();
	bool BenSpectrometer();
	bool TraceAverage();
	bool SaveTraces();
	bool MarcusScheduler();
	bool RoutineCalibration();  // placeholder, Tool TODO
	bool MatthewTransparency();
	
	std::vector<std::string> field_names;
	std::vector<std::string> criterions;
	std::vector<char> comparators;
	std::string error_ret;
	std::string query_string;
	int runnum=-1;
	int RunConfig=-1;
	
	// for identifying measurements stored in the database
	boost::posix_time::ptime last_measurement_time;
	int measurementnum=-1;
	
	// for identifying measurements stored in the database
	std::string last_measurement_time;
	int measurementnum=-1;
	
	bool get_ok;
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
};


#endif
