#ifndef MarcusScheduler_H
#define MarcusScheduler_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <pthread.h>
#include <sys/stat.h>
#include <unistd.h>

#include "Tool.h"

#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/chrono.hpp"

class MarcusScheduler: public Tool {
	
	public:
	MarcusScheduler();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	private:
	boost::posix_time::time_duration m_period;
	std::vector<std::string> commands;
	int current_command=-1;
	int command_step=0;
	std::string command_file;
	std::string measurement_name;
	int overwrite_saves=1;
	bool looping=false;
	int loop_count=0;
	
	std::string break_loop_flagfile_name="UserTools/MarcusScheduler/breakloop";
	
	std::map<std::string,int> LED_states{{"R",0}, {"G",0}, {"B",0}, {"White",0}, {"385",0}, {"260",0}, {"275",0}};
	const std::map<std::string,int> off_LED_states{{"R",0}, {"G",0}, {"B",0}, {"White",0}, {"385",0}, {"260",0}, {"275",0}};
	
	void WaitForDuration(std::string wait_string);
	bool check_break_loop();
	
};


#endif
