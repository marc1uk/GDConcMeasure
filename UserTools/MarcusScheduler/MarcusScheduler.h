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
	int loop_count=0;                   // flattened loop iteration counter, used for unique output file naming
	std::vector<int> loop_starts;       // start command number for loops
	std::vector<int> loop_iterations;   // total iterations for each loop
	std::vector<int> loop_counts;       // number performed so far
	std::string break_loop_flagfile_name;
	
	std::map<std::string,int> LED_states{{"R",0}, {"G",0}, {"B",0}, {"White",0}, {"385",0}, {"260",0}, {"275",0}};
	const std::map<std::string,int> off_LED_states{{"R",0}, {"G",0}, {"B",0}, {"White",0}, {"385",0}, {"260",0}, {"275",0}};
	
	
	void ReadCommandFile();
	void MainMenu();
	void PutSystemInSafeState();
	
	// handle all the automation commands
	void ProcessCommand(std::string& the_command);
	
	// individual handlers for automation commands
	void DoSave(std::string the_command);
	void DoAnalyse(std::string the_command);
	void DoPower(std::string the_command);
	void DoValves(std::string the_command);
	void DoMeasure(std::string the_command);
	void DoWait(std::string the_command);
	void SimpleWaitForDuration(std::string wait_string);
	void WaitForDuration(std::string wait_string);
	void StartLoop(std::string the_command);
	void EndLoop(std::string the_command);
	
	// helper functions
	bool check_break_loop(std::string& the_command);
	int check_for_break_file();
	
	// for logging
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;

};


#endif
