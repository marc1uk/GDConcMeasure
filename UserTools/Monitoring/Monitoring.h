#ifndef Monitoring_H
#define Monitoring_H

#include <string>
#include <iostream>
#include <chrono>

#include "Tool.h"
#include "DataModel.h"

class Monitoring: public Tool {
	
	public:
	Monitoring();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	private:
	std::chrono::time_point<std::chrono::high_resolution_clock> last_send;
	//std::chrono::duration<std::chrono::seconds> send_period_s;
	int send_period_s = 10; // seconds
	
	std::string m_configfile;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	int get_ok;
	
};


#endif
