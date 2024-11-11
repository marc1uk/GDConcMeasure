#ifndef CorrectStepChanges_H
#define CorrectStepChanges_H

#include <string>
#include <iostream>

#include "Tool.h"

class CorrectStepChanges: public Tool {
	public:
	
	CorrectStepChanges();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	private:
	
	double stddev_tolerance = 0.0001;
	double conc_change_tolerance = 0.003;
	int sample_size = 5;
	int ignored_pts = 1;
	
	std::deque<double> gdconcs_set1;
	std::deque<double> gdconcs_set2;
	std::deque<bool> step_vetoed;
	std::deque<bool> step_detected;
	std::deque<double> step_applied;
	std::deque<double> running_step;
	std::deque<std::string> lednames;
	std::deque<int> measurementnum;
	std::deque<std::string> timestamps;
	double new_conc;
	double running_change=0;
	int thismeasurementnum=-1;
	bool new_meas=false;
	
	std::deque<std::pair<double,double>> corrected_set1;
	std::deque<std::pair<double,double>> corrected_set2;
	
	double p_critical;
	
	int get_ok;
	int m_verbose=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	
};


#endif
