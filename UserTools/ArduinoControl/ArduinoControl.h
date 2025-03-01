#ifndef ArduinoControl_H
#define ArduinoControl_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "DataModel.h"
#include "ArduinoController.h"

class ArduinoControl: public Tool {
	
	public:
	
	ArduinoControl();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	private:
	ArduinoController controller;
	
	std::string m_configfile;
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	int get_ok;
	
};


#endif
