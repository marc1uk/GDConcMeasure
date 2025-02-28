#ifndef ArduinoController_H
#define ArduinoController_H

#include <string>
#include <iostream>
#include <array>

#include "Logging.h"
#include "Arduino.h"

using ToolFramework::Logging;

class ArduinoController {
	
	public:
	bool Configure(std::string com_port, int baud_rate, int verb, Logging* logger);
	
	bool SetState(const std::string& name, bool enable);                     // wrapper around Disable/Enable
	bool Disable(const std::string& name);                                   // S&R(N, "DISABLE"), N = LED275, White, RelayN, Deuterium, Tungsten,
	bool Enable(const std::string& name);                                    // S&R(N, "DISABLE"),     GAD_arm, Ref_arm, Lamp_shutter, Tube, Parallel
	bool Dark(); // turn off all lights                                      // S&R("DARK")
	bool ShutItDown();  // turn off lights, valves, close shutters           // S&R("OFF")
	bool GetLedTemp(double& temp);                                           // S&R("LED_TEMP")
	bool GetSolTemp(int sol, double& temp);                                  // S&R("SOLN_TEMP"), N=0-3
	bool GetSolTemps(std::array<double,3>& temps);                           // S&R("SOL_TEMPS")
	bool GetFlowRate(double& flow_rate);                                     // S&R("FLOW_SENSE")
	bool GetLeakStatus(int& leak_sensor_val);                                // S&R("LEAK_CHECK")
	bool Help(); // fetch and print valid commands from arduino              // S&R("HELP")
	// note Help not super useful as commands to ArduinoController don't match those on arduino. FIXME
	
	Arduino arduino;
	
	private:
	
	Logging* m_log = nullptr;
	template <typename T>
	void Log(T msg, int msg_verb, int class_verb){
		if(m_log){
			m_log->Log(msg, msg_verb,class_verb);
		} else if(msg_verb <= class_verb){  // probably unneeded but y'know
			if(msg_verb==0) std::cerr<<msg<<std::endl;
			else std::cout<<msg<<std::endl;
		}
	}
	
	std::string m_unique_name="ArduinoController";
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	int get_ok;
	
};

#endif
