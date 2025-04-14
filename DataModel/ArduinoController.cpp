#ifndef ARDUINO
#include "ArduinoController.h"
#include <vector>

bool ArduinoController::Configure(std::string device, int baud, int verb, Logging* logger){
	m_log = logger;
	verbosity = verb;
	return arduino.Configure(device, baud, verb, logger);
}

bool ArduinoController::Help(){
	// FIXME this isn't very helpful as commands used by the arduino are different to the keys
	// set into the CStore to trigger the corresponding functions. TODO align them
	std::string response;
	get_ok = arduino.SendAndReceive("HELP",response);
	std::cout<<response<<std::endl;
	return get_ok;
}

bool ArduinoController::SetState(const std::string& name, bool enable){
	Log(m_unique_name+" SetState("+name+", "+std::to_string(enable)+")",v_debug,verbosity);
	if(enable) return Enable(name);
	return Disable(name);
}

bool ArduinoController::Disable(const std::string& name){
	std::string resp;
	return arduino.SendAndReceive(name+" DISABLE", resp);
}

bool ArduinoController::Enable(const std::string& name){
	std::string resp;
	return arduino.SendAndReceive(name+" ENABLE", resp);
}

bool ArduinoController::Dark(){
	// turn off lights
	std::string resp;
	get_ok = arduino.SendAndReceive("DARK",resp);
	return get_ok; // && (resp=="Disabling all lights");
}

bool ArduinoController::ShutItDown(){
	
	/*
	// turn off lights
	std::map<std::string, bool> results;
	results["LED275"] = Disable("LED275");
	results["White"] = Disable("WHITE");
	results["Deuterium"] = Disable("DEUTERIUM");
	results["Tungsten"] = Disable("TUNGSTEN");
	// close shutters
	results["GAD"] = Disable("GAD_ARM");
	results["REF"] = Disable("REF_ARM");
	results["LAMP"] = Disable("LAMP_SHUTTER");
	// disable lamp DB15 control
	results["LAMP_DB15"] = Disable("LAMP_DB15");
	// disable other relays, whatever they may be used for
	results["Relay2"] = Disable("RELAY2"); // pump, perhaps.
	results["Relay3"] = Disable("RELAY3"); // unused.
	// close valves
	// FIXME ah! but while we're using the pump, we need at least one set of valves open!
	// for such situations comment out the below to not automatically disable them.
	results["PARALLEL"] = Disable("PARALLEL");
	results["TUBE"] = Disable("TUBE");
	
	bool allok=true;
	for(auto&& res : results){
		if(!res.second){
			allok = false;
			Log(m_unique_name+"::ShutItDown error disabling "+res.first,v_error,verbosity);
		}
	}
	*/
	
	std::string resp;
	bool allok = arduino.SendAndReceive("OFF",resp);
	
	return allok; // && resp=="Turning everything off");
}

bool ArduinoController::GetLedTemp(double& led_temp){
	
	led_temp = -999; // default
	
	std::string resp;
	bool ok = arduino.SendAndReceive("LED_TEMP",resp);
	if(!ok){
		Log(m_unique_name+"::GetLedTemp error '"+resp+"'",v_error,verbosity);
		return false;
	}
	
	// format is 'LED_TEMP: X' (so we know we're getting the right number)
	std::stringstream ss(resp);
	std::string tmp;
	ss >> tmp >> led_temp;
	ok = (!ss.fail() && (ss>>std::ws).eof());
	if(!ok || tmp!="LED_TEMP:"){
		Log(m_unique_name+"::GetLedTemp error '"+resp+"'",v_error,verbosity);
		return false;
	}
	
	return true;
}

bool ArduinoController::GetSolTemp(int sol_num, double& sol_temp){
	
	sol_temp=-999; // default
	
	// solenoid numbers should be 0-2
	if(sol_num <0 || sol_num > 2){
		Log(m_unique_name+"::GetSolTemp invalid solenoid num "+std::to_string(sol_num),v_error,verbosity);
		return false;
	}
	
	std::string key = "SOL"+std::to_string(sol_num)+"_TEMP";
	std::string resp;
	bool ok = arduino.SendAndReceive(key,resp);
	if(!ok){
		Log(m_unique_name+"::GetSolTemp error '"+resp+"'",v_error,verbosity);
		return false;
	}
	
	// response format is 'SOLN_TEMP: X' (so we know we're parsing the response string)
	std::string expected_key = "SOL"+std::to_string(sol_num)+"_TEMP:";
	std::stringstream ss(resp);
	std::string tmp;
	ss >> tmp >> sol_temp;
	ok = (!ss.fail() && (ss>>std::ws).eof());
	if(!ok || tmp!=expected_key){
		Log(m_unique_name+"::GetSolTemp error '"+resp+"'",v_error,verbosity);
		return false;
	}
	return true;
}

bool ArduinoController::GetSolTemps(std::array<double,3>& temps){
	
	for(int i=0; i<3; ++i) temps.at(i) = -999; // default
	
	// more likely we're gonna want all 3, so we have a command for that
	std::string resp;
	bool ok = arduino.SendAndReceive("SOL_TEMPS",resp);
	if(!ok){
		Log(m_unique_name+"::GetSolTemps error '"+resp+"'",v_error,verbosity);
		return false;
	}
	// response should be of the form: 'SOL_TEMPS: T1,T2,T3'
	std::stringstream ss(resp);
	std::string tmp;
	int temp_i=0;
	while(std::getline(ss, tmp, ',')){
		Log(m_unique_name+": next part: '"+tmp+"'",v_debug,verbosity);
		if(tmp=="SOL_TEMPS:") continue;
		if(temp_i>2){
			Log(m_unique_name+"::GetSolTemps too many temperatures in response? '"+resp+"'",v_error,verbosity);
			return false;
		}
		try{
			double temp = std::stod(tmp);
			temps.at(temp_i) = temp;
		} catch(...){
			Log(m_unique_name+"::GetSolTemps error; bad temp format '"+resp+"'",v_error,verbosity);
			return false;
		}
		++temp_i;
	}
	return true;
}

bool ArduinoController::GetFlowRate(double& flow_rate){
	
	flow_rate=-999; // default
	
	std::string resp;
	bool ok = arduino.SendAndReceive("FLOW_SENSE",resp);
	if(!ok){
		Log(m_unique_name+"::GetFlowStatus error '"+resp+"'",v_error,verbosity);
		return false;
	}
	// format is 'FLOW_SENSE: X' (so we know we're getting the right number)
	// X is a flow rate in revolutions/sec
	std::stringstream ss(resp);
	std::string tmp;
	ss >> tmp >> flow_rate;
	ok = (!ss.fail() && (ss>>std::ws).eof());
	if(!ok || tmp!="FLOW_SENSE:"){
		Log(m_unique_name+"::GetFlowRate error '"+resp+"'",v_error,verbosity);
		return false;
	}
	Log(m_unique_name+"::GetFlowRate sensor value: "+std::to_string(flow_rate),v_debug,verbosity);
	return true;
}

bool ArduinoController::GetLeakStatus(int& leak_sensor_val){
	std::string resp;
	bool ok = arduino.SendAndReceive("LEAK_CHECK",resp);
	if(!ok){
		Log(m_unique_name+"::GetLeakStatus error '"+resp+"'",v_error,verbosity);
		return false;
	}
	// format is a voltage. 1023 is totally dry, generally any appreciable amount
	// of water on any one sensor drops it to <450, and it bottoms out around 330
	// even if all sensors are soaked. So a threshold of 700 should be a
	// pretty definitive test of water
	std::stringstream ss(resp);
	std::string tmp;
	ss >> tmp >> leak_sensor_val;
	ok = (!ss.fail() && (ss>>std::ws).eof());
	bool good_header = (tmp=="LEAK_CHECK:");
	if(!ok || !good_header){
		Log(m_unique_name+"::GetLeakStatus error '"+resp+"'",v_error,verbosity);
		return false;
	}
	Log(m_unique_name+"::GetLeakStatus sensor value: "+std::to_string(leak_sensor_val),v_debug,verbosity);
	std::string status_msg = (leak_sensor_val<700) ? "WARNING: GOT WATER!" : "OK: no water";
	//std::cout<<"Leak_Status: "<<status_msg<<std::endl;
	return true;
}
#endif
