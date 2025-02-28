#include "ArduinoControl.h"
#include "rs232.h"

ArduinoControl::ArduinoControl():Tool(){}

bool ArduinoControl::Initialise(std::string configfile, DataModel &data){
	
	m_data= &data;
	
	/* - new method, Retrieve configuration options from the postgres database - */
	int RunConfig=-1;
	m_data->vars.Get("RunConfig",RunConfig);
	
	if(RunConfig>=0){
		std::string configtext;
		bool get_ok = m_data->postgres_helper.GetToolConfig(m_unique_name, configtext);
		if(!get_ok){
			Log(m_unique_name+" Failed to get Tool config from database!",v_error,verbosity);
			return false;
		}
		// parse the configuration to populate the m_variables Store.
		if(configtext!="") m_variables.Initialise(std::stringstream(configtext));
		
	}
	
	/* - old method, read config from local file - */
	if(configfile!="") m_variables.Initialise(configfile);
	
	m_variables.Get("verbosity",verbosity);
	//m_variables.Print();
	
	std::string com_port_str="";  // e.g. '/dev/ttyS0'
	get_ok = m_variables.Get("com_port", com_port_str);
	if(!get_ok){
		Log(m_unique_name+" No com port string specified!",v_error,verbosity);
		return false;
	}
	int baud_rate=19200;
	m_variables.Get("baud_rate", baud_rate);
	
	// configure
	get_ok = controller.Configure(com_port_str, baud_rate, verbosity, m_data->Log);
	if(!get_ok) return false;
	
	// connect to arduino
	get_ok = controller.arduino.Connect();
	if(!get_ok){
		Log(m_unique_name+" failed to connect to arduino!",v_error,verbosity);
		return false;
	}
	
	// put in datamodel for other tools (monitoring tool)
	m_data->arduino = &controller;
	
	return true;
}

bool ArduinoControl::Execute(){
	
	get_ok = true;
	
	// check for LED state changes
	// TODO sanity check that 'tmp' is "1" or "0"
	std::string tmp;
	if(m_data->CStore.Get("LED",tmp) && tmp=="Change"){
		// check Dark first as it turns everything off, then turn on what we want.
		if(m_data->CStore.Get("Dark",tmp)){
			std::string resp;
			get_ok = get_ok && controller.Dark();
			m_data->CStore.Remove("Dark");
		}
			
		if(m_data->CStore.Get("White",tmp)){
			get_ok = get_ok && controller.SetState("WHITE", (tmp=="1"));
			m_data->CStore.Remove("White");
		}
		
		if(m_data->CStore.Get("275_A",tmp)){
			get_ok = get_ok && controller.SetState("LED275", (tmp=="1"));
			m_data->CStore.Remove("275_A");
		}
		
		if(m_data->CStore.Get("Deuterium",tmp)){
			get_ok = get_ok && controller.SetState("DEUTERIUM", (tmp=="1"));
			m_data->CStore.Remove("Deuterium");
		}
		
		if(m_data->CStore.Get("Tungsten",tmp)){
			get_ok = get_ok && controller.SetState("TUNGSTEN", (tmp=="1"));
			m_data->CStore.Remove("Tungsten");
		}
		
		// generic Grove relay control
		for(int i=0; i<5; ++i){
			std::string key="Relay"+std::to_string(i);
			if(m_data->CStore.Get(key,tmp)){
				if(i<4) get_ok = get_ok && controller.SetState(key, (tmp=="1"));
				else Log(m_unique_name+" Relay4 received; relays are numbered 0-3",v_error, verbosity);
				m_data->CStore.Remove(key);
			}
		}
	}
	m_data->CStore.Remove("LED");
	
	if(m_data->CStore.Get("Shutter_gad",tmp)){
		get_ok = get_ok && controller.SetState("GAD_ARM", (tmp=="OPEN"));
		m_data->CStore.Remove("Shutter_gad");
	}
	
	if(m_data->CStore.Get("Shutter_ref",tmp)){
		get_ok = get_ok && controller.SetState("REF_ARM", (tmp=="OPEN"));
		m_data->CStore.Remove("Shutter_ref");
	}
	
	if(m_data->CStore.Get("Shutter_lamp",tmp)){
		get_ok = get_ok && controller.SetState("LAMP_SHUTTER", (tmp=="OPEN"));
		m_data->CStore.Remove("Shutter_lamp");
	}
	
	// this just enables/disables control via the DB15 connector on the back of the lamp
	// note not all electronics boxes support this; some have it hard-wired to 5V
	if(m_data->CStore.Get("Lamp_DB15",tmp)){
		get_ok = get_ok && controller.SetState("LAMP_DB15", (tmp=="ENABLE"));
		m_data->CStore.Remove("Lamp_DB15");
	}
	
	if(m_data->CStore.Get("Valve_gad",tmp)){
		get_ok = get_ok && controller.SetState("TUBE", (tmp=="OPEN"));
		m_data->CStore.Remove("Valve_gad");
	}
	
	if(m_data->CStore.Get("Valve_parallel",tmp)){
		get_ok = get_ok && controller.SetState("PARALLEL", (tmp=="OPEN"));
		m_data->CStore.Remove("Valve_parallel");
	}
	
	if(m_data->CStore.Get("Valve_pump",tmp)){
		get_ok = get_ok && controller.SetState("PUMP", (tmp=="OPEN"));
		m_data->CStore.Remove("Valve_pump");
	}
	
	if(m_data->CStore.Get("BEEP",tmp)){
		// TODO add support for beep patterns
		get_ok = get_ok && controller.arduino.SerialWrite("BEEP");
		m_data->CStore.Remove("BEEP");
	}
	
	return get_ok;
}


bool ArduinoControl::Finalise(){
	
	bool allok = controller.ShutItDown();
	allok = allok && controller.arduino.Disconnect();
	
	return allok;
	
}

