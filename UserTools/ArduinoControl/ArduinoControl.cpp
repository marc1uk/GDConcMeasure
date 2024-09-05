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
	
	// serial settings
	baud_rate=19200;
	m_variables.Get("baud_rate", baud_rate);
	std::string com_port_str="";  // e.g. '/dev/ttyS0'
	get_ok = m_variables.Get("com_port", com_port_str);
	if(!get_ok){
		Log(m_unique_name+" No com port string specified!",v_error,verbosity);
		return false;
	}
	com_port = DeviceNameToNumber(com_port_str);
	if(com_port<0){
		Log(m_unique_name+" invalid com port '"+com_port_str+"'",v_error,verbosity);
		return false;
	}
	Log(m_unique_name+" com port '"+com_port_str+"' mapped to port #"+std::to_string(com_port),v_debug,verbosity);
	get_ok = Connect();
	if(!get_ok){
		Log(m_unique_name+" failed to connect to arduino!",v_error,verbosity);
		return false;
	}
	
	return true;
}


bool ArduinoControl::Execute(){
	
	get_ok = true;
	
	// TODO arduino control is bizarrely slow. Either find out why,
	// or re-write this to combine the set of new states into one command
	// so that we can set everything with a single write.
	
	// check for LED state changes
	// TODO sanity check that 'tmp' is "1" or "0"
	std::string tmp;
	if(m_data->CStore.Get("LED",tmp) && tmp=="Change"){
		// check Dark first as it turns everything off, then turn on what we want.
		if(m_data->CStore.Get("Dark",tmp)){
			std::string resp;
			get_ok = SendAndReceive("OFF",resp);
			get_ok = get_ok & resp=="Disabling all lights";
			m_data->CStore.Remove("Dark");
		}
			
		if(m_data->CStore.Get("White",tmp)){
			get_ok = get_ok && SetState("WHITE", (tmp=="1"));
			m_data->CStore.Remove("White");
		}
		
		if(m_data->CStore.Get("275_A",tmp)){
			get_ok = get_ok && SetState("LED275", (tmp=="1"));
			m_data->CStore.Remove("275_A");
		}
		
		/* FIXME so we can analyse */
		if(m_data->CStore.Get("275_B",tmp)){
			get_ok = get_ok && SetState("RELAY3", (tmp=="1"));
			m_data->CStore.Remove("275_B");
		}
		
		if(m_data->CStore.Get("Deuterium",tmp)){
			get_ok = get_ok && SetState("DEUTERIUM", (tmp=="1"));
			m_data->CStore.Remove("Deuterium");
		}
		
		if(m_data->CStore.Get("Tungsten",tmp)){
			get_ok = get_ok && SetState("TUNGSTEN", (tmp=="1"));
			m_data->CStore.Remove("Tungsten");
		}
		
		// generic Grove relay control
		for(int i=0; i<5; ++i){
			std::string key="Relay"+std::to_string(i);
			if(m_data->CStore.Get(key,tmp)){
				if(i<4) get_ok = get_ok && SetState(key, (tmp=="1"));
				else Log(m_unique_name+" Relay4 received; relays are numbered 0-3",v_error, verbosity);
				m_data->CStore.Remove(key);
			}
		}
	}
	m_data->CStore.Remove("LED");
	
	if(m_data->CStore.Get("Shutter_gad",tmp)){
		get_ok = get_ok && SetState("GAD_ARM", (tmp=="OPEN"));
		m_data->CStore.Remove("Shutter_gad");
	}
	
	if(m_data->CStore.Get("Shutter_ref",tmp)){
		get_ok = get_ok && SetState("REF_ARM", (tmp=="OPEN"));
		m_data->CStore.Remove("Shutter_ref");
	}
	
	if(m_data->CStore.Get("Shutter_lamp",tmp)){
		get_ok = get_ok && SetState("LAMP_SHUTTER", (tmp=="OPEN"));
		m_data->CStore.Remove("Shutter_lamp");
	}
	
	if(m_data->CStore.Get("Valve_gad",tmp)){
		get_ok = get_ok && SetState("TUBE", (tmp=="OPEN"));
		m_data->CStore.Remove("Valve_gad");
	}
	
	if(m_data->CStore.Get("Valve_parallel",tmp)){
		get_ok = get_ok && SetState("PARALLEL", (tmp=="OPEN"));
		m_data->CStore.Remove("Valve_parallel");
	}
	
	// this just enables/disables control via the DB15 connector on the back of the lamp
	if(m_data->CStore.Get("Lamp_DB15",tmp)){
		get_ok = get_ok && SetState("LAMP_DB15", (tmp=="ENABLE"));
		m_data->CStore.Remove("Lamp_DB15");
	}
	
	return get_ok;
}


bool ArduinoControl::Finalise(){
	
	bool allok = ShutItDown();
	allok = allok && Disconnect();
	
	return allok;
	
}

bool ArduinoControl::SetState(const std::string& name, bool enable){
	Log(m_unique_name+" SetState("+name+", "+std::to_string(enable)+")",v_debug,verbosity);
	if(enable) return Enable(name);
	return Disable(name);
}

bool ArduinoControl::Disable(const std::string& name){
	std::string resp;
	return SendAndReceive(name+" DISABLE", resp);
}

bool ArduinoControl::Enable(const std::string& name){
	std::string resp;
	return SendAndReceive(name+" ENABLE", resp);
}

bool ArduinoControl::ShutItDown(){
	
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
	results["Relay2"] = Disable("RELAY2");
	results["Relay3"] = Disable("RELAY3");
	// close valves
	// FIXME ah! but while we're using the pump, we need at least one set of valves open!
	// for the time being, do not automatically disable them.
	//results["PARALLEL"] = Disable("PARALLEL");
	//results["TUBE"] = Disable("TUBE");
	
	bool allok=true;
	for(auto&& res : results){
		if(!res.second){
			allok = false;
			Log(m_unique_name+"::ShutItDown error disabling "+res.first,v_error,verbosity);
		}
	}
	
	return allok;
}

bool ArduinoControl::Connect(){
	
	// default arduino serial is 8 data bits, 1 stop bit, no parity, no flow control
	char mode[]={'8','N','1',0};
	int flowctl=0;
	
	if(RS232_OpenComport(com_port, baud_rate, mode, flowctl)){
		Log(m_unique_name+" Unable to open com port "+std::to_string(com_port),v_error,verbosity);
		return false;
	}
	
	// make sure nothing stray in the buffers
	// wait for arduino to do initial setup - this prints things we can't suppress
	// (from relay i2c library)
	std::this_thread::sleep_for(std::chrono::milliseconds(1000));
	RS232_flushRXTX(com_port);
	
	// check comms are working as expected
	std::string response;
	get_ok = SendAndReceive("Hello",response);
	if(!get_ok || response!="Hello!"){
		if(get_ok) Log(m_unique_name+" Error with initial comms check; expected response 'Hello!', got '"+response+"'",v_error,verbosity);
		return false;
	}
	
	// request arduino to connect to the grove relay controller over SPI
	get_ok = SendAndReceive("GROVE",response);
	if(!get_ok) return false;
	
	// turn off all lights, close all shutters and valves.
	//get_ok = ShutItDown();
	
	return get_ok;
	
}

bool ArduinoControl::Disconnect(){
	// tell arduino to close serial comms
	// (it'll open it again in 5 seconds and start listening for new connections)
	get_ok = SerialWrite("QUIT");
	if(!get_ok){
		Log(m_unique_name+" Error sending message 'QUIT'",v_error,verbosity);
		return false;
	}
	return true;
}

bool ArduinoControl::SendAndReceive(std::string msg, std::string& response, int timeout){
	Log(m_unique_name+" SendAndReceive of command '"+msg+"'",v_debug,verbosity);
	// send command
	//std::cout<<"Writing command"<<std::endl;
	get_ok = SerialWrite(msg);
	if(!get_ok){
		Log(m_unique_name+" Error sending message '"+msg+"'",v_error,verbosity);
		return false;
	}
	// read response
	//std::cout<<"Reading response..."<<std::endl;
	response = SerialRead(timeout);
	if(response.empty()){
		Log(m_unique_name+" Error: No response to command '"+msg+"'",v_error,verbosity);
		return false;
	}
	// as a standard, the arduino sketch returns strings starting with 'Err: ABC' when it fails to do something
	if(response.substr(0,4)=="Err:"){
		Log(m_unique_name+" Error: command '"+msg+"' returned '"+response+"'",v_error,verbosity);
		return false;
	}
	
	return true;
}

bool ArduinoControl::SerialWrite(std::string msg){
	if(msg.size()==0){
		Log(m_unique_name+" SerialWrite with empty message!",v_error,verbosity);
		return false;
	}
	if(msg.back()!='\n') msg.push_back('\n');
	get_ok = RS232_SendBuf(com_port, (unsigned char*)(msg.data()), msg.size());
	if(get_ok!=msg.size()){
		Log(m_unique_name+" Error writing message '"+msg+"'; sent "+std::to_string(get_ok)
		    +" bytes of "+std::to_string(msg.size()),v_error,verbosity);
		return false;
	}
	return true;
}

std::string ArduinoControl::SerialRead(int timeout_ms){
	int bytesread=0;
	const int maxbytes = 4096; // don't receive more than this many bytes
	std::string response;
	unsigned char buf[4096];
	atimer.reset();
	while(atimer.ms_elapsed()<timeout_ms){
		bytesread = RS232_PollComport(com_port, buf, 4096);
		if(bytesread>0){
			//printf("received %d bytes\n",bytesread);
			//for(int i=0; i<bytesread; ++i) printf("'%c' ",buf[i]);  // "%X: %c",buf[i],buf[i]);
			//printf("\n");
			response += std::string((const char*)(buf),bytesread);
		}
		// arduino ends SerialWrite with CRLF ('0x0D 0x0A') ('\n' == LF)
		if(bytesread>0 && (buf[bytesread-1]=='\n' || buf[bytesread-1]=='\0')){
			// received termination char
			//printf("got termchar!\n");
			break;
		} else if(response.size()>maxbytes){
			// don't just keep adding indefinitely.. something wrong.
			// instead return what we got and flush the remainder of the buffer.
			RS232_flushRXTX(com_port);
			Log(m_unique_name+" Warning! Too many chars in serial buffer! Flushing any remainder",v_warning,verbosity);
			break;
		}
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	}
	//printf("total response: '%s'\n",response.c_str());
	// trim whitespace
	while(response.size() && (isspace(response.back()) || response.back()=='\n')) response.pop_back();
	//printf("trimmed response: '%s'\n",response.c_str());
	return response;
}

int ArduinoControl::DeviceNameToNumber(std::string name){
	static const std::vector<std::string> devicenames{
	 "/dev/ttyS0","/dev/ttyS1","/dev/ttyS2","/dev/ttyS3","/dev/ttyS4","/dev/ttyS5",
	 "/dev/ttyS6","/dev/ttyS7","/dev/ttyS8","/dev/ttyS9","/dev/ttyS10","/dev/ttyS11",
	 "/dev/ttyS12","/dev/ttyS13","/dev/ttyS14","/dev/ttyS15","/dev/ttyUSB0",
	 "/dev/ttyUSB1","/dev/ttyUSB2","/dev/ttyUSB3","/dev/ttyUSB4","/dev/ttyUSB5",
	 "/dev/ttyAMA0","/dev/ttyAMA1","/dev/ttyACM0","/dev/ttyACM1",
	 "/dev/rfcomm0","/dev/rfcomm1","/dev/ircomm0","/dev/ircomm1",
	 "/dev/cuau0","/dev/cuau1","/dev/cuau2","/dev/cuau3",
	 "/dev/cuaU0","/dev/cuaU1","/dev/cuaU2","/dev/cuaU3"
	 };
	 auto it = std::find(devicenames.begin(),devicenames.end(),name);
	if(it!=devicenames.end()) return std::distance(devicenames.begin(),it);
	else return -1;
}
