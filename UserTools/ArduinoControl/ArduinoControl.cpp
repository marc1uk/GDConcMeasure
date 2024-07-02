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
	com_port=0; // i.e. /dev/ttyS0, equivalent to COM1 on windows
	baud_rate=9600;
	m_variables.Get("com_port", com_port);
	m_variables.Get("baud_rate", baud_rate);
	
	
	
	return true;
}


bool ArduinoControl::Execute(){
	
	get_ok = true;
	
	// check for LED state changes
	std::string tmp;
	if(m_data->CStore.Get("LED",tmp) && tmp=="Change"){
		
		if(m_data->CStore.Get("White",tmp)){
			get_ok = get_ok && SetState("White", (tmp=="1"));
		}
		
		if(m_data->CStore.Get("275",tmp)){
			get_ok = get_ok && SetState("LED275", (tmp=="1"));
		}
		
		if(m_data->CStore.Get("Deuterium",tmp)){
			get_ok = get_ok && SetState("D2", (tmp=="1"));
		}
		
		if(m_data->CStore.Get("Tungsten",tmp)){
			get_ok = get_ok && SetState("Tungsten", (tmp=="1"));
		}
		
	}
	
	if(m_data->CStore.Get("Shutter_gad",tmp)){
		get_ok = get_ok && SetState("GAD", (tmp=="CLOSE"));
	}
	
	if(m_data->CStore.Get("Shutter_ref",tmp)){
		get_ok = get_ok && SetState("REF", (tmp=="CLOSE"));
	}
	
	if(m_data->CStore.Get("Shutter_lamp",tmp)){
		get_ok = get_ok && SetState("LAMP", (tmp=="CLOSE"));
	}
	
	if(m_data->CStore.Get("Valve_inlet",tmp)){
		get_ok = get_ok && SetState("TUBE", (tmp=="CLOSE"));
	}
	
	if(m_data->CStore.Get("Valve_outlet",tmp)){
		get_ok = get_ok && SetState("PARALLEL", (tmp=="CLOSE"));
	}
	
	return get_ok;
}


bool ArduinoControl::Finalise(){
	
	bool allok = ShutItDown();
	allok = allok && Disconnect();
	
	return allok;
	
}

bool ArduinoControl::SetState(const std::string& name, bool enable){
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
	
	// turn off lights, valves, close shutters
	std::map<std::string, bool> results;
	results["LED275"] = Disable("LED275");
	results["White"] = Disable("White");
	results["GAD"] = Disable("GAD");
	results["REF"] = Disable("REF");
	results["Lamp"] = Disable("Lamp");
	
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
	// send command
	get_ok = SerialWrite(msg);
	if(!get_ok){
		Log(m_unique_name+" Error sending message '"+msg+"'",v_error,verbosity);
		return false;
	}
	// read response
	response = SerialRead();
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
			response += std::string((const char*)(buf),bytesread);
		}
		if(buf[bytesread-1]=='\n' || buf[bytesread-1]=='\0'){
			// received termination char
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
	return response;
}

