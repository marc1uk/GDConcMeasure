#ifndef ARDUINO
#include "Arduino.h"
#include "rs232.h"
#include <thread>
#include <chrono>
#include <vector>
#include <algorithm>

bool Arduino::Configure(std::string com_port_str, int baud, int verb, Logging* logger){
	
	m_log = logger;
	verbosity = verb;
	
	com_port = DeviceNameToNumber(com_port_str);
	if(com_port<0){
		Log(m_unique_name+" invalid com port '"+com_port_str+"'",v_error,verbosity);
		return false;
	}
	Log(m_unique_name+" com port '"+com_port_str+"' mapped to port #"+std::to_string(com_port),v_debug,verbosity);
	
	baud_rate = baud;
	return true;
}


bool Arduino::Connect(){
	
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
	
	connected = true; // must set before SendAndReceive can be called
	
	// check comms are working as expected
	std::string response;
	get_ok = SendAndReceive("HELLO",response);
	if(!get_ok){
		Log(m_unique_name+" Error with initial comms check!",v_error,verbosity);
			connected=false;
			return false;
	}
	if(response!="Hello!"){
		Log(m_unique_name+" Warning with initial comms check; expected response 'Hello!', got '"+response+"'",v_warning,verbosity);
		// sometimes happens after re-flashing or initial power up... TODO match and discard.
		// (actually this should already be done by the above RS232_FlushRXTX but ???)
		//connected=false;
		//return false;
	}
	
	// request arduino to connect to the grove relay controller over SPI
	get_ok = SendAndReceive("GROVE",response);
	if(!get_ok){
		connected=false;
		return false;
	}
	
	return get_ok;
	
}

bool Arduino::Disconnect(){
	
	if(!connected) return true;
	
	// tell arduino to close serial comms
	// (it'll open it again in 5 seconds and start listening for new connections)
	get_ok = SerialWrite("QUIT");
	if(!get_ok){
		Log(m_unique_name+" Error sending message 'QUIT'",v_error,verbosity);
		return false;
	}
	
	return true;
	
}

bool Arduino::SendAndReceive(std::string msg, std::string& response, int timeout){
	
	if(!connected){
		Log(m_unique_name+" SendAndReceive before connect!",v_error,verbosity);
		return false;
	}
	
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

bool Arduino::SerialWrite(std::string msg){
	
	if(!connected){
		Log(m_unique_name+" SerialWrite before connect!",v_error,verbosity);
		return false;
	}
	
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

std::string Arduino::SerialRead(int timeout_ms){
	
	if(!connected){
		Log(m_unique_name+" SerialRead before connect!",v_error,verbosity);
		return "";
	}
	
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

int Arduino::DeviceNameToNumber(std::string name){
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
#endif
