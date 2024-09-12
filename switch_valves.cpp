#include "turn_off_arduino.h"
#include "rs232.h"
#include <iostream>
#include <map>
#include <vector>
#include <thread>
#include <chrono>
#include <algorithm>

int get_ok=0;
int com_port=0;
int baud_rate=0;

bool Connect();
bool Disconnect();
bool Disable(const std::string& name);
bool Enable(const std::string& name);
bool SetState(const std::string& name, bool enable);
bool SendAndReceive(std::string msg, std::string& response, int timeout=3000);
bool SerialWrite(std::string msg);
std::string SerialRead(int timeout_ms=1000);
int DeviceNameToNumber(std::string name);

struct timer {
	public:
	void reset(){
		start = std::chrono::steady_clock::now();
	}
	unsigned long long ms_elapsed() const {
		return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
	}
	private:
	std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
} atimer;

int main(int argc, const char* argv[]){
	
	bool open_gad=false;
	if(argc<3){
		std::cout<<"Usage: "<<argv[0]<<" <com_port> (tube=1|parallel=0, default)"<<std::endl;
	} else {
		try {
			open_gad = std::stoi(std::string(argv[2]));
		} catch(std::exception& e){
			std::cerr<<"error parsing arg '"<<argv[2]
			         <<"' - use 1 to open GAD, 0 to open parallel"<<std::endl;
			return 1;
			}
	}
	
	// serial settings
	baud_rate=19200;
	std::string com_port_str="/dev/ttyACM0";
	if(argc>1){
		com_port_str = argv[1];
	}
	com_port = DeviceNameToNumber(com_port_str);
	if(com_port<0){
		std::cerr<<"invalid com port '"<<com_port_str<<"'"<<std::endl;
		return 2;
	}
	std::cout<<" com port '"<<com_port_str<<"' mapped to port #"<<std::to_string(com_port)<<std::endl;
	bool get_ok = Connect();
	if(!get_ok){
		std::cerr<<" failed to connect to arduino!"<<std::endl;
		return 3;
	}
	
	// close valves
	// TODO implement a function to query the current state of the valves
	// for now, take user argument
	std::map<std::string, bool> results;
	if(open_gad){
		results["TUBE"] = Enable("TUBE");
		results["PARALLEL"] = Disable("PARALLEL");
	} else {
		results["PARALLEL"] = Enable("PARALLEL");
		results["TUBE"] = Disable("TUBE");
	}
	
	bool allok=true;
	for(auto&& res : results){
		if(!res.second){
			allok = false;
			std::cerr<<"error switching "+res.first<<std::endl;
		}
	}
	
	return (allok==0) ? 0 : 4;
}

bool Connect(){
	
	// default arduino serial is 8 data bits, 1 stop bit, no parity, no flow control
	char mode[]={'8','N','1',0};
	int flowctl=0;
	
	if(RS232_OpenComport(com_port, baud_rate, mode, flowctl)){
		std::cerr<<" Unable to open com port "<<std::to_string(com_port)<<std::endl;
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
		if(get_ok) std::cerr<<"Error with initial comms check; expected response 'Hello!', got '"
		                    <<response<<"'"<<std::endl;
		return false;
	}
	
	// request arduino to connect to the grove relay controller over SPI
	get_ok = SendAndReceive("GROVE",response);
	if(!get_ok) return false;
	
	return get_ok;
	
}

bool Disconnect(){
	// tell arduino to close serial comms
	// (it'll open it again in 5 seconds and start listening for new connections)
	get_ok = SerialWrite("QUIT");
	if(!get_ok){
		std::cerr<<"Error sending message 'QUIT'"<<std::endl;
		return false;
	}
	return true;
}

bool SetState(const std::string& name, bool enable){
	std::cout<<"SetState("<<name<<", "<<std::to_string(enable)<<")"<<std::endl;
	if(enable) return Enable(name);
	return Disable(name);
}

bool Disable(const std::string& name){
	std::string resp;
	return SendAndReceive(name+" DISABLE", resp);
}

bool Enable(const std::string& name){
	std::string resp;
	return SendAndReceive(name+" ENABLE", resp);
}

bool SendAndReceive(std::string msg, std::string& response, int timeout){
	std::cout<<"SendAndReceive of command '"<<msg<<"'"<<std::endl;
	// send command
	get_ok = SerialWrite(msg);
	if(!get_ok){
		std::cerr<<" Error sending message '"<<msg<<"'"<<std::endl;
		return false;
	}
	// read response
	//std::cout<<"Reading response..."<<std::endl;
	response = SerialRead(timeout);
	if(response.empty()){
		std::cerr<<"Error: No response to command '"<<msg<<"'"<<std::endl;
		return false;
	}
	// as a standard, the arduino sketch returns strings starting with 'Err: ABC' when it fails to do something
	if(response.substr(0,4)=="Err:"){
		std::cerr<<"Error: command '"<<msg<<"' returned '"<<response<<"'"<<std::endl;
		return false;
	}
	
	return true;
}

bool SerialWrite(std::string msg){
	if(msg.size()==0){
		std::cerr<<"SerialWrite with empty message!"<<std::endl;
		return false;
	}
	if(msg.back()!='\n') msg.push_back('\n');
	get_ok = RS232_SendBuf(com_port, (unsigned char*)(msg.data()), msg.size());
	if(get_ok!=msg.size()){
		std::cerr<<"Error writing message '"<<msg<<"'; sent "<<get_ok
		         <<" bytes of "<<msg.size()<<std::endl;
		return false;
	}
	return true;
}

std::string SerialRead(int timeout_ms){
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
			std::cerr<<"Warning! Too many chars in serial buffer! Flushing any remainder"<<std::endl;
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

int DeviceNameToNumber(std::string name){
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
