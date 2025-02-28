#ifndef Arduino_H
#define Arduino_H

#include "Logging.h"

#include <string>
#include <iostream>

using ToolFramework::Logging;

class Arduino {
	
	public:
	bool Configure(std::string com_port_str, int baud, int verb, Logging* logger=nullptr);
	int DeviceNameToNumber(std::string str);
	bool Connect();
	bool Disconnect();
	
	std::string SerialRead(int timeout_ms=1000);
	bool SerialWrite(std::string msg);
	bool SendAndReceive(std::string msg, std::string& response, int timeout=3000);
	
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
	
	private:
	//timer atimer;
	int com_port;
	int baud_rate;
	bool connected=false;
	
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
	
	std::string m_unique_name="Arduino";
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	int get_ok;
	
};

#endif
