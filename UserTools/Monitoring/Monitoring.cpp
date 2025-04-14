#include "Monitoring.h"

#include <array>

Monitoring::Monitoring():Tool(){}


bool Monitoring::Initialise(std::string configfile, DataModel &data){
	
	InitialiseTool(data);
	m_configfile=configfile;
	InitialiseConfiguration(configfile);
	
	// Retrieve configuration options from the postgres database
	int RunConfig=-1;
	m_data->vars.Get("RunConfig",RunConfig);
	Log(m_unique_name+" getting RunConfig from ToolChain",v_debug,m_verbose);
	if(RunConfig>=0){
		std::string configtext;
		Log(m_unique_name+" getting tool configuration from database for runconfig "
		    +std::to_string(RunConfig),v_debug,m_verbose);
		bool get_ok = m_data->postgres_helper.GetToolConfig(m_unique_name, configtext);
		if(!get_ok){
			Log(m_unique_name+" Failed to get Tool config from database!",v_error,m_verbose);
			return false;
		}
		// parse the configuration to populate the m_variables Store.
		if(configtext!="") m_variables.Initialise(std::stringstream(configtext));
	}
	
	//overide with any variables in local file
	Log(m_unique_name+" reading local overrides from config file "+configfile,v_debug,m_verbose);
	if(configfile!="") m_variables.Initialise(configfile);
	
	//m_variables.Print();
	m_variables.Get("verbosity",m_verbose);
	
	
	if(!m_variables.Get("mon_period_s",send_period_s)) send_period_s=60;
	last_send = std::chrono::high_resolution_clock::now();
	
	if(!m_data->arduino){
		Log(m_unique_name+"::Initialise - Warning - no arduino controller found in datamodel!",v_error,m_verbose);
	}
	
	m_data->vars.Set("Status","Initialising");
	
	return true;
}


bool Monitoring::Execute(){
	
	m_data->vars.Set("Status","Running");
	
	auto secs_since_last_send = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - last_send);
	
	// can't do it just based on time because this may stall the toolchain while the UV LED is on - we need to use
	// the scheduler to synchronize. However we can limit the rate by skipping if not enough time has passed.
	std::string monitoring="";
	if(m_data->CStore.Get("Monitor", monitoring) /*&& (secs_since_last_send.count() > send_period_s)*/){
		
		Log(m_unique_name+" Sending monitoring Info!",v_debug,m_verbose);
		// TODO add finer granularity over rates of specific things
		// by value in 'monitoring'...
		
		// TODO any other non-arduino monitoring here?
		
		if(m_data->arduino){
			Log(m_unique_name+" Got arduino!",v_debug,m_verbose);
			
			//double led_temp
			//m_data->arduino->GetLedTemp(led_temp);  // currently not on WCTE GAD
			//m_data->monitoring_store.Set("led_temp",led_temp);
			
			double flow_rate;
			m_data->arduino->GetFlowRate(flow_rate);
			m_data->monitoring_store.Set("flow_rate",flow_rate);
			
			int leak_sense;
			m_data->arduino->GetLeakStatus(leak_sense);
			m_data->monitoring_store.Set("leak_check",leak_sense);
			if(leak_sense<700){
				m_data->services->SendAlarm("GAD LEAK SENSORS INDICATE WATER!",0,"GAD");
			}
			
			std::array<double,3> sol_temps;
			m_data->arduino->GetSolTemps(sol_temps);
			double max_T=0;
			for(int i=0; i<3; ++i){
				std::string key="sol_temp_"+std::to_string(i);
				m_data->monitoring_store.Set(key,sol_temps.at(i));
				if(sol_temps.at(i)>max_T) max_T=sol_temps.at(i);
			}
			if(max_T>70){
				std::string msg="GAD TEMP SENSORS INDICATE SOLENOID OVER-TEMP! ";
				for(int i=0; i<3; ++i){
					if(i!=0) msg+=", ";
					msg+="Solenoid "+std::to_string(i)+" = "+sol_temps.at(i)+"°C";
				}
				m_data->services->SendAlarm(msg,0,"GAD");
			}
			
		}
		
		std::string json="";
		m_data->monitoring_store>>json;
		get_ok = m_data->services->SendMonitoringData(json);
		
		if(!get_ok){
			Log(m_unique_name+" ERROR: SendMonitoringData failed!",v_error,m_verbose);
			return false;
		} else {
			Log(m_unique_name+" Sent monitoring data '"+json+"'",v_debug,m_verbose);
		}
		
		last_send = std::chrono::high_resolution_clock::now();
		
	}
	
	m_data->CStore.Remove("Monitor");
	
	return true;
}


bool Monitoring::Finalise(){
	
	m_data->vars.Set("Status","Stopped");
	return true;
}
