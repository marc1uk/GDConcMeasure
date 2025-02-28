#include "Monitoring.h"

#include <array>

Monitoring::Monitoring():Tool(){}


bool Monitoring::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	// TODO add finer granularity over rates of specific things
	m_variables.Get("mon_period_s",send_period_s);
	
	last_send = std::chrono::high_resolution_clock::now();
	
	m_data= &data;
	
	if(!m_data->arduino){
		Log(m_unique_name+"::Initialise - Warning - no arduino controller found in datamodel!",v_error,verbosity);
	}
	
	m_data->vars.Set("Status","Initialising");
	
	return true;
}


bool Monitoring::Execute(){
	
	m_data->vars.Set("Status","Running");
	
	auto secs_since_last_send = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - last_send);
	if(secs_since_last_send.count() > send_period_s){
		
		// TODO add non-arduino monitoring here.
		// e.g. latest concentration, transparency trace, LED intensities...
		// all the things we monitore in the DB
		
		
		if(m_data->arduino){
			
			//double led_temp
			//m_data->arduino->GetLedTemp(led_temp);  // currently not on WCTE GAD
			//m_data->monitoring_store.Set("led_temp",led_temp);
			
			double flow_rate;
			m_data->arduino->GetFlowRate(flow_rate);
			m_data->monitoring_store.Set("flow_rate",flow_rate);
			
			int leak_sense;
			m_data->arduino->GetLeakStatus(leak_sense);
			m_data->monitoring_store.Set("leak_check",leak_sense);
			
			std::array<double,3> sol_temps;
			m_data->arduino->GetSolTemps(sol_temps);
			for(int i=0; i<3; ++i){
				std::string key="sol_temp_"+std::to_string(i);
				m_data->monitoring_store.Set(key,sol_temps.at(i));
			}
			
		}
		
		std::string json="";
		m_data->monitoring_store>>json;
		m_data->services->SendMonitoringData(json, "GAD info"); // second arg seems to be optional?
		
	}
	
	return true;
}


bool Monitoring::Finalise(){
	
	m_data->vars.Set("Status","Stopped");
	return true;
}
