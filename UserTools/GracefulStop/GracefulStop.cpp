/* vim:set noexpandtab tabstop=4 wrap */
#include "GracefulStop.h"

GracefulStop::GracefulStop():Tool(){}

bool GracefulStop::gotStopSignal = false;

bool GracefulStop::Initialise(std::string configfile, DataModel &data){
	
	/////////////////// Useful header ///////////////////////
	if(configfile!="") m_variables.Initialise(configfile); //loading config file
	//m_variables.Print();
	
	m_data= &data; //assigning transient data pointer
	/////////////////////////////////////////////////////////////////
	
	// TODO potentially we could take from config file what type of signal we wish to look for.
	
	if(signal((int) SIGUSR1, GracefulStop::stopSignalHandler) == SIG_ERR){
		Log("GracefulStop Tool: Failed to setup signal handler!", v_error, verbosity);
		return false;
	}
	
	return true;
}

/**
* This function is registered as the signal handler. It gets called when the program receives the signal, and sets 'gotStopSignal' member to true.
* On Execute(), if this member is true, StopLoop will be set.
* @param[in] _ignored Not used here (normally the signal type) but required to match the signal handler function prototype.
*/
void GracefulStop::stopSignalHandler(int _ignored){
	// technically we could choose what to do based on the signal type passed, if we registered this function with multliple signals.
	gotStopSignal = true;
}

bool GracefulStop::Execute(){
	if(gotStopSignal){
		// check if we're in the middle of taking a measurement,
		// and allow it to finish if so
		std::string command_step="0";
		//std::cout<<"checking if we have command_step"<<std::endl;
		//get_ok = m_data->CStore.Has("command_step");
		//std::cout<<"returned "<<get_ok<<std::endl;
		//m_data->CStore.Print(true);
		//std::cout<<"getting command step"<<std::endl;
		get_ok = m_data->CStore.Get("command_step",command_step);
		//std::cout<<"graceful stop get_ok="<<get_ok<<std::endl;
		Log("GracefulStop Tool: Received SIGUSR1. Command step is "+command_step,v_error,verbosity);
		if(get_ok && command_step!="0"){
				Log(", allowing measurement to complete...",v_error,verbosity);
		} else {
			Log(", terminating ToolChain",v_error,verbosity);
			m_data->vars.Set("StopLoop",1);
			gotStopSignal=false;
		}
	}
	return true;
}

bool GracefulStop::Finalise(){
	return true;
}
