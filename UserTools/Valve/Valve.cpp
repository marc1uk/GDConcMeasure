#include "Valve.h"
#include "Algorithms.h"
#include <thread>
#include <chrono>

Valve::Valve():Tool(){}


bool Valve::Initialise(std::string configfile, DataModel &data){
  
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
  if(configfile!="")  m_variables.Initialise(configfile);
  
  //m_variables.Print();
  
  
  m_variables.Get("verbosity",verbosity);
  
  // get which valve we're controlling
  std::string type="";
  m_variables.Get("type",type);   // 'inlet', 'outlet' or 'pump'
  if(type!="inlet" && type!="outlet" && type!= "pump"){
    Log("Valve unrecognised type '"+type+"'",v_error,verbosity);
    return false;
  }
  // we'll look for a corresponding flag in the DataModel
  CStoreKey = "Valve_"+type;
  
  // get the pin connected to this valve
  m_valve_pin=-1;
  if (!m_variables.Get("valve_pin",m_valve_pin)){
    Log("Valve pin not set",v_error,verbosity);
    return false;
  }
  
  std::stringstream command;
  command<<"if [ ! -d /sys/class/gpio/gpio"<<m_valve_pin<<" ]; then echo \""
         <<m_valve_pin<<"\" > /sys/class/gpio/export; fi";
  std::string errmsg;
  int ok = SystemCall(command.str(), errmsg);
  if(ok!=0){
    Log("Valve::Initialise "+errmsg,0,0);
    return false;
  }
  
  command.str("");
  command<<"STATE=$(cat /sys/class/gpio/gpio"<<m_valve_pin<<"/direction); if [ \"${STATE}\" != \"out\" ]; "
         <<"then echo \"out\" > /sys/class/gpio/gpio"<<m_valve_pin<<"/direction; fi";
  ok = SystemCall(command.str(),errmsg);
  if(ok!=0){
    Log("Valve::Initialise "+errmsg,0,0);
    return false;
  }
  
  ok = ValveClose();
  if(not ok) return ok;
  
  return true;
}


bool Valve::Execute(){
  
  Log(CStoreKey+" Executing...",v_debug,verbosity);
  
  std::string Valve="";
  bool ok = true;
  
  if(m_data->CStore.Get(CStoreKey,Valve) && Valve!=valve){
    
    // set the CStore status to what the actual valve state is
    m_data->CStore.Set(CStoreKey,valve);
    
    if(Valve=="OPEN"){
      Log(CStoreKey+"::Execute got OPEN",v_debug,verbosity);
      ok = ValveOpen();
    }
    if(Valve=="CLOSE"){
      Log(CStoreKey+"::Execute got CLOSE",v_debug,verbosity);
      ok = ValveClose();
    }
    if(ok) m_data->CStore.Set(CStoreKey,valve);
    
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  }
  
  return ok;
}


bool Valve::Finalise(){
  
  bool ok = ValveClose();
  
  return ok;
}


bool Valve::ValveOpen(){
  
  Log("valve open",v_message,verbosity);
  if(m_valve_pin<0){
    Log(CStoreKey+"::ValveOpen invalid valve pin "+std::to_string(m_valve_pin),v_error,verbosity);
    return false;
  }
  
  std::stringstream command;
  command<<"echo \"1\" > /sys/class/gpio/gpio"<<m_valve_pin<<"/value";
  std::string errmsg;
  int ok = SystemCall(command.str(), errmsg);
  if(ok!=0){
    Log(CStoreKey+"::ValveOpen "+errmsg,0,0);
    return false;
  }
  valve="OPEN";
  return true;
}

bool Valve::ValveClose(){
  
  Log("valve close",v_message,verbosity);
  if(m_valve_pin<0){
    Log(CStoreKey+"::ValveClose invalid valve pin "+std::to_string(m_valve_pin),v_error,verbosity);
    return false;
  }
  
  std::stringstream command;
  command<<"echo \"0\" > /sys/class/gpio/gpio"<<m_valve_pin<<"/value";
  std::string errmsg;
  int ok = SystemCall(command.str(), errmsg);
  if(ok!=0){
    Log(CStoreKey+"::ValveClose "+errmsg,0,0);
    return false;
  }
  valve="CLOSE";
  return true;
}
