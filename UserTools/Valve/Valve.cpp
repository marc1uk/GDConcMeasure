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
  type="";
  m_variables.Get("type",type);   // 'inlet', 'outlet' or 'pump'
  if(type!="tube" && type!="inlet" && type!="outlet" && type!= "pump"){
    Log("Valve unrecognised type '"+type+"'",v_error,verbosity);
    return false;
  }
  // we'll look for a corresponding flag in the DataModel
  CStoreKey = "Valve_"+type;
  
  int ok;
  
  // old behaviour: 1 pin per valve
  if(type!="tube"){
    
    // get the pin connected to this valve
    m_valve_pin=-1;
    if(!m_variables.Get("valve_pin",m_valve_pin)){
      Log("Valve pin not set",v_error,verbosity);
      return false;
    }
    if(!ConfigurePin(m_valve_pin)) return false;
    
  } else {
    
    // tube type has two pins, but they both control inlet and outlet
    // one pin applies a switching voltage, the other applies only a holding voltage
    m_switching_valve_pin=-1;
    m_holding_valve_pin=-1;
    if( (!m_variables.Get("switching_pin",m_switching_valve_pin))
     || (!m_variables.Get("holding_pin",m_holding_valve_pin)) ){
      Log("Switching or holding pin not set",v_error,verbosity);
      return false;
    }
    if(!ConfigurePin(m_switching_pin)) return false;
    if(!ConfigurePin(m_holding_pin)) return false;
  
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
  
  int ok;
  if(type!="tube"){
    ok = SwitchPin(m_valve_pin, 1);
  } else {
    ok = SwitchPin(m_switching_pin, 1) &&
         SwitchPin(m_holding_pin, 1);
    std::this_thread::sleep_for(std::chrono::seconds(1));
    ok = ok && SwitchPin(m_switching_pin, 0);
  }
  
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
  
  int ok;
  if(type!="tube"){
    ok = SwitchPin(m_valve_pin, 0);
  } else {
    ok = SwitchPin(m_switching_pin, 0) &&
         SwitchPin(m_holding_pin, 0);
  }
  
  if(ok!=0){
    Log(CStoreKey+"::ValveClose "+errmsg,0,0);
    return false;
  }
  valve="CLOSE";
  return true;
}

bool Valve::ConfigurePin(int pin_num){
  std::stringstream command;
  std::string errmsg;
  // configure for gpio control
  command<<"if [ ! -d /sys/class/gpio/gpio"<<pin_num<<" ]; then echo \""
         <<m_valve_pin<<"\" > /sys/class/gpio/export; fi";
  ok = SystemCall(command.str(), errmsg);
  if(ok!=0){
    Log("Valve::ConfigurePin "+std::to_string(pin_num)+": "+errmsg,0,0);
    return false;
  }
  // configure for output
  command.str("");
  command<<"STATE=$(cat /sys/class/gpio/gpio"<<pin_num<<"/direction); if [ \"${STATE}\" != \"out\" ]; "
         <<"then echo \"out\" > /sys/class/gpio/gpio"<<pin_num<<"/direction; fi";
  ok = SystemCall(command.str(),errmsg);
  if(ok!=0){
    Log("Valve::ConfigurePin "+std::to_string(pin_num)+": "+errmsg,0,0);
    return false;
  }
  return true;
}

bool Valve::SwitchPin(int pin_num, int state){
  std::stringstream command;
  command<<"echo \""<<state<<"\" > /sys/class/gpio/gpio"<<pin_num<<"/value";
  std::string errmsg;
  return SystemCall(command.str(), errmsg);
}
  
