#include "BenPower.h"
#include "Algorithms.h"

BenPower::BenPower() : Tool(){}


bool BenPower::Initialise(std::string configfile, DataModel &data){
  
  m_data= &data; //assigning transient data pointer
  
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
  
  /*- old method, read configuration from local file - */
  if(configfile!="") m_variables.Initialise(configfile); //loading config file
  
  //m_variables.Print();
  
  
  m_variables.Get("verbosity",verbosity);
  
  if(!m_variables.Get("wake_pin",m_powerpin)){
    Log("no wake pin for powersupply set",v_error,verbosity);
    return false;
  }
  
  Log("power pin is "+std::to_string(m_powerpin),v_message,verbosity);
  
  std::stringstream command;
  command<<"if [ ! -d /sys/class/gpio/gpio"<<m_powerpin<<" ]; then echo \""
         <<m_powerpin<<"\" > /sys/class/gpio/export; fi";
  std::string errmsg;
  int ok = SystemCall(command.str(), errmsg);
  if(ok!=0){
    Log("BenPower::Initialise "+errmsg,0,0);
    return false;
  }
  
  command.str("");
  command<<"echo \"out\" > /sys/class/gpio/gpio"<<m_powerpin<<"/direction";
  ok = SystemCall(command.str(), errmsg);
  if(ok!=0){
    Log("BenPower::Initialise "+errmsg,0,0);
    return false;
  }
  
  ok = TurnOff();
  
  return ok;
  
}

bool BenPower::Execute(){
  
  Log("BenPower Executing...",v_debug,verbosity);
  
  std::string Power="";
  int ok = true;
  
  if(m_data->CStore.Get("Power",Power) && Power!=power){
    
    Log("BenPower setting power to "+Power,v_debug,verbosity);
    m_data->CStore.Set("Power",power);
    
    if(Power=="ON") ok = TurnOn();
    else if (Power=="OFF") ok = TurnOff();
    
    if(ok){
        m_data->CStore.Set("Power",power);
        Log("Power is now "+power,v_debug,verbosity);
    }
    
  }
  
  return ok;
}


bool BenPower::Finalise(){
  
  bool ok = TurnOff();
  
  return ok;
}

bool BenPower::TurnOn(){
  
  //write to GPIO
  Log("power on",v_message,verbosity);
  std::stringstream command;
  command<<"echo \"0\" > /sys/class/gpio/gpio"<<m_powerpin<<"/value";
  
  std::string errmsg;
  int retval=SystemCall(command.str(), errmsg);
  if(retval!=0){
    Log("BenPower::TurnOn "+errmsg,0,0);
    return false;
  }
  
  sleep(4);
  power ="ON";
  
  return true;
  
}

bool BenPower::TurnOff(){
  
  //write to GPIO
  Log("power off",v_message,verbosity);
  std::stringstream command;
  command<<"echo \"1\" > /sys/class/gpio/gpio"<<m_powerpin<<"/value";
  
  std::string errmsg;
  int retval = SystemCall(command.str(),errmsg);
  if(retval!=0){
    Log("BenPower::TurnOff "+errmsg,0,0);
    return false;
  }
  
  sleep(4);
  power="OFF";
  
  return true;
  
}
