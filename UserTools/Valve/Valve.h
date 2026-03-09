#ifndef Valve_H
#define Valve_H

#include <string>
#include <iostream>

#include "Tool.h"

class Valve: public Tool {


 public:

  Valve();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:

  int m_valve_pin;
  int m_switching_pin;
  int m_holding_pin;
  std::string valve;
  std::string CStoreKey;
  std::string type;
  
  bool ValveOpen();
  bool ValveClose();
  bool ConfigurePin(int pin_num);
  bool SwitchPin(int pin_num, int state);
  
  int verbosity=1;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;

};


#endif
