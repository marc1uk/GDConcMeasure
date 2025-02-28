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
  std::string valve;
  std::string CStoreKey;
  
  bool ValveOpen();
  bool ValveClose();

  int verbosity=1;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;

};


#endif
