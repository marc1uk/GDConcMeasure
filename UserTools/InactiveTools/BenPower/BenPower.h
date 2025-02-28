#ifndef BenPower_H
#define BenPower_H

#include <string>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include "Tool.h"

class BenPower: public Tool {
  
  
public:
  
  BenPower();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();
  
  bool TurnOn();
  bool TurnOff();
  
private:
  
  int m_powerpin;
  std::string power;
  
  int verbosity=1;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  
};


#endif
