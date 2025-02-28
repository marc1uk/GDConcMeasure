#ifndef StateTest_H
#define StateTest_H

#include <string>
#include <iostream>

#include "Tool.h"

class StateTest: public Tool {


 public:

  StateTest();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:





};


#endif
