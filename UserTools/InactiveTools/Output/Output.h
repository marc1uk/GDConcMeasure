#ifndef Output_H
#define Output_H

#include <string>
#include <iostream>

#include "Tool.h"

class Output: public Tool {


 public:

  Output();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:





};


#endif
