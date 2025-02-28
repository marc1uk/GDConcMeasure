#ifndef Increment_H
#define Increment_H

#include <string>
#include <iostream>

#include "Tool.h"

class Increment: public Tool {


 public:

  Increment();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:





};


#endif
