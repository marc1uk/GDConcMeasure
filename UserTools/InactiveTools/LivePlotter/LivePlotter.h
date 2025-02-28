#ifndef LivePlotter_H
#define LivePlotter_H

#include <string>
#include <iostream>

#include "Tool.h"

class LivePlotter: public Tool {
  
  public:
  
  LivePlotter();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();
  
  private:
  int verbosity=1;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  
};


#endif
