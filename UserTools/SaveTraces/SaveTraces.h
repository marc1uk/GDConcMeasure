#ifndef SaveTraces_H
#define SaveTraces_H

#include <string>
#include <iostream>

#include "Tool.h"

class SaveTraces: public Tool {


 public:

  SaveTraces();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:
 
  TFile* file=nullptr;
  std::string lastname="";
  int verbosity=1;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;




};


#endif
