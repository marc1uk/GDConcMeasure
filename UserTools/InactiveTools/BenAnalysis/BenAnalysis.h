#ifndef BenAnalysis_H
#define BenAnalysis_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

class BenAnalysis: public Tool {


 public:

  BenAnalysis();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:





};


#endif
