#ifndef BenWriter_H
#define BenWriter_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "TTree.h"
#include "TFile.h"
#include "TObject.h"

class BenWriter: public Tool {


 public:

  BenWriter();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:

  TFile* file;
  short m_month;

  void InitTTree(TTree* tree);

};


#endif
