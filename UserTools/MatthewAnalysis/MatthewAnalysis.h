#ifndef MatthewAnalysis_H
#define MatthewAnalysis_H

#include <string>
#include <iostream>

#include "Tool.h"

#include "TGraphErrors.h"
#include "TTree.h"

class MatthewAnalysis: public Tool {

 public:

  MatthewAnalysis();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:
  
  void RetrieveDarkSubPure();
  void RetrieveCalibrationCurve();
  bool ReadyToAnalyse();
  std::pair<TTree*, TTree*> SortTrees(); 
  double* FunctionalFit(TGraphErrors&);
  TGraphErrors RemoveRegion(TGraphErrors&);
  TGraphErrors PerformDarkSubtraction(TTree*,  TTree*);
  TGraphErrors Scale(TGraphErrors&, const double*, const double*);
  TGraphErrors CalculateAbsorbance( TGraphErrors&, TGraphErrors&);
  std::vector<double> FindPeakDifference(TGraphErrors&);
  std::vector<double> CalculateConcentration(const std::vector<double>);
  
};


#endif
