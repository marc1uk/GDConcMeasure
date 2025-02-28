#ifndef MatthewAnalysis_H
#define MatthewAnalysis_H

#include <string>
#include <iostream>
#include <utility>

#include "Tool.h"

#include "TGraphErrors.h"
#include "TFitResultPtr.h"
class TTree;
class TF1;

class MatthewAnalysis: public Tool {

 public:

  MatthewAnalysis();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:
  void MakeLivePlot();
  bool RetrieveDarkSubPure();
  bool RetrieveCalibrationCurve();
  bool RetrieveCalibrationCurveDB();
  bool ReadyToAnalyse();
  std::pair<TTree*, TTree*> FindDarkAndLEDTrees(); 
  double* FunctionalFit(TGraphErrors&);
  TGraphErrors RemoveRegion(TGraphErrors&);
  TGraphErrors PerformDarkSubtraction(TTree*, TTree*);
  TGraphErrors Scale(TGraphErrors&, const double*, const double*);
  TGraphErrors CalculateAbsorbance(TGraphErrors&, TGraphErrors&);
  std::pair<double,double> FindPeakDifference(TGraphErrors&);
  std::pair<double,double> CalculateConcentration(const std::pair<double,double>&);
  void ReInit();
  
  std::string led_name;
  
  TGraphErrors dark_subtracted_gd;
  TGraphErrors pure_scaled;
  TGraphErrors absorbance;
  TFitResultPtr fitresptr;
  
  TFitResultPtr lpf, rpf;
  TF1* left_peak = nullptr;
  TF1* right_peak = nullptr;
  
  int verbosity=1;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
  
};


#endif
