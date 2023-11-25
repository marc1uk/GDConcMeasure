#ifndef MatthewAnalysisStrikesBack_H
#define MatthewAnalysisStrikesBack_H

#include "TGraph.h"
#include "TTree.h"

#include <map>
#include <string>
#include <iostream>
#include <memory>

#include "Tool.h"
#include "gad_utils.h"

class MatthewAnalysisStrikesBack: public Tool {

 public:

  MatthewAnalysisStrikesBack();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

  private:

  std::map<std::string, LEDInfo> led_info_map;
  std::string current_led = "";
  TTree* led_tree_ptr = nullptr;
  TTree* dark_tree_ptr = nullptr;
  
  
  TGraph dark_subtracted_data_in;
  TGraph dark_subtracted_data_out;
  TGraph purefitgraph;
  TGraph absfitgraph;
  // yeah this is a bit messy
  TFitResult datafitres;
  TFitResult absfitres;
  TFitResultPtr datafitresp;
  TFitResultPtr absfitresp;
  
  // although it isn't *required*, I do still think it's helpful
  // to name your function parameters in the declaration....
  TGraph GetPure(const std::string& name);
  TGraph GetPureFromDB(const int&, const std::string&) const;
  TGraph GetHighConc(const std::string& name);
  TGraph GetHighConcFromDB(const int&, const std::string&) const;
  TF1* GetCalibrationCurve(const std::string&);
  TF1* GetCalibCurveFromDB(const int& calibID, const std::string& ledToAnalyse);
  void GetCurrentLED();
  void GetDarkAndLEDTrees();
  bool ReadyToAnalyse() const;
  bool UpdateDataModel(double metric, double conc);
  bool ReinitDataModel();
  void GetWavelengthIndices(const TGraph& g);
  void GetAbsRegion(const TGraph& datagraph, TGraph& absgraph);
  void GetSidebandRegion(const TGraph& datagraph, TGraph& absgraph);
  void GetRawMinMax(double& raw_min, double& raw_max);
  void GetDarkTraceParams(double& mean, double& sigma);
  
  // debugging
  TFile* fdebug = nullptr;
  int SaveDebug(const TObject* obj, const std::string& name);
  int SaveDebug(const TObject& obj, const std::string& name);
  
  int get_ok;
  int m_verbose=1;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
};


#endif
