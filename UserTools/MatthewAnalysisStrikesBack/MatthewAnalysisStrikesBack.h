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

struct LEDInfo {
  //  std::string name = "";
  std::string calib_fname = "";
  std::string calib_curve_name = "";
  int dark_offset = 0;
  //TGraph pure_ds;
  //TGraph high_conc_ds;
  //std::vector<double> initial_combined_fit_param_values;
  //std::vector<std::pair<double, double>> fit_param_ranges;#

  std::shared_ptr<CombinedGdPureFunc_DATA> combined_func;
  std::shared_ptr<AbsFunc> abs_func;
  FunctionalFit combined_fit;
  FunctionalFit absorbtion_fit;
  TF1* calibration_curve_ptr;  
};

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
  
  TGraph GetPure(const std::string& name);
  TGraph GetPureFromDB(const int&, const std::string&) const;
  TGraph GetHighConc(const std::string& name);
  TGraph GetHighConcFromDB(const int&, const std::string&) const;
  TF1* GetCalibrationCurve(const std::string&);
  void GetCurrentLED();
  void GetDarkAndLEDTrees();
  bool ReadyToAnalyse() const;
};


#endif
