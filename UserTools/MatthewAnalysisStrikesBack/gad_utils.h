#ifndef GAD_UTILS_H
#define GAD_UTILS_H

#include <iostream>
#include <memory>
#include <map>
#include <stdexcept>
#include <filesystem>
#include <algorithm>

#include "TApplication.h"
#include "TSystem.h"
#include <thread>
#include <chrono>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TCanvas.h"

[[maybe_unused]] static int number_of_points = 2068;
[[maybe_unused]] static int wave_min = 265;     // region suitable for fitting
[[maybe_unused]] static int wave_max = 290;     // (including absorption region)
[[maybe_unused]] static int wave_middle = 275;

[[maybe_unused]] static double abs_region_low = 270;
[[maybe_unused]] static double abs_region_high = 280;

[[maybe_unused]] static double sideband_region_low = 260;
[[maybe_unused]] static double sideband_region_high = 300;

[[maybe_unused]] static int abs_region_low_index = 0;
[[maybe_unused]] static int abs_region_high_index = 0;
[[maybe_unused]] static int abs_region_npoints = 0;

[[maybe_unused]] static int sideband_region_low_index = 0;
[[maybe_unused]] static int sideband_region_high_index = 0;
[[maybe_unused]] static int sideband_region_npoints = 0;


static std::string dark_name = "Dark";

class Func {
  // class with Evaluate method - replaces lambda based method from before
public:
  int n_fit = 0;
  double fit_min = 0;
  double fit_max = 0;
  virtual double Evaluate(double* x, double* p) = 0;
};

class SimplePureFunc : public Func {
public:
  TGraph pure_dark_subtracted;

  const int PURE_SCALING = n_fit++;
  const int PURE_FIRST_ORDER_CORRECTION = n_fit++;
  const int PURE_TRANSLATION = n_fit++;
  const int LINEAR_BACKGROUND_GRADIENT = n_fit++;
  const int LINEAR_BACKGROUND_OFFSET = n_fit++;

  SimplePureFunc(const TGraph&);
  double Evaluate(double*, double*);
};

class CombinedGdPureFunc : public Func {
public:
  TGraph pure_dark_subtracted;
  TGraph spec_abs;

  const int PURE_SCALING = n_fit++;
  const int PURE_FIRST_ORDER_CORRECTION = n_fit++;
  const int PURE_TRANSLATION = n_fit++;
  const int ABS_SCALING = n_fit++;
  const int ABS_TRANSLATION = n_fit++;
  const int FIRST_ORDER_BACKGROUND = n_fit++;
  const int ZEROTH_ORDER_BACKGROUND = n_fit++;
  
  CombinedGdPureFunc(const TGraph&, const TGraph&);
  double Evaluate(double*, double*);
};

class CombinedGdPureFunc_DATA : public Func {
public:
  TGraph pure_dark_subtracted;
  TGraph rat_abs;

  const int PURE_SCALING = n_fit++;
  //const int PURE_FIRST_ORDER_CORRECTION = n_fit++;
  //const int PURE_SECOND_ORDER_CORRECTION = n_fit++;
  const int PURE_TRANSLATION = n_fit++;
  const int PURE_STRETCH = n_fit++;
  const int ABS_SCALING = n_fit++;
  //const int ABS_TRANSLATION = n_fit++;
  //const int ABS_FIRST_ORDER_CORRECTION = n_fit++;
  //const int ABS_SECOND_ORDER_CORRECTION = n_fit++;
  //const int ABS_STRETCH = n_fit++;
  const int SECOND_ORDER_BACKGROUND = n_fit++;
  const int FIRST_ORDER_BACKGROUND = n_fit++;
  const int ZEROTH_ORDER_BACKGROUND = n_fit++;
  
  CombinedGdPureFunc_DATA(const TGraph&, const TGraph&);
  double Evaluate(double*, double*);
};


class AbsFunc : public Func {
public:
  TGraph abs_ds;

  const int ABS_SCALING = n_fit++;
  const int ABS_TRANSLATION = n_fit++;
  const int THIRD_BACKGROUND = n_fit++;
  const int SECOND_BACKGROUND = n_fit++;
  const int FIRST_BACKGROUND = n_fit++;
  const int ZEROTH_BACKGROUND = n_fit++;

  AbsFunc(const TGraph&);
  double Evaluate(double*, double*);
};

class FunctionalFit {
public:
  TF1 fit_funct;
  std::string fit_name = "";
  TGraph example;
  bool fitted = false;

  FunctionalFit(Func*, const std::string&);
  FunctionalFit() = default;
  void SetFitParameters(const std::vector<double>&);
  void SetFitParameterRanges(const std::vector<std::pair<double, double>>&);
  TFitResultPtr PerformFitOnData(TGraph, bool i = false);
  void SetExampleGraph(const TGraph&);
  double GetParameterValue(const int&) const;
  double GetChiSquared() const ;
  TGraph GetGraph() const;
  TGraph GetGraphOfJust(const std::vector<int>&) const;
  TGraph GetGraphExcluding(const std::vector<int>&) const;

};

class LEDInfo {
  public:
  LEDInfo();
  ~LEDInfo();
  LEDInfo(LEDInfo&& rhs); // move constructor
  LEDInfo(const LEDInfo& rhs) = delete;  // copy constructor
  LEDInfo& operator=(const LEDInfo&) = delete; // copy assignment
  // Matthew i don't know what you're up to with these shared pointers
  // so to be safe i'm gonna delete the copy constructor; feel free
  // to implement it if needed.
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

TGraph DarkSubtractFromTreePtrs(TTree*, TTree*, const int);
TGraph GetDarkSubtractFromFile(const std::string, const std::string, const int&);
TGraph RemoveRegion(const TGraph&, const double&, const double&);

template<class L>
TGraph BinaryOperation(const TGraph&, const TGraph&, const std::string&, const L&);

TGraph CalculateGdAbs(const TGraph&, const TGraph&);
TGraph PWRatio(const TGraph&, const TGraph&);
TGraph PWLogRatio(const TGraph&, const TGraph&);
TGraph PWPercentageDiff(const TGraph&, const TGraph&);
TGraph PWDifference(const TGraph&, const TGraph&);
TGraph PWMultiply(const TGraph&, const TGraph&);
double SuccessMetric(const TGraph&);

double GetPeakDiff(const TGraph&);
double GetPeakRatio(const TGraph&);
TGraph TrimGraph(const TGraph&, double l=wave_min, double h=wave_max);
TGraph Normalise(TGraph);
TGraph ZeroNegative(TGraph);
#endif 
