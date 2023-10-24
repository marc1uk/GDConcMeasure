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
#include "TCanvas.h"

[[maybe_unused]] static int number_of_points = 2068;
[[maybe_unused]] static int wave_min = 260;
[[maybe_unused]] static int wave_max = 290;

[[maybe_unused]] static double abs_region_low = 270;
[[maybe_unused]] static double abs_region_high = 280;

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
  TGraph pure_dark_subtracted;  
public:

  const int PURE_SCALING = n_fit++;
  const int PURE_FIRST_ORDER_CORRECTION = n_fit++;
  const int PURE_TRANSLATION = n_fit++;
  const int LINEAR_BACKGROUND_GRADIENT = n_fit++;
  const int LINEAR_BACKGROUND_OFFSET = n_fit++;

  SimplePureFunc(const TGraph&);
  double Evaluate(double*, double*);
};

class CombinedGdPureFunc : public Func {
  TGraph pure_dark_subtracted;
  TGraph spec_abs;
public:

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
  TGraph pure_dark_subtracted;
  TGraph rat_abs;
public:

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
  TGraph abs_ds;
public:

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
private:
  TF1 fit_funct;
  std::string fit_name = "";
  TGraph example;
  bool fitted = false;
  
public:

  FunctionalFit(Func*, const std::string&);
  FunctionalFit() = default;
  void SetFitParameters(const std::vector<double>&);
  void SetFitParameterRanges(const std::vector<std::pair<double, double>>&);
  void PerformFitOnData(TGraph, bool i = false);
  void SetExampleGraph(const TGraph&);
  double GetParameterValue(const int&) const;
  double GetChiSquared() const ;
  TGraph GetGraph() const;
  TGraph GetGraphOfJust(const std::vector<int>&) const;
  TGraph GetGraphExcluding(const std::vector<int>&) const;

};

TGraph DarkSubtractFromTreePtrs(TTree*, TTree*, const int);
TGraph GetDarkSubtractFromFile(const std::string, const std::string, const int&);
TGraph RemoveRegion(const TGraph&, const double&, const double&);

template<class L>
TGraph BinaryOperation(const TGraph&, const TGraph&, const std::string&, const L&);

TGraph CalculateGdAbs(const TGraph&, const TGraph&);
TGraph PWRatio(const TGraph&, const TGraph&);
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
