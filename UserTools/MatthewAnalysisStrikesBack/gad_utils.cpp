#include "gad_utils.h"


LEDInfo::LEDInfo() : calibration_curve_ptr(nullptr) {};

LEDInfo::~LEDInfo(){
  if(calibration_curve_ptr) delete calibration_curve_ptr;
}

LEDInfo::LEDInfo(LEDInfo&& rhs){
  calibration_curve_ptr = rhs.calibration_curve_ptr;
  rhs.calibration_curve_ptr = nullptr;
}

SimplePureFunc::SimplePureFunc(const TGraph& pds)
  : pure_dark_subtracted{pds} {
  fit_min = wave_min;
  fit_max = wave_max;
}

double SimplePureFunc::Evaluate(double* x, double* p){
  return (p[PURE_SCALING] + p[PURE_FIRST_ORDER_CORRECTION] * (x[0]-276) ) * (pure_dark_subtracted.Eval(x[0] - p[PURE_TRANSLATION])) +
          p[LINEAR_BACKGROUND_GRADIENT] * (x[0]-276)
          + p[LINEAR_BACKGROUND_OFFSET];
}


CombinedGdPureFunc::CombinedGdPureFunc(const TGraph& pds, const TGraph& a)
  : pure_dark_subtracted{pds}, spec_abs{a} {
  fit_min = wave_min;
  fit_max = wave_max;
}

double CombinedGdPureFunc::Evaluate(double* x, double* p){
  return (p[PURE_SCALING] + p[PURE_FIRST_ORDER_CORRECTION] * (x[0]-276)) * pure_dark_subtracted.Eval(x[0] - p[PURE_TRANSLATION]) +
    p[ABS_SCALING] * spec_abs.Eval(x[0] - p[ABS_TRANSLATION]) +
    p[FIRST_ORDER_BACKGROUND] * (x[0]-276) +
    p[ZEROTH_ORDER_BACKGROUND];
}

CombinedGdPureFunc_DATA::CombinedGdPureFunc_DATA(const TGraph& pds, const TGraph& a)
  : pure_dark_subtracted{pds}, rat_abs{a} {
  fit_min = wave_min;
  fit_max = wave_max;
}
 
double CombinedGdPureFunc_DATA::Evaluate(double* x, double* p){
  return p[PURE_SCALING]  // p[PURE_FIRST_ORDER_CORRECTION] * x[0]
    * pure_dark_subtracted.Eval(p[PURE_STRETCH] * (x[0] - 276) + 276 - p[PURE_TRANSLATION])
    * (std::max(0.0, (1-p[ABS_SCALING]*rat_abs.Eval(x[0]))) +
       p[SECOND_ORDER_BACKGROUND] *  (x[0]-276) * (x[0]-276) +
       p[FIRST_ORDER_BACKGROUND] * (x[0]-276) +
       p[ZEROTH_ORDER_BACKGROUND]);
//  return p[PURE_SCALING]  // p[PURE_FIRST_ORDER_CORRECTION] * x[0]
//    * pure_dark_subtracted.Eval(p[PURE_STRETCH] * (x[0] - 276) + 276 - p[PURE_TRANSLATION])
//    * ( 1 + p[SECOND_ORDER_BACKGROUND] *  (x[0]-276) * (x[0]-276)
//          + p[FIRST_ORDER_BACKGROUND] * (x[0]-276)
//          + p[ZEROTH_ORDER_BACKGROUND] )
//    * std::max(0.0, (1-p[ABS_SCALING]*rat_abs.Eval(x[0])));
      
      

}
  
// pure + abs that spec sees (0 outside of abs region) + pol1 - wiggled (non-trivial / impossible)

// pure * ratio absorbance (physics) + pol1 (wiggling all of this)
// how do we transform the pure: scale, shift, non-linear transform
// r = pure / gd => gd = pure / r

AbsFunc::AbsFunc(const TGraph& a)
  : abs_ds{a} {
  fit_min = abs_region_low;
  fit_max = abs_region_high;
}

double AbsFunc::Evaluate(double* x, double* p){
  return p[ABS_SCALING] * abs_ds.Eval(x[0] - p[ABS_TRANSLATION]) +
         p[THIRD_BACKGROUND]  * (x[0]-276) * (x[0]-276) * (x[0]-276) + 
         p[SECOND_BACKGROUND] * (x[0]-276) * (x[0]-276) + 
         p[FIRST_BACKGROUND]  * (x[0]-276) +
         p[ZEROTH_BACKGROUND];
}

FunctionalFit::FunctionalFit(Func* func_class_ptr, const std::string& fcn) : fit_name{fcn} {
  if (func_class_ptr->fit_min == 0 ||
      func_class_ptr->fit_max == 0 ||
      func_class_ptr->n_fit == 0){
    throw std::invalid_argument("FunctionalFit::FromFunctionClass: FIT REGION OR FIT DIMENSION NOT VALID!!!");
  }
  
  fit_funct = TF1(fit_name.c_str(),
                  func_class_ptr,
                  &Func::Evaluate,
                  func_class_ptr->fit_min, func_class_ptr->fit_max, func_class_ptr->n_fit,
                  fit_name.c_str(),"Evaluate");
  if (!fit_funct.IsValid()){
    throw std::runtime_error("FunctionalFit::FunctionalFit FIT FUNCTION IS NOT VALID!!!\n");
  }
}

void FunctionalFit::SetFitParameters(const std::vector<double>& pv){
  if (static_cast<int>(pv.size()) != fit_funct.GetNpar()){
    throw std::invalid_argument("FunctionalFit::SetFitParameters: WRONG NUMBER OF INITIAL FIT PARAMETER VALUES!!!\n"); 
  }
  fit_funct.SetParameters(pv.data());
  return;
}

void FunctionalFit::SetFitParameterRanges(const std::vector<std::pair<double, double>>& fpr){
  if (static_cast<int>(fpr.size()) != fit_funct.GetNpar()){
    throw std::invalid_argument("FunctionalFit::SetFitParameterRanges: WRONG NUMBER OF INITIAL FIT PARAMETER RANGES!!!\n"); 
  }
  for (int i = 0; i < fit_funct.GetNpar(); ++i){
    fit_funct.SetParLimits(i, fpr.at(i).first, fpr.at(i).second);
  }
}


TFitResultPtr FunctionalFit::PerformFitOnData(TGraph data, bool interactive){
  fit_funct.SetNpx(10000);
  TFitResultPtr res;
  
  if (interactive){
    int argc=0;
    char* argv=nullptr;
    static TApplication app("app", &argc, &argv);
    TCanvas c1("fiddle_canv", "fiddle", 1280, 1024);
    fit_funct.Draw();
    data.Draw("*L same");
    c1.Modified();
    c1.Update();
    gSystem->ProcessEvents();
    gPad->WaitPrimitive();
    res = data.Fit(&fit_funct, "NRS"); // keep this one
    while(gROOT->FindObject("fiddle_canv") != 0){
      c1.Modified();
      c1.Update();
      gSystem->ProcessEvents();
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
  }
  else {
//    TDirectory* cwd = gDirectory;
//    TFile* fsav=new TFile("fsav.root","RECREATE");
//    data.Write("absgraph_beingfit");
//    fit_funct.Write("fitfunc_beforefit");
    res = data.Fit(&fit_funct, "NRSQ"); // keep this one
//    fit_funct.Write("fitfunc_afterfit");
//    fsav->Close();
//    delete fsav;
//    cwd->cd();
  }
  
//  if (res->IsEmpty() || !res->IsValid() || res->Status() != 0){
//    throw std::runtime_error("FunctionalFit::PerformFitOnData: FIT NOT SUCCESSFUL!!!\n");
//  }
  
  fitted = true;
  return res;
}

void FunctionalFit::SetExampleGraph(const TGraph& e){
  if (e.GetN() == 0){  // != number_of_points){
    throw std::invalid_argument("FunctionalFit::SetExampleGraph: EXAMPLE HAS WRONG NUMBER OF POINTS!!!\n");
  }
  example = e;
}

double FunctionalFit::GetParameterValue(const int& p) const {
  if (!fitted){
    throw std::runtime_error("FunctionalFit::GetGraph - NOT YET FITTED!!!\n");
  }
  return fit_funct.GetParameter(p);
}

double FunctionalFit::GetChiSquared() const {
  if (!fitted){
    throw std::runtime_error("FunctionalFit::GetGraph - NOT YET FITTED!!!\n");
  }
  return fit_funct.GetChisquare();
}

TGraph FunctionalFit::GetGraph() const {
  if (!fitted){
    throw std::runtime_error("FunctionalFit::GetGraph - NOT YET FITTED!!!\n");
  }
  if (example.GetN() == 0){ // != number_of_points){
    throw std::runtime_error("FunctionalFit::GetGraph - NO EXAMPLE FOR POPULATING SET!!!\n");
  }
  TGraph result(example.GetN());
  for (int i = 0; i < example.GetN(); ++i){
    double x = 0, y = 0;
    example.GetPoint(i, x, y);
    result.SetPoint(i, x, fit_funct.Eval(x));
  }
  return result;
}
  
TGraph FunctionalFit::GetGraphOfJust(const std::vector<int>& p) const {
  if (!fitted){
    throw std::runtime_error("FunctionalFit::GetGraphOfJust - NOT YET FITTED!!!\n");
  }
  if (example.GetN() == 0){ // != number_of_points){
    throw std::runtime_error("FunctionalFit::GetGraph - NO EXAMPLE FOR POPULATING SET!!!\n");
  }
  for (const auto param : p){
    if (param > fit_funct.GetNpar() || param < 0){
      throw std::invalid_argument("FunctionalFit::GetGraphOfJust - INVALID PARAMETER CHOICE!!!\n");
    }
  }
  TF1 local_ff = fit_funct;
  for (int i = 0; i < fit_funct.GetNpar(); ++i){
    if (std::find(p.begin(), p.end(), i) == p.end()){
      local_ff.SetParameter(i, 0);
    }
  }
  TGraph result(example.GetN());
  for (int i = 0; i < example.GetN(); ++i){
    double x = 0, y = 0;
    example.GetPoint(i, x, y);
    result.SetPoint(i, x, local_ff.Eval(x));
  }
  return result;
}

TGraph FunctionalFit::GetGraphExcluding(const std::vector<int>& p) const {
  if (!fitted){
    throw std::runtime_error("FunctionalFit::GetGraphExcluding - NOT YET FITTED!!!\n");
  }
  if (example.GetN() == 0){ // != number_of_points){
    throw std::runtime_error("FunctionalFit::GetGraph - NO EXAMPLE FOR POPULATING SET!!!\n");
  }
  for (const auto param : p){
    if (param > fit_funct.GetNpar() || param < 0){
      throw std::invalid_argument("FunctionalFit::GetGraphExcluding - INVALID PARAMETER CHOICE!!!\n");
    }
  }
  TF1 local_ff = fit_funct;
  for (int i = 0; i < fit_funct.GetNpar(); ++i){
    if (std::find(p.begin(), p.end(), i) != p.end()){
      local_ff.SetParameter(i, 0);
    }
  }
  TGraph result(example.GetN());
  for (int i = 0; i < example.GetN(); ++i){
    double x = 0, y = 0;
    example.GetPoint(i, x, y);
    result.SetPoint(i, x, local_ff.Eval(x));
  }
  return result;
}

TGraph DarkSubtractFromTreePtrs(TTree* led_ptr, TTree* dark_ptr, const int dark_entry = -1){
  if (led_ptr == nullptr || dark_ptr ==  nullptr){
    throw std::invalid_argument("DarkSubtractFromTreePtrs: one or more input trees is nullptr!!!\n");
  }

  std::vector<double> led_values, dark_values, wavelengths;
  std::vector<double>* led_values_p = &led_values;
  std::vector<double>* dark_values_p = &dark_values;
  std::vector<double>* wavelengths_p = &wavelengths;
  
  const bool set_branch_success = 
    (led_ptr->SetBranchAddress("value", &led_values_p) == 0) &&
    (led_ptr->SetBranchAddress("wavelength", &wavelengths_p) == 0) &&
    (dark_ptr->SetBranchAddress("value", &dark_values_p) == 0);
  if (!set_branch_success){
    throw std::runtime_error("DarkSubtractFromTreePtrs: error when setting one or more branch addresses!!!\n");
  }
  
  const bool get_entry_success = 
    led_ptr->GetEntry(0) && 
    dark_entry == -1 ? dark_ptr->GetEntry(dark_ptr->GetEntries()-1) : dark_ptr->GetEntry(dark_entry); 
  if (!get_entry_success){
    throw std::runtime_error("DarkSubtractFromTreePtrs: failed to get led or dark entry!!!\n");
  }
  if (static_cast<int>(led_values.size()) != number_of_points ||
      static_cast<int>(dark_values.size()) != number_of_points ||
      static_cast<int>(wavelengths.size()) != number_of_points){
    throw std::runtime_error("DarkSubtractFromTreePtrs: retrieved vectors are incorrect size!!!\n");
  }
  
  TGraph result = TGraph(number_of_points);
  for (int i = 0; i < number_of_points; ++i){
    result.SetPoint(i, wavelengths.at(i), led_values.at(i) - dark_values.at(i));
  }
  return result;
}

TGraph GetDarkSubtractFromFile(const std::string fname, const std::string led_name, const int& dark_offset = -1){
  TFile* file_ptr = TFile::Open(fname.c_str(), "READ");
  if (!file_ptr || !file_ptr->IsOpen() || file_ptr->IsZombie()){
    throw std::invalid_argument("GetDarkSubtractFromFile: FAILED TO OPEN FILE" + fname + "!!!\n");
  }
  TGraph result =  DarkSubtractFromTreePtrs(static_cast<TTree*>(file_ptr->Get(led_name.c_str())),
                                            static_cast<TTree*>(file_ptr->Get(dark_name.c_str())),
                                            dark_offset);
  file_ptr->Close();
  return result;
}

TGraph RemoveRegion(const TGraph& g, const double& low, const double& high){
  /*
  if (g.GetN() != number_of_points){
    throw std::invalid_argument("RemoveRegion: INPUT TGraph HAS INCORRECT NUMBER OF POINTS!!!\n");
  }
  */
  int npts= g.GetN();
  if(npts==0){
    throw std::invalid_argument("RemoveRegion: EMPTY INPUT GRAPH!!!\n");
  }
  double* xarr = g.GetX();
  double* yarr = g.GetY();
  if (xarr[0] > low || xarr[npts-1] < high){
    throw std::invalid_argument("RemoveRegion: INVALID REGION BOUNDS!!!\n");
  }
  
  std::vector<double> xvec, yvec;
  for (int i = 0; i<npts; ++i){
    if( (xarr[i] < low) || (xarr[i] > high) ){
      xvec.push_back(xarr[i]);
      yvec.push_back(yarr[i]);
    }
  }
  TGraph result(xvec.size(), xvec.data(), yvec.data());
  return result;
}

template<class L>
TGraph BinaryOperation(const TGraph& a, const TGraph& b, const std::string& name, const L& op){
  if (a.GetN() != b.GetN()){
    throw std::invalid_argument((name +": INPUT TGraphs HAVE UNEQUAL NUMBER OF POINTS, a has "+a.GetN()+"  b has "+b.GetN()+"!!!\n"));
  }
  const int N = a.GetN();
  TGraph result(N);
  for (int i = 0; i < N; ++i){
    double x = 0, y = 0; 
    a.GetPoint(i, x, y);
    double xa = x;
    b.GetPoint(i, x, y);
    double xb = x;
    result.SetPoint(i, 0.5 * (xa + xb), op(a,b,i));
  }
  return result;
}

TGraph PWRatio(const TGraph& n, const TGraph& d){
  return BinaryOperation(n, d, "Point-Wise Ratio",
                         [](const TGraph& _n, const TGraph& _d, int _i){
                           double x = 0, y = 0;
                           _n.GetPoint(_i, x, y);
                           double ny = y;
                           _d.GetPoint(_i, x, y);
                           double dy = y;
                           return ny / dy;});
}

TGraph PWLogRatio(const TGraph& n, const TGraph& d){
  return BinaryOperation(n, d, "Point-Wise Ratio",
                         [](const TGraph& _n, const TGraph& _d, int _i){
                           double x = 0, y = 0;
                           _n.GetPoint(_i, x, y);
                           double ny = y;
                           _d.GetPoint(_i, x, y);
                           double dy = y;
                           return log10(ny / dy);});
}

TGraph PWMultiply(const TGraph& n, const TGraph& d){
  return BinaryOperation(n, d, "Point-Wise Multiply",
                         [](const TGraph& _n, const TGraph& _d, int _i){
                           double x = 0, y = 0;
                           _n.GetPoint(_i, x, y);
                           double ny = y;
                           _d.GetPoint(_i, x, y);
                           double dy = y;
                           return ny * dy;});
}


TGraph PWPercentageDiff(const TGraph& n, const TGraph& d){
  return BinaryOperation(n, d, "Point-Wise Difference",
                         [](const TGraph& _n, const TGraph& _d, int _i){
                           double x = 0, y = 0;
                           _n.GetPoint(_i, x, y);
                           double ny = y;
                           _d.GetPoint(_i, x, y);
                           double dy = y;
                           return (ny / dy) - 1;});

}

TGraph PWDifference(const TGraph& n, const TGraph& d){
  return BinaryOperation(n, d, "Point-Wise Difference",
                         [](const TGraph& _n, const TGraph& _d, int _i){
                           double x = 0, y = 0;
                           _n.GetPoint(_i, x, y);
                           double ny = y;
                           _d.GetPoint(_i, x, y);
                           double dy = y;
                           return ny - dy;});
}

double SuccessMetric(const TGraph& res){
  TGraph abs_w_cut(res.GetN());
  for (int i = 0; i < res.GetN(); ++i){
    double x = 0, y = 0;
    res.GetPoint(i, x, y);
      if (x > wave_min && x < wave_max){
        abs_w_cut.SetPoint(i, x, abs(y));
      }
    }
    return abs_w_cut.Integral();
}

double GetPeakDiff(const TGraph& abs){
  return abs.Eval(273) - abs.Eval(276);
}

double GetPeakRatio(const TGraph& abs){
  return abs.Eval(273) / abs.Eval(276);
}

TGraph TrimGraph(const TGraph& g, double l, double h){
  std::vector<double> xvec,yvec;
  double* xarr = g.GetX();
  double* yarr = g.GetY();
  for (int i = 0; i < g.GetN(); ++i){
    if( xarr[i] < l ) continue;
    if( xarr[i] > h ) break;
    xvec.push_back(xarr[i]);
    yvec.push_back(yarr[i]);
  }
  TGraph result(xvec.size(),xvec.data(),yvec.data());
  return result;
}

TGraph CalculateGdAbs(const TGraph& a, const TGraph& b){
  // a - pure, b - gd
  std::vector<double> vals;
  TGraph result = TGraph(a.GetN());
  for (int i = 0; i < a.GetN(); ++i){
    double xa = 0, ya = 0;
    a.GetPoint(i, xa, ya);
    if (xa > abs_region_high || xa < abs_region_low){
      result.SetPoint(i, xa, 0);
    }
    else {
      double xb = 0, yb = 0;
      b.GetPoint(i, xb, yb);
      vals.push_back(1 - (yb / ya));
      result.SetPoint(i, xa, 1 - (yb / ya));
    }
  }
  std::vector<double>::iterator scale_it = std::max_element(vals.begin(), vals.end());
  for (int i = 0; i < a.GetN(); ++i){
    double x = 0, y = 0;
    result.GetPoint(i, x, y);
    result.SetPoint(i, x, y / *scale_it);
  }
  return result;
}

/*
  I don't think this gets used at any point, basically this version of ROOT doesn't have a TGraph::Scale() method and after already 
  changing all my GetPointXs to GetPoints I think implementing my own Scale would force me to end my own life - so go lump it.
   - hey matthew you can implement TGraph::GetPointX(i) as TGraph::GetX()[i].
     TGraph::Scale(s) could be implemented as:
     auto scaler = [](TGraph& g, const double& s){ for(size_t i=0; i<g.GetN(); ++i) g.GetY()[i]*=s; return g; };
 */
// TGraph Normalise(TGraph a){
//   if (a.GetN() != a.GetN()){
//     std::invalid_argument("Normalise: INPUT TGraphs HAS WRONG NUMBER OF POINTS!!!\n");
//   }
//   double max = 0;
//   for (int i = 0; i < a.GetN(); ++i){
//     double x = 0, y = 0;
//     a.GetPoint(i, x, y);
//     if (y > max){max = y;}
//   }
//   a.Scale(1/max);
//   return a;
// }

TGraph ZeroNegative(TGraph a){
  TGraph result(a.GetN());
  for (int i = 0; i < a.GetN(); ++i){
    double x = 0, y = 0;
    a.GetPoint(i, x, y);
    result.SetPoint(i, x, abs(y) + 0.01);
  }
  return result;
}

