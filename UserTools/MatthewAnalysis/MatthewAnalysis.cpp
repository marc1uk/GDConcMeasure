#include "MatthewAnalysis.h"

#include <boost/algorithm/string.hpp>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TF1.h"

#include <map>
MatthewAnalysis::MatthewAnalysis():Tool(){}

bool MatthewAnalysis::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();
  m_data= &data;

  RetrieveCalibrationCurve();

  try {
    RetrieveDarkSubPure();
  } catch (std::exception &ex){
    std::cout << ex.what() << '\n';
    return false;
  }
    
  return true;
}

bool MatthewAnalysis::Execute(){

  if (ReadyToAnalyse()){
    
    std::pair<TTree*, TTree*> sorted_trees = SortTrees();
    TTree* dark_tree = sorted_trees.first, *led_tree = sorted_trees.second;
         
    TGraphErrors dark_subtracted_gd = PerformDarkSubtraction(led_tree, dark_tree);
    TGraphErrors dark_sub_gd_w_abs_region_removed = RemoveRegion(dark_subtracted_gd);
      
    double* fitting_parameters = FunctionalFit(dark_sub_gd_w_abs_region_removed);
    double* fitting_errors; //fix this
    
    TGraphErrors dark_subtracted_pure = *(m_data->dark_sub_pure);
    TGraphErrors pure_scaled = Scale(dark_subtracted_pure, fitting_parameters, fitting_errors);      
      
    TGraphErrors absorbance = CalculateAbsorbance(pure_scaled, dark_subtracted_gd);
    std::vector<double> peaks_w_errors = FindPeakDifference(absorbance);
    std::vector<double> concentration_w_errors = CalculateConcentration(peaks_w_errors);
    
    std::cout << "Peak difference from " << "ledName" << " is: " << peaks_w_errors.at(0) << " +/- " << peaks_w_errors.at(1) << '\n';
    std::cout << "Concentration from " << "ledName" << " is: " << concentration_w_errors.at(0) << " +/- " << concentration_w_errors.at(1) << '\n';
  }
  return true;
}
  
bool MatthewAnalysis::Finalise(){

  delete m_data->dark_sub_pure;
  delete m_data->pure_fct;
  delete m_data->calib_curve;
  
  return true;
}

void MatthewAnalysis::RetrieveCalibrationCurve(){
  const std::string calib_prefix = "calib_coeff";
  double calib_coefficients[6];
  for (auto i = 0; i < 6; ++i){m_variables.Get(calib_prefix + std::to_string(i), calib_coefficients[i]);}
  TF1* calib_curve = new TF1("calib", "pol 6", 0, 0.4);
  calib_curve->SetParameters(calib_coefficients);
  m_data->calib_curve = calib_curve;
}

void MatthewAnalysis::RetrieveDarkSubPure(){

  std::string pure_filename;
  m_variables.Get("pure_filename", pure_filename); 
  
  TFile pure_file(pure_filename.c_str(), "READ");
  if (pure_file.IsOpen()){
    TGraphErrors* dark_subtracted_pure = new TGraphErrors(2068);
    dark_subtracted_pure = (TGraphErrors*) pure_file.Get("Graph");
    m_data->dark_sub_pure = dark_subtracted_pure;


    
    const int wave_min = 260, wave_max = 300, numb_of_fitting_parameters = 3;
    TF1* pure_fct = new TF1("pure_fct",
			    [dark_subtracted_pure](double* x, double* par){
			      return (par[0] * dark_subtracted_pure->Eval(x[0])) + (par[1] * x[0]) + par[2];},
			    wave_min, wave_max, numb_of_fitting_parameters);
    m_data->pure_fct = pure_fct;
  }
  else {throw std::runtime_error(pure_filename + " could not be opened.\n");}
}

bool MatthewAnalysis::ReadyToAnalyse(){
  bool ready = false;
  std::string analyse;
  bool success = m_data->CStore.Get("Analyse", analyse);
  if (success && analyse == "Analyse"){
    m_data->CStore.Remove("Analyse");
    analyse.clear();
    ready = true;
  }
  return ready;
}

std::pair<TTree*, TTree*> MatthewAnalysis::SortTrees(){
  std::pair<TTree*, TTree*> result;

  std::string led_name;
  m_data->CStore.Get("ledToAnalyse",led_name);

  for (const auto [name, tree] : m_data->m_trees){
    if (name == led_name){result.second = tree;}
    else if (boost::iequals(name, "dark")){result.first = tree;}
  }    
  return result;
}

TGraphErrors MatthewAnalysis::PerformDarkSubtraction(TTree* ledTree, TTree* darkTree){
  std::vector<double> *led_values= nullptr, *led_wavelengths = nullptr, *led_errors = nullptr;
  std::vector<double>* dark_values = nullptr, *dark_wavelengths = nullptr, *dark_errors = nullptr; 
   
  ledTree->SetBranchAddress("value", &led_values);
  darkTree->SetBranchAddress("value", &dark_values);
  ledTree->SetBranchAddress("wavelength", &led_wavelengths);
  darkTree->SetBranchAddress("wavelength", &dark_wavelengths);
  ledTree->SetBranchAddress("error", &led_errors);
  darkTree->SetBranchAddress("error", &dark_errors);
   
  ledTree->GetEntry(0);
  darkTree->GetEntry(0);

  const int number_of_points = 2068;  
  TGraphErrors result(number_of_points);
  for (auto i = 0; i < number_of_points; ++i){
    result.SetPoint(i, led_wavelengths->at(i), led_values->at(i) - dark_values->at(i));
    result.SetPointError(i, 0, sqrt(pow(led_errors->at(i), 2) + pow(dark_errors->at(i),2)));
  }
  return result;
}

double* MatthewAnalysis::FunctionalFit(TGraphErrors& trace){
  TF1* pure_fct = m_data->pure_fct;
  trace.Fit("pure_fct", "R");
  double* fitting_parameters = pure_fct->GetParameters();
  return fitting_parameters;
}


TGraphErrors MatthewAnalysis::RemoveRegion(TGraphErrors& trace){  
  const auto lower_bound = 270, upper_bound = 280;
  double temp_x, temp_y, temp_err_x, temp_err_y; 
  
  TGraphErrors result(trace.GetN());
  for (auto i = 0; i < result.GetN(); ++i){
    trace.GetPoint(i, temp_x, temp_y);
    temp_err_x = trace.GetErrorX(i);
    temp_err_y = trace.GetErrorY(i);
    
    if (temp_x < lower_bound || temp_x > upper_bound){
      result.SetPoint(i, temp_x, temp_y);
      result.SetPointError(i, temp_err_x, temp_err_y);
    }
    else {result.RemovePoint(i);}
  }                 
  return result;
}
 
TGraphErrors MatthewAnalysis::Scale(TGraphErrors& pure_spectrum, const double* par, const double* par_errors){
  TGraphErrors result(pure_spectrum.GetN());
  
  double x_temp, y_temp;

  for (auto i = 0; i < pure_spectrum.GetN(); ++i){
    pure_spectrum.GetPoint(i, x_temp, y_temp);
    result.SetPoint(i, x_temp, par[0] * y_temp + par[1] * x_temp + par[2]);
    result.SetPointError(i, pure_spectrum.GetErrorX(i), sqrt(pow(par[0] * pure_spectrum.GetErrorY(i), 2) + pow(par[1] * pure_spectrum.GetErrorX(i),2)));
  }
  return result;
}
 
TGraphErrors MatthewAnalysis::CalculateAbsorbance(TGraphErrors& transmitted, TGraphErrors& received){
  TGraphErrors result(transmitted.GetN());
  double trans_x, trans_y, rec_x, rec_y;

  const float LOG10_OF_E = log10(exp(1));
  
  for (auto i = 0; i < transmitted.GetN(); ++i){
    transmitted.GetPoint(i, trans_x, trans_y);
    received.GetPoint(i, rec_x, rec_y);
    if (trans_x > 260 && trans_x < 300){ //get these from config too
      result.SetPoint(i, trans_x, log10(trans_y / rec_y ));
      result.SetPointError(i, transmitted.GetErrorX(i), LOG10_OF_E * sqrt(pow((transmitted.GetErrorY(i) / trans_y), 2) + pow(received.GetErrorY(i) / rec_y, 2)));
    }
    else {
      result.SetPoint(i, trans_x, 0);
    }
  }     
  return result;
}   
 
 std::vector<double> MatthewAnalysis::FindPeakDifference(TGraphErrors& absorbance){
   TF1 left_peak = TF1("leftPeak", "gaus", 272, 274), right_peak = TF1("rightPeak", "gaus", 275, 277);
   absorbance.Fit("leftPeak", "0R");
   absorbance.Fit("rightPeak", "0R");

   std::vector<double> result;
   result.push_back(left_peak.GetParameter(0) - right_peak.GetParameter(0));
   result.push_back(sqrt(pow(left_peak.GetParError(0),2) + pow(right_peak.GetParError(0),2)));
   return result;
}

std::vector<double> MatthewAnalysis::CalculateConcentration(const std::vector<double> peaks_info){
  std::vector<double> result;  
  TF1 calib_curve = *(m_data->calib_curve);
  result.push_back(calib_curve.GetX(peaks_info.at(0), 0.01, 0.21));
  result.push_back(calib_curve.GetX(peaks_info.at(0) + peaks_info.at(1), 0.01, 0.21) - calib_curve.GetX(peaks_info.at(0) - peaks_info.at(1), 0.01, 0.21));
  return result;
}
