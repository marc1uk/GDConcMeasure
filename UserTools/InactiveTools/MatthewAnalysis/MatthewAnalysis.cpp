#include "MatthewAnalysis.h"

#include <boost/algorithm/string.hpp>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"

#include <map>
#include <limits>

MatthewAnalysis::MatthewAnalysis():Tool(){}

bool MatthewAnalysis::Initialise(std::string configfile, DataModel &data){

  m_data= &data;
  
  /* - new method, Retrieve configuration options from the postgres database - */
  int RunConfig=-1;
  m_data->vars.Get("RunConfig",RunConfig);
  
  if(RunConfig>=0){
    std::string configtext;
    bool get_ok = m_data->postgres_helper.GetToolConfig(m_unique_name, configtext);
    if(!get_ok){
      Log(m_unique_name+" Failed to get Tool config from database!",v_error,verbosity);
      return false;
    }
    // parse the configuration to populate the m_variables Store.
    if(configtext!="") m_variables.Initialise(std::stringstream(configtext));
    
  }
  
  /* - old method, read config from local file - */
  if(configfile!="")  m_variables.Initialise(configfile);
  
  //m_variables.Print();
  
  
  m_variables.Get("verbosity",verbosity);
  
  
  bool success = RetrieveDarkSubPure();
  if(!success) return success;
  
  // get calibration cofficients from either local file (if present) or DB
  if(m_variables.Has("CalibCoeff0")){
    success = RetrieveCalibrationCurve();
  } else {
    success = RetrieveCalibrationCurveDB();
  }
  
  return success;
}

bool MatthewAnalysis::Execute(){
  
  Log("MatthewAnalysis Executing...",v_debug,verbosity);
  
  // check if we have an Analyse flag in the CStore indicating intensity data
  // is available for fitting
  if (ReadyToAnalyse()){
    
    Log("MatthewAnalysis Processing new measurement...",v_debug,verbosity);
    
    // reinit analysis results to prevent accidental carry-over
    Log("MatthewAnalysis reinitializing variables",v_debug,verbosity);
    ReInit();
    
    // get the TTrees containing LED-on and LED-off traces
    Log("MatthewAnalysis searching for TTrees",v_debug,verbosity);
    std::pair<TTree*, TTree*> dark_and_led_trees = FindDarkAndLEDTrees();
    TTree* dark_tree = dark_and_led_trees.first, *led_tree = dark_and_led_trees.second;
    if(dark_tree==nullptr || led_tree==nullptr){
      std::stringstream tmp;
      tmp << "MatthewAnalysis::FindDarkAndLEDTrees could not find signal ("<<led_tree
          << ") and/or dark ("<<dark_tree<<") TTree in DataModel";
      Log(tmp.str(),v_error,verbosity);
      return false;
    }
    
    
    // pull LED-on and LED-off traces from TTrees and produce TGraph of dark-subtracted trace
    Log("MatthewAnalysis performing dark subtraction",v_debug,verbosity);
    dark_subtracted_gd = PerformDarkSubtraction(led_tree, dark_tree);
    if(dark_subtracted_gd.GetN()==0){
      return false;  // failed to get data from TTrees
    }
    
    // Remove points in the absorption region 270-280nm TODO make this range configurable
    Log("MatthewAnalysis removing absorption region",v_debug,verbosity);
    TGraphErrors dark_sub_gd_w_abs_region_removed = RemoveRegion(dark_subtracted_gd);
    
    // taking the reference pure water trace as a function, fit this function, scaled
    // and with a linear background added, to the measured trace with the absorption region removed.
    // The fit paramters represent the scaling and linear background added.
    Log("MatthewAnalysis performing functional fit",v_debug,verbosity);
    double* fitting_parameters = FunctionalFit(dark_sub_gd_w_abs_region_removed);
    if(fitting_parameters==nullptr){
      return false;   // fit failed
    }
    double* fitting_errors; // TODO can be obtained from TFitResultPtr in FunctionalFit for example
    
    // scale pure water trace by fit function parameters to obtain a version of the pure water
    // spectrum that should match the LED-on spectrum in the lobes of the LED peak
    TGraphErrors dark_subtracted_pure = *(m_data->dark_sub_pure);
    Log("MatthewAnalysis scaling pure reference curve to match data",v_debug,verbosity);
    pure_scaled = Scale(dark_subtracted_pure, fitting_parameters, fitting_errors);
    
    // use the difference between the scaled pure water trace and the current LED-on trace
    // in the gd absorption-peak region to determine the amount of absorbance due to Gd.
    // absorbance is defined as the log_10(received/transmitted light)
    Log("MatthewAnalysis generating absorption curve",v_debug,verbosity);
    absorbance = CalculateAbsorbance(pure_scaled, dark_subtracted_gd);
    
    // fit the two main absorption peaks at 273 and 276 nm with gaussians,
    // then take difference in their amplitudes. As with the amplitudes of each peak,
    // this should be related to concentration.
    Log("MatthewAnalysis calculating difference in heights of absorption peaks",v_debug,verbosity);
    std::pair<double,double> peaksdiff_w_error = FindPeakDifference(absorbance);
    //std::vector<double> peaksdiff_w_error = {0.02 * m_data->concs_and_peakdiffs.size() + 0.005, 0};
    if(peaksdiff_w_error==std::pair<double,double>{0,0}){
      return false;   // failure fitting gaussians to the two absorbance peaks
    }
    
    // use calibration curve to map difference in height of the two main absorbance peaks to gd concentration
    Log("MatthewAnalysis mapping peak height difference to concentration",v_debug,verbosity);
    std::pair<double,double> concentration_w_error = CalculateConcentration(peaksdiff_w_error);
    
    // Inform downstream tools that a new measurement is available
    // maybe we could use the value to indicate if the data is good?
    Log("MatthewAnalysis passing out results",v_debug,verbosity);
    m_data->CStore.Set("NewMeasurement",true);
    
    // debug prints
    std::cout << "Peak difference from " << led_name << " is: " << peaksdiff_w_error.first
              << " +/- " << peaksdiff_w_error.second << std::endl;
    std::cout << "Concentration from " << led_name << " is: " << concentration_w_error.first
              << " +/- " << concentration_w_error.second << std::endl;
    
    // update datamodel
    m_data->CStore.Set("dark_subtracted_gd",reinterpret_cast<intptr_t>(&dark_subtracted_gd));
    m_data->CStore.Set("data_fitresptr",reinterpret_cast<intptr_t>(&fitresptr));
    m_data->CStore.Set("pure_scaled",reinterpret_cast<intptr_t>(&pure_scaled));
    m_data->CStore.Set("absorbance",reinterpret_cast<intptr_t>(&absorbance));
    m_data->CStore.Set("LeftPeakFitFunc",reinterpret_cast<intptr_t>(left_peak));
    m_data->CStore.Set("RightPeakFitFunc",reinterpret_cast<intptr_t>(right_peak));
    m_data->CStore.Set("LeftPeakFitResult",reinterpret_cast<intptr_t>(&lpf));
    m_data->CStore.Set("RightPeakFitResult",reinterpret_cast<intptr_t>(&rpf));
    m_data->concs_and_errs.push_back(concentration_w_error);
    m_data->peakdiffs_and_errs.push_back(peaksdiff_w_error);
    
  } else {
    // else no data to Analyse
    m_data->CStore.Remove("NewMeasurement");
  }
  
  return true;
}

bool MatthewAnalysis::Finalise(){
  
  if(m_data->dark_sub_pure) delete m_data->dark_sub_pure;
  m_data->dark_sub_pure = nullptr;
  if(m_data->pure_fct) delete m_data->pure_fct;
  m_data->pure_fct=nullptr;
  if(m_data->calib_curve) delete m_data->calib_curve;
  m_data->calib_curve=nullptr;
  
  return true;
}

void MatthewAnalysis::ReInit(){
  // Reinitialize measurement outputs, so we don't carry over
  // results from a previous fit if we bail early
  
  // clear TGraphs
  pure_scaled.Set(0);
  dark_subtracted_gd.Set(0);
  absorbance.Set(0);
  
  // remake TFitResultPtrs
  fitresptr = TFitResultPtr();
  lpf = TFitResultPtr();
  rpf = TFitResultPtr();
  
  // TF1s are owned by gROOT so to ensure deletion (for reinitialization) we create them on the heap
  if(left_peak) delete left_peak;
  if(right_peak) delete right_peak;
  
}

bool MatthewAnalysis::RetrieveCalibrationCurve(){
  // Constructs calibration curve from coefficients in the config file, then stores in DataModel
  
  const std::string calib_prefix = "CalibCoeff";
  const int pol_order = 6;
  double calib_coefficients[pol_order + 1];
  bool success = true;
  
  // retrieve calibration coefficients from config file
  for (auto i = 0; i < pol_order + 1; ++i){
    success *= m_variables.Get(calib_prefix + std::to_string(i), calib_coefficients[i]);
    //std::cout << calib_coefficients[i] << std::endl;
  }
  if(!success){
    Log("MatthewAnalysis::RetrieveCalibrationCurve did not find "+std::to_string(pol_order)
        +" calibration constants in config file!",v_error,verbosity);
    return false;
  }
  
  // construct 6th order polynomial to convert intensity into concentration
  TF1* calib_curve = new TF1("calib", "pol 6", 0, 0.4);
  calib_curve->SetParameters(calib_coefficients);
  
  if(!calib_curve->IsValid()){
    Log("MatthewAnalysis::RetrieveCalibrationCurve calibration curve is not valid!",v_error,verbosity);
    return false;
  }
  m_data->calib_curve = calib_curve;
  
  return success;
}

bool MatthewAnalysis::RetrieveCalibrationCurveDB(){
  // pull calibration version from config file
  
  // new calibration curve can be inserted with e.g.:
  // psql -U postgres -d "rundb" -c "INSERT INTO data (timestamp, run, tool, name, values) VALUES ('now()', '0', 'MatthewAnalysis', 'calibration_curve', '{\"version\":0, \"formula\":\"pol6\", \"params\":[0.0079581671, 2.7415760, -31.591756, 478.69924, 18682.891, -29982.759] }' );"
  
  int calib_version;
  bool get_ok = m_variables.Get("calib_version",calib_version);
  if(not get_ok){
    Log("MatthewAnalysis::RetrieveCalibrationCurve failed to find calibration curve "
        "version number in config file",v_error,verbosity);
    return false;
  }
  
  // put in datamodel for database
  m_data->CStore.Set("calib_version",calib_version);
  
  // use version number to lookup formula from database
  std::string formula_str="";
  std::string query_string = "SELECT values->'formula' FROM data WHERE name='calibration_curve'"
                             " AND values->'formula' IS NOT NULL AND values->'version' IS NOT NULL"
                             " AND values->'version'=" + m_data->postgres.pqxx_quote(calib_version);
  bool success = m_data->postgres.ExecuteQuery(query_string, formula_str);
  if(formula_str==""){
    Log("MatthewAnalysis::RetrieveCalibrationCurveDB obtained empty formula "
        "for calibration curve function",v_error,verbosity);
    return false;
  }
  // we must strip enclosing quotations or the TF1 constructor goes berzerk
  formula_str = formula_str.substr(1,formula_str.length()-2);
  
  // another query for the curve parameters
  query_string = "SELECT values->>'params' FROM data WHERE name='calibration_curve'"
                 " AND values->'params' IS NOT NULL' AND values->'version' IS NOT NULL"
                 " AND values->'version'=" + m_data->postgres.pqxx_quote(calib_version);
  std::string params_str;
  success &= m_data->postgres.ExecuteQuery(query_string, params_str);
  
  // check for errors
  if(!success){
    Log("MatthewAnalysis::RetrieveCalibrationCurve failed to retrieve calibration curve "
        "or parameters for version "+std::to_string(calib_version),v_error,verbosity);
    return false;
  }
  
  // the params string is a json array; i.e. '[val1, val2, val3...]'
  // strip the '[' and ']'
  params_str = params_str.substr(1,params_str.length()-2);
  // parse it
  std::stringstream ss(params_str);
  std::string tmp;
  std::vector<double> calib_coefficients;
  while(std::getline(ss,tmp,',')){
    char* endptr = &tmp[0];
    double nextval = strtod(tmp.c_str(),&endptr);
    if(endptr==&tmp[0]){
      Log("MatthewAnalysis::RetrieveCalibrationCurveDB failed to parse calibration parameters from string '"
          +params_str+"' for calibration curve version "+std::to_string(calib_version),v_error,verbosity);
      return false;
    }
    calib_coefficients.push_back(nextval);
  }
  if(calib_coefficients.size()==0){
    Log("MatthewAnalysis::RetrieveCalibrationCurve parsed no parameters from string '"
        +params_str+" for calibration curve version "+std::to_string(calib_version),v_error,verbosity);
    return false;
  }
  
  // construct a TF1 to convert intensity into concentration
  TF1* calib_curve = new TF1("calib", formula_str.c_str(), 0, 0.4);
  calib_curve->SetParameters(calib_coefficients.data());
  if(!calib_curve->IsValid()){
    Log("MatthewAnalysis::RetrieveCalibrationCurve calibration curve constructed from formula '"
        +formula_str+"' and parameters '"+params_str+"' is not valid!",v_error,verbosity);
    return false;
  }
  
  m_data->calib_curve = calib_curve;
  
  return success;
}

bool MatthewAnalysis::RetrieveDarkSubPure(){
  // Retrieves pure trace from .root file then stores it in the DataModel
  
  // pull the filename of the pure water reference trace from the config file
  std::string pure_filename;
  m_variables.Get("pure_filename", pure_filename);
  // put pure water reference filename into CStore for recording in db.
  // TODO should the pure water reference data be stored in the database, rather than a file?
  m_data->CStore.Set("pure_filename", pure_filename);
  
  // open the reference file
  TFile pure_file(pure_filename.c_str(), "READ");
  if (! pure_file.IsOpen()){
    Log("MatthewAnalysis::RetrieveDarkSubPure failed to open pure file '"
        +pure_filename+"'",v_error,verbosity);
    return false;
  }
  
  //  retrieve the pure water trace (LED spectrum) as a TGraph
  TGraphErrors* dark_subtracted_pure = (TGraphErrors*) pure_file.Get("Graph");
  if(dark_subtracted_pure==nullptr){
    Log("MatthewAnalysis::RetrieveDarkSubPure did not find 'Graph' in pure file '"
        +pure_filename+"'",v_error,verbosity);
    return false;
  }
  
  // minor optimization; if the data in the TGraph is sorted by x value
  //(which it should be for us), then setting the following option can speed up Eval() calls
  dark_subtracted_pure->SetBit(TGraph::kIsSortedX);
  
  m_data->dark_sub_pure = dark_subtracted_pure;
  
  // construct functional fit. We'll scale the reference pure water, dark subtracted intensity
  // and fit it to the measured data in the LED peak lobes.
  // This allows measurement of the absorption line depth in a way that is less sensitive
  // to variations in absolute scale of measured intensity (e.g. from LED fluctuations)
  const int wave_min = 260, wave_max = 300, numb_of_fitting_parameters = 3;
  TF1* pure_fct = new TF1("pure_fct",
          // TF1 wraps a c++ lambda function that evaluates the reference TGraph at a given wavelength
          [dark_subtracted_pure](double* x, double* par){
          return (par[0] * dark_subtracted_pure->Eval(x[0])) + (par[1] * x[0]) + par[2];},
       wave_min, wave_max, numb_of_fitting_parameters);
  m_data->pure_fct = pure_fct;
  
  return true;
}

bool MatthewAnalysis::ReadyToAnalyse(){
  // Checks if analyse flag has been set by scheduler.
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

std::pair<TTree*, TTree*> MatthewAnalysis::FindDarkAndLEDTrees(){
  // Finds the trees for both dark and led traces; then sorts them into a (dark,led)-ordered pair.
  
  std::pair<TTree*, TTree*> result{nullptr,nullptr};
  
  // this tool only extracts concentration from one trace (one LED) at a time;
  // get the name of the LED whose trace we're analysing
  bool get_ok = m_data->CStore.Get("ledToAnalyse",led_name);
  if(not get_ok){
    Log("MatthewAnalysis::FindDarkAndLEDTrees - Analyse flag found in CStore "
        "but no ledToAnalyse",v_error,verbosity);
    return result;
  }
  
  // loop over TTrees (traces) in DataModel
  for (std::pair<const std::string, TTree*>& atree : m_data->m_trees){
    if (atree.first == led_name){result.second = atree.second;}                  // signal trace
    else if (boost::iequals(atree.first, "dark")){result.first = atree.second;}  // dark trace
  }
  
  return result;
}

TGraphErrors MatthewAnalysis::PerformDarkSubtraction(TTree* ledTree, TTree* darkTree){
  // Takes the dark and led trees and returns a TGraphErrors of the dark subtracted trace
  
  // retrieve data from TTrees
  std::vector<double> *led_values= nullptr, *led_wavelengths = nullptr, *led_errors = nullptr;
  std::vector<double>* dark_values = nullptr, *dark_wavelengths = nullptr, *dark_errors = nullptr;
   
  Log("MatthewAnalysis setting branch addresses to retrieve data",v_debug,verbosity);
  bool failure=false;
  failure |= ((ledTree->SetBranchAddress("value", &led_values)) < 0);
  failure |= ((ledTree->SetBranchAddress("wavelength", &led_wavelengths)) < 0);
  failure |= ((ledTree->SetBranchAddress("error", &led_errors)) < 0);
  if(failure){
    Log("MatthewAnalysis failed to set branch addresses for "+led_name+" tree!",v_error,verbosity);
    return false;
  }
  failure |= ((darkTree->SetBranchAddress("value", &dark_values)) < 0);
  failure |= ((darkTree->SetBranchAddress("wavelength", &dark_wavelengths)) < 0);
  failure |= ((darkTree->SetBranchAddress("error", &dark_errors)) < 0);
  if(failure){
    Log("MatthewAnalysis failed to set branch addresses for dark tree!",v_error,verbosity);
    return false;
  }
  
  Log("MatthewAnalysis getting TTree entries",v_debug,verbosity);
  ledTree->GetEntry(0);
  darkTree->GetEntry(darkTree->GetEntries()-1);
  
  // check we got the data we expected
  const int number_of_points = 2068;    // expected trace size
  if( led_values==nullptr       || led_values->size()       != number_of_points ||
      led_wavelengths==nullptr  || led_wavelengths->size()  != number_of_points ||
      led_errors==nullptr       || led_errors->size()       != number_of_points ||
      dark_values==nullptr      || dark_values->size()      != number_of_points ||
      dark_wavelengths==nullptr || dark_wavelengths->size() != number_of_points ||
      dark_errors==nullptr      || dark_errors->size()      != number_of_points ){
      Log("MatthewAnalysis::PerformDarkSubtraction failed to get data from TTrees!",v_error,verbosity);
      return TGraphErrors{};
  }
  
  // construct a TGraphErrors of dark-subtracted data
  TGraphErrors result(number_of_points);
  for (int i = 0; i < number_of_points; ++i){
    result.SetPoint(i, led_wavelengths->at(i), led_values->at(i) - dark_values->at(i));
    result.SetPointError(i, 0, sqrt(pow(led_errors->at(i), 2) + pow(dark_errors->at(i),2)));
  }
  
  // for stability monitoring we'll record some characteristic information about the raw data
  // in the database. The dark trace should be pretty flat, so we'll histogram it,
  // fit it with a gaussian, and record the mean and sigma.
  TH1D* tmphist = new TH1D("tmphist","title",100,*std::min_element(dark_values->begin(), dark_values->end()),
                                                 *std::max_element(dark_values->begin(), dark_values->end()));
  for(int i=0; i<number_of_points; ++i){
    tmphist->Fill(dark_values->at(i));
  }
  std::string options = (verbosity<v_debug) ? "Q" : "";
  tmphist->Fit("gaus",options.c_str());
  double dark_mean = tmphist->GetFunction("gaus")->GetParameter(1);
  double dark_sigma = tmphist->GetFunction("gaus")->GetParameter(2);
  m_data->CStore.Set("dark_mean",dark_mean);
  m_data->CStore.Set("dark_sigma",dark_sigma);
  // for the raw LED-on trace we'll record the maximum and minimum value of the trace
  // within the absorption region of 270-280nm.
  double led_on_max = 0;
  double led_on_min = std::numeric_limits<double>::max();
  for(int i=0; i<number_of_points; ++i){
    if(led_wavelengths->at(i) > 270 && led_wavelengths->at(i) < 280){
      if(led_values->at(i) > led_on_max) led_on_max = led_values->at(i);
      if(led_values->at(i) < led_on_min) led_on_min = led_values->at(i);
    }
  }
  m_data->CStore.Set("raw_absregion_max",led_on_max);
  m_data->CStore.Set("raw_absregion_min",led_on_min);
  
  // cleanup
  delete tmphist;
  delete led_values;
  delete led_wavelengths;
  delete led_errors;
  delete dark_values;
  delete dark_wavelengths;
  delete dark_errors;
  
  return result;
}

TGraphErrors MatthewAnalysis::RemoveRegion(TGraphErrors& trace){
  // Returns a copy of the input trace with datapoints in the absorption region removed.
  
  const auto lower_bound = 270, upper_bound = 280;
  m_data->CStore.Set("gd_fit_range",std::pair<double,double>{lower_bound,upper_bound});
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

double* MatthewAnalysis::FunctionalFit(TGraphErrors& trace){
  // Fits the trace (with absorption region removed) to the pure water trace, and returns the fit parameters
  
  // get the pure-water reference function
  TF1* pure_fct = m_data->pure_fct;
  if(pure_fct==nullptr){
    Log("MatthewAnalysis::FunctionalFit obtained nullptr from DataModel for "
        "functional fit TF1",v_error,verbosity);
    return nullptr;
  }
  
  // do the fit
  std::string options = (verbosity<v_debug) ? "RSQ" : "RS";
  fitresptr = trace.Fit("pure_fct", options.c_str());
  if(fitresptr->IsEmpty() || !fitresptr->IsValid() || fitresptr->Status()!=0){
    Log("MatthewAnalysis::FunctionalFit fit status indicates fit failed!",v_error,verbosity);
    // see https://root.cern.ch/doc/master/classTH1.html#a7e7d34c91d5ebab4fc9bba3ca47dabdd
    // status section for an interpretation of int(res) values
    return nullptr;
  }
  
  // get fit parameters
  double* fitting_parameters = pure_fct->GetParameters();
  return fitting_parameters;
}
 
TGraphErrors MatthewAnalysis::Scale(TGraphErrors& pure_spectrum, const double* par, const double* par_errors){
  // Creates a TGraphErrors of the pure trace that has been scaled by parameters from a functional fit.
  
  TGraphErrors result(pure_spectrum.GetN());
  
  double x_temp, y_temp;
  
  // scale pure spectrum by fit function parameters
  for (auto i = 0; i < pure_spectrum.GetN(); ++i){
    pure_spectrum.GetPoint(i, x_temp, y_temp);
    result.SetPoint(i, x_temp, par[0] * y_temp + par[1] * x_temp + par[2]);
    result.SetPointError(i, pure_spectrum.GetErrorX(i), sqrt(pow(par[0] * pure_spectrum.GetErrorY(i), 2) + pow(par[1] * pure_spectrum.GetErrorX(i),2)));
  }
  
  return result;
}

TGraphErrors MatthewAnalysis::CalculateAbsorbance(TGraphErrors& transmitted, TGraphErrors& received){
  // Creates absorbance spectrum; ie A = log10(received_light / transmitted_light)
  
  TGraphErrors result(transmitted.GetN());
  double trans_x, trans_y, rec_x, rec_y;
  
  const float LOG10_OF_E = log10(exp(1));
  
  for (auto i = 0; i < transmitted.GetN(); ++i){
    transmitted.GetPoint(i, trans_x, trans_y);
    received.GetPoint(i, rec_x, rec_y);
    if (trans_x > 260 && trans_x < 300 && (trans_y/rec_y) > 0){
      result.SetPoint(i, trans_x, log10(trans_y / rec_y ));
      result.SetPointError(i, transmitted.GetErrorX(i), LOG10_OF_E * sqrt(pow((transmitted.GetErrorY(i) / trans_y), 2) + pow(received.GetErrorY(i) / rec_y, 2)));
    }
    else {
      result.SetPoint(i, trans_x, 0);
    }
  }
  return result;
}

std::pair<double,double> MatthewAnalysis::FindPeakDifference(TGraphErrors& absorbance_g){
  // Fit Gaussians to the two largest peaks in absorbance and calculate the difference in their heights
  // Since the height of both peaks are linearly related, taking their difference shares the same
  // relation to concentration as peak height, while helping to remove any remaining baseline offset
  
  // fit twin peaks  TODO make fit ranges / functions configurable
  // re-create the fit functions (probably overkill but ensures they're reset)
  left_peak  = new TF1("leftPeak",  "gaus", 272, 274);
  right_peak = new TF1("rightPeak", "gaus", 275, 277);
  
  std::string options = (verbosity<v_debug) ? "0RSQ" : "0RS";
  lpf = absorbance_g.Fit("leftPeak", options.c_str());
  rpf = absorbance_g.Fit("rightPeak", options.c_str());
  if(lpf->IsEmpty() || !lpf->IsValid() || int(lpf)!=0 ||
     rpf->IsEmpty() || !rpf->IsValid() || int(rpf)!=0){
    Log("MatthewAnalysis::FindPeakDifference fit to absorption peaks failed!",v_error,verbosity);
    return std::pair<double,double>{0,0};
  }
  
  // calculate difference and error
  std::pair<double,double> result;
  result.first = (left_peak->GetParameter(0) - right_peak->GetParameter(0));
  result.second = (sqrt(pow(left_peak->GetParError(0),2) + pow(right_peak->GetParError(0),2)));
  
  return result;
}

std::pair<double,double> MatthewAnalysis::CalculateConcentration(const std::pair<double,double>& peaks_info){
  // Uses the calibration curve to calculate the concentration corresponding to the measured peak difference
  
  // get calibration curve
  TF1* calib_curve = m_data->calib_curve;
  if(calib_curve==nullptr){
    Log("MatthewAnalysis::CalculateConcentration calibration curve is null!",v_error,verbosity);
    return std::pair<double,double>{0,0};
  }
  
  std::pair<double,double> result;
  
  // solve for concentration (x) from absorbance (y), with 0.01 < x < 0.21
  result.first= calib_curve->GetX(peaks_info.first, 0.01, 0.21);
  
  // estimate the error on measured concentration as the difference in concentrations
  // corresponding to measured absorbance +- the error on absorbance.
  result.second = (calib_curve->GetX(peaks_info.first + peaks_info.second, 0.01, 0.21) -
                   calib_curve->GetX(peaks_info.first - peaks_info.second, 0.01, 0.21));
  
  return result;
}
