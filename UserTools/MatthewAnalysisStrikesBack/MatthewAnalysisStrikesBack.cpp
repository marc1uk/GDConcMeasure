#include "MatthewAnalysisStrikesBack.h"

#include "gad_utils.h"
#include "TFile.h"
#include "TROOT.h"
//#include "TDirectory.h"
#include <memory>
#include <numeric>

MatthewAnalysisStrikesBack::MatthewAnalysisStrikesBack():Tool(){}

bool MatthewAnalysisStrikesBack::Initialise(std::string configfile, DataModel &data){
  
  if(configfile!="")  m_variables.Initialise(configfile);
  m_variables.Get("verbosity",m_verbose);
  int save_debug = 0;
  m_variables.Get("savedebug",save_debug);
  //m_variables.Print();
  
  m_data = &data;
  
  if(save_debug){
    TDirectory* fcur = gDirectory;
    fdebug = new TFile("NewMatthewAnalysisDebug.root","RECREATE");
    fcur->cd();
  }
  
  for (const auto& led_name : {"275_A", "275_B"}){
    Log(m_unique_name+"::Initialise - Getting info for LED: "+led_name, v_debug,m_verbose);
    LEDInfo& led_info = led_info_map[led_name];
    
    //retrieve reference pure and non-pure 
    Log(m_unique_name+" getting pure reference trace for LED "+led_name, v_debug,m_verbose);
    const TGraph pure_ds = TrimGraph(GetPure(led_name));
    SaveDebug(&pure_ds, std::string{"g_ref_pure_"}+led_name);
    Log(m_unique_name+" getting high concentration reference trace for LED "+led_name, v_debug,m_verbose);
    const TGraph high_conc_ds = TrimGraph(GetHighConc(led_name));
    SaveDebug(&high_conc_ds, std::string{"g_ref_highconc_"}+led_name);
    
    // get calibration curve
    Log(m_unique_name+" getting calibration curve for LED: "+led_name, v_debug,m_verbose);
    led_info.calibration_curve_ptr = GetCalibrationCurve(led_name);
    
    //create FunctionalFit object to do initial "simple", ie remove Gd region and fit sidebands
    std::shared_ptr<SimplePureFunc> simple_pure_func = std::make_shared<SimplePureFunc>(pure_ds);
    FunctionalFit simple_fit = FunctionalFit(simple_pure_func.get(), "simple");
    simple_fit.SetExampleGraph(pure_ds);  // i think this is only used for x values...
    
    // ROOT needs starting values and limits to stand any chance of doing the right thing
    // pure scaling, pure first order correction, pure translation, linear background gradient, constant y offset
    std::vector<double> pars{1,0,0,0,0};
    simple_fit.SetFitParameters(pars);
    std::vector<std::pair<double,double>> limits{ {0.01,100},{-1,1},{-10,10},{-10,10},{-1000,1000} };
    simple_fit.SetFitParameterRanges(limits);
    
    //fit non-pure with the simple fit and extract gd attenuation and ratio absorbance shapes
    Log(m_unique_name+" doing simple fit of pure reference on high conc reference for LED: "+led_name, v_debug,m_verbose);
    simple_fit.PerformFitOnData(RemoveRegion(high_conc_ds, abs_region_low, abs_region_high));
    
    const TGraph simple_fit_result = simple_fit.GetGraph();
    SaveDebug(&simple_fit_result, std::string{"g_ref_simple_fit_"}+led_name);
    Log(m_unique_name+" extracting gd shape for LED: "+led_name, v_debug,m_verbose);
    const TGraph gd_abs_attenuation = CalculateGdAbs(simple_fit_result, high_conc_ds);
    SaveDebug(&gd_abs_attenuation, std::string{"g_ref_atten_"}+led_name);
    const TGraph ratio_absorbance = PWRatio(simple_fit_result, high_conc_ds);
    SaveDebug(&ratio_absorbance, std::string{"g_ref_abs_"}+led_name);
    
    //create FunctionalFit object to do inital pure removal fit
    led_info.combined_func = std::make_unique<CombinedGdPureFunc_DATA>(pure_ds, gd_abs_attenuation);
    //led_info.combined_func = std::make_unique<CombinedGdPureFunc_DATA>(pure_ds, ratio_absorbance);
    // gd_abs_attenuation == (1-ratio_absorbance). CombinedGdPureFunc_DATA should use gd_abs_attenuation,
    // but shouldn't there be a log10 somewhere in here?
    led_info.combined_fit = FunctionalFit(led_info.combined_func.get(), "combined");
    led_info.combined_fit.SetExampleGraph(pure_ds);
    
    // y_scaling, x_translation, x_scaling, absorbance_scaling, 2nd order bg, 1st order bg, constant bg
    std::vector<std::string> parnames{"pure_scaling",
                                      "pure_translation",
                                      "pure_stretch",
                                      "abs_scaling",
                                      "second_order_background",
                                      "first_order_background",
                                      "zeroth_order_background"};
    for(int i=0; i<parnames.size(); ++i) led_info.combined_fit.fit_funct.SetParName(i, parnames.at(i).c_str());
    //pars = std::vector<double>{1,0,1,0.5,0,0,0};
    //pars = std::vector<double>{9,-0.25,0.8,0,0.001,-0.06,0.12}; << initial for set 1
    pars = std::vector<double>{20.222116, 0.065019, 0.809043, 0.000000, 0.001000, -0.048069, -0.310535}; // << last of set 1
    led_info.combined_fit.SetFitParameters(pars);
    limits = std::vector<std::pair<double,double>>{ {0.1,50},{-1,1},{0.7,1.1},{0,1.1},{0.002,0.002}, {-0.1,0.1}, {-1,1} };
    led_info.combined_fit.SetFitParameterRanges(limits);
    
    //create FunctionalFit object to do absorbance fit on data
    led_info.abs_func = std::make_shared<AbsFunc>(ratio_absorbance);
    led_info.absorbtion_fit = FunctionalFit(led_info.abs_func.get(), "abs");
    led_info.absorbtion_fit.SetExampleGraph(pure_ds);
    
    // y_scaling, x_translation, 3rd order bg, 2nd order bg, 1st order bg, constant bg
    //pars = std::vector<double>{1,0,0,0,0,0}; << initial for set 1
    pars = std::vector<double>{0.328834, -0.043217, 0.000000, 0.000000, 0.000966, 0.406963}; // << last of set 1
    led_info.absorbtion_fit.SetFitParameters(pars);
    limits = std::vector<std::pair<double,double>>{ {0.01,100},{-10,10},{-2,2},{10,10},{-100,100},{-1000,1000} };
    led_info.absorbtion_fit.SetFitParameterRanges(limits);
    
  }
  
  return true;
}

int MatthewAnalysisStrikesBack::SaveDebug(const TObject* obj, const std::string& name){
  if(fdebug==nullptr) return 1;
  TDirectory* fcur = gDirectory;
  fdebug->cd();
  int byteswritten = obj->Write(name.c_str());
  fcur->cd();
  return byteswritten;
}

bool MatthewAnalysisStrikesBack::Execute(){
  
  double* xv;
  
  // clear any previous results
  Log(m_unique_name+" resetting DataModel",v_debug,m_verbose);
  ReinitDataModel();
  
  // waits for analysis flag
  if (!ReadyToAnalyse()){
    // no data to Analyse
    // see if there's an old flag from this instance and remove it if so
    if(m_data->CStore.Get("NewMatthewAnalyse",current_led)){
      m_data->CStore.Remove("NewMatthewAnalyse");
    }
    return true;
  }
  
  // for debug
  static int measurementnum=0;
  get_ok = m_data->CStore.Get("dbmeasurementnum",measurementnum);
  if(!get_ok) ++measurementnum;

  // get the name of the name of the current led
  Log(m_unique_name+" getting LED to analyse", v_debug,m_verbose);
  GetCurrentLED();
  // from the currently taken traces, sort out the dark and led trees
  Log(m_unique_name+" getting dark & LED tree "+current_led, v_debug,m_verbose);
  GetDarkAndLEDTrees();

  double dark_mean=0, dark_sigma=0;
  Log(m_unique_name+" getting dark trace parameters", v_debug,m_verbose);
  GetDarkTraceParams(dark_mean, dark_sigma);
  m_data->CStore.Set("dark_mean",dark_mean);
  m_data->CStore.Set("dark_sigma",dark_sigma);
  
  double raw_min=0, raw_max=0;
  Log(m_unique_name+" getting led-on trace parameters", v_debug,m_verbose);
  GetRawMinMax(raw_min, raw_max);
  m_data->CStore.Set("raw_absregion_min",raw_min);
  m_data->CStore.Set("raw_absregion_max",raw_max);

  // calculate dark subtract from these traces
  Log(m_unique_name+" doing dark subtraction", v_debug,m_verbose);
  const TGraph darksubdata = DarkSubtractFromTreePtrs(led_tree_ptr, dark_tree_ptr);
  SaveDebug(&darksubdata, std::string{"darksub_"}+current_led+"_"+std::to_string(measurementnum));
  
  // update the datamodel
  Log(m_unique_name+" recording abs region",v_debug,m_verbose);
  GetAbsRegion(darksubdata, dark_subtracted_data_in);
  m_data->CStore.Set("dark_subtracted_data_in",reinterpret_cast<intptr_t>(&dark_subtracted_data_in));

  Log(m_unique_name+" recording sideband region",v_debug,m_verbose);
  GetSidebandRegion(darksubdata, dark_subtracted_data_out);
  m_data->CStore.Set("dark_subtracted_data_out",reinterpret_cast<intptr_t>(&dark_subtracted_data_out));
  
  const TGraph current_dark_sub = TrimGraph(darksubdata);
  SaveDebug(&current_dark_sub, std::string{"datatrimmed_"}+current_led+"_"+std::to_string(measurementnum));
  
  // retrieve the functional fit for this led
  LEDInfo& current_led_info = led_info_map.at(current_led);
  FunctionalFit curr_comb_fit = current_led_info.combined_fit;

  // perform the fit on this new data
  Log(m_unique_name+" fitting data", v_debug,m_verbose);
  
  // y_scaling, x_translation, x_scaling, absorbance_scaling, 2nd order bg, 1st order bg, constant bg
  //TFitResultPtr datafitresptr = curr_comb_fit.PerformFitOnData(current_dark_sub,true);
  // FIXME HACK: the fit's not working, so mask out the absorbance region and just fit the sidebands
  curr_comb_fit.fit_funct.FixParameter(3,0);
  TFitResultPtr datafitresptr = curr_comb_fit.PerformFitOnData(RemoveRegion(current_dark_sub, abs_region_low, abs_region_high),false);
  
  /* -- nah this still doesn't work
  // fix the pure parts
  std::vector<int> parstofix{0,1,2,4,5,6};
  for(int& pari : parstofix){
    curr_comb_fit.fit_funct.FixParameter(pari,curr_comb_fit.fit_funct.GetParameter(pari));
  }
  // release the absorbance and try to fit just that component
  curr_comb_fit.fit_funct.ReleaseParameter(3);
  curr_comb_fit.fit_funct.SetParLimits(3,0,1.1);
  TFitResultPtr datafitresptr = curr_comb_fit.PerformFitOnData(current_dark_sub,true);
  */
  
  // debug
  purefitgraph = curr_comb_fit.GetGraph();
  // FIXME this fit is baaaad.
  SaveDebug(&purefitgraph, std::string{"datafit_"}+current_led+"_"+std::to_string(measurementnum));
  
  // update the datamodel
  purefitgraph = curr_comb_fit.GetGraphExcluding({current_led_info.combined_func->ABS_SCALING});
  m_data->CStore.Set("purefit",reinterpret_cast<intptr_t>(&purefitgraph));
  SaveDebug(&purefitgraph, std::string{"purefit_"}+current_led+"_"+std::to_string(measurementnum));
  
  // doesn't this need a log10???
  //current_ratio_absorbtion = PWLogRatio(curr_comb_fit.GetGraphExcluding({current_led_info.combined_func->ABS_SCALING}), current_dark_sub);  -- the resulting concentrations from this are way off... did the calibration curve use log10?
  current_ratio_absorbtion = PWRatio(curr_comb_fit.GetGraphExcluding({current_led_info.combined_func->ABS_SCALING}), current_dark_sub);
  SaveDebug(&current_ratio_absorbtion, std::string{"ratioabs_"}+current_led+"_"+std::to_string(measurementnum));
  
  // update the datamodel
  m_data->CStore.Set("absorbance",reinterpret_cast<intptr_t>(&current_ratio_absorbtion));
  datafitresp = TFitResultPtr((TFitResult*)datafitresptr->Clone());
  intptr_t datafitresi = reinterpret_cast<intptr_t>(&datafitresp);
  m_data->CStore.Set("datafitresptr",datafitresi);
  
  // FIXME figure out a suitable metric of whether these fits worked or not
  int datafit_success =  !datafitresptr->IsEmpty() &&
                         datafitresptr->IsValid() &&
                         datafitresptr->Status()==0 &&
                         ((datafitresptr->Chi2()/datafitresptr->Ndf()) < 10.) &&
                         !TMath::IsNaN(datafitresptr->GetParams()[0]) &&
                         (datafitresptr->GetErrors()[0] < 0.5);
  m_data->CStore.Set("datafit_success",datafit_success);
  
  // now perform abs fit on ratio absorbance
  Log(m_unique_name+" fitting absorption curve", v_debug,m_verbose);
  FunctionalFit curr_abs_fit = current_led_info.absorbtion_fit;
  TFitResultPtr absfitresptr = curr_abs_fit.PerformFitOnData(current_ratio_absorbtion,false);

  // update the datamodel
  GetAbsRegion(curr_abs_fit.GetGraph(), absfitgraph);
  SaveDebug(&absfitgraph, std::string{"absfit_"}+current_led+"_"+std::to_string(measurementnum));
  m_data->CStore.Set("absfit",reinterpret_cast<intptr_t>(&absfitgraph));
  
  absfitresp = TFitResultPtr((TFitResult*)absfitresptr->Clone());
  intptr_t absfitresi = reinterpret_cast<intptr_t>(&absfitresp);
  m_data->CStore.Set("absfitresptr",absfitresi);
  
  int absfit_success =  !absfitresptr->IsEmpty() &&
                        absfitresptr->IsValid() &&
                        absfitresptr->Status()==0 &&
                        ((absfitresptr->Chi2()/absfitresptr->Ndf()) < 10.) &&
                        !TMath::IsNaN(absfitresptr->GetParams()[0]) &&
                        (absfitresptr->GetErrors()[0] < 0.5);
  m_data->CStore.Set("absfit_success",absfit_success);
  
  //extract metric from fit result, then get the concentration prediction from the calibration curve 
  double metric = curr_abs_fit.GetParameterValue(current_led_info.abs_func->ABS_SCALING);
  double conc_prediction = current_led_info.calibration_curve_ptr->GetX(metric);
  
  // debug: detect steps
  /*
  std::string rawfilename;
  get_ok = m_data->CStore.Get("Filename",rawfilename);
  static double last_conc_prediction = conc_prediction;
  if(std::abs(conc_prediction - last_conc_prediction)>0.002){
    Log(m_unique_name+" Step change detected! Concentration goes from "+std::to_string(last_conc_prediction)
       +" to "+std::to_string(conc_prediction)+" in file "+rawfilename+" at "+GetCurrentTimestamp(),
       v_error,m_verbose);
  }
  last_conc_prediction = conc_prediction;
  */
  
  // calculate suitable errors
  double metric_err = absfitresptr->GetErrors()[0];
  std::pair<double,double> metric_and_err{metric, metric_err};
  m_data->CStore.Set("metric_and_err",metric_and_err);
  
  TF1* calib_curve_func = led_info_map.at(current_led).calibration_curve_ptr;
  double gd_conc_err = metric_err * calib_curve_func->Derivative(metric);
  std::pair<double,double> conc_and_err{conc_prediction, gd_conc_err};
  m_data->CStore.Set("conc_and_err",conc_and_err);

  // Inform downstream tools that a new measurement is available
  // maybe we could use the value to indicate if the data is good?
  m_data->CStore.Set("NewMatthewAnalyse",current_led);
  Log(m_unique_name+" finished", v_debug,m_verbose);
  
  return true;
}

std::string MatthewAnalysisStrikesBack::GetCurrentTimestamp(){
  std::string timestamp;
  get_ok = m_data->CStore.Get("dbtimestamp",timestamp);
  if(!get_ok){
    timestamp=to_simple_string(m_data->measurment_time);
  }
  return timestamp;
}

bool MatthewAnalysisStrikesBack::Finalise(){
  if(fdebug){
    fdebug->Close();
    delete fdebug;
    fdebug=nullptr;
  }
  return true;
}

bool MatthewAnalysisStrikesBack::ReinitDataModel(){
  
  m_data->CStore.Remove("dark_subtracted_data_in");
  m_data->CStore.Remove("dark_subtracted_data_out");
  m_data->CStore.Remove("purefit");
  m_data->CStore.Remove("absorbance");
  m_data->CStore.Remove("absfit");
  m_data->CStore.Remove("dark_mean");
  m_data->CStore.Remove("dark_sigma");
  m_data->CStore.Remove("raw_absregion_min");
  m_data->CStore.Remove("raw_absregion_max");
  m_data->CStore.Remove("datafitresptr");
  m_data->CStore.Remove("absfitresptr");
  m_data->CStore.Remove("datafit_success");
  m_data->CStore.Remove("absfit_success");
  m_data->CStore.Remove("metric_and_err");
  m_data->CStore.Remove("conc_and_err");
  
  return true;
  
}

void MatthewAnalysisStrikesBack::GetWavelengthIndices(const TGraph& g){
  
  if(abs_region_npoints==0){
    
    double* wlarray = g.GetX();
    if(wlarray==nullptr){
      throw std::runtime_error("GetWavelengthIndices: null wavelengths array!!!\n");
    }
    
    for(int i=0; i<g.GetN(); ++i){
      if(wlarray[i] < sideband_region_low) ++sideband_region_low_index;
      if(wlarray[i] < abs_region_low) ++abs_region_low_index;
      if(wlarray[i] < abs_region_high) ++abs_region_high_index;
      if(wlarray[i] < sideband_region_high) ++sideband_region_high_index;
      else break;
    }
    abs_region_low = wlarray[abs_region_low_index];
    abs_region_high = wlarray[abs_region_high_index];
    abs_region_npoints = abs_region_high_index - abs_region_low_index;
    sideband_region_npoints = sideband_region_high_index - sideband_region_low_index;
  }
  
  return;
}

void MatthewAnalysisStrikesBack::GetAbsRegion(const TGraph& datagraph, TGraph& absgraph){
  
  double* x_in_vals = datagraph.GetX();
  double* y_in_vals = datagraph.GetY();
  if(x_in_vals==nullptr || y_in_vals==nullptr){
    throw std::invalid_argument("GetAbsRegion: x or y array nullptr!!!\n");
  }
  if(datagraph.GetN()==0){
    throw std::invalid_argument("GetAbsRegion: no datapoints!!!\n");
  }
  
  GetWavelengthIndices(datagraph);
  absgraph.Set(abs_region_npoints);
  double* x_out_vals = absgraph.GetX();
  double* y_out_vals = absgraph.GetY();
  
  if(datagraph.GetN()==number_of_points){
    // get the indices of the absorbance region
    for(int i=0, j=abs_region_low_index; i<abs_region_npoints; ++i, ++j){
      x_out_vals[i] = x_in_vals[j];
      y_out_vals[i] = y_in_vals[j];
      ++i;
    }
  } else {
    // tgraph is already trimmed
    for(int i=0, j=0; j<datagraph.GetN(); ++j){
      if(x_in_vals[j]>abs_region_high) break;
      if(x_in_vals[j]<abs_region_low) continue;
      x_out_vals[i] = x_in_vals[j];
      y_out_vals[i] = y_in_vals[j];
      ++i;
    }
  }
  
  return;
}

void MatthewAnalysisStrikesBack::GetSidebandRegion(const TGraph& datagraph, TGraph& absgraph){
  
  double* x_in_vals = datagraph.GetX();
  double* y_in_vals = datagraph.GetY();
  if(x_in_vals==nullptr || y_in_vals==nullptr){
    throw std::invalid_argument("GetSidebandRegion: x or y array nullptr!!!\n");
  }
  if(datagraph.GetN()==0){
    throw std::invalid_argument("GetAbsRegion: no datapoints!!!\n");
  }
  
  absgraph.Set(sideband_region_npoints - abs_region_npoints);
  double* x_out_vals = absgraph.GetX();
  double* y_out_vals = absgraph.GetY();
  
  if(datagraph.GetN()==number_of_points){
    // get the indices of the sideband region
    GetWavelengthIndices(datagraph);
    for(int i=0, j=abs_region_low_index; i<abs_region_npoints; ++i, ++j){
      x_out_vals[i] = x_in_vals[j];
      y_out_vals[i] = y_in_vals[j];
      ++i;
    }
  } else {
    // tgraph is already trimmed
    for(int i=0, j=0; j<datagraph.GetN(); ++j){
      if(x_in_vals[j]<sideband_region_high) break;
      if(x_in_vals[j]<sideband_region_low) continue;
      if((x_in_vals[j]>abs_region_low) && (x_in_vals[j]<abs_region_high)) continue;
      x_out_vals[i] = x_in_vals[j];
      y_out_vals[i] = y_in_vals[j];
      ++i;
    }
  }
  
  return;
}

void MatthewAnalysisStrikesBack::GetDarkTraceParams(double& mean, double& sigma){
  if (dark_tree_ptr ==  nullptr){
    throw std::invalid_argument("GetDarkTraceParams: dark tree is nullptr!!!\n");
  }
  
  std::vector<double> dark_values;
  std::vector<double>* dark_values_p = &dark_values;
  
  const bool set_branch_success = (dark_tree_ptr->SetBranchAddress("value", &dark_values_p) == 0);
  if (!set_branch_success){
    throw std::runtime_error("GetDarkTraceParams: error when setting branch address!!!\n");
  }
  
  const bool get_entry_success = dark_tree_ptr->GetEntry(dark_tree_ptr->GetEntries()-1);
  if (!get_entry_success){
    throw std::runtime_error("GetDarkTraceParams: failed to get entry!!!\n");
  }
  
  if (static_cast<int>(dark_values.size()) != number_of_points){
    throw std::runtime_error("GetDarkTraceParams: retrieved vector is incorrect size!!!\n");
  }
  
  mean = (std::accumulate(dark_values.begin(), dark_values.end(), 0))/number_of_points;
  sigma = TMath::StdDev(number_of_points, dark_values.data());
  
  return;
  
}

void MatthewAnalysisStrikesBack::GetRawMinMax(double& raw_min, double& raw_max){
  
  if (led_tree_ptr ==  nullptr){
    throw std::invalid_argument("GetRawMinMax: led tree is nullptr!!!\n");
  }
  
  std::vector<double> led_values, wavelengths;
  std::vector<double>* led_values_p = &led_values;
  std::vector<double>* wavelengths_p = &wavelengths;
  
  const bool set_branch_success = 
    (led_tree_ptr->SetBranchAddress("value", &led_values_p) == 0) &&
    (led_tree_ptr->SetBranchAddress("wavelength", &wavelengths_p) == 0);
  if (!set_branch_success){
    throw std::runtime_error("GetRawMinMax: error when setting branch addresses!!!\n");
  }
  
  const bool get_entry_success = led_tree_ptr->GetEntry(led_tree_ptr->GetEntries()-1);
  if (!get_entry_success){
    throw std::runtime_error("GetRawMinMax: failed to get entry!!!\n");
  }
  
  if ( static_cast<int>(led_values.size()) != number_of_points ||
       static_cast<int>(wavelengths.size()) != number_of_points ){
    throw std::runtime_error("GetRawMinMax: retrieved vectors of incorrect size!!!\n");
  }
  
  std::vector<double>::iterator abs_it_low = led_values.begin();
  std::vector<double>::iterator abs_it_high = led_values.begin();
  for(auto&& awl : wavelengths){
    if(awl < abs_region_low) ++abs_it_low;
    if(awl < abs_region_high) ++abs_it_high;
    else break;
  }
  
  raw_min = *std::min_element(abs_it_low, abs_it_high);
  raw_max = *std::max_element(abs_it_low, abs_it_high);
  
  return;
}



TGraph MatthewAnalysisStrikesBack::GetPure(const std::string& led_name) {
  // retrieve pure trace from filename specified in the config file, or else get from the database
  std::string pureID;  // for noting what was used in the database
  TGraph result;
  std::string pure_fname = "";
  int pure_offset = -1;
  bool ok = m_variables.Get("pure_fname", pure_fname) &&
    m_variables.Get("pure_offset", pure_offset);
  
  if (ok && pure_fname != "" && pure_offset != -1){
    Log(m_unique_name+" getting pure trace for led "+led_name+" from "+pure_fname
        +" entry "+std::to_string(pure_offset), v_debug,m_verbose);
    result = GetDarkSubtractFromFile(pure_fname, led_name, pure_offset);
    pureID = pure_fname + "::" + std::to_string(pure_offset);
  }
  else {
    int  pureref_ver = -1;
    ok = m_variables.Get("pureref_ver", pureref_ver);
    if (ok && pureref_ver != -1){
      Log(m_unique_name+" getting pure trace for led "+led_name
          +" from database entry "+std::to_string(pureref_ver), v_debug,m_verbose);
      result = GetPureFromDB(pureref_ver, led_name);
      pureID = std::to_string(pureref_ver);
    }
    else {
      throw std::runtime_error("MatthewAnalysisStrikesBack::GetPure - no pure specified!");
    }
  }
  
  std::string key = "purerefID_"+led_name;
  m_data->CStore.Set(key, pureID);
  
  return result;
}

TGraph MatthewAnalysisStrikesBack::GetHighConc(const std::string& led_name){
  // retrieve high conc trace from filename specified in the config file, or else get from the database 
  std::string highConcID;  // for noting what was used in the database
  TGraph result;
  std::string highconc_fname = "";
  int highconc_offset = -1;
  bool ok = m_variables.Get("highconc_fname", highconc_fname) &&
    m_variables.Get("highconc_offset", highconc_offset);
  
  if (ok && highconc_fname != "" && highconc_offset != -1){
    Log(m_unique_name+" getting high conc reference trace for led "+led_name+" from "+highconc_fname
        +" entry "+std::to_string(highconc_offset), v_debug,m_verbose);
    result = GetDarkSubtractFromFile(highconc_fname, led_name, highconc_offset);
    highConcID = highconc_fname + "::" + std::to_string(highconc_offset);
  }
  else {
    int  highconcref_ver = -1;
    ok = m_variables.Get("highconcref_ver", highconcref_ver);
    if (ok && highconcref_ver != -1){
      Log(m_unique_name+" getting high conc reference trace for led "+led_name
          +" from database entry "+std::to_string(highconcref_ver), v_debug,m_verbose);
      result = GetHighConcFromDB(highconcref_ver, led_name);
      highConcID = std::to_string(highconcref_ver);
    }
    else {
      throw std::runtime_error("MatthewAnalysisStrikesBack::GetHighconc - no highconc specified!");
    }
  }
  
  std::string key = "highconcrefID_"+led_name;
  m_data->CStore.Set(key, highConcID);
  
  return result;
}


TGraph MatthewAnalysisStrikesBack::GetPureFromDB(const int& pureref_ver, const std::string& ledToAnalyse) const {
  // build pure trace from values extracted from the database
  std::string query_string = "SELECT values->'yvals' FROM data WHERE name='pure_curve'"
    " AND values->'yvals' IS NOT NULL AND ledname IS NOT NULL"
    " AND values->'version' IS NOT NULL"
    " AND ledname="+m_data->postgres.pqxx_quote(ledToAnalyse)+
    " AND values->'version'="+ m_data->postgres.pqxx_quote(pureref_ver);
  std::string pureref_json="";
  bool get_ok = m_data->postgres.ExecuteQuery(query_string, pureref_json);
  if(!get_ok || pureref_json==""){
    throw std::runtime_error("MatthewAnalysisStrikesBack::GetPureFromDB: Failed to obtain y values for pure reference for led: " +
                             ledToAnalyse+" and pure version: " + std::to_string(pureref_ver));
  }
  // the values string is a json array; i.e. '[val1, val2, val3...]'
  // strip the '[' and ']'
  pureref_json = pureref_json.substr(1,pureref_json.length()-2);
  // parse it
  std::stringstream ss(pureref_json);
  std::string tmp;
  std::vector<double> pureref_yvals;
  while(std::getline(ss,tmp,',')){
    char* endptr = &tmp[0];
    double nextval = strtod(tmp.c_str(),&endptr);
    if(endptr==&tmp[0]){
    throw std::runtime_error("MatthewAnalysisStrikesBack::GetPureFromDB: Failed to parse y values for pure reference for led: " +
                             ledToAnalyse+" and pure version: " + std::to_string(pureref_ver));
    }
    pureref_yvals.push_back(nextval);
  }
  if(pureref_yvals.size()==0){
    throw std::runtime_error("MatthewAnalysisStrikesBack::GetPureFromDB: Parse zero  y values for pure reference for led: " +
                             ledToAnalyse+" and pure version: " + std::to_string(pureref_ver));
  }
  
  // repeat for the x-values
  query_string = "SELECT values->'xvals' FROM data WHERE name='pure_curve'"
    " AND values->'xvals' IS NOT NULL AND ledname IS NOT NULL"
    " AND values->'version' IS NOT NULL"
    " AND ledname="+m_data->postgres.pqxx_quote(ledToAnalyse)+
    " AND values->'version'="+ m_data->postgres.pqxx_quote(pureref_ver);
  pureref_json="";
  get_ok = m_data->postgres.ExecuteQuery(query_string, pureref_json);
  if(!get_ok || pureref_json==""){
    throw std::runtime_error("MatthewAnalysisStrikesBack::GetPureFromDB: Failed to obtain x values for pure reference for led: " +
                             ledToAnalyse+" and pure version: " + std::to_string(pureref_ver));
  }
  // the values string is a json array; i.e. '[val1, val2, val3...]'
  // strip the '[' and ']'
  pureref_json = pureref_json.substr(1,pureref_json.length()-2);
  // parse it
  ss.clear();
  ss.str(pureref_json);
  std::vector<double> pureref_xvals;
  while(std::getline(ss,tmp,',')){
    char* endptr = &tmp[0];
    double nextval = strtod(tmp.c_str(),&endptr);
    if(endptr==&tmp[0]){
      throw std::runtime_error("MatthewAnalysisStrikesBack::GetPureFromDB: Failed to parse x values for pure reference for led: " +
                               ledToAnalyse+" and pure version: " + std::to_string(pureref_ver));
    }
    pureref_xvals.push_back(nextval);
  }
  if(pureref_xvals.size()==0){
    throw std::runtime_error("MatthewAnalysisStrikesBack::GetPureFromDB: Parsed zero x values for pure reference for led: " +
                             ledToAnalyse+" and pure version: " + std::to_string(pureref_ver));
  }
  
  TGraph dark_subtracted_pure = TGraph(pureref_xvals.size(), pureref_xvals.data(), pureref_yvals.data());
  std::string purename="g_pureref_"+ledToAnalyse;
  dark_subtracted_pure.SetName(purename.c_str());
  dark_subtracted_pure.SetTitle(purename.c_str());
  return dark_subtracted_pure;
}

TGraph MatthewAnalysisStrikesBack::GetHighConcFromDB(const int& highconcref_ver, const std::string& ledToAnalyse) const {
  // build high conc from values extracted from database 
  std::string query_string = "SELECT values->'yvals' FROM data WHERE name='highconc_curve'"
    " AND values->'yvals' IS NOT NULL AND ledname IS NOT NULL"
    " AND values->'version' IS NOT NULL"
    " AND ledname="+m_data->postgres.pqxx_quote(ledToAnalyse)+
    " AND values->'version'="+ m_data->postgres.pqxx_quote(highconcref_ver);
  std::string highconcref_json="";
  bool get_ok = m_data->postgres.ExecuteQuery(query_string, highconcref_json);
  if(!get_ok || highconcref_json==""){
    throw std::runtime_error("MatthewAnalysisStrikesBack::GetHighconcFromDB: Failed to obtain y values for highconc reference for led: " +
                             ledToAnalyse+" and highconc version: " + std::to_string(highconcref_ver));
  }
  // the values string is a json array; i.e. '[val1, val2, val3...]'
  // strip the '[' and ']'
  highconcref_json = highconcref_json.substr(1,highconcref_json.length()-2);
  // parse it
  std::stringstream ss(highconcref_json);
  std::string tmp;
  std::vector<double> highconcref_yvals;
  while(std::getline(ss,tmp,',')){
    char* endptr = &tmp[0];
    double nextval = strtod(tmp.c_str(),&endptr);
    if(endptr==&tmp[0]){
    throw std::runtime_error("MatthewAnalysisStrikesBack::GetHighconcFromDB: Failed to parse y values for highconc reference for led: " +
                             ledToAnalyse+" and highconc version: " + std::to_string(highconcref_ver));
    }
    highconcref_yvals.push_back(nextval);
  }
  if(highconcref_yvals.size()==0){
    throw std::runtime_error("MatthewAnalysisStrikesBack::GetHighconcFromDB: Parse zero  y values for highconc reference for led: " +
                             ledToAnalyse+" and highconc version: " + std::to_string(highconcref_ver));
  }
  
  // repeat for the x-values
  query_string = "SELECT values->'xvals' FROM data WHERE name='highconc_curve'"
    " AND values->'xvals' IS NOT NULL AND ledname IS NOT NULL"
    " AND values->'version' IS NOT NULL"
    " AND ledname="+m_data->postgres.pqxx_quote(ledToAnalyse)+
    " AND values->'version'="+ m_data->postgres.pqxx_quote(highconcref_ver);
  highconcref_json="";
  get_ok = m_data->postgres.ExecuteQuery(query_string, highconcref_json);
  if(!get_ok || highconcref_json==""){
    throw std::runtime_error("MatthewAnalysisStrikesBack::GetHighConcFromDB: Failed to obtain x values for highconc reference for led: " +
                             ledToAnalyse+" and highconc version: " + std::to_string(highconcref_ver));
  }
  // the values string is a json array; i.e. '[val1, val2, val3...]'
  // strip the '[' and ']'
  highconcref_json = highconcref_json.substr(1,highconcref_json.length()-2);
  // parse it
  ss.clear();
  ss.str(highconcref_json);
  std::vector<double> highconcref_xvals;
  while(std::getline(ss,tmp,',')){
    char* endptr = &tmp[0];
    double nextval = strtod(tmp.c_str(),&endptr);
    if(endptr==&tmp[0]){
      throw std::runtime_error("MatthewAnalysisStrikesBack::GetHighconcFromDB: Failed to parse x values for highconc reference for led: " +
                               ledToAnalyse+" and highconc version: " + std::to_string(highconcref_ver));
    }
    highconcref_xvals.push_back(nextval);
  }
  if(highconcref_xvals.size()==0){
    throw std::runtime_error("MatthewAnalysisStrikesBack::GetHighconcFromDB: Parsed zero x values for highconc reference for led: " +
                             ledToAnalyse+" and highconc version: " + std::to_string(highconcref_ver));
  }
  
  TGraph dark_subtracted_highconc = TGraph(highconcref_xvals.size(), highconcref_xvals.data(), highconcref_yvals.data());
  std::string highconcname="g_highconcref_"+ledToAnalyse;
  dark_subtracted_highconc.SetName(highconcname.c_str());
  dark_subtracted_highconc.SetTitle(highconcname.c_str());
  return dark_subtracted_highconc;
}

bool MatthewAnalysisStrikesBack::ReadyToAnalyse() const {
  // checks the value of the analysis CStore variable 
  bool ready = false;
  std::string analyse="";
  m_data->CStore.Get("Analyse", analyse);
  if (analyse == "Analyse"){
    m_data->CStore.Remove("Analyse");
    ready = true;
  }
  return ready;
}

void MatthewAnalysisStrikesBack::GetCurrentLED(){
  // get name of current led - throw exception if not found, or cannot be analysed
  bool ok = m_data->CStore.Get("ledToAnalyse", current_led);
  if (!ok || current_led.empty()){
    throw std::invalid_argument("MatthewAnalysisStrikesBack::Execute - LED to analyse is blank!");
  }
  if ( led_info_map.count(current_led) != 1){
    throw std::invalid_argument("MatthewAnalysisStrikesBack::Execute - Cannot Analyse LED: "+current_led+" !");
  }
  return;
}

void MatthewAnalysisStrikesBack::GetDarkAndLEDTrees(){
  // from the current traces taken, sort the led and dark tree and store as class members
  for (const auto& [name, tree_ptr] : m_data->m_trees){
    if (name == current_led){led_tree_ptr = tree_ptr;}
    else if (boost::iequals(name, "dark")){dark_tree_ptr = tree_ptr;}
  }
  if(led_tree_ptr == nullptr){
    throw std::runtime_error("MatthewAnalysisStrikesBack::GetDarkAndLEDTrees: Failed to find led trace!");
  }
  if(dark_tree_ptr == nullptr){
    throw std::runtime_error("MatthewAnalysisStrikesBack::GetDarkAndLEDTrees: Failed to find dark trace!");
  }
  return;
}

TF1* MatthewAnalysisStrikesBack::GetCalibrationCurve(const std::string& led_name){
  // retrieve calibration curve from file whose name is specified on config file
  std::string calibcurveID;  // for storing in database
  TF1* result = nullptr;
  const std::string calib_f_str = "calib_curve_";
  std::string calib_fname = "";
  const std::string calib_v_str = "cal_tf1_name_";
  std::string calib_vname = "";
  bool ok = m_variables.Get(calib_f_str+led_name, calib_fname) &&
    m_variables.Get(calib_v_str+led_name, calib_vname);

  if (ok && !calib_fname.empty() && !calib_vname.empty()){
    Log(m_unique_name+" getting calibration curve for led "+led_name
      +" from file "+calib_fname+", key '"+calib_vname+"'", v_debug,m_verbose);
    TFile* f = TFile::Open(calib_fname.c_str());
    if (f->IsZombie()){
      throw std::runtime_error("MatthewAnalysisStrikesBack::GetCalibrationCurve: Couldn't open calibration curve file!");
    }
    result = static_cast<TF1*>(f->Get(calib_vname.c_str()));
    if (result == nullptr){
      throw std::runtime_error("MatthewAnalysisStrikesBack::GetCalibrationCurve: Couldn't retrieve calibration curve!");
    }
    calibcurveID = calib_fname+"::"+calib_vname;
    // it is very hard to get a clear picture (the relevant DOxygen page does not say)
    // but based on https://roottalk.root.cern.narkive.com/bjocjY0o/root-readobj-memory-leak
    // it seems we now own this TF1. We can (and should) close the file, and delete it later.
    f->Close();
  }
  else {
    // get from db
    int  calibcurve_ver = -1;
    ok = m_variables.Get("calibcurve_ver", calibcurve_ver);
    if (ok && calibcurve_ver != -1){
      Log(m_unique_name+" getting calibration curve for led "+led_name
          +" from database entry "+std::to_string(calibcurve_ver), v_debug,m_verbose);
      result = GetCalibCurveFromDB(calibcurve_ver, led_name);
      calibcurveID = std::to_string(calibcurve_ver);
    }
    else {
      throw std::runtime_error("MatthewAnalysisStrikesBack::GetCalibrationCurve - no calibcurve specified!");
    }
  }
  
  std::string key = "calibcurveID_"+led_name;
  m_data->CStore.Set(key, calibcurveID);
  
  return result;
}

TF1* MatthewAnalysisStrikesBack::GetCalibCurveFromDB(const int& calibID, const std::string& ledToAnalyse){
  // Constructs calibration curve from parameters in the database
  
  // new calibration curve can be inserted with e.g.:
  // psql -U postgres -d "rundb" -c "INSERT INTO data (timestamp, tool, ledname, name, values) VALUES ('now()', 'MatthewAnalysisStrikesBack', '275_A', 'calibration_curve', '{\"version\":0, \"formula\":\"pol6\", \"params\":[0.0079581671, 2.7415760, -31.591756, 478.69924, 18682.891, -29982.759] }' );"
  
  // use version number to lookup formula from database
  std::string formula_str="";
  std::string query_string = "SELECT values->>'formula' FROM data WHERE name='calibration_curve'"
                             " AND tool='MatthewAnalysisStrikesBack'"
                             " AND ledname="+m_data->postgres.pqxx_quote(ledToAnalyse)+
                             " AND values->'formula' IS NOT NULL"
                             " AND values->'version' IS NOT NULL"
                             " AND values->'version'="+m_data->postgres.pqxx_quote(calibID);
  get_ok = m_data->postgres.ExecuteQuery(query_string, formula_str);
  if(!get_ok || formula_str==""){
    throw std::runtime_error("MatthewAnalysisStrikesBack::GetCalibCurveFromDB - failed to find DB entry!");
  }
  // we must strip enclosing quotations or the TF1 constructor goes berzerk
  //formula_str = formula_str.substr(1,formula_str.length()-2);
  
  // another query for the curve parameters
  // when querying attributes from JSON fields, we need to explicitly ensure the checked
  // attributes exist, or exclude the entry from the search, or the query will fail.
  query_string = " SELECT values->'params' FROM data WHERE name='calibration_curve'"
                 " AND tool='MatthewAnalysisStrikesBack'"
                 " AND ledname="+m_data->postgres.pqxx_quote(ledToAnalyse)+
                 " AND values->'params' IS NOT NULL"
                 " AND values->'version' IS NOT NULL"
                 " AND values->'version'=" + m_data->postgres.pqxx_quote(calibID);
  std::string params_str;
  get_ok &= m_data->postgres.ExecuteQuery(query_string, params_str);
  
  // check for errors
  if(!get_ok){
    throw std::runtime_error("MatthewAnalysisStrikesBack::GetCalibCurveFromDB - failed to get parameters!");
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
        throw std::runtime_error("MatthewAnalysisStrikesBack::GetCalibCurveFromDB - failed to parse params!");
    }
    calib_coefficients.push_back(nextval);
  }
  if(calib_coefficients.size()==0){
    throw std::runtime_error("MatthewAnalysisStrikesBack::GetCalibCurveFromDB - no params!");
  }
  
  // construct a TF1 in the map from the formula and parameters given
  std::string calib_name = "f_cal_"+ledToAnalyse;
  TF1* cal_curve = new TF1(calib_name.c_str(), formula_str.c_str(), 0, 0.4);
  cal_curve->SetParameters(calib_coefficients.data());
  
  // check the TF1 constructed successfully
  if(!cal_curve->IsValid()){
    throw std::runtime_error("MatthewAnalysisStrikesBack::GetCalibCurveFromDB - TF1 invalid!");
  }
  
  return cal_curve;
}
