#include "MatthewAnalysisStrikesBack.h"

#include "gad_utils.h"
#include <memory>

MatthewAnalysisStrikesBack::MatthewAnalysisStrikesBack():Tool(){}

bool MatthewAnalysisStrikesBack::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data = &data;

  for (const auto& led_name : {"275_A", "275_B"}){
    Log(std::string("MatthewAnalysisStrikesBack:Initialise - Getting info for LED: ")+led_name, 0,0);
      
    LEDInfo led_info;

    //retrieve reference pure and non-pure 
    const TGraph pure_ds = GetPure(led_name);
    const TGraph high_conc_ds = GetHighConc(led_name);

    // get calibration curve
    led_info.calibration_curve_ptr = GetCalibrationCurve(led_name);
    
    //create FunctionalFit object to do initial "simple", ie remove Gd region and fit sidebands
    std::shared_ptr<SimplePureFunc> simple_pure_func = std::make_shared<SimplePureFunc>(pure_ds);
    FunctionalFit simple_fit = FunctionalFit(simple_pure_func.get(), "simple"); 
    simple_fit.SetExampleGraph(pure_ds);

    //fit non-pure with the simple fit and extract gd attenuation and ratio absorbance shapes
    simple_fit.PerformFitOnData(high_conc_ds);
    const TGraph simple_fit_result = simple_fit.GetGraph();
    const TGraph gd_abs_attenuation = CalculateGdAbs(simple_fit_result, high_conc_ds);
    const TGraph ratio_absorbance = PWRatio(simple_fit_result, high_conc_ds);

    //create FunctionalFit object to do absorbance fit on data
    led_info.abs_func = std::make_shared<AbsFunc>(ratio_absorbance);
    led_info.absorbtion_fit = FunctionalFit(led_info.abs_func.get(), "abs");
    led_info.absorbtion_fit.SetExampleGraph(pure_ds);

    //create FunctionalFit object to do inital pure removal fit
    led_info.combined_func = std::make_unique<CombinedGdPureFunc_DATA>(pure_ds, gd_abs_attenuation);
    led_info.combined_fit = FunctionalFit(led_info.combined_func.get(), "combined");
    led_info.combined_fit.SetExampleGraph(pure_ds);

    //add the above fitting objects to a map of leds
    led_info_map[led_name] = led_info;
  }
    
  return true;
}

bool MatthewAnalysisStrikesBack::Execute(){

  // waits for analysis flag
  if (!ReadyToAnalyse()){
    return true;
  }

  // get the name of the name of the current led
  GetCurrentLED();
  // from the currently taken traces, sort out the dark and led trees
  GetDarkAndLEDTrees();

  // calculate dark subtract from these traces
  const TGraph current_dark_sub = DarkSubtractFromTreePtrs(led_tree_ptr, dark_tree_ptr);

  // retrieve the functional fit for thios led
  LEDInfo current_led_info = led_info_map.at(current_led);
  FunctionalFit curr_comb_fit = current_led_info.combined_fit;

  // perform the fit on this new data
  curr_comb_fit.PerformFitOnData(current_dark_sub);
  const TGraph current_ratio_absorbtion = TrimGraph(PWRatio(curr_comb_fit.GetGraphExcluding({current_led_info.combined_func->ABS_SCALING}), current_dark_sub));

  // now perform abs fit on ratio absorbance
  FunctionalFit curr_abs_fit = current_led_info.absorbtion_fit;
  curr_abs_fit.PerformFitOnData(current_ratio_absorbtion);

  //extract metric from fit result, then get the concentration prediction from the calibration curve 
  double metric = curr_abs_fit.GetParameterValue(current_led_info.abs_func->ABS_SCALING);
  double conc_prediction = current_led_info.calibration_curve_ptr->GetX(metric);
  // where should I put this?
    
  return true;
}


bool MatthewAnalysisStrikesBack::Finalise(){

  return true;
}

TGraph MatthewAnalysisStrikesBack::GetPure(const std::string& led_name) {
  // retrieve pure trace from filename specified in the config file, or else get from the database 
  TGraph result;
  std::string pure_fname = "";
  int pure_offset = -1;
  bool ok = m_variables.Get("pure_fname", pure_fname) &&
    m_variables.Get("pure_offset", pure_offset);
  
  if (ok && pure_fname != "" && pure_offset != -1){
    result = GetDarkSubtractFromFile(pure_fname, led_name, pure_offset);
  }
  else {
    int  pureref_ver = -1;
    ok = m_variables.Get("pureref_ver", pureref_ver);
    if (ok && pureref_ver != -1){
      result = GetPureFromDB(pureref_ver, led_name);
    }
    else {
      throw std::runtime_error("MatthewAnalysisStrikesBack::GetPure - no pure specified!");
    }
  }
  return result;
}

TGraph MatthewAnalysisStrikesBack::GetHighConc(const std::string& led_name){
  // retrieve high conc trace from filename specified in the config file, or else get from the database 
  TGraph result;
  std::string highconc_fname = "";
  int highconc_offset = -1;
  bool ok = m_variables.Get("highconc_fname", highconc_fname) &&
    m_variables.Get("highconc_offset", highconc_offset);
  
  if (ok && highconc_fname != "" && highconc_offset != -1){
    result = GetDarkSubtractFromFile(highconc_fname, led_name, highconc_offset);
  }
  else {
    int  highconcref_ver = -1;
    ok = m_variables.Get("highconcref_ver", highconcref_ver);
    if (ok && highconcref_ver != -1){
      result = GetHighConcFromDB(highconcref_ver, led_name);
    }
    else {
      throw std::runtime_error("MatthewAnalysisStrikesBack::GetHighconc - no highconc specified!");
    }
  }
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
  m_data->CStore.Get("analyse_new", analyse);
  if (analyse == "analyse_new"){
    m_data->CStore.Remove("analyse_new");
    ready = true;
  }
  return ready;
}

void MatthewAnalysisStrikesBack::GetCurrentLED(){
  // get name of current led - throw exception if not found, or cannot be analysed
    bool ok = m_variables.Get("led", current_led);
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
  TF1* result = nullptr;
  const std::string calib_f_str = "calib_curve_";
  std::string calib_fname = "";
  const std::string calib_v_str = "cal_tf1_name_";
  std::string calib_vname = "";
  bool ok = m_variables.Get(calib_f_str+led_name, calib_fname) &&
    m_variables.Get(calib_v_str+led_name, calib_vname);

  if (ok && !calib_fname.empty() && !calib_vname.empty()){
    TFile* f = TFile::Open(calib_fname.c_str());
    if (f->IsZombie()){
      throw std::runtime_error("MatthewAnalysisStrikesBack::GetCalibrationCurve: Couldn't open calibration curve file!");
    }		     
    result = static_cast<TF1*>(f->Get(calib_vname.c_str()));
    if (result == nullptr){
      throw std::runtime_error("MatthewAnalysisStrikesBack::GetCalibrationCurve: Couldn't retrieve calibration curve!");
    }
    else {
      // get from db?
    }
  }
  return result;
}
