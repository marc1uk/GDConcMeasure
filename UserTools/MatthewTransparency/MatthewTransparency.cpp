#include "MatthewTransparency.h"
#include <sys/stat.h>

#include <memory>
#include <chrono>
#include <time.h>
#include <array>
#include <vector>
#include <string>
#include <sstream>

#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH2D.h"

MatthewTransparency::MatthewTransparency():Tool(){}

//struct MatthewTransparency::Transparency;

bool MatthewTransparency::Initialise(std::string configfile, DataModel &data){

  m_data= &data;
  
	// Retrieve configuration options from the postgres database
	int RunConfig=-1;
	m_data->vars.Get("RunConfig",RunConfig);
	Log(m_unique_name+" getting RunConfig from ToolChain",v_debug,verbosity);
	if(RunConfig>=0){
		std::string configtext;
		Log(m_unique_name+" getting tool configuration from database for runconfig "
		    +std::to_string(RunConfig),v_debug,verbosity);
		bool get_ok = m_data->postgres_helper.GetToolConfig(m_unique_name, configtext);
		if(!get_ok){
			Log(m_unique_name+" Failed to get Tool config from database!",v_error,verbosity);
			return false;
		}
		// parse the configuration to populate the m_variables Store.
		if(configtext!="") m_variables.Initialise(std::stringstream(configtext));
	}
  
  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  return true;
}

bool MatthewTransparency::Execute(){
  if (Ready()){
    m_data->CStore.Set("NewMatthewTransparency", true);

    //const std::string pure_file_name{"~/GDConcMeasure/data/2022/07/00117_12July22Calib_00000.root"}; // for debug
    
    int pure_ref_ver; //when DB support is ready
    m_variables.Get("pure_ref_ver", pure_ref_ver);
    m_data->CStore.Set("pure_refID_transparency", std::to_string(pure_ref_ver)); 
    
    std::array<double, N> dark_sub_sum{}, pure_dark_sub{}, pure_sum{}; 

    TTree *dark_tree_ptr = nullptr, *led_tree_ptr = nullptr;

    for(std::pair<const std::string, TTree*>& atree : m_data->m_trees){
      if (wavelengths.back() == 0){wavelengths = PopulateWavelength(atree.second);}
      if (atree.first == "Dark" || atree.first == "dark"){
          dark_tree_ptr = atree.second;
      }
    }

	for(std::pair<const std::string, TTree*>& atree : m_data->m_trees){
		if(atree.first != "Dark" && atree.first != "dark"){
			//pure_dark_sub = RetrievePureValuesFromFile(pure_file_name, atree.first); // for debug
			pure_dark_sub = RetrievePureValues(pure_ref_ver, atree.first);
			std::array<double, N> dark_sub = DarkSubtract(atree.second, dark_tree_ptr);
			for (auto i = 0; i < N; ++i){
				pure_sum.at(i) += pure_dark_sub.at(i);
				dark_sub_sum.at(i) += dark_sub.at(i);
			}
		}
	}
    
    const auto transparency_values = [=](){
				       std::array<double, N> result{};
				       for (auto i = 0; i < N; ++i){
					 //	 std::cout << i << std::endl;
					 // std::cout << dark_sub_sum.at(i) << std::endl;
					 // std::cout << pure_sum.at(i) << std::endl;
					 result.at(i) = dark_sub_sum.at(i) / pure_sum.at(i) ;
					 // std::cout << std::endl;
				       }
				       return result;
				     }();

    const Transparency t{transparency_values, GetDateTimeVec()};
    std::map<std::string, std::pair<double, double>> samples_map = CreateSamplesMap(t);
    m_data->CStore.Set("transparency_samples", samples_map);

    std::vector<double> temp_values(t.values.begin(), t.values.end());
    std::vector<double> temp_wav(wavelengths.begin(), wavelengths.end());
    m_data->CStore.Set("TransparencyTrace",temp_values);
    m_data->CStore.Set("TransparencyWls", temp_wav);
    
    SaveToMonthlyFile(t);
    AmmendSetOfTraces(t);
    MakeHeatMap(transparency_traces, wavelengths, "heat_map.root" );
    
    //    MakeQuickPlot(t.values, wavelengths, "trans_trans.root");
    // MakeQuickPlot(dark_sub_sum, wavelengths, "trans_darksub.root");
    //  MakeQuickPlot(pure_sum, wavelengths, "trans_pure.root");
    //    MakeDoublePlot(pure_sum, dark_sub_sum, wavelengths, "trans_double.root");
  } else {
  	// remove any previous flags, no new data
  	m_data->CStore.Remove("NewMatthewTransparency");
  }
      
  return true;
}

std::map<std::string, std::pair<double, double>> MatthewTransparency::CreateSamplesMap(const Transparency& t) const {
  std::map<std::string, std::pair<double, double>> result;
  const std::map<std::string, double> names_to_samples = {{"red", 620},
							  {"green", 500},
							  {"blue", 450}};

  for(const std::pair<const std::string, double>& asample : names_to_samples){
    const std::string name = asample.first;
    const double wav = asample.second;
    result.emplace(std::make_pair(name, std::make_pair(wav, GetValAtWavelength(t, wav))));
  }

  return result;
}

bool MatthewTransparency::Finalise(){
  //std::cout << "you did it Jimmy! YOU DID IT!" << std::endl;
  return true;
}

auto MatthewTransparency::PopulateWavelength(TTree* tree) const -> std::array<double, N> {
  // Takes a TTree* and returns an N sized vector with wavelength values, this is split from usual TTree reading since we only need to calculate wavelength vector once.
  
  std::vector<double>* w_ptr = nullptr;
  if(tree==nullptr){
  	std::cerr<<"PopulateWavelength with null tree"<<std::endl;
  	return std::array<double,N>{};
  }
  tree->SetBranchAddress("wavelength", &w_ptr);

  if (w_ptr == nullptr){
    std::cout << "Failed to set branch address to retrieve wavelength" << std::endl;
    return std::array<double,N>{};
  }
  else if (tree->GetEntriesFast() == 0){
    std::cout << "No entries in wavelength tree" << std::endl;
    return std::array<double,N>{};
  }

  const bool success = tree->GetEntry(0);
  const auto wav = [=](){
		     std::array<double, N> result;
		     for (auto i = 0; i < N; ++i){result.at(i) = w_ptr->at(i);}
		     return result;}();
  return wav;
}
 
auto MatthewTransparency::DarkSubtract(TTree* led, TTree* dark) const -> std::array<double, N> {
  // Takes two TTree* and returns an N sized vector containing dark subtracted intensity values of the LED that TTree* led corresponds to. 

  if (led == nullptr){
    std::cout << "led pointer is null in dark subtract" << std::endl;
  }

  if (dark == nullptr){
    std::cout << "dark pointer is null in dark subtract" << std::endl;
  }

  std::vector<double> *led_v = nullptr, *dark_v = nullptr;

  led->SetBranchAddress("value", &led_v);
  dark->SetBranchAddress("value", &dark_v);

  if (led_v == nullptr || dark_v == nullptr){
    std::cout << "Failed to set branch address in dark subtract calculation" << std::endl;
  }
  else if (led->GetEntriesFast() == 0 || dark->GetEntriesFast() == 0){
    std::cout << "No entries in dark or led tree" << std::endl;
  }
  led->GetEntry(0);
  dark->GetEntry(0);
  
  const auto dark_sub = [=](){
			  std::array<double, N> result;
			  for (auto i = 0; i < N; ++i){
			    result.at(i) = (led_v->at(i) - dark_v->at(i));
			  }
			  return result;
			}();
  
  return dark_sub;
}

auto MatthewTransparency::RetrievePureValues(const int pureref_ver, const std::string led) const -> std::array<double, N> {
  // Retrieve dark subtracted pure values from database for a given LED
  const std::string query_string{"SELECT values->'yvals' FROM data WHERE name='pure_curve'"
				 " AND values->'yvals' IS NOT NULL AND ledname IS NOT NULL"
				 " AND ledname="+m_data->postgres.pqxx_quote(led)+
				 " AND values->'version' IS NOT NULL"
				 " AND values->'version'="+ m_data->postgres.pqxx_quote(pureref_ver)};
  std::string pureref_json{""};

  m_data->postgres.ExecuteQuery(query_string, pureref_json);
  pureref_json = pureref_json.substr(1,pureref_json.length()-2);

  std::stringstream ss(pureref_json);
  std::string tmp;
  std::array<double, N> pureref_yvals{};
  auto i = 0;
  while(std::getline(ss,tmp,',')){
    char* endptr = &tmp[0];
    pureref_yvals.at(i) = strtod(tmp.c_str(), &endptr);
    ++i;
  }
  return pureref_yvals;  
}

auto MatthewTransparency::RetrievePureValuesFromFile(const std::string& fname, const std::string& led) const -> std::array<double, N> {
  // Until database integration is finished we'll just get the pure values from the first file in the calibration run. 
  std::unique_ptr<TFile> file(TFile::Open(fname.c_str()));
  TTree *led_tree = nullptr, *dark_tree = nullptr;
  if (file->IsZombie()){
    std::cout << "file not open" << std::endl;
  }
    
  led_tree = (TTree*) file->Get(led.c_str());

  if (led_tree == nullptr){
    std::cout << "couldn't get " << led << " tree"  << std::endl;
  }

  dark_tree = (TTree*) file->Get("Dark");
  
  if (dark_tree == nullptr){
    std::cout << "couldn't get dark tree" << std::endl;
  }
    
  const std::array<double, N> dark_sub = DarkSubtract(led_tree, dark_tree);
  
  return dark_sub;
}

void MatthewTransparency::SaveToMonthlyFile(const Transparency& t) {
  // Transparency values will be stored in a file that is newly generated once per month, 
  const std::array<double, N>* v_ptr = &(t.values);
  const std::array<int, 6>* t_ptr = &(t.date_time);

  const std::string transp_str{"transparency"}, val_str{"values"}, dt_str{"date_time"};
  const std::string fname{transp_str + "_" + t.date_time.at(0) + "_" + t.date_time.at(1) + ".root"};

  if (FileExists(fname)){
    std::unique_ptr<TFile> file(TFile::Open(fname.c_str(), "UPDATE"));
    if(!file){
        Log(m_unique_name+"::SaveToMonthlyFile unable to open file '"+fname+"'",v_error,verbosity);
        return;
    }
    TTree* tree = (TTree*) file->Get(transp_str.c_str());
    if(!tree){
        Log(m_unique_name+"::SaveToMonthlyFile failed to get tree '"+transp_str
           +"' from file '"+fname+"'",v_error,verbosity);
        return;
    }
    tree->SetBranchAddress(val_str.c_str(), &v_ptr);
    tree->SetBranchAddress(dt_str.c_str(), &t_ptr);
    tree->Fill();
    file->Write("*", TObject::kOverwrite);
  }
  else {
    std::unique_ptr<TFile> file(TFile::Open(fname.c_str(), "RECREATE"));
    if(!file){
        Log(m_unique_name+"::SaveToMonthlyFile failed to created new monthly file '"+fname+"'",v_error,verbosity);
        return;
    }
    std::unique_ptr<TTree> tree(new TTree(transp_str.c_str(), transp_str.c_str()));
    if(!tree){
        Log(m_unique_name+"::SaveToMonthlyFile failed to make new tree '"+transp_str
                   +"' for file '"+fname+"'",v_error,verbosity);
        return;
    }
    tree->Branch(val_str.c_str(), &v_ptr);
    tree->Branch(dt_str.c_str(), &t_ptr);
    tree->Fill();
    file->Write();
  }
  return;
}

bool MatthewTransparency::FileExists(std::string f) const {
  struct stat s;
  if (stat(f.c_str(), &s) == 0){
    return true;
  } else {
    return false;
  }
}
  
std::array<int, 6> MatthewTransparency::GetDateTimeVec() const{
  const auto now = std::chrono::system_clock::now();
  const std::time_t t = std::chrono::system_clock::to_time_t(now);
  tm local_tm = *localtime(&t);
  return {local_tm.tm_year+1900, local_tm.tm_mon+1, local_tm.tm_mday,
	  local_tm.tm_hour, local_tm.tm_min, local_tm.tm_sec};
}

bool MatthewTransparency::Ready() const {
  bool ready{false};
  std::string transparency{""};
  bool success = m_data->CStore.Get("Transparency", transparency);
  if (success && transparency == "Transparency"){
    m_data->CStore.Remove("Transparency");
    transparency.clear();
    ready = true;
  }
  return ready;
}

void MatthewTransparency::AmmendSetOfTraces(const Transparency& t) {
  
  constexpr int numb_measurements_in_month{2688};
  if (transparency_traces.size() ==  numb_measurements_in_month){
    transparency_traces.pop_back();
    transparency_times.pop_back();
  }
  transparency_traces.push_back(t.values);
  transparency_times.push_back(t.date_time);
}

double MatthewTransparency::GetValAtWavelength(const Transparency& t, const double wav) const {
  auto it  = std::find_if(wavelengths.begin(), wavelengths.end(), [wav](double v){return v >= wav;});
  return t.values.at(it - wavelengths.begin());
}

void MatthewTransparency::MakeQuickPlot(const std::array<double, N>& v, const std::array<double, N>& w, const std::string& f) const {
  TGraph g(v.size(), w.data(), v.data());
  TCanvas c1("c1", "c1", 2000, 2000);
  g.Draw();
  c1.SaveAs(("/home/pi/GDConcMeasure/" + f).c_str());
  
}

void MatthewTransparency::MakeDoublePlot(const std::array<double, N>& v1, const std::array<double, N>& v2, const std::array<double, N>& w, const std::string& f) const {
  TGraph g1(v1.size(), w.data(), v1.data());
  TGraph g2(v2.size(), w.data(), v2.data());
  TCanvas c1("c1", "c1", 2000, 2000);
  g1.Draw();
  g2.Draw("SAME");
  c1.SaveAs(("/home/pi/GDConcMeasure/" + f).c_str());
  
}

void MatthewTransparency::MakeHeatMap(const std::vector< std::array<double,N> >& traces, const std::array<double,N>& wav, const std::string& f) const {

  //Just get this from data model and make wavelength a member variable.
  TCanvas c1("c1", "c1", 2000, 2000);
  const int numb_of_meas = traces.size();
  int m{0};
  TH2D heat_map("heat_map", "Transparency Measurement Heat Map", numb_of_meas, 0, numb_of_meas-1, N, wav.front(), wav.back());
  for (const auto& trace : traces){
    for (auto i = 0; i < N; ++i){
      heat_map.SetBinContent(m, i, trace.at(i));
    }
    ++m;
  }

  std::string date_of_measure{""};
  for (auto j = 0; j < m; ++j){
    auto this_time = transparency_times.at(j);
    if (date_of_measure != std::to_string(this_time.at(2))+"-"+std::to_string(this_time.at(1))+"-"+std::to_string(this_time.at(0))){
      date_of_measure = std::to_string(this_time.at(2))+"-"+std::to_string(this_time.at(1))+"-"+std::to_string(this_time.at(0));
      heat_map.GetXaxis()->SetBinLabel(j+1, date_of_measure.c_str());
    }
    else {
      heat_map.GetXaxis()->SetBinLabel(j+1, "");
    }
  }

  heat_map.GetXaxis()->SetTitle("date of measurement");
  heat_map.GetYaxis()->SetTitle("wavelength / nm");  
  
  heat_map.Draw("COLZ");
  heat_map.SaveAs(("/home/pi/GDConcMeasure/obj_" + f).c_str());
  c1.SaveAs(("/home/pi/GDConcMeasure/" + f).c_str());
}

