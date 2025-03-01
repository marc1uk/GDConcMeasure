#include "ReturnOfTheMarcusAnalysis.h"
#include <stdexcept>
#include <strings.h> // strcasecmp

#include "TH1.h"
#include "TBufferJSON.h"

ReturnOfTheMarcusAnalysis::ReturnOfTheMarcusAnalysis():Tool(){}

bool ReturnOfTheMarcusAnalysis::Initialise(std::string configfile, DataModel &data){
	
	InitialiseTool(data);
	m_configfile=configfile;
	InitialiseConfiguration(configfile);
	
	verbosity=v_warning;
	
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
	
	//overide with any variables in local file
	Log(m_unique_name+" reading local overrides from config file "+configfile,v_debug,verbosity);
	if(configfile!="") m_variables.Initialise(configfile);
	
	//m_variables.Print();
	
	m_variables.Get("verbosity",verbosity);
	
	// for simplicity each instance of ReturnOfTheMarcusAnalysis only analyses one LED
	m_variables.Get("ledToAnalyse",ledToAnalyse);
	if(ledToAnalyse=="") throw std::runtime_error(m_unique_name+" has no LedToAnalyse!");
	Log(m_unique_name+" will analyse LED "+ledToAnalyse,v_debug,verbosity);
	
	// see if saving traces to ROOT file (debug)
	m_variables.Get("save_trees",save_trees);
	
	// get the transparency of pure water; this is a correction we need to apply to go from
	// the reference arm to what would be expected for the GAD arm with no Gd or contaminants
	GetPureWaterTransparency();
	
	// set up pointers for getting data from Trees
	wavelengthsp = &wavelengths;
	gad_valuesp= &gad_values;
	ref_valuesp= &ref_values;
	gad_darkp= &gad_dark;
	ref_darkp= &ref_dark;
	
	// pointers for writing data output to Trees
	ref_corr_valuesp = &ref_corr_values;
	absorbancesp = &absorbances;
	
	// probably not strictly necessary
	SetGraphTitles();
	
	return get_ok;
}


bool ReturnOfTheMarcusAnalysis::Execute(){
	
	Log(m_unique_name+" Executing...",v_debug,verbosity);
	
	// check if we have an Analyse flag in the CStore indicating intensity data
	// is available for fitting
	if (ReadyToAnalyse()){
		
		Log(m_unique_name+" Processing new measurement...",v_debug,verbosity);
		
		try {
			// reinit datamodel results to prevent accidental carry-over if we bail early
			Log(m_unique_name+" reinitializing variables",v_debug,verbosity);
			ReInit();
			
			// get pointers to the led-on and dark TTrees
			Log(m_unique_name+" getting input trees",v_debug,verbosity);
			GetTrees();
			
			// read in led-on and dark data entries for both GAD and reference arm measurements
			// and do dark subtraction
			Log(m_unique_name+" reading values",v_debug,verbosity);
			ReadValues();
			
			if(make_pureref){
				
				// generate pure water transparency from this measurement
				Log(m_unique_name+" generating pure reference",v_debug,verbosity);
				GeneratePureWaterTransparency();
				
			}
			
			// else {
				
				// calculate absorbance from log10(transmitted / received) light
				Log(m_unique_name+" calculating absorbance",v_debug,verbosity);
				CalculateAbsorbance();
				
			//}
			
			// place results into DataModel for storage
			Log(m_unique_name+" updating DataModel",v_debug,verbosity);
			UpdateDataModel();
			
			// Inform downstream tools that a new measurement is available
			// maybe we could use the value to indicate if the data is good?
			m_data->CStore.Set("NewMarcusAnalyse",ledToAnalyse);
			
		} catch(std::exception& e){
			Log(m_unique_name+" Error! Caught "+e.what(),v_error,verbosity);
			return false;
		}
		
	} else {
		// else no data to Analyse
		// see if there's an old flag from this instance and remove it if so
		std::string lastAnalyse="";
		if(m_data->CStore.Get("NewMarcusAnalyse",lastAnalyse) && lastAnalyse==ledToAnalyse){
			m_data->CStore.Remove("NewMarcusAnalyse");
		}
	}
	
	Log(m_unique_name+" done",v_debug,verbosity);
	
	return true;
}

bool ReturnOfTheMarcusAnalysis::Finalise(){
	if(outfile){
		outfile->Close();
		delete outfile;
		outfile=nullptr;
	}
	return true;
}

bool ReturnOfTheMarcusAnalysis::ReadyToAnalyse(){
	// Checks if analyse flag for our LED has been set by scheduler. Removes it if found.
	bool ready = false;
	std::string analyse="";
	std::string currentLED="";
	m_data->CStore.Get("Analyse", analyse);
	m_data->CStore.Get("ledToAnalyse", currentLED);
	if (analyse == "Analyse" && currentLED == ledToAnalyse){
		m_data->CStore.Remove("Analyse");
		ready = true;
	}
	
	return ready;
}

void ReturnOfTheMarcusAnalysis::SetGraphTitles(){
	
	std::vector<TGraph*>     graphs{ &g_ref,  &g_gad,  &g_ref_corr,  &g_abs, };
	std::vector<std::string> names { "g_ref", "g_gad", "g_ref_corr", "g_abs" };
	for(int i=0; i<graphs.size(); ++i){
		graphs.at(i)->SetName(names.at(i).c_str());
		graphs.at(i)->SetTitle(names.at(i).c_str());
	}
	
	return;
}

// -------------------------------------------------------------------------//

bool ReturnOfTheMarcusAnalysis::GetPureWaterTransparency(){
	// Retrieve pure water transparency from either DB or file
	// prioritize local filename if we have one
	std::string pureref_file;
	get_ok = m_variables.Get("pureref_file",pureref_file);
	if(get_ok){
		if(pureref_file=="create"){
			// create one from first measurement seen.
			// output file is hard-coded - bail if it exists rather than clobbering FIXME?
			TFile fout("PureWaterTransparency.root","OPEN");
			if(!fout.IsZombie()){
				throw std::runtime_error(m_unique_name+" pureref_file 'create' given,"
				      " please remove existing 'PureWaterTransparency.root'");
				fout.Close();
			}
			make_pureref=true;
			return true;
		}
		GetPureWaterTransparency(pureref_file);
	} else {
		// if no filename, see if we have a version number
		// for a database entry
		int pureref_ver=0;
		get_ok = m_variables.Get("pureref_ver",pureref_ver);
		if(!get_ok){
			throw std::runtime_error(m_unique_name+" No pure reference given!");
		}
		GetPureWaterTransparency(pureref_ver);
	}
	
	// set name and title
	std::string purename="g_pureref_"+ledToAnalyse;
	g_pure_absorbance.SetName(purename.c_str());
	g_pure_absorbance.SetTitle(purename.c_str());
	
	// also store a pointer to the graph for plotting on the webpage
	intptr_t puregraphp = reinterpret_cast<intptr_t>(&g_pure_absorbance);
	std::string key = "purerefData_"+ledToAnalyse;
	m_data->CStore.Set(key, puregraphp);
	
	Log(m_unique_name+" loaded pure reference trace of "+g_pure_absorbance.GetN()
	    +" points",v_debug,verbosity);
	
	return true;
}

bool ReturnOfTheMarcusAnalysis::GetPureWaterTransparency(int pureref_ver){
	// get reference trace representing absorption of pure water
	
	// new reference curve can be inserted with e.g.: FIXME update
	// psql -U postgres -d "rundb" -c "INSERT INTO data (timestamp, name, ledname, values) VALUES ('now()', 'pure_curve', '275_A', '{\"version\":0, \"xvals\":[200.0, ..., 800.0], \"yvals\":[0.0079581671, ..., -29982.759] }' );"
	
	// FIXME update
	std::string query_string = "SELECT values->'yvals' FROM data WHERE name='pure_transparency'"
	                           " AND ledname="+m_data->postgres.pqxx_quote(ledToAnalyse)+
	                           " AND values->'yvals' IS NOT NULL"
	                           " AND values->'version' IS NOT NULL"
	                           " AND values->'version'="+ m_data->postgres.pqxx_quote(pureref_ver);
	
	std::string pureref_json="";
	get_ok = m_data->postgres.ExecuteQuery(query_string, pureref_json);
	
	if(!get_ok || pureref_json==""){
		throw std::runtime_error(m_unique_name+" GetPureRefDB obtained empty y array for led "
		                        +ledToAnalyse+", version "+std::to_string(pureref_ver));
	}
	
	// the values string is a json array; i.e. '[val1, val2, val3...]'
	// first strip the '[' and ']' ...
	pureref_json = pureref_json.substr(1,pureref_json.length()-2);
	
	// then parse the remaining list of values
	std::stringstream ss(pureref_json);
	std::string tmp;
	std::vector<double> pureref_yvals;
	
	while(std::getline(ss,tmp,',')){
		char* endptr = &tmp[0];
		double nextval = strtod(tmp.c_str(),&endptr);
		if(endptr==&tmp[0]){
			throw std::runtime_error(m_unique_name+" GetPureRefDB failed to parse y array for led "
			                        +ledToAnalyse+" version "+std::to_string(pureref_ver));
		}
		pureref_yvals.push_back(nextval);
	}
	
	if(pureref_yvals.size()==0){
		throw std::runtime_error(m_unique_name+" GetPureRefDB parsed no y values for led "
		                         +ledToAnalyse+" version "+std::to_string(pureref_ver));
	}
	
	// repeat for the x-values
	query_string = "SELECT values->'xvals' FROM data WHERE name='pure_transparency'"
	               " AND ledname="+m_data->postgres.pqxx_quote(ledToAnalyse)+
	               " AND values->'xvals' IS NOT NULL"
	               " AND values->'version' IS NOT NULL"
	               " AND values->'version'="+ m_data->postgres.pqxx_quote(pureref_ver);
	
	pureref_json="";
	get_ok = m_data->postgres.ExecuteQuery(query_string, pureref_json);
	if(!get_ok || pureref_json==""){
		throw std::runtime_error(m_unique_name+" GetPureRefDB obtained empty x array for led "
		                        +ledToAnalyse+", version "+std::to_string(pureref_ver));
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
			throw std::runtime_error(m_unique_name+" GetPureRefDB failed to parse x array for led "
			                        +ledToAnalyse+" version "+std::to_string(pureref_ver));
		}
		pureref_xvals.push_back(nextval);
	}
	if(pureref_xvals.size()==0){
		throw std::runtime_error(m_unique_name+" GetPureRefDB parsed no x values for led "
		                         +ledToAnalyse+" version "+std::to_string(pureref_ver));
	}
	
	g_pure_absorbance = TGraph(pureref_xvals.size(), pureref_xvals.data(), pureref_yvals.data());
	
	// put the version number used in the CStore for later tools
	std::string key = "purerefID_"+ledToAnalyse;
	m_data->CStore.Set(key, std::to_string(pureref_ver));
	
	return true;
}

bool ReturnOfTheMarcusAnalysis::GetPureWaterTransparency(std::string filename){
	// get pure reference trace from a local file
	Log(m_unique_name+" loading pure reference trace from local file "+filename,v_debug,verbosity);
	
	TFile* puref = nullptr;
	try {
		puref = TFile::Open(filename.c_str(),"READ");
		
		if(puref==nullptr || puref->IsZombie()){
			throw std::runtime_error(m_unique_name+" Error opening pure reference file "+filename);
		}
		
		// returns number of bytes read
		get_ok = puref->ReadTObject(&g_pure_absorbance,"Graph");
		
		// N.B. DO NOT NORMALISE
		// that would introduce a scaling error such that measuring pure water would produce
		// results reporting a constant attenuation across all wavelengths
		// (although even this is not true - it's not uniform? why not?)
		
		if(get_ok<=0){
			throw std::runtime_error(m_unique_name+" failed to read pure reference TGraph 'Graph' from file "
				                     +filename);
		}
		
	} catch (std::exception& e){
		
		// attempt cleanup
		if(puref){
			puref->Close();
			delete puref;
		}
		
		// rethrow
		throw;
	}
	
	// put the version number used in the CStore for later tools
	std::string key = "purerefID_"+ledToAnalyse;
	m_data->CStore.Set(key, filename);
	
	return true;
}


// -------------------------------------------------------------------------//

bool ReturnOfTheMarcusAnalysis::GetTrees(){
	// Get the TTree pointers for the current dark and led traces from the DataModel
	Log(m_unique_name+" getting TTrees",v_debug,verbosity);
	led_tree = nullptr;
	dark_tree = nullptr;
	for(std::pair<const std::string, TTree*>& atree : m_data->m_trees){
		if (atree.first == ledToAnalyse) led_tree = atree.second;
		else if (strcasecmp(atree.first.c_str(), "dark")==0) dark_tree = atree.second;
		if(led_tree && dark_tree) break;
	}
	
	if(!led_tree) throw std::runtime_error(m_unique_name+" Failed to find led tree!");
	if(!dark_tree) throw std::runtime_error(m_unique_name+" Failed to find dark tree!");
	
	return bool(led_tree) && bool(dark_tree);
}

bool ReturnOfTheMarcusAnalysis::ReadBranch(TTree* tree, const std::string& branch, const size_t entry, std::vector<double>* values){
	get_ok = ((tree->SetBranchAddress(branch.c_str(), &values)) >= 0);
	if(!get_ok){
		throw std::runtime_error(m_unique_name+" failed to set address for tree "+tree->GetName()
		      +", branch "+branch);
	}
	// if SetBranchAddress(branch) succeeded we know tree->GetBranch(branch) will return a valid pointer.
	get_ok = tree->GetBranch(branch.c_str())->GetEntry(entry);
	if(get_ok<=0){
		throw std::runtime_error(m_unique_name+" failed to get entry "+std::to_string(entry)
		                        +" from tree "+tree->GetName()+", branch "+branch);
	}
	tree->GetBranch(branch.c_str())->ResetAddress();
	return true;
}

bool ReturnOfTheMarcusAnalysis::ReadValues(){
	// retrieve data for GAD and reference arms and do dark subtraction
	
	// We assume the measurement process is:
	// 1. measure dark
	// 2. measure ref arm
	// 3. measure dark
	// 4. measure gad arm.
	if(dark_tree->GetEntries()<1){
		Log(m_unique_name+" no entries in dark tree!",v_error,verbosity);
		return false;
	}
	if(led_tree->GetEntries()<2){
		Log(m_unique_name+" insufficient entries in '"+ledToAnalyse+"' tree!",v_error,verbosity);
		return false;
	}
	
	Log(m_unique_name+" reading branches",v_debug,verbosity);
	if(wavelengths.size()==0){
		ReadBranch(led_tree, "wavelength", led_tree->GetEntries()-1, wavelengthsp);
	}
	// we assume we did 'dark-ref-dark-gad', since this is what the MarcusScheduler currently does
	// if that changes we'll need to scan for the right dark entry for the reference arm
	// (i.e possibly not GetEntries()-2. the gad arm will always be the last one)
	ReadBranch(led_tree, "value", led_tree->GetEntries()-2, ref_valuesp);
	ReadBranch(led_tree, "value", led_tree->GetEntries()-1, gad_valuesp);
	ReadBranch(dark_tree, "value", dark_tree->GetEntries()-2, ref_darkp);
	ReadBranch(dark_tree, "value", dark_tree->GetEntries()-1, gad_darkp);
	
	if(g_ref.GetN()==0) g_ref.Set(wavelengths.size());
	if(g_gad.GetN()==0) g_gad.Set(wavelengths.size());
	
	// do dark subtraction
	Log(m_unique_name+" doing dark subtraction",v_debug,verbosity);
	try {
		for(size_t i=0; i<wavelengths.size(); ++i){
			gad_values.at(i) -= gad_dark.at(i);
			ref_values.at(i) -= ref_dark.at(i);
			
			if(TMath::IsNaN(gad_values.at(i)) || TMath::IsNaN(ref_values.at(i)) ||
			  !TMath::Finite(gad_values.at(i)) || !TMath::Finite(ref_values.at(i)) ){
				throw std::runtime_error(m_unique_name+" NaN value in trace point "+std::to_string(i)
				     +"gad: "+std::to_string(gad_values.at(i))+", ref: "+std::to_string(ref_values.at(i)));
			}
			
			g_gad.SetPoint(i, wavelengths.at(i), gad_values.at(i));
			g_ref.SetPoint(i, wavelengths.at(i), ref_values.at(i));
		}
	} catch(std::out_of_range& e){
		std::stringstream ss;
		ss << m_unique_name << " Caught " << e.what() << " doing dark subtraction!\n"
		   << "\twavelengths.size() = "+std::to_string(wavelengths.size())<<"\n"
		   << "\tGAD values.size() = "+std::to_string(gad_values.size())<<"\n"
		   << "\tref values.size() = "+std::to_string(ref_values.size())<<"\n"
		   << "\tGAD darks.size() = "+std::to_string(gad_dark.size())<<"\n"
		   << "\tref darks.size() = "+std::to_string(ref_dark.size())<<"\n";
		Log(ss.str(),v_error,verbosity);
		throw std::runtime_error(m_unique_name+" Error getting data from trees");
	}
	
	// for stability monitoring we'll record some characteristic information about the raw data
	// in the database. The dark trace should be pretty flat, so we'll histogram it,
	// fit it with a gaussian, and record the mean and sigma. - do this just for gad arm measurement.
	TH1D tmphist("tmphist","title",100,*std::min_element(gad_dark.begin(), gad_dark.end()),
		                               *std::max_element(gad_dark.begin(), gad_dark.end()));
	for(size_t i=0; i<gad_dark.size(); ++i){
		tmphist.Fill(gad_dark.at(i));
	}
	tmphist.Fit("gaus","Q");
	dark_mean = tmphist.GetFunction("gaus")->GetParameter(1);
	dark_sigma = tmphist.GetFunction("gaus")->GetParameter(2);
	
	// for the raw LED-on trace we'll record the maximum and minimum value of the trace
	// within the absorption region.
	ref_max = *std::max_element(ref_values.begin(), ref_values.end());
	ref_min = *std::min_element(ref_values.begin(), ref_values.end());
	
	gad_max = *std::max_element(gad_values.begin(), gad_values.end());
	gad_min = *std::min_element(gad_values.begin(), gad_values.end());
	
	Log(m_unique_name+" ref arm max: "+std::to_string(ref_max)
	   +", gad arm max: "+std::to_string(gad_max),v_debug,verbosity);
	
	return true;
}

// -------------------------------------------------------------------------//

bool ReturnOfTheMarcusAnalysis::CalculateAbsorbance(){
	// generate absorbance as log of ratio of (corrected) reference arm to GAD arm
	
	// TODO we could fit the reference arm data to the gad arm data in the sidebands.
	// This should be a no-op, but could detect changes
	// 1. from fluctuations in LED output between the two measurements (hopefully small)
	// 2. from variations in absorption down the gad arm since the 'water transparency' reference was taken - e.g. solarization of the fibres
	
	if(g_ref_corr.GetN()==0) g_ref_corr.Set(wavelengths.size());
	if(g_abs.GetN()==0) g_abs.Set(wavelengths.size());
	if(ref_corr_values.size()==0) ref_corr_values.resize(wavelengths.size());
	if(absorbances.size()==0) absorbances.resize(wavelengths.size());
	
	for(size_t i=0; i<wavelengths.size(); ++i){
		// correct reference arm values for water transparency (and other GAD optical path elements)
		// to obtain expected GAD arm measurement for pure water
		double ref_value_corr = ref_values.at(i) * g_pure_absorbance.GetY()[i];
		ref_corr_values[i]=ref_value_corr;
		g_ref_corr.SetPoint(i, wavelengths.at(i), ref_value_corr);
		// we'll get NaN if the argument to log10 is negative; i.e. if either value in the ratio is negative
		// while technically we cannot have negative light, after dark subtraction we can get negative vlaues.
		// if ref arm is <=0, call absorbance 0. If gad arm is <=0, set gad value to 1*
		// if gad arm is > ref arm, this is probably noise, so also set absorbance to 0
		// *a gad value of 0 means ref/gad is inf, so set to 1 ADC count.
		if(gad_values.at(i)<=0) gad_values.at(i)=1.;
		double absval=-1;
		if(ref_value_corr<=0) absval=1.; // assume no absorbance, since nothing to absorb? no generally appropriate value tbh
		else if(gad_values.at(i)>ref_value_corr) absval=1.; // FIXME sanity check that ref value is very small?
		else absval = ref_value_corr/gad_values.at(i); //log10(ref_value_corr/gad_values.at(i));
		if(TMath::IsNaN(absval) || !TMath::Finite(absval)){
			throw std::runtime_error(m_unique_name+" NaN absorbance value "+std::to_string(absval)
			                        +" for datapoint "+ std::to_string(i)
			                        +" from gad value "+std::to_string(gad_values.at(i))+", ref value: "
			                        +std::to_string(ref_values.at(i))+", corrected for pure water: "
			                        + std::to_string(ref_value_corr));
		}
		absorbances[i]=absval;
		g_abs.SetPoint(i, wavelengths.at(i), absval);
	}
	corrected_ref_max = *std::max_element(g_ref_corr.GetY(),g_ref_corr.GetY()+g_ref_corr.GetN());
	
	if(!save_trees) return true;
	static std::vector<double> puretranspvals(g_pure_absorbance.GetY(),g_pure_absorbance.GetY()+g_pure_absorbance.GetN());
	static std::vector<double>* puretranspvalsp = &puretranspvals;
	
	// save calculated traces to output tree if requested
	std::string filename;
	if(strcmp(gDirectory->GetFile()->GetOption(),"READ")==0){
		// if we're processing offline we won't have an open file being written
		// (perhaps we also should put debug info into an alternative file anyway?)
		Log("making rotma file",v_debug,verbosity);
		outfile = new TFile("rotma.root","RECREATE");
	} else if(outfile){
		outfile->cd();
	}
	if(!outtree){
		outtree = new TTree("rotma","Return of the Marcus Analysis");
		outtree->Branch("wavelength",&wavelengths);
		outtree->Branch("gad_values",&gad_values);
		outtree->Branch("ref_values",&ref_values);
		outtree->Branch("pure_transp",&puretranspvals);
		outtree->Branch("ref_corr",&ref_corr_values);
		outtree->Branch("abs",&absorbancesp);
	} else {
		outtree->SetBranchAddress("wavelength",&wavelengthsp);
		outtree->SetBranchAddress("gad_values",&gad_valuesp);
		outtree->SetBranchAddress("ref_values",&ref_valuesp);
		outtree->SetBranchAddress("pure_transp",&puretranspvalsp);
		outtree->SetBranchAddress("ref_corr",&ref_corr_valuesp);
		outtree->SetBranchAddress("abs",&absorbancesp);
	}
	outtree->Fill();
	if(outfile) outfile->Write("",TObject::kOverwrite);
	outtree->ResetBranchAddresses();
	
	if(m_data->CStore.Get("Filename",filename)){
		m_data->m_trees.emplace("rotma",outtree);
		// SaveTraces deletes all entries of m_data->m_trees when save is called
		// so we'll need to make a new one next Execute
		outtree=nullptr;
	}
	
	return true;
}

//==========================================================================//
// below functions could be moved to a new Tool
//==========================================================================//

// -------------------------------------------------------------------------//

void ReturnOfTheMarcusAnalysis::ReInit(){
	
	// Remove outputs from previous Executions, so we don't carry over results
	// from a previous fit if we bail early
	m_data->CStore.Remove("NewMarcusAnalyse");
	
	m_data->CStore.Remove("data_gad");
	m_data->CStore.Remove("data_ref_corrected");
	m_data->CStore.Remove("data_ref");
	
	m_data->CStore.Remove("absorbance_all");
	
	m_data->CStore.Remove("dark_mean");
	m_data->CStore.Remove("dark_sigma");
	m_data->CStore.Remove("raw_ref_max");
	m_data->CStore.Remove("corrected_ref_max");
	
	return;
}

void ReturnOfTheMarcusAnalysis::UpdateDataModel(){
	
	// static parameters from Initialization
	// since these don't change, don't bother re-setting them each time
//	m_data->CStore.Set("purerefID", filename);                                            // filename or DB version ID of pure water transparency trace
//	m_data->CStore.Set("purerefData", puregraphp);                                        // pointer to TGraph of pure water transparency
	
	// results from analysis
	m_data->CStore.Set("NewMarcusAnalysis",ledToAnalyse);                                  // A flag informing downstream tools that new results from this Tool are available
	
	m_data->CStore.Set("data_gad",reinterpret_cast<intptr_t>(&g_gad));                    // plot this
	m_data->CStore.Set("data_ref_corrected",reinterpret_cast<intptr_t>(&g_ref_corr));     // and this on webpage
	m_data->CStore.Set("data_ref",reinterpret_cast<intptr_t>(&g_ref));                    // this can be plotted but hidden by default
	
	m_data->CStore.Set("absorbance_all",reinterpret_cast<intptr_t>(&g_abs));              // full wl range
	
	m_data->CStore.Set("dark_mean",dark_mean);                                            // TODO maybe by storing dark info for both gad & ref,
	m_data->CStore.Set("dark_sigma",dark_sigma);                                          // we could tell if the spectrometer was warming up?
	m_data->CStore.Set("raw_ref_max",ref_max);                                            // led intensity down ref arm
	m_data->CStore.Set("corrected_ref_max",corrected_ref_max);                            // *expected* led intensity down gad arm
	m_data->CStore.Set("gad_max",gad_max);                                                // probably not terribly useful by itself
	
	// also put into monitoring store anything being sent for shift plots
	// e.g. absorbance trace, LED intensities, Gd concentration...
	m_data->monitoring_store.Set("ref_arm_"+ledToAnalyse+"_intensity",ref_max);

	// FIXME this is only useful for the moment until we add Gd, then we will need to do a fit
	// and extract the LED intensity from that. For now we have no Gd absorption so fit is unnecessary.
	m_data->monitoring_store.Set("gad_arm_"+ledToAnalyse+"_intensity",gad_max);
	
	// send absorbance graph
	std::string graph_json = TBufferJSON::ToJSON(&g_abs).Data();
	std::string graph_name = "absorbance_"+ledToAnalyse;
	m_data->services->SendROOTplot(graph_name, "AL", graph_json, true);
	
	return;
	
}

// -------------------------------------------------------------------------//

bool ReturnOfTheMarcusAnalysis::GeneratePureWaterTransparency(){
	// measure absorption along GAD vs reference arm
	// (water, fibres, lenses, other optical path elements, difference in source coupling, etc etc)
	// to obtain expected GAD arm measurement for pure water
	
	TGraph pure_water_transp(wavelengths.size());
	for(size_t i=0; i<wavelengths.size(); ++i){
		
		// we'll get NaN if the argument to log10 is negative; i.e. if either value in the ratio is negative
		// while technically we cannot have negative light, after dark subtraction we can get negative vlaues.
		// if ref arm is <=0, call absorbance 0. If gad arm is <=0, set gad value to 1*
		// if gad arm is > ref arm, this is probably noise, so also set absorbance to 0
		// *a gad value of 0 means ref/gad is inf, so set to 1 ADC count.
		double absval=-99;
		if(ref_values.at(i)<=0) ref_values.at(i)=1.;
		if(gad_values.at(i)>ref_values.at(i)) absval=0;
		else absval = gad_values.at(i) / ref_values.at(i); //log10(gad_values.at(i)/ref_values.at(i));
		if(TMath::IsNaN(absval) || !TMath::Finite(absval)){
			throw std::runtime_error(m_unique_name+" NaN absorbance value "+std::to_string(absval)
			                        +" for datapoint "+ std::to_string(i)
			                        +" from gad value "+std::to_string(gad_values.at(i))+", ref value: "
			                        +std::to_string(ref_values.at(i))
			                        +" in ReturnOfTheMarcusAnalysis::GeneratePureWaterTransparency");
		}
		
		pure_water_transp.SetPoint(i, wavelengths.at(i), absval);
	}
	// N.B. This should *not* be normalised!
	// (no scaling unless we have some way to externally measure the relative intensity of the LED during GAD and Ref arm
	//  measurements of this pure water measurement, in which case we can account for that here... but we don't. so we can't.)
	
	TFile fout("PureWaterTransparency.root","CREATE");
	if(fout.IsZombie()){
		throw std::runtime_error(m_unique_name+" transparency file exists!?");
	}
	pure_water_transp.Write();
	fout.Close();
	make_pureref = false;
	
	m_variables.Set("pureref_file","PureWaterTransparency.root");
	GetPureWaterTransparency();
	
	return true;
}



/*
// unrelated snippet: drawing with smoothing:
rotma->Draw("abs:wavelength","wavelength>400&&wavelength<700","L");
TGraph* g=(TGraph*)c1->GetListOfPrimitives()->At(3);
// or `auto g = (TGraph*)gPad->GetPrimitive("Graph");`
// or draw with "goff" and use `TGraph g(n,tree->GetV1(),tree->GetV2());`
TH1F hist("h","h",g->GetN(),g->GetXaxis()->GetXmin(),g->GetXaxis()->GetXmax());
for(int i=0; i<g->GetN(); ++i){ hist.SetBinContent(i,g->GetY()[i]); }
hist.SetLineColor(kSpring-1);
hist.Smooth(100);
hist.Draw("L");
*/
