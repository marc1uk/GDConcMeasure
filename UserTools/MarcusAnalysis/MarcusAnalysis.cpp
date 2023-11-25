#include "MarcusAnalysis.h"

MarcusAnalysis::MarcusAnalysis():Tool(){}


bool MarcusAnalysis::Initialise(std::string configfile, DataModel &data){
	
	m_data = &data;
	
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
	
	// basically to make handling configurations easier, each instance of MarcusAnalysis
	// will only analyse one LED - but with potentially several calibration curves
	// so, which LED is this Tool to analyse?
	m_variables.Get("ledToAnalyse",ledToAnalyse);
	if(ledToAnalyse==""){
		Log(m_unique_name+" has no ledToAnalyse!",v_error,verbosity);
		return false;
	}
	Log(m_unique_name+" will analyse LED "+ledToAnalyse,v_debug,verbosity);
	
	// get calibration cofficients from local file and/or DB
	get_ok = GetCalibrations();
	if(!get_ok) return false;
	
	// get pure water reference
	get_ok = GetPureTrace();
	if(!get_ok) return false;
	
	// construct a TF1 from the pure trace
	get_ok = GetPureScaledPlusExtras();
	if(!get_ok) return false;
	
	// construct a TF1 for fitting the absorption graph
	get_ok = GetFourGaussianFunc();
	if(!get_ok) return false;
	
	// init name and title of graphs. probably not technically needed.
	std::string gname="g_abs_"+ledToAnalyse;
	g_absorb.SetName(gname.c_str());
	g_absorb.SetTitle(gname.c_str());
	
	std::string pname="g_purefit_"+ledToAnalyse;
	g_purefit.SetName(pname.c_str());
	g_purefit.SetTitle(pname.c_str());
	
	return true;
}


bool MarcusAnalysis::Execute(){
	
	Log(m_unique_name+" Executing...",v_debug,verbosity);
	
	if(ledToAnalyse=="") return false;        // sanity check
	if(calib_curves.size()==0) return false;  // sanity check
	
	// check if we have an Analyse flag in the CStore indicating intensity data
	// is available for fitting
	if (ReadyToAnalyse()){
		
		Log(m_unique_name+" Processing new measurement...",v_debug,verbosity);
		
		// reinit datamodel results to prevent accidental carry-over if we bail early
		Log(m_unique_name+" reinitializing variables",v_debug,verbosity);
		ReInit();
		
		// get pointers to the led-on and dark TTrees
		Log(m_unique_name+" getting input data",v_debug,verbosity);
		get_ok = FindDarkAndLEDTrees();
		if(not get_ok) return false;
		
		// get led-on data
		Log(m_unique_name+" doing dark-subtraction",v_debug,verbosity);
		get_ok = GetDarkSubTrace();
		if(!get_ok) return false;
		
		// fit pure reference to data in sideband
		Log(m_unique_name+" fitting LED reference trace",v_debug,verbosity);
		get_ok = FitPure();
		if(!get_ok) return false;
		
		// calculate absorbance from log10(transmitted / received) light
		Log(m_unique_name+" calculating absorbance",v_debug,verbosity);
		get_ok = CalculateAbsorbance();
		if(not get_ok) return false;
		
		// fit absorption peaks to obtain difference and convert to concentration.
		// for each fitting method, calculate the difference in absorbtion peak heights
		// and convert to concentration. Store the results BoostStore map.
		// some of these fits may fail, but we won't abort at this point...
		for(auto&& method : calib_curves){
			Log(m_unique_name+" calculating concentration with method "+method.first,v_debug,verbosity);
			get_ok &= CalculateConcentration(method.first);
		}
		if(!get_ok){
			Log(m_unique_name+" Failed to extract concentration for some methods!",v_error,verbosity);
		}
		
		// place results into DataModel for storage
		Log(m_unique_name+" updating DataModel",v_debug,verbosity);
		UpdateDataModel();
		
		// Inform downstream tools that a new measurement is available
		// maybe we could use the value to indicate if the data is good?
		m_data->CStore.Set("NewMarcusAnalyse",ledToAnalyse);
		
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

bool MarcusAnalysis::ReadyToAnalyse(){
	// Checks if analyse flag has been set by scheduler,
	// and if the LED to analyse is the one associated with this instance of MarcusAnalysis
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

bool MarcusAnalysis::Finalise(){
	
	return true;
}

bool MarcusAnalysis::FindDarkAndLEDTrees(){
	// Finds the trees for both dark and led traces; then puts them into a (dark,light) pair.
	light_and_dark_trees = std::pair<TTree*, TTree*>{nullptr,nullptr};
	
	// loop over TTrees (traces) in DataModel and extract the dark and signal TTrees
	for (std::pair<const std::string, TTree*>& atree : m_data->m_trees){
		if (atree.first == ledToAnalyse){light_and_dark_trees.second = atree.second;}              // signal trace
		else if (boost::iequals(atree.first, "dark")){light_and_dark_trees.first = atree.second;}  // dark trace
	}
	
	if(light_and_dark_trees.first==nullptr){
		Log(m_unique_name+" Failed to find 'dark' tree in m_data->m_trees!",
		    v_error,verbosity);
		return false;
	}
	if(light_and_dark_trees.second==nullptr){
		Log(m_unique_name+" Failed to find led-on tree "+ledToAnalyse+" in m_data->m_trees!",
		    v_error,verbosity);
		return false;
	}
	
	return true;
}

void MarcusAnalysis::UpdateDataModel(){
	
	// update datamodel
	m_data->CStore.Set("dark_subtracted_data_in",reinterpret_cast<intptr_t>(&g_inband));
	m_data->CStore.Set("dark_subtracted_data_out",reinterpret_cast<intptr_t>(&g_sideband));
	m_data->CStore.Set("purefit",reinterpret_cast<intptr_t>(&g_purefit));
	m_data->CStore.Set("absorbance",reinterpret_cast<intptr_t>(&g_absorb));
	
	m_data->CStore.Set("purefitresptr",reinterpret_cast<intptr_t>(&purefitresptr));
	m_data->CStore.Set("results",reinterpret_cast<intptr_t>(&results));
	
	m_data->CStore.Set("dark_mean",dark_mean);
	m_data->CStore.Set("dark_sigma",dark_sigma);
	m_data->CStore.Set("raw_absregion_max",led_on_max);
	m_data->CStore.Set("raw_absregion_min",led_on_min);
	
}

void MarcusAnalysis::ReInit(){
	// Reinitialize measurement outputs, so we don't carry over
	// results from a previous fit if we bail early
	
	// clear TGraphs
	g_inband.Set(0);
	g_sideband.Set(0);
	g_other.Set(0);
	g_purefit.Set(0);
	g_absorb.Set(0);
	
	// clear characteristic values
	dark_mean = -1;
	dark_sigma = -1;
	led_on_max = -1;
	led_on_min = -1;
	
	// remake TFitResultPtrs
	purefitresptr = TFitResultPtr();
	
	// clear old fit result pointers
	resultsmap.clear();
	
	// clear absorbance fit and concentration results
	for(std::pair<const std::string, BoostStore>& aresult : results){
		//aresult.second.Delete();
		// actually don't remove everything, keep things like calib curve version
		// just clear the results from this fit
		aresult.second.Remove("peak_posns");
		aresult.second.Remove("peak_heights");
		aresult.second.Remove("peak_errs");
		aresult.second.Remove("fitresultptrs");
		aresult.second.Remove("conc_and_err");
	}
	
	UpdateDataModel();
}

// ==============================

bool MarcusAnalysis::GetPureTrace(){
	// Retrieve pure water reference trace
	// first see if we have a local filename
	std::string pureref_file;
	get_ok = m_variables.Get("pureref_file",pureref_file);
	if(get_ok){
		// if so read in the TGraph from that file
		get_ok = GetPureRef(pureref_file);
	} else {
		// if no filename, see if we have a version number
		// for a database entry
		int pureref_ver=0;
		get_ok = m_variables.Get("pureref_ver",pureref_ver);
		if(!get_ok){
			// neither given - no pure reference trace!
			Log(m_unique_name+" No pure reference given!",v_error,verbosity);
			return false;
		}
		get_ok = GetPureRefDB(pureref_ver);
	}
	return get_ok;
}

// ==============================

bool MarcusAnalysis::GetPureRefDB(int pureref_ver){
	// get reference trace representing the pure-water spectrum for this LED
	
	// new reference curve can be inserted with e.g.:
	// psql -U postgres -d "rundb" -c "INSERT INTO data (timestamp, name, ledname, values) VALUES ('now()', 'pure_curve', '275_A', '{\"version\":0, \"xvals\":[200.0, ..., 800.0], \"yvals\":[0.0079581671, ..., -29982.759] }' );"
	
	std::string query_string = "SELECT values->'yvals' FROM data WHERE name='pure_curve'"
	                           " AND ledname="+m_data->postgres.pqxx_quote(ledToAnalyse)+
	                           " AND values->'yvals' IS NOT NULL"
	                           " AND values->'version' IS NOT NULL"
	                           " AND values->'version'="+ m_data->postgres.pqxx_quote(pureref_ver);
	std::string pureref_json="";
	get_ok = m_data->postgres.ExecuteQuery(query_string, pureref_json);
	if(!get_ok || pureref_json==""){
		Log(m_unique_name+" GetPureRefDB obtained empty y array for led "+ledToAnalyse+
		    ", pure reference version "+std::to_string(pureref_ver),v_error,verbosity);
		return false;
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
			Log(m_unique_name+" GetPureRefDB failed to parse y values array for led "
			    +ledToAnalyse+" pure ref version "+std::to_string(pureref_ver),v_error,verbosity);
			return false;
		}
		pureref_yvals.push_back(nextval);
	}
	if(pureref_yvals.size()==0){
		Log(m_unique_name+" GetPureRefDB parsed no y values for led "
			    +ledToAnalyse+" pure ref version "+std::to_string(pureref_ver),v_error,verbosity);
		return false;
	}
	
	// repeat for the x-values
	query_string = "SELECT values->'xvals' FROM data WHERE name='pure_curve'"
	               " AND ledname="+m_data->postgres.pqxx_quote(ledToAnalyse)+
	               " AND values->'xvals' IS NOT NULL"
	               " AND values->'version' IS NOT NULL"
	               " AND values->'version'="+ m_data->postgres.pqxx_quote(pureref_ver);
	pureref_json="";
	get_ok = m_data->postgres.ExecuteQuery(query_string, pureref_json);
	if(!get_ok || pureref_json==""){
		Log(m_unique_name+" GetPureRefDB obtained empty x array for led "+ledToAnalyse+
		    ", pure reference version "+std::to_string(pureref_ver),v_error,verbosity);
		return false;
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
			Log(m_unique_name+" GetPureRefDB failed to parse y values array for led "
			    +ledToAnalyse+" pure ref version "+std::to_string(pureref_ver),v_error,verbosity);
			return false;
		}
		pureref_xvals.push_back(nextval);
	}
	if(pureref_xvals.size()==0){
		Log(m_unique_name+" GetPureRefDB parsed no y values for led "
			    +ledToAnalyse+" pure ref version "+std::to_string(pureref_ver),v_error,verbosity);
		return false;
	}
	
	dark_subtracted_pure = TGraph(pureref_xvals.size(), pureref_xvals.data(), pureref_yvals.data());
	std::string purename="g_pureref_"+ledToAnalyse;
	dark_subtracted_pure.SetName(purename.c_str());
	dark_subtracted_pure.SetTitle(purename.c_str());
	
	// put the version number used in the CStore for later tools
	std::string idkey = "purerefID_"+ledToAnalyse;
	m_data->CStore.Set(idkey, std::to_string(pureref_ver));
	
	// also store a pointer to the graph for plotting on the webpage
	std::string datakey = "purerefData_"+ledToAnalyse;
	intptr_t puregraphp = reinterpret_cast<intptr_t>(&dark_subtracted_pure);
	m_data->CStore.Set(datakey, puregraphp);
	
	
	return true;
}

// ==============================

bool MarcusAnalysis::GetPureRef(std::string filename){
	// get pure reference trace from a local file
	Log(m_unique_name+" loading pure reference trace from local file "+filename,v_debug,verbosity);

	
	TFile* puref = TFile::Open(filename.c_str());
	if(puref==nullptr || puref->IsZombie()){
		std::cerr<<"Error opening pure reference file "<<filename<<std::endl;
		return false;
	}
	TGraph* dark_subtracted_pure_p = (TGraph*)puref->Get("Graph");
	if(dark_subtracted_pure_p==nullptr){
		Log(m_unique_name+" failed to find pure reference TGraph 'Graph' "
		    "in file "+filename,v_error,verbosity);
		return false;
	}
	puref->Close();
	delete puref;
	dark_subtracted_pure = *dark_subtracted_pure_p;
	std::string purename="g_pureref_"+ledToAnalyse;
	dark_subtracted_pure.SetName(purename.c_str());
	dark_subtracted_pure.SetTitle(purename.c_str());
	
	// put the version number used in the CStore for later tools
	std::string idkey = "purerefID_"+ledToAnalyse;
	m_data->CStore.Set(idkey, filename);
	
	// also store a pointer to the graph for plotting on the webpage
	std::string datakey = "purerefData_"+ledToAnalyse;
	intptr_t puregraphp = reinterpret_cast<intptr_t>(&dark_subtracted_pure);
	m_data->CStore.Set(datakey, puregraphp);
	
	Log(m_unique_name+" loaded pure reference trace of "+dark_subtracted_pure.GetN()
	    +" points",v_debug,verbosity);
	
	return true;
}


// ==============================

bool MarcusAnalysis::GetCalibrations(){
	// m_variables should specify one or more calibration curves.
	// scan for calibration curves. Numbering must be consecutive.
	int method_i=0;
	get_ok = false;  // must be overridden by at least one successfully loaded calib function
	Log(m_unique_name+" getting calibration method...",v_debug,verbosity);
	while(true){
		
		// terminate on last method
		std::string methodname;
		if(!m_variables.Get("Method_"+std::to_string(method_i),methodname)) break;
		
		// method name is a key that will be used to look up fitting method.
		// methods should have corresponding code in CalculatePeakHeights
		// and a corresponding calibration curve in the database (or specified in config file)
		if(!(methodname=="raw" || methodname=="simple" || methodname=="complex")){
			Log(m_unique_name+" unknown absorption fitting method "+methodname,v_error,verbosity);
			get_ok = false;
			break;
		}
		Log(m_unique_name+" adding calibration method "+std::to_string(method_i)
		   +", "+methodname,v_debug,verbosity);
		
		// each abosorption fitting method yields a difference in peak heights;
		// we need a calibration curve to map this to concentration. For DB lookup
		// this can be an ID number, but the user may also override in the local
		// config file by specifying an ID of "Local"
		std::string calibID;
		if(!(m_variables.Get("Calib_"+std::to_string(method_i), calibID))){
			Log(m_unique_name+" no calibration ID for method "+methodname,v_error,verbosity);
			get_ok = false;
			break;
		}
		Log(m_unique_name+" method "+std::to_string(method_i)+" using calibration "+calibID,
		    v_debug,verbosity);
		
		// initialise results map
		results[methodname].Set("method",methodname);
		results[methodname].Set("led",ledToAnalyse);
		results[methodname].Set("calcurveID",calibID);
		
		// check for local override
		if(calibID == "Local"){
			// load function and parameters from local file
			Log(m_unique_name+" loading calibration parameters from local file",v_debug,verbosity);
			get_ok = GetCalCurve(methodname, method_i);
			if(not get_ok) break;
			++method_i;
			
		} else {
			// get the version id and query the curve parameters from the DB
			Log(m_unique_name+" loading calibration parameters from database",v_debug,verbosity);
			get_ok &= GetCalCurveDB(methodname, calibID);
			if(not get_ok) break;
			++method_i;
			
		}
		
		
	}
	if(method_i==0){
		Log(m_unique_name+" No calibration curves given!",v_error,verbosity);
	}
	Log(m_unique_name+" loaded "+std::to_string(method_i)+" calibration methods",v_debug,verbosity);
	
	// put in datamodel for database
	intptr_t calib_curves_ptr = reinterpret_cast<intptr_t>(&calib_curves);
	std::string calib_curves_key = "calib_curves_"+ledToAnalyse;
	m_data->CStore.Set(calib_curves_key,calib_curves_ptr);
	
	return get_ok;
}

bool MarcusAnalysis::GetCalCurve(std::string method, int method_i){
	// Constructs calibration curve from parameters in the config file.
	// each calibration curve maps peak height difference to concentration.
	
	std::string formula_str;
	std::string calfunc = "CalibFunc_"+std::to_string(method_i);
	get_ok = m_variables.Get(calfunc,formula_str);
	if(!get_ok){
		Log(m_unique_name+" no local calibration curve function given for method "+std::to_string(method_i),
		    v_error,verbosity);
		return false;
	}
	int npars=0;
	std::string calnparstr = "CalibNpar_"+std::to_string(method_i);
	get_ok = m_variables.Get(calnparstr,npars);
	if(!get_ok){
		Log(m_unique_name+" no local calibration num parameters given for method "+std::to_string(method_i),
		    v_error, verbosity);
		return false;
	}
	Log(m_unique_name+" calibration method "+std::to_string(method_i)+" corresponds to function "
	    +formula_str+" which has "+std::to_string(npars)+" parameters",v_debug,verbosity);
	
	// retrieve calibration coefficients from config file
	double calib_coefficients[npars];
	const std::string calib_prefix = "CalibPar_"+std::to_string(method_i)+"_";
	std::string params_str="";
	for (auto i = 0; i < npars; ++i){
		get_ok &= m_variables.Get(calib_prefix + std::to_string(i), calib_coefficients[i]);
		params_str += std::to_string(calib_coefficients[i]);
		if(i<(npars-1)) params_str += ", ";
		Log(m_unique_name+" parameter "+std::to_string(i)+" = "
		    +std::to_string(calib_coefficients[i]),v_debug,verbosity);
	}
	if(!get_ok){
		Log(m_unique_name+" did not find "+std::to_string(npars)
		    +" local calibration constants for method "+std::to_string(method_i),
		    v_error,verbosity);
		return false;
	}
	
	// construct a TF1 in the map from the formula and parameters given
	std::string calib_name = "f_cal_"+ledToAnalyse+"_"+method;
	calib_curves.insert({method, {calib_name.c_str(), formula_str.c_str(), 0, 0.25}});
	calib_curves.at(method).SetParameters(calib_coefficients);
	
	// check the TF1 constructed successfully
	if(!calib_curves.at(method).IsValid()){
		Log(m_unique_name+" local calibration curve constructed for method "+std::to_string(method_i)
		    +" from formula '"+formula_str+"' and parameters '"+params_str+"' is not valid!",
		    v_error,verbosity);
		calib_curves.erase(method);
		return false;
	}
	Log(m_unique_name+" constructed function successfully",v_debug,verbosity);
	
	results[method].Set("calcurveFunc",formula_str);
	results[method].Set("calcurvePars",params_str);
	
	return true;
}

bool MarcusAnalysis::GetCalCurveDB(std::string method, std::string calibID){
	// Constructs calibration curve from parameters in the database
	
	// new calibration curve can be inserted with e.g.:
	// psql -U postgres -d "rundb" -c "INSERT INTO data (timestamp, tool, ledname, name, values) VALUES ('now()', 'MarcusAnalysis', '275_A', 'calibration_curve', '{\"version\":0, \"formula\":\"pol6\", \"params\":[0.0079581671, 2.7415760, -31.591756, 478.69924, 18682.891, -29982.759] }' );"
	
	// use version number to lookup formula from database
	std::string formula_str="";
	std::string query_string = "SELECT values->>'formula' FROM data WHERE name='calibration_curve'"
	                           " AND tool='MarcusAnalysis'"
	                           " AND ledname="+m_data->postgres.pqxx_quote(ledToAnalyse)+
	                           " AND values->'formula' IS NOT NULL"
	                           " AND values->'method' IS NOT NULL AND values->'version' IS NOT NULL"
	                           " AND values->>'method'="+m_data->postgres.pqxx_quote(method)+
	                           " AND values->'version'="+m_data->postgres.pqxx_quote(calibID);
	get_ok = m_data->postgres.ExecuteQuery(query_string, formula_str);
	if(!get_ok || formula_str==""){
		Log(m_unique_name+" GetCalCurveDB failed to find calibration formula for method "
		    +method+", version "+calibID,v_error,verbosity);
		return false;
	}
	// we must strip enclosing quotations or the TF1 constructor goes berzerk
	//formula_str = formula_str.substr(1,formula_str.length()-2);
	
	// another query for the curve parameters
	// when querying attributes from JSON fields, we need to explicitly ensure the checked
	// attributes exist, or exclude the entry from the search, or the query will fail.
	query_string = " SELECT values->'params' FROM data WHERE name='calibration_curve'"
	               " AND tool='MarcusAnalysis'"
	               " AND ledname="+m_data->postgres.pqxx_quote(ledToAnalyse)+
	               " AND values->'params' IS NOT NULL"
	               " AND values->'method' IS NOT NULL AND values->'version' IS NOT NULL"
	               " AND values->>'method'="+m_data->postgres.pqxx_quote(method)+
	               " AND values->'version'=" + m_data->postgres.pqxx_quote(calibID);
	std::string params_str;
	get_ok &= m_data->postgres.ExecuteQuery(query_string, params_str);
	
	// check for errors
	if(!get_ok){
		Log(m_unique_name+" GetCalCurveDB failed to retrieve calibration curve parameters"
		    " for method "+method+", version "+calibID, v_error,verbosity);
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
			Log(m_unique_name+" GetCalCurveDB failed to parse calibration parameters from string '"
			    +params_str+"' for calibration curve method "+method+", version "+calibID,
			    v_error,verbosity);
			return false;
		}
		calib_coefficients.push_back(nextval);
	}
	if(calib_coefficients.size()==0){
		Log(m_unique_name+" GetCalCurveDB parsed no parameters from string '"
		    +params_str+" for calibration curve method "+method+", version "+calibID,
		    v_error,verbosity);
		return false;
	}
	
	// construct a TF1 in the map from the formula and parameters given
	std::string calib_name = "f_cal_"+ledToAnalyse+"_"+method;
	calib_curves.insert({method, {calib_name.c_str(), formula_str.c_str(), 0, 0.4}});
	calib_curves.at(method).SetParameters(calib_coefficients.data());
	
	// check the TF1 constructed successfully
	if(!calib_curves.at(method).IsValid()){
		Log(m_unique_name+" GetCalCurveDB invalid calibration curve for method "+method
		    +", version "+calibID+", constructed from formula '"+formula_str
		    +"' and parameters '"+params_str+"' is not valid!",v_error,verbosity);
		calib_curves.erase(method);
		return false;
	}
	
	return true;
}


bool MarcusAnalysis::GetPureScaledPlusExtras(){
	
	// construct functional fit.
	Log(m_unique_name+" constructing functional fit TF1 from pure reference trace",v_debug,verbosity);
	
	// We'll scale and add a linear background. We'll also limit the pure function
	// to the region in which we have light - no point fitting outside this region
	std::string name="f_purefit_"+ledToAnalyse;
	const int wave_min = 260, wave_max = 300, numb_of_fitting_parameters = 8;
	
	// pure_fct is a TF1 built from a C++ function representing a transformed version
	// of the pure reference data.
	// if we use a member as the c++ function it must be static, but then can only access static members
	// (which isn't suitable as pure reference is specific to this tool instance).
	// if in anonymous namespace, the resulting function is not a member, so cannot access any members
	// at all (including the pure reference data).
	// the easiest thing is a lambda that captures a pointer to the pure reference TGraph member
	Log(m_unique_name+" constructing functional TF1",v_debug,verbosity);
	TGraph* dark_subtracted_purep = &dark_subtracted_pure;
	pure_fct = new TF1(name.c_str(),
		[dark_subtracted_purep](double* x, double* par) -> double {
			// par [0] = y-scaling
			// par [1] = x-scaling
			// par [2] = x-offset
			// par [3] = y-offset
			// par [4] = linear baseline offset
			// par [5] = shoulder gaussian scaling
			// par [6] = shoulder gaussian centre, restricted to > 282nm (RH shoulder)
			// par [7] = shoulder gaussian spread
			double purepart = (par[0] * dark_subtracted_purep->Eval((par[1]*(x[0]-276))+276+par[2]));
			double linpart = (par[4] * (x[0]-276)) + par[3];
			double shoulderpart = par[5]*exp(-0.5*TMath::Sq((x[0]-282.-abs(par[6]))/par[7]));
			double retval = purepart + linpart + shoulderpart;
			
			return retval;
		},
		wave_min, wave_max, numb_of_fitting_parameters);
	
	
	// set default parameters
	//pure_fct->SetParameters(1.,1.,0.,0.,0.,0.,0.,10.);   // default
	// following from some random entry on 12th march, let's see if it restabilises
	// disabled after april2023 refurb
	//std::vector<double> LED_a_inits{0.749393, 1.002087, 0.269252, -10.00000, -1.254840, -195.183650, 0.00, 9.714654};
	//std::vector<double> LED_b_inits{2.245291, 1.041264, 0.599252, -9.999999, -0.864187, 85.906849, 0.00, 3.079667};
	std::vector<double> LED_a_inits{17.,0.9,3.4,1000.,-50.,-2000.,5.,20.};
	std::vector<double> LED_b_inits = LED_a_inits;
	//std::vector<double> LED_a_inits{18.24,0.79,1.08,-1725.,-189.6,-21385.,0.014,8.8};
	//std::vector<double> LED_b_inits{18.32,0.87,1.91,147.24,-78.59,-5837.,0.,8.88};
	
	pure_init_params = (ledToAnalyse=="275_A") ? LED_a_inits : LED_b_inits;
	pure_fct->SetParameters(pure_init_params.data());
	// set parameter limits
	pure_fct->SetParLimits(0,0,20);     // stretch y
	pure_fct->SetParLimits(1,0.5,1.2);  // stretch x
	pure_fct->SetParLimits(2,-20,20);   // x offset
	pure_fct->SetParLimits(3,-3000,3000);   // y offset
	pure_fct->SetParLimits(4,-500,10);   // linear baseline gradient
	pure_fct->SetParLimits(5,-30000,500); // shoulder gaussian amplitude XXX reduce
	pure_fct->SetParLimits(6,0,60);     // shoulder gaussian position - 282
	pure_fct->SetParLimits(7,0,50);     // shoulder gaussian width
	
	if(!pure_fct->IsValid()){
		Log(m_unique_name+" GetPureScaledPlusExtras failed to construct TF1 from pure reference data",
		    v_error,verbosity);
		return false;
	} else {
		Log(m_unique_name+" functional fit TF1 constructed",v_debug,verbosity);
	}
	
	return true;
	
}

bool MarcusAnalysis::FitTwoGaussians(std::pair<double,double>& peak_heights, std::pair<double,double>& peak_errs, std::pair<double,double>& peak_posns, std::vector<TFitResultPtr>& fitresptrs){
	// given the absorption graph, fit the two main peaks at 273 and 275nm
	// with gaussians, but only over a narrow region around the peak centre,
	// where a gaussian approximation is reasonable
	// TODO move to class so we don't have to construct them repeatedly on each execution
	std::string name1="f_gaus1_"+ledToAnalyse;
	std::string name2="f_gaus2_"+ledToAnalyse;
	TF1 gaus1(name1.c_str(),"gaus",272.5,273.5);
	TF1 gaus2(name2.c_str(),"gaus",275.2,276.2);
	
	double peakval = *std::max_element(g_absorb.GetY(),g_absorb.GetY()+g_absorb.GetN());
	
	// defaults and limits
	//gaus1.SetParameter(0,0.5);
	gaus1.SetParameter(0,peakval);
	gaus1.SetParameter(1,273);
	gaus1.SetParameter(2,0.6);
	gaus1.SetParLimits(0,0.,1.);
	gaus1.SetParLimits(1,272.75,273.25);
	gaus1.SetParLimits(2,0.3,1.);
	
	//gaus2.SetParameter(0,0.3);
	gaus2.SetParameter(0,peakval*0.1);
	gaus2.SetParameter(1,275.65);
	gaus2.SetParameter(2,0.55);
	gaus2.SetParLimits(0,0.,1.);
	gaus2.SetParLimits(1,275.3,276.0);
	gaus2.SetParLimits(2,0.3,1.);
	
	// fit the two gaussians
	fitresptrs[0] = g_absorb.Fit(name1.c_str(),"NRQS");
	fitresptrs[1] = g_absorb.Fit(name2.c_str(),"NRQS");
	// set names of the TFitResultPtrs so they retain info about which function they relate to
	fitresptrs[0]->SetName(name1.c_str());
	fitresptrs[1]->SetName(name2.c_str());
	
	double gausamp1 = gaus1.GetParameter(0);
	double gausamp2 = gaus2.GetParameter(0);
	
	double gaus1amperr = gaus1.GetParError(0);
	double gaus2amperr = gaus2.GetParError(0);
	//std::cout<<"simple peak fit gaus 1 has ampltiude "<<gausamp1<<"+-"<<gaus1amperr<<std::endl;
	//std::cout<<"simple peak fit gaus 2 has ampltiude "<<gausamp2<<"+-"<<gaus2amperr<<std::endl;
	
	double gaus1pos = gaus1.GetParameter(1);
	double gaus2pos = gaus2.GetParameter(1);
	
	peak_posns = std::pair<double,double>{gaus1pos,gaus2pos};
	peak_heights = std::pair<double,double>{gausamp1,gausamp2};
	peak_errs = std::pair<double,double>{gaus1amperr,gaus2amperr};
	
	get_ok = 1;
	if(fitresptrs[0]->IsEmpty() || !fitresptrs[0]->IsValid() || fitresptrs[0]->Status()!=0){
		std::string fitstat;
		fitstat += " IsEmpty=" + std::to_string(fitresptrs[0]->IsEmpty());
		fitstat += " IsValid=" + std::to_string(fitresptrs[0]->IsValid());
		fitstat += " Status=" + std::to_string(fitresptrs[0]->Status());
		Log(m_unique_name+" Two gaussian fit gaus1 fit failed: "+fitstat,v_error,verbosity);
		// this is not a robust check of a bad fit
		//get_ok = 0;
	}
	if(fitresptrs[1]->IsEmpty() || !fitresptrs[1]->IsValid() || fitresptrs[1]->Status()!=0){
		std::string fitstat;
		fitstat += " IsEmpty=" + std::to_string(fitresptrs[1]->IsEmpty());
		fitstat += " IsValid=" + std::to_string(fitresptrs[1]->IsValid());
		fitstat += " Status=" + std::to_string(fitresptrs[1]->Status());
		Log(m_unique_name+" Two gaussian fit gaus2 fit failed: "+fitstat,v_error,verbosity);
		// this is not a robust check of a bad fit
		//get_ok = 0;
	}
	
	if(TMath::IsNaN(peak_heights.first)||TMath::IsNaN(peak_heights.second)){
		Log(m_unique_name+" Two gaussian fit peak heights are NaN",v_error,verbosity);
		// not sure how the database deals with NaN: make life easier by just using 0
		if(TMath::IsNaN(peak_heights.first)) peak_heights.first=0;
		if(TMath::IsNaN(peak_heights.second)) peak_heights.second=0;
		// this IS a robust indicator of a bad fit!
		get_ok = false;
	}
	
	return get_ok;
}

bool MarcusAnalysis::GetFourGaussianFunc(){
	// to really have a reasonable fit to data we in fact need *four* gaussians.
	Log(m_unique_name+" constructing 4-gaussian fit",v_debug,verbosity);
	
	std::string tpname = "f_fourgaussianfunc_"+ledToAnalyse;
	if(fourgaussianfunc) delete fourgaussianfunc;
	fourgaussianfunc = new TF1(tpname.c_str()," [0]*exp(-0.5*TMath::Sq((x[0]-[1])/[2]))"
	                                          "+[0]*[3]*exp(-0.5*TMath::Sq((x[0]-[4])/[5]))"
	                                          "+[0]*[6]*exp(-0.5*TMath::Sq((x[0]-[7])/[8]))"
	                                          "+[0]*[9]*exp(-0.5*TMath::Sq((x[0]-[10])/[11]))",
	                                          270,277);
	
	// we may wish to make the amplitudes of all peaks relative to the amplitude of peak 1
	// by limiting the amplitudes 0-1, this ensures peak 1 is the largest, which helps
	// prevent misfits at low concentratiosn
	
	// par [0] = peak 1 (LH) amplitude
	// par [1] = peak 1 (LH) position
	// par [2] = peak 1 (LH) width
	
	// par [3] = amplitude of peak 1 RH shoulder (relative to peak 1)
	// par [4] = peak 1 RH shoulder position
	// par [5] = peak 1 RH shoulder width
	
	// par [6] = amplitude of peak 1 LH shoulder (relative to peak 1)
	// par [7] = peak 1 LH shoulder position
	// par [8] = peak 1 LH shoulder width
	
	// par [9]  = amplitude of peak 2
	// par [10] = peak 2 position
	// par [11] = peak 2 width
	fourgaussianfunc->SetParName(0,"peak 1 amp");
	fourgaussianfunc->SetParName(1,"peak 1 pos");
	fourgaussianfunc->SetParName(2,"peak 1 wid");
	
	fourgaussianfunc->SetParName(3,"RH shoulder amp");
	fourgaussianfunc->SetParName(4,"RH shoulder pos");
	fourgaussianfunc->SetParName(5,"RH shoulder wid");
	
	fourgaussianfunc->SetParName(6,"LH shoulder amp");
	fourgaussianfunc->SetParName(7,"LH shoulder pos");
	fourgaussianfunc->SetParName(8,"LH shoulder wid");
	
	fourgaussianfunc->SetParName(9,"peak 2 amp");
	fourgaussianfunc->SetParName(10,"peak 2 pos");
	fourgaussianfunc->SetParName(11,"peak 2 wid");
	
	// restrict ranges of position parameters to ensure
	// each gaussian fits the intended peak
	fourgaussianfunc->SetParLimits(1,272.5,273.5);       // peak 1
	fourgaussianfunc->SetParLimits(4,274.,275.);         // peak 1 RH shoulder
	fourgaussianfunc->SetParLimits(7,270,272.5);         // peak 1 LH shoulder
	fourgaussianfunc->SetParLimits(10,275.,276.);        // peak 2
	
	// restrict range of all widths to prevent crazy values
	fourgaussianfunc->SetParLimits(2,0.3,1.);
	fourgaussianfunc->SetParLimits(5,0.3,1.);
	fourgaussianfunc->SetParLimits(8,0.5,1.);
	fourgaussianfunc->SetParLimits(11,0.3,1.);
	
	// restrict ranges of amplitudes based on intended meaning.
	// peak 1 is 0-1 by definition of absorption,
	// others are 0-1 since they're defined as being relative to
	// (and are smaller than) peak 1
	fourgaussianfunc->SetParLimits(0,0.,1.);  // peak 1
	fourgaussianfunc->SetParLimits(3,0.,1.);  // peak 1
	fourgaussianfunc->SetParLimits(6,0.,1.);  // peak 1
	fourgaussianfunc->SetParLimits(9,0.,1.);  // peak 1
	
	if(!fourgaussianfunc->IsValid()){
		Log(m_unique_name+" Error constructing four gaussian fit!",v_error,verbosity);
		return false;
	} else {
		Log(m_unique_name+" 4-gaussian fit constructed",v_debug,verbosity);
	}
	
	return ReInitFourGaussians();
}


bool MarcusAnalysis::ReInitFourGaussians(){
	// set initial positions of the absorption peaks
	fourgaussianfunc->SetParameter(1,272.9);             // peak 1
	fourgaussianfunc->SetParameter(4,274.);              // peak 1 RH shoulder
	fourgaussianfunc->SetParameter(7,271.8);             // peak 1 LH shoulder
	fourgaussianfunc->SetParameter(10,275.66);           // peak 2
	
	// initial widths based on a good fit
	fourgaussianfunc->SetParameter(2,0.5);
	fourgaussianfunc->SetParameter(5,0.5);
	fourgaussianfunc->SetParameter(8,0.6);
	fourgaussianfunc->SetParameter(11,0.55);
	
	// initial amplitudes from a good fit to a concentration of 0.2% Gd.
	// n.b. the default values for amplitudes of peaks 1 and 2 set here
	// will be overrriden based on results from one of the simpler methods anyway.
	fourgaussianfunc->SetParameter(0,0.46);
	fourgaussianfunc->SetParameter(3,0.5);    // absolute amplitude ~0.23
	fourgaussianfunc->SetParameter(6,0.2);    // absolute amplitude ~0.09
	fourgaussianfunc->SetParameter(9,0.6);    // absolute amplitude ~0.27
	
	if(!fourgaussianfunc->IsValid()){
		Log(m_unique_name+" Error setting default parameter values for 4-gaussian fit",
		    v_error,verbosity);
		return false;
	}
	
	return true;
}

bool MarcusAnalysis::FitFourGaussians(std::pair<double,double>& peak_heights, std::pair<double,double>& peak_errs, std::pair<double,double>& peak_posns, std::vector<TFitResultPtr>& fitresptrs){
	ReInitFourGaussians();
	
	// initialise starting values of peak heights with the results
	// from the raw method, but coerced to a minimum of 0.
	peak_heights.first=std::max(0.,peak_heights.first);
	peak_heights.second=std::max(0.,peak_heights.second);
	
	fourgaussianfunc->SetParameter("peak 1 amp",peak_heights.first);
	fourgaussianfunc->SetParameter("peak 2 amp",peak_heights.second/peak_heights.first);
	fitresptrs[0] = g_absorb.Fit(fourgaussianfunc,"RQNS");
	fitresptrs[0]->SetName(fourgaussianfunc->GetName());
	if(fitresptrs[0]->IsEmpty() || !fitresptrs[0]->IsValid() || fitresptrs[0]->Status()!=0){
		std::string fitstat;
		fitstat += " IsEmpty=" + std::to_string(fitresptrs[0]->IsEmpty());
		fitstat += " IsValid=" + std::to_string(fitresptrs[0]->IsValid());
		fitstat += " Status=" + std::to_string(fitresptrs[0]->Status());
		Log(m_unique_name+" Four gaussian fit failed: "+fitstat,v_error,verbosity);
		// this is not a robust check of a bad fit
		//get_ok = false;
	}
	
	// get the height of the peaks as the maximum of the curve around the peak location
	peak_posns.first = fourgaussianfunc->GetMaximumX(272.5,273.5);
	peak_posns.second = fourgaussianfunc->GetMaximumX(275.,276.);
	peak_heights.first  = fourgaussianfunc->Eval(peak_posns.first);
	peak_heights.second  = fourgaussianfunc->Eval(peak_posns.second);
	if(TMath::IsNaN(peak_heights.first)||TMath::IsNaN(peak_heights.second)){
		Log(m_unique_name+" Four gaussian peak heights are NaN",v_error,verbosity);
		// not sure how the database deals with NaN: make life easier by just using 0
		if(TMath::IsNaN(peak_heights.first)) peak_heights.first=0;
		if(TMath::IsNaN(peak_heights.second)) peak_heights.second=0;
		// this IS a robust indicator of a bad fit!
		get_ok = false;
	}
	
	logmessage.str("");
	logmessage.clear();
	logmessage << m_unique_name<<" complex peaks at: "<<peak_posns.first<<", "<<peak_posns.second
	           <<" are "<<peak_heights.first<<" and "<<peak_heights.second
	           <<" with diff "<<peak_heights.first-peak_heights.second
	           <<" and ratio "<<peak_heights.first/peak_heights.second;
	Log(logmessage.str(),v_debug,verbosity);
	
	peak_errs = CalculateError(peak_posns.first, peak_posns.second);
	
	/*
	static int imagenum=0;
	imagenum++;
	TCanvas c_ttmp("c__tmp","c__tmp",1024,800);
	g_absorb.Draw("AXL*");
	fourgaussianfunc->Draw("same");
	c_ttmp.Modified();
	c_ttmp.Update();
	std::string filename="fourgausfit_"+std::to_string(imagenum)+".png";
	c_ttmp.SaveAs(filename.c_str());
	*/
	
	return get_ok;
}

std::pair<double,double> MarcusAnalysis::CalculateError(double peak1_pos, double peak2_pos){
	
	// ok so we have multiple (4) contributions to the function at each peak.
	// to estimate the error on the height of each peak, i'm gonna add in quadrature
	// the error from each contribution to that peak. The error from each contribution is
	// the error on that gaussian's amplitude, but scaled by it's contribution
	// to the function at that peak's position. e.g. if gaussian 4
	// (the peak2 RH shoulder) does not contribute to the amplitude at
	// peak 1, then its error would be scaled to 0. The contribution from the gaussian
	// centered on peak 1 would be scaled to 1, whereas the contribution to error
	// from the shoulders (peaks 2 and 3) would be something between.
	
	double totalerror1=0;
	double totalerror2=0;
	
	// we have 4 gaussians
	for(int i=0; i<4; ++i){
		
		// error on amplitude of this gaus
		double error_centre = fourgaussianfunc->GetParError(i*3);
		// note that amplitudes of gaussians 1-3 are relative to that of gaussian 0
		// so an error of 0.5 on gaussian 1 is actually only an error of 0.05
		// in the case that gaussian 0 peak height is 0.1
		if(i!=0) error_centre *= fourgaussianfunc->GetParameter(0);
		
		// get the relative contribution to the function from this gaussian
		// first get the height at it's centre...
		double amp_centre;
		if(i==0){
			amp_centre = fourgaussianfunc->GetParameter(0);
		} else {
			amp_centre = fourgaussianfunc->GetParameter(i*3)*fourgaussianfunc->GetParameter(0);
		}
		// now make a temporary gaus to obtain its height at peak 1
		TF1 agaus("agaus","[0]*exp(-0.5*TMath::Sq((x-[1])/[2]))",270., 300.);
		// copy over its parameters...
		for(int j=0; j<3; ++j){
			if(j==0 && i!=0){
				agaus.SetParameter(j,fourgaussianfunc->GetParameter(i*3+j)*fourgaussianfunc->GetParameter(0));
			} else {
				agaus.SetParameter(j,fourgaussianfunc->GetParameter(j+i*3));
			}
		}
		// evaluate at peak 1
		double amp_here = agaus.Eval(peak1_pos);
		double relative_amp = amp_here/amp_centre;
		// scale the error (on the ampltiude at the center) by the relative height
		// at peak 1, compared to at the centre
		double error_here = error_centre*relative_amp;
		
		/*
		std::cout<<"complex peak fit gaus "<<i<<" has ampltiude "
		         <<agaus.GetParameter(0)<<"+-"<<error_centre<<" scaled by "
		         <<relative_amp<<" to "<<error_here<<" at peak 1";
		*/
		
		totalerror1+= TMath::Sq(error_here);
		
		// repeat for contribution to peak2
		amp_here = agaus.Eval(peak2_pos);
		relative_amp = amp_here/amp_centre;
		// scale the error on the amplitude at the centre
		error_here = error_centre*relative_amp;
		
		//std::cout<<" and scaled by "<<relative_amp<<" to "<<error_here<<" at peak 2"<<std::endl;
		
		totalerror2+= TMath::Sq(error_here);
	}
	
	// take sqrt of totalerror1 and totalerror2 each to get quadrature sum of contributions...
	// but then take square to add the errors on the two peak amplitudes in quadrature
	//std::cout<<"total error at peak 1 "<<sqrt(totalerror1)<<" and at peak 2 "<<sqrt(totalerror2)<<std::endl;
	
	return std::pair<double,double>{sqrt(totalerror1),sqrt(totalerror2)};
}

////////////////////////

bool MarcusAnalysis::GetDarkSubTrace(){
	// Takes the dark and led trees and returns a TGraphErrors of the dark subtracted trace
	
	TTree* darkTree = light_and_dark_trees.first;
	TTree* ledTree = light_and_dark_trees.second;
	
	// sanity check
	if(darkTree->GetEntries()==0){
		Log(m_unique_name+" no entries in dark tree!",v_error,verbosity);
		return false;
	}
	if(ledTree->GetEntries()==0){
		Log(m_unique_name+" no entries in led-on tree!",v_error,verbosity);
		return false;
	}
	
	// retrieve data from TTrees
	std::vector<double> led_values, wavelengths, led_errors;
	std::vector<double> dark_values, dark_errors;
	
	std::vector<double>* led_valuesp= &led_values;
	std::vector<double>* wavelengthsp = &wavelengths;
	std::vector<double>* led_errorsp = &led_errors;
	std::vector<double>* dark_valuesp = &dark_values;
	std::vector<double>* dark_errorsp = &dark_errors;
	
	Log(m_unique_name+" setting branch addresses to retrieve data",v_debug,verbosity);
	bool failure=false;
	failure |= ((ledTree->SetBranchAddress("value", &led_valuesp)) < 0);
	failure |= ((ledTree->SetBranchAddress("wavelength", &wavelengthsp)) < 0);
	failure |= ((ledTree->SetBranchAddress("error", &led_errorsp)) < 0);
	if(failure){
		Log(m_unique_name+" failed to set branch addresses for "+ledToAnalyse+" tree!",v_error,verbosity);
		return false;
	}
	failure |= ((darkTree->SetBranchAddress("value", &dark_valuesp)) < 0);
	failure |= ((darkTree->SetBranchAddress("error", &dark_errorsp)) < 0);
	if(failure){
		Log(m_unique_name+" failed to set branch addresses for dark tree!",v_error,verbosity);
		return false;
	}
	
	Log(m_unique_name+" getting TTree entries",v_debug,verbosity);
	get_ok = ledTree->GetEntry(ledTree->GetEntries()-1);
	if(get_ok<=0){
		Log(m_unique_name+" failed to GetEntry of ledTree!",v_error,verbosity);
		return false;
	}
	get_ok = darkTree->GetEntry(darkTree->GetEntries()-1);
	if(get_ok<=0){
		Log(m_unique_name+" failed to GetEntry of darkTree!",v_error,verbosity);
		return false;
	}
	// after getting entries we should reset the branch addresses
	// as we just pointed them at local variables that will soon be out of scope
	ledTree->ResetBranchAddresses();
	darkTree->ResetBranchAddresses();
	
	// sanity check that we got the number of datapoints we expected
	const int number_of_points = 2068;  // expected trace size TODO maybe don't hard-code this...
	if( led_values.size()       != number_of_points ||
		wavelengths.size()      != number_of_points ||
		led_errors.size()       != number_of_points ||
		dark_values.size()      != number_of_points ||
		dark_errors.size()      != number_of_points ){
		Log(m_unique_name+" failed to get data from TTrees!",v_error,verbosity);
		return false;
	}
	
	// do dark subtraction
	for(int i=0; i<led_values.size(); ++i){
		led_values.at(i) -= dark_values.at(i);
	}
	
	// split into regions: in-band (absorption region)
	// sideband (fit region), and other (ignored region)
	// we do dark subtraction here while we're at it
	// n.b. other isn't used for now...
	std::vector<double> inband_values;
	std::vector<double> inband_errors;
	std::vector<double> inband_wls;
	std::vector<double> sideband_values;
	std::vector<double> sideband_errors;
	std::vector<double> sideband_wls;
	std::vector<double> other_values;
	std::vector<double> other_errors;
	std::vector<double> other_wls;
	sample_273=0;
	sample_276=0;
	for(int j=0; j<number_of_points; ++j){
		if(TMath::IsNaN(led_values.at(j))){
			Log("Nan value in data trace point "+std::to_string(j),v_error,verbosity);
		}
		if(wavelengths.at(j)>270 && wavelengths.at(j)<280){
			// in band
			inband_values.push_back(led_values.at(j));
			inband_wls.push_back(wavelengths.at(j));
			inband_errors.push_back(sqrt(led_errors[j]*led_errors[j] + dark_errors[j]*dark_errors[j]));
		} else if(wavelengths.at(j)>260 && wavelengths.at(j)<300){
			// side band (lobes)
			sideband_values.push_back(led_values.at(j));
			sideband_wls.push_back(wavelengths.at(j));
			sideband_errors.push_back(sqrt(led_errors[j]*led_errors[j] + dark_errors[j]*dark_errors[j]));
		} else {
			// other data (not used)
			other_values.push_back(led_values.at(j));
			other_wls.push_back(wavelengths.at(j));
			other_errors.push_back(sqrt(led_errors[j]*led_errors[j] + dark_errors[j]*dark_errors[j]));
		}
		
		// make a note of the index of the 273 and 276 peaks
		// within the raw in-band absorption graph arrays
		if(wavelengths.at(j)>272.8 && sample_273==0){
			sample_273 = inband_values.size()-1;
		} else if(wavelengths.at(j)>275.3 && sample_276==0){
			sample_276 = inband_values.size()-1;
		}
	}
	
	std::vector<double> inband_wl_errors(inband_wls.size(),0.5*(wavelengths[1]-wavelengths[0]));
	std::vector<double> sideband_wl_errors(sideband_wls.size(),0.5*(wavelengths[1]-wavelengths[0]));
	std::vector<double> other_wl_errors(other_wls.size(),0.5*(wavelengths[1]-wavelengths[0]));
	
	// construct TGraphErrors from dark-subtracted data in each region
	// i believe that ROOT does not manage a collection of TGraphs
	// (unless they're contour graphs that are part of a Draw call with "CONT, LIST" options)
	// https://root.cern.ch/root/htmldoc/guides/users-guide/ObjectOwnership.html#ownership-by-the-master-troot-object-groot
	// so we shouldn't have to worry about leaking memory just by replacing it.
	g_inband = TGraphErrors(inband_wls.size(), inband_wls.data(), inband_values.data(), inband_wl_errors.data(), inband_errors.data());
	g_sideband = TGraphErrors(sideband_wls.size(), sideband_wls.data(), sideband_values.data(), sideband_wl_errors.data(), sideband_errors.data());
	g_other = TGraphErrors(other_wls.size(), other_wls.data(), other_values.data(), other_wl_errors.data(), other_errors.data());
	// actually we seem to get better pure fits if we use no errors!
	g_sideband_noerrs = TGraph(sideband_wls.size(), sideband_wls.data(), sideband_values.data());
	
	std::string intitle = "g_inband_"+ledToAnalyse;
	std::string sidetitle = "g_sideband_"+ledToAnalyse;
	std::string othertitle = "g_other_"+ledToAnalyse;
	g_inband.SetTitle(intitle.c_str());
	g_inband.SetName(intitle.c_str());
	g_sideband.SetTitle(sidetitle.c_str());
	g_sideband.SetName(sidetitle.c_str());
	g_other.SetTitle(othertitle.c_str());
	g_other.SetName(othertitle.c_str());
	
	// for stability monitoring we'll record some characteristic information about the raw data
	// in the database. The dark trace should be pretty flat, so we'll histogram it,
	// fit it with a gaussian, and record the mean and sigma.
	TH1D tmphist("tmphist","title",100,*std::min_element(dark_values.begin(), dark_values.end()),
		                               *std::max_element(dark_values.begin(), dark_values.end()));
	for(int i=0; i<number_of_points; ++i){
		tmphist.Fill(dark_values.at(i));
	}
	//std::string options = (verbosity<v_debug) ? "Q" : "";
	tmphist.Fit("gaus","Q");
	dark_mean = tmphist.GetFunction("gaus")->GetParameter(1);
	dark_sigma = tmphist.GetFunction("gaus")->GetParameter(2);
	
	// for the raw LED-on trace we'll record the maximum and minimum value of the trace
	// within the absorption region.
	led_on_max = *std::max_element(inband_values.begin(), inband_values.end());
	led_on_min = *std::min_element(inband_values.begin(), inband_values.end());
	
	return true;
}

//////////////////////

bool MarcusAnalysis::FitPure(){
	
	// fit sideband with pure function to account for led power fluctuations and remove background
	pure_fct->SetParameters(pure_init_params.data());
	//purefitresptr = g_sideband.Fit(pure_fct,"RNMQS");
	purefitresptr = g_sideband_noerrs.Fit(pure_fct,"RNMQS");
	bool fit_ok = true;
	if(purefitresptr->IsEmpty() || !purefitresptr->IsValid() || purefitresptr->Status()!=0){
		std::string fitstat;
		fitstat += " IsEmpty=" + std::to_string(purefitresptr->IsEmpty());
		fitstat += " IsValid=" + std::to_string(purefitresptr->IsValid());
		fitstat += " Status=" + std::to_string(purefitresptr->Status());
		Log(m_unique_name+" pure fit to sideband region failed: "+fitstat,v_error,verbosity);
		// do not return false, these are not robust checks of a bad fit.
		//fit_ok = false;
	}
	
	// make a TGraph of the scaled pure fit for the website....
	double wlerr = g_inband.GetEX()[0];
	g_purefit.Set(g_inband.GetN()+g_sideband.GetN());
	for (int i = 0, j=0, k=0; k < g_inband.GetN()+g_sideband.GetN(); ++k){
		double next_inband_wl = g_inband.GetX()[i];
		double next_sideband_wl = g_sideband.GetX()[j];
		double next_wl;
		if(next_inband_wl < next_sideband_wl && i < g_inband.GetN() ){
			next_wl = next_inband_wl;
			++i;
		} else {
			next_wl = next_sideband_wl;
			++ j;
		}
		if(TMath::IsNaN(pure_fct->Eval(next_wl))){
			Log("Pure fit function eval to NaN",v_error,verbosity);
		}
		g_purefit.SetPoint(k, next_wl, pure_fct->Eval(next_wl));
		double fiterr = 0.00; // TODO FIXME how do we calculate this
		g_purefit.SetPointError(k, wlerr, fiterr);
	}
	
	/*
	TCanvas c_ttmp("c_ttmp","c_ttmp",1024,800);
	g_sideband.SetMarkerColor(kBlue);
	g_inband.SetMarkerColor(kRed);
	pure_fct->SetLineWidth(1);
	pure_fct->SetLineColor(kBlack);
	g_sideband.Draw("AX*");
	g_inband.Draw("same X*");
	pure_fct->Draw("same");
	c_ttmp.SaveAs("purefit.png");
	*/

	return fit_ok;
}

bool MarcusAnalysis::CalculateAbsorbance(){
	
	// calculate absorbance from ratio of fit to data in absorption region
	g_absorb.Set(g_inband.GetN());
	for(int i=0; i<g_inband.GetN(); ++i){
		double wlval, dataval;
		g_inband.GetPoint(i, wlval, dataval);
		double fitval = pure_fct->Eval(wlval);
		double absval = log10(fitval/dataval);
		if(TMath::IsNaN(absval)){
			logmessage.str("");
			logmessage << "absorbance calculation NaN value: fitval: "<<fitval<<", dataval: "<<dataval
			           <<" at wavelength "<<wlval<<std::endl;
			Log(logmessage.str(),v_error,verbosity);
			absval=0;
		}
		g_absorb.SetPoint(i, wlval, absval);
		
		// error on ratio is sqrt((Δa/a)²+(Δb/b)²)
		// https://www.statisticshowto.com/error-propagation/
		// error on logY(X) is (ΔX/X)*(1/ln(Y))
		// https://physics.stackexchange.com/questions/95254/the-error-of-the-natural-logarithm/95278#95278
		// at least in the regime where ΔX<<X. if outside this regime the non-linear
		// nature of logs means the errors will be asymmetric. In this case one can use
		// as a rough guide ΔY = log(X-ΔX) -> log(X+ΔX).
		double error_on_data = g_inband.GetEY()[i];
		// uhhh what's the error on the extrapolated fit value - i.e. transmitted intensity.
		// no idea. Let's assume it's of the same order as the data??? FIXME
		double error_on_fitval = error_on_data;
		double err_on_ratio = sqrt(TMath::Sq(error_on_data/dataval)
		                          +TMath::Sq(error_on_fitval/fitval));
		double err_on_ratio_over_ratio = err_on_ratio / (fitval/dataval);
		g_absorb.SetPointError(i, g_inband.GetEX()[0]/2., err_on_ratio_over_ratio*(1./log(10.)));
	}
	
	return true;
}

bool MarcusAnalysis::CalculateConcentration(std::string method){
	
	BoostStore& result = results[method];
	
	// extract the peak heights and errors
	Log(m_unique_name+" calculating peak heights for method "+method,v_debug,verbosity);
	get_ok = CalculatePeakHeights(method);
	result.Set("absfit_success", get_ok);
	
	// convert to gd concentration based on specified calibration curve
	Log(m_unique_name+" calculating concentration for method "+method,v_debug,verbosity);
	if(get_ok) get_ok = PeakDiffToConcentration(method);
	result.Set("conc_fit_success",get_ok);
	
	return get_ok;
}

bool MarcusAnalysis::CalculatePeakHeights(std::string method){
	
	BoostStore& result = results[method];
	std::pair<double,double> peak_posns;
	std::pair<double,double> peak_heights;
	std::pair<double,double> peak_errs;
	resultsmap[method] = std::vector<TFitResultPtr>{};
	std::vector<TFitResultPtr>& fitresptrs = resultsmap[method];
	get_ok = true;
	
	// raw peaks: just the graph value at the peak wavelengths
	if(method=="raw" || method=="complex"){
		g_absorb.GetPoint(sample_273, peak_posns.first, peak_heights.first);
		g_absorb.GetPoint(sample_276, peak_posns.second, peak_heights.second);
		if(TMath::IsNaN(peak_heights.first) || TMath::IsNaN(peak_heights.second)){
			logmessage <<"raw peak heights are NaN! samples: "<<sample_273<<", "<<sample_276
			           <<" with x: "<<peak_posns.first<<", "<<peak_posns.second<<" had values"
			           <<peak_heights.first<<" and "<<peak_heights.second<<std::endl
			           <<"absorption graph had "<<g_absorb.GetN()<<" points"<<std::endl;
			Log(logmessage.str(),v_error,verbosity);
		}
		
		// if peaks are negative, coerce to 0
		//peak_heights.first=std::max(0.,peak_heights.first);
		//peak_heights.second=std::max(0.,peak_heights.second);
		
		peak_errs = std::pair<double,double>{g_absorb.GetEY()[sample_273],g_absorb.GetEY()[sample_276]};
	}
	
	// simple peaks: fit the peaks with gaussians
	if(method=="simple"){
		fitresptrs.resize(2);
		get_ok = FitTwoGaussians(peak_heights, peak_errs, peak_posns, fitresptrs);
	}
	
	// complex peaks: fit entire region with a combination of 4 gaussians
	if(method=="complex"){
		fitresptrs.resize(1);
		get_ok = FitFourGaussians(peak_heights, peak_errs, peak_posns, fitresptrs);
	}
	
	Log(m_unique_name+" peak heights for "+method+" are "+
	    std::to_string(peak_heights.first)+", "+std::to_string(peak_heights.second),
	    v_debug,verbosity);
	
	result.Set("peak_posns", peak_posns);
	result.Set("peak_heights", peak_heights);
	result.Set("peak_errors", peak_errs);
	//result.Set("fitresultptrs",fitresptrs);
	// ughhh can't do that because TFitResultPtr isn't serialisable
	result.Set("fitresultptrs", reinterpret_cast<intptr_t>(&resultsmap[method]));
	
	return get_ok;
}

bool MarcusAnalysis::PeakDiffToConcentration(std::string method){
	
	BoostStore& result = results[method];
	std::pair<double,double> peak_heights;
	std::pair<double,double> peak_errors;
	get_ok = result.Get("peak_heights", peak_heights);
	get_ok &= result.Get("peak_errors", peak_errors);
	if(!get_ok){
		
		Log(m_unique_name+" Absorbance fit result did not contain peak heights and errors",
		    v_error,verbosity);
		
	} else {
		
		double peakheightdiff = peak_heights.first - peak_heights.second;
		
		// solve for concentration (x) from absorbance (y), with 0.01 < x < 0.21
		TF1& calib_curve = calib_curves.at(method);
		double conc = calib_curve.GetX(peakheightdiff, 0.001, 0.25);
		Log(m_unique_name+" concentration for method "+method+" is "+std::to_string(conc),v_debug,verbosity);
		
		// FIXME TODO
		// we have errors on the two peak heights but they may be correlated...
		// so calculation of the errors is probably going to require functions specialised
		// for the corresponding fit method
		// as a rough first estimate, assume uncorrelated
		double peakheightdifferr = TMath::Sqrt(TMath::Sq(peak_errors.first)+TMath::Sq(peak_errors.second));
		// error on concentration is error on height difference times gradient at that point.
		double conc_err = peakheightdifferr * calib_curve.Derivative(peakheightdiff);
		
		std::pair<double,double> conc_and_err{conc, conc_err};
		result.Set("conc_and_err", conc_and_err);
	}
	
	return get_ok;
}

