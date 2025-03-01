#include "ReturnOfTheMarcusAnalysisEpisode2.h"
#include <stdexcept>

ReturnOfTheMarcusAnalysisEpisode2::ReturnOfTheMarcusAnalysisEpisode2():Tool(){}

bool ReturnOfTheMarcusAnalysisEpisode2::Initialise(std::string configfile, DataModel &data){
	
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
	
	// for simplicity each instance of ReturnOfTheMarcusAnalysisEpisode2 only analyses one LED
	m_variables.Get("ledToAnalyse",ledToAnalyse);
	if(ledToAnalyse=="") throw std::runtime_error(m_unique_name+" has no LedToAnalyse!");
	Log(m_unique_name+" will analyse LED "+ledToAnalyse,v_debug,verbosity);
	
	// get reference Gd absorption shape
	GetAbsorptionRef();
	
	// turn reference graph into functional TF1 fit
	GetAbsFunc();
	
	// get calibration cofficients for converting absorbance to gd concentration
	GetCalibrationCurve();
	
	// set up pointers for getting data from Trees
	wavelengthsp = &wavelengths;
	gad_valuesp= &gad_values;
	ref_valuesp= &ref_values;
	gad_darkp= &gad_dark;
	ref_darkp= &ref_dark;
	
	// see if saving traces to ROOT file (debug)
	m_variables.Get("save_trees",save_trees);
	
	// pointers for writing data output to Trees
	ref_corr_valuesp = &ref_corr_values;
	absorbancesp = &absorbances;
	absfitvaluesp = &absfitvalues;
	
	// probably not strictly necessary
	SetGraphTitles();
	
	return get_ok;
}


bool ReturnOfTheMarcusAnalysisEpisode2::Execute(){
	
	Log(m_unique_name+" Executing...",v_debug,verbosity);
	
	// check if we have an Analyse flag in the CStore indicating intensity data
	// is available for fitting
	if (ReadyToAnalyse()){
		
		Log(m_unique_name+" Processing new measurement...",v_debug,verbosity);
		
		try {
			// reinit datamodel results to prevent accidental carry-over if we bail early
			Log(m_unique_name+" reinitializing variables",v_debug,verbosity);
			ReInit();
			
			// get absorbance data
			Log(m_unique_name+" getting absorbance data",v_debug,verbosity);
			GetAbsorbance();
			
			// fit absorbance trace with reference Gd absorbance shape
			Log(m_unique_name+" fitting absorbance",v_debug,verbosity);
			FitAbsorbance();
			
			// fit absorption peaks to obtain difference and convert to concentration.
			// for each fitting method, calculate the difference in absorbtion peak heights
			// and convert to concentration. Store the results BoostStore map.
			// some of these fits may fail, but we won't abort at this point...
			Log(m_unique_name+" calculating concentration",v_debug,verbosity);
			CalculateConcentration();
			
			// place results into DataModel for storage
			Log(m_unique_name+" updating DataModel",v_debug,verbosity);
			UpdateDataModel();
			
			// Inform downstream tools that a new measurement is available
			// maybe we could use the value to indicate if the data is good?
			m_data->CStore.Set("NewMarcusAnalyseEp2",ledToAnalyse);
			
		} catch(std::exception& e){
			Log(m_unique_name+" Error! Caught "+e.what(),v_error,verbosity);
			return false;
		}
		
	} else {
		// else no data to Analyse
		// see if there's an old flag from this instance and remove it if so
		std::string lastAnalyse="";
		if(m_data->CStore.Get("NewMarcusAnalyseEp2",lastAnalyse) && lastAnalyse==ledToAnalyse){
			m_data->CStore.Remove("NewMarcusAnalyseEp2");
		}
	}
	
	Log(m_unique_name+" done",v_debug,verbosity);
	
	return true;
}

bool ReturnOfTheMarcusAnalysisEpisode2::Finalise(){
	
	return true;
}

bool ReturnOfTheMarcusAnalysisEpisode2::ReadyToAnalyse(){
	
	// XXX XXX XXX NEEDS FIXING TODO XXX FIXME
	// Checks if analyse flag for our LED has been set by scheduler. Removes it if found.
	bool ready = false;
	std::string currentLED="";
	m_data->CStore.Get("NewMarcusAnalyse", currentLED);
	if (currentLED == ledToAnalyse){
		m_data->CStore.Remove("Analyse");
		ready = true;
	}
	
	return ready;
}

void ReturnOfTheMarcusAnalysisEpisode2::SetGraphTitles(){
	
	std::vector<TGraph*>     graphs{ &g_ref,  &g_gad,  &g_ref_corr,  &g_gadfit,  &g_abs,  &g_abs_gd,  &g_absfit  };
	std::vector<std::string> names { "g_ref", "g_gad", "g_ref_corr", "g_gadfit", "g_abs", "g_abs_gd", "g_absfit" };
	for(int i=0; i<graphs.size(); ++i){
		graphs.at(i)->SetName(names.at(i).c_str());
		graphs.at(i)->SetTitle(names.at(i).c_str());
	}
	
	return;
}

// -------------------------------------------------------------------------//

bool ReturnOfTheMarcusAnalysisEpisode2::GetAbsorbance(){
	return true;
}

//==========================================================================//
// below functions could be moved to a new Tool
//==========================================================================//

bool ReturnOfTheMarcusAnalysisEpisode2::GetAbsorptionRef(){
	// Retrieve reference absorption shape from either DB or file
	// prioritize local filename if we have one
	std::string absref_file;
	get_ok = m_variables.Get("absref_file",absref_file);
	if(get_ok){
		GetAbsorptionRef(absref_file);
	} else {
		// if no filename, see if we have a version number
		// for a database entry
		int absref_ver=0;
		get_ok = m_variables.Get("absref_ver",absref_ver);
		if(!get_ok){
			throw std::runtime_error(m_unique_name+" No absorption reference given!");
		}
		GetAbsorptionRef(absref_ver);
	}
	
	// set name and title
	std::string absname="g_absref_"+ledToAnalyse;
	g_absorption_ref.SetName(absname.c_str());
	g_absorption_ref.SetTitle(absname.c_str());
	
	// also store a pointer to the graph for plotting on the webpage
	intptr_t absrefgraphp = reinterpret_cast<intptr_t>(&g_absorption_ref);
	std::string key = "absrefData_"+ledToAnalyse;
	m_data->CStore.Set(key, absrefgraphp);
	
	Log(m_unique_name+" loaded Gd absorption reference trace of "+g_absorption_ref.GetN()
	    +" points",v_debug,verbosity);
	
	return true;
}

bool ReturnOfTheMarcusAnalysisEpisode2::GetAbsorptionRef(int absref_ver){
	// get reference Gd absorption trace
	
	// TODO make a datamodel function for getting curves from the database,
	// this functionality is duplicated in multiple tools
	
	// new reference curve can be inserted with e.g.: FIXME update
	// psql -U postgres -d "rundb" -c "INSERT INTO data (timestamp, name, ledname, values) VALUES ('now()', 'gd_abs_curve', '275_A', '{\"version\":0, \"xvals\":[200.0, ..., 800.0], \"yvals\":[0.0079581671, ..., -29982.759] }' );"
	
	// FIXME update
	std::string query_string = "SELECT values->'yvals' FROM data WHERE name='gd_abs_curve'"
	                           " AND ledname="+m_data->postgres.pqxx_quote(ledToAnalyse)+
	                           " AND values->'yvals' IS NOT NULL"
	                           " AND values->'version' IS NOT NULL"
	                           " AND values->'version'="+ m_data->postgres.pqxx_quote(absref_ver);
	
	std::string absref_json="";
	get_ok = m_data->postgres.ExecuteQuery(query_string, absref_json);
	
	if(!get_ok || absref_json==""){
		throw std::runtime_error(m_unique_name+" GetAbsorptionRef obtained empty y array for led "
		                        +ledToAnalyse+", version "+std::to_string(absref_ver));
	}
	
	// the values string is a json array; i.e. '[val1, val2, val3...]'
	// first strip the '[' and ']' ...
	absref_json = absref_json.substr(1,absref_json.length()-2);
	
	// then parse the remaining list of values
	std::stringstream ss(absref_json);
	std::string tmp;
	std::vector<double> absref_yvals;
	
	while(std::getline(ss,tmp,',')){
		char* endptr = &tmp[0];
		double nextval = strtod(tmp.c_str(),&endptr);
		if(endptr==&tmp[0]){
			throw std::runtime_error(m_unique_name+" GetAbsorptionRef failed to parse y array for led "
			                        +ledToAnalyse+" version "+std::to_string(absref_ver));
		}
		absref_yvals.push_back(nextval);
	}
	
	if(absref_yvals.size()==0){
		throw std::runtime_error(m_unique_name+" GetAbsorptionRef parsed no y values for led "
		                         +ledToAnalyse+" version "+std::to_string(absref_ver));
	}
	
	// repeat for the x-values
	query_string = "SELECT values->'xvals' FROM data WHERE name='absref_curve'"
	               " AND ledname="+m_data->postgres.pqxx_quote(ledToAnalyse)+
	               " AND values->'xvals' IS NOT NULL"
	               " AND values->'version' IS NOT NULL"
	               " AND values->'version'="+ m_data->postgres.pqxx_quote(absref_ver);
	
	absref_json="";
	get_ok = m_data->postgres.ExecuteQuery(query_string, absref_json);
	if(!get_ok || absref_json==""){
		throw std::runtime_error(m_unique_name+" GetAbsorptionRef obtained empty x array for led "
		                        +ledToAnalyse+", version "+std::to_string(absref_ver));
	}
	
	// the values string is a json array; i.e. '[val1, val2, val3...]'
	// strip the '[' and ']'
	absref_json = absref_json.substr(1,absref_json.length()-2);
	// parse it
	ss.clear();
	ss.str(absref_json);
	std::vector<double> absref_xvals;
	while(std::getline(ss,tmp,',')){
		char* endptr = &tmp[0];
		double nextval = strtod(tmp.c_str(),&endptr);
		if(endptr==&tmp[0]){
			throw std::runtime_error(m_unique_name+" GetAbsorptionRef failed to parse x array for led "
			                        +ledToAnalyse+" version "+std::to_string(absref_ver));
		}
		absref_xvals.push_back(nextval);
	}
	if(absref_xvals.size()==0){
		throw std::runtime_error(m_unique_name+" GetAbsorptionRef parsed no x values for led "
		                         +ledToAnalyse+" version "+std::to_string(absref_ver));
	}
	
	// put the version number used in the CStore for later tools
	std::string key="absrefID_"+ledToAnalyse;
	m_data->CStore.Set(key, std::to_string(absref_ver));
	
	return true;
}

bool ReturnOfTheMarcusAnalysisEpisode2::GetAbsorptionRef(std::string filename){
	// get abs reference trace from a local file
	Log(m_unique_name+" loading abs reference trace from local file "+filename,v_debug,verbosity);
	
	TFile* absf = nullptr;
	try {
		absf = TFile::Open(filename.c_str());
		
		if(absf==nullptr || absf->IsZombie()){
			throw std::runtime_error(m_unique_name+" Error opening abs reference file "+filename);
		}
		
		// returns number of bytes read
		get_ok = absf->ReadTObject(&g_absorption_ref,"Graph");
		
		if(get_ok<=0){
			throw std::runtime_error(m_unique_name+" failed to read abs reference TGraph 'Graph' from file "
				                     +filename);
		}
		
	} catch (std::exception& e){
		
		// attempt cleanup
		if(absf){
			absf->Close();
			delete absf;
		}
		
		// rethrow
		throw;
	}
	
	// put the version number used in the CStore for later tools
	std::string key="absrefID_"+ledToAnalyse;
	m_data->CStore.Set(key, filename);
	
	return true;
}

bool ReturnOfTheMarcusAnalysisEpisode2::GetAbsFunc(){
	
	// construct functional fit of reference gd absorption
	Log(m_unique_name+" constructing functional fit TF1 from reference Gd absorbance trace",v_debug,verbosity);
	
	// We'll scale it, and add a linear background to account for contaminants.
	// TODO we could potentially make this a pol2 or pol3 background,
	// but we should constrain it to being very small
	std::string name="f_absfit_"+ledToAnalyse;
	const int n_absfit_pars = 3;
	
	// for reasons explained in MarcusAnalysis, the easiest way to make a functional fit
	// is to make a lambda function that captures a pointer to the fitted TGraph member
	TGraph* g_abs_ref_p = &g_absorption_ref;
	abs_fct = new TF1(name.c_str(),
		[g_abs_ref_p](double* x, double* par) -> double {
			// par [0] = y-scaling
			// par [1] = baseline offset (c)
			// par [2] = baseline gradient (m)
			double abs = par[0]*g_abs_ref_p->Eval(x[0]);
			double baseline = par[2]*(x[0]-276) + par[1];
			return (abs + baseline);
		},
		ROI_min, ROI_max, n_absfit_pars);
	
	// set default parameters
	// TODO is it worth making these configuration parameters?
	// particularly if we don't reset them between Execute loops, probably not...
	std::vector<double> init_params{1,0,0};
	abs_fct->SetParameters(init_params.data());
	
	// set parameter limits
	abs_fct->SetParLimits(0,0,20);         // y scaling
	abs_fct->SetParLimits(1,-0.2,1.2);     // y offset
	abs_fct->SetParLimits(2,-20,20);       // linear baseline gradient
	
	if(!abs_fct->IsValid()){
		throw std::runtime_error(m_unique_name+" GetAbsFunc failed to construct valid TF1");
	}
	Log(m_unique_name+" functional fit TF1 constructed",v_debug,verbosity);
	
	return true;
}

//--------------------------------------------------------------------------//


bool ReturnOfTheMarcusAnalysisEpisode2::GetCalibrationCurve(){
	// get calibration curve from config variables, local file or filedatabase
	// (prioritised in that order)
	
	// look for config file specification
	get_ok = GetCalibrationCurveFromConfigs();
	if(!get_ok) get_ok = GetCalibrationCurveFromFile();
	if(!get_ok) get_ok = GetCalibrationCurveFromDB();
	
	// check we got a curve from somewhere
	if(!get_ok){
		throw std::runtime_error(m_unique_name+" Could not find any calibration curve!");
	}
	
	return false;
}

bool ReturnOfTheMarcusAnalysisEpisode2::GetCalibrationCurveFromConfigs(){
	// Constructs calibration curve from strings in the config file.
	
	// first TF1 formula
	std::string formula_str;
	get_ok = m_variables.Get("cal_tf1_formula",formula_str);
	if(!get_ok){
		// do not throw here, if this doesn't exist, assume we'll get calibration curve elsewhere
		Log(m_unique_name+" no local calibration curve function given",v_error,verbosity);
		return false;
	}
	
	// number of parameters
	int npars=0;
	get_ok = m_variables.Get("cal_n_par",npars);
	if(!get_ok){
		throw std::runtime_error(m_unique_name+" no local calibration num parameters given");
		return false;
	}
	
	// parameter values
	double calib_coefficients[npars];
	std::string params_str="";   // for recording in results DB
	for (auto i = 0; i < npars; ++i){
		get_ok = m_variables.Get("cal_par_" + std::to_string(i), calib_coefficients[i]);
		if(!get_ok){
			throw std::runtime_error(m_unique_name+" Missing local calibration curve parameter "+std::to_string(i));
		}
		params_str += std::to_string(calib_coefficients[i]);
		if(i<(npars-1)) params_str += ", ";
		Log(m_unique_name+" parameter "+std::to_string(i)+" = "
		    +std::to_string(calib_coefficients[i]),v_debug,verbosity);
	}
	
	// construct a TF1 from the formula and parameters given
	std::string calib_name = "f_cal_"+ledToAnalyse;
	calib_curve = TF1(calib_name.c_str(), formula_str.c_str(), 0, 0.25);
	calib_curve.SetParameters(calib_coefficients);
	
	// check the TF1 constructed successfully
	if(!calib_curve.IsValid()){
		throw std::runtime_error(m_unique_name+" calibration curve constructed from formula '"
		    +formula_str+"' and parameters '"+params_str+"' is not valid!");
	}
	Log(m_unique_name+" constructed calibration function from local configs successfully",v_debug,verbosity);
	
	std::string calcurvejson = "{ \"type\":\"configfile\",\"function\":\""+formula_str+"\",\"params\":\""+params_str+"\" }";
	m_data->CStore.Set("calcurve_"+ledToAnalyse,calcurvejson);
	
	return true;
}

bool ReturnOfTheMarcusAnalysisEpisode2::GetCalibrationCurveFromFile(){
	// Constructs calibration curve from file.
	
	// get filename
	std::string filename;
	get_ok = m_variables.Get("cal_curve_file",filename);
	if(!get_ok){
		// do not throw here, if this doesn't exist, assume we'll get calibration curve elsewhere
		Log(m_unique_name+" no calibration curve filename given",v_error,verbosity);
		return false;
	}
	
	TFile* calibf = nullptr;
	try {
		calibf = TFile::Open(filename.c_str());
		
		if(calibf==nullptr || calibf->IsZombie()){
			throw std::runtime_error(m_unique_name+" Error opening calibration curve file "+filename);
		}
		
		// returns number of bytes read
		get_ok = calibf->ReadTObject(&calib_curve,"CalCurve");
		
		if(get_ok<=0){
			throw std::runtime_error(m_unique_name+" failed to read calibration curve TF1 'CalCurve' from file "
				                     +filename);
		}
		
	} catch (std::exception& e){
		
		// attempt cleanup
		if(calibf){
			calibf->Close();
			delete calibf;
		}
		
		// rethrow
		throw;
	}
	
	// check the TF1 is valid
	if(!calib_curve.IsValid()){
		throw std::runtime_error(m_unique_name+" calibration curve retrieved from file '"
		                        +filename+"' is not valid!");
	}
	
	Log(m_unique_name+" constructed calibration function from local file '"+filename
	    +"' successfully",v_debug,verbosity);
	
	std::string calcurvejson = "{ \"type\":\"localfile\",\"filename\":\""+filename+"\" }";
	m_data->CStore.Set("calcurve_"+ledToAnalyse,calcurvejson);
	
	return true;
}

bool ReturnOfTheMarcusAnalysisEpisode2::GetCalibrationCurveFromDB(){
	// Constructs calibration curve from parameters in the database
	
	int calibID;
	get_ok = m_variables.Get("cal_curve_ver",calibID);
	if(!get_ok){
		// do not throw here, if this doesn't exist, assume we'll get calibration curve elsewhere
		Log(m_unique_name+" no DB version number for calibration curve given",v_error,verbosity);
		return false;
	}
	
	// TODO if version number is '-1', query latest version number.
	
	// new calibration curve can be inserted with e.g.:
	// psql -U postgres -d "rundb" -c "INSERT INTO data (timestamp, tool, ledname, name, values) VALUES ('now()', 'ReturnOfTheMarcusAnalysisEpisode2', '275_A', 'calibration_curve', '{\"version\":0, \"formula\":\"pol6\", \"params\":[0.0079581671, 2.7415760, -31.591756, 478.69924, 18682.891, -29982.759] }' );"
	
	// use version number to lookup formula from database
	std::string formula_str="";
	std::string query_string = "SELECT values->>'formula' FROM data WHERE name='calibration_curve'"
	                           " AND tool='ReturnOfTheMarcusAnalysisEpisode2'"
	                           " AND ledname="+m_data->postgres.pqxx_quote(ledToAnalyse)+
	                           " AND values->'formula' IS NOT NULL"
	                           " AND values->'version' IS NOT NULL"
	                           " AND values->'version'="+m_data->postgres.pqxx_quote(calibID);
	get_ok = m_data->postgres.ExecuteQuery(query_string, formula_str);
	if(!get_ok || formula_str==""){
		throw std::runtime_error(m_unique_name+" GetCalibrationCurveFromDB failed to find "
		                         "calibration curve version "+calibID);
	}
	// we must strip enclosing quotations or the TF1 constructor goes berzerk
	//formula_str = formula_str.substr(1,formula_str.length()-2);
	
	// another query for the curve parameters
	// when querying attributes from JSON fields, we need to explicitly ensure the checked
	// attributes exist, or exclude the entry from the search, or the query will fail.
	query_string = " SELECT values->'params' FROM data WHERE name='calibration_curve'"
	               " AND tool='ReturnOfTheMarcusAnalysisEpisode2'"
	               " AND ledname="+m_data->postgres.pqxx_quote(ledToAnalyse)+
	               " AND values->'params' IS NOT NULL"
	               " AND values->'version' IS NOT NULL"
	               " AND values->'version'=" + m_data->postgres.pqxx_quote(calibID);
	std::string params_str;
	get_ok &= m_data->postgres.ExecuteQuery(query_string, params_str);
	
	// check for errors
	if(!get_ok){
		throw std::runtime_error(m_unique_name+" GetCalibrationCurveFromDB failed to retrieve calibration"
		                         " curve parameters for version "+calibID);
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
			throw std::runtime_error(m_unique_name+" GetCalibrationCurveFromDB failed to parse calibration"
			                        " curve parameters for version "+calibID);
		}
		calib_coefficients.push_back(nextval);
	}
	if(calib_coefficients.size()==0){
		throw std::runtime_error(m_unique_name+" GetCalibrationCurveFromDB parsed no parameters from string '"
		                         +params_str+"' for calibration curve version "+calibID);
	}
	
	// construct a TF1 in the map from the formula and parameters given
	std::string calib_name = "f_cal_"+ledToAnalyse;
	calib_curve = TF1(calib_name.c_str(), formula_str.c_str(), 0, 0.4);
	calib_curve.SetParameters(calib_coefficients.data());
	
	// check the TF1 constructed successfully
	if(!calib_curve.IsValid()){
		throw std::runtime_error(m_unique_name+" GetCalibrationCurveFromDB invalid calibration curve"
		                         " for version "+calibID+", constructed from formula '"+formula_str
		                         +"' and parameters '"+params_str+"' is not valid!");
	}
	
	std::string calcurvejson = "{ \"type\":\"DB\",\"ID\":\""+std::to_string(calibID)+"\" }";
	m_data->CStore.Set("calcurve_"+ledToAnalyse,calcurvejson);
	
	return true;
}

// -------------------------------------------------------------------------//

bool ReturnOfTheMarcusAnalysisEpisode2::GetROI(){
	// get array indices corresponding to region of UV LED / Gd absorption
	Log(m_unique_name+" extracting ROI",v_debug,verbosity);
	
	npoints_all = wavelengths.size();
	npoints_gd = 0;
	int ROI_min_light=50; // FIXME make configurable
	for(int i=0; i<npoints_all; ++i){
		if(wavelengths.at(i)>ROI_max || (end_gd>start_gd && (gad_values.at(i)<ROI_min_light || ref_values.at(i)<ROI_min_light))) break;
		if(wavelengths.at(i)<ROI_min) continue;
		if(npoints_gd==0) start_gd=i;
		end_gd=i;
		++npoints_gd;
	}
	Log(m_unique_name+": ROI spans "+std::to_string(ROI_min)+" to "+std::to_string(ROI_max)
	    +" nm, corresponding to indices "+std::to_string(start_gd)+" to "+std::to_string(end_gd),
	    v_debug,verbosity);
	return true;
}

bool ReturnOfTheMarcusAnalysisEpisode2::FitAbsorbance(){
	// fit absorbance in UV region with reference shape and extract the scaling required
	
	// extract ROI
	// for absorbance fit to work best, we need to only fit the region where there is light
	// so do this on every execution, as our range of wavelengths for which this is true can shift
	GetROI();  // find indices of 260nm - 300nm range
	if(absfitvalues.size()==0){
		g_abs_gd.Set(npoints_gd);
		absfitvalues.resize(npoints_gd);
	}
	
	// extract subset of absorbance around Gd region
	for(size_t i=start_gd, j=0; i<end_gd; ++i, ++j){
		g_abs_gd.SetPoint(j, wavelengths.at(i), g_abs.GetY()[i]);
	}
	
	// initialise fit parameters (skip to carry over previous values)
	//abs_fct->SetParameters(absfunc_init_params.data());
	
	absfitresptr = g_abs_gd.Fit(abs_fct,"RNMQS"); // or make a new one and call it 'tmp'
	//absfitresptr = TFitResultPtr((TFitResult*)tmp->Clone());  // i don't know if Clone is required
	
	if(absfitresptr->IsEmpty() || !absfitresptr->IsValid() || absfitresptr->Status()!=0){
		std::string fitstat;
		fitstat += " IsEmpty=" + std::to_string(absfitresptr->IsEmpty());
		fitstat += " IsValid=" + std::to_string(absfitresptr->IsValid());
		fitstat += " Status=" + std::to_string(absfitresptr->Status());
		Log(m_unique_name+" Warning: abs fit status: "+fitstat,v_error,verbosity);
		// do not return false, these are not robust checks of a bad fit.
		// TODO implement a better check based on chi2
		// TODO based on past experience, do the fit multiple times
		//absfit_success = false;
	}
	
	absfit_success =  !absfitresptr->IsEmpty() &&
	                       absfitresptr->IsValid() &&
	                       absfitresptr->Status()==0 &&
	                     ((absfitresptr->Chi2()/absfitresptr->Ndf()) < 10.) &&
	                     !TMath::IsNaN(absfitresptr->GetParams()[0]) &&
	                      (absfitresptr->GetErrors()[0] < 0.5);
	
	// make a TGraph of the fit for the website....
	if(g_absfit.GetN()==0) g_absfit.Set(npoints_gd);
	for (int i = 0; i < npoints_gd; ++i){
		double next_wl = wavelengths[i+start_gd];
		double next_abs = abs_fct->Eval(next_wl);
		if(TMath::IsNaN(next_abs)){
			Log("abs fit function eval to NaN",v_error,verbosity);
			next_abs = 0;
		}
		absfitvalues[i] = next_abs;
		g_absfit.SetPoint(i, next_wl, next_abs);
	}
	
	/*
	TCanvas c_ttmp("c_ttmp","c_ttmp",1024,800);
	g_sideband.SetMarkerColor(kBlue);
	g_inband.SetMarkerColor(kRed);
	abs_fct->SetLineWidth(1);
	abs_fct->SetLineColor(kBlack);
	g_sideband.Draw("AX*");
	g_inband.Draw("same X*");
	abs_fct->Draw("same");
	c_ttmp.SaveAs("absfit.png");
	*/
	
	if(!save_trees) return absfit_success;
	
	TBranch* absbranch = outtree->Branch("absfit",&absfitvaluesp);
	absbranch->Fill();
	outtree->ResetBranchAddresses();
	
	return absfit_success;
}

// -------------------------------------------------------------------------//

bool ReturnOfTheMarcusAnalysisEpisode2::CalculateConcentration(){
	
	// convert to gd concentration based on specified calibration curve
	Log(m_unique_name+" calculating concentration",v_debug,verbosity);
	metric = absfitresptr->Parameter(0);
	gd_conc = calib_curve.GetX(metric);
	
	double metric_err = absfitresptr->GetErrors()[0];
	metric_and_err = std::pair<double,double>{metric, metric_err};
	
	double gd_conc_err = metric_err * calib_curve.Derivative(metric);
	conc_and_err = std::pair<double,double>{gd_conc, gd_conc_err};
	
	if(!save_trees) return get_ok;
	
	TBranch* bptr = nullptr;
	bptr = outtree->Branch("metric",&metric);
	bptr->Fill();
	bptr = outtree->Branch("gd_conc",&gd_conc);
	bptr->Fill();
	
	return get_ok;
}

// -------------------------------------------------------------------------//

void ReturnOfTheMarcusAnalysisEpisode2::ReInit(){
	
	// Remove outputs from previous Executions, so we don't carry over results
	// from a previous fit if we bail early
	m_data->CStore.Remove("NewMarcusAnalyse");
	
	m_data->CStore.Remove("data_gad");
	m_data->CStore.Remove("data_ref_corrected");
	m_data->CStore.Remove("data_ref");
	m_data->CStore.Remove("gad_fit");
	
	m_data->CStore.Remove("absorbance_all");
	m_data->CStore.Remove("absorbance_gd");
	m_data->CStore.Remove("absfit");
	
	m_data->CStore.Remove("absfitresptr");
	m_data->CStore.Remove("absfit_success");
	
	m_data->CStore.Remove("metric_and_err");
	m_data->CStore.Remove("conc_and_err");
	
	m_data->CStore.Remove("dark_mean");
	m_data->CStore.Remove("dark_sigma");
	m_data->CStore.Remove("raw_ref_max");
	m_data->CStore.Remove("corrected_ref_max");
	//m_data->CStore.Remove("gad_max");
	m_data->CStore.Remove("gad_fitted_max");
	
	return;
}

void ReturnOfTheMarcusAnalysisEpisode2::UpdateDataModel(){
	
	// static parameters from Initialization
	// since these don't change, don't bother re-setting them each time
//	m_data->CStore.Set("absrefID", filename);                                             // filename or DB version ID of reference absorption trace
//	m_data->CStore.Set("absrefData", absrefgraphp);                                       // pointer to TGraph of reference absorption trace
//	m_data->CStore.Set("purerefID", filename);                                            // filename or DB version ID of pure water transparency trace
//	m_data->CStore.Set("purerefData", puregraphp);                                        // pointer to TGraph of pure water transparency
	
	// results from analysis
	m_data->CStore.Set("NewMarcusAnalysis",ledToAnalyse);                                  // A flag informing downstream tools that new results from this Tool are available
	
	m_data->CStore.Set("data_gad",reinterpret_cast<intptr_t>(&g_gad));                    // plot this
	m_data->CStore.Set("data_ref_corrected",reinterpret_cast<intptr_t>(&g_ref_corr));     // and this on webpage
	m_data->CStore.Set("data_ref",reinterpret_cast<intptr_t>(&g_ref));                    // this can be plotted but hidden by default
	m_data->CStore.Set("gad_fit",reinterpret_cast<intptr_t>(&g_gadfit));                  // this can be plotted but hidden by default, used mainly for gad arm intensity extraction
	
	m_data->CStore.Set("absorbance_all",reinterpret_cast<intptr_t>(&g_abs));              // full wl range
	m_data->CStore.Set("absorbance_gd",reinterpret_cast<intptr_t>(&g_abs_gd));            // 260-300nm wl range
	m_data->CStore.Set("absfit",reinterpret_cast<intptr_t>(&g_absfit));                   // these two plotted in separate expansion
	
	m_data->CStore.Set("absfitresptr", reinterpret_cast<intptr_t>(&absfitresptr));        // results for DB
	m_data->CStore.Set("absfit_success",absfit_success);                                  // results for DB
	
	m_data->CStore.Set("metric_and_err",reinterpret_cast<intptr_t>(&metric_and_err));     // results for DB
	m_data->CStore.Set("conc_and_err",reinterpret_cast<intptr_t>(&conc_and_err));         // results for DB
	
	m_data->CStore.Set("dark_mean",dark_mean);                                            // TODO maybe by storing dark info for both gad & ref,
	m_data->CStore.Set("dark_sigma",dark_sigma);                                          // we could tell if the spectrometer was warming up?
	m_data->CStore.Set("raw_ref_max",ref_max);                                            // led intensity down ref arm
	m_data->CStore.Set("corrected_ref_max",corrected_ref_max);                            // *expected* led intensity down gad arm
//	m_data->CStore.Set("gad_max",gad_max);                                                // probably not terribly useful by itself
	m_data->CStore.Set("gad_fitted_max",gad_fitted_max);                                  // TODO *~measured* led intensity down gad arm
	                                                                                      // (obtained by fitting corrected ref to gad data in sidebands)
	                                                                                      // this could be used to measure e.g. solarization.
	
	return;
	
}
