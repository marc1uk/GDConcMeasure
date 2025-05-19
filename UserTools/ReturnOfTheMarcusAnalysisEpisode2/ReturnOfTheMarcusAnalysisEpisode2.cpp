#include "ReturnOfTheMarcusAnalysisEpisode2.h"
#include <stdexcept>
#include "TBufferJSON.h"

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
	
	// pointers for writing data output to Trees
	absfitvaluesp = &absfitvalues;
	bgfitvaluesp = &bgfitvalues;
	
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
			
			// get debug tree if requested
			if(save_trees && m_data->m_trees.count("rotma")){
				outtree = m_data->m_trees.at("rotma");
			} else {
				save_trees=false;
			}
			
			// get absorbance data
			Log(m_unique_name+" getting absorbance data",v_debug,verbosity);
			GetAbsorbance();
			
			// fit absorbance trace with reference Gd absorbance shape
			Log(m_unique_name+" fitting background absorbance",v_debug,verbosity);
			RemoveBackgroundAbsorbance();
			Log(m_unique_name+" fitting gd absorbance",v_debug,verbosity);
			FitAbsorbance(true);
			
			// fit absorption peaks to obtain difference and convert to concentration.
			// for each fitting method, calculate the difference in absorbtion peak heights
			// and convert to concentration. Store the results BoostStore map.
			// some of these fits may fail, but we won't abort at this point...
			Log(m_unique_name+" calculating concentration",v_debug,verbosity);
			CalculateConcentration();
			
			// place results into DataModel for storage
			Log(m_unique_name+" updating DataModel",v_debug,verbosity);
			UpdateDataModel();
			
			// write debug tree
			if(save_trees){
				outtree->Fill();
				TFile* outfile = outtree->GetCurrentFile();
				outfile->Write("",TObject::kOverwrite);
				outtree->ResetBranchAddresses();
			}
			
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
	
	// check if ReturnOfTheMarcusAnalysis tool has set flag indicating results
	bool ready = false;
	std::string currentLED="";
	m_data->CStore.Get("NewMarcusAnalyse", currentLED);
	if (currentLED == ledToAnalyse){
		ready = true;
	}
	
	return ready;
}

void ReturnOfTheMarcusAnalysisEpisode2::SetGraphTitles(){
	
	
	std::vector<TGraph*>     graphs{&g_abs_gd,  &g_absfit  };
	std::vector<std::string> names {"g_abs_gd", "g_absfit" };
	for(int i=0; i<graphs.size(); ++i){
		graphs.at(i)->SetName(names.at(i).c_str());
		graphs.at(i)->SetTitle(names.at(i).c_str());
	}
	
	return;
}

// -------------------------------------------------------------------------//

bool ReturnOfTheMarcusAnalysisEpisode2::GetAbsorbance(){
	
	intptr_t g_ptr=0;
	get_ok = m_data->CStore.Get("absorbance_all",g_ptr);
	TGraph* g_abs_ptr = reinterpret_cast<TGraph*>(g_ptr);
	g_abs_gd = TGraph(*g_abs_ptr);
	
	return get_ok;
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
	
	// current Gd absorption graph (ratio_abs_purev4_highv3.root) is taken from
	// EGADS since it seems to fit the data very well, and is expected to be Gd in very clean water.
	// however, for unknown reasons it is defined as purewater/gdloaded rather than vice versa.
	// convert it to the more sensible gdloaded/purewater (which should have a nice range of 0->1).
	// we also subtract off the baseline so that the range is from 0->-1, where 0 is no absorption.
	// Note that this baseline will be ~1 outside the Gd absorption peaks in EGADS
	// (since both purewater and gdloaded are measured in water), but will always be <1 in WCTE
	// (since the reference arm is fibre so has no water absorption)
	for(int i=0; i<g_absorption_ref.GetN(); ++i){
		g_absorption_ref.GetY()[i] = (1./g_absorption_ref.GetY()[i])-1.;
	}
	
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
	
	// We'll scale it, and add a background to account for contaminants.
	// during calibration a 3rd-order polynomial background was needed to extract clean Gd peaks,
	// but this could prove problematic when trying to fit both together....
	// it may be better to fit the background separately first, masking out the absorbance region
	// then subtract the background fit and fit the absorbance separately
	// (this is actually what was done to generate the calibration curves)
	std::string name="f_absfit_"+ledToAnalyse;
	const int n_absfit_pars = 1;
	
	// for reasons explained in MarcusAnalysis, the easiest way to make a functional fit
	// is to make a lambda function that captures a pointer to the fitted TGraph member
	TGraph* g_abs_ref_p = &g_absorption_ref;
	abs_fct = new TF1(name.c_str(),
		[g_abs_ref_p](double* x, double* par) -> double {
			/* old linear baseline
			// par [0] = y-scaling
			// par [1] = baseline offset (c)
			// par [2] = baseline gradient (m)
			double abs = par[0]*g_abs_ref_p->Eval(x[0]);
			double baseline = par[2]*(x[0]-276) + par[1];
			return (abs + baseline);
			*/
			double abs = par[0]*g_abs_ref_p->Eval(x[0]);
			/* -- yea looks like we need to split it
			// new, 3rd-order baseline
			// "[0]+[1]*(x-[2])+[3]*(x-[2])*(x-[2])+[4]*(x-[2])*(x-[2])*(x-[2])"
			double baseline = par[1]+
			                  par[2]*(x[0]-par[3])+
			                  par[4]*(x[0]-par[3])*(x[0]-par[3])+
			                  par[5]*(x[0]-par[3])*(x[0]-par[3])*(x[0]-par[3]);
			return (abs + baseline);
			*/
			return abs;
		},
		ROI_min, ROI_max, n_absfit_pars);
	
	// set initial parameters
	// TODO make these configuration parameters
	//std::vector<double> init_params{1,0,0};
	absfunc_init_params = std::vector<double>{0};
	abs_fct->SetParameters(absfunc_init_params.data());
	
	// a separate function for fitting baseline absorbance (absorbance of pure water + contaminants etc)
	name = "f_bgfit_"+ledToAnalyse;
	int n_bgfit_pars = 5;
	bg_abs_fct = new TF1(name.c_str(),"[0]+[1]*(x-[2])+[3]*(x-[2])*(x-[2])+[4]*(x-[2])*(x-[2])*(x-[2])", ROI_min, ROI_max);
	bgfunc_init_params = std::vector<double>{0.0603867,-0.00111141,265,0.000341765,-7.3633e-06};
	bg_abs_fct->SetParameters(bgfunc_init_params.data());
	
	// set parameter limits FIXME really should do this
	/*
	abs_fct->SetParLimits(0,0,20);         // y scaling
	abs_fct->SetParLimits(1,-0.2,1.2);     // y offset
	abs_fct->SetParLimits(2,-20,20);       // linear baseline gradient
	*/
	
	if(!abs_fct->IsValid()){
		throw std::runtime_error(m_unique_name+" GetAbsFunc failed to construct valid Gd TF1");
	}
	if(!bg_abs_fct->IsValid()){
		throw std::runtime_error(m_unique_name+" GetAbsFunc failed to construct valid background TF1");
	}
	Log(m_unique_name+" functional fit TF1 constructed",v_debug,verbosity);
	
	// note in CStore for recording in DB
	std::string formula_str = bg_abs_fct->GetExpFormula().Data(); // or bg_abs_fct->GetTitle();
	std::string params_str;
	for(size_t i=0; i<bg_abs_fct->GetNpar(); ++i){
		if(i!=0) params_str+=", ";
		params_str+=std::to_string(bg_abs_fct->GetParameter(i));
	}
	std::string absbgjson = "{ \"type\":\"configfile\",\"function\":\""+formula_str+"\",\"params\":\""+params_str+"\" }";
	m_data->CStore.Set("bgfitfunc_"+ledToAnalyse,absbgjson);
	
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
	// uh, for some reason it fails this test, even though it draws and evals just fine. :|
	/*
	if(!calib_curve.IsValid()){
		throw std::runtime_error(m_unique_name+" calibration curve retrieved from file '"
		                        +filename+"' is not valid!");
	}
	*/
	
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

bool ReturnOfTheMarcusAnalysisEpisode2::RemoveBackgroundAbsorbance(){
	// fit absorbance in UV excluding gd absorbance region with background poly, then subtract the fit across the whole UV ROI
	
	// initialise fit parameters (skip to carry over previous values)
	//bg_abs_fct->SetParameters(bgfunc_init_params.data());
	
	// make a TGraph with the gd absorbance region masked
	if(bg_indices.size()==0){
		for(size_t i=0; i<g_abs_gd.GetN(); ++i){
			if(g_abs_gd.GetX()[i]<269 || g_abs_gd.GetX()[i]>281) bg_indices.push_back(i);
			
		}
		g_abs_masked.Set(bg_indices.size());
	}
	for(size_t i=0; i<bg_indices.size(); ++i){
		g_abs_masked.SetPoint(i,g_abs_gd.GetX()[bg_indices.at(i)],g_abs_gd.GetY()[bg_indices.at(i)]);
	}
	
	// fit with background function
	bgfitresptr = g_abs_masked.Fit(bg_abs_fct,"RNMQS"); // or make a new one and call it 'tmp'
	//bgfitresptr = TFitResultPtr((TFitResult*)tmp->Clone());  // i don't know if Clone is required
	
	//g_abs_masked.SetName("g_abs_masked");
	//g_abs_masked.SaveAs("g_abs_masked.root");
	//bg_abs_fct->SaveAs("f_bg_fit.root");
	
	// record status of fit. probably redundant as none of these turned out to be reliable
	if(bgfitresptr->IsEmpty() || !bgfitresptr->IsValid() || bgfitresptr->Status()!=0){
		std::string fitstat;
		fitstat += " IsEmpty=" + std::to_string(bgfitresptr->IsEmpty());
		fitstat += " IsValid=" + std::to_string(bgfitresptr->IsValid());
		fitstat += " Status=" + std::to_string(bgfitresptr->Status());
		Log(m_unique_name+" Warning: abs bg fit status: "+fitstat,v_error,verbosity);
		// do not return false, these are not robust checks of a bad fit.
		// TODO implement a better check based on chi2
		// TODO based on past experience, do the fit multiple times
		//bgfit_success = false;
	}
	// we do it manually instead
	bgfit_success =  !bgfitresptr->IsEmpty() &&
	                  bgfitresptr->IsValid() &&
	                  bgfitresptr->Status()==0 &&
	                ((bgfitresptr->Chi2()/bgfitresptr->Ndf()) < 10.) &&
	    !TMath::IsNaN(bgfitresptr->GetParams()[0]) &&
	                 (bgfitresptr->GetErrors()[0] < 0.5);
	
	//  make a new background-subtracted TGraph of UV region
	if(g_bgfit.GetN()==0){
		g_bgfit.Set(g_abs_gd.GetN());
		bgfitvalues.resize(g_abs_gd.GetN());
		g_abs_bgrem.Set(g_abs_gd.GetN());
	}
	for (int i = 0; i < g_abs_gd.GetN(); ++i){
		double next_wl = g_abs_gd.GetX()[i];
		double next_bg = bg_abs_fct->Eval(next_wl);
		if(TMath::IsNaN(next_bg)){
			Log("abs fit function eval to NaN",v_error,verbosity);
			next_bg = 0;
		}
		bgfitvalues[i] = next_bg;
		g_bgfit.SetPoint(i, next_wl, next_bg);
		g_abs_bgrem.SetPoint(i, next_wl, g_abs_gd.GetY()[i] - next_bg);
		// we could 'reconstruct' the GAD arm data without Gd at this point
		//g_gad_fit.SetPoint(i,next_wl,next_bg*g_ref.GetY()[i]); // g_ref: reference arm graph
		// .... but there's no real reason to?
	}
	
	//g_abs_bgrem.SaveAs("g_abs_bgrem.root");
	
	if(!save_trees) return bgfit_success;
	
	// maybe save some stuff here for debug, or don't
	if(!outtree->GetBranch("bgfit")) outtree->Branch("bgfit",&bgfitvaluesp);
	else outtree->SetBranchAddress("bgfit",&bgfitvaluesp);
	
	return bgfit_success;
	
}

// -------------------------------------------------------------------------//


bool ReturnOfTheMarcusAnalysisEpisode2::FitAbsorbance(bool bgrem){
	// fit absorbance in UV region with reference shape and extract the scaling required
	
	// initialise fit parameters (skip to carry over previous values)
	//abs_fct->SetParameters(absfunc_init_params.data());
	
	TGraph& g_tofit = bgrem ? g_abs_bgrem : g_abs_gd;
	absfitresptr = g_tofit.Fit(abs_fct,"RNMQS"); // or make a new one and call it 'tmp'
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
	if(g_absfit.GetN()==0){
		g_absfit.Set(g_abs_gd.GetN());
		absfitvalues.resize(g_abs_gd.GetN());
	}
	for (int i = 0; i < g_abs_gd.GetN(); ++i){
		double next_wl = g_abs_gd.GetX()[i];
		double next_abs = abs_fct->Eval(next_wl);
		if(TMath::IsNaN(next_abs)){
			Log("abs fit function eval to NaN",v_error,verbosity);
			next_abs = 0;
		}
		absfitvalues[i] = next_abs;
		g_absfit.SetPoint(i, next_wl, next_abs);
	}
	
	if(!save_trees) return absfit_success;
	
	if(!outtree->GetBranch("absfit")) outtree->Branch("absfit",&absfitvaluesp);
	else outtree->SetBranchAddress("absfit",&absfitvaluesp);
	
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
	
	if(!outtree->GetBranch("metric")) outtree->Branch("metric",&metric);
	else outtree->SetBranchAddress("metric",&metric);
	if(!outtree->GetBranch("gd_conc")) outtree->Branch("gd_conc",&gd_conc);
	else outtree->SetBranchAddress("gd_conc",&gd_conc);
	
	return get_ok;
}

// -------------------------------------------------------------------------//

void ReturnOfTheMarcusAnalysisEpisode2::ReInit(){
	
	// Remove outputs from previous Executions, so we don't carry over results
	// from a previous fit if we bail early
	m_data->CStore.Remove("NewMarcusAnalyseEp2");
	
	m_data->CStore.Remove("absfit");
	m_data->CStore.Remove("absfitresptr");
	m_data->CStore.Remove("absfit_success");
	m_data->CStore.Remove("bgfit_success");
	
	m_data->CStore.Remove("metric_and_err");
	m_data->CStore.Remove("conc_and_err");
	
	return;
}

void ReturnOfTheMarcusAnalysisEpisode2::UpdateDataModel(){
	
	// A flag informing downstream tools that new results from this Tool are available
	m_data->CStore.Set("NewMarcusAnalyseEp2",ledToAnalyse);
	
	// for website
	m_data->CStore.Set("bgfit",reinterpret_cast<intptr_t>(&g_bgfit));
	m_data->CStore.Set("bg_rem_abs",reinterpret_cast<intptr_t>(&g_abs_bgrem));
	m_data->CStore.Set("absfit",reinterpret_cast<intptr_t>(&g_absfit));
	
	// results for DB
	m_data->CStore.Set("bgfitresptr", reinterpret_cast<intptr_t>(&bgfitresptr));
	m_data->CStore.Set("absfitresptr", reinterpret_cast<intptr_t>(&absfitresptr));
	m_data->CStore.Set("bgfit_success",bgfit_success);
	m_data->CStore.Set("absfit_success",absfit_success);
	m_data->CStore.Set("metric_and_err",metric_and_err);
	m_data->CStore.Set("conc_and_err",conc_and_err);
	
	// TODO
	//m_data->CStore.Set("gad_fitted_max",gad_fitted_max);
	// led intensity down gad arm, obtained by fitting corrected ref to gad data in sidebands
	// (this could be used to measure e.g. solarization?)
	
	// send stuff to WCTE DB
	m_data->monitoring_store.Set("gd_conc",conc_and_err.first);
	
	std::string graph_json, graph_name;
	graph_name = "g_abs_bgrem";
	graph_json = TBufferJSON::ToJSON(&g_abs_bgrem).Data();
	if(m_data->services) m_data->services->SendROOTplot(graph_name, "AL", graph_json, true);
	
	graph_name = "g_absfit";
	graph_json = TBufferJSON::ToJSON(&g_absfit).Data();
	if(m_data->services) m_data->services->SendROOTplot(graph_name, "AL", graph_json, true);
	
	graph_name = "g_bgfit";
	graph_json = TBufferJSON::ToJSON(&g_bgfit).Data();
	if(m_data->services) m_data->services->SendROOTplot(graph_name, "AL", graph_json, true);
	
	return;
	
}
