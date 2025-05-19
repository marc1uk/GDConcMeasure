#ifndef ReturnOfTheMarcusAnalysisEpisode2_H
#define ReturnOfTheMarcusAnalysisEpisode2_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "DataModel.h"
// &g_ref,  &g_gad,  &g_ref_corr,  &g_gadfit,  &g_abs
namespace {
	const double ROI_min = 260; // nm
	const double ROI_max = 300; // nm
	
}

class ReturnOfTheMarcusAnalysisEpisode2: public Tool {
	
	public:
	ReturnOfTheMarcusAnalysisEpisode2();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	
	private:
	
	// Initialise:
	bool GetCalibrationCurve();
	bool GetCalibrationCurveFromConfigs();
	bool GetCalibrationCurveFromFile();
	bool GetCalibrationCurveFromDB();
	
	bool GetAbsorptionRef();
	bool GetAbsorptionRef(int absref_ver);
	bool GetAbsorptionRef(std::string filename);
	bool GetAbsFunc();
	
	std::string ledToAnalyse;
	TGraph g_absorption_ref;     // gd absorption reference shape, scaled to fit data
	TGraph g_abs_gd;             // gd region only (from ReturnOfTheMarcusAnalysis Tool)
	TF1* bg_abs_fct;
	std::vector<double> bgfunc_init_params;
	TF1* abs_fct;
	std::vector<double> absfunc_init_params; // TODO read from config, add param limits
	TF1 calib_curve;
	
	// Execute:
	bool ReadyToAnalyse();
	void ReInit();
	void SetGraphTitles();
	bool GetAbsorbance();
	bool RemoveBackgroundAbsorbance();
	bool FitAbsorbance(bool bgrem);
	bool CalculateConcentration();
	void UpdateDataModel();
	
	// filled in RemoveBackgroundAbsorbance
	TGraph g_bgfit;
	TGraph g_abs_masked; // absorbance in UV region but with gd region masked (used to fit background)
	std::vector<size_t> bg_indices;
	TGraph g_abs_bgrem;
	TFitResultPtr bgfitresptr;
	bool bgfit_success = false; // our own metric as we can't trust the status of TFitResultPtr
	
	// filled in FitAbsorbance
	TGraph g_absfit;
	TFitResultPtr absfitresptr;
	bool absfit_success = false; // our own metric as we can't trust the status of TFitResultPtr
	double gad_fitted_max;
	
	// filled in CalculateConcentration
	double metric, gd_conc;
	std::pair<double,double> metric_and_err;
	std::pair<double,double> conc_and_err;
	
	// to save traces to output file (for debug, for now?)
	bool save_trees=true;   // let it follow ReturnOfTheMarcusAnalysis
	TTree* outtree=nullptr; // made in ReturnOfTheMarcusAnalysis
	// we add one new branch storing absorbance fit curve
	std::vector<double> bgfitvalues;
	std::vector<double>* bgfitvaluesp;
	std::vector<double> absfitvalues;
	std::vector<double>* absfitvaluesp=nullptr;
	
	std::string m_configfile;
	// for logging
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	int get_ok;
	
};


#endif
