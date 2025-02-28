#ifndef MarcusAnalysis_H
#define MarcusAnalysis_H

#include <string>
#include <iostream>
#include <sstream>

#include "Tool.h"
#include "DataModel.h"

class MarcusAnalysis: public Tool {
	
	public:
	
	MarcusAnalysis();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	private:
	
	bool GetPureRefDB(int pureref_ver);
	bool GetPureTrace();
	bool GetPureRef(std::string filename);
	static double PureScaledPlusExtras(double* x, double* par);
	bool GetPureScaledPlusExtras();
	bool GetCalibrations();
	bool GetCalCurve(std::string name, int method_i);
	bool GetCalCurveDB(std::string method, std::string calibID);
	bool GetFourGaussianFunc();
	
	void ReInit();
	bool ReInitFourGaussians();
	bool ReadyToAnalyse();
	bool FindDarkAndLEDTrees();
	bool GetDarkSubTrace();
	bool FitPure();
	bool CalculateAbsorbance();
	bool CalculateConcentration(std::string method);
	bool CalculatePeakHeights(std::string method);
	bool PeakDiffToConcentration(std::string method);
	void UpdateDataModel();
	
	// CalculatePeakHeights for simple method
	bool FitTwoGaussians(std::pair<double,double>& peak_heights, std::pair<double,double>& peak_errs, std::pair<double,double>& peak_posns, std::vector<TFitResultPtr>& fitresptrs);
	// CalculatePeakHeights for complex method
	
	bool FitFourGaussians(std::pair<double,double>& peak_heights, std::pair<double,double>& peak_errs, std::pair<double,double>& peak_posns, std::vector<TFitResultPtr>& fitresptrs);
	
	// calculating errors on peak heights for complex method is a tricky business
	// due to correlations between the fits
	std::pair<double,double> CalculateError(double peak1_pos, double peak2_pos);
	
	std::string ledToAnalyse="";
	
	// input data source
	std::pair<TTree*, TTree*> light_and_dark_trees{nullptr,nullptr};
	
	// input data
	TGraphErrors g_inband;
	TGraphErrors g_sideband;
	TGraphErrors g_other;
	TGraph g_sideband_noerrs;
	
	// pure reference used as part of fit, and saved to database
	TGraph dark_subtracted_pure;       // g_pureref_ledname
	
	// values that go into the db webpage table
	double dark_mean;
	double dark_sigma;
	double led_on_max;
	double led_on_min;
	
	// function built from pure reference, fit to data in sidebands
	TF1* pure_fct = nullptr;           // f_purefit_ledname
	std::vector<double> pure_init_params;
	
	// TGraph of fit result in over sideband+absorbance region,
	// only for convenience to store into database for webpage
	TGraphErrors g_purefit;              // g_purefit_ledname
	
	// pure fit result status, parameters and errors
	TFitResultPtr purefitresptr;
	
	// abosorption calculated from pure_fct and dark_subtracted_data
	TGraphErrors g_absorb;                  // g_abs_ledname
	
	// used by 'raw' fitting method - the data sample from the inband array
	// at which the 273 and 276 nm absorption peaks are centered
	int sample_273;
	int sample_276;
	
	// function fit to absorption data in complex method
	TF1* fourgaussianfunc = nullptr;
	
	// functions for absorbance to concentration conversion for each method
	std::map<std::string,TF1> calib_curves;
	
	// results of absorbance fits and conversion to concentration, for each method
	std::map<std::string, BoostStore> results;
	
	// also need to retain a map of TFitResultPtrs for each method
	std::map<std::string, std::vector<TFitResultPtr>> resultsmap;
	
	// for logging
	std::stringstream logmessage;
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	int get_ok;
	
};


#endif
