#ifndef ReturnOfTheMarcusAnalysis_H
#define ReturnOfTheMarcusAnalysis_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "DataModel.h"

struct measurement {
	
	int led_entry;                      // entry number in LED-on tree for gad arm measurement
	int ref_entry;                      // entry number in LED-on tree for ref arm measurement
	int ref_dark_entry;                 // entry number in Dark tree for ref arm measurement
	int gad_dark_entry;                 // entry number in Dark tree for gad arm measurement
	
	std::vector<double> values_ref;     // LED-on intensity values for reference arm measurement
	std::vector<double> dark_ref;       // Dark intensity values for reference arm measurement
	std::vector<double> values_gad;     // LED-on intensity values for gad arm measurement
	std::vector<double> dark_gad;       // Dark intensity values for gad arm measurement
	
	std::vector<double> darksub_ref;    // dark-subtracted reference arm intensities
	std::vector<double> darksub_gad;    // dark-subtracted gad arm intensities
	
	std::vector<double> corrected_ref;  // reference arm intensities corrected for pure water absorption
	std::vector<double> ratio_abs;      // ratio of gad arm intensity to corrected reference arm intensity
	std::vector<double> ratio_err;      // based on 1/sqrt(X) or something... TBD
	
	Short_t LED_yr, LED_mon, LED_dy, LED_hr, LED_mn, LED_sc;    // time the LED-on measurement was taken
	
};

class ReturnOfTheMarcusAnalysis: public Tool {
	
	public:
	ReturnOfTheMarcusAnalysis();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	
	private:
	
	// Initialise:
	bool GetPureWaterTransparency();
	bool GetPureWaterTransparency(int pureref_ver);
	bool GetPureWaterTransparency(std::string filename);
	
	// Execute:
	bool ReadyToAnalyse();
	void ReInit();
	void SetGraphTitles();
	bool GetTrees();
	bool ReadBranch(TTree* tree, const std::string& branch, const size_t entry, std::vector<double>* values);
	bool ReadValues();
	bool CalculateAbsorbance();
	void UpdateDataModel();
	
	std::string ledToAnalyse;
	
	// option to use first measurement to define water transparency
	bool make_pureref = false;
	bool GeneratePureWaterTransparency();
	
	// filled during initialise
	TGraph g_pure_absorbance;    // absorbance of pure water (ideally normalised)
	
	// filled in ReadValues
	TGraph g_ref;                // reference arm data
	TGraph g_gad;                // gad arm data
	
	// filled in CalculateAbsorbance
	TGraph g_ref_corr;           // ref arm corrected for water transparency
	TGraph g_abs;                // log(ratio) of corrected ref arm to gad arm
	
	// input trees
	TTree* led_tree = nullptr;
	TTree* dark_tree = nullptr;
	
	// for getting data from input trees
	std::vector<double> gad_values, ref_values, gad_dark, ref_dark, wavelengths;
	
	std::vector<double>* gad_valuesp= nullptr;
	std::vector<double>* ref_valuesp = nullptr;
	std::vector<double>* gad_darkp = nullptr;
	std::vector<double>* ref_darkp = nullptr;
	std::vector<double>* wavelengthsp = nullptr;
	
	// for random stats tracking
	double dark_mean, dark_sigma;
	double ref_max, corrected_ref_max, ref_min;
	double gad_max, gad_min;
	
	// for ROI
	bool init=true;
	double min_wl=200;
	double max_wl=800;
	int min_wl_index=-1;
	int max_wl_index=-1;
	
	// to save traces to output file (for debug, for now?)
	bool save_trees=false;
	TFile* outfile=nullptr; // used when putting debug in a separate file
	TTree* outtree=nullptr;
	// branches
	std::vector<double> ref_corr_values, absorbances;
	std::vector<double>* ref_corr_valuesp=nullptr;
	std::vector<double>* absorbancesp=nullptr;
	
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
