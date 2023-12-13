#ifndef plotter_h
#define plotter_h
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1.h"

#include <iostream>
#include <chrono>
#include <thread>


class Plotter{
	public:
	Plotter(TCanvas* canv, int verb=1) : c1(canv), verbosity(verb) {};
	
	int verbosity=1;
	TFile* rootfile=nullptr;
	TTree* ledtree=nullptr;
	TTree* darktree=nullptr;
	
	std::vector<double> values;
	std::vector<double> wavelengths;
	std::vector<double> errors;
	std::vector<double> values_dark;
	std::vector<double> errors_dark;
	std::vector<double> wavelength_errors;  // make this up based on 1/2 bin width.
	Short_t yr, mon, dy, hr, mn, sc;
	
	// think these need to remain in scope, for some reason
	std::vector<double>* values_p = &values;
	std::vector<double>* errors_p = &errors;
	std::vector<double>* values_dark_p = &values_dark;
	std::vector<double>* errors_dark_p = &errors_dark;
	std::vector<double>* wavelengths_p = &wavelengths;
	
	int n_datapoints;
	TGraphErrors* g_all = nullptr;
	TCanvas* c1=nullptr;
	
	int LoadFile(std::string file, std::string trace);
	int GetNextDarkEntry(int ledon_entry_num);
	TGraphErrors* GetSimpleGraph(int entryi=0);
	TMultiGraph* MakePlot(int entryi=0);
	
};

#endif
