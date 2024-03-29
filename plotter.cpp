#include "plotter.h"

int Plotter::LoadFile(std::string file, std::string trace){
	
	rootfile = TFile::Open(file.c_str(),"READ");
	if(rootfile==nullptr || rootfile->IsZombie()){
		std::cerr<<"no such file "<<file<<std::endl;
		return 0;
	}
	
	ledtree = (TTree*)rootfile->Get(trace.c_str());
	if(ledtree==nullptr || ledtree->GetEntries()==0){
		std::cerr<<"no led tree "<<trace<<" or it has no entries"<<std::endl;
		return 0;
	}
	darktree = (TTree*)rootfile->Get("Dark");
	if(darktree==nullptr || darktree->GetEntries()==0){
		std::cerr<<"no Dark tree or it has no entries, skipping dark subtraction"<<std::endl;
	}
	
	// get number of measurements
	int num_meas = ledtree->GetEntries();
	if(verbosity) std::cout<<"led tree had "<<num_meas<<" measurements"<<std::endl;
	
	if(verbosity) std::cout<<"setting led tree branch addresses"<<std::endl;
	ledtree->SetBranchAddress("value", &values_p);
	ledtree->SetBranchAddress("error", &errors_p);
	ledtree->SetBranchAddress("year",&yr);
	ledtree->SetBranchAddress("month",&mon);
	ledtree->SetBranchAddress("day",&dy);
	ledtree->SetBranchAddress("hour",&hr);
	ledtree->SetBranchAddress("min",&mn);
	ledtree->SetBranchAddress("sec",&sc);
	
	if(verbosity) std::cout<<"setting dark tree branch addresses"<<std::endl;
	darktree->SetBranchAddress("value", &values_dark_p);
	darktree->SetBranchAddress("wavelength", &wavelengths_p);
	darktree->SetBranchAddress("error", &errors_dark_p);
	darktree->SetBranchAddress("year",&yr);
	darktree->SetBranchAddress("month",&mon);
	darktree->SetBranchAddress("day",&dy);
	darktree->SetBranchAddress("hour",&hr);
	darktree->SetBranchAddress("min",&mn);
	darktree->SetBranchAddress("sec",&sc);
	
	// get num datapoints
	if(verbosity) std::cout<<"getting n datapoints"<<std::endl;
	darktree->GetEntry(0);
	n_datapoints = values_dark.size();
	
	// fill wavevelength errors using bin widths
	wavelength_errors.resize(n_datapoints);
	std::fill(wavelength_errors.begin(), wavelength_errors.end(), wavelengths.at(1)-wavelengths.at(0));
	
	return num_meas;
}

TMultiGraph* Plotter::MakePlot(int entryi){
	
	// plot entry 0 from led.
	ledtree->GetEntry(entryi);
	
	// Look up associated closest dark from dark tree.
	int dark_entry = GetNextDarkEntry(entryi);
	darktree->GetEntry(dark_entry);
	
	// do dark subtraction
	for(int j=0; j<n_datapoints; ++j){
		values.at(j) -= values_dark.at(j);
	}
	
		// split into in- and -out of absorbance bands
		std::vector<double> inband_values;
		std::vector<double> inband_wls;
		std::vector<double> sideband_values;
		std::vector<double> sideband_wls;
		std::vector<double> other_values;
		std::vector<double> other_wls;
		
		std::vector<double> inband_errs;
		std::vector<double> sideband_errs;
		
		//std::cout<<"values: [";
		for(int j=0; j<n_datapoints; ++j){
			//if(j>0) std::cout<<", ";
			//std::cout<<values.at(j);
			if(wavelengths.at(j)>270 && wavelengths.at(j)<280){
				// in band
				inband_values.push_back(values.at(j));
				inband_wls.push_back(wavelengths.at(j));
				inband_errs.push_back(errors.at(j));
			} else if(wavelengths.at(j)>260 && wavelengths.at(j)<300){
				// side band (lobes)
				sideband_values.push_back(values.at(j));
				sideband_wls.push_back(wavelengths.at(j));
				sideband_errs.push_back(errors.at(j));
			} else {
				other_values.push_back(values.at(j));
				other_wls.push_back(wavelengths.at(j));
			}
		}
		//std::cout<<"]"<<std::endl;
		
		// make TGraphErrors from data
		g_all = new TGraphErrors(n_datapoints, wavelengths.data(), values.data(), wavelength_errors.data(), errors.data());
		TGraph* g_inband = new TGraph(inband_values.size(), inband_wls.data(), inband_values.data());
		TGraph* g_sideband = new TGraph(sideband_values.size(), sideband_wls.data(), sideband_values.data());
		TGraph* g_other = new TGraph(other_values.size(), other_wls.data(), other_values.data());
		
		g_inband->SetLineColor(kRed);
		g_sideband->SetLineColor(kBlue);
		g_other->SetLineColor(kSpring+5);
		
		// PRINT OUT JUST TEMP
		/*
		std::cerr<<"{\"xerrs\": [";
		for(int i=0; i<inband_values.size(); ++i){
			std::cerr<<"0, ";
		}
		std::cerr<<"], \"xvals\": [";
		for(int i=0; i<inband_wls.size(); ++i){
			std::cerr<<inband_wls.at(i)<<", ";
		}
		std::cerr<<"], \"yerrs\": [";
		for(int i=0; i<inband_errs.size(); ++i){
			std::cerr<<inband_errs.at(i)<<", ";
		}
		std::cerr<<"], \"yvals\": [";
		for(int i=0; i<inband_values.size(); ++i){
			std::cerr<<inband_values.at(i)<<", ";
		}
		std::cerr<<"]}"<<std::endl;
		*/
		
		/*
		std::cerr<<"{\"xerrs\": [";
		for(int i=0; i<sideband_values.size(); ++i){
			std::cerr<<"0, ";
		}
		std::cerr<<"], \"xvals\": [";
		for(int i=0; i<sideband_wls.size(); ++i){
			std::cerr<<sideband_wls.at(i)<<", ";
		}
		std::cerr<<"], \"yerrs\": [";
		for(int i=0; i<sideband_errs.size(); ++i){
			std::cerr<<sideband_errs.at(i)<<", ";
		}
		std::cerr<<"], \"yvals\": [";
		for(int i=0; i<sideband_values.size(); ++i){
			std::cerr<<sideband_values.at(i)<<", ";
		}
		std::cerr<<"]}"<<std::endl;
		*/
	
	TMultiGraph* mg_all = new TMultiGraph("mg_all","mg_all");
	mg_all->Add(g_inband);
	mg_all->Add(g_sideband);
	mg_all->Add(g_other);
	
	if(c1){
		c1->cd();
		mg_all->Draw("AL");
	}
	
	return mg_all;
}


int Plotter::GetNextDarkEntry(int ledon_entry_num){
	
	// get led timestamp
	ledtree->GetEntry(ledon_entry_num);
	struct tm ledtime;
	ledtime.tm_year = yr - 1900;
	ledtime.tm_mon = mon - 1;
	ledtime.tm_mday = dy;
	ledtime.tm_hour = hr;
	ledtime.tm_min = mn;
	ledtime.tm_sec = sc;
	time_t ledtime_t = mktime(&ledtime);
	
	// loop over dark entries
	int darkentry=-1;
	for(int i=0; i<darktree->GetEntries(); ++i){
		darktree->GetEntry(i);
		
		struct tm darktime;
		darktime.tm_year = yr - 1900;
		darktime.tm_mon = mon - 1;
		darktime.tm_mday = dy;
		darktime.tm_hour = hr;
		darktime.tm_min = mn;
		darktime.tm_sec = sc;
		time_t darktime_t = mktime(&darktime);
		
		// difftime does [ time_a - time_b ]
		double numsecs = difftime(ledtime_t, darktime_t);
		// do we use the closest in time dark, or the last dark before the led-on?
		// what if someone decides to do the dark measurements after the led?
		if(numsecs<0 && darkentry>=0) break;
		darkentry=i;
	}
	
	return darkentry;
}

