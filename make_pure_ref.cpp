#include "TFile.h"
#include "TGraphErrors.h"
#include "plotter.h"
#include <iostream>

int main(int argc, const char** argv){
	if(argc<4){
		std::cerr<<"usage: "<<argv[0]<<" <file> <ledname> <pure/highconc> <date>"<<std::endl;
		return 0;
	}
	
	std::string fname = argv[1];
	std::string ledname = argv[2];
	std::string type = argv[3];
	if(type!="pure" && type!="highconc"){
		std::cerr<<"usage: "<<argv[0]<<" <file> <ledname> <pure/highconc> <date>"<<std::endl;
		return 0;
	}
	type+= "_curve";
	std::string date  = (argc>4) ? argv[4] : "now()";
	
	TFile* f = TFile::Open(fname.c_str());
	if(f==nullptr){
		std::cerr<<"couldn't open file "<<fname<<std::endl;
		return 1;
	}
	
	// see if it contains a graph, if so use that.
	TGraph* g = (TGraph*)f->Get("Graph");
	if(g==nullptr){
		std::cerr<<"no graph named 'Graph' in file "<<fname<<std::endl;
		// probably a data file
		Plotter myplotter(nullptr, 0);
		bool ok =  myplotter.LoadFile(fname, ledname);
		if(!ok){
			std::cerr<<"couldn't load trees from file "<<fname<<std::endl;
			return 2;
		}
		myplotter.MakePlot();
		g = myplotter.g_all;  // n.b. i guess we own this, plotter class leaks it ^_^;
		if(g==nullptr || g->GetN()==0){
			std::cerr<<"couldn't make data graph"<<std::endl;
			return 2;
		}
	}
	std::cerr<<"graph has "<<g->GetN()<<" points"<<std::endl;
	if(dynamic_cast<TGraphErrors*>(g)!=nullptr) std::cerr<<"graph has errors"<<std::endl;
	
	std::cerr<<"# This script only generates the psql command to run; it does not make the database entry.\n"
	         <<"# Please validate the below printout and run the command if everything looks good\n"<<std::endl;
	
	std::string psqlstring = "psql -c \"INSERT INTO data (timestamp, tool, ledname, name, values) VALUES "
	                         "( '"+date+"', 'MarcusAnalysis', '" +ledname + "', '"+type+"', '{ \\\"xvals\\\":[";
	for(int i=0; i<g->GetN(); ++i){
		psqlstring += std::to_string(g->GetX()[i]) + ",";
	}
	psqlstring.pop_back(); // remove trailing ','
	psqlstring += "], \\\"yvals\\\":[";
	for(int i=0; i<g->GetN(); ++i){
		psqlstring += std::to_string(g->GetY()[i]) + ",";
	}
	psqlstring.pop_back(); // remove trailing ','
	psqlstring += "]";
	
	TGraphErrors* ge = dynamic_cast<TGraphErrors*>(g);
	if(ge!=nullptr){
		// add errors if we have them
		psqlstring += ", \\\"xerrs\\\":[";
		for(int i=0; i<g->GetN(); ++i){
			psqlstring += std::to_string(g->GetEX()[i]) + ",";
		}
		psqlstring.pop_back(); // remove trailing ','
		psqlstring += "], \\\"yerrs\\\":[";
		for(int i=0; i<g->GetN(); ++i){
			psqlstring += std::to_string(g->GetEY()[i]) + ",";
		}
		psqlstring.pop_back(); // remove trailing ','
		psqlstring += "]";
	}
	psqlstring += ", \\\"version\\\":";   //+ version
	
	// inline next version number by querying for it
	std::string cmd = "set -o pipefail; psql -At -c \"SELECT MAX((values->>'version')::integer)+1 FROM data WHERE name='"+type+"' "
	                  "AND ledname='" + ledname + "'\" | grep -E '[0-9]+' | tr -d '\n'";
	std::cout<<psqlstring<<std::flush;
	int ret = system(cmd.c_str());
	if(ret!=0) std::cout<<"0"<<std::flush;  // in case no db entry
	std::cout<<"}' );\" "<<std::endl;
	
	return 0;
}
