#include <iostream>
#include <vector>
#include "TFile.h"
#include "TF1.h"

int main(int argc, char* argv[]){
	
	if(argc<3){
		std::cout<<"usage: "<<argv[0]<<" [file] [ledname] [date]"<<std::endl;
		std::cout<<"[file] should contain a TF1 with name corresponding to the led\n";
		std::cout<<"[date] should be in format e.g. '2023-08-12 12:03:02';\n"
		         <<"\tif not given 'now()' will be used"<<std::endl;
		return 0;
	}
	
	std::string filename = argv[1];
	std::string ledname = argv[2];
	std::string caldate="now()";
	if(argc>3){
		caldate=argv[3];
	}
	
	TFile* f = TFile::Open(filename.c_str(),"READ");
	if(f==nullptr || f->IsZombie()) return 1;
	TF1* t = dynamic_cast<TF1*>(f->Get(ledname.c_str()));
	if(t==nullptr) return 2;
	std::string formula=t->GetTitle(); // stores short formula, e.g. 'pol6'
	// but could technically be overridden by user...
	// could maybe use 'GetExpFormula' which stores the expanded version
	// ('p[0]+p[1]*x[1]+...')? might be more robust, but less user friendly
	
	std::vector<double> params(t->GetParameters(),t->GetParameters()+t->GetNpar());
	int ok=0;
	
	// this script just generates the psql command, it doesn't run it.
	std::cerr<<"# This script generates the appropriate psql commands, but does not run them.\n"
	         <<"# Please verify the following commands, and if acceptable, run them.\n"<<std::endl;
	
	std::string psqlstring = "psql -c \"INSERT INTO data (timestamp, name, tool, ledname, values ) "
	                           "VALUES ( '"+caldate+"', 'calibration_curve', 'MatthewAnalysisStrikesBack', '"+ledname+"', "
	                           "'{\\\"formula\\\":\\\"" + formula + "\\\", "
	                           "\\\"params\\\": [";
	for(int i=0; i<params.size(); ++i){
	        if(i>0) psqlstring += ",";
	        psqlstring += std::to_string(params.at(i));
	}
	psqlstring += "], \\\"version\\\": ";
	
	// we can run a psql command to get the next version number
	std::string cmd = "psql -At -c \"SELECT MAX((values->>'version')::integer)+1 FROM data "
	                  "WHERE name='calibration_curve' AND ledname='"+ledname+"' \""
	                  " | tr -d '\n' ";
	
	// we need to use a system call to insert it into our query
	std::cout<<psqlstring<<std::flush;
	system(cmd.c_str());
	std::cout<<"}') \";"<<std::endl;
	
	return 0;
}
