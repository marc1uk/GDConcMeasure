#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"

int main(int argc, char* argv[]){
	
	if(argc<3){
		std::cout<<"usage: "<<argv[0]<<" [ledname] [file] [date] [formula]"<<std::endl;
		std::cout<<"[file] should contain a TTree called 'params' with 3 branches: 'raw','simple' and 'complex',\n"
		         <<"\teach holding a vector<double> of parameter values"<<std::endl;
		std::cout<<"[date] should be in format e.g. '2023-08-12 12:03:02';\n"
		         <<"\tif not given 'now()' will be used"<<std::endl;
		std::cout<<"[formula] should be an expression that may be used to define a TF1\n"
		         <<"\tdefining the curve of concentration vs metric;\n"
		         <<"\tif not given 'pol6' will be used"<<std::endl;
		return 0;
	}
	
	std::string ledname = argv[1];
	std::string filename = argv[2];
	std::string caldate="now()";
	if(argc>3){
		caldate=argv[3];
	}
	std::string formula="pol6";
	if(argc>4){
		formula=argv[4];
	}
	
	TFile* f = TFile::Open(filename.c_str(),"READ");
	if(f==nullptr || f->IsZombie()) return 1;
	TTree* t = (TTree*)f->Get("params");
	if(t==nullptr) return 2;
	std::vector<double> raw_params;
	std::vector<double> simple_params;
	std::vector<double> complex_params;
	std::vector<double>* raw_paramsp=&raw_params;
	std::vector<double>* simple_paramsp=&simple_params;
	std::vector<double>* complex_paramsp=&complex_params;
	int ok=0;
	ok += t->SetBranchAddress("raw",&raw_paramsp);
	ok += t->SetBranchAddress("simple",&simple_paramsp);
	ok += t->SetBranchAddress("complex",&complex_paramsp);
	if(ok!=0){
		std::cerr<<"failed to set branch addresses"<<std::endl;
		return 3;
	}
	ok = t->GetEntry(0);
	if(ok<=0){
		std::cout<<"failed to get entry 0"<<std::endl;
		return 4;
	}
	
	// this script just generates the psql command, it doesn't run it.
	std::cerr<<"# This script generates the appropriate psql commands, but does not run them.\n"
	         <<"# Please verify the following commands, and if acceptable, run them.\n"<<std::endl;
	
	std::string psqlstring = "psql -c \"INSERT INTO data (timestamp, name, tool, ledname, values ) "
	                           "VALUES ( '"+caldate+"', 'calibration_curve', 'MarcusAnalysis', '"+ledname+"', "
	                           "'{\\\"formula\\\":\\\"" + formula + "\\\", \\\"method\\\":\\\"raw\\\", "
	                           "\\\"params\\\": [";
	for(int i=0; i<raw_params.size(); ++i){
	        if(i>0) psqlstring += ",";
	        psqlstring += std::to_string(raw_params.at(i));
	}
	psqlstring += "], \\\"version\\\": ";
	
	// we can run a psql command to get the next version number
	std::string cmd = "psql -At -c \"SELECT MAX((values->>'version')::integer)+1 FROM data "
	                  "WHERE name='calibration_curve' AND ledname='"+ledname+"' "
	                  "AND values->>'method'='raw'\" | tr -d '\n' ";
	
	// we need to use a system call to insert it into our query
	std::cout<<psqlstring<<std::flush;
	system(cmd.c_str());
	std::cout<<"}') \";"<<std::endl;
	
	// repeat it all for the other methods
	psqlstring = "psql -c \"INSERT INTO data (timestamp, name, tool, ledname, values ) "
	                           "VALUES ( '"+caldate+"', 'calibration_curve', 'MarcusAnalysis', '"+ledname+"', "
	                           "'{\\\"formula\\\":\\\"" + formula + "\\\", \\\"method\\\":\\\"simple\\\", "
	                           "\\\"params\\\": [";
	for(int i=0; i<simple_params.size(); ++i){
	        if(i>0) psqlstring += ",";
	        psqlstring += std::to_string(simple_params.at(i));
	}
	psqlstring += "], \\\"version\\\": ";
	cmd = "psql -At -c \"SELECT MAX((values->>'version')::integer)+1 FROM data "
	      "WHERE name='calibration_curve' AND ledname='"+ledname+"' AND values->>'method'='simple'\" | tr -d '\n' ";
	std::cout<<psqlstring<<std::flush;
	system(cmd.c_str());
	std::cout<<"}') \";"<<std::endl;
	
	// repeat it all for the other methods
	psqlstring = "psql -c \"INSERT INTO data (timestamp, name, tool, ledname, values ) "
	                           "VALUES ( '"+caldate+"', 'calibration_curve', 'MarcusAnalysis', '"+ledname+"', "
	                           "'{\\\"formula\\\":\\\"" + formula + "\\\", \\\"method\\\":\\\"complex\\\", "
	                           "\\\"params\\\": [";
	for(int i=0; i<complex_params.size(); ++i){
	        if(i>0) psqlstring += ",";
	        psqlstring += std::to_string(complex_params.at(i));
	}
	psqlstring += "], \\\"version\\\": ";
	cmd = "psql -At -c \"SELECT MAX((values->>'version')::integer)+1 FROM data "
	      "WHERE name='calibration_curve' AND ledname='"+ledname+"' AND values->>'method'='complex'\" | tr -d '\n' ";
	std::cout<<psqlstring<<std::flush;
	system(cmd.c_str());
	std::cout<<"}') \";"<<std::endl;
	
	return 0;
}
