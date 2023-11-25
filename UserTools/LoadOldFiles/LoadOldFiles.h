#ifndef LoadOldFiles_H
#define LoadOldFiles_H

#include <string>
#include <iostream>
#include <map>

#include "Tool.h"

struct oldfile{
	std::string filename;
	std::string treename;
	int runnum;
};

class LoadOldFiles: public Tool {


 public:

  LoadOldFiles();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

 private:
  int ParseFileList(std::string file_list, int max_files=-1);
  
  std::vector<oldfile> files;
  std::vector<oldfile>::iterator next_file_iter;
  TTree* darkTree=nullptr;
  TTree* darkTreeNew=nullptr;
  TTree* ledTree=nullptr;
  TFile* nextfile=nullptr;
  TFile* tmpfile=nullptr;
  int measurementnum=0;
  
  int get_ok = 0;
  int verbosity=1;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;

};


#endif
