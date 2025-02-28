#ifndef TraceAverage_H
#define TraceAverage_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include "Tool.h"
#include "DataModel.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TColor.h"
#include <math.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>

class TraceAverage: public Tool {


 public:

  TraceAverage();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

  static void pipeCloseHandler(int);
  static bool pipeclosed;

 private:

  bool InitTTree(TTree* tree);
  bool CheckFifo();
  
  bool livedraw=false;
  bool hold_max_plot=false;
  bool hold_max_range=false;
  bool plot_gd_region=false;
  double held_max=0;
  bool live_darksub=false;
  std::vector<double> darkvals;
  bool normalise_livedraw=false;
  TCanvas* cspec=nullptr;
  TGraphErrors* ge=nullptr;
  Color_t linecol = kRed;
  std::ofstream fifo;
  std::string fifoname;

  int verbosity=1;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;

};


#endif
