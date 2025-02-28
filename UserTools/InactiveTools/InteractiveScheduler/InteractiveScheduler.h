#ifndef InteractiveScheduler_H
#define InteractiveScheduler_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "Tool.h"
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/chrono.hpp"
#include <pthread.h>

//Thread_args(zmq::context_t* contextin, ){

//  context=contextin;
  
//}

class InteractiveScheduler: public Tool {


 public:

  InteractiveScheduler();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:
  
  boost::posix_time::time_duration m_period;
  std::vector<std::string> commands;
  int pos;
  std::string filename;
  std::string current_measurement;
  
};


#endif
