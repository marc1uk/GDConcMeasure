#ifndef WebScheduler_H
#define WebScheduler_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/chrono.hpp"
#include <pthread.h>

//Thread_args(zmq::context_t* contextin, ){

//  context=contextin;
  
//}

class WebScheduler: public Tool {


 public:

  WebScheduler();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:
  zmq::socket_t* sock;
  zmq::pollitem_t initems[1];
  
  boost::posix_time::time_duration m_period;
  boost::posix_time::ptime current, m_reference_time;
  BoostStore *tt;

  static void* Thread(void *arg);
  //  Thread_args* args;
  
};


#endif
