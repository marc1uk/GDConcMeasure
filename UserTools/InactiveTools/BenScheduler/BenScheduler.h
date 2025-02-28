#ifndef BenScheduler_H
#define BenScheduler_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/chrono.hpp"

class BenScheduler: public Tool {


 public:

  BenScheduler();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:

  boost::posix_time::time_duration m_period;
  boost::posix_time::ptime current, m_reference_time;


};


#endif
