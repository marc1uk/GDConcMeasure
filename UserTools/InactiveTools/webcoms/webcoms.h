#ifndef webcoms_H
#define webcoms_H

#include <string>
#include <iostream>
#include <sstream>
#include "Tool.h"

class webcoms: public Tool {


 public:

  webcoms();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();


 private:

  zmq::socket_t* sock;
  zmq::pollitem_t initems[1];


};


#endif
