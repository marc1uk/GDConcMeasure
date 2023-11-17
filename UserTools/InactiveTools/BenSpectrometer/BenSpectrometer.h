#ifndef BenSpectrometer_H
#define BenSpectrometer_H

#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "api/seabreezeapi/SeaBreezeAPI.h"
#include <unistd.h>

#include "Tool.h"

class BenSpectrometer: public Tool {
  
  
public:
  
  BenSpectrometer();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();
  
  bool EstablishUSB();
  bool EstablishUSB_try_catch();
  bool Connect();
  bool GetData();
  void RelinquishUSB();
  
private:
  
  std::string m_configfile;
  std::string SpectrometerModel="Blaze";
  std::string SpectrometerIP="192.168.50.4";
  int SpectrometerPort=57357;
  long*  device_ids = nullptr;
  long* spectrometer_ids = nullptr;
  int error;
  
  int intTime;
  int nTraces;
  bool connected;
  int wakedelay;
  
  int verbosity=1;
  int v_error=0;
  int v_warning=1;
  int v_message=2;
  int v_debug=3;
};


#endif
