#ifndef FakeSpectrometer_H
#define FakeSpectrometer_H

#include <string>
#include <iostream>
#include <random>

#include "Tool.h"

class FakeSpectrometer: public Tool {


 public:

  FakeSpectrometer();
  bool Initialise(std::string configfile,DataModel &data);
  bool Execute();
  bool Finalise();

  bool GetX();
  bool GetDark();
  bool GetPure();
  bool GetData(double cc);
  double Pure(int i);
  double PureFunc(int i);
  double Data(int i, double cc);
  double DataFunc(int i, double cc);
  double Absorb(int i, double cc);
  double Dark(int i = 0);
  double Error();
  double X(int i);

  void Wait();

 private:

  std::mt19937 mt;
  std::uniform_real_distribution<> ran, err;
  std::string m_configfile;

  double noise, pedestal, error;
  int nnn;

  double conc;

};


#endif
