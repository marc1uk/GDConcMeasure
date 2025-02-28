#include "Output.h"

Output::Output():Tool(){}


bool Output::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;

  return true;
}


bool Output::Execute(){

  return true;
}


bool Output::Finalise(){

  return true;
}
