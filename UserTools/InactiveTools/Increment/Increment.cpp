#include "Increment.h"

Increment::Increment():Tool(){}


bool Increment::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;

  return true;
}


bool Increment::Execute(){

  int a;
  std::cout<<"Increement?";
  std::cin>>a;
  std::cout<<std::endl;
  
  return true;
}


bool Increment::Finalise(){

  return true;
}
