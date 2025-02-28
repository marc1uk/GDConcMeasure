#include "BenScheduler.h"

BenScheduler::BenScheduler():Tool(){}


bool BenScheduler::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  double t=0.0;
  if(!m_variables.Get("wait",t)) t=1;
  m_period=boost::posix_time::time_duration(0, 0, t, 0);

  m_data->state= Analyse;
  
  return true;

    
}


bool BenScheduler::Execute(){

  std::cout<<"Last State is "<<m_data->state<<std::endl;
  boost::posix_time::ptime current(boost::posix_time::second_clock::local_time());
  boost::posix_time::time_duration lapse(m_period - (current - m_reference_time));
  
  if(m_data->state!=Sleep || lapse.is_negative()){
     
    m_data->state++;
     std::cout<<"New State is "<<m_data->state<<std::endl;
    if( m_data->state >= (int)State::Count) m_data->state=ReplaceWater;
    if( m_data->state==Sleep)  m_reference_time=current;
    else if( m_data->state==ReplaceWater) m_data->measurment_time=current;
    std::cout<<"Current State is "<<m_data->state<<std::endl;    
  }
  
  else usleep(100);
  
  return true;
}


bool BenScheduler::Finalise(){


  
  return true;
}
