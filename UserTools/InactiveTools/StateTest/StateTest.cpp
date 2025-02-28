#include "StateTest.h"

StateTest::StateTest():Tool(){}


bool StateTest::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  
  m_data->mode=state::idle;
  
  
  return true;
}


bool StateTest::Execute(){
  
  switch (m_data->mode){

  case state::idle:
    m_data->mode=state::power_up;
    break;
  case state::power_up:
    sleep(35);
    m_data->mode=state::init;
    break;
  case state::init:	//wake up, connect spectrometer on USB
    m_data->mode=state::calibration;
    break;
  case state::calibration:
    m_data->mode=state::calibration_done;
    break;
  case state::calibration_done:
    m_data->mode=state::measurement;
    break;
  case state::measurement:
    m_data->mode=state::measurement_done;
    break;
  case state::measurement_done:
    m_data->mode=state::finalise;
    break;
  case state::finalise:
    m_data->mode=state::power_down;
    break;
  case state::power_down:
    m_data->mode=state::idle;
    break;
  default:
    break;
  }

  std::cout<<"State is : "<< m_data->mode<<std::endl;
  
  
  
  return true;
}


bool StateTest::Finalise(){

  return true;
}
