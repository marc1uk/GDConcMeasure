#include "WebScheduler.h"

WebScheduler::WebScheduler():Tool(){}


bool WebScheduler::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;

  sock=new zmq::socket_t(*m_data->context, ZMQ_SUB);
  //  sock->bind("ipc:///tmp/WebInput");
  sock->bind("tcp://127.0.0.1:5555");
  sock->setsockopt(ZMQ_SUBSCRIBE, "", 0);
  
  initems[0].socket=*sock;
  initems[0].fd=0;
  initems[0].events=ZMQ_POLLIN;
  initems[0].revents=0;


  double t=0.0;
  if(!m_variables.Get("wait",t)) t=1;
  m_period=boost::posix_time::time_duration(0, 0, t, 0);
  
  //m_data->state= Analyse;

  m_data->state= Sleep;
  m_data->mode="Manual";
  
  return true;

    
}


bool WebScheduler::Execute(){

  boost::posix_time::ptime current(boost::posix_time::second_clock::local_time());
  m_data->measurment_time=current;
  
  zmq::poll(&initems[0], 1, 100);
  
  if ((initems[0].revents & ZMQ_POLLIN)){
    
    zmq::message_t msg;
    sock->recv(&msg);
    std::istringstream iss(static_cast<char*>(msg.data()));
    m_data->CStore.JsonParser(iss.str());
    
  }
  //m_data->CStore.Print();
  
  std::string tmp="";
  if(m_data->CStore.Get("Auto",tmp)){
    if(tmp=="Start") m_data->mode="Auto";
    else if(tmp=="Stop") m_data->mode="Manual";
    
  }
  
  if(m_data->mode=="Manual"){
    
  }
  /*
    else if(m_data->mode=="Auto"){    
    
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
    
  }
  */
  return true;
  
}


bool WebScheduler::Finalise(){

  delete sock;
  sock=0;
  
  return true;
}
