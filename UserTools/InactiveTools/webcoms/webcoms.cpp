#include "webcoms.h"

webcoms::webcoms():Tool(){}


bool webcoms::Initialise(std::string configfile, DataModel &data){

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
  
  return true;


}


bool webcoms::Execute(){

  zmq::poll(&initems[0], 1, 100);

  if ((initems[0].revents & ZMQ_POLLIN)){

    zmq::message_t msg;
    sock->recv(&msg);
    std::istringstream iss(static_cast<char*>(msg.data()));
    std::cout<<iss.str()<<std::endl;
  }

  
  return true;
}


bool webcoms::Finalise(){

  delete sock;
  sock=0;
  
  return true;
}
