#include "InteractiveScheduler.h"

InteractiveScheduler::InteractiveScheduler():Tool(){}


bool InteractiveScheduler::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;
  
  if(!m_variables.Get("command_file",filename)){

    std::cout<<"no conds file set"<<std::endl;
    return false;
  }

  std::string line;
  std::ifstream myfile (filename);
  if (myfile.is_open())
    {
      while ( getline (myfile,line) )
	{
	  commands.push_back(line);
	}
      myfile.close();
    }

  else{
    std::cout << "Unable to open file"<<std::endl; 
    return false;
  }
  
  pos=-1;
  current_measurement="";
  
  return true;

    
}


bool InteractiveScheduler::Execute(){

  if(pos==-1){
    std::string command;
    std::cout<<std::endl<<"Type Start, Stop, or measurment filename > ";
    std::cin>>command;
    if(command=="Start"){
      std::string tmp("{\"Auto\":\"Stop\",\"Power\":\"ON\"}");
      m_data->CStore.JsonParser(tmp);
    }
      
    else if(command=="Stop"){
      std::string tmp("{\"Auto\":\"Stop\",\"Power\":\"OFF\"}");
      m_data->CStore.JsonParser(tmp);
    }
    else if(command!=""){
      pos=0;
      m_data->CStore.JsonParser(commands.at(pos));
      current_measurement=command;
    }
  }
  
  else if(pos==commands.size()){

    std::stringstream tmp;
    tmp<<"{\"Auto\":\"Stop\",\"Filename\":\""<<current_measurement<<"\"}";
    m_data->CStore.JsonParser(tmp.str());
    pos=-1;
  }
  
  else {
    
    pos++;
    if(commands.at(pos)[0]=='{'){
      m_data->measurment_time= boost::posix_time::ptime(boost::posix_time::second_clock::local_time());
      m_data->CStore.JsonParser(commands.at(pos));
    }
    else{
      std::cout<<"in wait"<<std::endl;
      std::stringstream tmp(commands.at(pos));
      long wait;
      tmp>>wait;
      m_period=boost::posix_time::seconds(wait);
      boost::posix_time::ptime current(boost::posix_time::second_clock::local_time());
      boost::posix_time::ptime last(boost::posix_time::second_clock::local_time());
      boost::posix_time::time_duration lapse(m_period - (current - last));
      
      while(!lapse.is_negative()){
	std::cout<<"while test="<<!lapse.is_negative()<<std::endl;
	std::cout<<"lapse="<<lapse<<std::endl;
	sleep(1);
	current=boost::posix_time::second_clock::local_time();
	lapse=boost::posix_time::time_duration(m_period - (current - last));
      }
      
      
    }
    
  }
  
  
  return true;
  
}



bool InteractiveScheduler::Finalise(){
  
  commands.clear();
  
  return true;
}
