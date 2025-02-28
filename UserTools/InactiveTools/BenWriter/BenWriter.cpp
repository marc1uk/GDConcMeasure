#include "BenWriter.h"

BenWriter::BenWriter():Tool(){}


bool BenWriter::Initialise(std::string configfile, DataModel &data){
  
  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();
  
  m_data= &data;
  
  //  file = new TFile("date","RECREATE");
  
  //  boost::posix_time::ptime current(boost::posix_time::second_clock::local_time());
  //tm mytm=to_tm(current);
  //m_month=mytm.tm_mon;
  
  
  return true;
}


bool BenWriter::Execute(){
  
  if( m_data->state==PowerUP){
    tm mytm=to_tm(m_data->measurment_time);
    
    std::stringstream filename;
    filename<<"DataY"<<(mytm.tm_year + 1900)<<"M"<<(mytm.tm_mon +1)<<".root";
    file = new TFile(filename.str().c_str(),"UPDATE");
    
    for(State tmp=Dark; tmp<=All; tmp++){    
      std::stringstream name;
      name<<tmp;
      TTree* tree=0;
      tree =(TTree*)file->Get(name.str().c_str());
      if (!tree){
	std::cout<<"BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"<<std::endl;
	tree= new TTree(name.str().c_str(),name.str().c_str());   
	InitTTree(tree);
	//	tree->Write();
      }
      
      TTree* tmptree=m_data->GetTTree(name.str());
      if(tmptree){
	std::cout<<"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"<<std::endl;
	m_data->DeleteTTree(name.str());
      }
      m_data->AddTTree(name.str(),tree);
      
      
    }
  }
  
  else if( m_data->state==Analyse){

    std::cout<<"d1"<<std::endl;
    for( std::map<std::string,TTree*>::iterator it=m_data->m_trees.begin(); it!=m_data->m_trees.end(); it++){
 std::cout<<"d2"<<std::endl;
       it->second->Write("",TObject::kOverwrite);  
       std::cout<<"d3"<<std::endl;
    }

     std::cout<<"d4"<<std::endl;
    file->Save();
     std::cout<<"d5"<<std::endl;
    file->Close();
 std::cout<<"d6"<<std::endl;
    file=0;
 std::cout<<"d7"<<std::endl;
 
    for( std::map<std::string,TTree*>::iterator it=m_data->m_trees.begin(); it!=m_data->m_trees.end(); it++){
  std::cout<<"d8"<<std::endl;      
  //      delete it->second;
       std::cout<<"d9"<<std::endl;
      it->second=0;
       std::cout<<"d10"<<std::endl;
    }

     std::cout<<"d11"<<std::endl;
    m_data->m_trees.clear();
     std::cout<<"d12"<<std::endl;
  }
  
 std::cout<<"d13"<<std::endl;
  
  
  return true;
}


bool BenWriter::Finalise(){
  
  if (file!=0) file->Close();
  
  
    for( std::map<std::string,TTree*>::iterator it=m_data->m_trees.begin(); it!=m_data->m_trees.end(); it++){

      delete it->second;
      it->second=0;
    }

    
    m_data->m_trees.clear();
    
  
  return true;
}

 void BenWriter::InitTTree(TTree* tree){

    std::vector<double> value;
    std::vector<double> error;   
    std::vector<double> wavelength;
   
    short year;
    short month;
    short day;
    short hour; 
    short min;
    
    tree->Branch("value", &value);
    tree->Branch("error", &error);
    tree->Branch("wavelength", &wavelength);
    tree->Branch("year", &year, "year/S");
    tree->Branch("month", &month, "month/S");
    tree->Branch("day", &day, "day/S");
    tree->Branch("hour", &hour, "hour/S");
    tree->Branch("min", &min, "min/S");

    
 }
