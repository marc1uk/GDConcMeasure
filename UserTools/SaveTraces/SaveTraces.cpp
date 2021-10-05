#include "SaveTraces.h"

SaveTraces::SaveTraces():Tool(){}


bool SaveTraces::Initialise(std::string configfile, DataModel &data){

  if(configfile!="")  m_variables.Initialise(configfile);
  //m_variables.Print();

  m_data= &data;

  return true;
}


bool SaveTraces::Execute(){
  std::string save="";
  std::string name="";
  int overwrite=1;
    
  if(m_data->CStore.Get("Save",save) && save=="Save"){

    m_data->CStore.Print();
    save="";
    m_data->CStore.Set("Save",save);
    m_data->CStore.Get("Filename",name);
    m_data->CStore.Get("Overwrite",overwrite);

    std::cout<<"filename="<<name<<std::endl;
    std::string overwritestring = (overwrite==1) ? "RECREATE" : "UPDATE";
    // TODO if update, we need to open the file, retrieve the trees,
    // match them with the ones in the m_data->m_trees, then transfer our entries...
    TFile file(name.c_str(), "RECREATE");
    file.cd();
    
    for (std::map<std::string,TTree*>::iterator it=m_data->m_trees.begin(); it!=m_data->m_trees.end(); ++it){
      it->second->Write();
    }

    file.Save();
    file.Close();

  
  
  for (std::map<std::string,TTree*>::iterator it=m_data->m_trees.begin(); it!=m_data->m_trees.end(); ++it){
    delete it->second;
    it->second=0;
  }

  m_data->m_trees.clear();

    
  }

  
  return true;
}


bool SaveTraces::Finalise(){

  return true;
}
