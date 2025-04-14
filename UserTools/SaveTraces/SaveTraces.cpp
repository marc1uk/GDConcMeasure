#include "SaveTraces.h"

SaveTraces::SaveTraces():Tool(){}


bool SaveTraces::Initialise(std::string configfile, DataModel &data){
  
  InitialiseTool(data);
  m_configfile=configfile;
  InitialiseConfiguration(configfile);
  
  /* - new method, Retrieve configuration options from the postgres database - */
  int RunConfig=-1;
  m_data->vars.Get("RunConfig",RunConfig);
  
  if(RunConfig>=0){
    std::string configtext;
    bool get_ok = m_data->postgres_helper.GetToolConfig(m_unique_name, configtext);
    if(!get_ok){
      Log(m_unique_name+" Failed to get Tool config from database!",v_error,verbosity);
      return false;
    }
    // parse the configuration to populate the m_variables Store.
    if(configtext!="") m_variables.Initialise(std::stringstream(configtext));
    
  }
  
  /* - old method, read config from local file - */
  if(configfile!="")  m_variables.Initialise(configfile);
  
  //m_variables.Print();
  
  m_variables.Get("verbosity",verbosity);
  m_variables.Get("overwrite",overwrite);
  
  return true;
}


bool SaveTraces::Execute(){
  
  Log("SaveTraces Executing...",v_debug,verbosity);
  
  std::string save="";
  
  std::string name="";
  if(m_data->CStore.Get("Filename",name) && name!=lastname){
    Log(m_unique_name+" making new file '"+name+"'",v_debug,verbosity);
    
    std::string filemodestring = "CREATE";
    /*
    if(update) filemodestring= "UPDATE";
    // TODO if wanted to support 'update' mode we would need to open the file, retrieve the trees,
    // match them with the ones in the m_data->m_trees, then transfer our entries... or something?
    // basically update mode not supported for now, probably never will be.
    */
    if(overwrite) filemodestring="RECREATE";
    
    if(file){
      Log(m_unique_name+" deleting file that wasn't saved? '"+file->GetName()+"'",v_warning,verbosity);
      file->Close();
      delete file;
    }
    file = new TFile(name.c_str(), filemodestring.c_str());
    if(file==nullptr || file->IsZombie()){
      Log(m_unique_name+" Error making file '"+name+"'!",v_error,verbosity);
      if(file) file->Close();
      file=nullptr;
      return false;
    }
    file->cd();
    lastname=name;
  }
    
  if(m_data->CStore.Get("Save",save) && save=="Save"){
    
    save="";
    m_data->CStore.Set("Save",save);
    
    Log("SaveTraces saving to filename: '"+name+"'",v_message,verbosity);
    
    if(file){
      
      for (std::map<std::string,TTree*>::iterator it=m_data->m_trees.begin(); it!=m_data->m_trees.end(); ++it){
        it->second->Write();
      }
    
      file->Save();
      file->Close();
      delete file;
      file=nullptr;
      
      // this also deletes the trees as they're owned by the file
      m_data->m_trees.clear();
      
    } else {
      
      Log(m_unique_name+" Save called with no open file!!!",v_error,verbosity);
      for (std::map<std::string,TTree*>::iterator it=m_data->m_trees.begin(); it!=m_data->m_trees.end(); ++it){
        delete it->second;
        it->second=0;
      }
      m_data->m_trees.clear();
      return false;
    }
    
  }
  
  
  return true;
}


bool SaveTraces::Finalise(){
  
  return true;
}
