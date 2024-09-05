#include "SaveTraces.h"

SaveTraces::SaveTraces():Tool(){}


bool SaveTraces::Initialise(std::string configfile, DataModel &data){
  
  m_data= &data;
  
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
  
  return true;
}


bool SaveTraces::Execute(){
  
  Log("SaveTraces Executing...",v_debug,verbosity);
  
  std::string save="";
  int overwrite=1;
  
  std::string name="";
  if(m_data->CStore.Get("Filename",name) && name!=lastname){
    Log(m_unique_name+" making new file '"+name+"'",v_debug,verbosity);
    std::string overwritestring = (overwrite==1) ? "RECREATE" : "UPDATE";
    // TODO if update, we need to open the file, retrieve the trees,
    // match them with the ones in the m_data->m_trees, then transfer our entries...
    if(file){
    	Log(m_unique_name+" deleting file that wasn't saved? '"+file->GetName()+"'",v_warning,verbosity);
    	file->Close();
    	delete file;
    }
    file = new TFile(name.c_str(), "RECREATE");
    file->cd();
    lastname=name;
  }
    
  if(m_data->CStore.Get("Save",save) && save=="Save"){
    
    //m_data->CStore.Print();
    save="";
    m_data->CStore.Set("Save",save);
    //m_data->CStore.Get("Filename",name);
    m_data->CStore.Get("Overwrite",overwrite);
    
    Log("SaveTraces saving to filename: '"+name+"'",v_message,verbosity);
    
    
    for (std::map<std::string,TTree*>::iterator it=m_data->m_trees.begin(); it!=m_data->m_trees.end(); ++it){
      it->second->Write();
    }
    
    file->Save();
    file->Close();
    delete file;
    file=nullptr;
    
    /*  automatically deleted when file is closed now
    for (std::map<std::string,TTree*>::iterator it=m_data->m_trees.begin(); it!=m_data->m_trees.end(); ++it){
      delete it->second;
      it->second=0;
    }
    */
    
    m_data->m_trees.clear();
    
    
  }
  
  
  return true;
}


bool SaveTraces::Finalise(){
  
  return true;
}
