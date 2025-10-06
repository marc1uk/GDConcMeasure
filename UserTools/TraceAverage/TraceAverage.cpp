#include "TraceAverage.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TH1.h"

TraceAverage::TraceAverage():Tool(){}


bool TraceAverage::Initialise(std::string configfile, DataModel &data){
  
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
  
  livedraw = false;
  m_variables.Get("livedraw",livedraw);
  if(livedraw){
      Log("TraceAverage: live draw enabled, making TApp",v_debug,verbosity);
      // make a TApplication for live viewing the spectrum
	  int tapp_users=1;
	  if(m_data->CStore.Get("tapp_users",tapp_users)){
	  	++tapp_users;
	  	m_data->CStore.Set("tapp_users",tapp_users);
	  } else {
	    TApplication* tapp = new TApplication("TApp",0,0,0);
	    intptr_t tapp_p = reinterpret_cast<intptr_t>(tapp);
	    m_data->CStore.Set("tapp",tapp_p);
	    m_data->CStore.Set("tapp_users",tapp_users);
	  }
	  // make canvas to draw it on
          Log("TraceAverage: making live plot canvas",v_debug,verbosity);
	  cspec = new TCanvas("cspec","cspec",1200,700);
  }
  
  return true;
}


bool TraceAverage::Execute(){
  
  Log("TraceAverage Executing...",v_debug,verbosity);
  
  std::string measure;
  if(m_data->CStore.Get("Measure",measure) && measure == "Start"){
    
    Log("TraceAverage new data to process...",v_debug,verbosity);
    
    std::vector<double> value;
    std::vector<double> error;
    std::vector<double> wavelength=m_data->wavelength;
    
    // convert timestamp to string
    tm mytm=to_tm(m_data->measurment_time);
    short year=mytm.tm_year + 1900;
    short month=mytm.tm_mon + 1;
    short day=mytm.tm_mday;
    short hour=mytm.tm_hour;
    short min=mytm.tm_min;
    short sec=mytm.tm_sec;
    std::string date=to_simple_string(m_data->measurment_time);
    Log("TraceAverage timestamp for this measurement is " + date,v_debug,verbosity);
    
    // get trace name - we'll use this for the TTree name
    std::string name="";
    m_data->CStore.Get("Trace",name);
    if(name=="") name="test";
    
    // see if we have this TTree in the datamodel, or build one if not
    TTree* tree=0;
    tree=m_data->GetTTree(name);
    if(tree==0){
      tree=new TTree(name.c_str(),name.c_str());
      m_data->AddTTree(name,tree);
      bool ok = InitTTree(tree);
      //if(!ok) return false;  // XXX is it better that we at least make a graph for the website...?
    }
    
    // make local variables to hold measurement data
    std::vector<double>* a=&value;
    std::vector<double>* b=&error;
    std::vector<double>* c=&wavelength;
    bool err=false;
    // setup tree addresses
    err |= (tree->SetBranchAddress("value", &a) < 0);
    err |= (tree->SetBranchAddress("error", &b) < 0);
    err |= (tree->SetBranchAddress("wavelength", &c) < 0);
    err |= (tree->SetBranchAddress("year", &year) < 0);
    err |= (tree->SetBranchAddress("month", &month) < 0);
    err |= (tree->SetBranchAddress("day", &day) < 0);
    err |= (tree->SetBranchAddress("hour", &hour) < 0);
    err |= (tree->SetBranchAddress("min", &min) < 0);
    err |= (tree->SetBranchAddress("sec", &sec) < 0);
    if(err){
      Log("TraceAverage::Execute error setting branch addresses",0,0);
      // don't return, is it better that we at least make a graph for the website...?
    }
    
    // retrieve data from datamodel, averaging over spectrometer acquisitions
    double maxval=0;
    for(int i=0;i<m_data->wavelength.size(); i++){
      double valsum=0;
      double errorsum=0;
      
      for(int trace=0;trace<m_data->traceCollect.size(); trace++){
        
        valsum+=m_data->traceCollect.at(trace).at(i)/((double)m_data->traceCollect.size());
        if(m_data->traceCollect.at(trace).at(i)>maxval) maxval=m_data->traceCollect.at(trace).at(i);
      }
      for(int trace=0;trace<m_data->traceCollect.size(); trace++){
        
        errorsum+=((m_data->traceCollect.at(trace).at(i)-valsum)*(m_data->traceCollect.at(trace).at(i)-valsum))/(((double)m_data->traceCollect.size())-1.0);
      }
      
      errorsum =sqrt(errorsum)/sqrt(m_data->traceCollect.size());
      value.push_back(valsum);
      error.push_back(errorsum);
    }
    if((name=="Dark" || name=="dark") && maxval>2000){
      Log("TraceAverage tool found Dark trace with peak value of "+std::to_string(maxval)+" - could be LED stuck on! "
          "Terminating toolchain to power down as a precaution!",v_error,verbosity);
      m_data->vars.Set("StopLoop",1);
    }
    
    // fill TTree
    int nbytes = tree->Fill();
    if(nbytes<=0){
      Log("TraceAverage::Execute error calling TTree::Fill, returned "+std::to_string(nbytes),0,0);
      // don't return, is it better that we at least make a graph for the website...?
    }
    // we should reset the branch addresses as they are currently pointed at local
    // variables that will soon go out of scope.
    // if later tools call 'm_data->trees[this_tree]->GetEntry()'
    // without setting the branch addresses, it can segfault if we don't.
    tree->ResetBranchAddresses();
    
    // make a TGraph for the website
    TGraphErrors gr(value.size(),&wavelength[0],&value[0],0,&error[0]);
    gr.SetTitle("TGraphErrors Example");
    //gr.SetMarkerColor(4);
    //gr.SetMarkerStyle(21);
    
    // for live viewing
    if(livedraw){
	    if(ge) delete ge;
	    ge = new TGraphErrors(value.size(), wavelength.data(), value.data(), nullptr, error.data());
	    cspec->cd();
	    if(linecol==kRed) linecol=kBlue;
	    else linecol=kRed;
	    ge->SetLineColor(linecol);
	    if(*std::max_element(value.begin(), value.end())>maxvalue){
	    	maxvalue=*std::max_element(value.begin(), value.end());
	    }
	    ge->GetHistogram()->GetYaxis()->SetRangeUser(-50,maxvalue*1.1);
	    ge->Draw("AL");
	    cspec->Modified();
	    cspec->Update();
	    gSystem->ProcessEvents();
    }
    
    //TCanvas c1("c1","A Simple Graph with error bars",200,10,700,500);
    //gr.Draw("AP");
    //c1.SaveAs("/home/pi/WebServer/html/images/last.png");
    //gr.Write(tmp.str().c_str(),TObject::kOverwrite);
    
    m_data->CStore.Remove("Measure");
    m_data->CStore.Set("NewData",true);  // flag for downstream tools
    
  } else {
    m_data->CStore.Set("NewData",false);
  }
  
  
  return true;
}


bool TraceAverage::Finalise(){
  
  if(livedraw){
	  if(cspec) delete cspec;
	  cspec=nullptr;
	  
	  // TODO move this to a DataModel Method ('AddTAppUser / RemoveTAppUser')
	  int tapp_users=1;
	  m_data->CStore.Get("tapp_users",tapp_users);
	  --tapp_users;
	  m_data->CStore.Set("tapp_users",tapp_users);
	  if(tapp_users==0){
	    intptr_t tapp_p;
	    m_data->CStore.Get("tapp",tapp_p);
	    TApplication* tapp = reinterpret_cast<TApplication*>(tapp_p);
	    delete tapp;
	    tapp = nullptr;
	    tapp_p = 0;
	    m_data->CStore.Set("tapp",tapp_p);
	  }
  }
  
  return true;
}

bool TraceAverage::InitTTree(TTree* tree){
  
  std::vector<double> value;
  std::vector<double> error;
  std::vector<double> wavelength;
  
  short year;
  short month;
  short day;
  short hour;
  short min;
  short sec;
  
  bool bptr = true;
  bptr &= bool(tree->Branch("value", &value));
  bptr &= bool(tree->Branch("error", &error));
  bptr &= bool(tree->Branch("wavelength", &wavelength));
  bptr &= bool(tree->Branch("year", &year, "year/S"));
  bptr &= bool(tree->Branch("month", &month, "month/S"));
  bptr &= bool(tree->Branch("day", &day, "day/S"));
  bptr &= bool(tree->Branch("hour", &hour, "hour/S"));
  bptr &= bool(tree->Branch("min", &min, "min/S"));
  bptr &= bool(tree->Branch("sec", &sec, "sec/S"));
  
  if(!bptr){
    Log("TraceAverage::InitTTree error making TTree branches!",0,0);
    return false;
  }
  return true;
  
}
