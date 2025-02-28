#include "LivePlotter.h"
#include "TGraph.h"
#include "TCanvas.h"

LivePlotter::LivePlotter():Tool(){}


bool LivePlotter::Initialise(std::string configfile, DataModel &data){

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
  
  
  m_variables.Get("verbosity", verbosity);
  
  return true;
}


bool LivePlotter::Execute(){
  
  Log("LivePlotter Executing...",v_debug,verbosity);
  
  // Check if NewMeasurement flag has been set in the CStore by the MatthewAnalysis Tool.
  bool new_measurement;
  bool get_ok = m_data->CStore.Get("NewMeasurement",new_measurement);
  
  // if so, make a TGraph of measurements so far to publish on the website
  if(get_ok && new_measurement){
    
    // can only make a meaningful graph when we have more than one datapoint, really
    int num_of_meas = m_data->concs_and_errs.size();
    if (num_of_meas > 1){
      
      // concentration vs difference in amplitudes of main absorbance peaks
      // the points on this graph show the sampled range of our calibration plot
      TGraphErrors live_plot(num_of_meas);
      live_plot.SetTitle("Concentration vs Peak Diff;Gd Concentration [%];"
                         "(Absorbance at 273nm) - (Absorbance at 276nm) [no units]");
      
      // fill data
      for (int i = 0; i < num_of_meas; ++i){
        live_plot.SetPoint(i, m_data->concs_and_errs.at(i).first, m_data->peakdiffs_and_errs.at(i).first);
        live_plot.SetPointError(i, m_data->concs_and_errs.at(i).second, m_data->peakdiffs_and_errs.at(i).second);
      }
      
      // draw and save
      TCanvas c1("c1", "c1", 1000, 1000);
      live_plot.Draw();
      c1.SaveAs("/var/www/html/images/live_conc.png");
      
      // TODO add other live plots, e.g. concentration vs time, or vs measurement.
      // or peak height 1&2 vs measurement etc
      
    }
  }
  
  return true;
}


bool LivePlotter::Finalise(){
  
  return true;
}
