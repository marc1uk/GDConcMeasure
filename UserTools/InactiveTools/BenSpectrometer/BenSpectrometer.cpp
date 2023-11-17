#include "BenSpectrometer.h"

BenSpectrometer::BenSpectrometer() : Tool(){}

bool BenSpectrometer::Initialise(std::string configfile, DataModel &data){
  
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
  
  
  // get Tool configuration
  m_variables.Get("traces", nTraces);
  m_variables.Get("integration", intTime);
  m_variables.Get("verbosity", verbosity);
  m_variables.Get("SpectrometerModel",SpectrometerModel);
  m_variables.Get("SpectrometerIP",SpectrometerIP);
  m_variables.Get("SpectrometerPort",SpectrometerPort);
  m_variables.Get("wakedelay",wakedelay);
  
  // default state not connected
  connected=false;
  return true;
}


bool BenSpectrometer::Execute(){
  
  Log("BenSpectrometer Executing...",v_debug,verbosity);
  
  // check if we're being asked to connect
  std::string connect="";
  if(m_data->CStore.Get("Connect",connect)){
    if(connected){
      Log("spectrometer already connected",v_message,verbosity);
    } else {
      
      Log("spectrometer on",v_message,verbosity);
      
      Log("BenSpectrometer sleeping for "+std::to_string(wakedelay)+"s after power-up "
          +"to allow spectrometer to wake up",v_warning,verbosity);
      sleep(wakedelay);
      Log("attempting to connect to spectrometer...",v_message,verbosity);
      // make 50(!) attempts to connect
      // to be honest this seems pointless; i've never seen it connect on any attempt other than the first.
      // if the first fails, all others similarly fail.
      /*
      int tries=0;
      while(!EstablishUSB()){
        tries++;
        Log("Attempt "+std::to_string(tries)+"/50",v_warning,verbosity);
        if(tries>=50){
          Log("couldnt establish connetion to spectrometer",v_error,verbosity);
          return false;
        }
        sleep(1);
      } */
      bool ok = Connect();
      if(not ok){
        Log("BenSpectrometer::Execute failed to connect to spectrometer",v_error,verbosity);
        return false;
      } else {
        Log("BenSpectrometer::Execute successfully connected to spectrometer",v_message,verbosity);
      }
      
      // update our internal connection status
      connected=true;
      
      // update connection status in DataModel for display on website
      m_data->CStore.Set("SpectrometerConnected",true);
    }
  }
  
  // we'll also be disconnected if power went down
  std::string Power="";
  if(m_data->CStore.Get("Power",Power) && Power==false){
      connected=false;
      m_data->CStore.Set("SpectrometerConnected",false);
  }
  
  // check if we have a flag in the DataModel to start a measurement
  std::string tmp="";
  if(m_data->CStore.Get("Measure",tmp) && tmp=="Start"){
    
    Log("got measure",v_message,verbosity);
    // read measurement from spectrometer
    bool ok = GetData();
    if(not ok){
      // should we stop the measurement here, so that downstream Tools
      // do not see the 'Measure' flag and try to process non-existent data?
      m_data->CStore.Remove("Measure");
      return ok;
    }
  
  }
  
  
  return true;
}


bool BenSpectrometer::Finalise(){
  
  //RelinquishUSB();
  return true;
}


bool BenSpectrometer::EstablishUSB_try_catch(){
  // unused
  
  error = 0;
  char nameBuffer[80];
  int flag=1;
  sbapi_initialize();
  sbapi_probe_devices();
  //sbapi_add_TCPIPv4_device_location(SpectrometerModel.c_str(), SpectrometerIP.c_str(), SpectrometerPort);
  
  if (sbapi_get_number_of_device_ids()!=1) return false;
  device_ids = (long *)calloc(1, sizeof(long));
  sbapi_get_device_ids(device_ids, 1);
  flag *= sbapi_get_device_type(device_ids[0], &error, nameBuffer, 79);
  //flag *= sbapi_open_device(device_ids[0], &error);
  int trial = 0;
  while (trial < 4){
    try {
      Log("Connecting spectrometer, attempt "+std::to_string(trial),v_message,verbosity);
      sbapi_open_device(device_ids[0], &error);
      return true;
    }
    catch (const std::exception &e) {
      Log("Can't connect to spectrometer",v_error,verbosity);
      Log("Attempting, trial: "+std::to_string(++trial),v_error,verbosity);
    }
  }
  
  sbapi_open_device(device_ids[0], &error);
  
  if(flag == 0)  return false;
  
  spectrometer_ids = (long *)calloc(1, sizeof(long));
  sbapi_get_spectrometer_features(device_ids[0], &error, spectrometer_ids, 1); /// this may not be needed
  
  return true;
}

bool BenSpectrometer::EstablishUSB(){
  // attempt to open a connection to the spectrometer
  
  error = 0;
  char nameBuffer[80];
  int flag=1;
  sbapi_initialize();
  sbapi_probe_devices();
  sbapi_add_TCPIPv4_device_location(const_cast<char*>(SpectrometerModel.c_str()),
                                    const_cast<char*>(SpectrometerIP.c_str()),
                                    SpectrometerPort);
  if (sbapi_get_number_of_device_ids()!=1) return false;
  device_ids = (long *)calloc(1, sizeof(long));
  sbapi_get_device_ids(device_ids, 1);
  flag *= sbapi_get_device_type(device_ids[0], &error, nameBuffer, 79);
  Log("namebuffer="+std::string(nameBuffer),v_debug,verbosity);
  Log("error="+std::to_string(error),((error==0) ? v_debug : v_error),verbosity);
  //flag *= sbapi_open_device(device_ids[0], &error);
  sleep(30);
  sbapi_open_device(device_ids[0], &error);
  Log("error="+std::to_string(error),((error==0) ? v_debug : v_error),verbosity);
  if(flag ==0 || error!=0) return false;
  spectrometer_ids = (long *)calloc(1, sizeof(long));
  sbapi_get_spectrometer_features(device_ids[0], &error, spectrometer_ids, 1); /// this may not be needed
  if((*spectrometer_ids)==0) return false;
  Log("error="+std::to_string(error),((error==0) ? v_debug : v_error),verbosity);
  Log("spectrometer_ids="+std::to_string(*spectrometer_ids),v_debug,verbosity);
  Log("device_ids[0]="+std::to_string(device_ids[0]),v_debug,verbosity);
  Log("established",v_debug,verbosity);
  return true;
}

bool BenSpectrometer::Connect(){
  
  error = 0;
  
  Log("BenSpectrometer Connect: initializing - ensuring availability of resources (no success check)",v_debug,verbosity);
  /* 
  This should be called prior to any other sbapi_call.
  The API may recover gracefully if this is not called, but future releases may assume this is called first.
  This should be called synchronously – a single thread should call this.
  */
  sbapi_initialize();  // returns void
  
  Log("BenSpectrometer Connect: probing devices - performing scan for recognised "
      "devices on all buses that support autodetection",v_debug,verbosity);
  /*
  Seems like TCP does not support autodetection, and since we use add_TCPIPv4_device,
  this scan is probably redundant. Incidentally sbapi_probe_devices still does not show
  the device even after calling sbapi_add_TCPIPv4_device_location.
  */
  int found_devices = sbapi_probe_devices();
  Log("BenSpectrometer Connect: sbapi_probe_devices found "+std::to_string(found_devices)+" devices",v_debug,verbosity);
  
  Log("BenSpectrometer Connect: going to connect to model "+SpectrometerModel+" on IP "+SpectrometerIP+" port "
      +std::to_string(SpectrometerPort),v_debug,verbosity);
  
  Log("BenSpectrometer Connect: registering spectrometer IP address as a device location",v_debug,verbosity);
  int ok = sbapi_add_TCPIPv4_device_location(const_cast<char*>(SpectrometerModel.c_str()),
                                         const_cast<char*>(SpectrometerIP.c_str()),
                                         SpectrometerPort);
  Log(std::string("BenSpectrometer Connect: sbapi_add_TCPIPv4_device returned ")+((ok==0) ? "success" : "failure"),v_debug,verbosity);
  
  Log("BenSpectrometer Connect: getting number of available devices",v_debug,verbosity);
  /*
     This returns the number of devices either found by probe_devices, or explicitly specified
     by sbapi_add_***_device_*. There appears to be no actual check that such specified devices
     are indeed accessible, so this is likely to always return 1.
  */
  int num_available_devices = sbapi_get_number_of_device_ids();
  Log("BenSpectrometer Connect: number of available devices = "+std::to_string(num_available_devices),v_debug,verbosity);
  
  if (num_available_devices<1){
    Log("BenSpectrometer Connect: giving up as we found no devices",v_error,verbosity);
    return false;
  }
  
  long* my_device_ids = new long[num_available_devices];
  Log("BenSpectrometer Connect: getting device ids for available devices",v_debug,verbosity);
  /*
     populates up to N entries in my_device_ids, where N is the second argument, and returns
     the number of device ids actually obtained. Note that it may return non-zero device ids
     for manually added devices, regardless of whether they are actually present.
  */
  int num_device_ids = sbapi_get_device_ids(my_device_ids, num_available_devices);
  Log("BenSpectrometer Connect: got "+std::to_string(num_device_ids)+" device ids:",v_debug,verbosity);
  for(int i=0; i<num_device_ids; ++i){
     Log("BenSpectrometer Connect: device "+std::to_string(i)+" has id "+std::to_string(my_device_ids[i]),v_debug,verbosity);
  }
  
  int max_chars = 79;
  char* nameBuffer = new char[max_chars];
  Log("BenSpectrometer Connect: getting names of devices",v_debug,verbosity);
  /*
     again this seems to return no error, and even returns a valid name (in uppercase,
     i.e. different to that given to sbapi_add_TCPIPv4_device_location), even when the
     device is not present!
  */
  for(int i=0; i<num_device_ids; ++i){
    memset(nameBuffer, 0, max_chars);
    nameBuffer[0]= '\0';
    int nbytes = sbapi_get_device_type(my_device_ids[i], &error, nameBuffer, max_chars);
    if(error==0){
       Log("BenSpectrometer Connect: no error ",v_debug,verbosity);
    } else {
      Log("BenSpectrometer Connect: got error code "+std::to_string(error)
             +sbapi_get_error_string(error)+" when getting type for device "+std::to_string(i),v_error,verbosity);
    }
    Log("BenSpectrometer Connect: device type was "+std::to_string(nbytes)+" chars long: '"+nameBuffer+"'",v_debug,verbosity);
  }
  
  Log("BenSpectrometer Connect: sleep for 30 seconds...",v_debug,verbosity);
  sleep(30);
  Log("BenSpectrometer Connect: resuming",v_debug,verbosity);
  
  long device_id = my_device_ids[0];
  Log("BenSpectrometer Connect: opening first device, id "+std::to_string(device_id),v_debug,verbosity);
  /*
     it is only at this point that a connection to the device is attempted, and we can
     really tell if the device is actually present and functioning.
  */
  error=0;
  ok = sbapi_open_device(device_id, &error);
  if(ok==0){
    Log("BenSpectrometer Connect: opened ok!",v_debug,verbosity);
  } else {
    Log("BenSpectrometer Connect: failed to open device! error was "+std::to_string(error)
        +" = "+sbapi_get_error_string(error),v_error,verbosity);
    return false;
  }
  
  Log("BenSpectrometer Connect: getting device "+std::to_string(device_id)+" number of features",v_debug,verbosity);
  error=0;
  int max_feature_codes = sbapi_get_number_of_spectrometer_features(device_id, &error);
  if(error!=0){
    Log("BenSpectrometer Connect: error "+std::to_string(error)+" = "+sbapi_get_error_string(error)+
        " getting number of features for device "+std::to_string(device_id),v_error,verbosity);
    return false;
  }
    
  long* features_codes = new long[max_feature_codes];
  error=0;
  /*
     Somewhat oddly named, this gets 'instance ids', where "each instance refers to a single
     optical bench". We need to specify an instance id when we retrieve data.
     For reference we expect 2 "features" for our Blaze, and we use the first.
     If the spectrometer is not connected, this returns no features, but also no error!
  */
  int num_features = sbapi_get_spectrometer_features(device_id, &error, features_codes, max_feature_codes);
  std::string msgg="BenSpectrometer Connect: got "+std::to_string(num_features)+" feature codes";
  if(max_feature_codes==num_features){
    msgg+=" (limited - there may be more)";
  }
  Log(msgg,v_debug,verbosity);
  for(int i=0; i<num_features; ++i){
        Log("BenSpectrometer Connect: feature "+std::to_string(i)+" had code "
            +std::to_string(features_codes[i]),v_debug,verbosity);
  }
  if(error!=0){
        Log("BenSpectrometer Connect: returned error code "+std::to_string(error)+" = "
            +sbapi_get_error_string(error),v_error,verbosity);
  }
  
  long feature_id = features_codes[0];
  Log("BenSpectrometer Connect: feature_id is "+std::to_string(feature_id),v_debug,verbosity);
  
  Log("BenSpectrometer Connect: consider this established",v_debug,verbosity);
  
  device_ids = (long *)calloc(1, sizeof(long));
  device_ids[0]=device_id;
  spectrometer_ids = (long *)calloc(1, sizeof(long));
  spectrometer_ids[0]=feature_id;
  
  Log("BenSpectrometer Connect: cleanup",v_debug,verbosity);
  delete[] my_device_ids;
  delete[] features_codes;
  
  Log("BenSpectrometer Connect: returning",v_debug,verbosity);
  return true;
}

void BenSpectrometer::RelinquishUSB(){
  
  // close session, if we have one open
  if(device_ids[0]>0){
    sbapi_close_device(device_ids[0], &error);
    if(error!=0){
      Log(std::string("BenSpectrometer::RelinquishUSB - sbapi_close_device returned error code ")
          +std::to_string(error)+std::string(", ")+sbapi_get_error_string(error),v_error,verbosity);
    }
  }
  
  // free resources
  free(device_ids);
  free(spectrometer_ids);
  sbapi_shutdown();   // no success checks available
  
}

bool BenSpectrometer::GetData(){
  
  if(device_ids==nullptr || device_ids[0]<0){
    Log("BenSpectrometer::GetData has no spectrometer open!",v_error,verbosity);
    return false;
  }
  /*
   * take N spectrometer traces and stores result in DataModel
   */
  
  error = 0;
  double *data           = 0;
  double *wavelength = 0;
  
  // query amount of memory required to store a spectrum
  int size = sbapi_spectrometer_get_formatted_spectrum_length(device_ids[0],spectrometer_ids[0], &error);
  if(size==0 || error!=0){
    Log("BenSpectrometer::GetData failed getting trace length! sbapi_spectrometer_get_formatted_spectrum_length"
        " returned "+std::to_string(size)+", error code "+std::to_string(error),v_error,verbosity);
    return false;
  }
  Log("size="+std::to_string(size),v_message,verbosity);
  
  // allocate storage
  data           = (double *)calloc(size, sizeof(double));
  wavelength = (double *)calloc(size, sizeof(double));
  
  // configure measurement integration time
  sbapi_spectrometer_set_integration_time_micros(device_ids[0], spectrometer_ids[0],&error, intTime);
  if(error){
    Log("BenSpectrometer::GetData sbapi_spectrometer_set_integration_time_micros returned error "
       +std::to_string(error)+std::string(": ")+sbapi_get_error_string(error),v_error,verbosity);
    return 0;
  }
  
  // clear old data in DataModel
  m_data->traceCollect.clear();
  
  // loop over N traces (we'll average them later to reduce statistical uncertainty)
  for(int i = 0; i < nTraces; i++){
    
    // retrieve this spectrum
    int ret = sbapi_spectrometer_get_formatted_spectrum(device_ids[0], spectrometer_ids[0],&error, data, size);
    if(ret==0 || error!=0){
      Log("BenSpectrometer::GetData failed getting data! sbapi_spectrometer_get_formatted_spectrum returned "
          +std::to_string(ret)+", error code "+std::to_string(error),v_error,verbosity);
      return false;
    }
    
    // convert to vector
    std::vector<double> trace(size);
    memcpy(&trace[0], &data[0], size*sizeof(double));
    
    // add this trace to datamodel vector
    m_data->traceCollect.push_back(trace);
  }
  
  // wavelengths are always the same so we only need one copy in the datamodel
  int ret = sbapi_spectrometer_get_wavelengths(device_ids[0], spectrometer_ids[0],
                                     &error, wavelength, size);
  if(ret==0 || error!=0){
    Log("BenSpectrometer::GetData failed getting wavevelengths! sbapi_spectrometer_get_wavelengths returned "
        +std::to_string(ret)+", error code "+std::to_string(error),v_error,verbosity);
    return false;
  }
  m_data->wavelength.resize(size);
  memcpy(&(m_data->wavelength)[0], &wavelength[0], size*sizeof(double));
  
  // free temporaries for current trace
  free(data);
  free(wavelength);
  
  return true;
}


