#include "BenSpectrometer.h"

BenSpectrometer::BenSpectrometer() : Tool(){}

bool BenSpectrometer::Initialise(std::string configfile, DataModel &data){
  
  if(configfile!=""){
    //m_configfile = configfile;
    m_variables.Initialise(configfile);
  }
	//m_variables.Print();
  
  m_data = &data;
  
  m_variables.Get("traces", nTraces);
  m_variables.Get("integration", intTime);

  power="OFF";
  return true;
}


bool BenSpectrometer::Execute(){

  std::string Power="";

  if(m_data->CStore.Get("Power",Power) && Power!=power){

    if(Power=="ON"){
      power=Power;
      std::cout<<"spectromiter on"<<std::endl;
      int tries=0;
      while(!EstablishUSB()){
	tries++;
	std::cout<<"Attempt "<<tries<<":50"<<std::endl;
	if(tries>=50){
	  std::cout<<"couldnt establish connetion to spectromiter"<<std::endl;
	  return false;
	}
	sleep(1);
      }
      
    }
    else if (Power=="OFF") Power=power;

  }

  std::string tmp="";

  if(m_data->CStore.Get("Measure",tmp) && tmp=="Start"){
    std::cout<<"got measure"<<std::endl;
    GetData();

  }
  
  //old
  /*
  if(m_data->state==PowerUP){

    int tries=0;
    while(!EstablishUSB()){
      tries++;
      std::cout<<"Attempt "<<tries<<":50"<<std::endl;
      if(tries>=50){
	std::cout<<"couldnt establish connetion to spectromiter"<<std::endl;
	return false;
      }
      sleep(1);
    }
    
  }
  else if(m_data->state==Analyse) RelinquishUSB();
  else if(m_data->state>PowerUP && m_data->state<Analyse) GetData();
  */

  //old old
  /*
    switch (m_data->mode)
    {
    case state::init:	//wake up, connect spectrometer on USB
			//try and catch statement
			//Initialise(m_configfile, *m_data);
			if(!EstablishUSB())
			return false;
			//GetData();	//measure dark noise on wake up
			break;
			//if (!m_data->calibrationComplete || m_data->calibrationName != m_data->concentrationName)
			//	GetData();
			//break;
			case state::calibration:	//take dark cause LED is off
			case state::measurement:
			break;
			case state::take_dark:
			case state::take_spectrum:
			GetData();
			break;
			case state::finalise:
			//Finalise();
			RelinquishUSB();
			break;
			default:
			break;
			}
  */
  return true;
}


bool BenSpectrometer::Finalise(){
  
  return true;
}


bool BenSpectrometer::EstablishUSB_try_catch(){
  
  error = 0;
  char nameBuffer[80];
  int flag=1;
  sbapi_initialize();
  sbapi_probe_devices();
  // sbapi_add_TCPIPv4_device_location("Blaze", "192.168.1.151", 57357);
  //sbapi_add_TCPIPv4_device_location("Jaz", "192.168.1.150", 7654);
  
  if (sbapi_get_number_of_device_ids()!=1) return false;
  device_ids = (long *)calloc(1, sizeof(long));
  sbapi_get_device_ids(device_ids, 1);
  flag *= sbapi_get_device_type(device_ids[0], &error, nameBuffer, 79);
  // flag *= sbapi_open_device(device_ids[0], &error);
  int trial = 0;
  while (trial < 4){
    try {
      std::cout << "Connecting spectrometer, attempt " << trial << std::endl;
      sbapi_open_device(device_ids[0], &error);
      return true;
    }
    catch (const std::exception &e) {
      std::cout << "Can't connect to spectrometer\n";
      std::cout << "Attempting, trial: " << ++trial << std::endl;
    }
  }
  
  sbapi_open_device(device_ids[0], &error);
  
  if(flag == 0)  return false;
  
  spectrometer_ids = (long *)calloc(1, sizeof(long));
  sbapi_get_spectrometer_features(device_ids[0], &error, spectrometer_ids, 1); /// this may not be needed
}

bool BenSpectrometer::EstablishUSB(){
  
  error = 0;
  char nameBuffer[80];
  int flag=1;
  sbapi_initialize();
  sbapi_probe_devices();
  sbapi_add_TCPIPv4_device_location("Blaze", "192.168.50.2", 57357);
  //  sbapi_add_TCPIPv4_device_location("Jaz", "192.168.1.150", 7654);
  //sbapi_add_TCPIPv4_device_location("Jaz", "192.168.50.2", 7654);
  if (sbapi_get_number_of_device_ids()!=1) return false;
  device_ids = (long *)calloc(1, sizeof(long));
  sbapi_get_device_ids(device_ids, 1);
  flag *= sbapi_get_device_type(device_ids[0], &error, nameBuffer, 79);
  std::cout<<"namebuffer="<<nameBuffer<<std::endl;
  std::cout<<"error="<<error<<std::endl;
  //flag *= sbapi_open_device(device_ids[0], &error);
  sleep(30);
  sbapi_open_device(device_ids[0], &error);
  std::cout<<"error="<<error<<std::endl;
  if(flag ==0) return false;
  spectrometer_ids = (long *)calloc(1, sizeof(long));
  std::cout<<"spectrometer_ids="<<spectrometer_ids<<std::endl;
  std::cout<<"device_ids[0]="<<device_ids[0]<<std::endl;
  sbapi_get_spectrometer_features(device_ids[0], &error, spectrometer_ids, 1); /// this may not be needed
  std::cout<<"error="<<error<<std::endl;
  std::cout<<"spectrometer_ids="<<spectrometer_ids<<std::endl;
  std::cout<<"device_ids[0]="<<device_ids[0]<<std::endl;
  std::cout<<"established"<<std::endl;
  return true;
}

void BenSpectrometer::RelinquishUSB(){
  
  sbapi_close_device(device_ids[0], &error);
  free(device_ids);
  free(spectrometer_ids);
  sbapi_shutdown();
  
}

bool BenSpectrometer::GetData(){
  
  /*
   * take measurement with LED off
   * this should be background to following measurements
   * save in vTraceCollect all traces
   * save in xAxis the wavelength of x-axis
   */
  error = 0;
  // spectrometer_ids = 0;
  //  std::vector<std::vector<double> > datavec;
  double *data	   = 0;
  double *wavelength = 0;
  
  int size = sbapi_spectrometer_get_formatted_spectrum_length(device_ids[0],spectrometer_ids[0], &error);
  std::cout<<"size="<<size<<std::endl;
  sbapi_spectrometer_set_integration_time_micros(device_ids[0], spectrometer_ids[0],&error, intTime);
  data	   = (double *)calloc(size, sizeof(double));
  wavelength = (double *)calloc(size, sizeof(double));
  
  m_data->traceCollect.clear(); 
  for(int i = 0; i < nTraces; i++){
    
    sbapi_spectrometer_get_formatted_spectrum(device_ids[0], spectrometer_ids[0],&error, data, size);
    std::vector<double> trace(size);
    memcpy(&trace[0], &data[0], size*sizeof(double));
    m_data->traceCollect.push_back(trace);
  }
  
  sbapi_spectrometer_get_wavelengths(device_ids[0], spectrometer_ids[0],
				     &error, wavelength, size);
  m_data->wavelength.resize(size);
  memcpy(&(m_data->wavelength)[0], &wavelength[0], size*sizeof(double));
  
  free(data);
  free(wavelength);
  
  return true;
}


