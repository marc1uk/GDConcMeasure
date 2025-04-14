#ifndef DATAMODEL_H
#define DATAMODEL_H

#include <map>
#include <string>
#include <vector>
#include <utility>

#include "TTree.h"

#include "SlowControlCollection.h"
#include "Store.h"
#include "BoostStore.h"
#include "Logging.h"
#include "DAQDataModelBase.h"

#include "TGraphErrors.h"
#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "GdTree.h"
#include "Postgres.h"
#include "PGHelper.h"

#include "ArduinoController.h"

#include <zmq.hpp>

using namespace ToolFramework;

enum State {
	    ReplaceWater,
	    Sleep,
	    PowerUP,
	    Dark,
	    //Dark1,
	    //Dark2,
	    //LEDW,//
	    //LED385L,//
	    DLED260J,
	    LED260J,
	    DLED275J,  
	    LED275J,
	    //LEDR,//
	    //LEDG,//
	    //LEDB,//
	    //Dark3,
	    DAll,
	    All,
	    Analyse,
	    Count,
	    LEDR,
	    LEDG,
	    LEDB
};

inline State& operator++(State& eDOW, int)  // <--- note -- must be a reference
{
	const int i = static_cast<int>(eDOW)+1;
	eDOW = static_cast<State>((i) % ((int)State::Count));
	return eDOW;
}

//enumeration types (both scoped and unscoped) can have overloaded operators
std::ostream& operator<<(std::ostream& os, State s);
/*enum state
{
	idle,
	power_up,
	power_down,
	init,
	calibration,
	calibration_done,
	measurement,
	measurement_done,
	take_dark,
	take_spectrum,
	turn_off_led,
	analyse,
	change_water,
	settle_water,
	manual_on,
	manual_off,
	finalise
	};*/

class DataModel : public DAQDataModelBase {
  
public:
  
  DataModel();
  
  TTree* GetTTree(std::string name);
  void AddTTree(std::string name, TTree *tree);
  void DeleteTTree(std::string name);
  
  //	GdTree* GetGdTree(std::string name);
  //		void AddGdTree(std::string name, GdTree *tree);
  //		void DeleteGdTree(std::string name);
  //		int SizeGdTree();
  
  //Store vars;
  //BoostStore CStore;
  //std::map<std::string,BoostStore*> Stores;
  
  //Logging *Log;
  
  //zmq::context_t* context;
  
  State state;
  std::string mode;
  std::vector<std::vector<double> > traceCollect;
  std::vector<double> wavelength;	//can be retrieved from seabreeze
  boost::posix_time::ptime measurment_time;

  TGraphErrors* dark_sub_pure = nullptr;
  TF1* pure_fct = nullptr;
  TF1* calib_curve = nullptr;

  std::vector<std::pair<double,double>> concs_and_errs;
  std::vector<std::pair<double,double>> peakdiffs_and_errs;

  //  std::array<double, 2068> current_transparency_values;
  //  std::array<int, 6> current_transparency_date_time;

  //  std::vector<std::array<double, 2068>> transparency_traces;
  //std::vector<std::array<int, 6>> transparency_times;
  
  /*
    state mode;	//state for scheduler
    
    bool isCalibrationTool;
    bool isMeasurementTool;
    bool depleteWater;
    bool circulateWater;
    
    //Calibration manager
    //created and destroyed in the class
    //
    //	used globally
    std::string calibrationName;
    bool isCalibrated;	//true if there is calibratio file!
    bool calibrationDone;	//true if calibration is finished
    //bool calibrationComplete;	//true if calibration can be completed
    
    std::string calibrationTime;	//timestamp for calibration
    //
    //	used in Analysis
    std::string concentrationName;	//name of tree for concentration measurmenet
    //function accepts absorbance in and gives concentration
    TF1 *concentrationFunction;
    TF2 *concentrationFuncStat;	//y variable is error on absorbance
    TF1 *concentrationFuncSyst;
    TFitResultPtr concentrationFit;	//fit result
    std::string fitFuncName;	//name of fit result
    
    double gdconc, gd_err;	//value of gd concentration and error
    //
    //	used in Writing
    std::string calibrationFile;
    std::string calibrationBase;
    
    
    //Calibration manager
    //created and destroyed in the class
    //
    //	used globally
    bool measurementDone;	//true if the measurement is finished
    //
    //	used in Writing & used in analysis
    std::string measurementFile;
    std::string measurementName;
    std::string measurementBase;
    
    //LEDmanager
    //string used by calibration and measurement managers
    bool isLEDon;		//true if led on
    std::string turnOnLED;	//led to be turned on
    
    //Pump
    //
    //override scheduler
    bool turnOnPump;	//must change water
    
    
    //filled by Spectrometer, for Analysis
    //
    std::vector<std::vector<double> > traceCollect;
    std::vector<double> wavelength;	//can be retrieved from seabreeze
    
    std::map<std::string, GdTree*> m_gdtrees; 
  */    
    
  std::map<std::string, TTree*> m_trees;
  
  // database manager
  Postgres postgres;        // manages interface with the database
  PGHelper postgres_helper; // provides user methods to perform routine database operations
  
  ArduinoController* arduino=nullptr;
  
  Store monitoring_store; // put variables in here to send to the monitoring DB by Monitoring tool
  
private:
  
};


#endif
