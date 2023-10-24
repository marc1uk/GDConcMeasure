#include "../Unity.cpp"

Tool* Factory(std::string tool){
Tool* ret=0;

// if (tool=="Type") tool=new Type;
if (tool=="DummyTool") ret=new DummyTool;

//if (tool=="Scheduler") ret=new Scheduler;
//if (tool=="LEDManager") ret=new LEDManager;
//if (tool=="Power") ret=new Power;
//if (tool=="Spectrometer") ret=new Spectrometer;
//if (tool=="Pump") ret=new Pump;
//if (tool=="Calibration") ret=new Calibration;
//if (tool=="Measurement") ret=new Measurement;
//if (tool=="Analysis") ret=new Analysis;
//if (tool=="Output") ret=new Output;
//if (tool=="Writer") ret=new Writer;
//if (tool=="CalibrationManager") ret=new CalibrationManager;
//if (tool=="MeasurementManager") ret=new MeasurementManager;
//if (tool=="StateTest") ret=new StateTest;
//if (tool=="FakeSpectrometer") ret=new FakeSpectrometer;
if (tool=="BenScheduler") ret=new BenScheduler;
if (tool=="Valve") ret=new Valve;
if (tool=="BenPower") ret=new BenPower;
if (tool=="BenSpectrometer") ret=new BenSpectrometer;
if (tool=="TraceAverage") ret=new TraceAverage;
if (tool=="Increment") ret=new Increment;
if (tool=="BenLED") ret=new BenLED;
if (tool=="BenWriter") ret=new BenWriter;
if (tool=="BenAnalysis") ret=new BenAnalysis;
if (tool=="webcoms") ret=new webcoms;
if (tool=="WebScheduler") ret=new WebScheduler;
if (tool=="SaveTraces") ret=new SaveTraces;
if (tool=="InteractiveScheduler") ret=new InteractiveScheduler;
if (tool=="MarcusScheduler") ret=new MarcusScheduler;
if (tool=="GracefulStop") ret=new GracefulStop;
if (tool=="MatthewAnalysis") ret=new MatthewAnalysis;
if (tool=="LivePlotter") ret=new LivePlotter;
if (tool=="SaveToDB") ret=new SaveToDB;
if (tool=="PGTool") ret=new PGTool;
if (tool=="LoadOldFiles") ret=new LoadOldFiles;
if (tool=="MarcusAnalysis") ret=new MarcusAnalysis;
if (tool=="MatthewTransparency") ret=new MatthewTransparency;
 if(tool=="MatthewAnalysisStrikesBack") ret=new MatthewAnalysisStrikesBack;
return ret;
}

