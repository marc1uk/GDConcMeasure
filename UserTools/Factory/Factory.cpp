#include "Factory.h"
#include "Unity.h"

Tool* Factory(std::string tool) {
Tool* ret=0;

// if (tool=="Type") tool=new Type;
if (tool=="ArduinoControl") ret=new ArduinoControl;
if (tool=="DummyTool") ret=new DummyTool;
if (tool=="GracefulStop") ret=new GracefulStop;
if (tool=="LoadOldFiles") ret=new LoadOldFiles;
if (tool=="MarcusScheduler") ret=new MarcusScheduler;
if (tool=="MatthewAnalysisStrikesBack") ret=new MatthewAnalysisStrikesBack;
if (tool=="MatthewTransparency") ret=new MatthewTransparency;
if (tool=="Monitoring") ret=new Monitoring;
if (tool=="ReturnOfTheMarcusAnalysisEpisode2") ret=new ReturnOfTheMarcusAnalysisEpisode2;
if (tool=="ReturnOfTheMarcusAnalysis") ret=new ReturnOfTheMarcusAnalysis;
if (tool=="SaveToDB") ret=new SaveToDB;
if (tool=="SaveTraces") ret=new SaveTraces;
if (tool=="TraceAverage") ret=new TraceAverage;
if (tool=="BenSpectrometer") ret=new BenSpectrometer;
return ret;
}
