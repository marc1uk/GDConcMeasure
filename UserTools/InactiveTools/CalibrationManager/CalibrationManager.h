#ifndef CalibrationManager_H
#define CalibrationManager_H

#include <string>
#include <iostream>

#include <vector>
#include <map>
#include <algorithm>

#include "TList.h"
#include "TF1.h"
#include "TF2.h"

#include "Tool.h"

class CalibrationManager: public Tool
{
	public:

		CalibrationManager();
		bool Initialise(std::string configfile,DataModel &data);
		bool Execute();
		bool Finalise();

		std::vector<std::string> CalibrationList();
		std::vector<std::string> LoadCalibration(std::vector<std::string> &cList);
		void Load(TFile &f, std::string name, std::string type);
		void Create(std::string name);
		bool IsUpdate(std::string name, int timeUpdate);
		bool Calibrate();

	private:

		std::string m_configfile;
		int verbose;

		bool skipCalib;
		std::string calibFile;		//file in which calibration list is
		std::string concFuncName;	//name of concentration function
		std::string statFuncName;	//name of uncertainity function
		std::string systFuncName;	//name of uncertainity function
		std::string concTreeName;

		std::vector<std::string> calibList, allList;
		std::vector<std::string>::iterator ic;

		std::map<std::string, int> timeUpdate;
};


#endif
