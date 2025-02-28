#ifndef Writer_H
#define Writer_H

#include <string>
#include <iostream>

#include "Tool.h"
#include "TTree.h"
#include "TFile.h"

class Writer: public Tool {


	public:

		Writer();
		bool Initialise(std::string configfile,DataModel &data);
		bool Execute();
		bool Finalise();

		void WriteCalibration();
		void WriteMeasurement();

		void Loop();

	private:





};


#endif
