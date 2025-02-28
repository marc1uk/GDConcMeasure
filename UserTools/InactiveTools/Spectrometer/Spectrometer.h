#ifndef Spectrometer_H
#define Spectrometer_H

#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "api/seabreezeapi/SeaBreezeAPI.h"
#include <unistd.h>

#include "Tool.h"

class Spectrometer: public Tool {


	public:

		Spectrometer();
		bool Initialise(std::string configfile,DataModel &data);
		bool Execute();
		bool Finalise();

		bool EstablishUSB();
		bool EstablishUSB_try_catch();
		bool GetData();
		void RelinquishUSB();
		
	private:

		std::string m_configfile;
		int verbose;
		long*  device_ids;
		long* spectrometer_ids;
		int error;

		int intTime;
		int nTraces;

};


#endif
