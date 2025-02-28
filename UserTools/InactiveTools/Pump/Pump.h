#ifndef Pump_H
#define Pump_H

#include <string>
#include <iostream>

#include "Tool.h"

class Pump: public Tool {


	public:

		Pump();
		bool Initialise(std::string configfile,DataModel &data);
		bool Execute();
		bool Finalise();

		void TurnOn();
		void TurnOff();

	private:

		int pumppin;



};


#endif
