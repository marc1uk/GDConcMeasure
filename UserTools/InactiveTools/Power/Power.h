#ifndef Power_H
#define Power_H

#include <string>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include "Tool.h"

class Power: public Tool {


	public:

		Power();
		bool Initialise(std::string configfile,DataModel &data);
		bool Execute();
		bool Finalise();

		void TurnOn();
		void TurnOff();

	private:

		int powerpin;



};


#endif
