#ifndef LEDManager_H
#define LEDManager_H

#include <string>
#include <iostream>
#include <errno.h>

#include "Tool.h"

#include "wiringPiI2C.h"

enum LED
{
	led_260 = 1,
	led_275 = 2,
	led_385 = 4,
	led_405 = 8,
	led_W7  = 16,
	led_R   = 32,
	led_G   = 64,
	led_B   = 128,
};

class LEDManager: public Tool
{
	public:
	
		LEDManager();
		bool Initialise(std::string configfile,DataModel &data);
		bool Execute();
		bool Finalise();

		void Configure();
		bool EstablishI2C();
		void MapLED();
		bool SetPWMfreq(double freq);
		
		bool TurnOn();
		bool TurnOff();
		bool TurnOffAndSleep();
		bool IsSleeping();
		bool WakeUpDriver();
		bool SleepDriver();

		bool TurnLEDArray(unsigned int ledON);
		bool TurnLEDon(std::string ledName);
		bool TurnLEDoff(std::string ledName);

		bool Read(int reg, int &data);
		bool Write(int reg, int data);
		bool ReadAI(int reg, int num_reg, std::vector<int> &block);
		bool WriteAI(int reg, const std::vector<int> &block);
	
	
	private:
		std::string m_configfile;
		//verbosity
		int verbose;
		//set up by configure
		std::string wiringLED, configLED;
		unsigned int resolution;
		double fVin, fDelay, frequencyPWM;
		unsigned long lastTime;

		//counting and flags
		int file_descript;
		unsigned int ledONstate;

		//map between led name and its "binary" value
		std::map<std::string, unsigned int> mLED_name;
		std::map<std::string, unsigned int>::iterator iLname;

		//map between led name and its gpIO pin
		std::map<std::string, unsigned char> mLED_chan;
		std::map<std::string, unsigned char>::iterator iLchan;

		//map between led name and its gpIO pin
		std::map<std::string, double> mLED_duty;
		std::map<std::string, double>::iterator iLduty;
};

#endif
