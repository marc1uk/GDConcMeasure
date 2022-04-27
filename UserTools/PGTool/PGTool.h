/* vim:set noexpandtab tabstop=4 wrap filetype=cpp */
#ifndef PGTool_H
#define PGTool_H

#include <string>
#include <iostream>

#include "Tool.h"

/**
* \class PGTool
*
* This is a demonstratoin Tool to show how to retrieve configuration variables from the Run Database. *
* $Author: M.O'Flaherty $
*/
class PGTool: public Tool {
	public:
	bool Initialise(std::string configfile,DataModel &data); ///< Initialise Function for setting up Tool resources. @param configfile The path and name of the dynamic configuration file to read in. @param data A reference to the transient data class used to pass information between Tools.
	bool Execute(); ///< Execute function used to perform Tool purpose.
	bool Finalise(); ///< Finalise function used to clean up resources.
	
	private:
	std::string toolName;
	
	// verbosity levels: if 'verbosity' < this level, the message type will be logged.
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage;
	int get_ok;
};


#endif
