#ifndef PGHELPER_H
#define PGHELPER_H

#include "Postgres.h"
#include "Store.h"
#include "BoostStore.h"

using namespace ToolFramework;

class DataModel;

class PGHelper{
	public:
	PGHelper(DataModel* m_data_in=nullptr);
	~PGHelper(){m_data=nullptr;}
	void SetDataModel(DataModel* m_data_in);
	void SetVerbosity(int verb);
	
	bool GetCurrentRun(int& runnum, int* runconfig=nullptr, std::string* err=nullptr);
	bool GetRunConfig(int& runconfig_in, int* runnum_in=nullptr, std::string* err=nullptr);
	bool GetToolsConfig(std::string& toolsconfig, int* runnum=nullptr, int* runconfig_in=nullptr, std::string* err=nullptr);
	bool ParseToolsConfig(std::string toolsconfig, BoostStore& configStore, std::string* err); ///< parse the toolsconfig string into a Store where key is the unique toolname and value is a std::pair< std::string ToolClassName, int ToolConfigFileVersion >
	bool GetToolConfig(std::string className, int version, std::string& toolconfig, std::string* err); ///< Get configuration for a given tool given its class name and the configuration version
	bool GetToolConfig(std::string toolname, std::string& toolconfig, int* runconfig_in=nullptr, int* runnum_in=nullptr, std::string* err=nullptr); ///< Get the configuration for the given instance of a tool in a particular toolchain, under the unique name 'toolname'
	
	int InsertToolConfig(Store config, std::string toolname, std::string author, std::string description, std::string* err);
	
	private:
	DataModel* m_data;
	
	int verbosity=1;
	int v_error=0;
	int v_warning=1;
	int v_message=2;
	int v_debug=3;
	std::string logmessage;
	int get_ok;
	
};


#endif
