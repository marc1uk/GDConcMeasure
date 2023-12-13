#include "LoadOldFiles.h"

LoadOldFiles::LoadOldFiles():Tool(){}

bool LoadOldFiles::Initialise(std::string configfile, DataModel &data){
	m_data= &data;
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_variables.Get("verbosity",verbosity);
	
	std::string file_list="";
	m_variables.Get("file_list",file_list);
	int max_files=-1;
	m_variables.Get("max_files",max_files);
	int n_files = ParseFileList(file_list, max_files);
	if(n_files==0){
		Log(m_unique_name+" No files in files list "+file_list,v_error,verbosity);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	
	next_file_iter = files.begin();
	
	// we need a temporary file because we can't define TTrees without putting them in a file...
	tmpfile = new TFile("/tmp/tmpfile.root","RECREATE");
	
	return true;
}


bool LoadOldFiles::Execute(){
	
	bool got_file=false;
	std::string filename;
	std::string treename;
	while(next_file_iter!=files.end()){
		filename = next_file_iter->filename;
		treename = next_file_iter->treename;
		int runnum = next_file_iter->runnum;
		Log("LoadOldFiles: processing file "+filename+", tree "+treename+", run "+std::to_string(runnum),v_debug,verbosity);
		
		// set run number to use in database
		m_data->CStore.Set("dbrunnum",runnum);
		
		++next_file_iter;
		if(next_file_iter==files.end()){
			Log("LoadOldFiles: last file, setting StopLoop",v_message,verbosity);
			m_data->vars.Set("StopLoop",1);
		}
		
		// try to open the file
		Log("LoadOldFiles: opening file",v_debug,verbosity);
		if(nextfile) nextfile->Close();
		if(darkTreeNew) delete darkTreeNew;
		nextfile = TFile::Open(filename.c_str(),"READ");
		if(nextfile==nullptr || nextfile->IsZombie()){
			Log(std::string("LoadOldFiles failed to open file '")+filename+"'",v_error,verbosity);
			continue;
		}
	
		// get the led-on tree
		Log("LoadOldFiles: getting ledtree",v_debug,verbosity);
		ledTree = (TTree*)nextfile->Get(treename.c_str());
		if(ledTree==nullptr){
			Log("LoadOldFiles failed to find specified tree '"+treename+"' in file '"+filename+"'",v_error,verbosity);
			continue;
		}
		
		// under normal operation, each measurement results in a TTree with one entry... probably not very efficient.
		// the MatthewAnalysis tool analyses led Tree entry 0, so warn if there are others
		Log("LoadOldFiles: checking one entry",v_debug,verbosity);
		if(ledTree->GetEntries()==0){
			Log("LoadOldFiles found no entries in tree '"+treename+"' in file '"+filename+"!",v_warning,verbosity);
			continue;
		} else if(ledTree->GetEntries()>1){
			Log("LoadOldFiles found more than one entry in tree '"+treename
				+"' in file '"+filename+"! Only entry 0 will be analysed!",v_warning,verbosity);
		}
		
		got_file=true;
		break;
	}
	
	// check we managed to load a file successfully
	if(!got_file) return false;
	
	// arbitrary measurement number to uniquely identify a measurement
	// we should probably retrieve this from the database when re-analysing old data
	// so that we can match the new results to the old results.... TODO
	++measurementnum;
	get_ok = m_data->CStore.Get("dbmeasurementnum",measurementnum);
	
	// while each LED gets saved to a different file, the dark traces for all LED measurements are currently saved
	// in one common file. It is also the case that several unused dark traces are taken between measurements
	// to warm up the spectrometer. The result of this is that the MatthewAnalysis takes the last (most recent)
	// entry from the Dark tree for its dark subtraction, but other dark traces may be taken after this.
	// When reading from the file, we therefore need to be careful about picking the right Dark trace.
	// To do this, we'll find the most recent dark by timestamp
	Log("LoadOldFiles: setting branch addresses",v_debug,verbosity);
	Short_t yr, mon, dy, hr, mn, sc;
	ledTree->SetBranchAddress("year",&yr);
	ledTree->SetBranchAddress("month",&mon);
	ledTree->SetBranchAddress("day",&dy);
	ledTree->SetBranchAddress("hour",&hr);
	ledTree->SetBranchAddress("min",&mn);
	ledTree->SetBranchAddress("sec",&sc);
	Log("LoadOldFiles: getting entry 0",v_debug,verbosity);
	ledTree->GetEntry(0);
	struct tm ledtime;
	ledtime.tm_year = yr - 1900;
	ledtime.tm_mon = mon - 1;
	ledtime.tm_mday = dy;
	ledtime.tm_hour = hr;
	ledtime.tm_min = mn;
	ledtime.tm_sec = sc;
	time_t ledtime_t = mktime(&ledtime);
	// for some reason on the first toolchain Execute (perhaps just the first mktime call)
	// this produces a timestamp one hour later than the passed ledtime, and even changes
	// the ledtime.tm_hour to be one hour later. re-set the hour and regenerate.
	ledtime.tm_hour = hr;
	ledtime_t = mktime(&ledtime);
	
	// we can override the default "now()" timestamp passed to postgres by providing it explicitly
	char dbtimestamp[20];
	snprintf(dbtimestamp, 20, "%04d-%02d-%02d %02d:%02d:%02d", yr, mon, dy, hr, mn, sc);
	std::string dbtimestampstring(dbtimestamp);
	m_data->CStore.Set("dbtimestamp",dbtimestampstring);
	
	// get the dark tree
	Log("LoadOldFiles: getting dark tree",v_debug,verbosity);
	darkTree = (TTree*)nextfile->Get("Dark");
	if(darkTree==nullptr){
		// hack because of inconsistencies in old data
		darkTree = (TTree*)nextfile->Get("dark");
	}
	if(darkTree==nullptr){
		Log(std::string("LoadOldFiles failed to find dark tree in file '")+filename+"'",v_error,verbosity);
		return false;
	}
	
	Log("LoadOldFiles: setting branch addresses",v_debug,verbosity);
	darkTree->SetBranchAddress("year",&yr);
	darkTree->SetBranchAddress("month",&mon);
	darkTree->SetBranchAddress("day",&dy);
	darkTree->SetBranchAddress("hour",&hr);
	darkTree->SetBranchAddress("min",&mn);
	darkTree->SetBranchAddress("sec",&sc);
	
	// loop over dark entries
	int darkentry=0;
	Log("LoadOldFiles: scanning "+std::to_string(darkTree->GetEntries())+" entries for closest timestamp",v_debug,verbosity);
	for(int i=0; i<darkTree->GetEntries(); ++i){
		darkTree->GetEntry(i);
		
		struct tm darktime;
		darktime.tm_year = yr - 1900;
		darktime.tm_mon = mon - 1;
		darktime.tm_mday = dy;
		darktime.tm_hour = hr;
		darktime.tm_min = mn;
		darktime.tm_sec = sc;
		time_t darktime_t = mktime(&darktime);
		
		// difftime does [ time_a - time_b ]
		double numsecs = difftime(ledtime_t, darktime_t);
		// do we use the closest in time dark, or the last dark before the led-on?
		// what if someone decides to do the dark measurements after the led?
		if(numsecs<0) break;
		darkentry=i;
	}
	Log("LoadOldFiles: dark entry number = "+std::to_string(darkentry),v_debug,verbosity);
	
	// make a new TTree with just the dark entry we're intersted in
	// because TTrees must attach to a file, attach it to a temporary file
	Log("LoadOldFiles: cloning tree",v_debug,verbosity);
	auto curdir = gDirectory;
	tmpfile->cd();
	// copy the branch structure and addresses from the input tree
	darkTreeNew = darkTree->CloneTree(0);
	if(darkTreeNew==nullptr){
		Log("LoadOldFiles: cloned tree is null!",v_error,verbosity);
		return false;
	}
	curdir->cd();
	
	// load the desired entry to copy
	Log("LoadOldFiles: getting desired entry",v_debug,verbosity);
	darkTree->GetEntry(darkentry);
	// write it to the clone tree
	Log("LoadOldFiles: filling copy",v_debug,verbosity);
	darkTreeNew->Fill();
	
	//m_data->m_trees is a std::map<std::string, TTree*>
	// we need to put the 'dark' tree and the 'led' tree in it,
	// where 'led' tree name must match ledToAnalyse
	Log("LoadOldFiles: setting names",v_debug,verbosity);
	m_data->m_trees["dark"] = darkTreeNew;
	m_data->m_trees[treename] = ledTree;
	
	std::pair<int,int> treeentrynums{0,darkentry}; // light then dark
	get_ok = m_data->CStore.Get("dbtreeentries",treeentrynums);
	
	m_data->CStore.Set("ledToAnalyse",treename);
	m_data->CStore.Set("Filename",filename);
	
	// trigger the MatthewAnalysis
	std::string analyse = "Analyse";
	m_data->CStore.Set("Analyse", analyse);
	
	return true;
}


bool LoadOldFiles::Finalise(){
	
	if(tmpfile){
		tmpfile->Close();
		delete tmpfile;
	}
	if(nextfile){
		nextfile->Close();
	}
	return true;
}

int LoadOldFiles::ParseFileList(std::string file_list, int max_files){
	std::ifstream files_list(file_list.c_str());
	if(!files_list.is_open()){
		Log(std::string("LoadOldFiles Failed to load file list '")+file_list+"'",v_error,verbosity);
		return 0;
	}
	std::string line;
	std::stringstream linestream;
	std::string filename;
	std::string treename;
	int runnum;
	oldfile anoldfile;
	while(getline(files_list,line)){
		if(line.empty()) continue;
		if(line[0]=='#') continue;
		linestream.clear();
		linestream.str(line);
		if(!(linestream >> anoldfile.filename >> anoldfile.treename >> anoldfile.runnum)){
			Log("LoadOldFiles failure parsing line "+line,v_error,verbosity);
			continue;
		}
		files.push_back(anoldfile);
		if(max_files>0 && files.size()==max_files) break;
	}
	files_list.close();
	
	return files.size();
}
