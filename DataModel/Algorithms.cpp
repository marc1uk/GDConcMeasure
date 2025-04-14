#include "Algorithms.h"
#include <stdio.h>
#include <sys/wait.h>
#include <locale>
#include <errno.h>
#include <string.h>
#include <sys/stat.h> // for stat() test to see if file or folder

int SystemCall(std::string cmd, std::string& retstring){
  // execute command, capture return status and any output
  
  int retstatus;                     // exit code of command
  retstring="";
  
  // gotcha: our particular command redirects stdout,
  // so if we append 2>&1 we get stderr into our output file too!
  // fortunately we can put redirects whereever we like,
  // so insert it before anything else
  cmd.insert(0,"2>&1; ");
  cmd.append(";");
  
  // fork, execute cmd with '/bin/sh -c' and open a read pipe connected to outout (r)
  FILE * stream = popen(cmd.c_str(), "r");
  // check the fork succeeded and we got a pipe to read from
  if(stream){
    // read any return
    unsigned int bufsize=255;
    char buffer[bufsize]; // temporary buffer for reading return
    while(!feof(stream)){
      // if we're very paranoid we can check for read error with (ferror(stream)!=0)
      if (fgets(buffer, bufsize, stream) != NULL) retstring.append(buffer);
    }
    // pop off trailing newline
    if(retstring.size() && std::isspace(retstring.back())) retstring.pop_back();
    // close the pipe, and capture command return value
    int stat = pclose(stream);
    retstatus = (WIFEXITED(stat)) ? WEXITSTATUS(stat) : WTERMSIG(stat);
  } else {
    retstring = "SystemCall popen returned nullptr for command '"+cmd
               +"'', error is: "+strerror(errno);
    return -1;
  }
  if(retstatus!=0){
    retstring="SystemCall with command "+cmd+" failed with error code "
        +std::to_string(retstatus)+", stderr returned '"+retstring+"'";
  }
  
  return retstatus;
  
}


bool CheckPath(std::string path, std::string& type){
	struct stat s;
	if(stat(path.c_str(),&s)==0){
		if(s.st_mode & S_IFDIR){        // mask to extract if it's a directory?? how does this work?
			type="d";  //it's a directory
			return true;
		} else if(s.st_mode & S_IFREG){ // mask to check if it's a file??
			type="f"; //it's a file
			return true;
		} else {
			// exists, but neither file nor directory?
			type="???";
			return false;
			//assert(false&&"Check input path: stat says it's neither file nor directory..?");
		}
	} else {
		// does not exist - could be a pattern, e.g. "/path/to/rootfiles_*.root"
		type="none";
		return false;
	}
	return false;
}
