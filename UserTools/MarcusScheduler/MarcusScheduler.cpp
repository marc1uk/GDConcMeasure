#include "MarcusScheduler.h"

MarcusScheduler::MarcusScheduler():Tool(){}

bool MarcusScheduler::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	
	if(!m_variables.Get("command_file",command_file)){
		std::cout<<"no command_file specified in config file!"<<std::endl;
		return false;
	}
	m_variables.Get("overwrite_saves",overwrite_saves);
	
	std::string line;
	std::ifstream myfile (command_file);
	if (myfile.is_open()){
		while(getline(myfile,line)){
			if(line.size()==0) continue;
			if(line[0]=='#') continue;
			commands.push_back(line);
			if(line.substr(0,line.find(' '))=="loop") looping=true; // note this for later
		}
		myfile.close();
	} else {
		std::cout << "Unable to open file"<<command_file<<std::endl;
		return false;
	}
	
	current_command=-1;
	command_step=0;
	
	return true;
}


bool MarcusScheduler::Execute(){
	
	// startup, do not start the command process until the user types 'Start'
	try{
	if(current_command==-1){
		std::string command;
		std::cout<<std::endl<<"Type Start to power up, Stop to power down, or anything else to begin automation: ";
		std::cin>>command;
		if(command=="Start"){
			// bring power up
			std::string tmp("{\"Auto\":\"Stop\",\"Power\":\"ON\"}");
			m_data->CStore.JsonParser(tmp);
		} else if(command=="Stop"){
			// ensure power is off
			std::string tmp("{\"Auto\":\"Stop\",\"Power\":\"OFF\"}");
			m_data->CStore.JsonParser(tmp);
		} else if(command=="Quit"){
			std::string tmp("{\"Auto\":\"Stop\",\"Power\":\"OFF\"}");
			m_data->CStore.JsonParser(tmp);
			m_data->vars.Set("StopLoop",1);
		} else if(command!=""){
			// on anything else, load the first command to begin automation
			current_command=0;
		}
	}
	
	// if we've performed all the automation steps, go back to the beginning
	else if(current_command==commands.size()){
		std::cout<<"Reached end of command list, returning to main menu"<<std::endl;
		std::stringstream tmp;
		tmp<<"{\"Auto\":\"Stop\"}";
		m_data->CStore.JsonParser(tmp.str());
		current_command=-1;
	}
	// otherwise, we're doing an automation step!
	else {
		// if we've completed the last command, move to the next one
		std::string the_command = commands.at(current_command);
		std::cout<<"Processing command "<<current_command<<" step "<<command_step<<std::endl;
		std::cout<<"current command is "<<the_command<<std::endl;
		
		// to allow infinite loops we support the presence of a 'loop' command
		// but we would also like a mechanism to break such a loop
		// we'll do this by checking for existence of a flag file, i guess
		if(check_break_loop()){
			// hold off until we've finished our current measurement
			// if we're in the middle of one
			if(command_step==0){
				// once we're no longer in the middle of a measurement,
				// search for the loop point
				std::string a_command=the_command;
				int loop_point=current_command;
				// scan forward until we find the 'loop' command, or the end of the command chain
				while(loop_point<(commands.size()-1) && a_command.substr(0,a_command.find(' '))!="loop"){
					a_command = commands.at(++loop_point);
				}
				if(loop_point==(commands.size()-1)){
					if(a_command.substr(0,a_command.find(' '))=="loop"){
						// if the loop command is the last one in the command file,
						// return to the main menu
						std::cout<<"No further actions beyond loop point. "
							 <<"Returning to Main Menu"<<std::endl;
				                std::stringstream tmp;
				                tmp<<"{\"Auto\":\"Stop\"}";
				                m_data->CStore.JsonParser(tmp.str());
				                current_command=-1;
						the_command="dummy"; // avoid triggering anything else in this Execute()
					} else {
						// if we didn't find a 'loop' command, ignore the flag file
						std::cout<<"Did not find a 'loop' command, ignoring break loop"<<std::endl;
						std::string syscmd="rm -f "+break_loop_flagfile_name;
						system(syscmd.c_str());
					}
				// if we found a 'loop' command and there are other commands after it,
				// break the loop and carry on from the next command
				} else {
					std::cout<<"loop point at command "<<loop_point
						 <<" of "<<commands.size()<<" commands, "
						 <<"advancing to command "<<++loop_point<<std::endl;
					current_command=loop_point;
					std::string syscmd="rm -f "+break_loop_flagfile_name;
					system(syscmd.c_str());
					the_command="dummy"; // avoid triggering anything else in this Execute()
				}
			} // else hold off until the current measurement completes
		}

		// parse the command
		// if it's a delay
		if(the_command.substr(0,4)=="wait"){
			std::cout<<"Performing wait"<<std::endl;
			// strip off the wait
			std::string wait_time = the_command.substr(5,std::string::npos);
			// strip off any trailing comment
			wait_time = wait_time.substr(0,wait_time.find('#'));
			WaitForDuration(wait_time);
			++current_command;
		// if it's a measurement
		} else if(the_command.substr(0,7)=="measure"){
			std::cout<<"Performing Measurement step "<<command_step<<std::endl;
			// this is needed by the TraceAverage tool
			m_data->measurment_time= boost::posix_time::ptime(boost::posix_time::second_clock::local_time());
			std::string json_string="";
			switch (command_step){
				case 0:
					{
					// first thing, ensure all LEDs are off
					// ====================================
					json_string = "{\"R\":\"0\",\"G\":\"0\",\"B\":\"0\",\"White\":\"0\",\"385\":\"0\",\"260\":\"0\",\"275\":\"0\",\"LED\":\"Change\"}";
					++command_step;
					//break;  // XXX
					}
				case 1:
					{
					// now take a dark measurement
					// ===========================
					// if we have multiple measurements with a dark each, they need to be named accordingly
					// parse the set of LEDs in this measurement to append to the trace name
					std::string led_list = the_command.substr(8,the_command.find('#')-8);
					while(led_list.back()==' ') led_list.pop_back();              // trim trailing spaces
					std::replace( led_list.begin(), led_list.end(), ' ', '_' );   // replace spaces with underscores
					measurement_name = "Dark_"+led_list;
					json_string = "{\"Measure\":\"Start\",\"Trace\":\""+measurement_name+"\"}";
					++command_step;
					//break;  // XXX
					}
				case 2:
					// XXX taking a dark for _every_ measurement may be excessive. Disable the 'break' on case 0 and 1 to let it run straight through to 2.
					// all these do is set json_string and measurement_name, the former is overwritten and the latter is useable.
					{
					// then turn on the desired LEDs
					// =============================
					// these should be listed in order after the 'measure' string
					std::cout<<"will parse led list"<<std::endl;
					// initialize LED_states to all off
					LED_states = off_LED_states;
					// obtain the list of LEDs starting after 'measure' and ending at first '#' to string any trailing comments
					std::string led_list = the_command.substr(8,the_command.find('#')-8);
					while(led_list.back()==' ') led_list.pop_back();
					// split up the list of LEDs and enable all those that appear and are known
					std::cout<<"parsing LED list: \""<<led_list<<"\""<<std::endl;
					int last_space=0;
					while(true){
						int next_space = led_list.find(' ',last_space);
						std::string next_led_name = led_list.substr(last_space, next_space-last_space);
						std::cout<<"Setting LED "<<next_led_name<<" to ON for next measurement"<<std::endl;
						// check this is a valid LED name
						// compatibility
						if(LED_states.count(next_led_name)){
							LED_states.at(next_led_name) = 1; // enable this LED
						} else if(next_led_name!="Dark" && next_led_name!="None"){
							std::cerr<<"Unknown LED "<<next_led_name<<std::endl;
						}
						if(next_space!=std::string::npos){
							last_space = led_list.find_first_not_of(' ',next_space+1);
						} else {
							break;
						}
					}
					std::cout<<"building measurement name"<<std::endl;
					// add the LED states to the JSON string
					json_string="{";
					for(const std::pair<std::string,int>& an_led_state : LED_states){
						json_string += "\""+an_led_state.first+"\":\""+an_led_state.second+"\",";
					}
					measurement_name = measurement_name.substr(5,std::string::npos);  // use same measurement name as for Dark, but remove initial 'Dark_'
					// add the LED change command
					json_string += "\"LED\":\"Change\"}";
					++command_step;
					break;
					}
				case 3:
					{
					// now take the light measurement
					// ==============================
					json_string = "{\"Measure\":\"Start\",\"Trace\":\""+measurement_name+"\"}";
					++command_step;
					break;
					}
				case 4:
					{
					// finally, disable the LEDs. VERY IMPORTANT.
					// ==========================================
					json_string = "{\"R\":\"0\",\"G\":\"0\",\"B\":\"0\",\"White\":\"0\",\"385\":\"0\",\"260\":\"0\",\"275\":\"0\",\"LED\":\"Change\"}";
					// time to move to the next command
					command_step=0;
					++current_command;
					break;
					}
				default:
					std::cerr<<"measurement process step "<<command_step<<" unknown! How did we get here?!"<<std::endl;
			}
			std::cout<<"Queuing action: "<<json_string<<std::endl;
			// queue the action for this step of the measurement process
			m_data->CStore.JsonParser(json_string);
		} else if(the_command.substr(0,7) == "analyse"){
		  std::cout << "Analysing measurements\n";
		  std::string ledToAnalyse = the_command.substr(8, std::string::npos);
		  std::string json_string = "{\"Analyse\":\"Analyse\",\"ledToAnalyse\":\"" + ledToAnalyse + "\"}";
		  std::cout<<"Queuing action: "<<json_string<<'\n';
		  m_data->CStore.JsonParser(json_string);	 
		  ++current_command;
		} else if(the_command.substr(0,4)=="save"){
		        std::cout<<"Saving measurements"<<std::endl;
			// get the filename to save it as
			std::string filename = the_command.substr(5,std::string::npos);
			// trim comments & trailing whitespace
			filename = filename.substr(0,filename.find('#'));
			filename = filename.substr(0,filename.find(' '));
			std::string json_string;
			if(!looping){
				json_string = "{\"Save\":\"Save\",\"Filename\":\""+filename+"\"}";
			} else {
				std::string filenamesub = filename.substr(0,filename.find_last_of('.'));
				std::string filename_indexed = filenamesub+"_"+std::to_string(loop_count)+".root";
				json_string = "{\"Save\":\"Save\",\"Filename\":\""+filename_indexed+"\"}";
			}
			//std::string json_string = "{\"Save\":\"Save\",\"Filename\":\""+filename
			//	+"\,\"Overwrite\",\""+std::to_string(overwrite_saves)+"\"}";
			
			std::cout<<"Queuing action: "<<json_string<<std::endl;
			// queue the action for this step of the measurement process
			m_data->CStore.JsonParser(json_string);
			++current_command;
		} else if(the_command.substr(0,5)=="valve"){
			std::string open_or_close = the_command.substr(6,std::string::npos);
			open_or_close = open_or_close.substr(0,open_or_close.find(' '));
			std::cout<<"valve open or close is: \""<<open_or_close<<"\""<<std::endl;
			std::string json_string;
			if(open_or_close=="open"){
				json_string="{\"Valve\":\"OPEN\"}";
			} else if(open_or_close=="close"){
				json_string="{\"Valve\":\"CLOSE\"}";
			} else {
				std::cerr<<"Unknown valve state: '"<<open_or_close<<"'"<<std::endl;
			}
			m_data->CStore.JsonParser(json_string);
			++current_command;
		} else if(the_command.substr(0,the_command.find(' '))=="quit"){
                        std::string json_string("{\"Auto\":\"Stop\",\"Power\":\"OFF\"}");
                        m_data->CStore.JsonParser(json_string);
                        m_data->vars.Set("StopLoop",1);
		} else if(the_command.substr(0,the_command.find(' '))=="loop"){
			// for command files containing a series of steps that should be infinitely repeated
			++loop_count;
			current_command=0;
		} else if(the_command.substr(0,the_command.find(' '))=="dummy"){
			// dummy command, do nothing
		} else if(the_command[0]=='{'){
			// handle json style manual commands i guess
			m_data->CStore.JsonParser(the_command);
			++current_command;
		} else {
			std::cout<<"ERROR: Unknown command action: "<<the_command<<std::endl;
			++current_command;
		}
	}
	}
	catch(std::exception& e){
		std::cerr<<"MarcusScheduler caught exception "<<e.what()<<" during command "
			 <<current_command<<" of "<<commands.size()<<" at step "<<command_step<<std::endl;
		if(current_command<commands.size()){
			std::cerr<<"Corresponding command is: "<<commands.at(current_command)<<std::endl;
		}
		std::cerr<<"Aborting this command!"<<std::endl;
		// ensure system is in a safe state - turn off all LEDs, close the pump valves
		std::string json_string = "{\"R\":\"0\",\"G\":\"0\",\"B\":\"0\",\"White\":\"0\",\"385\":\"0\",\"260\":\"0\",\"275\":\"0\",\"LED\":\"Change\",\"Valve\":\"CLOSE\"}";
		command_step=0;
		++current_command;
	}
	
	// for the GracefulStop tool to be able to terminate infinite loops without interrupting
	// a mutli-loop measurement process, it needs to know when we are in the middle of a measurement
	m_data->CStore.Set("command_step",std::to_string(command_step));
	
	return true;
	
}

void MarcusScheduler::WaitForDuration(std::string wait_string){
	printf("Waiting for %s seconds\n",wait_string.c_str());
	std::stringstream tmp(wait_string);
	long wait;
	tmp>>wait;
	m_period=boost::posix_time::seconds(wait);
	boost::posix_time::ptime current(boost::posix_time::second_clock::local_time());
	boost::posix_time::ptime last(boost::posix_time::second_clock::local_time());
	boost::posix_time::time_duration lapse(m_period - (current - last));
	
	char lapse_buffer[100];
	int string_length = wait_string.length();
	std::string string_length_string = std::to_string(string_length);
	std::string lapse_buffer_format = "%0"+string_length_string+"ld";
	printf(lapse_buffer_format.c_str(),lapse.total_seconds());
	fflush(stdout);
	while(!lapse.is_negative()){
		// format the remaining time as a fixed width string
		long remaining_secs = lapse.total_seconds();
		string_length = snprintf(lapse_buffer, 100, lapse_buffer_format.c_str(), remaining_secs);
		// erase the current lapsed time
		for(int i=0; i<string_length; ++i) printf("\b");
		// overwrite it with the new lapsed time
		printf(lapse_buffer_format.c_str(),remaining_secs);
		fflush(stdout);
		sleep(1);
		current=boost::posix_time::second_clock::local_time();
		lapse=boost::posix_time::time_duration(m_period - (current - last));
	}
}


bool MarcusScheduler::Finalise(){
	
	commands.clear();
	
	return true;
}

bool MarcusScheduler::check_break_loop(){
  struct stat buffer;
  return (stat (break_loop_flagfile_name.c_str(), &buffer) == 0);
}
