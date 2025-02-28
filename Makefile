ToolDAQPath=ToolDAQ
ToolFrameworkCore=$(ToolDAQPath)/ToolFrameworkCore
ToolDAQFramework=$(ToolDAQPath)/ToolDAQFramework
SOURCEDIR=`pwd`

CXXFLAGS=  -fPIC -Wpedantic -Wall -std=c++14 -Wno-psabi -Wno-comment #-Wno-unused -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept  -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef #-Werror -Wold-style-cast 
#-fsanitize=address -O1 -fno-omit-frame-pointer -g

ifeq ($(MAKECMDGOALS),debug)
CXXFLAGS+= -O0 -g -lSegFault -rdynamic -DDEBUG
else
CXXFLAGS+= -O3
endif

#ToolDAQFrameworkInclude = $(ToolDAQFrameworkDIR)/include
#ToolDAQFrameworkLib = $(ToolDAQFrameworkDIR)/lib

ZMQLib= -L $(ToolDAQPath)/zeromq-4.0.7/lib -lzmq
ZMQInclude= -isystem $(ToolDAQPath)/zeromq-4.0.7/include/

RootLib = `root-config --libs`
RootInclude = -isystem `root-config --incdir`

BoostLib= -L $(ToolDAQPath)/boost_1_66_0/install/lib -lboost_date_time -lboost_serialization -lboost_iostreams -lboost_system
BoostInclude= -isystem $(ToolDAQPath)/boost_1_66_0/install/include

# need to check postgres install path - /usr/pgsql-12/lib...?
PostgresLib= -L $(ToolDAQPath)/libpqxx-6.4.7/install/lib -lpqxx -lpq
PostgresInclude= -isystem $(ToolDAQPath)/libpqxx-6.4.7/install/include

# simple c++ serial comms library
SerialLib = -L $(ToolDAQPath)/SerialCpp -lserial
SerialInclude = -isystem $(ToolDAQPath)/SerialCpp

# pi wiring library for gpio
#WiringPiLib = -L $(ToolDAQPath)/WiringPi/wiringPi -lwiringPi
#WiringPiInclude = -I $(ToolDAQPath)/WiringPi/wiringPi

# spectrometer
SeaBreezeLib = -L $(ToolDAQPath)/seabreeze-3.0.11/SeaBreeze/lib/ -lseabreeze
SeaBreezeInclude = -isystem ToolDAQ/seabreeze-3.0.11/SeaBreeze/include

# 64-bit location
LIBUSBPATH:=/usr/lib/aarch64-linux-gnu
ifeq ($(wildcard $(LIBUSBPATH)/.),)
  # if not present try 32-bit location
  LIBUSBPATH:=/usr/lib/arm-linux-gnueabihf
endif

DataModelInclude = $(RootInclude) $(SerialInclude) $(PostgresInclude)
DataModelLib = $(RootLib) $(SerialLib) $(PostgresLib)

MyToolsInclude = $(SeaBreezeInclude) $(WiringPiInclude)
MyToolsLib = $(WiringPiLib) $(SeaBreezeLib) -L$(LIBUSBPATH) -lusb

Includes=-I $(ToolFrameworkCore)/include/ -I $(ToolDAQFramework)/include/ -I $(SOURCEDIR)/include/ $(ZMQInclude) $(BoostInclude) $(DataModelInclude)
ToolLibraries = $(patsubst %, lib/%, $(filter lib%, $(subst /, , $(wildcard UserTools/*/*.so))))
LIBRARIES=lib/libDataModel.so lib/libMyTools.so $(ToolLibraries)
DataModelHEADERS:=$(patsubst %.h, include/%.h, $(filter %.h, $(subst /, ,$(wildcard DataModel/*.h))))
MyToolHEADERS:=$(patsubst %.h, include/%.h, $(filter %.h, $(subst /, ,$(wildcard UserTools/*/*.h) $(wildcard UserTools/*.h))))
ToolLibs = $(patsubst %.so, %, $(patsubst lib%, -l%,$(filter lib%, $(subst /, , $(wildcard UserTools/*/*.so)))))
AlreadyCompiled = $(wildcard UserTools/$(filter-out %.so UserTools , $(subst /, ,$(wildcard UserTools/*/*.so)))/*.cpp)
SOURCEFILES:=$(patsubst %.cpp, %.o,  $(filter-out $(AlreadyCompiled), $(wildcard */*.cpp) $(wildcard */*/*.cpp)))
Libs=-L $(SOURCEDIR)/lib/ -lDataModel -L $(ToolDAQFramework)/lib/ -lToolDAQChain -lDAQDataModelBase  -lDAQLogging -lServiceDiscovery -lDAQStore -L $(ToolFrameworkCore)/lib/ -lToolChain -lMyTools -lDataModelBase -lLogging -lStore -lpthread  $(ToolLibs) -L $(ToolDAQFramework)/lib/ -lToolDAQChain -lDAQDataModelBase  -lDAQLogging -lServiceDiscovery -lDAQStore $(ZMQLib) $(BoostLib)

#.SECONDARY: $(%.o)

all: $(DataModelHEADERS) $(MyToolHEADERS) $(SOURCEFILES) $(LIBRARIES) GAD_ToolChain NodeDaemon RemoteControl

debug: all

GAD_ToolChain: src/main.o $(LIBRARIES) $(DataModelHEADERS) $(MyToolHEADERS) | $(SOURCEFILES)
	@echo -e "\e[38;5;11m\n*************** Making " $@ " ****************\e[0m"
	g++  $(CXXFLAGS) $< -o $@ $(Includes) $(Libs) $(DataModelInclude) $(DataModelLib) $(MyToolsInclude) $(MyToolsLib) 

include/%.h:
	@echo -e "\e[38;5;87m\n*************** sym linking headers ****************\e[0m"
	ln -s  `pwd`/$(filter %$(strip $(patsubst include/%.h, /%.h, $@)), $(wildcard DataModel/*.h) $(wildcard UserTools/*/*.h) $(wildcard UserTools/*.h)) $@

src/%.o :  src/%.cpp   
	@echo -e "\e[38;5;214m\n*************** Making " $@ "****************\e[0m"
	g++ $(CXXFLAGS) -c $< -o $@ $(Includes)

UserTools/Factory/Factory.o :  UserTools/Factory/Factory.cpp  $(DataModelHEADERS) $(MyToolHEADERS)
	@echo -e "\e[38;5;214m\n*************** Making " $@ "****************\e[0m"
	g++ $(CXXFLAGS) -c $< -o $@ $(Includes) $(DataModelInclude) $(ToolsInclude)

UserTools/%.o :  UserTools/%.cpp  $(DataModelHEADERS) UserTools/%.h
	@echo -e "\e[38;5;214m\n*************** Making " $@ "****************\e[0m"
	g++ $(CXXFLAGS) -c $< -o $@ $(Includes) $(DataModelInclude) $(ToolsInclude)

DataModel/%.o : DataModel/%.cpp DataModel/%.h  $(DataModelHEADERS)
	@echo -e "\e[38;5;214m\n*************** Making " $@ "****************\e[0m"
	g++ $(CXXFLAGS) -c $< -o $@ $(Includes) $(DataModelInclude)

lib/libDataModel.so: $(patsubst %.cpp, %.o , $(wildcard DataModel/*.cpp)) |   $(DataModelHEADERS)
	@echo -e "\e[38;5;201m\n*************** Making " $@ "****************\e[0m"
	g++ $(CXXFLAGS) --shared $^ -o $@ $(Includes) $(DataModelInclude)

lib/libMyTools.so: $(patsubst %.cpp, %.o , $(filter-out $(AlreadyCompiled), $(wildcard UserTools/*/*.cpp))) |   $(DataModelHEADERS) $(MyToolHEADERS)
	@echo -e "\e[38;5;201m\n*************** Making " $@ "****************\e[0m"
	g++ $(CXXFLAGS) --shared $^ -o $@ $(Includes) $(DataModelInclude) $(MyToolsInclude)

lib/%.so:
	@echo -e "\e[38;5;87m\n*************** sym linking Tool libs ****************\e[0m"
	ln -s `pwd`/$(filter %$(strip $(patsubst lib/%.so, /%.so ,$@)), $(wildcard UserTools/*/*.so)) $@

NodeDaemon: $(ToolDAQFramework)/NodeDaemon
	@echo -e "\e[38;5;87m\n*************** sym linking " $@ " ****************\e[0m"
	ln -s $(ToolDAQFramework)/NodeDaemon ./

RemoteControl: $(ToolDAQFramework)/RemoteControl
	@echo -e "\e[38;5;87m\n*************** sym linking " $@ " ****************\e[0m"
	ln -s $(ToolDAQFramework)/RemoteControl ./

clean:
	@echo -e "\e[38;5;201m\n*************** Cleaning up ****************\e[0m"
	rm -f */*/*.o
	rm -f */*.o
	rm -f include/*.h
	rm -f lib/*.so
	rm -rf main
	rm -rf NodeDaemon
	rm -rf RemoteControl

make_pureref_DB_entry: make_pure_ref.cpp plotter.h plotter.cpp
	g++ -g -std=c++11 $< plotter.cpp -I./ `root-config --cflags --libs` -o $@

make_calcurve_DB_entry: make_cal_curve.cpp
	g++ -g -std=c++11 $^ `root-config --cflags --libs` -o $@
        
make_calcurve_DB_entry2: make_cal_curve2.cpp
	g++ -g -std=c++11 $^ `root-config --cflags --libs` -o $@
