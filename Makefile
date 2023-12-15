ToolDAQPath=ToolDAQ

CPPFLAGS= -Wno-psabi -std=c++14

ZMQLib= -L $(ToolDAQPath)/zeromq-4.0.7/lib -lzmq
ZMQInclude= -I $(ToolDAQPath)/zeromq-4.0.7/include/

RootLib = `root-config --libs`
RootInclude = `root-config --cflags`

BoostLib= -L $(ToolDAQPath)/boost_1_66_0/install/lib -lboost_date_time -lboost_serialization -lboost_iostreams -lboost_system
BoostInclude= -I $(ToolDAQPath)/boost_1_66_0/install/include

# need to check postgres install path - /usr/pgsql-12/lib...?
PostgresLib= -L $(ToolDAQPath)/libpqxx-6.4.7/install/lib -lpqxx -lpq
PostgresInclude= -I $(ToolDAQPath)/libpqxx-6.4.7/install/include

DataModelInclude= $(RootInclude) -I  ToolDAQ/eigen-3.3.7/
DataModelLib=  $(RootLib)

# 64-bit location
LIBUSBPATH:=/usr/lib/aarch64-linux-gnu
ifeq ($(wildcard $(LIBUSBPATH)/.),)
  # if not present try 32-bit location
  LIBUSBPATH:=/usr/lib/arm-linux-gnueabihf
endif

MyToolsInclude= $(RootInclude) -isystem ToolDAQ/seabreeze-3.0.11/SeaBreeze/include/ $(PostgresInclude) -I $(ToolDAQPath)/WiringPi/wiringPi/
MyToolsLib= $(RootLib) -L $(ToolDAQPath)/WiringPi/wiringPi/ -lwiringPi -L ToolDAQ/seabreeze-3.0.11/SeaBreeze/lib/ -lseabreeze $(PostgresLib) -L$(LIBUSBPATH) -lusb

all: lib/libStore.so lib/libLogging.so lib/libDataModel.so include/Tool.h lib/libMyTools.so lib/libServiceDiscovery.so lib/libToolChain.so GAD_ToolChain RemoteControl  NodeDaemon

GAD_ToolChain: src/main.cpp | lib/libMyTools.so lib/libStore.so lib/libLogging.so lib/libToolChain.so lib/libDataModel.so lib/libServiceDiscovery.so
	@echo -e "\n*************** Making " $@ "****************"
	g++ -g $(CPPFLAGS) src/main.cpp -o $@ -I include -L lib -lStore -lMyTools -lToolChain -lDataModel -lLogging -lServiceDiscovery -lpthread $(DataModelInclude) $(DataModelib) $(MyToolsInclude)  $(MyToolsLib) $(ZMQLib) $(ZMQInclude)  $(BoostLib) $(BoostInclude)


lib/libStore.so: $(ToolDAQPath)/ToolDAQFramework/src/Store/*
	cd $(ToolDAQPath)/ToolDAQFramework && make lib/libStore.so
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/src/Store/*.h include/
	cp $(ToolDAQPath)/ToolDAQFramework/lib/libStore.so lib/
	#g++ -g -fPIC -shared  -I include $(ToolDAQPath)/ToolDAQFramework/src/Store/*.cpp -o lib/libStore.so $(BoostLib) $(BoostInclude)


include/Tool.h:  $(ToolDAQPath)/ToolDAQFramework/src/Tool/Tool.h
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/src/Tool/Tool.h include/


lib/libToolChain.so: $(ToolDAQPath)/ToolDAQFramework/src/ToolChain/* | lib/libLogging.so lib/libStore.so lib/libMyTools.so lib/libServiceDiscovery.so lib/libLogging.so lib/libDataModel.so
	@echo -e "\n*************** Making " $@ "****************"
	#cp $(ToolDAQPath)/ToolDAQFramework/UserTools/Factory/*.h include/
	cp $(ToolDAQPath)/ToolDAQFramework/src/ToolChain/*.h include/
	g++ -g $(CPPFLAGS) -fPIC -shared $(ToolDAQPath)/ToolDAQFramework/src/ToolChain/ToolChain.cpp -I include -lpthread -L lib -lStore -lDataModel -lServiceDiscovery -lLogging -lMyTools -o lib/libToolChain.so $(DataModelInclude) $(DataModelib) $(ZMQLib) $(ZMQInclude) $(MyToolsInclude) $(MyToolsLib) $(BoostLib) $(BoostInclude)


clean:
	@echo -e "\n*************** Cleaning up ****************"
	rm -f include/*.h
	rm -f lib/*.so
	rm -f GAD_ToolChain
	rm -f RemoteControl
	rm -f NodeDaemon

lib/libDataModel.so: DataModel/* lib/libLogging.so | lib/libStore.so
	@echo -e "\n*************** Making " $@ "****************"
	cp DataModel/*.h include/
	g++ -g $(CPPFLAGS) -fPIC -shared DataModel/*.cpp -I include -L lib -lStore  -lLogging  -o lib/libDataModel.so $(DataModelInclude) $(DataModelLib) $(ZMQLib) $(ZMQInclude)  $(BoostLib) $(BoostInclude) $(PostgresInclude) $(PostgresLib)

lib/libMyTools.so: UserTools/*/* UserTools/* | include/Tool.h lib/libDataModel.so lib/libLogging.so lib/libStore.so
	@echo -e "\n*************** Making " $@ "****************"
	cp UserTools/*/*.h include/
	cp UserTools/Factory/*.h include/
	g++ -g $(CPPFLAGS) -fPIC -shared  UserTools/Factory/Factory.cpp -I include -L lib -lStore -lDataModel -lLogging -o lib/libMyTools.so $(MyToolsInclude) $(MyToolsLib) $(DataModelInclude) $(DataModelib) $(ZMQLib) $(ZMQInclude) $(BoostLib) $(BoostInclude)

RemoteControl:
	cd $(ToolDAQPath)/ToolDAQFramework/ && make RemoteControl
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/RemoteControl ./

NodeDaemon:
	cd $(ToolDAQPath)/ToolDAQFramework/ && make NodeDaemon
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/NodeDaemon ./

lib/libServiceDiscovery.so: $(ToolDAQPath)/ToolDAQFramework/src/ServiceDiscovery/* | lib/libStore.so
	cd $(ToolDAQPath)/ToolDAQFramework && make lib/libServiceDiscovery.so
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/src/ServiceDiscovery/ServiceDiscovery.h include/
	cp $(ToolDAQPath)/ToolDAQFramework/lib/libServiceDiscovery.so lib/
	#g++ -shared -fPIC -I include $(ToolDAQPath)/ToolDAQFramework/src/ServiceDiscovery/ServiceDiscovery.cpp -o lib/libServiceDiscovery.so -L lib/ -lStore  $(ZMQInclude) $(ZMQLib) $(BoostLib) $(BoostInclude)

lib/libLogging.so:  $(ToolDAQPath)/ToolDAQFramework/src/Logging/* | lib/libStore.so
	cd $(ToolDAQPath)/ToolDAQFramework && make lib/libLogging.so
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/src/Logging/Logging.h include/
	cp $(ToolDAQPath)/ToolDAQFramework/lib/libLogging.so lib/
	#g++ -shared -fPIC -I include $(ToolDAQPath)/ToolDAQFramework/src/Logging/Logging.cpp -o lib/libLogging.so -L lib/ -lStore $(ZMQInclude) $(ZMQLib) $(BoostLib) $(BoostInclude)

update:
	@echo -e "\n*************** Updating ****************"
	cd $(ToolDAQPath)/ToolDAQFramework; git pull
	cd $(ToolDAQPath)/zeromq-4.0.7; git pull
	git pull

make_pureref_DB_entry: make_pure_ref.cpp
	g++ -std=c++11 $^ `root-config --cflags --libs` -o $@

make_calcurve_DB_entry: make_cal_curve.cpp
	g++ -std=c++11 $^ `root-config --cflags --libs` -o $@
