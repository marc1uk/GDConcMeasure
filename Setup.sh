#!/bin/bash

#Application path location of applicaiton (i.e. this script)
ToolDAQapp="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ -d ${ToolDAQapp}/ToolDAQ/root-6.28.06 ]; then
	. ${ToolDAQapp}/ToolDAQ/root-6.28.06/install/bin/thisroot.sh
else
	. ${ToolDAQapp}/ToolDAQ/root-6.14.06/install/bin/thisroot.sh
fi

export LD_LIBRARY_PATH=${ToolDAQapp}/lib:${ToolDAQapp}/ToolDAQ/zeromq-4.0.7/lib:${ToolDAQapp}/ToolDAQ/boost_1_66_0/install/lib:${ToolDAQapp}/ToolDAQ/seabreeze-3.0.11/SeaBreeze/lib:${ToolDAQapp}/ToolDAQ/libpqxx-6.4.7/install/lib:$LD_LIBRARY_PATH

if [ -f /var/lib/postgresql/setup_db.sh ]; then
	. /var/lib/postgresql/setup_db.sh
else
	. ${ToolDAQapp}/setup_db.sh
fi
