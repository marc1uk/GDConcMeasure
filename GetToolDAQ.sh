#!/bin/bash

# Saves the directory where this script is located to the ToolDAQapp
# variable. This method isn't foolproof. See
# https://stackoverflow.com/a/246128/4081973 if you need something more robust
# for edge cases (e.g., you're calling the script using symlinks).
TOPDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

init=1
tooldaq=1
boostflag=1
zmq=1
final=1
rootflag=1
eigen=0
seabreeze=1
libpqxx=1
serialcomms=1

# basic essentials
sudo dnf install -y git gcc g++ cmake binutils

while [ ! $# -eq 0 ]
do
    case "$1" in
	--help | -h)
	    helpmenu
	    exit
	    ;;

	--with_root | -r)
	    echo "Installing ToolDAQ with root"
	    rootflag=1 
	    ;;
	
	--no_boost | -b)
            echo "Installing ToolDAQ without boost"
            boostflag=0
	    ;;
	
	--no_init )
	     echo "Installing ToolDAQ without creating ToolDAQ Folder"
	    init=0;
	    ;;

	--no_zmq )
            echo "Installing ToolDAQ without zmq"
            zmq=0
            ;;

	--no_tooldaq )
	    echo "Installing dependancies without ToolDAQ"
	    tooldaq=0
	    ;;

	--no_final )
            echo "Installing ToolDAQ without compiling ToolAnalysis"
            final=0
            ;;

	--no_serial )
	    echo "Installing without RS232 library"
	    serialcomms=0
	    ;;

	--ToolDAQ_ZMQ )
            echo "Installing ToolDAQ & ZMQ"
	    boostflag=0
	    rootflag=0
	    final=0
            ;;

	--Boost )
            echo "Installing Boost"
	    init=0
	    tooldaq=0
	    zmq=0
	    final=0
	    rootflag=0
            ;;

	--Root )
            echo "Installing Root"
	    init=0
	    tooldaq=0
	    boostflag=0
	    zmq=0
	    final=0
            ;;
	
	--Serial )
	    echo "installing RS232 library"
	    init=0
	    tooldaq=0
	    boostflag=0
	    zmq=0
	    final=0
	    serialcomms=1
	
	--Final )
            echo "Compiling ToolDAQ"
	    init=0
	    tooldaq=0
	    boostflag=0
	    rootflag=0
	    zmq=0
            ;;

    esac
    shift
done

if [ $init -eq 1 ]
then
    
    mkdir ${TOPDIR}/ToolDAQ
fi

if [ $tooldaq -eq 1 ]
then
    cd ${TOPDIR}/ToolDAQ
    #git clone https://github.com/ToolDAQ/ToolDAQFramework.git
    git clone https://github.com/marc1uk/ToolDAQFramework.git
    git checkout GAD
fi

if [ $zmq -eq 1 ]
then
    cd ${TOPDIR}/ToolDAQ
    git clone https://github.com/ToolDAQ/zeromq-4.0.7.git
    
    cd zeromq-4.0.7
    
    ./configure --prefix=`pwd`
    make -j8
    make install
    
    export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
    
    cd ../
fi

if [ $eigen -eq 1 ]
then
    cd ${TOPDIR}/ToolDAQ
    wget https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
    tar -zxf eigen-3.3.7.tar.gz
    rm eigen-3.3.7.tar.gz
fi

if [ $seabreeze -eq 1 ]
then
    sudo dnf --enablerepo=crb install libusb-devel
    cd ${TOPDIR}/ToolDAQ
    git clone git@github.com:GDconcentration/SeaBreeze.git seabreeze-3.0.11
    cd seabreeze-3.0.11/SeaBreeze
    make
fi

if [ $libpqxx -eq 1 ]
then
    sudo dnf install -y postgresql postgresql-server libpq
    sudo dnf --enablerepo=crb install libpq-devel
    cd ${TOPDIR}/ToolDAQ
    wget https://github.com/jtv/libpqxx/archive/refs/tags/6.4.7.tar.gz
    tar -xzf 6.4.7.tar.gz
    rm 6.4.7.tar.gz
    cd libpqxx-6.4.7
    mkdir install
    ./configure --disable-documentation --enable-shared --prefix=$PWD/install
    make
    make install
    cd ..
fi

if [ $boostflag -eq 1 ]
then
    
    cd ${TOPDIR}/ToolDAQ
    git clone https://github.com/ToolDAQ/boost_1_66_0.git
    #wget http://downloads.sourceforge.net/project/boost/boost/1.66.0/boost_1_66_0.tar.gz
    
    #tar zxf boost_1_66_0.tar.gz
    #rm -rf boost_1_66_0.tar.gz
     
    cd boost_1_66_0
    
    mkdir install
    
    ./bootstrap.sh --prefix=`pwd`/install/  > /dev/null 2>/dev/null
    ./b2 install iostreams
    
    export LD_LIBRARY_PATH=`pwd`/install/lib:$LD_LIBRARY_PATH
    cd ../
fi


if [ $rootflag -eq 1 ]
then
    sudo dnf install -y libX11-devel libXpm-devel libXft-devel libXext-devel python python-devel openssl-devel
    cd ${TOPDIR}/ToolDAQ
    wget https://root.cern.ch/download/root_v6.14.06.source.tar.gz
    tar zxf root_v6.14.06.source.tar.gz
    rm -rf root_v6.14.06.source.tar.gz
    cd root-6.14.06
    mkdir install
    cd install
    cmake ../ -Dcxx14=OFF -Dcxx11=ON -Dxml=ON -Dmt=ON -Dmathmore=ON -Dx11=ON -Dimt=ON -DCMAKE_BUILD_TYPE=RelWithDebInfo -Dbuiltin_gsl=ON -Dfftw3=ON -DCMAKE_SHARED_LINKER_FLAGS='-latomic'
    # note latomic flag required by gcc8!
    # https://root-forum.cern.ch/t/root-6-18-04-fails-to-compile-on-raspbian-10-with-solution/36514
    make -j4
    make install
    source bin/thisroot.sh
    
    cd ../
    
fi

if [ $serialflag -eq 1 ]; then
	cd ${TOPDIR}/ToolDAQ
	git clone git@github.com:marc1uk/SerialCpp.git
	cd SerialCpp
	make
	cd ../
fi

cd ${TOPDIR}

if [ $final -eq 1 ]
then
    
    echo "current directory"
    echo `pwd`
    make clean
    make 
    
    export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
fi
