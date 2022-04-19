#!/bin/bash

init=1
tooldaq=1
boostflag=1
zmq=1
final=1
rootflag=0
eigen=1
seabreeze=1
libpqxx=1

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
    
    mkdir ToolDAQ
fi

cd ToolDAQ

if [ $tooldaq -eq 1 ]
then

git clone https://github.com/ToolDAQ/ToolDAQFramework.git
fi

if [ $zmq -eq 1 ]
then
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
    wget https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
    tar -zxf eigen-3.3.7.tar.gz
    rm eigen-3.3.7.tar.gz
fi

if [ $seabreeze -eq 1 ]
then
    echo "please unpack seabreeze sourcefiles into ToolDAQ folder" # TODO add to git repo
fi

if [ $libpqxx -eq 1 ]
then
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

cd ../

if [ $final -eq 1 ]
then
    
    echo "current directory"
    echo `pwd`
    make clean
    make 
    
    export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
fi
