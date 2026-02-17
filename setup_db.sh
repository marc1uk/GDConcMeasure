#!/bin/bash
export PATH=/opt/postgresql11/install/bin:$PATH
export LD_LIBRARY_PATH=/opt/postgresql11/install/lib:$LD_LIBRARY_PATH
export PGUSER=postgres
export PG_COLOR=always
export PGDATABASE=rundb
#export PGHOST=/var/run/postgresql
export PGHOST=192.168.50.100
export PGPORT=5432
