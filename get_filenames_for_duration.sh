#!/bin/bash
cd /home/gad/GDConcMeasure
. Setup.sh
if [ $# -lt 2 ]; then
	echo "usage: $0 [start time] [end time] [LED (275_A)] [CRITERIA]"
	echo "timestamps should be in format: YYYY-MM-DD (HH:MM:SS)"
	echo "(times are optional)"
	echo "[CRITERIA] may be any other psql criteria(s)"
else
	START="'"${1}"'"
	END="'"${2}"'"
	LED="'"${3:-275_A}"'"
	CRITERIA=''
	if [ $# -gt 3 ]; then
		for i in `seq 4 $#`; do
			CRITERIA='AND ${!i} '
		done
	fi
	# debug override
	CRITERIA='run>10000'
	psql -c "SELECT timestamp,values->'rawfile' FROM data WHERE name='rawfile' AND ledname=${LED} AND timestamp > '2023-12-22' AND timestamp < '2023-12-26' ${CRITERIA} ORDER BY timestamp ASC"
fi

