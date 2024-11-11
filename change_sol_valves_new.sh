#!/usr/bin/env bash

if [ ! -d /sys/class/gpio/gpio15 ]; then
	echo "exporting 15"
	echo "15" > /sys/class/gpio/export
fi
if [ ! -d /sys/class/gpio/gpio18 ]; then
	echo "exporting 18"
	echo "18" > /sys/class/gpio/export
fi

dir=$(cat /sys/class/gpio/gpio15/direction)
if [ "${dir}" != "out" ]; then
	echo "setting direction 15"
	echo "out" > /sys/class/gpio/gpio15/direction
fi
dir=$(cat /sys/class/gpio/gpio18/direction)
if [ "${dir}" != "out" ]; then
	echo "setting direction 18"
	echo "out" > /sys/class/gpio/gpio18/direction
fi

# if power is off, set valves to 0 before powering on
powerison=1
if [ ! -d /sys/class/gpio/gpio4 ]; then
	powerison=0
	echo "power not exported"
elif [ $(cat /sys/class/gpio/gpio4/direction) != "out" ]; then
	echo "power not set to output"
	powerison=0
elif [ $(cat /sys/class/gpio/gpio4/value) != "0" ]; then
	echo "power state 0"
	powerison=0
fi

if [ $powerison -eq 0 ]; then
	echo "power does not appear to be on, therefore valves closed"
	# switching valve
	echo "0" > /sys/class/gpio/gpio15/value
	# holding valve
	echo "0" > /sys/class/gpio/gpio18/value
	
	# power must be on to power the values
	echo "powering on"
	ret=$(/home/pi/poweron.sh)
	echo "poweron.sh returned $ret"
	if [ ! ret == 0 ]; then
		"powerup error"
		exit 1
	fi
fi


if grep -Fxq "1" /sys/class/gpio/gpio18/value || grep -Fxq "1" /sys/class/gpio/gpio15/value; then
	echo "valves open, closing"
	echo 0 > /sys/class/gpio/gpio15/value
	echo 0 > /sys/class/gpio/gpio18/value
	echo "valves closed"
elif grep -Fxq "0" /sys/class/gpio/gpio18/value && grep -Fxq "0" /sys/class/gpio/gpio15/value; then
	echo "valves closed, switching..."
	# switching (gpio15)
	echo 1 > /sys/class/gpio/gpio15/value
	# holding (gpio18)
	echo 1 > /sys/class/gpio/gpio18/value
	echo "sleep 1..."
	sleep 1
	echo "switch to holding only"
	echo 0 > /sys/class/gpio/gpio15/value
	echo "valves opened"
fi  
