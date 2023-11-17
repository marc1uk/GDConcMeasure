#!/usr/bin/env bash

if [ ! -d /sys/class/gpio/gpio15 ]; then
	echo "15" > /sys/class/gpio/export
fi
if [ ! -d /sys/class/gpio/gpio15 ]; then
	echo "18" > /sys/class/gpio/export
fi

dir=$(cat /sys/class/gpio/gpio15/direction)
if [ "${dir}" != "out" ]; then
	echo "out" > /sys/class/gpio/gpio15/direction
fi
dir=$(cat /sys/class/gpio/gpio18/direction)
if [ "${dir}" != "out" ]; then
	echo "out" > /sys/class/gpio/gpio18/direction
fi

# to stagger the values, we should ensure that if power is off,
# that both valves are set to 0 before powering on
powerison=1
if [ ! -d /sys/class/gpio/gpio4 ]; then
	powerison=0
elif [ $(cat /sys/class/gpio/gpio4/direction) != "out" ]; then
	powerison=0
elif [ $(cat /sys/class/gpio/gpio4/value) != "0" ]; then
	powerison=0
fi

if [ $powerison -eq 0 ]; then
	echo "power does not appear to be on, therefore valves closed"
	echo "setting valve state to closed to prevent simultaneous closing on powerup"
	echo "0" > /sys/class/gpio/gpio15/value
	# just in case
	sleep 1
	echo "0" > /sys/class/gpio/gpio18/value
fi

# power must be on to power the values
echo "ensuring power is on"
/home/pi/poweron.sh

if grep -Fxq "0" /sys/class/gpio/gpio18/value && grep -Fxq "0" /sys/class/gpio/gpio15/value; then
    echo 1 > /sys/class/gpio/gpio18/value
    sleep 1
    echo 1 > /sys/class/gpio/gpio15/value
    echo "valves opened"
elif grep -Fxq "1" /sys/class/gpio/gpio18/value && grep -Fxq "1" /sys/class/gpio/gpio15/value; then
    echo 0 > /sys/class/gpio/gpio18/value
    sleep 1
    echo 0 > /sys/class/gpio/gpio15/value
    echo "valve closed"
else
    echo "valve state mismatch"
fi  
    
