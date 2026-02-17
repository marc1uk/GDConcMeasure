#!/bin/bash
echo "4" | sudo tee /sys/class/gpio/export
echo "15" | sudo tee /sys/class/gpio/export
echo "17" | sudo tee /sys/class/gpio/export
echo "18" | sudo tee /sys/class/gpio/export

echo "out" | sudo tee /sys/class/gpio/gpio4/direction
echo "out" | sudo tee /sys/class/gpio/gpio15/direction
echo "out" | sudo tee /sys/class/gpio/gpio17/direction
echo "out" | sudo tee /sys/class/gpio/gpio18/direction

echo "1" | sudo tee /sys/class/gpio/gpio4/value
echo "0" | sudo tee /sys/class/gpio/gpio18/value
echo "0" | sudo tee /sys/class/gpio/gpio17/value
echo "0" | sudo tee /sys/class/gpio/gpio15/value
