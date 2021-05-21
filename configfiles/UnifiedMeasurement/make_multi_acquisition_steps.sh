#!/bin/bash
# one initial set before any mixing, should be pure DI at this point
echo "measure Dark" > measurement_commands
echo "measure 275" >> measurement_commands
echo "save mixing_settling_tests_0.root" >> measurement_commands
for i in `seq 1 10`; do
	# each set of single_acquisition_steps.txt runs the pump for 60 seconds
	# then takes 17 pairs of measurements, plotting the transmission over around 30 mins
	echo "# acquisition $i" >> measurement_commands
	cat single_acquisition_steps.txt >> measurement_commands
	echo "save mixing_settling_tests_$i.root" >> measurement_commands
done
# after 10 such measurements, we'll have done 10 minutes of mixing and 300 minutes of waiting
# quit, ensuring everything is turned off and powered down
echo "quit" >> measurement_commands
