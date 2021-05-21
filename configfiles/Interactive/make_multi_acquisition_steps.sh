#!/bin/bash
cat single_acquisition_steps.txt > multi_acquisition_steps.txt
echo -n '{"Save":"Save","Filename":"' >> multi_acquisition_steps.txt
echo -n "before_mix_0.root" >> multi_acquisition_steps.txt
echo '"}' >> multi_acquisition_steps.txt
echo '{"Valve":"OPEN"}' >> multi_acquisition_steps.txt
echo '60' >> multi_acquisition_steps.txt
echo '{"Valve":"CLOSE"}' >> multi_acquisition_steps.txt
echo '1800' >> multi_acquisition_steps.txt
for i in `seq 1 120`; do
	cat single_acquisition_steps.txt >> multi_acquisition_steps.txt
	echo -n '{"Save":"Save","Filename":"' >> multi_acquisition_steps.txt
	echo -n "mix_${i}.root" >> multi_acquisition_steps.txt
	echo '"}' >> multi_acquisition_steps.txt
	echo '{"Valve":"OPEN"}' >> multi_acquisition_steps.txt
	echo '60' >> multi_acquisition_steps.txt
	echo '{"Valve":"CLOSE"}' >> multi_acquisition_steps.txt
	echo '1800' >> multi_acquisition_steps.txt
done
