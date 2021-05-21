#!/bin/bash

for j in `seq 1 28`; do  #2week loop
    for i in `seq 1 48`; do   #12 hr loop
	cat single_acquisition_steps.txt >> multi_acquisition_steps_ben.txt
	echo 900 >> multi_acquisition_steps_ben.txt
    done
    echo -n '{"Save":"Save","Filename":"' >> multi_acquisition_steps_ben.txt
    echo -n "stability_${j}.root" >> multi_acquisition_steps_ben.txt
    echo '"}' >> multi_acquisition_steps_ben.txt
done
    
