#!/bin/bash

CONTENTS="
#name           channel   voltage
LEDW            1         3.2
LEDRGBAnnode    12        2.0
LEDR            13        2.0
LEDG            14        2.0
LEDB            15        2.0
LED385L         8         3.5
#
# Note: the PMW output to the 275nm LEDs is upscaled from 5V to 12V line voltage,
# while maintaining the same duty cycle. So the actual voltage applied to the LEDs
# is the voltage given here multiplied by (12/5).
LED275J_A       7         3.3
LED275J_B       11        3.3'
"

psql -c "INSERT INTO auxiliary ( created, name, version, description, author, contents ) VALUES ( 'now()', 'LEDwiring', 0, 'first EGADS set with database recording', 'Marcus OFlaherty', '${CONTENTS}' );"

CONTENTS="
start
power on
wait 900
start_loop
setFile 05SeptEGADS
valve outlet open
valve inlet open
measure Dark
measure Dark
measure Dark
measure Dark
measure Dark
valve inlet close
valve outlet close
wait 30
measure Dark
measure Dark
measure Dark
measure Dark
measure Dark
measure 275_B
analyse 275_B
measure Dark
measure 275_A
analyse 275_A

measure Dark
measure R G B
measure Dark
measure White 385

transparency

save
loop
quit
"

psql -c "INSERT INTO auxiliary ( created, name, version, description, author, contents ) VALUES ( 'now()', 'measurement_commands', 1, 'first EGADS set with database recording', 'Marcus OFlaherty', '${CONTENTS}' );"

CONTENTS="
myGracefulStop         GracefulStop            0
myMarcusScheduler      MarcusScheduler         0
myBenPower             BenPower                0
myPump                 Valve                   0
myValveInlet           Valve                   1
myValveOutlet          Valve                   2
BenLED                 BenLED                  0
BenSpectrometer        BenSpectrometer         0
TraceAverage           TraceAverage            0
MarcusAnalysisLEDA     MarcusAnalysis          0
MarcusAnalysisLEDB     MarcusAnalysis          1
SaveTraces             SaveTraces              0
MatthewTransparency    MatthewTransparency     0
mySaveToDB             SaveToDB                0
"

psql -c "INSERT INTO runconfigs ( runconfig, created, description, author, toolsconfig ) VALUES ( 1, 'now()', 'first EGADS set with database recording', 'Marcus OFlaherty', '${CONTENTS}' )"

CONTENTS="
verbosity 1
ledToAnalyse 275_A
pureref_ver 0         # DB pure trace ID
#pureref_file ./pureDarkSubtracted_275_A_julcal2.root
# may have several methods for fitting absorption curve
# all will be used and their corresponding gd concentrations saved
Method_0 raw          # absorption graph fitting method
Calib_0 0             # DB calibration curve ID
Method_1 simple
Calib_1 0
Method_2 complex
Calib_2 0

# alternatively one may specify a local calibration curve
# for example (not an actual calibration curve)
#Calib_0 Local
#CalibFunc_0 pol6
#CalibNpar_0 7
#CalibPar_0_0 0.000325653
#CalibPar_0_1 2.9966
#CalibPar_0_2 -11.0959
#CalibPar_0_3 12.3771
#CalibPar_0_4 54.1008
#CalibPar_0_5 -183.526
#CalibPar_0_6 128.346

#Calib_1 Local
#CalibFunc_1 pol6
#CalibNpar_1 7
#CalibPar_1_0 0.000510926
#CalibPar_1_1 2.75227
#CalibPar_1_2 -12.2684
#CalibPar_1_3 21.3031
#CalibPar_1_4 59.946
#CalibPar_1_5 -385.552
#CalibPar_1_6 650.412

#Calib_2 Local
#CalibFunc_2 pol6
#CalibNpar_2 7
#CalibPar_2_0 0.000505721
#CalibPar_2_1 2.79712
#CalibPar_2_2 -12.3976
#CalibPar_2_3 19.9358
#CalibPar_2_4 61.6582
#CalibPar_2_5 -352.829
#CalibPar_2_6 563.938
"

psql -c "INSERT INTO configfiles ( tool, version, created, author, description, contents ) VALUES ( 'MarcusAnalysis', 0, 'now()', 'Marcus OFlaherty', 'First version, LED A', '${CONTENTS}' ); "

CONTENTS="
verbosity 1
ledToAnalyse 275_B
pureref_ver 0         # DB pure trace ID
#pureref_file ./pureDarkSubtracted_275_B_julcal2.root
# may have several methods for fitting absorption curve
# all will be used and their corresponding gd concentrations saved
Method_0 raw          # absorption graph fitting method
Calib_0 0             # DB calibration curve ID
Method_1 simple
Calib_1 0
Method_2 complex
Calib_2 0
"

psql -c "INSERT INTO configfiles ( tool, version, created, author, description, contents ) VALUES ( 'MarcusAnalysis', 1, 'now()', 'Marcus OFlaherty', 'First version, LED A', '${CONTENTS}' ); "


