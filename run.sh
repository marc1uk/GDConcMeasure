#!/bin/bash
sudo /bin/bash -c "cd /home/pi/GDConcMeasure && . Setup.sh && ./GAD_ToolChain configfiles/${1:-UnifiedMeasurementDB}/ToolChainConfig"
