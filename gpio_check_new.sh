#!/bin/bash
. /home/pi/GDConcMeasure/Setup.sh

APPLICATION_NAME="GAD_ToolChain"
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# we can't store the webhook in public version control
SAFETYNETWEBHOOK=$(cat /home/pi/safety_net/safety_net_webhook.txt);

# sanity check
if [ -z "${SAFETYNETWEBHOOK}" ]; then
        echo "Webhook safety net not defined!"
        exit 1;
fi

OUTPINS=(4 15 17 18)
for PIN in "${OUTPINS[@]}"; do
	if [ ! -e "/sys/class/gpio/gpio${PIN}" ]; then
		sudo echo "${PIN}" > /sys/class/gpio/export
	fi
	if [ "$(cat /sys/class/gpio/gpio${PIN}/direction)" != "out" ]; then
		sudo echo "out" > /sys/class/gpio/gpio${PIN}/direction
	fi
done

# to monitor whether the valves are still open,
# intermittently check the valve states, and retain the last N readings in the database
# first check we have these files, as they do need to be remade on boot so may not yet exist
HOLDINGVALVESTATUS="$(cat /sys/class/gpio/gpio18/value 2> /dev/null)"
SWITCHINGVALVESTATUS="$(cat /sys/class/gpio/gpio15/value 2> /dev/null)"
PUMPSTATUS="$(cat /sys/class/gpio/gpio17/value 2> /dev/null)"
POWERSTATUS="$(cat /sys/class/gpio/gpio4/value 2>/dev/null)"
JSON="{ \"switching_valve\":${SWITCHINGVALVESTATUS}, \"holding_valve\":${HOLDINGVALVESTATUS}, \"pump\":${PUMPSTATUS}, \"power\":${POWERSTATUS} }"
#echo "JSON is '${JSON}'"
psql -U postgres -d rundb -c "INSERT INTO webpage ( name, timestamp, values ) VALUES ( 'gpio_status', 'NOW()', '${JSON}' )"

# the valves should remain open, with only the holding output (gpio 18) active
# let's sample once a minute, for fidelity. In that case to store the last 24 hours
# we need (24*60) = 1440 samples. Samples older than this we'll delete.
psql -U postgres -d rundb -c "DELETE FROM webpage WHERE name='gpio_status' AND timestamp < now()-'24 hours'::interval"

#####################################################
# uncomment this to disable the valve safety check! #
#####################################################
#exit 0

# safety check that the holding output is enabled
if [ ${POWERSTATUS} -eq 0 ] && [ ${HOLDINGVALVESTATUS} -ne 1 ]; then
	echo "valve holding output is not enabled!"
	curl -X POST -H 'Content-type: application/json' --data '{"text":" :warning: :warning: :warning: GAD VALVE HOLDING OUTPUT IS NOT OPEN! :warning: :warning: :warning:"}' ${SAFETYNETWEBHOOK}
	
	# should we open? perhaps not. let the user do it, just in case
	#echo "1" > /sys/class/gpio/gpio15/value
	#echo "1" > /sys/class/gpio/gpio18/value
	#sleep 1
	#echo "0" > /sys/class/gpio/gpio15/value
fi

# safety check that the switching valve is not open for more than 5 minutes on the trot
let MINSENABLED=$(psql -d rundb -U postgres -At -c "SELECT SUM((values->'switching_valve')::integer) FROM webpage WHERE name='gpio_status' AND timestamp>(now() - '6 minutes'::interval)")

if [ ${SWITCHINGVALVESTATUS} -eq 1 ] && [ ${MINSENABLED} -ge 4 ]; then
	echo "valve switching output has been enabled for more than 5 minutes! Switching to holding!"
	curl -X POST -H 'Content-type: application/json' --data '{"text":" :warning: :warning: :warning: GAD VALVE SWITCHING OUTPUT HAS BEEN ENABLED FOR THE LAST 5 CONSECUTIVE MINUTES! SWITCHING TO HOLDING! :warning: :warning: :warning:"}' ${SAFETYNETWEBHOOK}
	# switch to holding
	echo "1" > /sys/class/gpio/gpio18/value
	sleep 1
	echo "0" > /sys/class/gpio/gpio15/value
fi
 
 
