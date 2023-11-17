#!/bin/bash
. /home/pi/GDConcMeasure/setup_db.sh

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

# to monitor whether the GAD toolchain is actually turning the valves on and off,
# intermittently check the valve states, and retain the last 500 readings in the database
# first check we have these files, as they do need to be remade on boot so may not yet exist
INVALVESTATUS="$(cat /sys/class/gpio/gpio18/value 2> /dev/null)"
OUTVALVESTATUS="$(cat /sys/class/gpio/gpio15/value 2> /dev/null)"
PUMPSTATUS="$(cat /sys/class/gpio/gpio17/value 2> /dev/null)"
POWERSTATUS="$(cat /sys/class/gpio/gpio4/value 2>/dev/null)"
JSON="{ \"invalve\":${INVALVESTATUS}, \"outvalve\":${OUTVALVESTATUS}, \"pump\":${PUMPSTATUS}, \"power\":${POWERSTATUS} }"
#echo "JSON is '${JSON}'"
psql -U postgres -d rundb -c "INSERT INTO webpage ( name, timestamp, values ) VALUES ( 'gpio_status', 'NOW()', '${JSON}' )"

# the valves currently stay open for long enough to perform 5 dark measurements - about 4 minutes
# we need to sample at least once every 4 minutes to catch this happening, and properly
# recognise that we are indeed regularly turning the valves on and off.
# let's actually do it once a minute, for fidelity. In that case to store the last 24 hours
# we need (24*60) = 1440 samples. Samples older than this we'll delete.
psql -U postgres -d rundb -c "DELETE FROM webpage WHERE name='gpio_status' AND timestamp < now()-'24 hours'::interval"

#####################################################
# uncomment this to disable the valve safety check! #
#####################################################
#exit 0

# while we're doing this we can do some safety checks -
# check that the valves have not been open for the last 30 consecutive measurements
let INSOPEN=$(psql -d rundb -U postgres -At -c "SELECT SUM((values->'invalve')::integer) FROM webpage WHERE name='gpio_status' AND timestamp>(now() - '15 minutes'::interval)")
let OUTSOPEN=$(psql -d rundb -U postgres -At -c "SELECT SUM((values->'outvalve')::integer) FROM webpage WHERE name='gpio_status' AND timestamp>(now() - '15 minutes'::interval)")
if [ ${INSOPEN} -eq 15 ]; then
	echo "in valve has been open for the last 15 consecutive minutes! Closing valve!"
	curl -X POST -H 'Content-type: application/json' --data '{"text":" :warning: :warning: :warning: GAD INPUT VALVE HAS BEEN OPEN FOR THE LAST 15 CONSECUTIVE MINUTES! CLOSING THE VALVE! :warning: :warning: :warning:"}' ${SAFETYNETWEBHOOK}
	echo "0" > /sys/class/gpio/gpio18/value
fi
if [ ${OUTSOPEN} -eq 15 ]; then
	echo "out valve has been open for the last 15 consecutive minutes! Closing valve!"
	curl -X POST -H 'Content-type: application/json' --data '{"text":" :warning: :warning: :warning: GAD OUTPUT VALVE HAS BEEN OPEN FOR THE LAST 15 CONSECUTIVE MINUTES! CLOSING THE VALVE! :warning: :warning: :warning:"}' ${SAFETYNETWEBHOOK}
	echo "0" > /sys/class/gpio/gpio15/value
fi
 
 
