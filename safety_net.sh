#! /usr/bin/env bash

APPLICATION_NAME="GAD_ToolChain"
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SAFETYNETWEBHOOK=$(cat ${SCRIPTDIR}/safety_net_webhook.txt);

# sanity check
if [ -z "${SAFETYNETWEBHOOK}" ]; then
	echo "Webhook safety net not defined!"
	exit 1;
fi

# check the GAD toolchain is running
# initial sanity check that the safety net control file exists and contains either 'ON' or 'OFF'
if ! grep -wq "OFF\|ON" ${SCRIPTDIR}/safety_net_active.txt; then
    curl -X POST -H 'Content-type: application/json' --data '{"text":" :warning: :warning: :warning: GAD SAFETY NET ACTIVE FLAG IS NEITHER ON OR OFF - CHECK `'"${SCRIPTDIR}/safety_net_active.txt"'` FILE :warning: :warning: :warning:"}' ${SAFETYNETWEBHOOK}
    sudo wall "safety_net_active.txt file should contain 'OFF' or 'ON'"

# assuming this is the case, if safety net is not intentionally disabled...
elif grep -wq "ON" ${SCRIPTDIR}/safety_net_active.txt; then
    # ... then try to find a running GAD_ToolChain process
    if [ $(ps aux | grep $APPLICATION_NAME | grep -v sudo | grep -v gdb | grep -v grep | grep -v defunct | wc -l) -lt 1 ]; then
        # if no GAD_Toolchain process was found, ensure the valves are disabled
        # first check the power file exists; without this the valves cannot be enabled
        if [ "$(cat /sys/class/gpio/gpio4/value 2>/dev/null)" == "0" ]; then
            # power file exists and power is enabled, so valves could be enabled.
            # Ensure that they're set to disabled.
            # This disables both valves (18 and 15), the pump (17), and the PSU power (4).
            echo "0" > /sys/class/gpio/gpio17/value
            echo "0" > /sys/class/gpio/gpio15/value
            echo "0" > /sys/class/gpio/gpio18/value
            echo "1" > /sys/class/gpio/gpio4/value
            # notify slack we powered off the peripherals
            curl -X POST -H 'Content-type: application/json' --data '{"text":" :warning: :warning: :warning: GAD POWER SUPPLY WAS ON WITHOUT TOOLCHAIN RUNNING -  POWER NOW OFF :warning: :warning: :warning:"}' ${SAFETYNETWEBHOOK}
            echo "no toolchain running but power on, powering off"
        fi
        # notify slack that there's no toolchain running
        curl -X POST -H 'Content-type: application/json' --data '{"text":" :warning: :warning: :warning: GADToolChain IS NOT RUNNING :warning: :warning: :warning:"}' ${SAFETYNETWEBHOOK}
	sudo wall "no toolchain running" #  that's fine for this pi
    fi # else GAD_ToolChain is running
fi # else safety net is set to OFF - should we log a warning?

