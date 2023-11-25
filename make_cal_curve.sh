#!/bin/bash
./make_calcurve_DB_entry
echo '#!/bin/bash' > test.sh
#./make_calcurve_DB_entry 275_B fitParams_275_A_apr2023cal.root '2023-04-17 09:00:00' 'pol6' >> test.sh
./make_calcurve_DB_entry ${@} >> test.sh
chmod +x test.sh
echo "check contents of test.sh and execute if all looks good"
