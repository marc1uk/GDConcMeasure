#!/bin/bash
echo '#!/bin/bash' > test.sh
#./make_pureref_DB_entry pureDarkSubtracted_275_B_apr2023cal.root 275_B >> test.sh
./make_pureref_DB_entry ${@} >> test.sh
chmod +x test.sh
echo "check contents of test.sh and execute if all looks good"
