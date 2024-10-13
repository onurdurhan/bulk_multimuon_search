#!/bin/bash
echo "setting up the environment"
source /cvmfs/sndlhc.cern.ch/SNDLHC-2024/June25/setUp.sh
source /afs/cern.ch/work/o/onur/SNDLHCSOFT_2024/config.sh
partition=$(basename $1 | cut -d'-' -f2 | cut -d'.' -f1)
python /afs/cern.ch/user/o/onur/bulk_multimuon_search/run_mu3_engine.py filter -f "$EOSSHIP/$1" -r $2 -y $3 -t $4 -o filtered_time_$partition.root
python /afs/cern.ch/user/o/onur/bulk_multimuon_search/run_mu3_engine.py search -f filtered_time_$partition.root -r $2 -y $3 -t $4 -o mu3_search_run_00$2_$partition.root
DIR="/eos/experiment/sndlhc/users/odurhan/multi_muon_search/$3/$2"
# Check if the run directory exists
if ! xrdfs $EOSSHIP stat "$DIR" > /dev/null 2>&1; then
  # If DIR doesn't exist, create it using xrdfs
  echo "Directory $DIR does not exist. Creating it now using XRootD."
  xrdfs $EOSSHIP mkdir "$DIR"
  echo "Directory $DIR created."
else
  echo "Directory $DIR already exists."
fi
xrdcp mu3_search_run_00$2_$partition.root "$EOSSHIP/$DIR"
