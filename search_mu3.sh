#!/bin/bash
echo "setting up the environment"
source /cvmfs/sndlhc.cern.ch/SNDLHC-2024/June25/setUp.sh
source /afs/cern.ch/work/o/onur/SNDLHCSOFT_2024/config.sh
path_to_data=/eos/experiment/sndlhc/convertedData/physics/2022
partition=$(basename $1 | cut -d'-' -f2 | cut -d'.' -f1)
python /afs/cern.ch/user/o/onur/bulk_multimuon_search/run_mu3_engine.py filter -f $1 -r $2 -p $path_to_data -o filtered_time_$partition.root
python /afs/cern.ch/user/o/onur/bulk_multimuon_search/run_mu3_engine.py search -f filtered_time_$partition.root -r $2 -p $path_to_data -o mu3_search_run_00$2_$partition.root
xrdcp  mu3_search_run_00$2_$partition.root /eos/user/o/onur/multi_muon_search/$2
