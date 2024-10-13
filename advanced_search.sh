#!/bin/bash

# Set the value of the year variable
runNumber=$1
year=$2
target=$3

# Create the base HTCondor submission file content
cat <<EOL > job.sub
executable  = play_args.sh
arguments   = \$(file) $runNumber $year $target
output      = output/onur.\$(ClusterId).\$(ProcId).out
error       = error/onur.\$(ClusterId).\$(ProcId).err
log         = log/onur.log
+JobFlavour = "longlunch"
Transfer_output_files = ""
max_retries = 3
EOL

# Add conditional queue lines based on the value of year
if [[ $year -eq 2022 ]]; then
    echo "queue file matching /eos/experiment/sndlhc/convertedData/physics/$year/run_00$runNumber/sndsw*.root" >> job.sub
elif [[ $year -eq 2023 ]]; then
    echo "queue file matching /eos/experiment/sndlhc/convertedData/physics/${year}_reprocess/run_00$runNumber/sndsw*.root" >> job.sub
elif [[ $year -eq 2024 ]]; then
    echo "queue file matching /eos/experiment/sndlhc/convertedData/physics/$year/run_$target/run_00$runNumber/sndsw*.root" >> job.sub
else
    echo "Unsupported year: $year. No jobs will be queued."
    exit 1  # Exit the script with a non-zero status
fi

# Submit the job to HTCondor
condor_submit job.sub
#cat job.sub
