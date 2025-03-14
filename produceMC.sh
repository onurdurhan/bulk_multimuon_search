#!/bin/bash
source /cvmfs/sndlhc.cern.ch/SNDLHC-2024/June25/setUp.sh
source /afs/cern.ch/work/o/onur/SNDLHCSOFT_2024/config.sh
cd $(pwd)
#xrdcp root://eosuser.cern.ch//eos/user/o/onur/FLUKA/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr7GeVloss.root .
#xrdcp root://eosuser.cern.ch//eos/user/o/onur/FLUKA/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_30GeVcut.root .
xrdcp root://eosuser.cern.ch//eos/user/o/onur/FLUKA/LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr6.4GeVloss_with_reco_tracks.root .
FILE=/eos/experiment/sndlhc/users/odurhan/MuonicTridentMC_6.4GeVloss_MolasseRock_with_reco_tracks/sndLHC.Ntuple-TGeant4_boost100.0_LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_6.4GeVloss_with_reco_tracks-digCPP-$1.root
if ! xrdfs root://eosuser.cern.ch stat "$FILE" > /dev/null 2>&1; then
    echo "File does not exits. Running the simulation for $FILE"
    python $SNDSW_ROOT/shipLHC/run_simSND.py --boostFactor 100 --Ntuple -f LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr6.4GeVloss_with_reco_tracks.root -y 2022
    python $SNDSW_ROOT/shipLHC/run_digiSND.py -f sndLHC.Ntuple-TGeant4_boost100.0.root -g geofile_full.Ntuple-TGeant4_boost100.0.root -cpp
    xrdcp sndLHC.Ntuple-TGeant4_boost100.0_digCPP.root $FILE
else
    echo "echo $FILE already exists. Nothing to do."
fi

GEO_FILE=/eos/experiment/sndlhc/users/odurhan/MuonicTridentMC_6.4GeVloss_MolasseRock_with_reco_tracks/geofile_full.Ntuple-TGeant4_boost100.0.root
if ! xrdfs root://eosuser.cern.ch stat "$GEO_FILE"> /dev/null 2>&1; then
    echo "geo file does not exist, copying"
    xrdcp geofile_full.Ntuple-TGeant4_boost100.0.root $EOSSHIP/$GEO_FILE 
else
    echo "geo file is present. nothing to do"
fi

#FILTERED_FILE=/eos/experiment/sndlhc/users/odurhan/flushed_passing_mu3_withScifi_MolasseRock/sndLHC.Ntuple-TGeant4_boost100.0_LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_7GeVloss-digCPP-$1.root
#python /afs/cern.ch/user/o/onur/bulk_multimuon_search/MultiMuonMC.py -f $EOSSHIP/$FILE -g $EOSSHIP/$GEO_FILE -j $1 -t scifi
xrdcp histos*.root root://eosuser.cern.ch//eos/user/o/onur/multi_muon_histos_Molasse
#if ! xrdfs $EOSSHIP stat "$FILTERED_FILE" > /dev/null 2>&1; then
#    echo "filtered file does not exist, filtering"
#    python /afs/cern.ch/user/o/onur/bulk_multimuon_search/filter_MC.py -f $EOSSHIP/$FILE -g $EOSSHIP/$GEO_FILE -j $1 -t scifi
#    xrdcp filtered_clusters.root $EOSSHIP/$FILTERED_FILE
#    xrdcp data_weighted_*.json /eos/user/o/onur/dicts/scifi/weighted
#    xrdcp data_unweighted_*.json /eos/user/o/onur/dicts/scifi/unweighted
#else
#    echo "$FILTERED_FILE already exists. Nothing to do"
#fi
#xrdcp data_weighted_*.json /eos/user/o/onur/dicts_Molasse/scifi/weighted
#xrdcp data_unweighted_*.json /eos/user/o/onur/dicts_Molasse/scifi/unweighted
