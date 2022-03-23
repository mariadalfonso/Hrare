#!/bin/sh

#source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.24.06/x86_64-centos7-gcc48-opt/bin/thisroot.sh

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc900
scramv1 project CMSSW CMSSW_12_2_0 # cmsrel is an alias not on the workers
cd CMSSW_12_2_0/src/
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
cd ../..

python3 skim.py --year=$1 --era=$2 --PDType=$3 --SkimType=$4 --whichJob=$5

status=$?

rm -f functions_cc*

if [ $status -eq 0 ]; then
  echo "SUCCESS"

else
  echo "FAILURE"
  filename=("toResubmit_"$1"_"$2"_"$3"_"$4".txt")
  echo 'source ./skim_condor.sh '$1' '$2' '$3' '$4' '$5'' >> $filename

fi
