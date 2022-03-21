#!/bin/sh

source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.24.06/x86_64-centos7-gcc48-opt/bin/thisroot.sh

python3 skim.py --year=$1 --era=$2 --PDType=$3 --SkimType=$4 --whichJob=$5

status=$?

rm -f functions_cc*

if [ $status -eq 0 ]; then
  echo "SUCCESS"

else
  echo "FAILURE"

fi
