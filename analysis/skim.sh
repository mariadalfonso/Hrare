#!/bin/sh

echo "hostname"
hostname
whoami

echo $PWD

#source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.24.06/x86_64-centos7-gcc48-opt/bin/thisroot.sh

python3 skim5.py --year=$1 --era=$2 --PDType=$3 --SkimType=$4 --whichJob=$5
#python3 skim5.py --year=2018 --era="A" --PDType="EGamma" --SkimType="VH" --whichJob=101

status=$?

if [ $status -eq 0 ]; then
  echo "SUCCESS"

else
  echo "FAILURE"

fi
