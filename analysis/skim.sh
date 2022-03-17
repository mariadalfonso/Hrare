#!/bin/sh

python3 skim.py --year=2018 --era="A" --PDType="EGamma" --SkimType="VH"

status=$?

if [ $status -eq 0 ]; then
  echo "SUCCESS"

else
  echo "FAILURE"

fi
