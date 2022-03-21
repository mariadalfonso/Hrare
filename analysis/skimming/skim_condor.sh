#!/bin/sh

USERPROXY=`id -u`
echo ${USERPROXY}
#voms-proxy-init --voms cms --valid 168:00

alias cmsvoms='voms-proxy-init -rfc -voms cms; cp /tmp/x509up_u'$(id -u)' ~/'

cmsvoms
echo ${cmsvoms}
ls -ltra ~/x509up_u${USERPROXY} 

echo 'WHOAMI?'
whoami

line=1

while [ $line -le 5 ]
do

echo 'PROCESSING' $line

#set -- $line
whichYear=$1
whichEra=$2
whichPD=$3
whichSkim=$4
whichJob=$5

whichYear=2018
whichEra="A"
whichPD="EGamma"
whichSkim="VH" 
whichJob=4


if [ ! -d "/scratch/submit/cms/mariadlf/Hrare/SKIMS/D01/VBF/2018/" ]; then
  echo "creating output folders" /scratch/submit/cms/mariadlf/Hrare/SKIMS/D01/
  mkdir -p /scratch/submit/cms/mariadlf/Hrare/SKIMS/D01/VBF/2018/
  mkdir -p /scratch/submit/cms/mariadlf/Hrare/SKIMS/D01/VH/2018/

fi

cat << EOF > submit
Universe   = vanilla
Executable = skim.sh
Arguments  = ${whichYear} ${whichEra} ${whichPD} ${whichSkim} ${whichJob}
RequestMemory = 6000
RequestCpus = 1
RequestDisk = DiskUsage
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = ""
Log    = logs/simple_skim_${whichSample}_${whichJob}.log
Output = logs/simple_skim_${whichSample}_${whichJob}.out
Error  = logs/simple_skim_${whichSample}_${whichJob}.error
transfer_input_files = skim.py, ../utilsHrare.py, ../functions.cc, ../config/selection.json, ../config/skimDB.json
use_x509userproxy = True
x509userproxy = /home/submit/mariadlf/x509up_u${USERPROXY}
Requirements = ((BOSCOGroup == "bosco_cms" && BOSCOCluster == "ce03.cmsaf.mit.edu") || (BOSCOCluster == "t3serv008.mit.edu")) && (Machine != "t3btch070.mit.edu") && (Machine != "t3desk014.mit.edu")
+REQUIRED_OS = "rhel7"
Queue
EOF

line=$(( $line + 1 ))

condor_submit submit

done 

rm -f submit
