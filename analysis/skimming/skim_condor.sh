#!/bin/sh

USERPROXY=`id -u`
echo ${USERPROXY}

alias cmsvoms='voms-proxy-init -rfc -voms cms --valid 168:00 -pwstdin < $HOME/.grid-cert-passphrase; cp /tmp/x509up_u'$(id -u)' ~/'

#cmsvoms
echo ${cmsvoms}
ls -ltra ~/x509up_u${USERPROXY} 

line=1

filename=("toResubmit_"$1"_"$2"_"$3"_"$4".txt")
touch $filename

while [ $line -le 1 ]
do

#set -- $line
whichYear=$1
whichEra=$2
whichPD=$3
whichSkim=$4
whichJob=$5

echo 'PROCESSING' $whichJob

fOutDir=("/scratch/submit/cms/mariadlf/Hrare/SKIMS/D01/"$whichSkim"/"$whichYear"/"$whichPD"+Run"$whichEra)
echo $fOutDir

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
Log    = logs/simple_skim_${whichSample}_${whichJob}.log
Output = logs/simple_skim_${whichSample}_${whichJob}.out
Error  = logs/simple_skim_${whichSample}_${whichJob}.error
transfer_input_files = skim.py, ../utilsHrare.py, ../functions.cc, ../config/selection.json, ../config/skimDB.json
transfer_output_remaps = "out_${whichJob}.root = ${fOutDir}/out_${whichJob}.root"
use_x509userproxy = True
x509userproxy = /home/submit/mariadlf/x509up_u${USERPROXY}
+AccountingGroup = "analysis.mariadlf"
+DESIRED_Sites        = "T2_AT_Vienna,T2_BE_IIHE,T2_BE_UCL,T2_BR_SPRACE,T2_BR_UERJ,T2_CH_CERN,T2_CH_CERN_AI,T2_CH_CERN_HLT,T2_CH_CERN_Wigner,T2_CH_CSCS,T2_CH_CSCS_HPC,T2_CN_Beijing,T2_DE_DESY,T2_DE_RWTH,T2_EE_Estonia,T2_ES_CIEMAT,T2_ES_IFCA,T2_FI_HIP,T2_FR_CCIN2P3,T2_FR_GRIF_IRFU,T2_FR_GRIF_LLR,T2_FR_IPHC,T2_GR_Ioannina,T2_HU_Budapest,T2_IN_TIFR,T2_IT_Bari,T2_IT_Legnaro,T2_IT_Pisa,T2_IT_Rome,T2_KR_KISTI,T2_MY_SIFIR,T2_MY_UPM_BIRUNI,T2_PK_NCP,T2_PL_Swierk,T2_PL_Warsaw,T2_PT_NCG_Lisbon,T2_RU_IHEP,T2_RU_INR,T2_RU_ITEP,T2_RU_JINR,T2_RU_PNPI,T2_RU_SINP,T2_TH_CUNSTDA,T2_TR_METU,T2_TW_NCHC,T2_UA_KIPT,T2_UK_London_IC,T2_UK_SGrid_Bristol,T2_UK_SGrid_RALPP,T2_US_Caltech,T2_US_Florida,T2_US_MIT,T2_US_Nebraska,T2_US_Purdue,T2_US_UCSD,T2_US_Vanderbilt,T2_US_Wisconsin,T3_CH_CERN_CAF,T3_CH_CERN_DOMA,T3_CH_CERN_HelixNebula,T3_CH_CERN_HelixNebula_REHA,T3_CH_CMSAtHome,T3_CH_Volunteer,T3_US_HEPCloud,T3_US_NERSC,T3_US_OSG,T3_US_PSC,T3_US_SDSC"
#Requirements = ((BOSCOGroup == "bosco_cms" && BOSCOCluster == "ce03.cmsaf.mit.edu") || (BOSCOCluster == "t3serv008.mit.edu")) && (Machine != "t3btch070.mit.edu") && (Machine != "t3desk014.mit.edu")
+REQUIRED_OS = "rhel7"
Queue
EOF

line=$(( $line + 1 ))

condor_submit submit

done 

rm -f submit
