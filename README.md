# Hrare

scram p CMSSW_10_6_27

cd CMSSW_10_6_27/src/

cmsenv

git clone git@github.com:mariadalfonso/Hrare.git --branch main

# nanoV9 GT and campaigns for data and MC

https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/Releases/NanoAODv9

mc: 2016 HIPM : 106X_mcRun2_asymptotic_preVFP_v11
mc: 2016 : 106X_mcRun2_asymptotic_v17
mc: 2017 : 106X_mc2017_realistic_v9
mc: 2018 : 106X_upgrade2018_realistic_v16_L1v1
data 2016/17/18:  106X_dataRun2_v35

# cmsDriver Options
--customise=Hrare/NanoAOD/nano_cff.nanoAOD_customizeMesons
