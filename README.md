# Hrare

scram p CMSSW_10_6_26

cd CMSSW_10_6_26/src/

cmsenv

git clone git@github.com:mariadalfonso/Hrare.git --branch main

# nanoV9 GT and campaigns for data and MC

https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/Releases/NanoAODv9

# cmsDriver Options
--customise=Hrare/NanoAOD/nano_cff.nanoAOD_customizeMesons
