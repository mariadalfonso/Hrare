# Hrare

scram p CMSSW_10_6_27

cd CMSSW_10_6_27/src/

cmsenv

git clone git@github.com:mariadalfonso/Hrare.git --branch D02

# nanoV9 GT and campaigns for data and MC

https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/Releases/NanoAODv9

| Input dataset | GT | Era |
| ------ | ------ | ------ |
| RunIISummer20UL16MiniAODAPV MC | 106X_mcRun2_asymptotic_preVFP_v11 | Run2_2016_HIPM,run2_nanoAOD_106Xv2 |
| RunIISummer20UL16MiniAOD MC | 106X_mcRun2_asymptotic_v17 | Run2_2016,run2_nanoAOD_106Xv2 |
| RunIISummer20UL17MiniAOD MC | 106X_mc2017_realistic_v9 | Run2_2017,run2_nanoAOD_106Xv2 |
| RunIISummer20UL18MiniAOD MC | 106X_upgrade2018_realistic_v16_L1v1 | Run2_2018,run2_nanoAOD_106Xv2 |
| HIPM_UL2016 data (Summer19/20) | 106X_dataRun2_v35 | Run2_2016_HIPM,run2_nanoAOD_106Xv2 |
| 2016 UL data (Summer19/20) | 106X_dataRun2_v35 | Run2_2016,run2_nanoAOD_106Xv2 |
| 2017 UL data (Summer19/20) | 106X_dataRun2_v35 | Run2_2017,run2_nanoAOD_106Xv2 |
| 2018 UL data (Summer19/20) | 106X_dataRun2_v35 | Run2_2018,run2_nanoAOD_106Xv2 |

# cmsDriver Options
--customise=Hrare/NanoAOD/nano_cff.nanoAOD_customizeMesons


# current branches
http://dalfonso.web.cern.ch/dalfonso/Hrare/QCD_Pt_15to30_size_report.html#phi
