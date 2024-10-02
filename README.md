# Hrare

*==================  below is for Run3 =================


scram p CMSSW_13_0_10

cd CMSSW_13_0_10/src/

cmsenv

git clone git@github.com:mariadalfonso/Hrare.git --branch D05Run3

# nanoV12 GT and campaigns for data and MC

https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/Releases/NanoAODv12

| Input dataset | GT | Era |
| ------ | ------ | ------ |
| Run3Summer22 | 130X_mcRun3_2022_realistic_v5 | Run3 |
| Run3Summer22EE | 130X_mcRun3_2022_realistic_postEE_v6 | Run3 |
| Run3Summer23 | 130X_mcRun3_2023_realistic_v14 | Run3_2023 |
| Run3Summer23BPix | 130X_mcRun3_2023_realistic_postBPix_v2 | Run3_2023 |

| 2022-CDE | 130X_dataRun3_Prompt_v3 | Run3 |
| 2022-FG |  130X_dataRun3_PromptAnalysis_v1  | Run3 |
| 2023-BCD | 130X_dataRun3_PromptAnalysis_v1 | Run3 |

* PPD pointers
* Run2022  https://docs.google.com/presentation/d/1F4ndU7DBcyvrEEyLfYqb29NGkBPs20EAnBxe_l7AEII/edit
* Run2023  https://docs.google.com/presentation/d/1TjPem5jX0fzqvTGl271_nQFoVBabsrdrO0i8Qo1uD5E/edit#slide=id.g289f499aa6b_2_58
* see status of production here https://pdmv-pages.web.cern.ch/run_3_data/full_table.html

# cmsDriver Options
--customise=Hrare/NanoAOD/nano_cff.nanoAOD_customizeMesonsRun3


*==================  below is for Run2 =================


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
