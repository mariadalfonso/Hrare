#!/bin/bash

#cardDIR=$1
#wsDIR=$2

cardDIR="DATACARDS"
wsDIR="WS"
resultFile="resultsOCT12.txt"

#cardDIR="DATACARDSmva"
#wsDIR="WSmva"
#resultFile="resultsOCT12_MVA.txt"

##########

echo $cardDIR
echo $wsDIR

for meson in "Phi" "Rho";
#for meson in "Rho";
do
    echo $meson
    
    python bwsHrare.py --whichMeson=_$meson\Cat --whichCat=GFcat --inputFileSig=$wsDIR/Signal_GFcat__$meson\Cat_2018_workspace.root --inputFileBKG=$wsDIR/Bkg_GFcat__$meson\Cat_2018_workspace.root --output=$cardDIR/workspace_STAT_$meson\_GFcat_2018.root --datCardName=$cardDIR/datacard_STAT_$meson\_GFcat_2018.txt

    python bwsHrare.py --whichMeson=_$meson\Cat --whichCat=VBFcatlow --inputFileSig=$wsDIR/Signal_VBFcatlow__$meson\Cat_2018_workspace.root --inputFileBKG=$wsDIR/Bkg_VBFcatlow__$meson\Cat_2018_workspace.root --output=$cardDIR/workspace_STAT_$meson\_VBFcatlow_2018.root --datCardName=$cardDIR/datacard_STAT_$meson\_VBFcatlow_2018.txt

    python bwsHrare.py --whichMeson=_$meson\Cat --whichCat=VBFcat --inputFileSig=$wsDIR/Signal_VBFcat__$meson\Cat_Run2_workspace.root --inputFileBKG=$wsDIR/Bkg_VBFcat__$meson\Cat_Run2_workspace.root --output=$cardDIR/workspace_STAT_$meson\_VBFcat_Run2.root --datCardName=$cardDIR/datacard_STAT_$meson\_VBFcat_Run2.txt

    python bwsHrare.py --whichMeson=_$meson\Cat --whichCat=Zinvcat --inputFileSig=$wsDIR/Signal_Zinvcat__$meson\Cat_2018_workspace.root --inputFileBKG=$wsDIR/Bkg_Zinvcat__$meson\Cat_2018_workspace.root --output=$cardDIR/workspace_STAT_$meson\_Zinvcat_2018.root --datCardName=$cardDIR/datacard_STAT_$meson\_Zinvcat_2018.txt

##

#python bwsHrare.py --whichCat=Wcat --inputFileSig=$wsDIR/Signal_Wcat__$meson\Cat_2018_workspace.root --inputFileBKG=$wsDIR/Bkg_Wcat__$meson\Cat_2018_workspace.root --output=$cardDIR/workspace_STAT_$meson\Cat_Wcat_2018.root --datCardName=$cardDIR/datacard_STAT_$meson\_Wcat_2018.txt

#python bwsHrare.py --whichCat=Zcat --inputFileSig=$wsDIR/Signal_Zcat__$meson\Cat_2018_workspace.root --inputFileBKG=$wsDIR/Bkg_Zcat__$meson\Cat_2018_workspace.root --output=$cardDIR/workspace_STAT_$meson\Cat_Zcat_2018.root --datCardName=$cardDIR/datacard_STAT_$meson\_Zcat_2018.txt

##########
done

echo $resultFile

combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Phi_GFcat_2018.txt  > $resultFile
combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Rho_GFcat_2018.txt  >> $resultFile

combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Phi_VBFcatlow_2018.txt  >> $resultFile
combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Rho_VBFcatlow_2018.txt  >> $resultFile

combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Phi_VBFcat_Run2.txt  >> $resultFile
combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Rho_VBFcat_Run2.txt  >> $resultFile

combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Phi_Zinvcat_2018.txt >> $resultFile
combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Rho_Zinvcat_2018.txt >> $resultFile

##

#combineCards.py Zinvcat=$cardDIR/datacard_STAT_Rho_Zinvcat_2018.txt VBFcat=$cardDIR/datacard_STAT_Rho_VBFcat_Run2.txt VBFcatlow=$cardDIR/datacard_STAT_Rho_VBFcatlow_2018.txt > $cardDIR/datacard_STAT_Rho_comb.txt
combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Rho_comb.txt

#combineCards.py Zinvcat=$cardDIR/datacard_STAT_Phi_Zinvcat_2018.txt VBFcat=$cardDIR/datacard_STAT_Phi_VBFcat_Run2.txt VBFcatlow=$cardDIR/datacard_STAT_Phi_VBFcatlow_2018.txt > $cardDIR/datacard_STAT_Phi_comb.txt
combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Phi_comb.txt

exit

#combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Phi_Wcat_2018.txt >> $resultFile
#combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Rho_Wcat_2018.txt >> $resultFile

#combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Phi_Zcat_2018.txt >> $resultFile
#combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Rho_Zcat_2018.txt >> $resultFile

#combineCards.py Zinvcat=$cardDIR/datacard_STAT_Rho_Zinvcat_2018.txt Wcat=$cardDIR/datacard_STAT_Rho_Wcat_2018.txt Zcat=$cardDIR/datacard_STAT_Rho_Zcat_2018.txt VBFcat=$cardDIR/datacard_STAT_Rho_VBFcat_2018.txt VBFcatlow=$cardDIR/datacard_STAT_Rho_VBFcatlow_2018.txt > $cardDIR/datacard_STAT_Rho_comb.txt
#combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Rho_comb.txt

#combineCards.py Zinvcat=$cardDIR/datacard_STAT_Phi_Zinvcat_2018.txt Wcat=$cardDIR/datacard_STAT_Phi_Wcat_2018.txt Zcat=$cardDIR/datacard_STAT_Phi_Zcat_2018.txt VBFcat=$cardDIR/datacard_STAT_Phi_VBFcat_2018.txt VBFcatlow=$cardDIR/datacard_STAT_Phi_VBFcatlow_2018.txt > $cardDIR/datacard_STAT_Phi_comb.txt
#combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Phi_comb.txt

