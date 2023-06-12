#!/bin/bash

#cardDIR=$1
#wsDIR=$2

wsDIR="WS_MARCH20"
cardDIR="workspaces"
resultDir="workspaces"
resultFile="results_MARCH20.txt"
rm -rf DATACARDS_MARCH20

#wsDIR="WSmva_MARCH20"
#cardDIR="workspaces"
#resultDir="workspaces"
#resultFile="results_MVA_MARCH20.txt"
#rm -rf DATACARDSmva_MARCH20

##########

echo $cardDIR
echo $wsDIR

for meson in "Phi" "Rho";

do
    echo $meson
    
    python bwsHrare.py --whichMeson=_$meson\Cat --whichCat=GFcat --inputFileSig=$wsDIR/Signal_GFcat__$meson\Cat_2018_workspace.root --inputFileBKG=$wsDIR/Bkg_GFcat__$meson\Cat_2018_workspace.root --output=$cardDIR/workspace_$meson\_GFcat_2018.root --datCardName=$cardDIR/datacard_$meson\_GFcat_2018.txt

    python bwsHrare.py --whichMeson=_$meson\Cat --whichCat=VBFcatlow --inputFileSig=$wsDIR/Signal_VBFcatlow__$meson\Cat_2018_workspace.root --inputFileBKG=$wsDIR/Bkg_VBFcatlow__$meson\Cat_2018_workspace.root --output=$cardDIR/workspace_$meson\_VBFcatlow_2018.root --datCardName=$cardDIR/datacard_$meson\_VBFcatlow_2018.txt

    python bwsHrare.py --whichMeson=_$meson\Cat --whichCat=VBFcat --inputFileSig=$wsDIR/Signal_VBFcat__$meson\Cat_Run2_workspace.root --inputFileBKG=$wsDIR/Bkg_VBFcat__$meson\Cat_Run2_workspace.root --output=$cardDIR/workspace_$meson\_VBFcat_Run2.root --datCardName=$cardDIR/datacard_$meson\_VBFcat_Run2.txt

#    python bwsHrare.py --whichMeson=_$meson\Cat --whichCat=Zinvcat --inputFileSig=$wsDIR/Signal_Zinvcat__$meson\Cat_2018_workspace.root --inputFileBKG=$wsDIR/Bkg_Zinvcat__$meson\Cat_2018_workspace.root --output=$cardDIR/workspace_$meson\_Zinvcat_2018.root --datCardName=$cardDIR/datacard_$meson\_Zinvcat_2018.txt

#    python bwsHrare.py --whichMeson=_$meson\Cat --whichCat=Wcat --inputFileSig=$wsDIR/Signal_Wcat__$meson\Cat_Run2_workspace.root --inputFileBKG=$wsDIR/Bkg_Wcat__$meson\Cat_Run2_workspace.root --output=$cardDIR/workspace_$meson\Cat_Wcat_Run2.root --datCardName=$cardDIR/datacard_$meson\_Wcat_Run2.txt

#    python bwsHrare.py --whichMeson=_$meson\Cat --whichCat=Zcat --inputFileSig=$wsDIR/Signal_Zcat__$meson\Cat_Run2_workspace.root --inputFileBKG=$wsDIR/Bkg_Zcat__$meson\Cat_Run2_workspace.root --output=$cardDIR/workspace_$meson\Cat_Zcat_Run2.root --datCardName=$cardDIR/datacard_$meson\_Zcat_Run2.txt

    python bwsHrare.py --whichMeson=_$meson\Cat --whichCat=Vcat --inputFileSig=$wsDIR/Signal_Vcat__$meson\Cat_Run2_workspace.root --inputFileBKG=$wsDIR/Bkg_Vcat__$meson\Cat_Run2_workspace.root --output=$cardDIR/workspace_$meson\Cat_Vcat_Run2.root --datCardName=$cardDIR/datacard_$meson\_Vcat_Run2.txt

##########
done

echo $resultFile

echo ' **** GFcat ****' > $resultFile

combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Phi_GFcat_2018.txt -n PhiGFcat --run expected >> $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Rho_GFcat_2018.txt -n RhoGFcat --run expected >> $resultFile

echo ' **** VBFcatlow ****' >> $resultFile

combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Phi_VBFcatlow_2018.txt -n PhiVBFcatlow --run expected >> $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Rho_VBFcatlow_2018.txt -n RhoVBFcatlow --run expected >> $resultFile

echo ' **** VBFcat ****' >> $resultFile

combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Phi_VBFcat_Run2.txt -n PhiVBFcat --run expected >> $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Rho_VBFcat_Run2.txt -n RhoVBFcat --run expected >> $resultFile

#echo ' **** Zinvcat ****' >> $resultFile

#combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Phi_Zinvcat_2018.txt -n PhiZinvcat --run expected >> $resultFile
#combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Rho_Zinvcat_2018.txt -n RhoZinvcat --run expected >> $resultFile

#combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Phi_Wcat_Run2.txt -n PhiWcat --run expected >> $resultFile
#combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Rho_Wcat_Run2.txt -n RhoWcat --run expected >> $resultFile

#combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Phi_Zcat_Run2.txt -n PhiZcat --run expected >> $resultFile
#combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Rho_Zcat_Run2.txt -n RhoZcat --run expected  >> $resultFile

echo ' **** Vcat ****' >> $resultFile

combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Phi_Vcat_Run2.txt -n PhiVcat --run expected >> $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Rho_Vcat_Run2.txt -n RhoVcat --run expected >> $resultFile

mv higgsCombine*.AsymptoticLimits.mH125.root $resultDir

#cp -r workspaces DATACARDSmva_MARCH20
cp -r workspaces DATACARDS_MARCH20

##

exit


cd $cardDIR/

combineCards.py Wcat=datacard_Rho_Wcat_Run2.txt Zcat=datacard_Rho_Zcat_Run2.txt Zinvcat=datacard_Rho_Zinvcat_2018.txt VBFcat=datacard_Rho_VBFcat_Run2.txt VBFcatlow=datacard_Rho_VBFcatlow_2018.txt > $cardDIR/datacard_Rho_comb.txt
combineCards.py Wcat=datacard_Phi_Wcat_Run2.txt Zcat=datacard_Phi_Zcat_Run2.txt Zinvcat=datacard_Phi_Zinvcat_2018.txt VBFcat=datacard_Phi_VBFcat_Run2.txt VBFcatlow=datacard_Phi_VBFcatlow_2018.txt > $cardDIR/datacard_Phi_comb.txt

cd ..

combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Phi_comb.txt -n PhiCombcat --run expected >> $resultFile
combine -M AsymptoticLimits -m 125 -t -1 $cardDIR/datacard_Rho_comb.txt -n RhoCombcat --run expected >> $resultFile

exit

#combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Phi_Wcat_2018.txt >> $resultFile
#combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Rho_Wcat_2018.txt >> $resultFile

#combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Phi_Zcat_2018.txt >> $resultFile
#combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Rho_Zcat_2018.txt >> $resultFile

#combineCards.py Zinvcat=$cardDIR/datacard_STAT_Rho_Zinvcat_2018.txt Wcat=$cardDIR/datacard_STAT_Rho_Wcat_2018.txt Zcat=$cardDIR/datacard_STAT_Rho_Zcat_2018.txt VBFcat=$cardDIR/datacard_STAT_Rho_VBFcat_2018.txt VBFcatlow=$cardDIR/datacard_STAT_Rho_VBFcatlow_2018.txt > $cardDIR/datacard_STAT_Rho_comb.txt
#combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Rho_comb.txt

#combineCards.py Zinvcat=$cardDIR/datacard_STAT_Phi_Zinvcat_2018.txt Wcat=$cardDIR/datacard_STAT_Phi_Wcat_2018.txt Zcat=$cardDIR/datacard_STAT_Phi_Zcat_2018.txt VBFcat=$cardDIR/datacard_STAT_Phi_VBFcat_2018.txt VBFcatlow=$cardDIR/datacard_STAT_Phi_VBFcatlow_2018.txt > $cardDIR/datacard_STAT_Phi_comb.txt
#combine -M AsymptoticLimits -t -1 $cardDIR/datacard_STAT_Phi_comb.txt

