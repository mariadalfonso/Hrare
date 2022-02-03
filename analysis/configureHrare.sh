#!/bin/bash

for year in "2018";
#for year in "12016" "2016" "2017" "2018";
do

#for type in "isPhiCat" "isRhoCat";
for type in "isPhiCat"	
do

#for cat in "isVBFtag";
#for cat in "isZtag";
for cat in "isWtag";
do

TMPDATFILE="output_${cat}_${type}_${year}.txt"
rm $TMPDATFILE
echo " ------------------------------------------ " > $TMPDATFILE


if [ $cat = "isVBFtag" ]; then

python3 VGammaMeson_cat.py $cat $type 6 $year >> $TMPDATFILE 
python3 VGammaMeson_cat.py $cat $type 7 $year >> $TMPDATFILE
python3 VGammaMeson_cat.py $cat $type 8 $year >> $TMPDATFILE
python3 VGammaMeson_cat.py $cat $type 9 $year >> $TMPDATFILE

python3 VGammaMeson_cat.py $cat $type -31 $year >> $TMPDATFILE  # EGammaA 2018
python3 VGammaMeson_cat.py $cat $type -32 $year >> $TMPDATFILE  # EGammaB 2018
python3 VGammaMeson_cat.py $cat $type -33 $year >> $TMPDATFILE  # EGammaC 2018
python3 VGammaMeson_cat.py $cat $type -34 $year >> $TMPDATFILE  # EGammaD 2018

# signal 12; QCD 20-25

fi

if [ $cat = "isZtag" ]; then

#python3 VGammaMeson_cat.py $cat $type 10 $year >> $TMPDATFILE   # signal ZLL
python3 VGammaMeson_cat.py $cat $type 1  $year >> $TMPDATFILE   # ZG
python3 VGammaMeson_cat.py $cat $type 0  $year >> $TMPDATFILE   # DYjets
python3 VGammaMeson_cat.py $cat $type 34 $year >> $TMPDATFILE   # DYjets 0J
python3 VGammaMeson_cat.py $cat $type 35 $year >> $TMPDATFILE   # DYjets 1J
python3 VGammaMeson_cat.py $cat $type 36 $year >> $TMPDATFILE   # DYjets 2J

#python3 VGammaMeson_cat.py $cat $type 31 $year >> $TMPDATFILE   # Wjets 0J
#python3 VGammaMeson_cat.py $cat $type 32 $year >> $TMPDATFILE   # Wjets 1J
#python3 VGammaMeson_cat.py $cat $type 33 $year >> $TMPDATFILE   # Wjets 2J
#python3 VGammaMeson_cat.py $cat $type 2  $year >> $TMPDATFILE   # WG
#python3 VGammaMeson_cat.py $cat $type 3  $year >> $TMPDATFILE   # Wjets
#python3 VGammaMeson_cat.py $cat $type 4  $year >> $TMPDATFILE   # TT2l
#python3 VGammaMeson_cat.py $cat $type 5  $year >> $TMPDATFILE   # TT1l

python3 VGammaMeson_cat.py $cat $type -1 $year >> $TMPDATFILE  # SingleMuonA 2018
python3 VGammaMeson_cat.py $cat $type -2 $year >> $TMPDATFILE  # SingleMuonB 2018
python3 VGammaMeson_cat.py $cat $type -3 $year >> $TMPDATFILE  # SingleMuonC 2018
python3 VGammaMeson_cat.py $cat $type -4 $year >> $TMPDATFILE  # SingleMuonD 2018    
    

fi

if [ $cat = "isWtag" ]; then

##python3 VGammaMeson_cat.py $cat $type 11  >> $TMPDATFILE  # signal WLNU
python3 VGammaMeson_cat.py $cat $type 2  $year >> $TMPDATFILE  # WG
python3 VGammaMeson_cat.py $cat $type 3  $year >> $TMPDATFILE  # Wjets
python3 VGammaMeson_cat.py $cat $type 31 $year >> $TMPDATFILE  # Wjets 0J
python3 VGammaMeson_cat.py $cat $type 32 $year >> $TMPDATFILE  # Wjets 1J
python3 VGammaMeson_cat.py $cat $type 33 $year >> $TMPDATFILE  # Wjets 2J
python3 VGammaMeson_cat.py $cat $type 1  $year >> $TMPDATFILE  # ZG
python3 VGammaMeson_cat.py $cat $type 0  $year >> $TMPDATFILE  # DYjets
python3 VGammaMeson_cat.py $cat $type 34 $year >> $TMPDATFILE  # DYjets 0J
python3 VGammaMeson_cat.py $cat $type 35 $year >> $TMPDATFILE  # DYjets 1J
python3 VGammaMeson_cat.py $cat $type 36 $year >> $TMPDATFILE  # DYjets 2J
python3 VGammaMeson_cat.py $cat $type 4  $year >> $TMPDATFILE  # TT2l
python3 VGammaMeson_cat.py $cat $type 5  $year >> $TMPDATFILE  # TT1l  (BROKEN SAMPLE - 2018)

python3 VGammaMeson_cat.py $cat $type -1 $year >> $TMPDATFILE  # SingleMuonA 2018
python3 VGammaMeson_cat.py $cat $type -2 $year >> $TMPDATFILE  # SingleMuonB 2018
python3 VGammaMeson_cat.py $cat $type -3 $year >> $TMPDATFILE  # SingleMuonC 2018
python3 VGammaMeson_cat.py $cat $type -4 $year >> $TMPDATFILE  # SingleMuonD 2018    
    
fi

done
done
done

exit
