#!/bin/bash

# -1 is data, 2 is signal, 1 is BKG

isMC=-1

for year in "12016" "2016" "2017" "2018";
do

for type in "isPhiCat" "isRhoCat";
do

#for cat in "isVBFtag";
#for cat in "isZinvtag";
#for cat in "isWtag";
for cat in "isZtag";
do

TMPDATFILE="output_${cat}_${type}_${year}.txt"
rm $TMPDATFILE
echo " ------------------------------------------ " > $TMPDATFILE


if [ $cat = "isVBFtag" ]; then

    if [[ "$isMC" == '2' ]]; then        
	python3 VGammaMeson_cat.py $cat $type 1010 $year #phi
	python3 VGammaMeson_cat.py $cat $type 1020 $year #rho
    fi	
    if [[ "$isMC" == '1' ]]; then    
	python3 VGammaMeson_cat.py $cat $type 6 $year >> $TMPDATFILE 
	python3 VGammaMeson_cat.py $cat $type 7 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 8 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 9 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 20 $year >> $TMPDATFILE 
	python3 VGammaMeson_cat.py $cat $type 21 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 22 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 23 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 24 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 25 $year >> $TMPDATFILE    
    fi
    if [[ "$isMC" == '-1' ]]; then        
	python3 VGammaMeson_cat.py $cat $type -31 $year >> $TMPDATFILE  # EGammaA 2018
	python3 VGammaMeson_cat.py $cat $type -32 $year >> $TMPDATFILE  # EGammaB 2018
	python3 VGammaMeson_cat.py $cat $type -33 $year >> $TMPDATFILE  # EGammaC 2018
	python3 VGammaMeson_cat.py $cat $type -34 $year >> $TMPDATFILE  # EGammaD 2018
    fi
fi

if [ $cat = "isZinvtag" ]; then

    if [[ "$isMC" == '2' ]]; then
	python3 VGammaMeson_cat.py $cat $type 1015  $year >> $TMPDATFILE   # H phi
	python3 VGammaMeson_cat.py $cat $type 1025  $year >> $TMPDATFILE   # H rho	
    fi
  
    if [[ "$isMC" == '1' ]]; then
        python3 VGammaMeson_cat.py $cat $type 2  $year >> $TMPDATFILE   # WG
	python3 VGammaMeson_cat.py $cat $type 31 $year >> $TMPDATFILE   # Wjets 0J
	python3 VGammaMeson_cat.py $cat $type 32 $year >> $TMPDATFILE   # Wjets 1J
	python3 VGammaMeson_cat.py $cat $type 33 $year >> $TMPDATFILE   # Wjets 2J	
    fi
fi


if [ $cat = "isZtag" ]; then

    if [[ "$isMC" == '2' ]]; then
        python3 VGammaMeson_cat.py $cat $type 1013 $year >> $TMPDATFILE  # Z-phi
        python3 VGammaMeson_cat.py $cat $type 1023 $year >> $TMPDATFILE   # Z-omega
    fi

    if [[ "$isMC" == '1' ]]; then
	python3 VGammaMeson_cat.py $cat $type 1  $year >> $TMPDATFILE   # ZG
	#python3 VGammaMeson_cat.py $cat $type 0  $year >> $TMPDATFILE   # DYjets
	python3 VGammaMeson_cat.py $cat $type 34 $year >> $TMPDATFILE   # DYjets 0J
	python3 VGammaMeson_cat.py $cat $type 35 $year >> $TMPDATFILE   # DYjets 1J
	python3 VGammaMeson_cat.py $cat $type 36 $year >> $TMPDATFILE   # DYjets 2J

	##python3 VGammaMeson_cat.py $cat $type 31 $year >> $TMPDATFILE   # Wjets 0J
	##python3 VGammaMeson_cat.py $cat $type 32 $year >> $TMPDATFILE   # Wjets 1J
	##python3 VGammaMeson_cat.py $cat $type 33 $year >> $TMPDATFILE   # Wjets 2J
	##python3 VGammaMeson_cat.py $cat $type 2  $year >> $TMPDATFILE   # WG
	##python3 VGammaMeson_cat.py $cat $type 3  $year >> $TMPDATFILE   # Wjets
	##python3 VGammaMeson_cat.py $cat $type 4  $year >> $TMPDATFILE   # TT2l
	##python3 VGammaMeson_cat.py $cat $type 5  $year >> $TMPDATFILE   # TT1l
    fi

    if [[ "$isMC" == '-1' ]]; then
	python3 VGammaMeson_cat.py $cat $type -1 $year >> $TMPDATFILE  # SingleMuonA 2018
	python3 VGammaMeson_cat.py $cat $type -2 $year >> $TMPDATFILE  # SingleMuonB 2018
	python3 VGammaMeson_cat.py $cat $type -3 $year >> $TMPDATFILE  # SingleMuonC 2018
	python3 VGammaMeson_cat.py $cat $type -4 $year >> $TMPDATFILE  # SingleMuonD 2018

	python3 VGammaMeson_cat.py $cat $type -11 $year >> $TMPDATFILE  # DoubleMuonA 2018
	python3 VGammaMeson_cat.py $cat $type -12 $year >> $TMPDATFILE  # DoubleMuonB 2018
	python3 VGammaMeson_cat.py $cat $type -13 $year >> $TMPDATFILE  # DoubleMuonC 2018
	python3 VGammaMeson_cat.py $cat $type -14 $year >> $TMPDATFILE  # DoubleMuonD 2018

	python3 VGammaMeson_cat.py $cat $type -21 $year >> $TMPDATFILE  # MuonEGA 2018
	python3 VGammaMeson_cat.py $cat $type -22 $year >> $TMPDATFILE  # MuonEGB 2018
	python3 VGammaMeson_cat.py $cat $type -23 $year >> $TMPDATFILE  # MuonEGC 2018
	python3 VGammaMeson_cat.py $cat $type -24 $year >> $TMPDATFILE  # MuonEGD 2018

	python3 VGammaMeson_cat.py $cat $type -31 $year >> $TMPDATFILE  # EGammaA 2018
	python3 VGammaMeson_cat.py $cat $type -32 $year >> $TMPDATFILE  # EGammaB 2018
	python3 VGammaMeson_cat.py $cat $type -33 $year >> $TMPDATFILE  # EGammaC 2018
	python3 VGammaMeson_cat.py $cat $type -34 $year >> $TMPDATFILE  # EGammaD 2018
	fi	
    fi

if [ $cat = "isWtag" ]; then

    if [[ "$isMC" == '-1' ]]; then
	python3 VGammaMeson_cat.py $cat $type -1 $year >> $TMPDATFILE  # SingleMuonA 2018
	python3 VGammaMeson_cat.py $cat $type -2 $year >> $TMPDATFILE  # SingleMuonB 2018
	python3 VGammaMeson_cat.py $cat $type -3 $year >> $TMPDATFILE  # SingleMuonC 2018
	python3 VGammaMeson_cat.py $cat $type -4 $year >> $TMPDATFILE  # SingleMuonD 2018
#	python3 VGammaMeson_cat.py $cat $type -5 $year >> $TMPDATFILE  # SingleMuonD 2018
#	python3 VGammaMeson_cat.py $cat $type -6 $year >> $TMPDATFILE  # SingleMuonD 2018
#	python3 VGammaMeson_cat.py $cat $type -7 $year >> $TMPDATFILE  # SingleMuonD 2018
#	python3 VGammaMeson_cat.py $cat $type -8 $year >> $TMPDATFILE  # SingleMuonD 2018

	python3 VGammaMeson_cat.py $cat $type -11 $year >> $TMPDATFILE  # DoubleMuonA 2018
	python3 VGammaMeson_cat.py $cat $type -12 $year >> $TMPDATFILE  # DoubleMuonB 2018
	python3 VGammaMeson_cat.py $cat $type -13 $year >> $TMPDATFILE  # DoubleMuonC 2018
     	python3 VGammaMeson_cat.py $cat $type -14 $year >> $TMPDATFILE  # DoubleMuonD 2018
#	python3 VGammaMeson_cat.py $cat $type -15 $year >> $TMPDATFILE  # DoubleMuonD 2018
#	python3 VGammaMeson_cat.py $cat $type -16 $year >> $TMPDATFILE  # DoubleMuonD 2018
#	python3 VGammaMeson_cat.py $cat $type -17 $year >> $TMPDATFILE  # DoubleMuonD 2018
#	python3 VGammaMeson_cat.py $cat $type -18 $year >> $TMPDATFILE  # DoubleMuonD 2018

	python3 VGammaMeson_cat.py $cat $type -21 $year >> $TMPDATFILE  # MuonEGA 2018
	python3 VGammaMeson_cat.py $cat $type -22 $year >> $TMPDATFILE  # MuonEGB 2018
	python3 VGammaMeson_cat.py $cat $type -23 $year >> $TMPDATFILE  # MuonEGC 2018
	python3 VGammaMeson_cat.py $cat $type -24 $year >> $TMPDATFILE  # MuonEGD 2018
#	python3 VGammaMeson_cat.py $cat $type -25 $year >> $TMPDATFILE  # MuonEGD 2018
#	python3 VGammaMeson_cat.py $cat $type -26 $year >> $TMPDATFILE  # MuonEGD 2018
#	python3 VGammaMeson_cat.py $cat $type -27 $year >> $TMPDATFILE  # MuonEGD 2018
#	python3 VGammaMeson_cat.py $cat $type -28 $year >> $TMPDATFILE  # MuonEGD 2018	

	python3 VGammaMeson_cat.py $cat $type -31 $year >> $TMPDATFILE  # EGammaA 2018
	python3 VGammaMeson_cat.py $cat $type -32 $year >> $TMPDATFILE  # EGammaB 2018
	python3 VGammaMeson_cat.py $cat $type -33 $year >> $TMPDATFILE  # EGammaC 2018
	python3 VGammaMeson_cat.py $cat $type -34 $year >> $TMPDATFILE  # EGammaD 2018

#	python3 VGammaMeson_cat.py $cat $type -41 $year >> $TMPDATFILE  # DoubleEG 2017
#	python3 VGammaMeson_cat.py $cat $type -42 $year >> $TMPDATFILE  # DoubleEG 2017
#	python3 VGammaMeson_cat.py $cat $type -43 $year >> $TMPDATFILE  # DoubleEG 2017
#	python3 VGammaMeson_cat.py $cat $type -44 $year >> $TMPDATFILE  # DoubleEG 2017
#	python3 VGammaMeson_cat.py $cat $type -45 $year >> $TMPDATFILE  # DoubleEG 2017
#	python3 VGammaMeson_cat.py $cat $type -46 $year >> $TMPDATFILE  # DoubleEG 2017
#	python3 VGammaMeson_cat.py $cat $type -47 $year >> $TMPDATFILE  # DoubleEG 2017
#	python3 VGammaMeson_cat.py $cat $type -48 $year >> $TMPDATFILE  # DoubleEG 2017

#	python3 VGammaMeson_cat.py $cat $type -51 $year >> $TMPDATFILE  # SingleEle 2017-16
#	python3 VGammaMeson_cat.py $cat $type -52 $year >> $TMPDATFILE  # SingleEle 2017-16
#	python3 VGammaMeson_cat.py $cat $type -53 $year >> $TMPDATFILE  # SingleEle 2017-16
#	python3 VGammaMeson_cat.py $cat $type -54 $year >> $TMPDATFILE  # SingleEle 2017-16
#	python3 VGammaMeson_cat.py $cat $type -55 $year >> $TMPDATFILE  # SingleEle 2017-16
#	python3 VGammaMeson_cat.py $cat $type -56 $year >> $TMPDATFILE  # SingleEle 2017-16
#	python3 VGammaMeson_cat.py $cat $type -57 $year >> $TMPDATFILE  # SingleEle 2017-16
#	python3 VGammaMeson_cat.py $cat $type -58 $year >> $TMPDATFILE  # SingleEle 2017-16
    fi

    if [[ "$isMC" == '2' ]]; then
	python3 VGammaMeson_cat.py $cat $type 1011 $year >> $TMPDATFILE #w+ Phi
	python3 VGammaMeson_cat.py $cat $type 1012 $year >> $TMPDATFILE #w-
	python3 VGammaMeson_cat.py $cat $type 1013 $year >> $TMPDATFILE #zh
	python3 VGammaMeson_cat.py $cat $type 1021 $year >> $TMPDATFILE #w+ RHo
	python3 VGammaMeson_cat.py $cat $type 1022 $year >> $TMPDATFILE #w-
	python3 VGammaMeson_cat.py $cat $type 1023 $year >> $TMPDATFILE #zh
    fi

    if [[ "$isMC" == '1' ]]; then
	python3 VGammaMeson_cat.py $cat $type 2  $year >> $TMPDATFILE  # WG
##	python3 VGammaMeson_cat.py $cat $type 3  $year >> $TMPDATFILE  # Wjets
	python3 VGammaMeson_cat.py $cat $type 31 $year >> $TMPDATFILE  # Wjets 0J
	python3 VGammaMeson_cat.py $cat $type 32 $year >> $TMPDATFILE  # Wjets 1J
	python3 VGammaMeson_cat.py $cat $type 33 $year >> $TMPDATFILE  # Wjets 2J
	python3 VGammaMeson_cat.py $cat $type 1  $year >> $TMPDATFILE  # ZG
##	python3 VGammaMeson_cat.py $cat $type 0  $year >> $TMPDATFILE  # DYjets
	python3 VGammaMeson_cat.py $cat $type 34 $year >> $TMPDATFILE  # DYjets 0J
	python3 VGammaMeson_cat.py $cat $type 35 $year >> $TMPDATFILE  # DYjets 1J
	python3 VGammaMeson_cat.py $cat $type 36 $year >> $TMPDATFILE  # DYjets 2J
	python3 VGammaMeson_cat.py $cat $type 4  $year >> $TMPDATFILE  # TT2l
#	python3 VGammaMeson_cat.py $cat $type 5  $year >> $TMPDATFILE  # TT1l  (BROKEN SAMPLE - 2018)
    fi

fi

done
done
done
