#!/bin/bash

#SS
#python3 VGammaMeson_cat.py isGFtag isPhiCat -66 2018
#python3 VGammaMeson_cat.py isGFtag isPhiCat -65 2018

#python3 VGammaMeson_cat.py isVBFtag isPhiCat -66 2018
#python3 VGammaMeson_cat.py isVBFtag isPhiCat -65 2018

#python3 VGammaMeson_cat.py isVBFtaglow isPhiCat -66 2018
#python3 VGammaMeson_cat.py isVBFtaglow isPhiCat -65 2018

#python3 VGammaMeson_cat.py isZinvtag isPhiCat -66 2018
#python3 VGammaMeson_cat.py isZinvtag isPhiCat -65 2018

#exit

#-1 is data, 2 is signal, 1 is BKG

isMC=1

for year in "12016" "22016" "2017" "2018";

do

for type in "isRhoCat" "isPhiCat";

do

for cat in "isGFtag" "isVBFtag" "isVBFtaglow" "isZinvtag" "isZtag" "isWtag";
do

TMPDATFILE="output_${cat}_${type}_${year}.txt"
rm $TMPDATFILE
echo " ------------------------------------------ " > $TMPDATFILE

if [ $cat = "isGFtag" ]; then

    if [[ "$isMC" == '2' ]]; then
	python3 VGammaMeson_cat.py $cat $type 1017 $year >> $TMPDATFILE #phi
	python3 VGammaMeson_cat.py $cat $type 1027 $year >> $TMPDATFILE #rho
	python3 VGammaMeson_cat.py $cat $type 1010 $year >> $TMPDATFILE # VBF phi
	python3 VGammaMeson_cat.py $cat $type 1020 $year >> $TMPDATFILE # VBF rho
# 	 python3 VGammaMeson_cat.py $cat $type 1015 $year >> $TMPDATFILE # Zinv phi
#	 python3 VGammaMeson_cat.py $cat $type 1025 $year >> $TMPDATFILE # Zinv rho
#        python3 VGammaMeson_cat.py $cat $type 1016 $year >> $TMPDATFILE # ggZinv phi
#	 python3 VGammaMeson_cat.py $cat $type 1026 $year >> $TMPDATFILE # ggZinv rho
#        python3 VGammaMeson_cat.py $cat $type 1011 $year >> $TMPDATFILE #w+ Phi
#        python3 VGammaMeson_cat.py $cat $type 1012 $year >> $TMPDATFILE #w-
#        python3 VGammaMeson_cat.py $cat $type 1021 $year >> $TMPDATFILE #w+ Rho
#        python3 VGammaMeson_cat.py $cat $type 1022 $year >> $TMPDATFILE #w-
    fi

    if [[ "$isMC" == '1' ]]; then
#	python3 VGammaMeson_cat.py $cat $type 6 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 7 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 8 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 9 $year >> $TMPDATFILE

	python3 VGammaMeson_cat.py $cat $type 10 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 11 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 12 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 13 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 14 $year >> $TMPDATFILE

#	python3 VGammaMeson_cat.py $cat $type 20 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 21 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 22 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 23 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 24 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 25 $year >> $TMPDATFILE

#	python3 VGammaMeson_cat.py $cat $type 45 $year >> $TMPDATFILE  &  # WGamma had
#	python3 VGammaMeson_cat.py $cat $type 46 $year >> $TMPDATFILE  &  # ZGamma had
#	python3 VGammaMeson_cat.py $cat $type 47 $year >> $TMPDATFILE  &  # TTGamma

#	python3 VGammaMeson_cat.py $cat $type 31 $year >> $TMPDATFILE  &  # Wjets 0J
#	python3 VGammaMeson_cat.py $cat $type 32 $year >> $TMPDATFILE  &  # Wjets 1J
#	python3 VGammaMeson_cat.py $cat $type 33 $year >> $TMPDATFILE  &  # Wjets 2J

    fi

    if [[ "$isMC" == '-1' ]]; then
	python3 VGammaMeson_cat.py $cat $type -62 $year >> $TMPDATFILE  # EGammaB 2018
	python3 VGammaMeson_cat.py $cat $type -63 $year >> $TMPDATFILE  # EGammaC 2018
	python3 VGammaMeson_cat.py $cat $type -64 $year >> $TMPDATFILE  # EGammaD 2018
    fi

fi

if [ $cat = "isVBFtaglow" ]; then

    if [[ "$isMC" == '2' ]]; then
	python3 VGammaMeson_cat.py $cat $type 1010 $year >> $TMPDATFILE #phi
	python3 VGammaMeson_cat.py $cat $type 1020 $year >> $TMPDATFILE #rho
	python3 VGammaMeson_cat.py $cat $type 1017 $year >> $TMPDATFILE #phi
	python3 VGammaMeson_cat.py $cat $type 1027 $year >> $TMPDATFILE #rho
    fi

    if [[ "$isMC" == '1' ]]; then

#	python3 VGammaMeson_cat.py $cat $type 45 $year >> $TMPDATFILE #WG
#	python3 VGammaMeson_cat.py $cat $type 46 $year >> $TMPDATFILE #ZG
#	python3 VGammaMeson_cat.py $cat $type 47 $year >> $TMPDATFILE #TTG

#	python3 VGammaMeson_cat.py $cat $type 6 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 7 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 8 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 9 $year >> $TMPDATFILE

	python3 VGammaMeson_cat.py $cat $type 10 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 11 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 12 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 13 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 14 $year >> $TMPDATFILE

	python3 VGammaMeson_cat.py $cat $type 15 $year >> $TMPDATFILE #VBFGamma

#	python3 VGammaMeson_cat.py $cat $type 20 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 21 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 22 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 23 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 24 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 25 $year >> $TMPDATFILE
    fi

    if [[ "$isMC" == '-1' ]]; then
	python3 VGammaMeson_cat.py $cat $type -62 $year >> $TMPDATFILE  # EGammaB 2018
	python3 VGammaMeson_cat.py $cat $type -63 $year >> $TMPDATFILE  # EGammaC 2018
	python3 VGammaMeson_cat.py $cat $type -64 $year >> $TMPDATFILE  # EGammaD 2018
    fi

fi

if [ $cat = "isVBFtag" ]; then

    if [[ "$isMC" == '2' ]]; then
	python3 VGammaMeson_cat.py $cat $type 1010 $year >> $TMPDATFILE #phi
	python3 VGammaMeson_cat.py $cat $type 1020 $year >> $TMPDATFILE #rho
	python3 VGammaMeson_cat.py $cat $type 1017 $year >> $TMPDATFILE #phi
	python3 VGammaMeson_cat.py $cat $type 1027 $year >> $TMPDATFILE #rho
    fi

    if [[ "$isMC" == '1' ]]; then

#	python3 VGammaMeson_cat.py $cat $type 45 $year >> $TMPDATFILE #WGamma had
#	python3 VGammaMeson_cat.py $cat $type 46 $year >> $TMPDATFILE #ZGamma had
#	python3 VGammaMeson_cat.py $cat $type 47 $year >> $TMPDATFILE #TTG

#	python3 VGammaMeson_cat.py $cat $type 6 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 7 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 8 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 9 $year >> $TMPDATFILE

	python3 VGammaMeson_cat.py $cat $type 10 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 11 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 12 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 13 $year >> $TMPDATFILE
	python3 VGammaMeson_cat.py $cat $type 14 $year >> $TMPDATFILE

	python3 VGammaMeson_cat.py $cat $type 15 $year >> $TMPDATFILE #VBFGamma

#	python3 VGammaMeson_cat.py $cat $type 20 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 21 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 22 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 23 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 24 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 25 $year >> $TMPDATFILE
    fi

    if [[ "$isMC" == '-1' ]]; then

	if [ $year = "2018" ] ; then
	    python3 VGammaMeson_cat.py $cat $type -62 $year >> $TMPDATFILE  # TauA 2018
            python3 VGammaMeson_cat.py $cat $type -63 $year >> $TMPDATFILE  # TauB 2018
            python3 VGammaMeson_cat.py $cat $type -64 $year >> $TMPDATFILE  # TauB 2018

	    python3 VGammaMeson_cat.py $cat $type -31 $year >> $TMPDATFILE  # EGammaA 2018
	    python3 VGammaMeson_cat.py $cat $type -32 $year >> $TMPDATFILE  # EGammaB 2018
	    python3 VGammaMeson_cat.py $cat $type -33 $year >> $TMPDATFILE  # EGammaC 2018
	    python3 VGammaMeson_cat.py $cat $type -34 $year >> $TMPDATFILE  # EGammaD 2018
	fi

	if [ $year = "2017" ]; then
	   python3 VGammaMeson_cat.py $cat $type -76 $year >> $TMPDATFILE  # SinglePhoton RunF
	fi

	if [ $year = "12016" ]; then
	   python3 VGammaMeson_cat.py $cat $type -81 $year >> $TMPDATFILE  # SinglePhoton RunB1
	   python3 VGammaMeson_cat.py $cat $type -82 $year >> $TMPDATFILE  # SinglePhoton RunB2
	   python3 VGammaMeson_cat.py $cat $type -83 $year >> $TMPDATFILE  # SinglePhoton RunC
	   python3 VGammaMeson_cat.py $cat $type -84 $year >> $TMPDATFILE  # SinglePhoton RunD
	   python3 VGammaMeson_cat.py $cat $type -85 $year >> $TMPDATFILE  # SinglePhoton RunE
	   python3 VGammaMeson_cat.py $cat $type -86 $year >> $TMPDATFILE  # SinglePhoton RunF
	fi
    fi
fi

if [ $cat = "isZinvtag" ]; then

    if [[ "$isMC" == '1' ]]; then

#       python3 VGammaMeson_cat.py $cat $type 6 $year >> $TMPDATFILE
#       python3 VGammaMeson_cat.py $cat $type 7 $year >> $TMPDATFILE
#       python3 VGammaMeson_cat.py $cat $type 8 $year >> $TMPDATFILE
#       python3 VGammaMeson_cat.py $cat $type 9 $year >> $TMPDATFILE

#       python3 VGammaMeson_cat.py $cat $type 10 $year >> $TMPDATFILE
#       python3 VGammaMeson_cat.py $cat $type 11 $year >> $TMPDATFILE
#       python3 VGammaMeson_cat.py $cat $type 12 $year >> $TMPDATFILE
#       python3 VGammaMeson_cat.py $cat $type 13 $year >> $TMPDATFILE
#       python3 VGammaMeson_cat.py $cat $type 14 $year >> $TMPDATFILE

#	python3 VGammaMeson_cat.py $cat $type 37  $year >> $TMPDATFILE   # Zinv
#	python3 VGammaMeson_cat.py $cat $type 38  $year >> $TMPDATFILE   # Zinv
#	python3 VGammaMeson_cat.py $cat $type 39  $year >> $TMPDATFILE   # Zinv
#	python3 VGammaMeson_cat.py $cat $type 40  $year >> $TMPDATFILE   # Zinv
#	python3 VGammaMeson_cat.py $cat $type 41  $year >> $TMPDATFILE   # Zinv
#	python3 VGammaMeson_cat.py $cat $type 42  $year >> $TMPDATFILE   # Zinv
#	python3 VGammaMeson_cat.py $cat $type 43  $year >> $TMPDATFILE   # Zinv
#	python3 VGammaMeson_cat.py $cat $type 44  $year >> $TMPDATFILE   # Zinv

       python3 VGammaMeson_cat.py $cat $type 2  $year >> $TMPDATFILE   # WG to lnu
       python3 VGammaMeson_cat.py $cat $type 45 $year >> $TMPDATFILE   # WG had

       python3 VGammaMeson_cat.py $cat $type 31 $year >> $TMPDATFILE   # Wjets 0J
       python3 VGammaMeson_cat.py $cat $type 32 $year >> $TMPDATFILE   # Wjets 1J
       python3 VGammaMeson_cat.py $cat $type 33 $year >> $TMPDATFILE   # Wjets 2J

       python3 VGammaMeson_cat.py $cat $type 46 $year >> $TMPDATFILE   # ZG had
       python3 VGammaMeson_cat.py $cat $type 1  $year >> $TMPDATFILE   # ZG to ll
       python3 VGammaMeson_cat.py $cat $type 48 $year >> $TMPDATFILE   # ZG to nn

       python3 VGammaMeson_cat.py $cat $type 47  $year >> $TMPDATFILE   # TTG
#       python3 VGammaMeson_cat.py $cat $type 15  $year >> $TMPDATFILE   # G+jets VBF
#	python3 VGammaMeson_cat.py $cat $type 4  $year >> $TMPDATFILE   # TTG 1l
#	python3 VGammaMeson_cat.py $cat $type 5  $year >> $TMPDATFILE   # TTG 2l

#	python3 VGammaMeson_cat.py $cat $type 20 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 21 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 22 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 23 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 24 $year >> $TMPDATFILE
#	python3 VGammaMeson_cat.py $cat $type 25 $year >> $TMPDATFILE

    fi

    if [[ "$isMC" == '2' ]]; then
	python3 VGammaMeson_cat.py $cat $type 1015 $year >> $TMPDATFILE # Zinv phi
	python3 VGammaMeson_cat.py $cat $type 1025 $year >> $TMPDATFILE # Zinv rho
	python3 VGammaMeson_cat.py $cat $type 1016 $year >> $TMPDATFILE # ggZinv phi
	python3 VGammaMeson_cat.py $cat $type 1026 $year >> $TMPDATFILE # ggZinv rho

	python3 VGammaMeson_cat.py $cat $type 1011 $year >> $TMPDATFILE #w+ Phi
	python3 VGammaMeson_cat.py $cat $type 1012 $year >> $TMPDATFILE #w-
	python3 VGammaMeson_cat.py $cat $type 1021 $year >> $TMPDATFILE #w+ Rho
	python3 VGammaMeson_cat.py $cat $type 1022 $year >> $TMPDATFILE #w-

        python3 VGammaMeson_cat.py $cat $type 1013 $year >> $TMPDATFILE  # Z-phi
        python3 VGammaMeson_cat.py $cat $type 1023 $year >> $TMPDATFILE  # Z-rho
        python3 VGammaMeson_cat.py $cat $type 1014 $year >> $TMPDATFILE  # ggZ-phi
        python3 VGammaMeson_cat.py $cat $type 1024 $year >> $TMPDATFILE  # ggZ-rho
    fi
    if [[ "$isMC" == '-1' ]]; then
	python3 VGammaMeson_cat.py $cat $type -62 $year >> $TMPDATFILE  # TauA 2018
	python3 VGammaMeson_cat.py $cat $type -63 $year >> $TMPDATFILE  # TauB 2018
	python3 VGammaMeson_cat.py $cat $type -64 $year >> $TMPDATFILE  # TauB 2018
    fi
fi


if [ $cat = "isZtag" ]; then


    if [[ "$isMC" == '2' ]]; then

        python3 VGammaMeson_cat.py $cat $type 1013 $year >> $TMPDATFILE  # Z-phi
        python3 VGammaMeson_cat.py $cat $type 1023 $year >> $TMPDATFILE  # Z-rho
        python3 VGammaMeson_cat.py $cat $type 1014 $year >> $TMPDATFILE  # ggZ-phi
        python3 VGammaMeson_cat.py $cat $type 1024 $year >> $TMPDATFILE  # ggZ-rho

    fi

    if [[ "$isMC" == '1' ]]; then

##	python3 VGammaMeson_cat.py $cat $type 1  $year >> $TMPDATFILE   # ZG
	python3 VGammaMeson_cat.py $cat $type 0  $year >> $TMPDATFILE   # DYjets
#	python3 VGammaMeson_cat.py $cat $type 34 $year >> $TMPDATFILE   # DYjets 0J
#	python3 VGammaMeson_cat.py $cat $type 35 $year >> $TMPDATFILE   # DYjets 1J
#	python3 VGammaMeson_cat.py $cat $type 36 $year >> $TMPDATFILE   # DYjets 2J

	##python3 VGammaMeson_cat.py $cat $type 31 $year >> $TMPDATFILE   # Wjets 0J
	##python3 VGammaMeson_cat.py $cat $type 32 $year >> $TMPDATFILE   # Wjets 1J
	##python3 VGammaMeson_cat.py $cat $type 33 $year >> $TMPDATFILE   # Wjets 2J
	##python3 VGammaMeson_cat.py $cat $type 2  $year >> $TMPDATFILE   # WG
	##python3 VGammaMeson_cat.py $cat $type 3  $year >> $TMPDATFILE   # Wjets
	##python3 VGammaMeson_cat.py $cat $type 4  $year >> $TMPDATFILE   # TT2l
	##python3 VGammaMeson_cat.py $cat $type 5  $year >> $TMPDATFILE   # TT1l
        ##python3 VGammaMeson_cat.py $cat $type 47 $year >> $TMPDATFILE   #TTG

    fi

    if [[ "$isMC" == '-1' ]]; then

	if [ $year = "2018" ]; then
	    python3 VGammaMeson_cat.py $cat $type -1 $year >> $TMPDATFILE  # SingleMuonA 2018
	    python3 VGammaMeson_cat.py $cat $type -2 $year >> $TMPDATFILE  # SingleMuonB 2018
	    python3 VGammaMeson_cat.py $cat $type -3 $year >> $TMPDATFILE  # SingleMuonC 2018
	    python3 VGammaMeson_cat.py $cat $type -4 $year >> $TMPDATFILE  # SingleMuonD 2018

	    python3 VGammaMeson_cat.py $cat $type -11 $year >> $TMPDATFILE  # DoubleMuonA 2018
	    python3 VGammaMeson_cat.py $cat $type -12 $year >> $TMPDATFILE  # DoubleMuonB 2018
	    python3 VGammaMeson_cat.py $cat $type -13 $year >> $TMPDATFILE  # DoubleMuonC 2018
	    python3 VGammaMeson_cat.py $cat $type -14 $year >> $TMPDATFILE  # DoubleMuonD 2018

	    python3 VGammaMeson_cat.py $cat $type -31 $year >> $TMPDATFILE  # EGammaA 2018
	    python3 VGammaMeson_cat.py $cat $type -32 $year >> $TMPDATFILE  # EGammaB 2018
	    python3 VGammaMeson_cat.py $cat $type -33 $year >> $TMPDATFILE  # EGammaC 2018
	    python3 VGammaMeson_cat.py $cat $type -34 $year >> $TMPDATFILE  # EGammaD 2018

	    python3 VGammaMeson_cat.py $cat $type -21 $year >> $TMPDATFILE  # MuonEGA 2018
	    python3 VGammaMeson_cat.py $cat $type -22 $year >> $TMPDATFILE  # MuonEGB 2018
	    python3 VGammaMeson_cat.py $cat $type -23 $year >> $TMPDATFILE  # MuonEGC 2018
	    python3 VGammaMeson_cat.py $cat $type -24 $year >> $TMPDATFILE  # MuonEGD 2018
	fi

	if [ $year = "2017" ]; then

	    python3 VGammaMeson_cat.py $cat $type -2 $year >> $TMPDATFILE  # SingleMuonB
	    python3 VGammaMeson_cat.py $cat $type -3 $year >> $TMPDATFILE  # SingleMuonC
	    python3 VGammaMeson_cat.py $cat $type -4 $year >> $TMPDATFILE  # SingleMuonD
	    python3 VGammaMeson_cat.py $cat $type -5 $year >> $TMPDATFILE  # SingleMuonE
            python3 VGammaMeson_cat.py $cat $type -6 $year >> $TMPDATFILE  # SingleMuonF

	    python3 VGammaMeson_cat.py $cat $type -12 $year >> $TMPDATFILE  # DoubleMuonB
	    python3 VGammaMeson_cat.py $cat $type -13 $year >> $TMPDATFILE  # DoubleMuonC
	    python3 VGammaMeson_cat.py $cat $type -14 $year >> $TMPDATFILE  # DoubleMuonD
	    python3 VGammaMeson_cat.py $cat $type -15 $year >> $TMPDATFILE  # DoubleMuonE
	    python3 VGammaMeson_cat.py $cat $type -16 $year >> $TMPDATFILE  # DoubleMuonF

	    python3 VGammaMeson_cat.py $cat $type -42 $year >> $TMPDATFILE  # DoubleEG 2017
	    python3 VGammaMeson_cat.py $cat $type -43 $year >> $TMPDATFILE  # DoubleEG 2017
	    python3 VGammaMeson_cat.py $cat $type -44 $year >> $TMPDATFILE  # DoubleEG 2017
	    python3 VGammaMeson_cat.py $cat $type -45 $year >> $TMPDATFILE  # DoubleEG 2017
	    python3 VGammaMeson_cat.py $cat $type -46 $year >> $TMPDATFILE  # DoubleEG 2017

	    python3 VGammaMeson_cat.py $cat $type -52 $year >> $TMPDATFILE  # SingleEle 2017-16
	    python3 VGammaMeson_cat.py $cat $type -53 $year >> $TMPDATFILE  # SingleEle 2017-16
	    python3 VGammaMeson_cat.py $cat $type -54 $year >> $TMPDATFILE  # SingleEle 2017-16
	    python3 VGammaMeson_cat.py $cat $type -55 $year >> $TMPDATFILE  # SingleEle 2017-16
	    python3 VGammaMeson_cat.py $cat $type -56 $year >> $TMPDATFILE  # SingleEle 2017-16

	    python3 VGammaMeson_cat.py $cat $type -22 $year >> $TMPDATFILE  # MuonEGB 2017
	    python3 VGammaMeson_cat.py $cat $type -23 $year >> $TMPDATFILE  # MuonEGC 2017
	    python3 VGammaMeson_cat.py $cat $type -24 $year >> $TMPDATFILE  # MuonEGD 2017
	    python3 VGammaMeson_cat.py $cat $type -25 $year >> $TMPDATFILE  # MuonEGE 2017
	    python3 VGammaMeson_cat.py $cat $type -26 $year >> $TMPDATFILE  # MuonEGF 2017
	fi

	if [ $year = "22016" ]; then

	    python3 VGammaMeson_cat.py $cat $type -6 $year >> $TMPDATFILE  # SingleMuonF
	    python3 VGammaMeson_cat.py $cat $type -7 $year >> $TMPDATFILE  # SingleMuonG
	    python3 VGammaMeson_cat.py $cat $type -8 $year >> $TMPDATFILE  # SingleMuonH

	    python3 VGammaMeson_cat.py $cat $type -16 $year >> $TMPDATFILE  # DoubleMuonF
	    python3 VGammaMeson_cat.py $cat $type -17 $year >> $TMPDATFILE  # DoubleMuonG
	    python3 VGammaMeson_cat.py $cat $type -18 $year >> $TMPDATFILE  # DoubleMuonH

	    python3 VGammaMeson_cat.py $cat $type -26 $year >> $TMPDATFILE  # DoubleEGF
	    python3 VGammaMeson_cat.py $cat $type -27 $year >> $TMPDATFILE  # DoubleEGG
	    python3 VGammaMeson_cat.py $cat $type -28 $year >> $TMPDATFILE  # DoubleEGH

	    python3 VGammaMeson_cat.py $cat $type -46 $year >> $TMPDATFILE  # DoubleElectronF
	    python3 VGammaMeson_cat.py $cat $type -47 $year >> $TMPDATFILE  # DoubleElectronG
	    python3 VGammaMeson_cat.py $cat $type -48 $year >> $TMPDATFILE  # DoubleElectronH

	    python3 VGammaMeson_cat.py $cat $type -56 $year >> $TMPDATFILE  # SingleElectronF
	    python3 VGammaMeson_cat.py $cat $type -57 $year >> $TMPDATFILE  # SingleElectronG
	    python3 VGammaMeson_cat.py $cat $type -58 $year >> $TMPDATFILE  # SingleElectronH

	fi	

	if [ $year = "12016" ]; then

	    python3 VGammaMeson_cat.py $cat $type -1 $year >> $TMPDATFILE  # SingleMuon B-ver1
	    python3 VGammaMeson_cat.py $cat $type -2 $year >> $TMPDATFILE  # SingleMuon B-ver2
	    python3 VGammaMeson_cat.py $cat $type -3 $year >> $TMPDATFILE  # SingleMuon C
	    python3 VGammaMeson_cat.py $cat $type -4 $year >> $TMPDATFILE  # SingleMuon D
	    python3 VGammaMeson_cat.py $cat $type -5 $year >> $TMPDATFILE  # SingleMuon E
	    python3 VGammaMeson_cat.py $cat $type -6 $year >> $TMPDATFILE  # SingleMuon F

	    python3 VGammaMeson_cat.py $cat $type -11 $year >> $TMPDATFILE  # DoubleMuon B-ver1
	    python3 VGammaMeson_cat.py $cat $type -12 $year >> $TMPDATFILE  # DoubleMuon B-ver2
	    python3 VGammaMeson_cat.py $cat $type -13 $year >> $TMPDATFILE  # DoubleMuon C
	    python3 VGammaMeson_cat.py $cat $type -14 $year >> $TMPDATFILE  # DoubleMuon D
	    python3 VGammaMeson_cat.py $cat $type -15 $year >> $TMPDATFILE  # DoubleMuon E
	    python3 VGammaMeson_cat.py $cat $type -16 $year >> $TMPDATFILE  # DoubleMuon F

	    python3 VGammaMeson_cat.py $cat $type -21 $year >> $TMPDATFILE  # MuonEG B-ver1
	    python3 VGammaMeson_cat.py $cat $type -22 $year >> $TMPDATFILE  # MuonEG B-ver2
	    python3 VGammaMeson_cat.py $cat $type -23 $year >> $TMPDATFILE  # MuonEG C
	    python3 VGammaMeson_cat.py $cat $type -24 $year >> $TMPDATFILE  # MuonEG D
	    python3 VGammaMeson_cat.py $cat $type -25 $year >> $TMPDATFILE  # MuonEG E
	    python3 VGammaMeson_cat.py $cat $type -26 $year >> $TMPDATFILE  # MuonEG F

	    python3 VGammaMeson_cat.py $cat $type -41 $year >> $TMPDATFILE  # DoubleEG B-ver1
	    python3 VGammaMeson_cat.py $cat $type -42 $year >> $TMPDATFILE  # DoubleEG B-ver2
	    python3 VGammaMeson_cat.py $cat $type -43 $year >> $TMPDATFILE  # DoubleEG C
	    python3 VGammaMeson_cat.py $cat $type -44 $year >> $TMPDATFILE  # DoubleEG D
	    python3 VGammaMeson_cat.py $cat $type -45 $year >> $TMPDATFILE  # DoubleEG E
	    python3 VGammaMeson_cat.py $cat $type -46 $year >> $TMPDATFILE  # DoubleEG F

	    python3 VGammaMeson_cat.py $cat $type -51 $year >> $TMPDATFILE  # SingleEle B-ver1
	    python3 VGammaMeson_cat.py $cat $type -52 $year >> $TMPDATFILE  # SingleEle B-ver2
	    python3 VGammaMeson_cat.py $cat $type -53 $year >> $TMPDATFILE  # SingleEle C
	    python3 VGammaMeson_cat.py $cat $type -54 $year >> $TMPDATFILE  # SingleEle D
	    python3 VGammaMeson_cat.py $cat $type -55 $year >> $TMPDATFILE  # SingleEle E
	    python3 VGammaMeson_cat.py $cat $type -56 $year >> $TMPDATFILE  # SingleEle F
	fi
    fi
fi

if [ $cat = "isWtag" ]; then

    if [[ "$isMC" == '-1' ]]; then

	if [ $year = "2018" ]; then
	    python3 VGammaMeson_cat.py $cat $type -1 $year >> $TMPDATFILE  # SingleMuonA 2018
	    python3 VGammaMeson_cat.py $cat $type -2 $year >> $TMPDATFILE  # SingleMuonB 2018
	    python3 VGammaMeson_cat.py $cat $type -3 $year >> $TMPDATFILE  # SingleMuonC 2018
	    python3 VGammaMeson_cat.py $cat $type -4 $year >> $TMPDATFILE  # SingleMuonD 2018

	    python3 VGammaMeson_cat.py $cat $type -11 $year >> $TMPDATFILE  # DoubleMuonA 2018
	    python3 VGammaMeson_cat.py $cat $type -12 $year >> $TMPDATFILE  # DoubleMuonB 2018
	    python3 VGammaMeson_cat.py $cat $type -13 $year >> $TMPDATFILE  # DoubleMuonC 2018
	    python3 VGammaMeson_cat.py $cat $type -14 $year >> $TMPDATFILE  # DoubleMuonD 2018

	    python3 VGammaMeson_cat.py $cat $type -31 $year >> $TMPDATFILE  # EGammaA 2018
	    python3 VGammaMeson_cat.py $cat $type -32 $year >> $TMPDATFILE  # EGammaB 2018
	    python3 VGammaMeson_cat.py $cat $type -33 $year >> $TMPDATFILE  # EGammaC 2018
	    python3 VGammaMeson_cat.py $cat $type -34 $year >> $TMPDATFILE  # EGammaD 2018

	    python3 VGammaMeson_cat.py $cat $type -21 $year >> $TMPDATFILE  # MuonEGA 2018
	    python3 VGammaMeson_cat.py $cat $type -22 $year >> $TMPDATFILE  # MuonEGB 2018
	    python3 VGammaMeson_cat.py $cat $type -23 $year >> $TMPDATFILE  # MuonEGC 2018
	    python3 VGammaMeson_cat.py $cat $type -24 $year >> $TMPDATFILE  # MuonEGD 2018
	fi

	if [ $year = "2017" ]; then
	    python3 VGammaMeson_cat.py $cat $type -2 $year >> $TMPDATFILE  # SingleMuonB
	    python3 VGammaMeson_cat.py $cat $type -3 $year >> $TMPDATFILE  # SingleMuonC
	    python3 VGammaMeson_cat.py $cat $type -4 $year >> $TMPDATFILE  # SingleMuonD
	    python3 VGammaMeson_cat.py $cat $type -5 $year >> $TMPDATFILE  # SingleMuonE
	    python3 VGammaMeson_cat.py $cat $type -6 $year >> $TMPDATFILE  # SingleMuonF

	    python3 VGammaMeson_cat.py $cat $type -12 $year >> $TMPDATFILE  # DoubleMuonB
	    python3 VGammaMeson_cat.py $cat $type -13 $year >> $TMPDATFILE  # DoubleMuonC
	    python3 VGammaMeson_cat.py $cat $type -14 $year >> $TMPDATFILE  # DoubleMuonD
	    python3 VGammaMeson_cat.py $cat $type -15 $year >> $TMPDATFILE  # DoubleMuonE
	    python3 VGammaMeson_cat.py $cat $type -16 $year >> $TMPDATFILE  # DoubleMuonF

	    python3 VGammaMeson_cat.py $cat $type -42 $year >> $TMPDATFILE  # DoubleEG 2017
	    python3 VGammaMeson_cat.py $cat $type -43 $year >> $TMPDATFILE  # DoubleEG 2017
	    python3 VGammaMeson_cat.py $cat $type -44 $year >> $TMPDATFILE  # DoubleEG 2017
	    python3 VGammaMeson_cat.py $cat $type -45 $year >> $TMPDATFILE  # DoubleEG 2017
	    python3 VGammaMeson_cat.py $cat $type -46 $year >> $TMPDATFILE  # DoubleEG 2017

	    python3 VGammaMeson_cat.py $cat $type -52 $year >> $TMPDATFILE  # SingleEle 2017-16
	    python3 VGammaMeson_cat.py $cat $type -53 $year >> $TMPDATFILE  # SingleEle 2017-16
	    python3 VGammaMeson_cat.py $cat $type -54 $year >> $TMPDATFILE  # SingleEle 2017-16
	    python3 VGammaMeson_cat.py $cat $type -55 $year >> $TMPDATFILE  # SingleEle 2017-16
	    python3 VGammaMeson_cat.py $cat $type -56 $year >> $TMPDATFILE  # SingleEle 2017-16

	    python3 VGammaMeson_cat.py $cat $type -22 $year >> $TMPDATFILE  # MuonEGB 2017
	    python3 VGammaMeson_cat.py $cat $type -23 $year >> $TMPDATFILE  # MuonEGC 2017
	    python3 VGammaMeson_cat.py $cat $type -24 $year >> $TMPDATFILE  # MuonEGD 2017
	    python3 VGammaMeson_cat.py $cat $type -25 $year >> $TMPDATFILE  # MuonEGE 2017
	    python3 VGammaMeson_cat.py $cat $type -26 $year >> $TMPDATFILE  # MuonEGF 2017
	fi

	if [ $year = "22016" ]; then
	    python3 VGammaMeson_cat.py $cat $type -6 $year >> $TMPDATFILE  # SingleMuonF
	    python3 VGammaMeson_cat.py $cat $type -7 $year >> $TMPDATFILE  # SingleMuonG
	    python3 VGammaMeson_cat.py $cat $type -8 $year >> $TMPDATFILE  # SingleMuonH

	    python3 VGammaMeson_cat.py $cat $type -16 $year >> $TMPDATFILE  # DoubleMuonF
	    python3 VGammaMeson_cat.py $cat $type -17 $year >> $TMPDATFILE  # DoubleMuonG
	    python3 VGammaMeson_cat.py $cat $type -18 $year >> $TMPDATFILE  # DoubleMuonH

            python3 VGammaMeson_cat.py $cat $type -26 $year >> $TMPDATFILE  # MuonEGF
	    python3 VGammaMeson_cat.py $cat $type -27 $year >> $TMPDATFILE  # MuonEGG
	    python3 VGammaMeson_cat.py $cat $type -28 $year >> $TMPDATFILE  # MuonEGH

	    python3 VGammaMeson_cat.py $cat $type -46 $year >> $TMPDATFILE  # DoubleEGF
	    python3 VGammaMeson_cat.py $cat $type -47 $year >> $TMPDATFILE  # DoubleEGG
	    python3 VGammaMeson_cat.py $cat $type -48 $year >> $TMPDATFILE  # DoubleEGH

	    python3 VGammaMeson_cat.py $cat $type -56 $year >> $TMPDATFILE  # SingleElectronF
	    python3 VGammaMeson_cat.py $cat $type -57 $year >> $TMPDATFILE  # SingleElectronG
	    python3 VGammaMeson_cat.py $cat $type -58 $year >> $TMPDATFILE  # SingleElectronH
	fi

	if [ $year = "12016" ]; then
	    python3 VGammaMeson_cat.py $cat $type -1 $year >> $TMPDATFILE  # SingleMuon B-ver1
	    python3 VGammaMeson_cat.py $cat $type -2 $year >> $TMPDATFILE  # SingleMuon B-ver2
	    python3 VGammaMeson_cat.py $cat $type -3 $year >> $TMPDATFILE  # SingleMuon C
	    python3 VGammaMeson_cat.py $cat $type -4 $year >> $TMPDATFILE  # SingleMuon D
	    python3 VGammaMeson_cat.py $cat $type -5 $year >> $TMPDATFILE  # SingleMuon E
	    python3 VGammaMeson_cat.py $cat $type -6 $year >> $TMPDATFILE  # SingleMuon F

	    python3 VGammaMeson_cat.py $cat $type -11 $year >> $TMPDATFILE  # DoubleMuon B-ver1
	    python3 VGammaMeson_cat.py $cat $type -12 $year >> $TMPDATFILE  # DoubleMuon B-ver2
	    python3 VGammaMeson_cat.py $cat $type -13 $year >> $TMPDATFILE  # DoubleMuon C
	    python3 VGammaMeson_cat.py $cat $type -14 $year >> $TMPDATFILE  # DoubleMuon D
	    python3 VGammaMeson_cat.py $cat $type -15 $year >> $TMPDATFILE  # DoubleMuon E
	    python3 VGammaMeson_cat.py $cat $type -16 $year >> $TMPDATFILE  # DoubleMuon F

	    python3 VGammaMeson_cat.py $cat $type -21 $year >> $TMPDATFILE  # MuonEG B-ver1
	    python3 VGammaMeson_cat.py $cat $type -22 $year >> $TMPDATFILE  # MuonEG B-ver2
	    python3 VGammaMeson_cat.py $cat $type -23 $year >> $TMPDATFILE  # MuonEG C
	    python3 VGammaMeson_cat.py $cat $type -24 $year >> $TMPDATFILE  # MuonEG D
	    python3 VGammaMeson_cat.py $cat $type -25 $year >> $TMPDATFILE  # MuonEG E
	    python3 VGammaMeson_cat.py $cat $type -26 $year >> $TMPDATFILE  # MuonEG F

	    python3 VGammaMeson_cat.py $cat $type -42 $year >> $TMPDATFILE  # DoubleEG B-ver1
	    python3 VGammaMeson_cat.py $cat $type -42 $year >> $TMPDATFILE  # DoubleEG B-ver2
	    python3 VGammaMeson_cat.py $cat $type -43 $year >> $TMPDATFILE  # DoubleEG C
	    python3 VGammaMeson_cat.py $cat $type -44 $year >> $TMPDATFILE  # DoubleEG D
	    python3 VGammaMeson_cat.py $cat $type -45 $year >> $TMPDATFILE  # DoubleEG E
	    python3 VGammaMeson_cat.py $cat $type -46 $year >> $TMPDATFILE  # DoubleEG F

	    python3 VGammaMeson_cat.py $cat $type -52 $year >> $TMPDATFILE  # SingleEle B-ver1
	    python3 VGammaMeson_cat.py $cat $type -53 $year >> $TMPDATFILE  # SingleEle B-ver2
	    python3 VGammaMeson_cat.py $cat $type -54 $year >> $TMPDATFILE  # SingleEle C
	    python3 VGammaMeson_cat.py $cat $type -55 $year >> $TMPDATFILE  # SingleEle D
	    python3 VGammaMeson_cat.py $cat $type -56 $year >> $TMPDATFILE  # SingleEle E
	    python3 VGammaMeson_cat.py $cat $type -56 $year >> $TMPDATFILE  # SingleEle F
	fi
    fi

    if [[ "$isMC" == '2' ]]; then

	python3 VGammaMeson_cat.py $cat $type 1011 $year >> $TMPDATFILE #w+ Phi
	python3 VGammaMeson_cat.py $cat $type 1012 $year >> $TMPDATFILE #w-
	python3 VGammaMeson_cat.py $cat $type 1013 $year >> $TMPDATFILE #zh
	python3 VGammaMeson_cat.py $cat $type 1014 $year >> $TMPDATFILE #ggzh
	python3 VGammaMeson_cat.py $cat $type 1021 $year >> $TMPDATFILE #w+ RHo
	python3 VGammaMeson_cat.py $cat $type 1022 $year >> $TMPDATFILE #w-
	python3 VGammaMeson_cat.py $cat $type 1023 $year >> $TMPDATFILE #zh
	python3 VGammaMeson_cat.py $cat $type 1024 $year >> $TMPDATFILE #ggzh

    fi

    if [[ "$isMC" == '1' ]]; then

#	python3 VGammaMeson_cat.py $cat $type 2  $year >> $TMPDATFILE  # WG
#	python3 VGammaMeson_cat.py $cat $type 1  $year >> $TMPDATFILE  # ZG
#	python3 VGammaMeson_cat.py $cat $type 47 $year >> $TMPDATFILE  #TTG

#	python3 VGammaMeson_cat.py $cat $type 3  $year >> $TMPDATFILE  # Wjets
#	python3 VGammaMeson_cat.py $cat $type 32 $year >> $TMPDATFILE # Wjets 1J
#	python3 VGammaMeson_cat.py $cat $type 33 $year >> $TMPDATFILE # Wjets 2J
#	python3 VGammaMeson_cat.py $cat $type 31 $year >> $TMPDATFILE # Wjets 0J
	python3 VGammaMeson_cat.py $cat $type 0  $year >> $TMPDATFILE # DYjets
#	python3 VGammaMeson_cat.py $cat $type 34 $year >> $TMPDATFILE # DYjets 0J
#	python3 VGammaMeson_cat.py $cat $type 35 $year >> $TMPDATFILE # DYjets 1J
#	python3 VGammaMeson_cat.py $cat $type 36 $year >> $TMPDATFILE # DYjets 2J
#	python3 VGammaMeson_cat.py $cat $type 4  $year >> $TMPDATFILE  # TT2l
#	python3 VGammaMeson_cat.py $cat $type 5  $year >> $TMPDATFILE  # TT1l  (BROKEN SAMPLE - 2018)

    fi
fi

done
done
done
