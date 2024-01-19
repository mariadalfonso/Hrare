## set environment

#wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

#bash Miniforge3-Linux-x86_64.sh

conda create --name myenv

conda install -c conda-forge root

python3 -m pip install --user --no-binary=correctionlib correctionlib

python3 -m pip install dask-jobqueue

## run RDF on notebook
testDistCluster.ipynb

## run distributed test
python testDistLocalCluster.py

## test standalone xrootd
python testXROOTD.py

## set up VOMS
voms-proxy-init  --voms cms -valid 198:0

## notes
#1. access to cvmfs is needed i.e. /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/

#2. we can activate the environment with this: conda activate myenv

#3. some librarty is pre-compiled i.e. g++ $(correction config --cflags --ldflags --rpath) mysf.cc -shared -fPIC -o mysf.so 
or python utilsHrare.py pre-compile the functions_cc.so and functions_cc_ACLiC_dict_rdict.pcm

## standalone root
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch

source $VO_CMS_SW_DIR/cmsset_default.sh

source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.28.02/x86_64-centos8-gcc85-opt/bin/thisroot.sh
