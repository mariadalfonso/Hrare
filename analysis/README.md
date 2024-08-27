## set environment with	singularity

> APPTAINER_BIND="/home/submit,/work/submit,/data/submit,/scratch/submit,/cvmfs" singularity shell /cvmfs/unpacked.cern.ch/registry.hub.docker.com/rootproject/root:latest


## set environment with conda

> #wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

> #bash Miniforge3-Linux-x86_64.sh

> conda create --name myenv

> conda activate myenv

> conda install -c conda-forge root

> python3 -m pip install --user --no-binary=correctionlib correctionlib

> python3 -m pip install distributed

> python3 -m pip install dask-jobqueue

> python3 -m pip install dask_gateway

## run RDF on notebook
testDistCluster.ipynb

## run distributed test
python testDistCluster.py


## set up VOMS
voms-proxy-init  --voms cms -valid 198:0

## notes
#1. access to cvmfs is needed i.e. /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/

#2. we can activate the environment with this: conda activate myenv

## standalone root
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch

source $VO_CMS_SW_DIR/cmsset_default.sh

source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.28.02/x86_64-centos8-gcc85-opt/bin/thisroot.sh
