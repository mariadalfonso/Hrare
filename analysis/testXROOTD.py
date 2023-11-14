import ROOT
import json
import os
import sys
from XRootD import client
import glob
from subprocess import call,check_output

def findDIR(directory,useXROOTD=False):

    print(directory)
    counter = 0
    rootFiles = ROOT.vector('string')()
    maxFiles = 1000000000

    if(useXROOTD == True and "/data/submit/cms" in directory):
        xrd = "root://submit50.mit.edu/"
        xrdpath = directory.replace("/data/submit/cms","")
        f = check_output(['xrdfs', f'{xrd}', 'ls', xrdpath]).decode(sys.stdout.encoding)
        stringFiles = f.split()
        for e in range(len(stringFiles)):
            filePath = xrd + stringFiles[e]
            if "failed/" in filePath: continue
            if "log/" in filePath: continue
            if ".txt" in filePath: continue
            counter+=1
            if(counter > maxFiles): break
            rootFiles.push_back(filePath)
    else:
        filename = ("{}".format(directory)+"/*.root")
        for filenames in glob.glob(filename):
            counter+=1
            if(counter > maxFiles): break
            rootFiles.push_back(filenames)

    print(len(rootFiles))
    return rootFiles


useXROOTD=True

ROOT.ROOT.EnableImplicitMT()
RDataFrame = ROOT.RDataFrame

BARRELphotons = "(Photon_pt>30 and Photon_isScEtaEB and Photon_mvaID_WP90)"
ENDCAPphotons = "(Photon_pt>30 and Photon_isScEtaEE and Photon_mvaID_WP80)"

year = 2018
GOODphotons = "({0} or {1}) and Photon_pt>38 and Photon_electronVeto and abs(Photon_eta)<2.1".format(BARRELphotons,ENDCAPphotons) #90-80

def analysis():
    files = findDIR("/data/submit/cms/store/user/mariadlf/nano/D02/GluGlu_HToRhoGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3+MINIAODSIM/",useXROOTD)

    dfINI = RDataFrame("Events", files)
    print(GOODphotons)

    df = (dfINI
          .Filter("nPhoton>0 and PV_npvsGood>0","photon from nano >0 and PV_npvsGood > 0")
          .Define("w","{}".format(1))
          .Define("goodPhotons", "{}".format(GOODphotons))
          .Define("nGoodPhotons","Sum(goodPhotons)*1.0f")
          .Filter("Sum(goodPhotons)>0", "At least one good Photon")
          .Define("goodPhotons_pt", "Photon_pt[goodPhotons][0]")
          .Define("goodPhotons_eta", "Photon_eta[goodPhotons][0]")
          )

    ptsum = df.Sum("Muon_pt")
    res = ptsum.GetValue()
    print('sum of muons = ',res)

# In a Python script the Dask client needs to be initalized in a context
# Jupyter notebooks / Python session don't need this
if __name__ == "__main__":

    analysis()
