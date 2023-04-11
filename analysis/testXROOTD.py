import ROOT
import json
import os
from XRootD import client
import glob

def findDIR(directory,useXROOTD=False):

    print(directory)
    counter = 0
    rootFiles = ROOT.vector('string')()
    maxFiles = 1000000000

    if(useXROOTD == True and "/data/submit/cms" in directory):
        print(directory)
        fs = client.FileSystem('root://submit50.mit.edu/')
        lsst = fs.dirlist(directory.replace("/data/submit/cms",""))
        for e in lsst[1]:
            filePath = os.path.join(directory.replace("/data/submit/cms","root://submit50.mit.edu/"),e.name)
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
    files = findDIR("/data/submit/cms/store/user/mariadlf/nano/D01/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8+RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1+MINIAODSIM/",useXROOTD)

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
