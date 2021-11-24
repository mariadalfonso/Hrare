import ROOT
import os
import json
from subprocess import call,check_output
import fnmatch
import math

from utilsHrare import findDIR

TRIGGERvbf = "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3 or HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3 or HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50"

TRIGGERincl = "HLT_Photon35_TwoProngs35"

def groupFiles(fIns, group):

    filesPerGroup = math.ceil(len(fIns)/group)
    ret = []
    for i in range(0, group):
    
        a = i*filesPerGroup
        b = (i+1)*filesPerGroup
        subFiles = fIns[a:b]
        ret.append(subFiles)
    
    return ret

if __name__ == "__main__":

    if True: ### skim TTToSemiLeptonic
    
        group = 10
        fOutDir = "/work/submit/mariadlf/SKIMS/D00/GJets_DR-0p4_HT-100To200_Summer20UL18"
        files =findDIR("/mnt/T2_US_MIT/hadoop/cms/store/user/paus/nanohr/D00/GJets_DR-0p4_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM")
        print(len(files))

        groupedFiles = groupFiles(files, group)

        for i, group in enumerate(groupedFiles):
            
            fOutName = "%s/output_%d.root" % (fOutDir, i)
            print("Create %s" % fOutName)

            rdf = ROOT.ROOT.RDataFrame("Events", group).Define("trigger","{}".format(TRIGGERvbf)+" or {}".format(TRIGGERincl)).Filter("trigger>0").Snapshot("Events", fOutName)

            del rdf
