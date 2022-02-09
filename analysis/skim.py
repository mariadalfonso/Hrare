import ROOT
import os
import json
from subprocess import call,check_output
import fnmatch
import math
import json

ROOT.ROOT.EnableImplicitMT()
ROOT.ROOT.EnableThreadSafety()

from utilsHrare import findDIR, findMany
from utilsHrare import getTriggerFromJson, pickTRG

with open("config/selection.json") as jsonFile:
    jsonObject = json.load(jsonFile)
    jsonFile.close()

overall = jsonObject['triggers']

#TRIGGERvbf = "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3 or HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3 or HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50"
#TRIGGERincl = "HLT_Photon35_TwoProngs35"

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

    if True:

        type="VH"
        datasetNumber = -2
        year = 2018
        era = "B"
        PDType =  "SingleMuon"
        PRESELECTION = "(trigger>0 and nPhoton>0)"
        isVBF = False
        isW = True
        isZ = True

        fOutDir = "/work/submit/mariadlf/SKIMS/D01/"+type+"/"+str(year)+"/"+PDType+"+Run"+era
        fInDir = "/mnt/T2_US_MIT/hadoop/cms/store/user/paus/nanohr/D01/"
        datasetExp = (PDType+"+Run"+str(year)+era+"*")
        files = findMany(fInDir, datasetExp)
        print("number of INPUT files = ",len(files))

        TRIGGER=pickTRG(overall,year,PDType,isVBF,isW,isZ)
        print(TRIGGER)

        group = math.ceil((len(files)/1)) # one per files otherwise need to handle with hadd
        print("number of output files = ", group)
        groupedFiles = groupFiles(files, group)

        for i, group in enumerate(groupedFiles):
            
            fOutName = "%s/output_%s_%d.root" % (fOutDir, PDType, i)
            print("Create %s" % fOutName)

            regexToDrop = "^(?!.*omega).*$"
            rdf = ROOT.ROOT.RDataFrame("Events", group).Define("trigger","{}".format(TRIGGER)).Filter("{}".format(PRESELECTION)).Snapshot("Events", fOutName, regexToDrop)

            del rdf
