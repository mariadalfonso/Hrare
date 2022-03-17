import ROOT
import os, sys, getopt, json, time, subprocess, socket
from subprocess import call,check_output
import fnmatch
import math

ROOT.ROOT.EnableImplicitMT()
ROOT.ROOT.EnableThreadSafety()

from utilsHrare import findMany
from utilsHrare import pickTRG, getSkimsFromJson

with open("config/selection.json") as jsonFile:
    jsonObject = json.load(jsonFile)
    jsonFile.close()

overall = jsonObject['triggers']

with open("config/skimDB.json") as jsonFile:
    skimObject = json.load(jsonFile)
    jsonFile.close()

skims = skimObject['skims']

def groupFiles(fIns, group):

    filesPerGroup = math.ceil(len(fIns)/group)
    ret = []
    for i in range(0, group):
    
        a = i*filesPerGroup
        b = (i+1)*filesPerGroup
        subFiles = fIns[a:b]
        print(subFiles)
        ret.append(subFiles)
    
    return ret

if __name__ == "__main__":

    whichJob = -1

    valid = ['year=', "era=", "PDType=", "SkimType=","whichJob="]
    opts, args = getopt.getopt(sys.argv[1:], "", valid)

    for opt, arg in opts:
        if opt == "--year":
            year = int(arg)
        if opt == "--era":
            era = str(arg)
        if opt == "--PDType":
            PDType = str(arg)
        if opt == "--SkimType":
            skimType = str(arg)
        if opt == "--whichJob":
            whichJob = int(arg)

    print("year=",year," era=",era," PDType=",PDType," skimType=",skimType," whichJob=",whichJob)

    if True:
        isVBF = False
        isW = False
        isZ = False

        if skimType== "VBF": isVBF=True
        if skimType== "VH": isW=True
        if skimType== "VH": isZ=True

        fOutDir = "/scratch/submit/cms/mariadlf/Hrare/SKIMS/D01/"+skimType+"/"+str(year)+"/"+PDType+"+Run"+era
        fInDir = "/mnt/T2_US_MIT/hadoop/cms/store/user/paus/nanohr/D01/"
        datasetExp = (str(PDType)+"+Run"+str(year)+str(era)+"*")
        files = findMany(fInDir, datasetExp)
        print("number of INPUT files = ",len(files))

        TRIGGER=pickTRG(overall,year,PDType,isVBF,isW,isZ)
        print("TRIGGER=",TRIGGER)

        PRESELECTION=getSkimsFromJson(skims, skimType)
        print("PRESELECTION",PRESELECTION)

        group = math.ceil((len(files)/1))
        print("number of output files = ", group)
        groupedFiles = groupFiles(files, group)

        for i, group in enumerate(groupedFiles):
            passJob = whichJob == -1 or whichJob == i
            if(passJob == False): continue

            fOutName = "%s/output_%s_%d.root" % (fOutDir, PDType, i)
            print("Create %s" % fOutName)

            regexToDrop = "^(?!.*omega).*$"
            rdf = ROOT.ROOT.RDataFrame("Events", group).Define("trigger","{}".format(TRIGGER)).Filter("{}".format(PRESELECTION)).Snapshot("Events", fOutName, regexToDrop)

            del rdf
