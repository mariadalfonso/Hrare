import ROOT
import os, sys, getopt, json, time, subprocess, socket
from subprocess import call,check_output
import fnmatch
import math

ROOT.ROOT.EnableImplicitMT()
ROOT.ROOT.EnableThreadSafety()

sys.path.insert(0, '/home/submit/mariadlf/Hrare/CMSSW_10_6_27/src/Hrare/analysis')

#from utilsHrare import findManyXRDFS, findManyClient, readListFromFile
from utilsHrare import readListFromFile
from utilsHrare import pickTRG, getSkimsFromJson

with open("/home/submit/mariadlf/Hrare/CMSSW_10_6_27/src/Hrare/analysis/config/selection.json") as jsonFile:
#with open("selection.json") as jsonFile:
    jsonObject = json.load(jsonFile)
    jsonFile.close()

with open("/home/submit/mariadlf/Hrare/CMSSW_10_6_27/src/Hrare/analysis/config/trigger.json") as trgJsonFile:
    trgObject = json.load(trgJsonFile)
    trgJsonFile.close()

overall = trgObject['triggers']

with open("/home/submit/mariadlf/Hrare/CMSSW_10_6_27/src/Hrare/analysis/config/skimDB.json") as jsonFile:
#with open("skimDB.json") as jsonFile:
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
#        print(subFiles)
        ret.append(subFiles)
    
    return ret

if __name__ == "__main__":

    whichJob = -1

    valid = ['year=', "era=", "PDType=", "SkimType=","whichJob=","whichFile="]
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
        if opt == "--whichFile":
            whichFile = str(arg)

    print("year=",year," era=",era," PDType=",PDType," skimType=",skimType," whichJob=",whichJob," whichFile=",whichFile)

    if True:
        isVBF = False
        isW = False
        isZ = False
        isZinv = False

        if skimType== "VBF": isVBF=True
        if skimType== "Zinv": isZinv=True
        if skimType== "VH": isW=True
        if skimType== "VH": isZ=True

        checkFile = "catalog/files_"+str(PDType)+"_"+str(era)+"_"+str(year)+".txt"
        print("checkFile",checkFile)
        print("whichFile",whichFile)

        if checkFile!=whichFile:
            print("FILE MISMATCH")

        fOutDir = "/scratch/submit/cms/mariadlf/Hrare/newSKIMS/D01/"+skimType+"/"+str(year)+"/"+PDType+"+Run"+era
        fInDir = "/mnt/T2_US_MIT/hadoop/cms/store/user/paus/nanohr/D01/"
#        fInDir = "xrdfs root://xrootd.cmsaf.mit.edu ls /store/user/paus/nanohr/D01"
        datasetExp = ("/"+str(PDType)+"+Run"+str(year)+str(era)+"*")
        files = readListFromFile(whichFile)
        print("number of INPUT files = ",len(files))

        TRIGGER=pickTRG(overall,year,PDType,isVBF,isW,isZ,isZinv)
        print("TRIGGER=",TRIGGER)

        PRESELECTION=getSkimsFromJson(skims, skimType)
        print("PRESELECTION",PRESELECTION)

        group = math.ceil((len(files)/1))
        print("number of output files = ", group)
        groupedFiles = groupFiles(files, group)

        for i, group in enumerate(groupedFiles):
            passJob = whichJob == -1 or whichJob == i
            if(passJob == False): continue

            fOutName = "%s/output_%s_%d_%d.root" % (fOutDir, PDType, year, i)
            print("Create %s" % fOutName)
            print("input file %s"  % group)

            regexToDrop = "^(?!.*omega).*$"
            rdf = ROOT.ROOT.RDataFrame("Events", group).Define("trigger","{}".format(TRIGGER)).Filter("{}".format(PRESELECTION)).Snapshot("Events", fOutName, regexToDrop)

            del rdf
