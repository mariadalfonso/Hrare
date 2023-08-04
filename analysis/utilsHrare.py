import ROOT
import os
import json
import fnmatch
import glob
from XRootD import client

if "/home/submit/mariadlf/Hrare/CMSSW_10_6_27_new/src/Hrare/analysis/functions.so" not in ROOT.gSystem.GetLibraries():
    ROOT.gSystem.CompileMacro("/home/submit/mariadlf/Hrare/CMSSW_10_6_27_new/src/Hrare/analysis/functions.cc","k")

import correctionlib
correctionlib.register_pyroot_binding()

def loadPolarizationTree(mc,year):
    ROOT.gInterpreter.ProcessLine('initPol(%d,%d);'% (mc,year))

def loadCorrectionSet(year):
    ROOT.gInterpreter.Load("config/mysf.so")
    ROOT.gInterpreter.Declare('#include "config/mysf.h"')
    ROOT.gInterpreter.ProcessLine('auto corr_sf = MyCorrections(%d);' % (year))
    ROOT.gInterpreter.Declare('''
        #ifndef MYFUN
        #define MYFUN
        Vec_f computeJECuncertainties(MyCorrections corrSFs, Vec_f jet_pt, Vec_f jet_eta){
        Vec_f new_jet_delta(jet_pt.size(), 1.0);
        int type = 0;
        for (unsigned int idx = 0; idx < jet_pt.size(); ++idx) new_jet_delta[idx] = corrSFs.eval_jesUnc(jet_eta[idx], jet_pt[idx], type );
        return new_jet_delta;
        }
        #endif
        ''')

def getSkimsFromJson(overall, type ):

    for selection in overall:
        if selection['type'] == type : return selection['PRESELECTION']

def getTriggerFromJson(overall, type, year ):

    for trigger in overall:
        if trigger['name'] == type and trigger['year'] == year: return trigger['definition']

def getMesonFromJson(overall, type, cat ):

    for meson in overall:
        if meson['name'] == type and meson['type'] == cat: return meson['definition']

def getMVAFromJson(overall, type, cat ):

    for meson in overall:
        if meson['name'] == type and meson['type'] == cat: return meson['file']

def loadJSON(fIn):

    if not os.path.isfile(fIn):
        print("JSON file %s does not exist" % fIn)
        return

    if not hasattr(ROOT, "jsonMap"):
        print("jsonMap not found in ROOT dict")
        return

    info = json.load(open(fIn))
    print("JSON file %s loaded" % fIn)
    for k,v in info.items():

        vec = ROOT.std.vector["std::pair<unsigned int, unsigned int>"]()
        for combo in v:
            pair = ROOT.std.pair["unsigned int", "unsigned int"](*[int(c) for c in combo])
            vec.push_back(pair)
            ROOT.jsonMap[int(k)] = vec

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

def readListFromFile(filename):

    rootFiles = ROOT.vector('string')()
    with open(filename, "r") as f:
        for item in f:
            rootFiles.push_back(item.rstrip())

    return rootFiles

def findMany(basedir, regex):

    if basedir[-1] == "/": basedir = basedir[:-1]
    regex = basedir + "/" + regex

    print("regex=",regex)
    print("basedir=",basedir)

    counter = 0

    rootFiles = ROOT.vector('string')()

    for root, directories, filenames in os.walk(basedir):

        for f in filenames:

            filePath = os.path.join(os.path.abspath(root), f)
            if "failed/" in filePath: continue
            if "log/" in filePath: continue
            if fnmatch.fnmatch(filePath, regex): rootFiles.push_back(filePath)
            if fnmatch.fnmatch(filePath, regex): counter+=1
            #            if counter>10: break
            #            if counter>250: break
    return rootFiles

def getMClist(year, sampleNOW):

    files = findDIR("{}".format(SwitchSample(sampleNOW,year)[0]))
    return files

def getDATAlist(argument,year):

    dirJson = "config"

    if(year == 2018):
        loadJSON("{}/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt".format(dirJson))
    if(year == 2017):
        loadJSON("{}/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt".format(dirJson))
    if(year == 22016 or year == 12016):
        loadJSON("{}/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt".format(dirJson))

    dirT2 = "/mnt/T2_US_MIT/hadoop/cms/store/user/paus/nanohr/D02/"

# FIXME: this build toooo many lists

    switch = {
        -1: (("SingleMuon+Run"+str(year)+"A*"),"SingleMuon"),
        -2: (("SingleMuon+Run"+str(year)+"B*"),"SingleMuon"),
        -3: (("SingleMuon+Run"+str(year)+"C*"),"SingleMuon"),
        -4: (("SingleMuon+Run"+str(year)+"D*"),"SingleMuon"),
        #
        -11: (("DoubleMuon+Run"+str(year)+"A*"),"DoubleMuon"),
        -12: (("DoubleMuon+Run"+str(year)+"B*"),"DoubleMuon"),
        -13: (("DoubleMuon+Run"+str(year)+"C*"),"DoubleMuon"),
        -14: (("DoubleMuon+Run"+str(year)+"D*"),"DoubleMuon"),
        #
        -21: (("MuonEG+Run"+str(year)+"A*"),"MuonEG"),
        -22: (("MuonEG+Run"+str(year)+"B*"),"MuonEG"),
        -23: (("MuonEG+Run"+str(year)+"C*"),"MuonEG"),
        -24: (("MuonEG+Run"+str(year)+"D*"),"MuonEG"),
        #
        -31: (("EGamma+Run"+str(year)+"A*"),"EGamma"),
        -32: (("EGamma+Run"+str(year)+"B*"),"EGamma"),
        -33: (("EGamma+Run"+str(year)+"C*"),"EGamma"),
        -34: (("EGamma+Run"+str(year)+"D*"),"EGamma"),
        #  SinglePhoton (for the 2017 VBF)
        -52: (("SinglePhoton+Run"+str(year)+"B*"),"SinglePhoton"),
        -53: (("SinglePhoton+Run"+str(year)+"C*"),"SinglePhoton"),
        -54: (("SinglePhoton+Run"+str(year)+"D*"),"SinglePhoton"),
        -55: (("SinglePhoton+Run"+str(year)+"E*"),"SinglePhoton"),
        -56: (("SinglePhoton+Run"+str(year)+"F*"),"SinglePhoton"),
        # Tau
        -61: (("Tau+Run"+str(year)+"A*"),"Tau"),
        -62: (("Tau+Run"+str(year)+"B*"),"Tau"),
        -63: (("Tau+Run"+str(year)+"C*"),"Tau"),
        -64: (("Tau+Run"+str(year)+"D*"),"Tau"),
    }

    tmp_pair = switch.get(argument, "regex, PDtype")
    pair = (findMany(dirT2,tmp_pair[0]),tmp_pair[1])

    return pair

def getSkims(argument,year,category):

    dirJson = "config"

    if(year == 2018):
        loadJSON("{}/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt".format(dirJson))
    if(year == 2017):
        loadJSON("{}/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt".format(dirJson))
    if(year == 22016 or year == 12016):
        loadJSON("{}/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt".format(dirJson))

    dirScratch02 = "/scratch/submit/cms/mariadlf/Hrare/newSKIMS/D02"
    dirScratchSS = "/scratch/submit/cms/mariadlf/Hrare/SStau"

    switch = {
        -1: (("SingleMuon+Run"+"A"),"SingleMuon"),
        -2: (("SingleMuon+Run"+"B"),"SingleMuon"),
        -3: (("SingleMuon+Run"+"C"),"SingleMuon"),
        -4: (("SingleMuon+Run"+"D"),"SingleMuon"),
        -5: (("SingleMuon+Run"+"E"),"SingleMuon"),
        -6: (("SingleMuon+Run"+"F"),"SingleMuon"),
        -7: (("SingleMuon+Run"+"G"),"SingleMuon"),
        -8: (("SingleMuon+Run"+"H"),"SingleMuon"),
        #
        -11: (("DoubleMuon+Run"+"A"),"DoubleMuon"),
        -12: (("DoubleMuon+Run"+"B"),"DoubleMuon"),
        -13: (("DoubleMuon+Run"+"C"),"DoubleMuon"),
        -14: (("DoubleMuon+Run"+"D"),"DoubleMuon"),
        -15: (("DoubleMuon+Run"+"E"),"DoubleMuon"),
        -16: (("DoubleMuon+Run"+"F"),"DoubleMuon"),
        -17: (("DoubleMuon+Run"+"G"),"DoubleMuon"),
        -18: (("DoubleMuon+Run"+"H"),"DoubleMuon"),
        #
        -21: (("MuonEG+Run"+"A"),"MuonEG"),
        -22: (("MuonEG+Run"+"B"),"MuonEG"),
        -23: (("MuonEG+Run"+"C"),"MuonEG"),
        -24: (("MuonEG+Run"+"D"),"MuonEG"),
        -25: (("MuonEG+Run"+"E"),"MuonEG"),
        -26: (("MuonEG+Run"+"F"),"MuonEG"),
        -27: (("MuonEG+Run"+"G"),"MuonEG"),
        -28: (("MuonEG+Run"+"H"),"MuonEG"),
        #
        -31: (("EGamma+Run"+"A"),"EGamma"),
        -32: (("EGamma+Run"+"B"),"EGamma"),
        -33: (("EGamma+Run"+"C"),"EGamma"),
        -34: (("EGamma+Run"+"D"),"EGamma"),
        -35: (("EGamma+Run"+"E"),"EGamma"),
        -36: (("EGamma+Run"+"F"),"EGamma"),
        -37: (("EGamma+Run"+"G"),"EGamma"),
        -38: (("EGamma+Run"+"H"),"EGamma"),
        #
        -41: (("DoubleEG+Run"+"A"),"DoubleEG"),
        -42: (("DoubleEG+Run"+"B"),"DoubleEG"),
        -43: (("DoubleEG+Run"+"C"),"DoubleEG"),
        -44: (("DoubleEG+Run"+"D"),"DoubleEG"),
        -45: (("DoubleEG+Run"+"E"),"DoubleEG"),
        -46: (("DoubleEG+Run"+"F"),"DoubleEG"),
        -47: (("DoubleEG+Run"+"G"),"DoubleEG"),
        -48: (("DoubleEG+Run"+"H"),"DoubleEG"),
        #
        -51: (("SingleElectron+Run"+"A"),"SingleElectron"),
        -52: (("SingleElectron+Run"+"B"),"SingleElectron"),
        -53: (("SingleElectron+Run"+"C"),"SingleElectron"),
        -54: (("SingleElectron+Run"+"D"),"SingleElectron"),
        -55: (("SingleElectron+Run"+"E"),"SingleElectron"),
        -56: (("SingleElectron+Run"+"F"),"SingleElectron"),
        -57: (("SingleElectron+Run"+"G"),"SingleElectron"),
        -58: (("SingleElectron+Run"+"H"),"SingleElectron"),
        #
        -61: (("Tau+Run"+"A"),"Tau"),
        -62: (("Tau+Run"+"B"),"Tau"),
        -63: (("Tau+Run"+"C"),"Tau"),
        -64: (("Tau+Run"+"D"),"Tau"),
        -65: (("Tau"+str(year)+"C*"),"Tau"), # SS
        -66: (("Tau"+str(year)+"D*"),"Tau"), # SS
        #
#        -71: (("SinglePhoton+Run"+"A"),"SinglePhoton"),
#        -72: (("SinglePhoton+Run"+"B"),"SinglePhoton"),
#        -73: (("SinglePhoton+Run"+"C"),"SinglePhoton"),
#        -74: (("SinglePhoton+Run"+"D"),"SinglePhoton"),
#        -75: (("SinglePhoton+Run"+"E"),"SinglePhoton"),
        -76: (("SinglePhoton+Run"+"F"),"SinglePhoton"), # only this used 
#        -77: (("SinglePhoton+Run"+"G"),"SinglePhoton"),
#        -78: (("SinglePhoton+Run"+"H"),"SinglePhoton"),
        #
    }

    print(category)

    if(year == 12016):
        switch = {
            -1: (("SingleMuon+Run"+"B-ver1"),"SingleMuon"),
            -2: (("SingleMuon+Run"+"B-ver2"),"SingleMuon"),
            -3: (("SingleMuon+Run"+"C"),"SingleMuon"),
            -4: (("SingleMuon+Run"+"D"),"SingleMuon"),
            -5: (("SingleMuon+Run"+"E"),"SingleMuon"),
            -6: (("SingleMuon+Run"+"F"),"SingleMuon"),

            -11: (("DoubleMuon+Run"+"B-ver1"),"DoubleMuon"),
            -12: (("DoubleMuon+Run"+"B-ver2"),"DoubleMuon"),
            -13: (("DoubleMuon+Run"+"C"),"DoubleMuon"),
            -14: (("DoubleMuon+Run"+"D"),"DoubleMuon"),
            -15: (("DoubleMuon+Run"+"E"),"DoubleMuon"),
            -16: (("DoubleMuon+Run"+"F"),"DoubleMuon"),

            -21: (("MuonEG+Run"+"B-ver1"),"MuonEG"),
            -22: (("MuonEG+Run"+"B-ver2"),"MuonEG"),
            -23: (("MuonEG+Run"+"C"),"MuonEG"),
            -24: (("MuonEG+Run"+"D"),"MuonEG"),
            -25: (("MuonEG+Run"+"E"),"MuonEG"),
            -26: (("MuonEG+Run"+"F"),"MuonEG"),

            -41: (("DoubleEG+Run"+"B-ver1"),"DoubleEG"),
            -42: (("DoubleEG+Run"+"B-ver2"),"DoubleEG"),
            -43: (("DoubleEG+Run"+"C"),"DoubleEG"),
            -44: (("DoubleEG+Run"+"D"),"DoubleEG"),
            -45: (("DoubleEG+Run"+"E"),"DoubleEG"),
            -46: (("DoubleEG+Run"+"F"),"DoubleEG"),

            -51: (("SingleElectron+Run"+"B-ver1"),"SingleElectron"),
            -52: (("SingleElectron+Run"+"B-ver2"),"SingleElectron"),
            -53: (("SingleElectron+Run"+"C"),"SingleElectron"),
            -54: (("SingleElectron+Run"+"D"),"SingleElectron"),
            -55: (("SingleElectron+Run"+"E"),"SingleElectron"),
            -56: (("SingleElectron+Run"+"F"),"SingleElectron"),

            -81: (("SinglePhoton+Run"+"B-ver1-HIPM"),"SinglePhoton"),
            -82: (("SinglePhoton+Run"+"B-ver2-HIPM"),"SinglePhoton"),
            -83: (("SinglePhoton+Run"+"C-HIPM"),"SinglePhoton"),
            -84: (("SinglePhoton+Run"+"D-HIPM"),"SinglePhoton"),
            -85: (("SinglePhoton+Run"+"E-HIPM"),"SinglePhoton"),
            -86: (("SinglePhoton+Run"+"F-HIPM"),"SinglePhoton"),
        }

    tmp_pair = switch.get(argument, "regex, PDtype")
    finalDir = dirScratch02+"/"+category+"/"+str(year)+"/"+tmp_pair[0]
    if (category=="VBF" and (argument == -62 or argument == -63 or argument == -64) ): finalDir = dirScratch02+"/"+"Zinv"+"/"+str(year)+"/"+tmp_pair[0]
#    if (argument == -65 or argument == -66): finalDir = dirScratchSS+"/"+tmp_pair[0]
    print(finalDir)
    pair = (findDIR(finalDir),tmp_pair[1])
    return pair

def SwitchSample(argument,year):

    #https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
    # cross section from  https://cms-gen-dev.cern.ch/xsdb

    ####----------  relevant disks
    dirT2 = "/mnt/T2_US_MIT/hadoop/cms/store/user/paus/nanohr/D01/"
    dirGluster = "/data/submit/cms/store/user/mariadlf/nano/D01/"
    dirScratch = "/scratch/submit/cms/mariadlf/Hrare/SKIMS/D01/"
    ####----------
    dirLocal = "/work/submit/mariadlf/Hrare/OCT14/"
    dirLocalNEW = "/work/submit/mariadlf/Hrare/D01/2018/"
    dirLocalNEW2 = "/work/submit/mariadlf/Hrare/D02/2018/"
    dirLocalD01 = "/data/submit/mariadlf/HrareSIG/D01/2018/"
    dirLocalTEST = "/work/submit/mariadlf/Hrare/TEST/2018/"
    dirLocalTEST2 = "/work/submit/mariadlf/Hrare/TESTtrackFIX/2018/"

    ####----------
    dirScratch02 = "/scratch/submit/cms/mariadlf/Hrare/newSKIMS/D02/"
    dirGluster02 = "/data/submit/cms/store/user/mariadlf/nano/D02/"
    dirT202 = "/mnt/T2_US_MIT/hadoop/cms/store/user/paus/nanohr/D02/"

    campaign = ""
    if(year == 2018): campaign = "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1*"
    if(year == 2017): campaign = "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9*"
    if(year == 22016): campaign = "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17*"
    if(year == 12016): campaign = "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11*"

    campaignFIX = ""
    if(year == 2018): campaignFIX = "RunIISummer20UL18MiniAODv2-4cores5k_106X_upgrade2018_realistic_v16_L1v1*"
    if(year == 2017): campaignFIX = "RunIISummer20UL17MiniAODv2-4cores5k_106X_mc2017_realistic_v9*"
    if(year == 22016): campaignFIX = "RunIISummer20UL16MiniAODv2-4cores5k_106X_mcRun2_asymptotic_v17*"
    if(year == 12016): campaignFIX = "RunIISummer20UL16MiniAODAPVv2-4cores5k_106X_mcRun2_asymptotic_preVFP_v11*"

    switch = {
        ## SIGNAL
        1010: (dirGluster02+"VBF_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaign,3781.7*0.49), #NNLO
        1011: (dirGluster02+"WplusH_WToLNu_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaign,3*94.26*0.49), #xsec = 3*9.426E-02 (xsec*Wl) * BR(Hphigamma)=1 BR(phi->kk)=0.49
        1012: (dirGluster02+"WminusH_WToLNu_HToRhoGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaign,3*59.83*0.49), #xsec = 3*5.983E-02 (xsecWl) * BR(Hphigamma)=1 BR(phi->kk)=0.49
        1013: (dirGluster02+"ZH_ZToLL_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaign,3*(29.82 - 4.14)*0.49), #xsec = 3*9.426E-02 (xsec*Zll) * BR(Hphigamma)=1 BR(phi->kk)=0.49
        1014: (dirGluster02+"ggZH_ZToLL_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaign,3*4.14*0.49), #xsec = 3*9.426E-02 (xsec*Zll) * BR(Hphigamma)=1 BR(phi->kk)=0.49
        1015: (dirGluster+"ZH_ZToNuNu_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaign,(177.62 - 24.57)*0.49), #xsec = 3*9.426E-02 (xsec*Zll) * BR(Hphigamma)=1 BR(phi->kk)=0.49
        1016: (dirGluster+"ggZH_ZToNuNu_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaign,24.57*0.49), #xsec = 3*9.426E-02 (xsec*Zll) * BR(Hphigamma)=1 BR(phi->kk)=0.49
        1017: (dirGluster02+"GluGlu_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaign,48580*0.49), #xsec = 3*9.426E-02 (xsec*ggH) * BR(Hphigamma)=1 BR(phi->kk)=0.49
        1018: (dirGluster02+"TTHtoPhiG_M-125_TuneCP5_13TeV_powheg-pythia8+"+campaign,505.2*0.49),
        #
        1020: (dirGluster02+"VBF_HToRhoGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaign,3781.7), # xsec = 4pb * BR(Hrhogamma)=1 BR(rho->pipi)=1
        1021: (dirGluster02+"WplusH_WToLNu_HToRhoGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaign,3*94.26), #xsec = 3*9.426E-02 (Wl) * BR(Hrhogamma)=1 BR(rho->pipi)=1
        1022: (dirGluster02+"WminusH_WToLNu_HToRhoGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaign,3*59.83), #xsec = 3*5.983E-02 (Wl) * BR(Hrhogamma)=1 BR(rho->pipi)=1
        1023: (dirGluster02+"ZH_ZToLL_HToRhoGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaign,3*(29.82 - 4.14)), #xsec = 3*9.426E-02 (xsec*Wl) * BR(Hrhogamma)=1 BR(rho->pipi)=1
        1024: (dirGluster02+"ggZH_ZToLL_HToRhoGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaign,3*4.14), #xsec = 3*9.426E-02 (xsec*Wl) * BR(Hrhogamma)=1 BR(rho->pipi)=1
        1025: (dirGluster+"ZH_ZToNuNu_HToRhoGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaign,(177.62 - 24.57)), #xsec = 3*9.426E-02 (xsec*Wl) * BR(Hrhogamma)=1 BR(rho->pipi)=1
        1026: (dirGluster+"ggZH_ZToNuNu_HToRhoGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaign,24.57), #xsec = 3*9.426E-02 (xsec*Wl) * BR(Hrhogamma)=1 BR(rho->pipi)=1
        1027: (dirGluster02+"GluGlu_HToRhoGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaign,48580), #xsec = 3*9.426E-02 (xsec*ggH) * BR(Hrhogamma)=1 BR(rho->pipi)=1
        1028: (dirGluster02+"TTHtoRhoG_M-125_TuneCP5_13TeV_powheg-pythia8+"+campaign,505.2),
        #
        1030: (dirGluster02+"VBF_HtoK0starG_M-125_TuneCP5_PSWeights_13TeV_powheg-pythia8+"+campaign,3781.7),
        1031: (dirGluster02+"WplusH_WToLNu_HtoK0starG_M-125_TuneCP5_PSWeights_13TeV_powheg-pythia8+"+campaign,3*94.26),
        1032: (dirGluster02+"WminusH_WtoLNu_HtoK0starG_M-125_TuneCP5_PSWeights_13TeV_powheg-pythia8+"+campaign,3*59.83),
        1033: (dirGluster02+"ZH_ZToLL_HtoK0starG_M-125_TuneCP5_PSWeights_13TeV_powheg-pythia8+"+campaign,3*(29.82 - 4.14)),
        1034: (dirGluster02+"GluGluZH_ZToLL_HtoK0starG_M-125_TuneCP5_PSWeights_13TeV_powheg-pythia8+"+campaign,3*4.14),
        1035: (dirGluster+"ZH_ZToNuNu_HtoK0starG_M-125_TuneCP5_PSWeights_13TeV_powheg-pythia8+"+campaign,(177.62 - 24.57)),
        1036: (dirGluster+"ggZH_ZToNuNu_HtoK0starG_M-125_TuneCP5_PSWeights_13TeV_powheg-pythia8+"+campaign,24.57),
        1037: (dirGluster02+"GluGluHtoK0starG_M-125_TuneCP5_PSWeights_13TeV_powheg-pythia8+"+campaign,48580),
        1038: (dirGluster02+"TTH_HtoK0starG_M-125_TuneCP5_PSWeights_13TeV_powheg-pythia8+"+campaign,505.2),
        #
        1039: (dirLocalNEW2+"ggh-hK0Stargamma-powheg_2"+"/NANOAOD_02",48580), #xsec = 3*9.426E-02 (xsec*ggH) * BR(HkstaraGamma)=1
        #
        1040: (dirLocalNEW2+"ggh-homegagamma-powheg"+"/NANOAOD_03_test5",48580*0.892),
        1041: (dirLocalNEW2+"ggh-hphipipipi0gamma-powheg"+"/NANOAOD_03_test5",48580*0.153),
        1042: (dirLocalNEW2+"ggh-hD0StarKmPiPPi0gamma-powheg"+"/NANOAOD_03_test5",48580*0.14),
        1043: (dirLocalNEW2+"ggh-hD0Stargamma-powheg"+"/NANOAOD_03_test5",48580*0.0389),
        #
        1054: (dirLocalNEW+"vbf-hphiKLKSgamma-powheg"+"/NANOAOD_01",3781.7*0.24), # xsec = 4pb * BR(Hphigamma)=1 BR(phi->kLkS)=0.24
        #
        1053: (dirT2+"ZH_HToJPsiG_JPsiToMuMu_TuneCP5_13TeV-madgraph-pythia8+"+campaign,6067*1000), #check xSEC
        # TESTS
        1050: (dirLocalTEST+"ggh-hrhogamma-powheg"+"/NANOAOD_01",46870), #xsec = 3*9.426E-02 (xsec*ggH) * BR(Hrhogamma)=1 BR(rho->pipi)=1
        1051: (dirLocalTEST2+"ggh-hrhogamma-powheg"+"/NANOAOD_01",48580), #xsec = 3*9.426E-02 (xsec*ggH) * BR(Hrhogamma)=1 BR(rho->pipi)=1
        1052: (dirLocalTEST2+"ggh-hphigamma-powheg"+"/NANOAOD_01",48580*0.49), #xsec = 3*9.426E-02 (xsec*ggH) * BR(Hrhogamma)=1 BR(rho->pipi)=1
        216: (dirLocal+"ZLLphigamma_pythia8_genFix",0.10*1000), #xsec=1pb * 0.101 (Zll) * BR(Hphigamma)=1
#        10: (dirLocalNEW+"zh-hphigamma-powheg",0.10*1000), #xsec=1pb * 0.101 (Zll) * BR(Hphigamma)=1
###### new set
#        13: (dirLocalNEW+"zh_Phigamma",0.10*1000), #xsec=1pb * 0.101 (Zll) * BR(Hphigamma)=1
#        14: (dirLocalNEW+"zh_Rhogamma",0.10*1000), #xsec=1pb * 0.101 (Zll) * BR(Hphigamma)=1
#        15: (dirLocalNEW+"zh_Omegagamma",0.10*1000), #xsec=1pb * 0.101 (Zll) * BR(Hphigamma)=1
#        16: (dirLocalNEW+"zh_KLKSgamma",0.10*1000), #xsec=1pb * 0.101 (Zll) * BR(Hphigamma)=1
#        1000: (dirT2+"DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaign,6067*1000), #NNLO
        0: (dirGluster02+"DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaign,6067*1000), #NNLO (LO was 5398.0)
#        0: (dirT2+"DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaign,6067*1000), #NNLO (LO was 5398.0)
        1: (dirGluster02+"ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign, 51.1*1000), #LO
        2: (dirGluster02+"WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign, 191.3*1000), #LO
        #        2: (dirT2+"WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM", 412.7*1000), #LO
        3: (dirGluster02+"WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaign,61526.7*1000), #NNLO (LO was 53870.0)
        4: (dirT2+"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8+"+campaign,88.2*1000), #NNLO
        5: (dirT2+"TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8+"+campaign,365.3452*1000), #NNLO
        ##
        9: (dirGluster02+"VBFGamma_5f_TuneCP5_DipoleRecoil_13TeV-madgraph-pythia8+"+campaign,21.09*1000), #LO
        ##
        10: (dirGluster02+"GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaign,18540.0*1000*1.26), #LO *1.26
        11: (dirGluster02+"GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaignFIX,8644.0*1000*1.26), #LO *1.26
        12: (dirGluster02+"GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaign,2183.0*1000*1.26), #LO *1.26
        13: (dirGluster02+"GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaign,260.2*1000*1.26), #LO *1.26
        14: (dirGluster02+"GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaign,86.58*1000*1.26), #LO *1.26
        #
        15: (dirGluster02+"GJets_DR-0p4_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaign,15750.0*1000*1.26), #LO *1.26
        16: (dirGluster02+"GJets_DR-0p4_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaign,5034*1000*1.26), #LO *1.26
        17: (dirGluster02+"GJets_DR-0p4_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaign,1129*1000*1.26), #LO *1.26
        18: (dirGluster02+"GJets_DR-0p4_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaign,126.2*1000*1.26), #LO *1.26
        19: (dirGluster02+"GJets_DR-0p4_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaign,41.31*1000*1.26), #LO *1.26
        ####
        20: (dirT2+"QCD_Pt-30to50_EMEnriched_TuneCP5_13TeV-pythia8+"+campaign,6447000.0*1000*1.26), #LO *1.26
        21: (dirT2+"QCD_Pt-50to80_EMEnriched_TuneCP5_13TeV-pythia8+"+campaign,1988000.0*1000*1.26), #LO *1.26
        22: (dirT2+"QCD_Pt-80to120_EMEnriched_TuneCP5_13TeV-pythia8+"+campaign,367500.0*1000*1.26), #LO *1.26
        23: (dirT2+"QCD_Pt-120to170_EMEnriched_TuneCP5_13TeV-pythia8+"+campaign,66590.0*1000*1.26), #LO *1.26
        24: (dirT2+"QCD_Pt-170to300_EMEnriched_TuneCP5_13TeV-pythia8+"+campaign,16620.0*1000*1.26), #LO *1.26
        25: (dirT2+"QCD_Pt-300toInf_EMEnriched_TuneCP5_13TeV-pythia8+"+campaign,1104.0*1000*1.26), #LO *1.26

        31: (dirGluster02+"WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,53330.0*1000), #LO
        32: (dirGluster02+"WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,8875.0*1000), #LO
        33: (dirGluster02+"WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,3338.0*1000), #LO

        34: (dirGluster02+"DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,5129.0*1000), #LO
        35: (dirGluster02+"DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,951.5*1000), #LO
        36: (dirGluster02+"DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,361.4*1000), #LO

        37: (dirT2+"Z1JetsToNuNu_M-50_LHEFilterPtZ-50To150_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,579.9*1000),
        38: (dirT2+"Z1JetsToNuNu_M-50_LHEFilterPtZ-150To250_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,17.42*1000),
        39: (dirT2+"Z1JetsToNuNu_M-50_LHEFilterPtZ-250To400_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,1.987*1000),
        40: (dirT2+"Z1JetsToNuNu_M-50_LHEFilterPtZ-400ToInf_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,0.2179*1000),
        41: (dirT2+"Z2JetsToNuNu_M-50_LHEFilterPtZ-50To150_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,312.5*1000),
        42: (dirT2+"Z2JetsToNuNu_M-50_LHEFilterPtZ-150To250_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,30*1000),
        43: (dirT2+"Z2JetsToNuNu_M-50_LHEFilterPtZ-250To400_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,4.98*1000),
        44: (dirT2+"Z2JetsToNuNu_M-50_LHEFilterPtZ-400ToInf_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,0.8165*1000),

        45: (dirT2+"WGammaToJJGamma_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,8.635*1000),
        46: (dirT2+"ZGammaToJJGamma_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,4.144*1000),
        47: (dirGluster02+"TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8+"+campaign,3.757*1000),
        48: (dirT2+"ZGTo2NuG_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,30.11*1000),

        ##
        49: (dirT2+"WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8+"+campaign,11.09*1000),
        #50: (dirT2+"WWTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,51.65*1000),
        51: (dirT2+"WZTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,9.119*1000),
        #52: (dirT2+"WZTo1L3Nu_4f_TuneCP5_13TeV-amcatnloFXFX-pythia+"+campaign,3.414*1000),
        #53: (dirT2+"WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,5.213*1000),
        #54: (dirT2+"WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,6.419*1000),
        #55: (dirT2+"ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8+"+campaign,0.9738*1000),
        #56: (dirT2+"ZZTo2Nu2Q_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,4.545*1000),
        #57: (dirT2+"ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8+"+campaign,3.676*1000),
        #58: (dirT2+"ZZTo4L_TuneCP5_13TeV_powheg_pythia8+"+campaign,1.325*1000),
        #WWTo4Q_4f
        #ZZTo4Q_5f
    }

    return switch.get(argument, "BKGdefault, xsecDefault")

def pickTRG(overall,year,PDType,isVBF,isW,isZ,isZinv):

    TRIGGER=''
    if(year == 2018 and isVBF):
        if (PDType== "EGamma"): TRIGGER=getTriggerFromJson(overall, "isVBF", year)
        elif (PDType== "Tau"): TRIGGER="{0} and not {1}".format(getTriggerFromJson(overall, "isZinv", year),getTriggerFromJson(overall, "isVBF", year))
        elif (PDType== "NULL"): TRIGGER=getTriggerFromJson(overall, "isVBFor", year)   # MC seems the same
    if(year == 2018 and isZinv):
        if (PDType== "Tau"): TRIGGER=getTriggerFromJson(overall, "isZinv", year)
        elif (PDType== "NULL"): TRIGGER=getTriggerFromJson(overall, "isZinv", year)    # MC seems the same
    if(year == 2017 and isVBF):
        if (PDType== "SinglePhoton"): TRIGGER=getTriggerFromJson(overall, "isVBF", year)
        elif (PDType== "NULL"): TRIGGER=getTriggerFromJson(overall, "isVBF", year)     # MC seems the same
    if((year == 12016 or year == 22016) and isVBF):
        if (PDType== "SinglePhoton"): TRIGGER=getTriggerFromJson(overall, "isVBF", year)
        elif (PDType== "NULL"): TRIGGER=getTriggerFromJson(overall, "isVBF", year)     # MC seems the same

    if(year == 2018 and (isW or isZ)):
        if (PDType == "DoubleMuon"): TRIGGER = "{0}".format(getTriggerFromJson(overall, "diMU", year))
        elif (PDType == "SingleMuon"): TRIGGER = "{0} and not {1}".format(getTriggerFromJson(overall, "oneMU", year),getTriggerFromJson(overall, "diMU", year))
        elif (PDType == "EGamma"): TRIGGER = "({0} or {1}) and not {2} and not {3}".format(getTriggerFromJson(overall, "oneEL", year),getTriggerFromJson(overall, "diEL", year),getTriggerFromJson(overall, "oneMU", year),getTriggerFromJson(overall, "diMU", year))
        elif(PDType == "MuonEG"): TRIGGER = "{0} and not {1} and not {2} and not {3}".format(getTriggerFromJson(overall, "muG", year),getTriggerFromJson(overall, "oneEL", year),getTriggerFromJson(overall, "diEL", year),getTriggerFromJson(overall, "oneMU", year),getTriggerFromJson(overall, "diMU", year))
        elif(year == 2018):
            TRIGGER = "{0} or {1} or {2} or {3} or {4}".format(getTriggerFromJson(overall, "oneEL", year),getTriggerFromJson(overall, "diEL", year),getTriggerFromJson(overall, "oneMU", year),getTriggerFromJson(overall, "diMU", year),getTriggerFromJson(overall, "muG", year))
        else:
            print("PROBLEM with triggers!!!")

    if((year == 2017 or year == 22016 or year == 12016) and (isW or isZ)):
        if (PDType == "DoubleMuon"): TRIGGER = "{0}".format(getTriggerFromJson(overall, "diMU", year))
        elif (PDType == "SingleMuon"): TRIGGER = "{0} and not {1}".format(getTriggerFromJson(overall, "oneMU", year),getTriggerFromJson(overall, "diMU", year))
        elif (PDType == "DoubleEG"): TRIGGER = "{0} and not {1} and not {2}".format(getTriggerFromJson(overall, "diEL", year),getTriggerFromJson(overall, "oneMU", year),getTriggerFromJson(overall, "diMU", year))
        elif (PDType == "SingleElectron"): TRIGGER = "{0} and not {1} and not {2} and not {3}".format(getTriggerFromJson(overall, "oneEL", year),getTriggerFromJson(overall, "diEL", year),getTriggerFromJson(overall, "oneMU", year),getTriggerFromJson(overall, "diMU", year))
        elif(PDType == "MuonEG"): TRIGGER = "{0} and not {1} and not {2} and not {3}".format(getTriggerFromJson(overall, "muG", year),getTriggerFromJson(overall, "oneEL", year),getTriggerFromJson(overall, "diEL", year),getTriggerFromJson(overall, "oneMU", year),getTriggerFromJson(overall, "diMU", year))
        elif(year == 2017 or year == 22016 or year == 12016):
            TRIGGER = "{0} or {1} or {2} or {3} or {4}".format(getTriggerFromJson(overall, "oneEL", year),getTriggerFromJson(overall, "diEL", year),getTriggerFromJson(overall, "oneMU", year),getTriggerFromJson(overall, "diMU", year),getTriggerFromJson(overall, "muG", year))
        else:
            print("PROBLEM with triggers!!!")

    return TRIGGER

def computeWeigths(df, files, sampleNOW, year, isMC):

    nevents = df.Count().GetValue()
    print("%s entries in the dataset" %nevents)

    if not isMC:
        return 1.
    else:
        rdf = ROOT.RDataFrame("Runs", files)
        genEventSumWeight = rdf.Sum("genEventSumw").GetValue()
        genEventSumNoWeight = rdf.Sum("genEventCount").GetValue()

        weight = (SwitchSample(sampleNOW,year)[1] / genEventSumWeight)

        lumiEq = (genEventSumNoWeight / SwitchSample(sampleNOW,year)[1])
        print("lumi equivalent fb %s" %lumiEq)

        return weight

