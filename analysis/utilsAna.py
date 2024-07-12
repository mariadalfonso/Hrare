import ROOT
import glob
from subprocess import check_output
import sys

def loadUserCode():
    print('loadUserCode()')
    ROOT.gSystem.AddDynamicPath("./.")
    ROOT.gROOT.ProcessLine(".include ./.")
    ROOT.gInterpreter.AddIncludePath("./.")
#    ROOT.gInterpreter.Declare('#include "./config/functions.h"')
    ROOT.gInterpreter.ProcessLine('#include "./config/functions.h"')

def loadSFhisto(mc,year):
    ROOT.gInterpreter.ProcessLine('initTrigSF();')
    ROOT.gInterpreter.ProcessLine('initIsoSF();')

def loadPolarizationTree(mc,year):
    ROOT.gInterpreter.ProcessLine('initPol(%d,%d,%d);'% (mc,year,ROOT.GetThreadPoolSize()))

def loadtmva_helper():
    ROOT.gInterpreter.Declare('#include "tmva_helper_xml.h"')

def loadCorrectionSet(year):
    print('loadCorrectionSet()')
    import correctionlib
    correctionlib.register_pyroot_binding()
    ROOT.gInterpreter.Declare('#include "./config/sfCorrLib.h"')
#    ROOT.gInterpreter.Declare('#include "config/mysf.h"')
#    ROOT.gInterpreter.Load("config/mysf.so")
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
        '''
    )

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

def readDataQuality(year):
    dirJson = "./config"
    if(year == 2018):
        loadJSON("{}/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt".format(dirJson))
    if(year == 2017):
        loadJSON("{}/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt".format(dirJson))
    if(year == 22016 or year == 12016):
        loadJSON("{}/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt".format(dirJson))

## I have two here
def findDIR(directory,useXROOTD=False):

    print('HELLO')
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

def listDir(isT2,useXcache,xrdpath):
    print(xrdpath)
    if isT2:
        if useXcache: xrd = "root://bost-cms-xcache01.lhcone.es.net:1094/"
        else: xrd = "root://xrootd.cmsaf.mit.edu/"
    else:
        #for /data/submit/
        xrd = "root://submit50.mit.edu/"
    f = check_output(['xrdfs', f'{xrd}', 'ls', xrdpath]).decode(sys.stdout.encoding)

    rootFiles = ROOT.vector('string')()
    stringFiles = f.split()

    for e in range(len(stringFiles)):
        filePath = xrd + stringFiles[e]
        rootFiles.push_back(filePath)

    print(len(rootFiles))
    return rootFiles

def BuildDictJpsiCC():

    # cross section from  https://cms-gen-dev.cern.ch/xsdb  (normally in /fb)
    useXcache=False
    isT2=False
    dirNameSig="/data/submit/cms/store/user/mariadlf/nano/GluGluH_HJPsiCC/NANOAOD_test3"
    dirNameData="/scratch/submit/cms/mariadlf/Hrare/newSKIMS/D04/BPH/2018/"
    dirNameBKG="/data/submit/kyoon/D04/"
    dirNameBKG2="/data/submit/mariadlf/D04/"
    campaignv1="RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM"
    campaignv2="RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM"

    thisdict = {
        -1: (findDIR(dirNameData+"Charmonium+RunA"),-1),
        -2: (findDIR(dirNameData+"Charmonium+RunB"),-1),
        -3: (findDIR(dirNameData+"Charmonium+RunC"),-1),
        -4: (findDIR(dirNameData+"Charmonium+RunD"),-1),
        10: (findDIR(dirNameBKG+"JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8+"+campaignv2),868900.0*1000*0.0593), # BR for Jpsi to MuMu 5.93% (PHYTHIA was forced all with JpsiMuMu) + from cards filterEff = 0.26 and total xsection (30560000.0)
        11: (findDIR(dirNameBKG+"BToJpsi_JPsiToMuMu_BMuonFilter_HardQCD_TuneCP5_13TeV-pythia8-evtgen+"+campaignv1),146300.0*1000), # EvtJen has exclusive processes
#        12: (findDIR(dirNameBKG2+"BuToJpsiK_BMuonFilter_TuneCP5_13TeV-pythia8-evtgen+"+campaignv2),21930000.0*1000),
#        12: (findDIR(dirNameBKG2+"BcToJPsiMuMu_inclusive_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen+"+campaignv1),4207000.0*1000),        # there are 3 datasets with extentions
#        13: (findDIR(dirNameBKG2+"LambdaBToJpsiLambda_JpsiToMuMu_TuneCP5_13TeV-pythia8-evtgen+"+campaignv1),2148000.0*1000),
        1000: (findDIR(dirNameSig),1.*0.0593)
    }

    #BToJpsi
    #https://github.com/cms-data/GeneratorInterface-EvtGenInterface/blob/master/incl_BtoJpsi_mumu.dec
    #BPH-RunIISummer20UL18GEN-00015-fragment_JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8.py
    #BPH-RunIISummer20UL18GEN-00195-fragment_BToJpsi_JPsiToMuMu_BMuonFilter_HardQCD_TuneCP5_13TeV-pythia8-evtgen.py

    return thisdict

def BuildDict():

    #useXROOTD=False
    useXcache=False
    isT2=True
    if isT2:
        dirName="/store/user/paus/nanohr/D02/"
    else:
        # at CERN
        dirName="/store/user/mariadlf/nano/D02/"

    campaignv3 = "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3+MINIAODSIM/"
    campaignv1 = "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM/"
    campaignv2 = "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM/"
    campaignFIX = "RunIISummer20UL18MiniAODv2-4cores5k_106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM/"

    thisdict = {
        1010: (listDir(isT2,useXcache,dirName+"VBF_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaignv1),3781.7*0.49),
        1017: (listDir(isT2,useXcache,dirName+"GluGlu_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaignv3),48580*0.49),
        10: (listDir(isT2,useXcache,dirName+"GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaignv2),18540.0*1000*1.26), #LO *1.26
        11: (listDir(isT2,useXcache,dirName+"GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaignFIX),8644.0*1000*1.26), #LO *1.26
        12: (listDir(isT2,useXcache,dirName+"GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaignv2),2183.0*1000*1.26), #LO *1.26
        13: (listDir(isT2,useXcache,dirName+"GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaignv2),260.2*1000*1.26), #LO *1.26
        14: (listDir(isT2,useXcache,dirName+"GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaignv2),86.58*1000*1.26) #LO *1.26
    }

    return thisdict

def SwitchSample(thisdict,argument):

    return thisdict.get(argument, "BKGdefault, xsecDefault")

def pickTRG(overall,year,PDType,isVBF,isW,isZ,isZinv,isBPH):

    TRIGGER=''
    if(year == 2018 and isVBF):
        if (PDType== "EGamma"): TRIGGER=getTriggerFromJson(overall, "isVBF", year)
        elif (PDType== "Tau"): TRIGGER="{0} and not {1}".format(getTriggerFromJson(overall, "isZinv", year),getTriggerFromJson(overall, "isVBF", year))
        elif (PDType== "NULL"): TRIGGER=getTriggerFromJson(overall, "isVBFor", year)   # MC seems the same

    if(year == 2018 and isZinv):
        if (PDType== "Tau"): TRIGGER=getTriggerFromJson(overall, "isZinv", year)
        elif (PDType== "NULL"): TRIGGER=getTriggerFromJson(overall, "isZinv", year)    # MC seems the same

    if(year == 2018 and isBPH):
        if (PDType== "Charmonium"): TRIGGER=getTriggerFromJson(overall, "isBPH", year)
        elif (PDType== "NULL"): TRIGGER=getTriggerFromJson(overall, "isBPH", year)    # MC seems the same

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

def getTriggerFromJson(overall, type, year ):

    for trigger in overall:
        if trigger['name'] == type and trigger['year'] == year: return trigger['definition']

def getMesonFromJson(overall, type, cat ):

    for meson in overall:
        if meson['name'] == type and meson['type'] == cat: return meson['definition']

def getMVAFromJson(overall, type, cat ):

    for meson in overall:
        if meson['name'] == type and meson['type'] == cat: return meson['file']


def getSkimsFromJson(overall, type ):

    for selection in overall:
        if selection['type'] == type : return selection['PRESELECTION']
