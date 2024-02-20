import ROOT
import glob
from subprocess import check_output
import sys

def loadUserCode():
    print('loadUserCode()')
    ROOT.gSystem.AddDynamicPath("./.")
    ROOT.gROOT.ProcessLine(".include ./.")
    ROOT.gInterpreter.AddIncludePath("./.")
    ROOT.gInterpreter.Declare('#include "config/functions.h"')

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
    ROOT.gInterpreter.Declare('#include "config/sfCorrLib.h"')
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

def listDir(isT2,xrdpath):
    print(xrdpath)
    if isT2:
        xrd = "root://xrootd.cmsaf.mit.edu/"
    else:
        xrd = "root://submit50.mit.edu/"
    f = check_output(['xrdfs', f'{xrd}', 'ls', xrdpath]).decode(sys.stdout.encoding)

    rootFiles = ROOT.vector('string')()
    stringFiles = f.split()

    for e in range(len(stringFiles)):
        filePath = xrd + stringFiles[e]
        rootFiles.push_back(filePath)

    print(len(rootFiles))
    return rootFiles


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

