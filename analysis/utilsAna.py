import ROOT
import glob
from subprocess import check_output
import sys
import json
import os

lumisMgamma={
    '12016': 19.52, #APV #(B-F for 2016 pre)
    '22016': 16.80, #postVFP
    '2016': 35.9,
    '2017': 41.5,
    '12017': 7.7, #(F for 2017)
    '2018': 59.70,
    '12018': 39.54,
    'all': 86.92,      #19.52 + 7.7 + 59.70
    '12022':7.98, # C-D
    '22022':26.67, # E, F, G
    '12023':17.794, #C
    '22023':9.451, #D
}

lumisJpsiCC={
    '12016': 19.52, #APV #(B-F for 2016 pre)
    '22016': 16.80, #postVFP
    '2017': 41.5,
    '2018': 59.70,
    '12022':7.98, # C-D
    '22022':26.67, # E, F, G
    '12023':17.794, #C
    '22023':9.451, #D
}

xsecRun2={
    'ggH':48580,
    'VBFH':3781.7,
    'Z':61216030, #(6067*1000*10.09) ## 6067 is the NNLO DYJetsToLL_M-50 xsec 10.09 is to correct the BR ToLL,
}

#https://arxiv.org/pdf/2402.09955
xsecRun3={
    'ggH':52230, # 0.4 × σ(13 TeV) + 0.6 × σ(14 TeV)
    'VBFH':4078,
    'Z':61216030, # same as Run2 for now
}

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

def loadtmvahelper():
    print('loadtmvahelper()')
    ROOT.gInterpreter.ProcessLine('#include "./config/tmva_helper_xml.h"')

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

#    print('JSON=',fIn)

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
    print('HELLO readDataQuality', year)
    dirJson = "./config"
    if(str(year) == '2018'):
        loadJSON("{}/cert/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt".format(dirJson))
    if(str(year) == '2017'):
        loadJSON("{}/cert/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt".format(dirJson))
    if(str(year) == '22016' or year == '12016'):
        loadJSON("{}/cert/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt".format(dirJson))
    if(str(year) == '12022'):
        loadJSON("{}/cert/Cert_Collisions2022_355100_357900_eraBCD_Golden.json".format(dirJson))
    if(str(year) == '22022'):
        loadJSON("{}/cert/Cert_Collisions2022_359022_362760_eraEFG_Golden.json".format(dirJson))
    if(str(year) == '12023'):
        loadJSON("{}/cert/Cert_Collisions2023_366403_369802_eraBC_Golden.json".format(dirJson))
    if(str(year) == '22023'):
        loadJSON("{}/cert/Cert_Collisions2023_369803_370790_eraD_Golden.json".format(dirJson))

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

    xcache = "root://cms-xcache.rcac.purdue.edu:1094/";
#    xcache = "root://bost-cms-xcache01.lhcone.es.net:1094/";

    rootFiles = ROOT.vector('string')()
    print(xrdpath)
    if isT2:
        xrd = "root://xrootd.cmsaf.mit.edu/"
    else:
        xrd = "root://submit50.mit.edu/"
    f = check_output(['xrdfs', f'{xrd}', 'ls', xrdpath]).decode(sys.stdout.encoding)

    stringFiles = f.split()
    for e in range(len(stringFiles)):
        filePath = xrd + stringFiles[e]
        if ".root" not in filePath: continue
        filePath = filePath.replace("root://xrootd.cmsaf.mit.edu/",xcache)
        rootFiles.push_back(filePath)
    print(len(rootFiles))
    return rootFiles

def BuildDictJpsiCC(year):

    # cross section from  https://cms-gen-dev.cern.ch/xsdb  (normally in /fb)
    useXcache=False
    isT2=False
    dirNameSig="/ceph/submit/data/group/cms/store/user/mariadlf/nano/"
    #    dirNameSig="/data/submit/cms/store/user/mariadlf/nano/GluGluH_HJPsiCC/NANOAOD_D04"
    dirNameData="/scratch/submit/cms/mariadlf/Hrare/newSKIMS/D04/BPH/"
    dirNameData2="/scratch/submit/cms/mariadlf/Hrare/newSKIMS/D04-MoreTRG/BPH/"
    dirNameBKG="/ceph/submit/data/user/k/kyoon/D04/"
    dirNameBKG2="/data/submit/mariadlf/D04/"
    dirNameBKG3="/ceph/submit/data/user/m/mariadlf/Hrare_psiCC/D05/"
#    campaignv1="RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM"
#    campaignv2="RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM"
    campaign = ""
    if(year == '2018'): campaign = "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1*"
    if(year == '2017'): campaign = "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9*"
    if(year == '22016'): campaign = "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17*"
    if(year == '12016'): campaign = "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11*"
    if(year == '12022'): campaign = "Run3Summer22MiniAODv4-validDigi_130X_mcRun3_2022_realistic_v5*"

    thisdict = {
        12: (findDIR(dirNameBKG3+"InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8+"+campaign),55960000000.0*1000),
        10: (findDIR(dirNameBKG+"JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8+"+campaign),868900.0*1000*0.0593), # BR for Jpsi to MuMu 5.93% (PHYTHIA was forced all with JpsiMuMu) + from cards filterEff = 0.26 and total xsection (30560000.0)
        11: (findDIR(dirNameBKG+"BToJpsi_JPsiToMuMu_BMuonFilter_HardQCD_TuneCP5_13TeV-pythia8-evtgen+"+campaign),146300.0*1000), # EvtJen has exclusive processes
#        12: (findDIR(dirNameBKG2+"BuToJpsiK_BMuonFilter_TuneCP5_13TeV-pythia8-evtgen+"+campaignv2),21930000.0*1000),
#        12: (findDIR(dirNameBKG2+"BcToJPsiMuMu_inclusive_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen+"+campaignv1),4207000.0*1000),        # there are 3 datasets with extentions
#        13: (findDIR(dirNameBKG2+"LambdaBToJpsiLambda_JpsiToMuMu_TuneCP5_13TeV-pythia8-evtgen+"+campaignv1),2148000.0*1000),
        1000: (findDIR(dirNameSig+"GluGluH_HJPsiCC/NANOAOD_D04"),xsecRun2['ggH']*2*0.00001*0.0593), # 2 x 10^-5 is the BR from the theorists
        1001: (findDIR(dirNameSig+"ZJets_ZJPsiCC/NANOAOD_D04"),xsecRun2['Z']*7.3*0.00001*0.0593) # 0.0593 is the BR for rhe Jpsi->mumu, 7.3 x 10^-5 is the BR
    }

    if(year == '12016'):
        dict_ = {
            -1: (findDIR(dirNameData2+str(year)+"/Charmonium+RunB-ver1"),-1),
            -2: (findDIR(dirNameData2+str(year)+"/Charmonium+RunB-ver2"),-1),
            -3: (findDIR(dirNameData2+str(year)+"/Charmonium+RunC"),-1),
            -4: (findDIR(dirNameData2+str(year)+"/Charmonium+RunD"),-1),
            -5: (findDIR(dirNameData2+str(year)+"/Charmonium+RunE"),-1),
            -6: (findDIR(dirNameData2+str(year)+"/Charmonium+RunF"),-1),
        }
        thisdict.update(dict_)

    if(year == '22016'):
        dict_ = {
            -6: (findDIR(dirNameData2+str(year)+"/Charmonium+RunF"),-1),
            -7: (findDIR(dirNameData2+str(year)+"/Charmonium+RunG"),-1),
            -8: (findDIR(dirNameData2+str(year)+"/Charmonium+RunH"),-1),
            }
        thisdict.update(dict_)

    if(year == '2018' or year == '2017'):
        dict_ = {
            -1: (findDIR(dirNameData2+str(year)+"/Charmonium+RunA"),-1),
            -2: (findDIR(dirNameData2+str(year)+"/Charmonium+RunB"),-1),
            -3: (findDIR(dirNameData2+str(year)+"/Charmonium+RunC"),-1),
            -4: (findDIR(dirNameData2+str(year)+"/Charmonium+RunD"),-1),
            -5: (findDIR(dirNameData2+str(year)+"/Charmonium+RunE"),-1),
            -6: (findDIR(dirNameData2+str(year)+"/Charmonium+RunF"),-1),
        }
        thisdict.update(dict_)

    if(year == '12022'):
        dict_ = {
            -13: (findDIR(dirNameData2+str(year)+"/ParkingDoubleMuonLowMass*+RunC"),-1),
            -14: (findDIR(dirNameData2+str(year)+"/ParkingDoubleMuonLowMass*+RunD"),-1),
        }
        thisdict.update(dict_)

    if(year == '22022'):
        dict_ = {
            -15: (findDIR(dirNameData2+str(year)+"/ParkingDoubleMuonLowMass*+RunE"),-1),
            -16: (findDIR(dirNameData2+str(year)+"/ParkingDoubleMuonLowMass*+RunF"),-1),
            -17: (findDIR(dirNameData2+str(year)+"/ParkingDoubleMuonLowMass*+RunG"),-1),
        }
        thisdict.update(dict_)

    if(year == '12023'):
        dict_ = {
            -21: (findDIR(dirNameData2+str(year)+"/ParkingDoubleMuonLowMass*v1+RunC"),-1),
            -22: (findDIR(dirNameData2+str(year)+"/ParkingDoubleMuonLowMass*v2+RunC"),-1),
            -23: (findDIR(dirNameData2+str(year)+"/ParkingDoubleMuonLowMass*v3+RunC"),-1),
            -24: (findDIR(dirNameData2+str(year)+"/ParkingDoubleMuonLowMass*v4+RunC"),-1),
        }
        thisdict.update(dict_)

    if(year == '22023'):
        dict_ = {
            -31: (findDIR(dirNameData2+str(year)+"/ParkingDoubleMuonLowMass*v1+RunD"),-1),
            -32: (findDIR(dirNameData2+str(year)+"/ParkingDoubleMuonLowMass*v2+RunD"),-1),
        }
        thisdict.update(dict_)

    #BToJpsi
    #https://github.com/cms-data/GeneratorInterface-EvtGenInterface/blob/master/incl_BtoJpsi_mumu.dec
    #BPH-RunIISummer20UL18GEN-00015-fragment_JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8.py
    #BPH-RunIISummer20UL18GEN-00195-fragment_BToJpsi_JPsiToMuMu_BMuonFilter_HardQCD_TuneCP5_13TeV-pythia8-evtgen.py

    return thisdict

def BuildDictMgamma():

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
        1010: (listDir(isT2,dirName+"VBF_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaignv1),xsecRun2['VBFH']*0.49),
        1017: (listDir(isT2,dirName+"GluGlu_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaignv3),xsecRun2['ggH']*0.49),
        10: (listDir(isT2,dirName+"GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaignv2),18540.0*1000*1.26), #LO *1.26
        11: (listDir(isT2,dirName+"GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaignFIX),8644.0*1000*1.26), #LO *1.26
        12: (listDir(isT2,dirName+"GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaignv2),2183.0*1000*1.26), #LO *1.26
        13: (listDir(isT2,dirName+"GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaignv2),260.2*1000*1.26), #LO *1.26
        14: (listDir(isT2,dirName+"GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaignv2),86.58*1000*1.26), #LO *1.26
        -62: (listDir(isT2,dirName+"Tau+Run2018B-UL2018_MiniAODv2-v2+MINIAOD"),-1),
        -63: (listDir(isT2,dirName+"Tau+Run2018C-UL2018_MiniAODv2-v1+MINIAOD"),-1),
        -64: (listDir(isT2,dirName+"Tau+Run2018D-UL2018_MiniAODv2-v1+MINIAOD"),-1)
    }

    return thisdict

def BuildDictMgammaRun3(year):

    # Run3 xsection from https://xsecdb-xsdb-official.app.cern.ch/xsdb/
    dirCeph = '/ceph/submit/data/user/m/mariadlf/Hrare_gammaM/D05/2024/'
    dirCephGroup = '/ceph/submit/data/group/cms/store/user/mariadlf/D05/'
    dirNameDataSkims = '/scratch/submit/cms/mariadlf/Hrare/newSKIMS/D05/Zinv/'

    isT2=True
    dirName = "/store/user/paus/nanohr/D05/"
    campaign = ""
    #if(year == '2018'): campaign = "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1*"
    #if(year == '2017'): campaign = "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9*"
    #if(year == '22016'): campaign = "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17*"
    #if(year == '12016'): campaign = "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11*"

    if(str(year) == '12022'): campaign = "Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_*"
    if(str(year) == '22022'): campaign = "Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_*"
    if(str(year) == '12023'): campaign = "Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14*"
    if(str(year) == '22023'): campaign = "Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix_v*"

    print(campaign,'campaign')
    sampleName = 'G-4Jets'
    if (str(year) == '12023' or str(year) == '22023'): sampleName = 'GJ-4Jets'

    # need to be the xsection
    thisdict = {
        1010: (findDIR(dirCeph+"VBFHphigamma"),xsecRun3['VBFH']*0.49),
        1020: (findDIR(dirCeph+"VBFHrhogamma"),xsecRun3['VBFH']),
        1030: (findDIR(dirCeph+"VBFHkstgamma"),xsecRun3['VBFH']*(2./3)),
        1040: (findDIR(dirCeph+"VBFHDstgamma"),xsecRun3['VBFH']*0.0389),
        #
        1017: (findDIR(dirCeph+"ggHphigamma"),xsecRun3['ggH']*0.49),
        1027: (findDIR(dirCeph+"ggHrhogamma"),xsecRun3['ggH']),
        1037: (findDIR(dirCeph+"ggHkstgamma"),xsecRun3['ggH']*(2./3)),
        1047: (findDIR(dirCeph+"ggHDstgamma"),xsecRun3['ggH']*0.0389),
        #
        1019: (findDIR(dirCeph+"Zphigamma_MadGraph"),xsecRun3['Z']*0.49),
        1029: (findDIR(dirCeph+"Zrhogamma_MadGraph"),xsecRun3['Z']),
        1039: (findDIR(dirCeph+"Zkstgamma"),xsecRun3['Z']*(2./3)),
        1049: (findDIR(dirCeph+"ZDstgamma"),xsecRun3['Z']*0.0389),
        9: (findDIR(dirCephGroup+sampleName+"_HT-40to70_TuneCP5_13p6TeV_madgraphMLM-pythia8+"+campaign),15310*1000*1.26), #LO *1.26
        10: (findDIR(dirCephGroup+sampleName+"_HT-70to100_TuneCP5_13p6TeV_madgraphMLM-pythia8+"+campaign),8111.0*1000*1.26), #LO *1.26
        11: (findDIR(dirCephGroup+sampleName+"_HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8+"+campaign),7373.0*1000*1.26), #LO *1.26
        12: (findDIR(dirCephGroup+sampleName+"_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8+"+campaign),1539.0*1000*1.26), #LO *1.26
        13: (findDIR(dirCephGroup+sampleName+"_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8+"+campaign),167.7*1000*1.26), #LO *1.26
        14: (findDIR(dirCephGroup+sampleName+"_HT-600_TuneCP5_13p6TeV_madgraphMLM-pythia8+"+campaign),54.48*1000*1.26) #LO *1.26
    }

    if(str(year) == '12022'):
        dict_ = {
            -61: (findDIR(dirNameDataSkims+str(year)+"/Tau+RunC"),-1),
            -62: (findDIR(dirNameDataSkims+str(year)+"/Tau+RunD"),-1),
        }
        thisdict.update(dict_)

    if(str(year) == '22022'):
        dict_ = {
            -63: (findDIR(dirNameDataSkims+str(year)+"/Tau+RunE"),-1),
            -64: (findDIR(dirNameDataSkims+str(year)+"/Tau+RunF"),-1),
            -65: (findDIR(dirNameDataSkims+str(year)+"/Tau+RunG"),-1),
        }
        thisdict.update(dict_)

    if(str(year) == '12023'):
        dict_ = {
            -61: (findDIR(dirNameDataSkims+str(year)+"/Tau+RunCv1"),-1),
            -62: (findDIR(dirNameDataSkims+str(year)+"/Tau+RunCv2"),-1),
            -63: (findDIR(dirNameDataSkims+str(year)+"/Tau+RunCv3"),-1),
            -64: (findDIR(dirNameDataSkims+str(year)+"/Tau+RunCv4"),-1),
        }
        thisdict.update(dict_)

    if(str(year) == '22023'):
        dict_ = {
            -65: (findDIR(dirNameDataSkims+str(year)+"/Tau+RunDv1"),-1),
            -66: (findDIR(dirNameDataSkims+str(year)+"/Tau+RunDv2"),-1),
        }
        thisdict.update(dict_)

    return thisdict

def SwitchSample(thisdict,argument):

    return thisdict.get(argument, "BKGdefault, xsecDefault")

def computeWeigths(rdf,xsec):

    genEventSumWeight = rdf.Sum("genEventSumw").GetValue()
    genEventSumNoWeight = rdf.Sum("genEventCount").GetValue()

    print('genEventSumWeight',genEventSumWeight)
    print('genEventSumNoWeight',genEventSumNoWeight)

    weight = xsec / genEventSumWeight

    lumiEq = (genEventSumNoWeight / xsec)
    print("lumi equivalent fb %s" %lumiEq)

    # neglecting the lumi
    #    weight = (1./genEventSumWeight)
    print('weight',weight)
    return weight


def pickTRG(overall,year,PDType,isVBF,isW,isZ,isZinv,isBPH):

    TRIGGER=''
    if(year == 2018 and isVBF):
        if (PDType== "EGamma"): TRIGGER=getTriggerFromJson(overall, "isVBF", year)
        elif (PDType== "Tau"): TRIGGER="{0} and not {1}".format(getTriggerFromJson(overall, "isZinv", year),getTriggerFromJson(overall, "isVBF", year))
        elif (PDType== "NULL"): TRIGGER=getTriggerFromJson(overall, "isVBFor", year)   # MC seems the same

    if((year == 2018 or year == 12022 or year == 22022 or year == 12023 or year == 22023) and isZinv):
        print('hello')
        if (PDType== "Tau"): TRIGGER=getTriggerFromJson(overall, "isZinv", year)
        elif (PDType== "NULL"): TRIGGER=getTriggerFromJson(overall, "isZinv", year)    # MC seems the same
        print('PDType=',PDType)
        print('year=',year)

    if(isBPH):
        if (PDType== "Charmonium" or PDType=="ParkingDoubleMuonLowMass"): TRIGGER=getTriggerFromJson(overall, "isBPH", year)
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
