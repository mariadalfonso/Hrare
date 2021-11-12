import ROOT
import os
import json
from subprocess import call,check_output

if "/functions.so" not in ROOT.gSystem.GetLibraries():
    ROOT.gSystem.CompileMacro("functions.cc","k")

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

def findDataset(name):

    DASclient = "dasgoclient -query '%(query)s'"
    cmd= DASclient%{'query':'file dataset=%s'%name}
    print(cmd)
    check_output(cmd,shell=True)
    fileList=[ 'root://xrootd-cms.infn.it//'+x for x in check_output(cmd,shell=True).split() ]

    files_ROOT = ROOT.vector('string')()
    for f in fileList: files_ROOT.push_back(f)

    return files_ROOT

def findDIR(directory):

    print(directory)

    counter = 0
    rootFiles = ROOT.vector('string')()
    for root, directories, filenames in os.walk(directory):
        for f in filenames:

            counter+=1
            filePath = os.path.join(os.path.abspath(root), f)
            if "failed/" in filePath: continue
            if "log/" in filePath: continue
            rootFiles.push_back(filePath)
#            if counter>100: break
#            if counter>50: break
#            if counter>5: break

    return rootFiles

def concatenate(result, tmp1):
    for f in tmp1:
        result.push_back(f)


def getMClist(sampleNOW):

    files = findDIR("{}".format(SwitchSample(sampleNOW)[0]))
    if sampleNOW==0:
            files1 = findDIR("{}".format(SwitchSample(1000)[0]))
            concatenate(files, files1)
    return files

def getDATAlist():

    files1 = findDIR("/mnt/T2_US_MIT/hadoop/cms/store/user/paus/nanohr/D00/SingleMuon+Run2018B-UL2018_MiniAODv2-v2+MINIAOD")
    files2 = findDIR("/mnt/T2_US_MIT/hadoop/cms/store/user/paus/nanohr/D00/SingleMuon+Run2018C-UL2018_MiniAODv2-v2+MINIAOD")
    files3 = findDIR("/mnt/T2_US_MIT/hadoop/cms/store/user/paus/nanohr/D00/SingleMuon+Run2018D-UL2018_MiniAODv2-v3+MINIAOD")
    files4 = findDIR("/mnt/T2_US_MIT/hadoop/cms/store/user/paus/nanohr/D00/SingleMuon+Run2018A-UL2018_MiniAODv2-v2+MINIAOD")

    files = ROOT.vector('string')()
    concatenate(files, files1)
    concatenate(files, files2)
    concatenate(files, files3)
    concatenate(files, files4)

    return files

def SwitchSample(argument):

    dirT2 = "/mnt/T2_US_MIT/hadoop/cms/store/user/paus/nanohr/D00/"
    dirLocal = "/work/submit/mariadlf/Hrare/OCT14/"

    switch = {
        1000: (dirT2+"DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM",6067*1000), #NNLO
        0: (dirT2+"DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1_ext1-v1+MINIAODSIM",6067*1000), #NNLO
        1: (dirT2+"ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM", 51.1*1000), #LO
        2: (dirT2+"WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM", 412.7*1000), #LO 
        3: (dirT2+"WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM",53870.0*1000), #LO
        4: (dirT2+"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM",88.2*1000), #NNLO
        5: (dirT2+"TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM",365.3452*1000), #NNLO
        ##
        6: (dirT2+"GJets_DR-0p4_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM",5031*1000), #LO
        7: (dirT2+"GJets_DR-0p4_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM",1126*1000), #LO
        8: (dirT2+"GJets_DR-0p4_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM",124.3*1000), #LO
        9: (dirT2+"GJets_DR-0p4_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM",40.76*1000), #LO
        ##
        10: (dirLocal+"ZLLphigamma_pythia8_genFix",0.10*1000), #xsec=1pb * 0.101 (Zll) * BR(Hphigamma)=1
        11: (dirLocal+"WLNUphigamma_pythia8",0.66*1000), #xsec = 2pb * 0.33 (Wl) * BR(Hphigamma)=1
        12: (dirLocal+"vbf-hphigamma-powheg",2.*1000) # xsec = 4pb * BR(Hphigamma)=1 BR(phi->kk)=0.49

    }
    return switch.get(argument, "BKGdefault, xsecDefault")


def plot(h,filename,doLogX,color):

   ROOT.gStyle.SetOptStat(1);
   ROOT.gStyle.SetTextFont(42)
   c = ROOT.TCanvas("c", "", 800, 700)
   if doLogX: c.SetLogx();
#  c.SetLogy()

   h.SetTitle("")
   h.GetXaxis().SetTitleSize(0.04)
   h.GetYaxis().SetTitleSize(0.04)
   h.SetLineColor(color)

   h.Draw()

   label = ROOT.TLatex(); label.SetNDC(True)
   label.DrawLatex(0.175, 0.740, "#eta")
   label.DrawLatex(0.205, 0.775, "#rho,#omega")
   label.DrawLatex(0.270, 0.740, "#phi")
   label.SetTextSize(0.040); label.DrawLatex(0.100, 0.920, "#bf{CMS Simulation}")
   label.SetTextSize(0.030); label.DrawLatex(0.630, 0.920, "#sqrt{s} = 13 TeV, L_{int} = X fb^{-1}")

   print("saving file: {} ".format(filename))

   c.SaveAs(filename)
