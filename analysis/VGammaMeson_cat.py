import ROOT
import os
import sys
import json

ROOT.ROOT.EnableImplicitMT()
from utilsHrare import getMClist, getDATAlist, getSkims, computeWeigths
from utilsAna import getMesonFromJson, pickTRG, getMVAFromJson, getTriggerFromJson
from utilsAna import loadCorrectionSet, loadPolarizationTree, loadSFhisto, loadUserCode, loadtmva_helper
import tmva_helper_xml

doSyst = True
doPol = True
doMVA = True
if sys.argv[1]=='isZtag': doMVA = False
if sys.argv[1]=='isWtag': doMVA = False
#if sys.argv[2]=='isOmega3PiCat': doMVA = False
#if sys.argv[2]=='isPhi3PiCat': doMVA = False
#if sys.argv[2]=='isD0StarCat': doMVA = False
if sys.argv[2]=='isD0Pi0StarCat': doMVA = False
if sys.argv[2]=='isOmega3PiCat': doPol = False
if sys.argv[2]=='isPhi3PiCat': doPol = False
if sys.argv[2]=='isD0StarCat': doPol = False
if sys.argv[2]=='isD0Pi0StarCat': doPol = False
doPlot = True
doTrigger = False
doMesonMassSB = False

useD03=False
if sys.argv[2]=='isOmega3PiCat' or sys.argv[2]=='isPhi3PiCat' or sys.argv[2]=='isD0StarCat' or sys.argv[2]=='isD0Pi0StarCat': useD03=True

isGF = False
isZinv = False
isZ = False
isW = False
isVBF = False
isVBFlow = False

isPhiCat = "false"
isRhoCat = "false"
isK0StarCat = "false"
isOmega3PiCat = "false"
isPhi3PiCat = "false"
isD0StarCat = "false"
isD0Pi0StarCat = "false"

if sys.argv[1]=='isGFtag': isGF = True
if sys.argv[1]=='isZinvtag': isZinv = True
if sys.argv[1]=='isZtag': isZ = True
if sys.argv[1]=='isWtag': isW = True
if sys.argv[1]=='isVBFtag': isVBF = True
if sys.argv[1]=='isVBFtaglow': isVBFlow = True

if sys.argv[2]=='isPhiCat': isPhiCat = "true"
if sys.argv[2]=='isRhoCat': isRhoCat = "true"
if sys.argv[2]=='isK0StarCat': isK0StarCat = "true"
if sys.argv[2]=='isOmega3PiCat': isOmega3PiCat = "true"
if sys.argv[2]=='isPhi3PiCat': isPhi3PiCat = "true"
if sys.argv[2]=='isD0StarCat': isD0StarCat = "true"
if sys.argv[2]=='isD0Pi0StarCat': isD0Pi0StarCat = "true"

if sys.argv[4]=='2018': year = 2018
if sys.argv[4]=='2017': year = 2017
if sys.argv[4]=='22016': year = 22016 #F-H
if sys.argv[4]=='12016': year = 12016 #B-F

lumis={
    '12016': 19.52, #APV #(B-F for 2016 pre)
    '22016': 16.80, #postVFP
    '2016': 35.9,
    '2017': 41.5,
    '12017': 7.7, #(F for 2017)
    '2018': 59.70,
    '12018': 39.54,
    'all': 86.92,      #19.52 + 7.7 + 59.70
}

DEEP_B_LOOSE={
    '2018': 0.1208,
    '2017': 0.1355,
    '22016': 0.1918,
    '12016': 0.2027,
}

DEEP_B_MEDIUM={
    '2018': 0.4148,
    '2017': 0.4506,
    '22016': 0.5847,
    '12016': 0.6001,
}

DEEP_B_TIGHT={
    '2018': 0.7665,
    '2017': 0.7738,
    '22016': 0.8767,
    '12016': 0.8819,
}

#$$$$
#$$$$
#$$$$

PRESELECTION = "(nPhoton>0 && (nphi>0 or nrho>0 or nK0Star>0 or nomega>0) && PV_npvsGood>0)"
CLEAN_LepMes = "{}".format("(Sum(goodMeson)>0 and isMuorEle==1) ? deltaR(Muon_eta[goodMuons][0], Muon_phi[goodMuons][0], goodMeson_eta[index_pair[0]], goodMeson_phi[index_pair[0]]):(Sum(goodMeson)>0 and isMuorEle==2) ? deltaR(Electron_eta[goodElectrons][0], Electron_phi[goodElectrons][0], goodMeson_eta[index_pair[0]], goodMeson_phi[index_pair[0]]): -999")

CLEAN_JetMes = "{}".format("Sum(goodMeson)>0 ? std::min(deltaR(Jet_eta[goodJets][0], Jet_phi[goodJets][0], goodMeson_eta[index_pair[0]], goodMeson_phi[index_pair[0]]),deltaR(Jet_eta[goodJets][1], Jet_phi[goodJets][1], goodMeson_eta[index_pair[0]], goodMeson_phi[index_pair[0]])):-999")

CLEAN_JetPH = "{}".format("Sum(goodPhotons)>0 ? std::min(deltaR(Jet_eta[goodJets][0], Jet_phi[goodJets][0], goodPhotons_eta[index_pair[1]], goodPhotons_phi[index_pair[1]]),deltaR(Jet_eta[goodJets][1], Jet_phi[goodJets][1], goodPhotons_eta[index_pair[1]], goodPhotons_phi[index_pair[1]])):-999")

with open("/home/submit/mariadlf/Hrare/CMSSW_10_6_27_new/src/Hrare/analysis/config/selection.json") as jsonFile:
    jsonObject = json.load(jsonFile)
    jsonFile.close()

with open("/home/submit/mariadlf/Hrare/CMSSW_10_6_27_new/src/Hrare/analysis/config/trigger.json") as trgJsonFile:
    trgObject = json.load(trgJsonFile)
    trgJsonFile.close()

GOODJETS = jsonObject['GOODJETS']
LOOSEmuons = jsonObject['LOOSEmuons']
LOOSEelectrons = jsonObject['LOOSEelectrons']
GOODMUON = jsonObject['GOODMUON']
GOODELE = jsonObject['GOODELE']
JSON = jsonObject['JSON']
BARRELphotons = jsonObject['BARRELphotons']
ENDCAPphotons = jsonObject['ENDCAPphotons']
ENDCAPphotonsLoose = jsonObject['ENDCAPphotonsLoose']
photonsLoose = jsonObject['photonsLOOSE']

METFLAG = jsonObject['METFLAG']

MVA = jsonObject['MVAweights']
TRIGGERS = trgObject['triggers']
mesons = jsonObject['mesons']

#$$$$
#$$$$
#$$$$

def selectionTAG(df, doSyst, isData):

    if isZ:
        dftag = (df.Define("goodMuons","{}".format(GOODMUON)+" and Muon_mediumId and Muon_pfRelIso04_all < 0.25") # iso same as the loose
                 .Define("ele_mask", "cleaningMask(Photon_electronIdx[goodPhotons],nElectron)")
                 .Define("goodElectrons","{}".format(GOODELE)+" and Electron_mvaFall17V2Iso_WP90") # medium
                 .Define("looseMu","{}".format(LOOSEmuons))
                 .Define("looseEle","{}".format(LOOSEelectrons))
                 .Filter("Sum(looseEle)+Sum(looseMu)==2", "at least two muons or electrons, and no extra loose leptons")
                 .Define("isMuorEle","(Sum(looseMu)==2 and Sum(goodMuons)>=1)?1:(Sum(looseEle)==2 and Sum(goodElectrons)>=1)?2 :0")
                 .Filter("isMuorEle>0","at least one leading lepton")
                 .Define("Z_veto1", "Sum(looseEle)==2 ? Minv2(Electron_pt[looseEle][0], Electron_eta[looseEle][0], Electron_phi[looseEle][0], Electron_mass[looseEle][0], photon_pt, goodPhotons_eta[index_pair[1]], goodPhotons_phi[index_pair[1]]).first: -1")
                 .Define("Z_veto2", "Sum(looseEle)==2 ? Minv2(Electron_pt[looseEle][1], Electron_eta[looseEle][1], Electron_phi[looseEle][1], Electron_mass[looseEle][1], photon_pt, goodPhotons_eta[index_pair[1]], goodPhotons_phi[index_pair[1]]).first: -1")
                 .Filter("abs(Z_veto1-91) > 5 and abs(Z_veto2-91) > 5","kill the Z recontructed as gamma + electron")
                 .Define("V_mass", "(Sum(looseMu)==2 and Sum(Muon_charge[looseMu])==0)? Minv(Muon_pt[looseMu], Muon_eta[looseMu], Muon_phi[looseMu], Muon_mass[looseMu]) : (Sum(looseEle)==2 and Sum(Electron_charge[looseEle])==0) ? Minv(Electron_pt[looseEle], Electron_eta[looseEle], Electron_phi[looseEle], Electron_mass[looseEle]): 0.")
                 .Filter("(V_mass>(91-10) and V_mass<(91+15))","At least one good Z")
                 .Define("Visr_mass", "(Sum(looseMu)==2 and Sum(Muon_charge[looseMu])==0)? Minv3(Muon_pt[looseMu], Muon_eta[looseMu], Muon_phi[looseMu], Muon_mass[looseMu], photon_pt, goodPhotons_eta[index_pair[1]], goodPhotons_phi[index_pair[1]]) : (Sum(looseEle)==2 and Sum(Electron_charge[looseEle])==0) ? Minv3(Electron_pt[looseEle], Electron_eta[looseEle], Electron_phi[looseEle], Electron_mass[looseEle], photon_pt, goodPhotons_eta[index_pair[1]], goodPhotons_phi[index_pair[1]]): 0.")
                 .Filter("Visr_mass>91+5")
#                 .Define("Mu1_hasTriggerMatch", "hasTriggerMatch(Muon_eta[goodMuons][0], Muon_phi[goodMuons][0], TrigObj_eta, TrigObj_phi)")
#                 .Define("Mu2_hasTriggerMatch", "hasTriggerMatch(Muon_eta[goodMuons][1], Muon_phi[goodMuons][1], TrigObj_eta, TrigObj_phi)")
#                 .Define("Ele1_hasTriggerMatch", "Sum(goodElectrons)>1 ? hasTriggerMatch(Electron_eta[goodElectrons][0], Electron_phi[goodElectrons][0], TrigObj_eta, TrigObj_phi) : 0")
#                 .Define("Ele2_hasTriggerMatch", "Sum(goodElectrons)>1 ? hasTriggerMatch(Electron_eta[goodElectrons][1], Electron_phi[goodElectrons][1], TrigObj_eta, TrigObj_phi) : 0")
#                 .Filter("trigger>0 and ((Mu1_hasTriggerMatch and Muon_pt[goodMuons][0]>26) or (Mu2_hasTriggerMatch and Muon_pt[goodMuons][1]>26))","pass trigger")
                 .Define("LeadingLeptonPt", "isMuorEle==1 ? Muon_pt[looseMu][0]: isMuorEle==2 ? Electron_pt[looseEle][0] :0")
                 .Define("SubLeadingLeptonPt", "isMuorEle==1 ? Muon_pt[looseMu][1]: isMuorEle==2 ? Electron_pt[looseEle][1] :0")
                 .Define("LeadingLeptonEta", "isMuorEle==1 ? Muon_eta[looseMu][0]: isMuorEle==2 ? Electron_eta[looseEle][0] :0")
                 .Define("SubLeadingLeptonEta", "isMuorEle==1 ? Muon_eta[looseMu][1]: isMuorEle==2 ? Electron_eta[looseEle][1] :0")
                 )
        return dftag

    if isW:

        dftag = (df.Define("goodMuons","{}".format(GOODMUON)+" and Muon_tightId and Muon_pfRelIso04_all < 0.15") ## tight
                 .Define("ele_mask", "cleaningMask(Photon_electronIdx[goodPhotons],nElectron)")
                 .Define("goodElectrons","{}".format(GOODELE)+" and Electron_mvaFall17V2Iso_WP80 and Electron_pt>30") ## tight
                 .Define("vetoEle","{}".format(LOOSEelectrons))
                 .Define("vetoMu","{}".format(LOOSEmuons))
                 .Filter("(Sum(goodMuons)+Sum(goodElectrons))==1 and (Sum(vetoEle)+Sum(vetoMu))==1","one lepton")
                 .Define("isMuorEle","Sum(goodMuons)==1?1: Sum(goodElectrons)==1?2 :0")
                 .Define("V_mass","Sum(goodMuons)>0 ? mt(Muon_pt[goodMuons][0], Muon_phi[goodMuons][0], DeepMETResolutionTune_pt, DeepMETResolutionTune_phi) : mt(Electron_pt[goodElectrons][0], Electron_phi[goodElectrons][0], DeepMETResolutionTune_pt, DeepMETResolutionTune_phi)")
                 .Define("Z_veto", "Sum(goodElectrons)==1 ? Minv2(Electron_pt[goodElectrons][0], Electron_eta[goodElectrons][0], Electron_phi[goodElectrons][0], Electron_mass[goodElectrons][0], photon_pt,goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]]).first: -1")
                 .Filter("abs(Z_veto-91) > 10","kill the Z recontructed as gamma + electron")
#                 .Filter("DeepMETResolutionTune_pt>15","MET>15")
#                 .Filter("V_mass>15","MT>15")
#                 .Define("trigger","{}".format(TRIGGER))
#                 .Define("Mu1_hasTriggerMatch", "hasTriggerMatch(Muon_eta[goodMuons][0], Muon_phi[goodMuons][0], TrigObj_eta, TrigObj_phi)")
#                 .Filter("trigger>0 and Mu1_hasTriggerMatch and Muon_pt[goodMuons][0]>26","pass trigger")
                 .Define("deltaLepMeson","{}".format(CLEAN_LepMes))
                 .Filter("deltaLepMeson>0.5","kill the muon reconstructed as meson")
                 .Define("dPhiGammaMET","abs(deltaPhi(goodPhotons_phi[index_pair[1]], DeepMETResolutionTune_phi))")
                 .Define("dPhiMesonMET","abs(deltaPhi(goodMeson_phi[index_pair[0]], DeepMETResolutionTune_phi))")
                 .Define("dPhiLeptonMET","isMuorEle==1 ? abs(deltaPhi(Muon_phi[goodMuons][0], DeepMETResolutionTune_phi)): isMuorEle==2 ? abs(deltaPhi(Electron_phi[goodElectrons][0], DeepMETResolutionTune_phi)): 999.")
                 .Define("LeadingLeptonPt", "isMuorEle==1 ? Muon_pt[goodMuons][0]: isMuorEle==2 ? Electron_pt[goodElectrons][0] :0")
                 .Define("LeadingLeptonEta", "isMuorEle==1 ? Muon_eta[goodMuons][0]: isMuorEle==2 ? Electron_eta[goodElectrons][0] :0")
                 )
        return dftag

    if isVBF or isVBFlow:

        VBFcut = "mJJ>300 and dEtaJJ>3 and Y1Y2<0 and Jet_pt[goodJets][0]>40"
        if isVBFlow: VBFcut = "mJJ>250 and dEtaJJ>3. and Y1Y2<0"

# tight means less PU
## default is the medium (Jet_puId == 2)
## use the tight PU id (Jet_puId == 1) for 2017-2018: note those are swapped for the 2016
## use the tight PU id (Jet_puId == 4) for 2016-12016
## https://github.com/cms-nanoAOD/cmssw/issues/583
## flag means passlooseID*4+passmediumID*2+passtightID*1.
## https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID#miniAOD_and_nanoAOD

        PUjetID = "true"
#        if year == 2018 or year == 2017: PUjetID = "(((Jet_puId & 1) and abs(Jet_eta)>2.75) or ((Jet_puId & 2) and abs(Jet_eta)<=2.75))"
#        if year == 2016 or year == 12016: PUjetID = "(((Jet_puId & 4) and abs(Jet_eta)>2.75 or ((Jet_puId & 2) and abs(Jet_eta)<=2.75))"

        if (doSyst and isData == "false"):
            df = df.Define("Jet_delta",'computeJECuncertainties(corr_sf, Jet_pt, Jet_eta)').Vary("Jet_pt", "ROOT::RVec<ROOT::RVecF>{Jet_pt*(1-Jet_delta),Jet_pt*(1+Jet_delta)}", variationTags=["dn","up"], variationName = "JetSYST")

        dftag = (df.Define("goodJets","{}".format(GOODJETS)+" and {}".format(PUjetID))
                 .Define("nGoodJets","Sum(goodJets)*1.0f").Filter("Sum(goodJets)>1","two jets for VBF")
                 .Define("mJJ","Minv(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets])")
                 .Define("dEtaJJ","abs(Jet_eta[goodJets][0] - Jet_eta[goodJets][1])")
                 .Define("dPhiJJ","abs(deltaPhi(Jet_phi[goodJets][0],Jet_phi[goodJets][1]))")
                 .Define("Y1Y2","Jet_eta[goodJets][0]*Jet_eta[goodJets][1]")
                 .Filter("{}".format(VBFcut),"Filter on MJJ , Deta, Y1Y2")
                 .Define("ele_mask", "cleaningMask(Photon_electronIdx[goodPhotons],nElectron)")
                 .Define("vetoEle","{}".format(LOOSEelectrons))
                 .Define("vetoMu","{}".format(LOOSEmuons))
                 .Filter("(Sum(vetoEle)+Sum(vetoMu))==0", "no leptons")
#                 .Define("trigger","{}".format(TRIGGER))
#                 .Filter("trigger>0", "pass triggers")
                 .Define("SoftActivityJetNjets5F","SoftActivityJetNjets5*1.0f")
                 .Define("jet1Pt","Jet_pt[goodJets][0]")
                 .Define("jet2Pt","Jet_pt[goodJets][1]")
                 .Define("jet1Eta","Jet_eta[goodJets][0]")
                 .Define("jet2Eta","Jet_eta[goodJets][1]")
                 .Define("jet1hfsigmaPhiPhi","Jet_hfsigmaPhiPhi[goodJets][0]")
                 .Define("jet1hfsigmaEtaEta","Jet_hfsigmaEtaEta[goodJets][0]")
                 .Define("jet2hfsigmaPhiPhi","Jet_hfsigmaPhiPhi[goodJets][1]")
                 .Define("jet2hfsigmaEtaEta","Jet_hfsigmaEtaEta[goodJets][1]")
                 .Define("deltaJetMeson","{}".format(CLEAN_JetMes))
                 .Define("deltaJetPhoton","{}".format(CLEAN_JetPH))
                 .Define("zepVar","compute_jet_HiggsVars_var(Jet_pt[goodJets],Jet_eta[goodJets],Jet_phi[goodJets],Jet_mass[goodJets], photon_pt,goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]], goodMeson_pt[index_pair[0]],goodMeson_eta[index_pair[0]], goodMeson_phi[index_pair[0]], goodMeson_mass[index_pair[0]], 0)")
                 .Define("detaHigJet1","compute_jet_HiggsVars_var(Jet_pt[goodJets],Jet_eta[goodJets],Jet_phi[goodJets],Jet_mass[goodJets], photon_pt,goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]], goodMeson_pt[index_pair[0]],goodMeson_eta[index_pair[0]], goodMeson_phi[index_pair[0]], goodMeson_mass[index_pair[0]], 1)")
                 .Define("detaHigJet2","compute_jet_HiggsVars_var(Jet_pt[goodJets],Jet_eta[goodJets],Jet_phi[goodJets],Jet_mass[goodJets], photon_pt,goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]], goodMeson_pt[index_pair[0]],goodMeson_eta[index_pair[0]], goodMeson_phi[index_pair[0]], goodMeson_mass[index_pair[0]], 2)")
#                 .Filter("DeepMETResolutionTune_pt<75","DeepMETResolutionTune_pt<75") # not doing Zinv as nominal
                 )
        return dftag

    if isZinv:
        dftag = (df.Define("ele_mask", "cleaningMask(Photon_electronIdx[goodPhotons],nElectron)")
                 .Define("vetoEle","{}".format(LOOSEelectrons))
                 .Define("vetoMu","{}".format(LOOSEmuons))
                 .Filter("(Sum(vetoEle)+Sum(vetoMu))==0", "no leptons")
#                 .Define("trigger","{}".format(TRIGGER))
#                 .Filter("trigger>0", "pass triggers")
                 .Filter("DeepMETResolutionTune_pt>50","MET>50")
                 .Define("metFilter","{}".format(METFLAG))
                 .Filter("metFilter", "pass METfilter")
                 .Define("dPhiGammaMET","abs(deltaPhi(goodPhotons_phi[index_pair[1]], DeepMETResolutionTune_phi))")
                 .Define("dPhiMesonMET","abs(deltaPhi(goodMeson_phi[index_pair[0]], DeepMETResolutionTune_phi))")
                 .Define("ptRatioMEThiggs","abs(DeepMETResolutionTune_pt-HCandPT)/HCandPT")
#                 .Filter("ptRatioMEThiggs<0.8","ptRatioMEThiggs<0.8")
                 .Filter("dPhiGammaMET>1","dPhiGammaMET>1")
                 .Filter("dPhiMesonMET>1","dPhiMesonMET>1")
                 ##
                 .Define("goodJets","{}".format(GOODJETS))
                 .Define("nGoodJets","Sum(goodJets)*1.0f")
                 .Define("bjet", "Jet_btagDeepB[goodJets] > {}".format(DEEP_B_MEDIUM['2018']))
                 .Define("nbtag", "Sum(bjet)*1.0f")
                 .Define("PV_npvsGoodF","PV_npvsGood*1.0f")
                 .Define("SoftActivityJetNjets5F","SoftActivityJetNjets5*1.0f")
                 )
        return dftag

    if isGF:

        if (doSyst and isData == "false"):
            df = df.Define("Jet_delta",'computeJECuncertainties(corr_sf, Jet_pt, Jet_eta)').Vary("Jet_pt", "ROOT::RVec<ROOT::RVecF>{Jet_pt*(1-Jet_delta),Jet_pt*(1+Jet_delta)}", variationTags=["dn","up"], variationName = "JetSYST")

        dftag = (df.Define("ele_mask", "cleaningMask(Photon_electronIdx[goodPhotons],nElectron)")
                 .Define("vetoEle","{}".format(LOOSEelectrons))
                 .Define("vetoMu","{}".format(LOOSEmuons))
                 .Filter("(Sum(vetoEle)+Sum(vetoMu))==0", "no leptons")
                 #                 .Define("trigger","{}".format(TRIGGER))
                 #                 .Filter("trigger>0", "pass triggers")
                 #.Filter("DeepMETResolutionTune_pt<75","DeepMETResolutionTune_pt<75") # not doing Zinv as nominal
                 .Define("goodJets","{}".format(GOODJETS))
                 .Define("nGoodJets","Sum(goodJets)*1.0f")
                 .Define("SoftActivityJetNjets5F","SoftActivityJetNjets5*1.0f")
#                 .Filter("Sum(goodJets)<2 or (Sum(goodJets)>=2 and Jet_pt[goodJets][0]<30) or (Sum(goodJets)>=2 and Jet_eta[goodJets][0]*Jet_eta[goodJets][1]>0) or (Sum(goodJets)>=2 and Jet_eta[goodJets][0]*Jet_eta[goodJets][1]<0 and abs(Jet_eta[goodJets][0] - Jet_eta[goodJets][1])<3 )","0 or 1 jets (pt20, |eta|<4.7) or >=2 with dEta<3")
                 .Define("mJJ","Sum(goodJets)>=2?Minv(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets]):0.")
                 .Filter("Sum(goodJets)<2 or (Sum(goodJets)>=2 and Jet_pt[goodJets][0]<30) or (Sum(goodJets)>=2 and Jet_eta[goodJets][0]*Jet_eta[goodJets][1]>0) or (Sum(goodJets)>=2 and Jet_eta[goodJets][0]*Jet_eta[goodJets][1]<0 and abs(Jet_eta[goodJets][0] - Jet_eta[goodJets][1])<3 ) or (Sum(goodJets)>=2 and mJJ<300) ","0 or 1 jets (pt20, |eta|<4.7) or >=2 with dEta<3 or >=2 with mJJ<300")
                 )
        return dftag

def dfGammaMeson(df,PDType):

    TRIGGER=pickTRG(TRIGGERS,year,PDType,isVBF,isW,isZ,(isZinv or isVBFlow or isGF))

    GOODphotons = ""
    if(isGF): GOODphotons = "({0} or {1}) and Photon_pt>38 and Photon_electronVeto and abs(Photon_eta)<2.1 and Photon_mvaID_WP80".format(BARRELphotons,ENDCAPphotons) #90-80
    if(isVBF): GOODphotons = "{} and Photon_pt>75 and Photon_electronVeto".format(BARRELphotons) #90
    if(isVBFlow): GOODphotons = "({0} or {1}) and Photon_pt>38 and Photon_pt<75 and Photon_electronVeto and abs(Photon_eta)<2.1 and Photon_mvaID_WP80".format(BARRELphotons,ENDCAPphotons) #90-80
    if(isZinv): GOODphotons = "({0} or {1}) and Photon_pt>38 and abs(Photon_eta)<2.1 and Photon_electronVeto".format(BARRELphotons,ENDCAPphotons) #90-80
    if(isW or isZ): GOODphotons = "({0} or {1}) and (Photon_pixelSeed == false)".format(BARRELphotons,ENDCAPphotonsLoose) #90-90
    print("PHOTONS = ", GOODphotons)

    dfOBJ = (df.Filter("nPhoton>0 and PV_npvsGood>0","photon from nano >0 and PV_npvsGood > 0")
             .Define("triggerAna","{}".format(TRIGGER))
             .Filter("triggerAna>0", "pass triggers")  ## comment while doing trigger studies
             .Define("loosePhotons", "{}".format(photonsLoose))
             .Define("nPhotonsVeto","Sum(loosePhotons)")
             .Define("goodPhotons", "{}".format(GOODphotons))
             .Define("nGoodPhotons","Sum(goodPhotons)*1.0f")
             .Filter("Sum(goodPhotons)>0", "At least one good Photon")
             .Define("goodPhotons_pt", "Photon_pt[goodPhotons]")
             .Define("goodPhotons_eta", "Photon_eta[goodPhotons]")
             .Define("goodPhotons_phi", "Photon_phi[goodPhotons]")
             .Define("goodPhotons_pfRelIso03_all", "Photon_pfRelIso03_all[goodPhotons]")
             .Define("goodPhotons_pfRelIso03_chg", "Photon_pfRelIso03_chg[goodPhotons]")
             .Define("goodPhotons_hoe", "Photon_hoe[goodPhotons]")
#             .Define("goodPhotons_pixelSeed", "Photon_pixelSeed[goodPhotons]")
#             .Define("goodPhotons_r9", "Photon_r9[goodPhotons]")
#             .Define("goodPhotons_sieie", "Photon_sieie[goodPhotons]")
#             .Define("goodPhotons_x_calo", "Photon_x_calo[goodPhotons]")
#             .Define("goodPhotons_y_calo", "Photon_y_calo[goodPhotons]")
#             .Define("goodPhotons_z_calo", "Photon_z_calo[goodPhotons]")
             .Define("goodPhotons_mvaID", "Photon_mvaID[goodPhotons]")
             .Define("goodPhotons_energyErr", "Photon_energyErr[goodPhotons]")
             .Define("goodPhotons_isScEtaEB", "Photon_isScEtaEB[goodPhotons]")
             .Define("jet_mask", "cleaningMask(Photon_jetIdx[goodPhotons],nJet)")
             )

#    if useD03:
#        dfOBJ.Define("goodPhotons_x_calo", "Photon_x_calo[goodPhotons]").Define("goodPhotons_y_calo", "Photon_y_calo[goodPhotons]").Define("goodPhotons_z_calo", "Photon_z_calo[goodPhotons]")

    return dfOBJ

def dfHiggsCand(df,isData):

    GOODPHI = ""
    if doMesonMassSB:
        if(isVBF): GOODPHI = "{}".format(getMesonFromJson(mesons, "isVBF", "isPhiCatMassSB"))
        if(isVBFlow): GOODPHI = "{}".format(getMesonFromJson(mesons, "isVBFlow" , "isPhiCatMassSB"))
        if(isZinv or isGF): GOODPHI = "{}".format(getMesonFromJson(mesons, "isZinv", "isPhiCatMassSB"))
        if(isW or isZ): GOODPHI = "{}".format(getMesonFromJson(mesons, "VH", "isPhiCatMassSB"))
    else:
        if(isVBF): GOODPHI = "{}".format(getMesonFromJson(mesons, "isVBF", "isPhiCat"))
        if(isVBFlow): GOODPHI = "{}".format(getMesonFromJson(mesons, "isVBFlow" , "isPhiCat"))
        if(isZinv or isGF): GOODPHI = "{}".format(getMesonFromJson(mesons, "isZinv", "isPhiCat"))
        if(isW or isZ): GOODPHI = "{}".format(getMesonFromJson(mesons, "VH", "isPhiCat"))

    GOODRHO = ""
    if doMesonMassSB:
        if(isVBF): GOODRHO = "{}".format(getMesonFromJson(mesons, "isVBF", "isRhoCatMassSB"))
        if(isVBFlow): GOODRHO = "{}".format(getMesonFromJson(mesons, "isVBFlow" , "isRhoCatMassSB"))
        if(isZinv or isGF): GOODRHO = "{}".format(getMesonFromJson(mesons, "isZinv", "isRhoCatMassSB"))
        if(isW or isZ): GOODRHO = "{}".format(getMesonFromJson(mesons, "VH", "isRhoCatMassSB"))
    else:
        if(isVBF): GOODRHO = "{}".format(getMesonFromJson(mesons, "isVBF", "isRhoCat"))
        if(isVBFlow): GOODRHO = "{}".format(getMesonFromJson(mesons, "isVBFlow" , "isRhoCat"))
        if(isZinv or isGF): GOODRHO = "{}".format(getMesonFromJson(mesons, "isZinv", "isRhoCat"))
        if(isW or isZ): GOODRHO = "{}".format(getMesonFromJson(mesons, "VH", "isRhoCat"))

    GOODK0STAR = ""
    if doMesonMassSB:
        if(isVBF): GOODK0STAR = "{}".format(getMesonFromJson(mesons, "isVBF", "isK0StarCatMassSB"))
        if(isVBFlow): GOODK0STAR = "{}".format(getMesonFromJson(mesons, "isVBFlow" , "isK0StarCatMassSB"))
        if(isZinv or isGF): GOODK0STAR = "{}".format(getMesonFromJson(mesons, "isZinv", "isK0StarCatMassSB"))
        if(isW or isZ): GOODK0STAR = "{}".format(getMesonFromJson(mesons, "VH", "isK0StarCatMassSB"))
    else:
        if(isVBF): GOODK0STAR = "{}".format(getMesonFromJson(mesons, "isVBF", "isK0StarCat"))
        if(isVBFlow): GOODK0STAR = "{}".format(getMesonFromJson(mesons, "isVBFlow" , "isK0StarCat"))
        if(isZinv or isGF): GOODK0STAR = "{}".format(getMesonFromJson(mesons, "isZinv", "isK0StarCat"))
        if(isW or isZ): GOODK0STAR = "{}".format(getMesonFromJson(mesons, "VH", "isK0StarCat"))

    GOODOMEGA3Pi = ""
    if(isVBF): GOODOMEGA3Pi = "{}".format(getMesonFromJson(mesons, "isVBF", "isOmega3PiCat"))
    if(isVBFlow): GOODOMEGA3Pi = "{}".format(getMesonFromJson(mesons, "isVBFlow" , "isOmega3PiCat"))
    if(isZinv or isGF): GOODOMEGA3Pi = "{}".format(getMesonFromJson(mesons, "isZinv", "isOmega3PiCat"))
    if(isW or isZ): GOODOMEGA3Pi = "{}".format(getMesonFromJson(mesons, "VH", "isOmega3PiCat"))

    GOODPHI3Pi = ""
    if(isVBF): GOODPHI3Pi = "{}".format(getMesonFromJson(mesons, "isVBF", "isPhi3PiCat"))
    if(isVBFlow): GOODPHI3Pi = "{}".format(getMesonFromJson(mesons, "isVBFlow" , "isPhi3PiCat"))
    if(isZinv or isGF): GOODPHI3Pi = "{}".format(getMesonFromJson(mesons, "isZinv", "isPhi3PiCat"))
    if(isW or isZ): GOODPHI3Pi = "{}".format(getMesonFromJson(mesons, "VH", "isPhi3PiCat"))

    GOODD0STAR = ""
    if(isVBF): GOODD0STAR = "{}".format(getMesonFromJson(mesons, "isVBF", "isD0StarCat"))
    if(isVBFlow): GOODD0STAR = "{}".format(getMesonFromJson(mesons, "isVBFlow" , "isD0StarCat"))
    if(isZinv or isGF): GOODD0STAR = "{}".format(getMesonFromJson(mesons, "isZinv", "isD0StarCat"))
    if(isW or isZ): GOODD0STAR = "{}".format(getMesonFromJson(mesons, "VH", "isD0StarCat"))

    GOODD0PI0STAR = ""
    if(isVBF): GOODD0PI0STAR = "{}".format(getMesonFromJson(mesons, "isVBF", "isD0Pi0StarCat"))
    if(isVBFlow): GOODD0PI0STAR = "{}".format(getMesonFromJson(mesons, "isVBFlow" , "isD0Pi0StarCat"))
    if(isZinv or isGF): GOODD0PI0STAR = "{}".format(getMesonFromJson(mesons, "isZinv", "isD0Pi0StarCat"))
    if(isW or isZ): GOODD0PI0STAR = "{}".format(getMesonFromJson(mesons, "VH", "isD0Pi0StarCat"))

    if(isPhiCat=="true"):

        print("PHI = ", GOODPHI)
        dfbase = (df.Filter("nphi>0").Define("goodMeson","({}".format(GOODPHI)+" && {}".format(isPhiCat)+")")
                  .Filter("Sum(goodMeson)>0", "one good Phi (ptPhi, validfit, ptTracks)")
                  .Define("goodMeson_pt", "phi_kin_pt[goodMeson]")
                  .Define("goodMeson_eta", "phi_kin_eta[goodMeson]")
                  .Define("goodMeson_phi", "phi_kin_phi[goodMeson]")
                  .Define("goodMeson_mass", "phi_kin_mass[goodMeson]")
                  .Define("goodMeson_iso", "phi_iso[goodMeson]")
                  .Define("goodMeson_isoNeuHad", "phi_isoNeuHad[goodMeson]")
                  .Define("goodMeson_vtx_chi2dof", "phi_kin_vtx_chi2dof[goodMeson]")
                  .Define("goodMeson_vtx_prob", "phi_kin_vtx_prob[goodMeson]")
                  .Define("goodMeson_sipPV", "phi_kin_sipPV[goodMeson]")
                  .Define("goodMeson_bestVtx_idx", "phi_bestVtx_idx[goodMeson]")
                  .Define("goodMeson_bestVtx_X", "phi_bestVtx_X[goodMeson]")
                  .Define("goodMeson_bestVtx_Y", "phi_bestVtx_Y[goodMeson]")
                  .Define("goodMeson_bestVtx_Z", "phi_bestVtx_Z[goodMeson]")
                  .Define("goodMeson_massErr", "phi_kin_massErr[goodMeson]")
                  .Define("goodMeson_trk1_pt", "phi_trk1_pt[goodMeson]")
                  .Define("goodMeson_trk2_pt", "phi_trk2_pt[goodMeson]")
                  .Define("goodMeson_trk1_eta", "phi_trk1_eta[goodMeson]")
                  .Define("goodMeson_trk2_eta", "phi_trk2_eta[goodMeson]")
                  .Define("goodMeson_DR","DeltaR(phi_trk1_eta[goodMeson],phi_trk2_eta[goodMeson],phi_trk1_phi[goodMeson],phi_trk2_phi[goodMeson])")
                  .Define("wrongMeson","({}".format(GOODRHO)+")")
                  .Define("wrongMeson_pt","Sum(wrongMeson) > 0 ? rho_kin_pt[wrongMeson]: ROOT::VecOps::RVec<float>(0.f)")
                  .Define("wrongMeson2","({}".format(GOODK0STAR)+")")
                  .Define("wrongMeson2_pt","Sum(wrongMeson2) > 0 ? K0Star_kin_pt[wrongMeson2]: ROOT::VecOps::RVec<float>(0.f)")
#                  .Define("wrongMeson2_pt", "ROOT::VecOps::RVec<float>(0.f)" )
                  )

    if(isRhoCat=="true"):

        print("RHO = ", GOODRHO)
        dfbase = (df.Filter("nrho>0").Define("goodMeson","({}".format(GOODRHO)+" && {}".format(isRhoCat)+")")
                  .Filter("Sum(goodMeson)>0", "one good Rho (ptPhi, validfit, ptTracks)")
                  .Define("goodMeson_pt", "rho_kin_pt[goodMeson]")
                  .Define("goodMeson_eta", "rho_kin_eta[goodMeson]")
                  .Define("goodMeson_phi", "rho_kin_phi[goodMeson]")
                  .Define("goodMeson_iso", "rho_iso[goodMeson]")
                  .Define("goodMeson_isoPho", "rho_isoPho[goodMeson]")
                  .Define("goodMeson_isoNeuHad", "rho_isoNeuHad[goodMeson]")
                  .Define("goodMeson_mass", "rho_kin_mass[goodMeson]")
                  .Define("goodMeson_vtx_chi2dof", "rho_kin_vtx_chi2dof[goodMeson]")
                  .Define("goodMeson_vtx_prob", "rho_kin_vtx_prob[goodMeson]")
                  .Define("goodMeson_sipPV", "rho_kin_sipPV[goodMeson]")
                  .Define("goodMeson_bestVtx_idx", "rho_bestVtx_idx[goodMeson]")
                  .Define("goodMeson_bestVtx_X", "rho_bestVtx_X[goodMeson]")
                  .Define("goodMeson_bestVtx_Y", "rho_bestVtx_Y[goodMeson]")
                  .Define("goodMeson_bestVtx_Z", "rho_bestVtx_Z[goodMeson]")
                  .Define("goodMeson_massErr", "rho_kin_massErr[goodMeson]")
                  .Define("goodMeson_trk1_pt", "rho_trk1_pt[goodMeson]")
                  .Define("goodMeson_trk2_pt", "rho_trk2_pt[goodMeson]")
                  .Define("goodMeson_trk1_eta", "rho_trk1_eta[goodMeson]")
                  .Define("goodMeson_trk2_eta", "rho_trk2_eta[goodMeson]")
                  .Define("goodMeson_DR","DeltaR(rho_trk1_eta[goodMeson],rho_trk2_eta[goodMeson],rho_trk1_phi[goodMeson],rho_trk2_phi[goodMeson])")
                  .Define("wrongMeson","({}".format(GOODPHI)+")")
                  .Define("wrongMeson_pt","Sum(wrongMeson) > 0 ? phi_kin_pt[wrongMeson]: ROOT::VecOps::RVec<float>(0.f)")
                  .Define("wrongMeson2","({}".format(GOODK0STAR)+")")
                  .Define("wrongMeson2_pt","Sum(wrongMeson2) > 0 ? K0Star_kin_pt[wrongMeson2]: ROOT::VecOps::RVec<float>(0.f)")
#                  .Define("wrongMeson2_pt", "ROOT::VecOps::RVec<float>(0.f)" )
                  )

    if(isK0StarCat=="true"):

        print("K0Star = ", GOODK0STAR)
        dfbase = (df.Filter("nK0Star>0")
                  .Define("K0Star_mask", "cleaningPair(event,nK0Star,K0Star_kaon_phi,K0Star_kaon_eta,K0Star_kaon_pt,K0Star_pion_phi,K0Star_pion_eta,K0Star_pion_pt,K0Star_kin_mass)")
                  .Define("goodMeson","({}".format(GOODK0STAR)+" && {}".format(isK0StarCat)+" && K0Star_mask)")
                  .Filter("Sum(goodMeson)>0", "one good K0Star (ptPhi, validfit, ptTracks)")
                  .Define("goodMeson_pt", "K0Star_kin_pt[goodMeson]")
                  .Define("goodMeson_eta", "K0Star_kin_eta[goodMeson]")
                  .Define("goodMeson_phi", "K0Star_kin_phi[goodMeson]")
                  .Define("goodMeson_iso", "K0Star_iso[goodMeson]")
                  .Define("goodMeson_isoNeuHad", "K0Star_isoNeuHad[goodMeson]")
                  .Define("goodMeson_mass", "K0Star_kin_mass[goodMeson]")
                  .Define("goodMeson_vtx_chi2dof", "K0Star_kin_vtx_chi2dof[goodMeson]")
                  .Define("goodMeson_vtx_prob", "K0Star_kin_vtx_prob[goodMeson]")
                  .Define("goodMeson_sipPV", "K0Star_kin_sipPV[goodMeson]")
                  .Define("goodMeson_bestVtx_idx", "K0Star_bestVtx_idx[goodMeson]")
                  .Define("goodMeson_bestVtx_X", "K0Star_bestVtx_X[goodMeson]")
                  .Define("goodMeson_bestVtx_Y", "K0Star_bestVtx_Y[goodMeson]")
                  .Define("goodMeson_bestVtx_Z", "K0Star_bestVtx_Z[goodMeson]")
                  .Define("goodMeson_massErr", "K0Star_kin_massErr[goodMeson]")
                  .Define("goodMeson_trk1_pt", "K0Star_pion_pt[goodMeson]")
                  .Define("goodMeson_trk2_pt", "K0Star_kaon_pt[goodMeson]")
                  .Define("goodMeson_trk1_eta", "K0Star_pion_eta[goodMeson]")
                  .Define("goodMeson_trk2_eta", "K0Star_kaon_eta[goodMeson]")
                  .Define("goodMeson_DR","DeltaR(K0Star_kaon_eta[goodMeson],K0Star_pion_eta[goodMeson],K0Star_kaon_phi[goodMeson],K0Star_pion_phi[goodMeson])")
                  .Define("wrongMeson","({}".format(GOODRHO)+")")
                  .Define("wrongMeson_pt","Sum(wrongMeson) > 0 ? rho_kin_pt[wrongMeson]: ROOT::VecOps::RVec<float>(0.f)")
                  .Define("wrongMeson2","({}".format(GOODPHI)+")")
                  .Define("wrongMeson2_pt","Sum(wrongMeson2) > 0 ? phi_kin_pt[wrongMeson2]: ROOT::VecOps::RVec<float>(0.f)")
                  )

    if(isOmega3PiCat=="true"):

        print("OMEGA = ", GOODOMEGA3Pi)
        dfbase = (df.Filter("nomega>0")
                  .Define("goodMeson_corr_mass", "Sum(omega_Nphotons)==1 ? sum2BodyMass(omega_kin_pt, omega_kin_eta, omega_kin_phi, omega_kin_mass, omega_photon_pt, omega_photon_eta, omega_photon_phi) : omega_Nbody_mass")
                  .Define("goodMeson","({}".format(GOODOMEGA3Pi)+" && abs(goodMeson_corr_mass-0.782)<0.1 && {}".format(isOmega3PiCat)+")")
                  .Filter("Sum(goodMeson)>0", "one good Omega (ptPhi, validfit, ptTracks)")
                  .Define("goodMeson_mass", "goodMeson_corr_mass[goodMeson]")
                  .Define("goodMeson_mass_RAW", "omega_Nbody_mass[goodMeson]")
                  .Define("goodMeson_pt", "omega_Nbody_pt[goodMeson]")
                  .Define("goodMeson_eta", "omega_Nbody_eta[goodMeson]")
                  .Define("goodMeson_phi", "omega_Nbody_phi[goodMeson]")
                  .Define("goodMeson_iso", "omega_iso[goodMeson]")
                  .Define("goodMeson_vtx_chi2dof", "omega_kin_vtx_chi2dof[goodMeson]")
                  .Define("goodMeson_vtx_prob", "omega_kin_vtx_prob[goodMeson]")
                  .Define("goodMeson_sipPV", "omega_kin_sipPV[goodMeson]")
                  .Define("goodMeson_bestVtx_idx", "omega_bestVtx_idx[goodMeson]")
                  .Define("goodMeson_bestVtx_X", "omega_bestVtx_X[goodMeson]")
                  .Define("goodMeson_bestVtx_Y", "omega_bestVtx_Y[goodMeson]")
                  .Define("goodMeson_bestVtx_Z", "omega_bestVtx_Z[goodMeson]")
                  .Define("goodMeson_massErr", "omega_kin_massErr[goodMeson]")
                  .Define("goodMeson_charged_mass", "omega_kin_mass[goodMeson]")
                  .Define("goodMeson_charged_pt", "omega_kin_pt[goodMeson]")
                  .Define("goodMeson_charged_eta", "omega_kin_eta[goodMeson]")
                  .Define("goodMeson_charged_phi", "omega_kin_phi[goodMeson]")
                  .Define("goodMeson_trk1_pt", "omega_trk1_pt[goodMeson]")
                  .Define("goodMeson_trk2_pt", "omega_trk2_pt[goodMeson]")
                  .Define("goodMeson_trk1_eta", "omega_trk1_eta[goodMeson]")
                  .Define("goodMeson_trk2_eta", "omega_trk2_eta[goodMeson]")
                  .Define("goodMeson_DR","DeltaR(omega_trk1_eta[goodMeson],omega_trk2_eta[goodMeson],omega_trk1_phi[goodMeson],omega_trk2_phi[goodMeson])")
                  .Define("goodMeson_photon_pt", "(Sum(omega_Nphotons[goodMeson]) >= 1) ? omega_photon_pt[goodMeson] : Vec_f {0.0}")
                  .Define("goodMeson_photon2_pt", "(Sum(omega_Nphotons[goodMeson]) >= 2) ? omega_photon2_pt[goodMeson] : Vec_f {0.0}")
                  .Define("goodMeson_DR_trkPairPh","(Sum(omega_Nphotons[goodMeson]) >= 1) ? DeltaR(omega_kin_eta[goodMeson],omega_photon_eta[goodMeson],omega_kin_phi[goodMeson],omega_photon_phi[goodMeson]) : Vec_f {0.0}")
                  .Define("goodMeson_DR_trkPairPh2","(Sum(omega_Nphotons[goodMeson]) >= 2) ? DeltaR(omega_kin_eta[goodMeson],omega_photon2_eta[goodMeson],omega_kin_phi[goodMeson],omega_photon2_phi[goodMeson]) : Vec_f {0.0}")
                  .Define("goodMeson_dEta_trkPairPh","abs(omega_kin_eta[goodMeson] - omega_photon_eta[goodMeson])")
#                  .Define("wrongMeson","({}".format(GOODRHO)+")")
#                  .Define("wrongMeson_pt","Sum(wrongMeson) > 0 ? rho_kin_pt[wrongMeson]: ROOT::VecOps::RVec<float>(0.f)")
                  .Define("wrongMeson_pt", "ROOT::VecOps::RVec<float>(0.f)" )
                  .Define("wrongMeson2_pt", "ROOT::VecOps::RVec<float>(0.f)" )
                  )

    if(isPhi3PiCat=="true"):

        print("Phi3Pi = ", GOODPHI3Pi)
        dfbase = (df.Filter("nomega>0")
                  .Define("goodMeson_corr_mass", "Sum(omega_Nphotons)==1 ? sum2BodyMass(omega_kin_pt, omega_kin_eta, omega_kin_phi, omega_kin_mass, omega_photon_pt, omega_photon_eta, omega_photon_phi) : omega_Nbody_mass")
                  .Define("goodMeson","({}".format(GOODPHI3Pi)+" && abs(goodMeson_corr_mass-1.02)<0.1 && {}".format(isPhi3PiCat)+")")
                  .Filter("Sum(goodMeson)>0", "one good Omega mass)")
                  .Define("goodMeson_mass", "goodMeson_corr_mass[goodMeson]")
                  .Define("goodMeson_mass_RAW", "omega_Nbody_mass[goodMeson]")
                  .Define("goodMeson_pt", "omega_Nbody_pt[goodMeson]")
                  .Define("goodMeson_eta", "omega_Nbody_eta[goodMeson]")
                  .Define("goodMeson_phi", "omega_Nbody_phi[goodMeson]")
                  #
                  .Define("goodMeson_iso", "omega_iso[goodMeson]")
                  .Define("goodMeson_vtx_chi2dof", "omega_kin_vtx_chi2dof[goodMeson]")
                  .Define("goodMeson_vtx_prob", "omega_kin_vtx_prob[goodMeson]")
                  .Define("goodMeson_sipPV", "omega_kin_sipPV[goodMeson]")
                  .Define("goodMeson_bestVtx_idx", "omega_bestVtx_idx[goodMeson]")
                  .Define("goodMeson_bestVtx_X", "omega_bestVtx_X[goodMeson]")
                  .Define("goodMeson_bestVtx_Y", "omega_bestVtx_Y[goodMeson]")
                  .Define("goodMeson_bestVtx_Z", "omega_bestVtx_Z[goodMeson]")
                  .Define("goodMeson_massErr", "omega_kin_massErr[goodMeson]")
                  .Define("goodMeson_charged_mass", "omega_kin_mass[goodMeson]")
                  .Define("goodMeson_charged_pt", "omega_kin_pt[goodMeson]")
                  .Define("goodMeson_charged_eta", "omega_kin_eta[goodMeson]")
                  .Define("goodMeson_charged_phi", "omega_kin_phi[goodMeson]")
                  .Define("goodMeson_trk1_pt", "omega_trk1_pt[goodMeson]")
                  .Define("goodMeson_trk2_pt", "omega_trk2_pt[goodMeson]")
                  .Define("goodMeson_trk1_eta", "omega_trk1_eta[goodMeson]")
                  .Define("goodMeson_trk2_eta", "omega_trk2_eta[goodMeson]")
                  .Define("goodMeson_DR","DeltaR(omega_trk1_eta[goodMeson],omega_trk2_eta[goodMeson],omega_trk1_phi[goodMeson],omega_trk2_phi[goodMeson])")
                  .Define("goodMeson_photon_pt", "(Sum(omega_Nphotons[goodMeson]) >= 1) ? omega_photon_pt[goodMeson] : Vec_f {0.0}")
                  .Define("goodMeson_photon2_pt", "(Sum(omega_Nphotons[goodMeson]) >= 2) ? omega_photon2_pt[goodMeson] : Vec_f {0.0}")
                  .Define("goodMeson_DR_trkPairPh","(Sum(omega_Nphotons[goodMeson]) >= 1) ? DeltaR(omega_kin_eta[goodMeson],omega_photon_eta[goodMeson],omega_kin_phi[goodMeson],omega_photon_phi[goodMeson]) : Vec_f {0.0}")
                  .Define("goodMeson_DR_trkPairPh2","(Sum(omega_Nphotons[goodMeson]) >= 2) ? DeltaR(omega_kin_eta[goodMeson],omega_photon2_eta[goodMeson],omega_kin_phi[goodMeson],omega_photon2_phi[goodMeson]) : Vec_f {0.0}")
                  .Define("goodMeson_dEta_trkPairPh","abs(omega_kin_eta[goodMeson] - omega_photon_eta[goodMeson])")
#                  .Define("wrongMeson","({}".format(GOODRHO)+")")
#                  .Define("wrongMeson_pt","Sum(wrongMeson) > 0 ? rho_kin_pt[wrongMeson]: ROOT::VecOps::RVec<float>(0.f)")
                  .Define("wrongMeson_pt", "ROOT::VecOps::RVec<float>(0.f)" )
                  .Define("wrongMeson2_pt", "ROOT::VecOps::RVec<float>(0.f)" )
                  )

    if(isD0StarCat=="true"):

        print("D0Star = ", GOODD0STAR)
        dfbase = (df.Filter("nd0").Define("goodMeson","({}".format(GOODD0STAR)+" && {}".format(isD0StarCat)+")")
                  .Filter("Sum(goodMeson)>0", "one good D0Star (ptPhi, validfit, ptTracks)")
                  .Define("goodMeson_mass", "d0_d0Star_Nbody_mass[goodMeson]")
                  .Define("goodMeson_pt", "d0_d0Star_Nbody_pt[goodMeson]")
                  .Define("goodMeson_eta", "d0_d0Star_Nbody_eta[goodMeson]")
                  .Define("goodMeson_phi", "d0_d0Star_Nbody_phi[goodMeson]")
                  #
                  .Define("goodMeson_iso", "d0_iso[goodMeson]")
                  .Define("goodMeson_vtx_chi2dof", "d0_kin_vtx_chi2dof[goodMeson]")
                  .Define("goodMeson_vtx_prob", "d0_kin_vtx_prob[goodMeson]")
                  .Define("goodMeson_sipPV", "d0_kin_sipPV[goodMeson]")
                  .Define("goodMeson_bestVtx_idx", "d0_bestVtx_idx[goodMeson]")
                  .Define("goodMeson_bestVtx_X", "d0_bestVtx_X[goodMeson]")
                  .Define("goodMeson_bestVtx_Y", "d0_bestVtx_Y[goodMeson]")
                  .Define("goodMeson_bestVtx_Z", "d0_bestVtx_Z[goodMeson]")
                  .Define("goodMeson_massErr", "d0_kin_massErr[goodMeson]")
                  .Define("goodMeson_charged_mass", "d0_kin_mass[goodMeson]")
                  .Define("goodMeson_charged_pt", "d0_kin_pt[goodMeson]")
                  .Define("goodMeson_charged_eta", "d0_kin_eta[goodMeson]")
                  .Define("goodMeson_charged_phi", "d0_kin_phi[goodMeson]")
                  .Define("goodMeson_trk1_pt", "d0_pion_pt[goodMeson]")
                  .Define("goodMeson_trk2_pt", "d0_kaon_pt[goodMeson]")
                  .Define("goodMeson_trk1_eta", "d0_pion_eta[goodMeson]")
                  .Define("goodMeson_trk2_eta", "d0_kaon_eta[goodMeson]")
                  .Define("goodMeson_DR","DeltaR(d0_pion_eta[goodMeson],d0_kaon_eta[goodMeson],d0_pion_phi[goodMeson],d0_kaon_phi[goodMeson])")
                  .Define("goodMeson_photon_pt", "(Sum(d0_d0Star_Nphotons[goodMeson]) >= 1) ? d0_d0Star_photon_pt[goodMeson] : Vec_f {0.0}")
                  .Define("goodMeson_photon2_pt", "(Sum(d0_d0Star_Nphotons[goodMeson]) >= 2) ? d0_d0Star_photon2_pt[goodMeson] : Vec_f {0.0}")
                  .Define("goodMeson_DR_trkPairPh","(Sum(d0_d0Star_Nphotons[goodMeson]) >= 1) ? DeltaR(d0_kin_eta[goodMeson],d0_d0Star_photon_eta[goodMeson],d0_kin_phi[goodMeson],d0_d0Star_photon_phi[goodMeson]) : Vec_f {0.0}")
                  .Define("goodMeson_DR_trkPairPh2","(Sum(d0_d0Star_Nphotons[goodMeson]) >= 2) ? DeltaR(d0_kin_eta[goodMeson],d0_d0Star_photon2_eta[goodMeson],d0_kin_phi[goodMeson],d0_d0Star_photon2_phi[goodMeson]) : Vec_f {0.0}")
                  .Define("goodMeson_dEta_trkPairPh","abs(d0_kin_eta[goodMeson] - d0_d0Star_photon_eta[goodMeson])")
#                  .Define("wrongMeson","({}".format(GOODRHO)+")")
#                  .Define("wrongMeson_pt","Sum(wrongMeson) > 0 ? rho_kin_pt[wrongMeson]: ROOT::VecOps::RVec<float>(0.f)")
                  .Define("wrongMeson_pt", "ROOT::VecOps::RVec<float>(0.f)" )
                  .Define("wrongMeson2_pt", "ROOT::VecOps::RVec<float>(0.f)" )
                  )

    if(isD0Pi0StarCat=="true"):

        print("D0Star = ", GOODD0PI0STAR)
        dfbase = (df.Filter("nd0pi0")
                  .Define("goodMeson_corr_mass", "Sum(d0pi0_d0Star_Nphotons)>0 ? sum2BodyMass(d0pi0_kin_pt, d0pi0_kin_eta, d0pi0_kin_phi, d0pi0_kin_mass, d0pi0_d0Star_photon_pt, d0pi0_d0Star_photon_eta, d0pi0_d0Star_photon_phi) : d0pi0_d0Star_Nbody_mass")
                  .Define("goodMeson","({}".format(GOODD0PI0STAR)+" && {}".format(isD0Pi0StarCat)+")")
                  .Filter("Sum(goodMeson)>0", "one good D0Star (ptPhi, validfit, ptTracks)")
#                  .Define("goodMeson_mass", "Sum(d0pi0_d0Star_Nphotons[goodMeson])==1 ? sum2BodyMass(d0pi0_d0Star_kin_pt[goodMeson], d0pi0_d0Star_kin_eta[goodMeson], d0pi0_d0Star_kin_phi[goodMeson], d0pi0_d0Star_kin_mass[goodMeson], d0pi0_d0Star_photon_pt[goodMeson], d0pi0_d0Star_photon_eta[goodMeson], d0pi0_d0Star_photon_phi[goodMeson]) : d0pi0_d0Star_Nbody_mass[goodMeson]")
                  .Define("goodMeson_mass", "goodMeson_corr_mass[goodMeson]")
                  .Define("goodMeson_mass_RAW", "d0pi0_d0Star_Nbody_mass[goodMeson]")
                  .Define("goodMeson_pt", "d0pi0_d0Star_Nbody_pt[goodMeson]")
                  .Define("goodMeson_eta", "d0pi0_d0Star_Nbody_eta[goodMeson]")
                  .Define("goodMeson_phi", "d0pi0_d0Star_Nbody_phi[goodMeson]")
                  #
                  .Define("goodMeson_iso", "d0pi0_iso[goodMeson]")
                  .Define("goodMeson_vtx_chi2dof", "d0pi0_kin_vtx_chi2dof[goodMeson]")
                  .Define("goodMeson_vtx_prob", "d0pi0_kin_vtx_prob[goodMeson]")
                  .Define("goodMeson_sipPV", "d0pi0_kin_sipPV[goodMeson]")
                  .Define("goodMeson_bestVtx_idx", "d0pi0_bestVtx_idx[goodMeson]")
                  .Define("goodMeson_bestVtx_X", "d0pi0_bestVtx_X[goodMeson]")
                  .Define("goodMeson_bestVtx_Y", "d0pi0_bestVtx_Y[goodMeson]")
                  .Define("goodMeson_bestVtx_Z", "d0pi0_bestVtx_Z[goodMeson]")
                  .Define("goodMeson_massErr", "d0pi0_kin_massErr[goodMeson]")
                  .Define("goodMeson_charged_mass", "d0pi0_kin_mass[goodMeson]")
                  .Define("goodMeson_charged_pt", "d0pi0_kin_pt[goodMeson]")
                  .Define("goodMeson_charged_eta", "d0pi0_kin_eta[goodMeson]")
                  .Define("goodMeson_charged_phi", "d0pi0_kin_phi[goodMeson]")
                  .Define("goodMeson_trk1_pt", "d0pi0_pion_pt[goodMeson]")
                  .Define("goodMeson_trk2_pt", "d0pi0_kaon_pt[goodMeson]")
                  .Define("goodMeson_trk1_eta", "d0pi0_pion_eta[goodMeson]")
                  .Define("goodMeson_trk2_eta", "d0pi0_kaon_eta[goodMeson]")
                  .Define("goodMeson_DR","DeltaR(d0pi0_pion_eta[goodMeson],d0pi0_kaon_eta[goodMeson],d0pi0_pion_phi[goodMeson],d0pi0_kaon_phi[goodMeson])")
                  .Define("goodMeson_photon_pt", "(Sum(d0pi0_d0Star_Nphotons[goodMeson]) >= 1) ? d0pi0_d0Star_photon_pt[goodMeson] : Vec_f {0.0}")
                  .Define("goodMeson_photon2_pt", "(Sum(d0pi0_d0Star_Nphotons[goodMeson]) >= 2) ? d0pi0_d0Star_photon2_pt[goodMeson] : Vec_f {0.0}")
                  .Define("goodMeson_DR_trkPairPh","(Sum(d0pi0_d0Star_Nphotons[goodMeson]) >= 1) ? DeltaR(d0pi0_kin_eta[goodMeson],d0pi0_d0Star_photon_eta[goodMeson],d0pi0_kin_phi[goodMeson],d0pi0_d0Star_photon_phi[goodMeson]) : Vec_f {0.0}")
                  .Define("goodMeson_DR_trkPairPh2","(Sum(d0pi0_d0Star_Nphotons[goodMeson]) >= 2) ? DeltaR(d0pi0_kin_eta[goodMeson],d0pi0_d0Star_photon2_eta[goodMeson],d0pi0_kin_phi[goodMeson],d0pi0_d0Star_photon2_phi[goodMeson]) : Vec_f {0.0}")
                  .Define("goodMeson_dEta_trkPairPh","abs(d0pi0_kin_eta[goodMeson] - d0pi0_d0Star_photon_eta[goodMeson])")
#                  .Define("wrongMeson","({}".format(GOODRHO)+")")
#                  .Define("wrongMeson_pt","Sum(wrongMeson) > 0 ? rho_kin_pt[wrongMeson]: ROOT::VecOps::RVec<float>(0.f)")
                  .Define("wrongMeson_pt", "ROOT::VecOps::RVec<float>(0.f)" )
                  .Define("wrongMeson2_pt", "ROOT::VecOps::RVec<float>(0.f)" )
                  )

    genMatchPDFNum='-1'
    if isPhiCat=="true": genMatchPDFNum='333'
    if isRhoCat=="true": genMatchPDFNum='113'
    if isK0StarCat=="true": genMatchPDFNum='313'
    if isOmega3PiCat=="true": genMatchPDFNum='223'
    if isPhi3PiCat=="true": genMatchPDFNum='333'

    dfHig = (dfbase.Define("index_pair","HiggsCandFromRECO(goodMeson_pt, goodMeson_eta, goodMeson_phi, goodMeson_mass, goodMeson_trk1_pt, goodMeson_trk2_pt, wrongMeson_pt, wrongMeson2_pt, goodPhotons_pt, goodPhotons_eta, goodPhotons_phi)").Filter("index_pair[0]!= -1", "at least a good meson candidate")
             #               .Define("jet_mask", "cleaningMask(Photon_jetIdx[goodPhotons],nJet)")
             .Define("jet_mask2", "cleaningJetFromOBJ(Jet_eta, Jet_phi, goodMeson_eta[index_pair[0]], goodMeson_phi[index_pair[0]])")
             #.Define("meson_pt", "(index_pair[0]!= -1) ? goodMeson_pt[index_pair[0]]: 0.f")
             .Define("photon_pt", "(index_pair[1]!= -1) ? goodPhotons_pt[index_pair[1]]: 0.f")
             .Define("meson_isoNeuHad", "(index_pair[0]!= -1) ? goodMeson_isoNeuHad[index_pair[0]]: 0.f")
             )

    if(isOmega3PiCat=="true" or isPhi3PiCat=="true" or isD0StarCat=="true" or isD0Pi0StarCat=="true"): dfHig_ = callMVAregress(dfHig,isVBF,isVBFlow,isGF,isZinv)
    if(isRhoCat=="true" or isPhiCat=="true" or isK0StarCat=="true"): dfHig_ = dfHig.Define("meson_pt", "(index_pair[0]!= -1) ? goodMeson_pt[index_pair[0]]: 0.f")

    dfFinal = (dfHig_.Define("HCandMass", "compute_HiggsVars_var(meson_pt,goodMeson_eta[index_pair[0]],goodMeson_phi[index_pair[0]],goodMeson_mass[index_pair[0]],photon_pt,goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]],0)")
               .Define("HCandPT",   "compute_HiggsVars_var(meson_pt,goodMeson_eta[index_pair[0]],goodMeson_phi[index_pair[0]],goodMeson_mass[index_pair[0]],photon_pt,goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]],1)")
               .Define("HCandPHI",   "compute_HiggsVars_var(meson_pt,goodMeson_eta[index_pair[0]],goodMeson_phi[index_pair[0]],goodMeson_mass[index_pair[0]],photon_pt,goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]],2)")
               .Define("HCandMass_noregress", "compute_HiggsVars_var(goodMeson_pt[index_pair[0]],goodMeson_eta[index_pair[0]],goodMeson_phi[index_pair[0]],goodMeson_mass[index_pair[0]],photon_pt,goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]],0)")
#               .Define("HCandMass_VTXcorr", "compute_HiggsVars_var_VtxCorr(meson_pt,goodMeson_eta[index_pair[0]],goodMeson_phi[index_pair[0]],goodMeson_mass[index_pair[0]],goodMeson_bestVtx_X[index_pair[0]],goodMeson_bestVtx_Y[index_pair[0]], goodMeson_bestVtx_Z[index_pair[0]],photon_pt,goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]],goodPhotons_x_calo[index_pair[1]],goodPhotons_y_calo[index_pair[1]],goodPhotons_z_calo[index_pair[1]],0)")
               .Define("dPhiGammaMesonCand","abs(deltaPhi(goodPhotons_phi[index_pair[1]], goodMeson_phi[index_pair[0]]))")
               .Define("dEtaGammaMesonCand","abs(goodPhotons_eta[index_pair[1]] - goodMeson_eta[index_pair[0]])")
               .Define("dRGammaMesonCand","deltaR(goodPhotons_eta[index_pair[1]], goodPhotons_phi[index_pair[1]], goodMeson_eta[index_pair[0]], goodMeson_phi[index_pair[0]])")
               .Define("sigmaHCandMass_Rel2","(goodPhotons_energyErr[index_pair[1]]*goodPhotons_energyErr[index_pair[1]])/(goodPhotons_pt[index_pair[1]]*goodPhotons_pt[index_pair[1]]) + (goodMeson_massErr[index_pair[0]]*goodMeson_massErr[index_pair[0]])/(goodMeson_mass[index_pair[0]]*goodMeson_mass[index_pair[0]])")
               .Define("classify","topology(goodPhotons_eta[index_pair[1]], goodMeson_eta[index_pair[0]])")
               )

    if (isData == "false"):
        dfFinal = dfFinal.Define("indexMatch","genMatch(goodMeson_pt[index_pair[0]],goodMeson_eta[index_pair[0]],goodMeson_phi[index_pair[0]],goodMeson_mass[index_pair[0]],GenPart_eta,GenPart_phi,GenPart_pdgId,GenPart_genPartIdxMother,{0})".format(genMatchPDFNum))

    return dfFinal

def dfwithSYST(df,year):


    if df.HasColumn("Photon_dEsigmaUp") and df.HasColumn("Photon_dEsigmaDown"):

        dfVaryPh=(df
                  # ONLY FOR MC
                  #               .DefinePerSample("goodPhotons_dEsigmaUp", 'rdfsampleinfo_.Contains("+Run") ? 0.f : Photon_dEsigmaUp[goodPhotons]')
                  #               .DefinePerSample("goodPhotons_dEsigmaDown", 'rdfsampleinfo_.Contains("+Run") ? 0.f : Photon_dEsigmaDown[goodPhotons]')
                  #               .DefinePerSample("goodPhotons_dEsigmaUp", "mc>=0 ? Photon_dEsigmaUp[goodPhotons] : Photon_dEscaleUp[goodPhotons]")
                  #               .DefinePerSample("goodPhotons_dEsigmaDown", "mc>=0 ? Photon_dEsigmaDown[goodPhotons] : Photon_dEscaleDown[goodPhotons]")
                  #               .DefinePerSample("photon_dEsigmaUp",'rdfsampleinfo_.Contains("+Run") ? 0.0f : (1.f+Photon_dEsigmaUp[goodPhotons][index_pair[1]])')
                  #               .DefinePerSample("photon_dEsigmaDown",'rdfsampleinfo_.Contains("+Run") ? 0.0f : (1.f+Photon_dEsigmaDown[goodPhotons][index_pair[1]])')
                  #               .Define("photon_dEsigmaUp",'rdfsampleinfo_.Contains("+Run") ? 0.0f : (1.f+goodPhotons_dEsigmaUp[index_pair[1]])')
                  #               .Define("photon_dEsigmaDown",'rdfsampleinfo_.Contains("+Run") ? 0.0f : (1.f+goodPhotons_dEsigmaDown[index_pair[1]])')
                  # ONLY FOR MC
                  .Define("goodPhotons_dEsigmaUp", 'Photon_dEsigmaUp[goodPhotons]')
                  .Define("goodPhotons_dEsigmaDown", 'Photon_dEsigmaDown[goodPhotons]')
                  .Define("photon_dEsigmaUp",'(1.f+goodPhotons_dEsigmaUp[index_pair[1]])')
                  .Define("photon_dEsigmaDown",'(1.f+goodPhotons_dEsigmaDown[index_pair[1]])')
                  .Vary("photon_pt", "ROOT::RVecF{photon_pt*photon_dEsigmaDown,photon_pt*photon_dEsigmaUp}", variationTags=["dn","up"], variationName = "PhotonSYST")
                  ##
                  )
        dfVary = dfVaryPh
    else:
        dfVary = df

    photonIDyear=year
    if year==12016: photonIDyear = '2016preVFP'
    if year==22016: photonIDyear = '2016postVFP'

    photonIDwp="wp90"
    if isGF or isVBFlow: photonIDwp="wp80"

    muonIDyear=year
    if year==2018: muonIDyear = '2018_UL'
    if year==2017: muonIDyear = '2017_UL'
    if year==12016: muonIDyear = '2016preVFP_UL'
    if year==22016: muonIDyear = '2016postVFP_UL'

    ###"(index_pair[1]!= -1 && goodPhotons_isScEtaEB[index_pair[1]]>0) ? wp90 : wp80"

    dfFinal_withSF = (dfVary
                      .Define("SFphoton_ID_Nom",'corr_sf.eval_photonSF("{0}", "sf", "{1}", goodPhotons_eta[index_pair[1]], photon_pt)'.format(photonIDyear,photonIDwp))
                      .Define("SFphoton_ID_Up",'corr_sf.eval_photonSF("{0}", "sfup", "{1}" , goodPhotons_eta[index_pair[1]], photon_pt)'.format(photonIDyear,photonIDwp))
                      .Define("SFphoton_ID_Dn",'corr_sf.eval_photonSF("{0}", "sfdown", "{1}" , goodPhotons_eta[index_pair[1]], photon_pt)'.format(photonIDyear,photonIDwp))
                      #                      .Define("SFphoton_PixVeto_Nom",'corr_sf.eval_photonPixVetoSF("{0}", "sf", "{1}", goodPhotons_eta[index_pair[1]], photon_pt)'.format(photonIDyear,"MVA"))
                      #                      .Define("SFphoton_PixVeto_Up",'corr_sf.eval_photonPixVetoSF("{0}", "sfup", "{1}" , goodPhotons_eta[index_pair[1]], photon_pt)'.format(photonIDyear,"MVA"))
                      #                      .Define("SFphoton_PixVeto_Dn",'corr_sf.eval_photonPixVetoSF("{0}", "sfdown", "{1}" , goodPhotons_eta[index_pair[1]], photon_pt)'.format(photonIDyear,"MVA"))
                      .Define("SFphoton_PixVeto_Nom","1.f")
                      ##
                      ## c=a*b Dc = aDa + b*db + DaDb
                      .Define("SFmeson1_reco_Nom",'corr_sf.eval_muonTRKSF("{0}", "sf", goodMeson_trk1_eta[index_pair[0]], goodMeson_trk1_pt[index_pair[0]])'.format(muonIDyear))
                      .Define("meson1_reco_DeltaUp",'abs(SFmeson1_reco_Nom-corr_sf.eval_muonTRKSF("{0}", "systup", goodMeson_trk1_eta[index_pair[0]], goodMeson_trk1_pt[index_pair[0]]))'.format(muonIDyear))
                      .Define("meson1_reco_DeltaDn",'abs(SFmeson1_reco_Nom-corr_sf.eval_muonTRKSF("{0}", "systdown", goodMeson_trk1_eta[index_pair[0]], goodMeson_trk1_pt[index_pair[0]]))'.format(muonIDyear))
                      .Define("SFmeson2_reco_Nom",'corr_sf.eval_muonTRKSF("{0}", "sf", goodMeson_trk2_eta[index_pair[0]], goodMeson_trk2_pt[index_pair[0]])'.format(muonIDyear))
                      .Define("meson2_reco_DeltaUp",'abs(SFmeson2_reco_Nom-corr_sf.eval_muonTRKSF("{0}", "systup", goodMeson_trk2_eta[index_pair[0]], goodMeson_trk2_pt[index_pair[0]]))'.format(muonIDyear))
                      .Define("meson2_reco_DeltaDn",'abs(SFmeson2_reco_Nom-corr_sf.eval_muonTRKSF("{0}", "systdown", goodMeson_trk2_eta[index_pair[0]], goodMeson_trk2_pt[index_pair[0]]))'.format(muonIDyear))
                      .Define("SFmeson_reco_Nom","SFmeson1_reco_Nom*SFmeson2_reco_Nom")
                      .Define("SFmeson_reco_Up","SFmeson_reco_Nom+SFmeson1_reco_Nom*meson2_reco_DeltaUp+SFmeson2_reco_Nom*meson1_reco_DeltaUp")
                      .Define("SFmeson_reco_Dn","SFmeson_reco_Nom-SFmeson1_reco_Nom*meson2_reco_DeltaDn-SFmeson2_reco_Nom*meson1_reco_DeltaDn")
                      ##
                      .Define("SFpu_Nom",'corr_sf.eval_puSF(Pileup_nTrueInt,"nominal")')
                      .Define("SFpu_Up",'corr_sf.eval_puSF(Pileup_nTrueInt,"up")')
                      .Define("SFpu_Dn",'corr_sf.eval_puSF(Pileup_nTrueInt,"down")')
                      #
                      .Define("idx_nom_up_down", "indices(3)")
                      #
                      )


    if isGF or isVBFlow:
        dfFinal_withSF_2 = (dfFinal_withSF
                            .Define("photonTrigSF_Nom",'SF_HLT_leg(photon_pt,22,0)')
                            .Define("photonTrigSF_Up",'SF_HLT_leg(photon_pt,22,1)')
                            .Define("photonTrigSF_Dn",'SF_HLT_leg(photon_pt,22,-1)')
                            .Define("tauTrigSF_Nom",'SF_HLT_leg(meson_pt,15,0)')
                            .Define("tauTrigSF_Up",'SF_HLT_leg(meson_pt,15,1)')
                            .Define("tauTrigSF_Dn",'SF_HLT_leg(meson_pt,15,-1)')
                            .Define("chIso_Nom",'SF_chIso(meson_pt,goodMeson_eta[index_pair[0]],0)')
                            .Define("chIso_Up",'SF_chIso(meson_pt,goodMeson_eta[index_pair[0]],1)')
                            .Define("chIso_Dn",'SF_chIso(meson_pt,goodMeson_eta[index_pair[0]],-1)')
                            )

    if isZ or isW:
        dfFinal_withSF_2 = (dfFinal_withSF
                            .Define("SFmuon_ISO_Nom",'isMuorEle==1 ? corr_sf.eval_muonISOSF("{0}", "sf", LeadingLeptonEta, LeadingLeptonPt, "T"):1'.format(muonIDyear))
                            .Define("SFmuon_ISO_Up",'isMuorEle==1 ? corr_sf.eval_muonISOSF("{0}", "systup", LeadingLeptonEta, LeadingLeptonPt, "T"):1'.format(muonIDyear))
                            .Define("SFmuon_ISO_Dn",'isMuorEle==1 ? corr_sf.eval_muonISOSF("{0}", "systdown", LeadingLeptonEta, LeadingLeptonPt, "T"):1'.format(muonIDyear))
                            .Define("SFmuon_ID_Nom",'isMuorEle==1 ? corr_sf.eval_muonIDSF("{0}", "sf", LeadingLeptonEta, LeadingLeptonPt, "T"):1'.format(muonIDyear))
                            .Define("SFmuon_ID_Up",'isMuorEle==1 ? corr_sf.eval_muonIDSF("{0}", "systup", LeadingLeptonEta, LeadingLeptonPt, "T"):1'.format(muonIDyear))
                            .Define("SFmuon_ID_Dn",'isMuorEle==1 ? corr_sf.eval_muonIDSF("{0}", "systdown", LeadingLeptonEta, LeadingLeptonPt, "T"):1'.format(muonIDyear))
                            #
                            .Define("SFelectron_ID_Nom",'isMuorEle==2 ? corr_sf.eval_electronSF("{0}", "sf", "{1}", LeadingLeptonEta, LeadingLeptonPt):1'.format(photonIDyear,"wp80iso"))
                            .Define("SFelectron_ID_Up",'isMuorEle==2 ? corr_sf.eval_electronSF("{0}", "sfup", "{1}", LeadingLeptonEta, LeadingLeptonPt):1'.format(photonIDyear,"wp80iso"))
                            .Define("SFelectron_ID_Dn",'isMuorEle==2 ? corr_sf.eval_electronSF("{0}", "sfdown", "{1}", LeadingLeptonEta, LeadingLeptonPt):1'.format(photonIDyear,"wp80iso"))
                            .Define("SFelectron_reco_Nom",'isMuorEle==2 ? corr_sf.eval_electronSF("{0}", "sf", "{1}", LeadingLeptonEta, LeadingLeptonPt):1'.format(photonIDyear,"RecoAbove20"))
                            .Define("SFelectron_reco_Up",'isMuorEle==2 ? corr_sf.eval_electronSF("{0}", "sfup", "{1}", LeadingLeptonEta, LeadingLeptonPt):1'.format(photonIDyear,"RecoAbove20"))
                            .Define("SFelectron_reco_Dn",'isMuorEle==2 ? corr_sf.eval_electronSF("{0}", "sfdown", "{1}", LeadingLeptonEta, LeadingLeptonPt):1'.format(photonIDyear,"RecoAbove20"))
                            #
                            )


    if isZ:
        dfFinal_withSF_3 = (dfFinal_withSF_2
                            .Define("SFmuon2_ISO_Nom",'isMuorEle==1 ? corr_sf.eval_muonISOSF("{0}", "sf", SubLeadingLeptonEta, SubLeadingLeptonPt, "L"):1'.format(muonIDyear))
                            .Define("SFmuon2_ISO_Up",'isMuorEle==1 ? corr_sf.eval_muonISOSF("{0}", "systup", SubLeadingLeptonEta, SubLeadingLeptonPt, "L"):1'.format(muonIDyear))
                            .Define("SFmuon2_ISO_Dn",'isMuorEle==1 ? corr_sf.eval_muonISOSF("{0}", "systdown", SubLeadingLeptonEta, SubLeadingLeptonPt, "L"):1'.format(muonIDyear))
                            .Define("SFmuon2_ID_Nom",'isMuorEle==1 ? corr_sf.eval_muonIDSF("{0}", "sf", SubLeadingLeptonEta, SubLeadingLeptonPt, "M"):1'.format(muonIDyear))
                            .Define("SFmuon2_ID_Up",'isMuorEle==1 ? corr_sf.eval_muonIDSF("{0}", "systup", SubLeadingLeptonEta, SubLeadingLeptonPt, "M"):1'.format(muonIDyear))
                            .Define("SFmuon2_ID_Dn",'isMuorEle==1 ? corr_sf.eval_muonIDSF("{0}", "systdown", SubLeadingLeptonEta, SubLeadingLeptonPt, "M"):1'.format(muonIDyear))
                            #
#                            .Define("SFelectron2_ID_Nom",'isMuorEle==2 ? corr_sf.eval_electronSF("{0}", "sf", "{1}", SubLeadingLeptonEta, SubLeadingLeptonPt):1'.format(photonIDyear,"wp90iso"))
#                            .Define("SFelectron2_ID_Up",'isMuorEle==2 ? corr_sf.eval_electronSF("{0}", "sfup", "{1}", SubLeadingLeptonEta, SubLeadingLeptonPt):1'.format(photonIDyear,"wp90iso"))
#                            .Define("SFelectron2_ID_Dn",'isMuorEle==2 ? corr_sf.eval_electronSF("{0}", "sfdown", "{1}", SubLeadingLeptonEta, SubLeadingLeptonPt):1'.format(photonIDyear,"wp90iso"))
#                            .Define("SFelectron2_reco_Nom",'isMuorEle==2 ? corr_sf.eval_electronSF("{0}", "sf", "{1}", SubLeadingLeptonEta, SubLeadingLeptonPt):1'.format(photonIDyear,"RecoAbove20"))
#                            .Define("SFelectron2_reco_Up",'isMuorEle==2 ? corr_sf.eval_electronSF("{0}", "sfup", "{1}", SubLeadingLeptonEta, SubLeadingLeptonPt):1'.format(photonIDyear,"RecoAbove20"))
#                            .Define("SFelectron2_reco_Dn",'isMuorEle==2 ? corr_sf.eval_electronSF("{0}", "sfdown", "{1}", SubLeadingLeptonEta, SubLeadingLeptonPt):1'.format(photonIDyear,"RecoAbove20"))
                            #
                            )

    #####
    #####


    if isZ: dfFinal = (dfFinal_withSF_3
                       .Define("w_allSF", "w*SFpu_Nom*L1PreFiringWeight_Nom*SFphoton_ID_Nom*SFphoton_PixVeto_Nom*SFmeson_reco_Nom*SFelectron_ID_Nom*SFelectron_reco_Nom*SFmuon_ID_Nom*SFmuon_ISO_Nom*SFmuon2_ISO_Nom*SFmuon2_ID_Nom")
                       #*SFelectron2_reco_Nom*SFelectron2_ID_Nom")
                       .Define("muo2ID_weights", "NomUpDownVar(SFmuon2_ID_Nom, SFmuon2_ID_Up, SFmuon2_ID_Dn, w_allSF)")
                       .Define("muo2ISO_weights", "NomUpDownVar(SFmuon2_ISO_Nom, SFmuon2_ISO_Up, SFmuon2_ISO_Dn, w_allSF)")
#                       .Define("ele2ID_weights", "NomUpDownVar(SFelectron2_ID_Nom, SFelectron2_ID_Up, SFelectron2_ID_Dn, w_allSF)")
#                       .Define("ele2Reco_weights", "NomUpDownVar(SFelectron2_reco_Nom, SFelectron2_reco_Up, SFelectron2_reco_Dn, w_allSF)")
                       #
                       .Define("muoID_weights", "NomUpDownVar(SFmuon_ID_Nom, SFmuon_ID_Up, SFmuon_ID_Dn, w_allSF)")
                       .Define("muoISO_weights", "NomUpDownVar(SFmuon_ISO_Nom, SFmuon_ISO_Up, SFmuon_ISO_Dn, w_allSF)")
                       .Define("eleID_weights", "NomUpDownVar(SFelectron_ID_Nom, SFelectron_ID_Up, SFelectron_ID_Dn, w_allSF)")
                       .Define("eleReco_weights", "NomUpDownVar(SFelectron_reco_Nom, SFelectron_reco_Up, SFelectron_reco_Dn, w_allSF)")
                       #
                       .Define("PS_isr_weights", "NomUpDownVar(1., PSWeight[0], PSWeight[2], w_allSF)")
                       .Define("PS_fsr_weights", "NomUpDownVar(1., PSWeight[1], PSWeight[3], w_allSF)")
                       .Define("L1PreFiring_weights", "NomUpDownVar(L1PreFiringWeight_Nom, L1PreFiringWeight_Up, L1PreFiringWeight_Dn, w_allSF)")
                       .Define("pu_weights", "NomUpDownVar(SFpu_Nom, SFpu_Up, SFpu_Dn, w_allSF)")
                       .Define("phoID_weights", "NomUpDownVar(SFphoton_ID_Nom, SFphoton_ID_Up, SFphoton_ID_Dn, w_allSF)")
                       .Define("mesonRECO_weights", "NomUpDownVar(SFmeson_reco_Nom, SFmeson_reco_Up, SFmeson_reco_Dn, w_allSF)")
                       )

    elif isW: dfFinal = (dfFinal_withSF_2
                         .Define("w_allSF", "w*SFpu_Nom*L1PreFiringWeight_Nom*SFphoton_ID_Nom*SFphoton_PixVeto_Nom*SFmeson_reco_Nom*SFelectron_ID_Nom*SFelectron_reco_Nom*SFmuon_ID_Nom*SFmuon_ISO_Nom")
                         .Define("muoID_weights", "NomUpDownVar(SFmuon_ID_Nom, SFmuon_ID_Up, SFmuon_ID_Dn, w_allSF)")
                         .Define("muoISO_weights", "NomUpDownVar(SFmuon_ISO_Nom, SFmuon_ISO_Up, SFmuon_ISO_Dn, w_allSF)")
                         .Define("eleID_weights", "NomUpDownVar(SFelectron_ID_Nom, SFelectron_ID_Up, SFelectron_ID_Dn, w_allSF)")
                         .Define("eleReco_weights", "NomUpDownVar(SFelectron_reco_Nom, SFelectron_reco_Up, SFelectron_reco_Dn, w_allSF)")
                         #
                         .Define("PS_isr_weights", "NomUpDownVar(1., PSWeight[0], PSWeight[2], w_allSF)")
                         .Define("PS_fsr_weights", "NomUpDownVar(1., PSWeight[1], PSWeight[3], w_allSF)")
                         .Define("L1PreFiring_weights", "NomUpDownVar(L1PreFiringWeight_Nom, L1PreFiringWeight_Up, L1PreFiringWeight_Dn, w_allSF)")
                         .Define("pu_weights", "NomUpDownVar(SFpu_Nom, SFpu_Up, SFpu_Dn, w_allSF)")
                         .Define("phoID_weights", "NomUpDownVar(SFphoton_ID_Nom, SFphoton_ID_Up, SFphoton_ID_Dn, w_allSF)")
                         .Define("mesonRECO_weights", "NomUpDownVar(SFmeson_reco_Nom, SFmeson_reco_Up, SFmeson_reco_Dn, w_allSF)")
                         )

    elif isGF or isVBFlow: dfFinal = (dfFinal_withSF_2
                                      .Define("w_allSF", "w*SFpu_Nom*L1PreFiringWeight_Nom*SFphoton_ID_Nom*SFphoton_PixVeto_Nom*SFmeson_reco_Nom*photonTrigSF_Nom*tauTrigSF_Nom*chIso_Nom")
#                                      .Define("w_allSF", "w*SFpu_Nom*L1PreFiringWeight_Nom*SFphoton_ID_Nom*SFphoton_PixVeto_Nom*SFmeson_reco_Nom")
                                      .Define("phTrigSF_weights", "NomUpDownVar(photonTrigSF_Nom, photonTrigSF_Up, photonTrigSF_Dn, w_allSF)")
                                      .Define("tauTrigSF_weights", "NomUpDownVar(tauTrigSF_Nom, tauTrigSF_Up, tauTrigSF_Dn, w_allSF)")
                                      .Define("chIsoSF_weights", "NomUpDownVar(chIso_Nom, chIso_Up, chIso_Dn, w_allSF)")
                                      .Define("L1PreFiring_weights", "NomUpDownVar(L1PreFiringWeight_Nom, L1PreFiringWeight_Up, L1PreFiringWeight_Dn, w_allSF)")
                                      .Define("PS_isr_weights", "NomUpDownVar(1., PSWeight[0], PSWeight[2], w_allSF)")
                                      .Define("PS_fsr_weights", "NomUpDownVar(1., PSWeight[1], PSWeight[3], w_allSF)")
                                      .Define("pu_weights", "NomUpDownVar(SFpu_Nom, SFpu_Up, SFpu_Dn, w_allSF)")
                                      .Define("phoID_weights", "NomUpDownVar(SFphoton_ID_Nom, SFphoton_ID_Up, SFphoton_ID_Dn, w_allSF)")
                                      .Define("mesonRECO_weights", "NomUpDownVar(SFmeson_reco_Nom, SFmeson_reco_Up, SFmeson_reco_Dn, w_allSF)")
                                      )

    else: dfFinal = (dfFinal_withSF
                     .Define("w_allSF", "w*SFpu_Nom*L1PreFiringWeight_Nom*SFphoton_ID_Nom*SFphoton_PixVeto_Nom*SFmeson_reco_Nom")
                     .Define("L1PreFiring_weights", "NomUpDownVar(L1PreFiringWeight_Nom, L1PreFiringWeight_Up, L1PreFiringWeight_Dn, w_allSF)")
                     .Define("PS_isr_weights", "NomUpDownVar(1., PSWeight[0], PSWeight[2], w_allSF)")
                     .Define("PS_fsr_weights", "NomUpDownVar(1., PSWeight[1], PSWeight[3], w_allSF)")
                     .Define("pu_weights", "NomUpDownVar(SFpu_Nom, SFpu_Up, SFpu_Dn, w_allSF)")
                     .Define("phoID_weights", "NomUpDownVar(SFphoton_ID_Nom, SFphoton_ID_Up, SFphoton_ID_Dn, w_allSF)")
                     .Define("mesonRECO_weights", "NomUpDownVar(SFmeson_reco_Nom, SFmeson_reco_Up, SFmeson_reco_Dn, w_allSF)")
                     )

    return dfFinal

def dfCommon(df,year,isData,mc,sumw,isVBF,isVBFlow,isGF,isZinv):

    lumi = 1.
    weight = "{0}".format(1.)
    if mc>=0: weight = "{0}*genWeight*{1}".format(lumi,sumw)

    lumiIntegrated = 1.
    print('isData = ',isData)
    if (isData == "false"):
        if((isVBF or isW or isZ) and year == 2018): lumiIntegrated = lumis['2018']
        if((isW or isZ) and year == 2017): lumiIntegrated = lumis['2017']
        if((isVBF) and year == 2017): lumiIntegrated = lumis['12017']
        if((isVBF or isW or isZ) and year == 12016): lumiIntegrated = lumis['12016']
        if((isW or isZ) and year == 22016): lumiIntegrated = lumis['22016']
        if((isVBFlow or isGF or isZinv) and year == 2018): lumiIntegrated = lumis['12018']
        print('lumiIntegrated=',lumiIntegrated, ' year=',year)

    dfComm = (df
              .Define("mc","{}".format(mc))
              .Define("isData","{}".format(isData))
              .Define("applyJson","{}".format(JSON)).Filter("applyJson","pass JSON")
              .Define("w","{}".format(weight))
              .Define("wraw","{}".format(weight))
              .Define("lumiIntegrated","{}".format(lumiIntegrated))
              .Filter("PV_npvsGood>0","one good PV")
              )

    return dfComm

def callPolarization(dfcandtag, mc):


    if False:
        dfNew = (dfcandtag.Define('theta', 'getPhiPolarizationAngle(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass,event)')
                 .Define('isHPhiKaon', 'isHPhiKaon(GenPart_pdgId, GenPart_genPartIdxMother)')
                 .Define('isHPhiPhiKaon', 'isHPhiPhiKaon(GenPart_pdgId, GenPart_genPartIdxMother)')
                 .Define('fourvectors', 'getHPhiKaonFourVectors(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)')
                 .Define('ditrack_pt', 'fourvectors[0]')
                 .Define('ditrack_eta', 'fourvectors[1]')
                 .Define('ditrack_phi', 'fourvectors[2]')
                 .Define('ditrack_mass', 'fourvectors[3]')
                 .Define('genphi_pt', 'fourvectors[4]')
                 .Define('genphi_eta', 'fourvectors[5]')
                 .Define('genphi_phi', 'fourvectors[6]')
                 .Define('genphi_mass', 'fourvectors[7]')
                 )
    if (mc >= 1010 and mc <= 1038):
        dfNew = (dfcandtag.Define('theta',"getPolAngle(event,rdfslot_)")
                 )
    else:
        dfNew = dfcandtag.Define("theta","1.f")

    return dfNew

def callMVAregress(df,isVBF,isVBFlow,isGF,isZinv):

    MVAweights = ""
    if(isGF): MVAweights = "{}".format(getMVAFromJson(MVA, "isGF" , sys.argv[2] ))
    print(MVAweights)

    tmva_helper = tmva_helper_xml.TMVAHelperXML(MVAweights)
    print(tmva_helper.variables)

    dfWithMVA = (df.Define("goodMeson_photon1_pt_input_pred", "(index_pair[0]!= -1) ? goodMeson_photon_pt[index_pair[0]] : 0.f")
                 .Define("goodMeson_photon2_pt_input_pred", "(index_pair[0]!= -1) ? goodMeson_photon2_pt[index_pair[0]] : 0.f")
                 .Define("goodMeson_ditrk_mass_input_pred", "(index_pair[0]!= -1) ? goodMeson_charged_mass[index_pair[0]] : 0.f")
                 .Define("goodMeson_mass_input_pred", "(index_pair[0]!= -1) ? goodMeson_mass[index_pair[0]] : 0.f")
                 .Define("goodMeson_ditrk_eta_input_pred", "(index_pair[0]!= -1) ? goodMeson_charged_eta[index_pair[0]] : 0.f")
                 .Define("goodMeson_photon1_DR_input_pred", "(index_pair[0]!= -1) ? goodMeson_DR_trkPairPh[index_pair[0]] : 0.f")
                 .Define("goodMeson_photon2_DR_input_pred", "(index_pair[0]!= -1) ? goodMeson_DR_trkPairPh2[index_pair[0]] : 0.f")
                 .Define("goodMeson_DR_input_pred", "(index_pair[0]!= -1) ? goodMeson_DR[index_pair[0]] : 0.f")
                 .Define("delta_eta_goodMeson_ditrk_goodPhoton_input_pred", "(index_pair[0]!= -1 and index_pair[1]!= -1 ) ? abs(goodMeson_charged_eta[index_pair[0]] - goodPhotons_eta[index_pair[1]]) : 0.f")
                 .Define("delta_phi_goodMeson_ditrk_goodPhoton_input_pred", "(index_pair[0]!= -1 and index_pair[1]!= -1 ) ? abs(deltaPhi(goodMeson_charged_phi[index_pair[0]], goodPhotons_phi[index_pair[1]])) : 0.f")
                 .Define("var0_input_pred", "(goodPhotons_pt[index_pair[1]] != 0) ? goodMeson_pt[index_pair[0]]/goodPhotons_pt[index_pair[1]] : 0.f")
                 .Define("var1_input_pred", "(goodMeson_charged_pt[index_pair[1]] != 0) ? goodMeson_pt[index_pair[0]]/goodMeson_charged_pt[index_pair[1]] : 0.f")
                )

    dfWithMVA = tmva_helper.run_inference(dfWithMVA,"meson_regress_pt_sf")

    dfWithMVA = dfWithMVA.Define("meson_pt", "meson_regress_pt_sf[0]*goodMeson_pt[index_pair[0]]")

    return dfWithMVA

def callMVA(df,isVBF,isVBFlow,isGF,isZinv):

    MVAweights = ""
    if(isGF): MVAweights = "{}".format(getMVAFromJson(MVA, "isGF" , sys.argv[2] ))
    if(isVBF): MVAweights = "{}".format(getMVAFromJson(MVA, "isVBF" , sys.argv[2] ))
    if(isVBFlow): MVAweights = "{}".format(getMVAFromJson(MVA, "isVBFlow" , sys.argv[2] ))
    if(isZinv): MVAweights = "{}".format(getMVAFromJson(MVA, "isZinv" , sys.argv[2] ))
    print(MVAweights)

    tmva_helper = tmva_helper_xml.TMVAHelperXML(MVAweights)
    print(tmva_helper.variables)

    NVar = "10"
    if(isZinv): NVar = "10"
    print('NVAR=',NVar)

#    s ='''
#    TMVA::Experimental::RReader model("{0}");
#    computeModel = TMVA::Experimental::Compute<{1}, float>(model);
#    '''
#    print(s.format(MVAweights,NVar))
#    ROOT.gInterpreter.ProcessLine(s.format(MVAweights,NVar))
#    variables = ROOT.model.GetVariableNames()

    dfWithMVA = (df.Define("HCandPT__div_sqrtHCandMass", "(HCandMass>0) ? HCandPT/sqrt(HCandMass): 0.f")
               .Define("HCandPT__div_HCandMass", "(HCandMass>0) ? HCandPT/HCandMass: 0.f")
               #"goodPhotons_pt__div_HCandPT"
               .Define("photon_pt__div_HCandPT", "(index_pair[1]!= -1 && HCandPT>0 ) ? photon_pt/HCandPT: 0.f")
               .Define("photon_pt__div_HCandMass", "(index_pair[1]!= -1 && HCandMass>0) ? photon_pt/HCandMass: 0.f")
               #"goodPhotons_eta"
               .Define("photon_eta","(index_pair[1]!= -1) ? goodPhotons_eta[index_pair[1]]: 0.f")
               #"goodPhotons_mvaID"
               .Define("photon_mvaID","(index_pair[1]!= -1) ? goodPhotons_mvaID[index_pair[1]]: 0.f")
               #"goodPhotons_pfRelIso03_all"
               .Define("photon_pfRelIso03_all","(index_pair[1]!= -1) ? goodPhotons_pfRelIso03_all[index_pair[1]]: 0.f")
               #"goodPhotons_energyErr"
               .Define("photon_energyErr","(index_pair[1]!= -1) ? goodPhotons_energyErr[index_pair[1]]: 0.f")
               #
               .Define("meson_DR", "(index_pair[0]!= -1) ? goodMeson_DR[index_pair[0]]: 0.f")
               #goodMeson_DR__times_sqrtHCandMass
               .Define("meson_DR__times_sqrtHCandMass", "(index_pair[0]!= -1 && HCandMass>0) ? goodMeson_DR[index_pair[0]]*sqrt(HCandMass): 0.f")
               #"goodMeson_pt__div_HCandPT"
               .Define("meson_pt__div_HCandPT", "(index_pair[0]!= -1 && HCandPT>0) ? goodMeson_pt[index_pair[0]]/HCandPT: 0.f")
               .Define("meson_pt__div_HCandMass", "(index_pair[0]!= -1 && HCandMass>0) ? goodMeson_pt[index_pair[0]]/HCandMass: 0.f")
               #
               .Define("meson_mass","(index_pair[0]!= -1) ? goodMeson_mass[index_pair[0]]: 0.f")
               .Define("meson_massErr","(index_pair[0]!= -1) ? goodMeson_massErr[index_pair[0]]: 0.f")
               .Define("meson_iso","(index_pair[0]!= -1) ? goodMeson_iso[index_pair[0]]: 0.f")
               .Define("meson_sipPV","(index_pair[0]!= -1) ? goodMeson_sipPV[index_pair[0]]: 0.f")
               .Define("meson_trk1_eta","(index_pair[0]!= -1) ? goodMeson_trk1_eta[index_pair[0]]: 0.f")
               .Define("meson_trk2_eta","(index_pair[0]!= -1) ? goodMeson_trk2_eta[index_pair[0]]: 0.f")
               .Define("meson_trk1_pt","(index_pair[0]!= -1) ? goodMeson_trk1_pt[index_pair[0]]: 0.f")
               .Define("meson_trk2_pt","(index_pair[0]!= -1) ? goodMeson_trk2_pt[index_pair[0]]: 0.f")
               .Define("photon_pt__div_sqrtHCandMass", "(index_pair[1]!= -1 && HCandMass>0) ? goodPhotons_pt[index_pair[1]]/sqrt(HCandMass): 0.f")
               .Define("meson_pt__div_sqrtHCandMass", "(index_pair[0]!= -1 && HCandMass>0) ? goodMeson_pt[index_pair[0]]/sqrt(HCandMass): 0.f")
               .Define("meson_vtx_prob","(index_pair[0]!= -1) ? goodMeson_vtx_prob[index_pair[0]]: 0.f")
               ## for ggH
               .Define("dPhiGammaMesonCand__div_sqrtHCandMass","(HCandMass>0) ? dPhiGammaMesonCand/sqrt(HCandMass): 0.f")
               ## for VBF
               .Define("dEtaGammaMesonCand__div_HCandMass","(HCandMass>0) ? dEtaGammaMesonCand/HCandMass: 0.f")
               # both GGH and VBF
               .Define("dEtaGammaMesonCand__div_sqrtHCandMass","(HCandMass>0) ? dEtaGammaMesonCand/sqrt(HCandMass): 0.f")
#               .Define("MVAdisc", ROOT.computeModel, list(variables))
                )

    dfWithMVA = tmva_helper.run_inference(dfWithMVA,"MVAdisc")

    return dfWithMVA

def DefineContent(branchList,isData):

    for branchName in [
            "HCandMass",
            "HCandMass_noregress",
#            "HCandMass_VTXcorr",
            "HCandPT",
            "index_pair",
            "meson_pt",
            "photon_pt",
            "sigmaHCandMass_Rel2",
            #
            "goodPhotons_pt",
            "goodPhotons_eta",
            "goodPhotons_pfRelIso03_all",
            "goodPhotons_hoe",
#            "goodPhotons_r9",
#            "goodPhotons_sieie",
            "goodPhotons_mvaID",
            "goodPhotons_energyErr",
            "nPhoton",
            "nGoodPhotons",
            "nPhotonsVeto",
            #
            "SoftActivityJetNjets5",
            "DeepMETResolutionTune_pt",
            "DeepMETResolutionTune_phi",
            "dPhiGammaMesonCand",
            "dEtaGammaMesonCand",
            "classify",
            #
            "triggerAna",
            # these below are for the 2018
#            "L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3",
#            "L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3",
#            "L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3",
#            "L1_SingleIsoEG32er2p1",
#            "L1_SingleIsoEG30er2p1",
#            "L1_SingleIsoEG28er2p1",
#            "L1_SingleLooseIsoEG30er1p5",
#            "L1_SingleLooseIsoEG28er1p5",
            #
            "w",
            "wraw",
            "w_allSF",
            "mc",
            "PV_npvsGood",
            "run",
            "luminosityBlock",
            "event",
            "lumiIntegrated",
    ]:
        branchList.push_back(branchName)

    if doPol:
        for branchName in [
                "theta",
        ]:
            branchList.push_back(branchName)

    if (isData == "false"):
        for branchName in [
            "indexMatch",
        ]:
            branchList.push_back(branchName)

    if (doSyst and isData == "false"):
        for branchName in [
                "L1PreFiringWeight_Nom",
                "L1PreFiringWeight_Up",
                "L1PreFiringWeight_Dn",
                "SFphoton_ID_Nom",
                "SFphoton_ID_Up",
                "SFphoton_ID_Dn",
                "SFpu_Nom",
                "SFpu_Up",
                "SFpu_Dn",
#                "photonTrigSF_Nom",
#                "photonTrigSF_Up",
#                "photonTrigSF_Dn",
#                "tauTrigSF_Nom",
#                "tauTrigSF_Up",
#                "tauTrigSF_Dn",
        ]:
            branchList.push_back(branchName)

    if (doSyst and isData == "false" and (isW or isZ)):
        for branchName in [
                "SFelectron_ID_Nom",
                "SFelectron_ID_Up",
                "SFelectron_ID_Dn",
                "SFmuon_ID_Nom",
                "SFmuon_ID_Up",
                "SFmuon_ID_Dn",
        ]:
            branchList.push_back(branchName)

    for branchName in [
            "goodMeson",
            "goodMeson_DR",
            "goodMeson_mass",
            "goodMeson_massErr",
            "goodMeson_pt",
            "goodMeson_iso",
            "goodMeson_trk1_pt",
            "goodMeson_trk2_pt",
            "goodMeson_trk1_eta",
            "goodMeson_trk2_eta",
            "goodMeson_vtx_chi2dof",
            "goodMeson_vtx_prob",
            "goodMeson_sipPV",
            "goodMeson_bestVtx_Z",
            "goodMeson_bestVtx_X",
            "goodMeson_bestVtx_Y",
            "goodMeson_bestVtx_idx",
    ]:
        branchList.push_back(branchName)

    if isZ or isW:
        for branchName in [
                "V_mass",
                "isMuorEle",
                "LeadingLeptonPt",
        ]:
            branchList.push_back(branchName)

    if isZ:
        for branchName in [
                "Z_veto1",
                "Z_veto2",
                "Visr_mass",
                "SubLeadingLeptonPt",
        ]:
            branchList.push_back(branchName)

    if isW:
        for branchName in [
                "deltaLepMeson",
                "dPhiGammaMET",
                "dPhiMesonMET",
                "Z_veto",
        ]:
            branchList.push_back(branchName)

    if isZinv:
        for branchName in [
                "dPhiGammaMET",
                "dPhiMesonMET",
                "ptRatioMEThiggs",
                "HCandPHI",
                "nbtag",
        ]:
            branchList.push_back(branchName)

    if (isGF or isVBF or isVBFlow or isZinv) and doMVA and (isRhoCat=="true" or isPhiCat=="true" or isK0StarCat=="true"):
        for branchName in [
                "MVAdisc",
        ]:
            branchList.push_back(branchName)

    if isGF and doMVA and (isOmega3PiCat=="true" or isPhi3PiCat=="true" or isD0Pi0StarCat=="true" or isD0StarCat=="true"):
        for branchName in [
                "meson_regress_pt_sf",
        ]:
            branchList.push_back(branchName)

    if isGF or isZinv:
        for branchName in [
                "nGoodJets",
        ]:
            branchList.push_back(branchName)

    if isVBF or isVBFlow:
        for branchName in [
                "mJJ",
                "nGoodJets",
                "dEtaJJ",
                "dPhiJJ",
                "Y1Y2",
                "deltaJetMeson",
                "deltaJetPhoton",
                "jet1Pt",
                "jet2Pt",
                "jet1Eta",
                "jet2Eta",
                "jet1hfsigmaPhiPhi",
                "jet1hfsigmaEtaEta",
                "jet2hfsigmaPhiPhi",
                "jet2hfsigmaEtaEta",
                "zepVar",
                "detaHigJet1",
                "detaHigJet2"
        ]:
            branchList.push_back(branchName)

    return branchList

def analysis(df,year,mc,sumw,isData,PDType):

    if doTrigger:
        dfCom = dfCommon(df,year,isData,mc,sumw,isVBF,isVBFlow,isGF,isZinv)
        dfOBJ = dfGammaMeson(dfCom,PDType)
        dfbase = dfHiggsCand(dfOBJ,isData)
#        dfbase = dfOBJ
        dfFINAL = selectionTAG(dfbase,doSyst,isData)
    else:
        dfCom = dfCommon(df,year,isData,mc,sumw,isVBF,isVBFlow,isGF,isZinv)
        dfOBJ = dfGammaMeson(dfCom,PDType)
        dfbase = dfHiggsCand(dfOBJ,isData)
        dfcandtag = selectionTAG(dfbase,doSyst,isData)
        dfpol_ = callPolarization(dfcandtag,mc)
        if (doSyst and isData == "false"):
            dfpreFINAL = dfwithSYST(dfpol_,year)
        else:
            dfpreFINAL = dfpol_.Define("w_allSF", "w")

        if doMVA and (isRhoCat=="true" or isPhiCat=="true" or isK0StarCat=="true"):
            dfFINAL = callMVA(dfpreFINAL,isVBF,isVBFlow,isGF,isZinv)
        else: dfFINAL = dfpreFINAL

    branchList = ROOT.vector('string')()

    if doTrigger:
        for branchName in [
                "HCandMass",
                "meson_pt",
                "photon_pt",
                #
                "mJJ",
                "nGoodJets",
                "dPhiGammaMesonCand",
                "dEtaGammaMesonCand",
                "dRGammaMesonCand",
                #
                "w",
                "wraw",
                "mc",
                "PV_npvsGood",
                "run",
                "luminosityBlock",
                "event",
                "lumiIntegrated",
        ]:
            branchList.push_back(branchName)

        for branchName in [
                "triggerAna",
                "TrigObj_filterBits",
                "TrigObj_eta",
                "TrigObj_pt",
                "TrigObj_id",
        ]:
            branchList.push_back(branchName)

        for branchName in [
                "index_pair",
                "goodPhotons_pt",
                "goodPhotons_eta",
                "goodPhotons_phi",
                "goodPhotons_pfRelIso03_all",
#                "goodPhotons_pixelSeed",
#                "goodPhotons_r9",
#                "goodPhotons_sieie",
                "goodPhotons_hoe",
                "goodPhotons_mvaID",

                "goodMeson_mass",
                "goodMeson_eta",
                "goodMeson_iso",
                "goodMeson_phi",
                "goodMeson_pt",
                "goodMeson_trk1_pt",
                "goodMeson_trk2_pt",
                "goodMeson_trk1_eta",
                "goodMeson_trk2_eta",
#                "goodMeson_vtx_chi2dof",
#                "goodMeson_vtx_prob",
#                "goodMeson_bestVtx_idx",
#                "goodMeson_bestVtx_Z",
#                "goodMeson_bestVtx_X",
#                "goodMeson_bestVtx_Y",

        ]:
            branchList.push_back(branchName)

    else:

        branchList = DefineContent(branchList,isData)


        if(isOmega3PiCat=="true" or isPhi3PiCat=="true" or isD0Pi0StarCat=="true"):
            for branchName in [
                    "goodMeson_mass_RAW",
                    "goodMeson_photon_pt",
                    "goodMeson_DR_trkPairPh",
                    "goodMeson_dEta_trkPairPh",
                    "goodMeson_charged_mass",
                    "goodMeson_charged_pt"
            ]:
                branchList.push_back(branchName)

        if(isD0StarCat=="true"):
            for branchName in [
                    "goodMeson_charged_mass",
                    "goodMeson_charged_pt"
            ]:
                branchList.push_back(branchName)

        if(isRhoCat=="true" or isK0StarCat=="true" or isPhiCat=="true"):
            for branchName in [
#                    "goodMeson_isoPho",
                    "goodMeson_isoNeuHad",
                    "meson_isoNeuHad",
            ]:
                branchList.push_back(branchName)

        if(isK0StarCat=="true"):
            for branchName in [
                    "K0Star_mask",
            ]:
                branchList.push_back(branchName)

    catM = ""
    if(isPhiCat=="true"): catM = "PhiCat"
    if(isRhoCat=="true"): catM = "RhoCat"
    if(isK0StarCat=="true"): catM = "K0StarCat"
    if(isOmega3PiCat=="true"): catM = "Omega3PiCat"
    if(isPhi3PiCat=="true"): catM = "Phi3PiCat"
    if(isD0StarCat=="true"): catM = "D0StarCat"
    if(isD0Pi0StarCat=="true"): catM = "D0Pi0StarCat"
    catTag = ""
    if isZ: catTag = "Zcat"
    if isZinv: catTag = "Zinvcat"
    if isW: catTag = "Wcat"
    if isVBF: catTag = "VBFcat"
    if isVBFlow: catTag = "VBFcatlow"
    if isGF: catTag = "GFcat"

    if True:
        outputFile = "/work/submit/mariadlf/NOV1/{0}/outname_mc{1}_{2}_{3}_{0}.root".format(year,mc,catTag,catM)
        if sys.argv[2]=='isOmega3PiCat' or sys.argv[2]=='isPhi3PiCat' or sys.argv[2]=='isD0StarCat' or sys.argv[2]=='isD0Pi0StarCat':
            outputFile = "/work/submit/mariadlf/NOVMarti/{0}/outname_mc{1}_{2}_{3}_{0}.root".format(year,mc,catTag,catM)
        print(outputFile)
        snapshotOptions = ROOT.RDF.RSnapshotOptions()
        snapshotOptions.fCompressionAlgorithm = ROOT.kLZ4
        snapshot_tdf = dfFINAL.Snapshot("events", outputFile, branchList, snapshotOptions)
        print("snapshot_tdf DONE")
        print(outputFile)

    if False:
        print("---------------- SUMMARY -------------")
        ## this doens't work with the negative weights
        report = dfFINAL.Report()
        report.Print()

    if doPlot and doSyst and isData == "false" and mc>1000:
        print("---------------- PLOTTING with SYST -------------")
        hists = {
            #        "Z_mass":     {"name":"Z_mass","title":"Di Muon mass; m_{#mu^{+}#mu^{-}} (GeV);N_{Events}","bin":500,"xmin":70,"xmax":120},
#            "V_mass":     {"name":"V_mass","title":"transverse mass; m_{T}(#mu^{+} MET} (GeV);N_{Events}","bin":80,"xmin":40,"xmax":120},
            "HCandMass":  {"name":"HCandMass","title":"H mass;m_{k^{+}k^{-}#gamma} (GeV);N_{Events}","bin":70,"xmin":100,"xmax":170},
#            "phi_num":    {"name":"nphi","title":"Phi N;N {k^{+}k^{-}} (GeV);N_{Events}","bin":10,"xmin":0.,"xmax":10.},
#            "Phi_mass":   {"name":"phi_kin_mass","title":"Phi mass;m_{k^{+}k^{-}} (GeV);N_{Events}","bin":200,"xmin":0.95,"xmax":1.15},
#            "Phi_pt":     {"name":"phi_kin_pt","title":"Phi pt ;p^{T}_{k^{+}k^{-}} (GeV);N_{Events}","bin":1000,"xmin":0.25,"xmax":50.25},
#            "Phi_gen_mass":   {"name":"phi_gen_mass","title":"Phi gen mass;m_{k^{+}k^{-}} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":10.},
#            "Phi_mass_err":   {"name":"phi_kin_massErr","title":"Phi mass error;m_{k^{+}k^{-}} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":0.5},
#            "Phi_kin_vtx_chi2dof":   {"name":"phi_kin_vtx_chi2dof","title":"Phi vtx_chi2dof;m_{k^{+}k^{-}} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":5.0},
        }

        histos = []
        for h in hists:

            # 1D is for nom only
            model = (hists[h]["name"], hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"])
#            h1d = dfFINAL.Histo1D(model, hists[h]["name"], "w")
            h1d = dfFINAL.Histo1D(model, hists[h]["name"], "w_allSF")
            histos.append(h1d)
#            h1d_noSF = dfFINAL.Histo1D(model, hists[h]["name"], "w")
#            histos.append(h1d_noSF)

            ## to use the SYST that change the variable
            hx = ROOT.RDF.Experimental.VariationsFor(h1d)
            hx["PhotonSYST:dn"].SetName(hists[h]["name"]+":PhotonSYST:dn")
            histos.append(hx["PhotonSYST:dn"])
            hx["PhotonSYST:up"].SetName(hists[h]["name"]+":PhotonSYST:up")
            histos.append(hx["PhotonSYST:up"])

            if isGF or isVBF or isVBFlow:
                ## to use the SYST that change the variable
                hx = ROOT.RDF.Experimental.VariationsFor(h1d)
                hx["JetSYST:dn"].SetName(hists[h]["name"]+":JetSYST:dn")
                histos.append(hx["JetSYST:dn"])
                hx["JetSYST:up"].SetName(hists[h]["name"]+":JetSYST:up")
                histos.append(hx["JetSYST:up"])

            ## those that change the weights only
            # 2D is for nom, up, down
            model2d_pu = (hists[h]["name"]+":PU", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
            histos.append(dfFINAL.Histo2D(model2d_pu, hists[h]["name"], "idx_nom_up_down", "pu_weights"))
            model2d_L1 = (hists[h]["name"]+":L1", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
            histos.append(dfFINAL.Histo2D(model2d_L1, hists[h]["name"], "idx_nom_up_down", "L1PreFiring_weights"))
            model2d_phoID = (hists[h]["name"]+":phoID", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
            histos.append(dfFINAL.Histo2D(model2d_phoID, hists[h]["name"], "idx_nom_up_down", "phoID_weights"))
            model2d_mesonRECO = (hists[h]["name"]+":mesonRECO", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
            histos.append(dfFINAL.Histo2D(model2d_mesonRECO, hists[h]["name"], "idx_nom_up_down", "mesonRECO_weights"))
            if isGF or isVBFlow:
                model2d_mesonChISO = (hists[h]["name"]+":mesonChISO", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
                histos.append(dfFINAL.Histo2D(model2d_mesonChISO, hists[h]["name"], "idx_nom_up_down", "chIsoSF_weights"))
                model2d_phoTrig = (hists[h]["name"]+":phoTrig", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
                histos.append(dfFINAL.Histo2D(model2d_phoTrig, hists[h]["name"], "idx_nom_up_down", "phTrigSF_weights"))
                model2d_tauTrig = (hists[h]["name"]+":tauTrig", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
                histos.append(dfFINAL.Histo2D(model2d_tauTrig, hists[h]["name"], "idx_nom_up_down", "tauTrigSF_weights"))
            if isW or isZ:
                model2d_eleID = (hists[h]["name"]+":eleID", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
                histos.append(dfFINAL.Histo2D(model2d_eleID, hists[h]["name"], "idx_nom_up_down", "eleID_weights"))
                model2d_eleRECO = (hists[h]["name"]+":eleRECO", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
                histos.append(dfFINAL.Histo2D(model2d_eleRECO, hists[h]["name"], "idx_nom_up_down", "eleReco_weights"))
                model2d_muoID = (hists[h]["name"]+":muoID", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
                histos.append(dfFINAL.Histo2D(model2d_muoID, hists[h]["name"], "idx_nom_up_down", "muoID_weights"))
                model2d_muoISO = (hists[h]["name"]+":muoISO", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
                histos.append(dfFINAL.Histo2D(model2d_muoISO, hists[h]["name"], "idx_nom_up_down", "muoISO_weights"))
            if isZ:
                model2d_muo2ID = (hists[h]["name"]+":muo2ID", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
                histos.append(dfFINAL.Histo2D(model2d_muo2ID, hists[h]["name"], "idx_nom_up_down", "muo2ID_weights"))
                model2d_muo2ISO = (hists[h]["name"]+":muo2ISO", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
                histos.append(dfFINAL.Histo2D(model2d_muo2ISO, hists[h]["name"], "idx_nom_up_down", "muo2ISO_weights"))
#                model2d_ele2ID = (hists[h]["name"]+":ele2ID", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
#                histos.append(dfFINAL.Histo2D(model2d_ele2ID, hists[h]["name"], "idx_nom_up_down", "ele2ID_weights"))
            model2d_PSwei = (hists[h]["name"]+":PSisr", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
            histos.append(dfFINAL.Histo2D(model2d_PSwei, hists[h]["name"], "idx_nom_up_down", "PS_isr_weights"))
            model2d_PSwei = (hists[h]["name"]+":PSfsr", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
            histos.append(dfFINAL.Histo2D(model2d_PSwei, hists[h]["name"], "idx_nom_up_down", "PS_fsr_weights"))

#        evtcounts = []
#        evtcount = dfFINAL.Count()
#        evtcounts.append(evtcount)
#        ROOT.ROOT.RDF.RunGraphs(evtcounts)

        outputFileHisto = "/work/submit/mariadlf/OCT_SYSTvar/{0}/histoname_mc{1}_{2}_{3}_{0}_wSF.root".format(year,mc,catTag,catM,year)
        print(outputFileHisto)
        myfile = ROOT.TFile(outputFileHisto,"RECREATE")

        for h in histos:
            h.Write()
        myfile.Close()
        myfile.Write()

def readMCSample(year,sampleNOW,useD03):

    files = getMClist(year,sampleNOW,useD03)
    print(len(files))
    #local
    df = ROOT.RDataFrame("Events", files)
    rdf = ROOT.RDataFrame("Runs", files)

    sumW = computeWeigths(df, rdf, sampleNOW, year, True, useD03)
    if (doPol and sampleNOW>1000 and sampleNOW<1039): loadPolarizationTree(sampleNOW,year)
    if doSyst: loadSFhisto(sampleNOW,year)
    loadCorrectionSet(year)
    analysis(df,year,sampleNOW,sumW,"false","NULL")

def readDataSample(year,datasetNumber):

    pair = getDATAlist(datasetNumber,year)
    files = pair[0]
    PDType = pair[1]
    print(len(files))
    print(PDType)

    #local
    df = ROOT.RDataFrame("Events", files)
    nevents = df.Count().GetValue()
    print("%s entries in the dataset" %nevents)

    analysis(df,year,datasetNumber, 1. ,"true",PDType)

def readDataSkims(datasetNumber,year,category,useD03):

    print("enum",datasetNumber)
    print("year",year)
    print("cat",category)
    if (category=="isZtag" or category=="isWtag"):
        pair = getSkims(datasetNumber,year,"VH",useD03)
    elif category=="isVBFtag":
        pair = getSkims(datasetNumber,year,"VBF",useD03)
    if (category=="isZinvtag" or category=="isVBFtaglow" or category=="isGFtag"):
        pair = getSkims(datasetNumber,year,"Zinv",useD03)

    files = pair[0]
    PDType = pair[1]
    print(len(files))
    print(PDType)

    #local
    df = ROOT.RDataFrame("Events", files)
    nevents = df.Count().GetValue()
    print("%s entries in the dataset" %nevents)

    analysis(df,year,datasetNumber,1.,"true",PDType)
    print("***ANALYSIS DONE ***")

def runTest():

    df = ROOT.RDataFrame("Events", "root://eoscms.cern.ch//eos/cms//store/group/phys_higgs/HiggsExo/dalfonso/Hrare/D01/vbf-hphigamma-powheg/NANOAOD_01/step7_VBS_Phigamma_8.root")

    w=1.
    nevents = df.Count().GetValue()
    print("%s entries in the dataset" %nevents)

    sampleNOW=-1
    analysis(df,-1,w,"false")

   
if __name__ == "__main__":

    loadUserCode()
    loadtmva_helper()

#    runTest()
#    to run: python3 -i VGammaMeson_cat.py isVBFtag isPhiCat 12 2018
    print(int(sys.argv[3]))

    if ( sys.argv[1]=="isVBFtag" and int(sys.argv[3]) in [ -31, -32, -33, -34, -76, -81, -82, -83, -84, -85, -86, -62, -63, -64, -65, -66]):
        readDataSkims(int(sys.argv[3]),int(sys.argv[4]),sys.argv[1],False) # skims VBF

    elif ( (sys.argv[1]=="isZinvtag" or sys.argv[1]=="isVBFtaglow" or sys.argv[1]=="isGFtag") and int(sys.argv[3]) in [-62, -63, -64, -65, -66]):
        readDataSkims(int(sys.argv[3]),int(sys.argv[4]),sys.argv[1],useD03) # skims Tau

    elif ( (sys.argv[1]=="isWtag" or sys.argv[1]=="isZtag") and int(sys.argv[3]) in [-1, -2, -3, -4, -5, -6, -7, -8]):
        readDataSkims(int(sys.argv[3]),int(sys.argv[4]),sys.argv[1],False) # skims singleMu
    elif ( (sys.argv[1]=="isWtag" or sys.argv[1]=="isZtag") and int(sys.argv[3]) in [-11, -12, -13, -14, -15, -16, -17, -18]):
        readDataSkims(int(sys.argv[3]),int(sys.argv[4]),sys.argv[1],False) # skims doubleMu
    elif ( (sys.argv[1]=="isWtag" or sys.argv[1]=="isZtag") and int(sys.argv[3]) in [-21, -22, -23, -24, -25, -26, -27, -28]):
        readDataSkims(int(sys.argv[3]),int(sys.argv[4]),sys.argv[1],False) # skims MuEG
    elif ( (sys.argv[1]=="isWtag" or sys.argv[1]=="isZtag") and int(sys.argv[3]) in [-31, -32, -33, -34, -35, -36, -37, -38]):
        readDataSkims(int(sys.argv[3]),int(sys.argv[4]),sys.argv[1],False) # skims EG
    elif ( (sys.argv[1]=="isWtag" or sys.argv[1]=="isZtag") and int(sys.argv[3]) in [-41, -42, -43, -44, -45, -46, -47, -48]):
        readDataSkims(int(sys.argv[3]),int(sys.argv[4]),sys.argv[1],False) # skims DoubleEG
    elif ( (sys.argv[1]=="isWtag" or sys.argv[1]=="isZtag") and int(sys.argv[3]) in [-51, -52, -53, -54, -55, -56, -57, -58]):
        readDataSkims(int(sys.argv[3]),int(sys.argv[4]),sys.argv[1],False) # skims SingleElectron
    elif(int(sys.argv[3]) < 0):
        readDataSample(int(sys.argv[4]),int(sys.argv[3]) )  # DATA
    else: readMCSample(int(sys.argv[4]),int(sys.argv[3]),useD03) # to switch sample
