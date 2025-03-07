import ROOT
import json
import sys

from utilsAna import loadUserCode
from utilsAna import BuildDictJpsiCC, SwitchSample
from utilsAna import readDataQuality
from utilsAna import computeWeigths
from utilsAna import lumisJpsiCC

from datetime import datetime

ROOT.gROOT.SetBatch()
ROOT.ROOT.EnableImplicitMT()

loadUserCode()
import helper_tmva # need various definitions

with open("./config/qualityData.json") as qualityJsonFile:
    qualObject = json.load(qualityJsonFile)
    qualityJsonFile.close()

JSON = qualObject['JSON']

'''
chain = ROOT.TChain("Events")
chain.Add("/data/submit/cms/store/user/mariadlf/nano/GluGluH_HJPsiCC/NANOAOD_test3/*.root")
chainR = ROOT.TChain("Runs")
chainR.Add("/data/submit/cms/store/user/mariadlf/nano/GluGluH_HJPsiCC/NANOAOD_test3/*.root")
sumw = computeWeigths(chainR)
'''

#print('sumw',sumw)
#weight = "{0}*genWeight*{1}".format(lumi,10000000)
#weight = "{0}*genWeight".format(lumi)

#lumi = 59. #2018
#year = 2018
#TRIGGER = "HLT_Dimuon20_Jpsi_Barrel_Seagulls or HLT_Dimuon25_Jpsi"
#TRIGGER2017 = "HLT_Dimuon25_Jpsi"
#TRIGGER2018 = "(HLT_Dimuon25_Jpsi or HLT_Dimuon25_Jpsi_noCorrL1)"
#TRIGGER3 = "(HLT_Dimuon25_Jpsi or HLT_Dimuon25_Jpsi_noCorrL1 or HLT_Dimuon0_Jpsi_L1_NoOS or HLT_Dimuon0_Jpsi_NoVertexing_NoOS or HLT_Dimuon0_Jpsi_NoVertexing)"

#jetDef='Jet'
jetDef='JetNoMuFromJpsi'

tagger='Deep'
tagger='DeepFlav'
#if year=='12022': tagger='PNet'

#regress='cRegCorr'
#regress='PNetRegPtRawCorrNeutrino'
regress='{}_cRegCorr[goodJets][index_CloseFar[1]]'.format(jetDef)
regressClose='{}_cRegCorr[goodJets][index_CloseFar[0]]'.format(jetDef)

regress='1.0'
regressClose='1.0'

#GOODJETS = "(Jet_pt>30 and abs(Jet_eta)<2.5 and Jet_btag{}CvL>-1)".format(tagger)
GOODJETS = "({0}_pt>20 and abs({0}_eta)<2.5)".format(jetDef)

#GOODJPSI = "Jpsi_kin_pt>20 and Jpsi_kin_sipPV<1.5" #Jpsi_kin_sipPV<2 reduce the BToJpsi (IP 3D significance) # at skim level already vtx_prob > 0.001

def callMVARegression(df,target):

    fileName = "JPSImva/HJpsiCC_regression_model_"+target+".root"
    print(fileName)
    tmva_helper = helper_tmva.TMVAHelperXGB(fileName, "regress_model")
    print(tmva_helper.variables)

    newCol="regr_"+target
    print(newCol)
    dfWithMVA = tmva_helper.run_inference(df,newCol)

    return dfWithMVA


def callMVAclassification(df):

    tmva_helper = helper_tmva.TMVAHelperXGB("JPSImva/HJpsiCC_classification_model.root", "bdt_model")
    print(tmva_helper.variables)

    dfWithMVA = tmva_helper.run_inference(df,"discrMVA")

    return dfWithMVA

def addPfCands(df,mc):

    df = (df.Define("pfcand_close_pt","buildConstVar(Jpsi_muon1_pt[0],Jpsi_muon2_pt[0],PFCand_pt,JetPFCands_candIdx,JetPFCands_jetIdx,index_CloseFar[0],-1)")
          .Define("pfcand_close_eta","buildConstVar(Jpsi_muon1_eta[0],Jpsi_muon2_eta[0],PFCand_eta,JetPFCands_candIdx,JetPFCands_jetIdx,index_CloseFar[0],-1)")
          .Define("pfcand_close_phi","buildConstVar(Jpsi_muon1_phi[0],Jpsi_muon2_phi[0],PFCand_phi,JetPFCands_candIdx,JetPFCands_jetIdx,index_CloseFar[0],-1)")
          .Define("pfcand_close_dxy","buildConstVar(Jpsi_muon1_dxy[0],Jpsi_muon2_dxy[0],PFCand_dxy,JetPFCands_candIdx,JetPFCands_jetIdx,index_CloseFar[0],-1)")
          .Define("pfcand_close_dz","buildConstVar(Jpsi_muon1_dz[0],Jpsi_muon2_dz[0],PFCand_dz,JetPFCands_candIdx,JetPFCands_jetIdx,index_CloseFar[0],-1)")
          ## add the dzsig = dz dzerror (same for dx)
          ## need to fix those  pdg for muons for the muons
          .Define("pfcand_close_isMu","buildConstPDG(Jpsi_muon1_charge[0],Jpsi_muon2_charge[0],PFCand_pdgId,JetPFCands_candIdx,JetPFCands_jetIdx,index_CloseFar[0],13)")
          .Define("pfcand_close_isEl","buildConstPDG(Jpsi_muon1_charge[0],Jpsi_muon2_charge[0],PFCand_pdgId,JetPFCands_candIdx,JetPFCands_jetIdx,index_CloseFar[0],11)")
          .Define("pfcand_close_isChargedHad","buildConstPDG(Jpsi_muon1_charge[0],Jpsi_muon2_charge[0],PFCand_pdgId,JetPFCands_candIdx,JetPFCands_jetIdx,index_CloseFar[0],211)")
          .Define("pfcand_close_isGamma","buildConstPDG(Jpsi_muon1_charge[0],Jpsi_muon2_charge[0],PFCand_pdgId,JetPFCands_candIdx,JetPFCands_jetIdx,index_CloseFar[0],22)")
          .Define("pfcand_close_isNeutralHad","buildConstPDG(Jpsi_muon1_charge[0],Jpsi_muon2_charge[0],PFCand_pdgId,JetPFCands_candIdx,JetPFCands_jetIdx,index_CloseFar[0],130)")
          #
          .Define("jet_npfcand_close","buildConstN(event,1,1,PFCand_pt[JetPFCands_candIdx],JetPFCands_jetIdx,index_CloseFar[0])")
          # need to fix those btag for the muons
          .Define("pfcand_close_btagEtaRel","buildConstVarJet(0.f,0.f,JetPFCands_btagEtaRel,JetPFCands_jetIdx,index_CloseFar[0])")
          .Define("pfcand_close_btagPtRatio","buildConstVarJet(0.f,0.f,JetPFCands_btagPtRatio,JetPFCands_jetIdx,index_CloseFar[0])")
          .Define("pfcand_close_btagPParRatio","buildConstVarJet(0.f,0.f,JetPFCands_btagPParRatio,JetPFCands_jetIdx,index_CloseFar[0])")
          #
          #https://github.com/cms-sw/cmssw/blob/862b9fd82f99aae6090c6abe750de9f312432d45/DataFormats/BTauReco/interface/TaggingVariable.h#L22
          #trackEtaRel_ = reco::btau::etaRel(jetDir, candidate->momentum());
          #https://github.com/cms-sw/cmssw/blob/862b9fd82f99aae6090c6abe750de9f312432d45/RecoBTag/FeatureTools/src/TrackInfoBuilder.cc#L19
          #TVector3 jetDir3(jetDir.x(), jetDir.y(), jetDir.z());
          #TVector3 trackMom3(candidate->momentum().x(), candidate->momentum().y(), candidate->momentum().z());
          #trackPtRel_ = trackMom3.Perp(jetDir3);
          #trackPParRatio_ = jetDir.Dot(candidate->momentum()) / candidate->p();
          #
          .Define("pfcand_close_btagSip3dSig","buildConstVarJet(Jpsi_muon1_Sip3dSig[0],Jpsi_muon2_Sip3dSig[0],JetPFCands_btagSip3dSig,JetPFCands_jetIdx,index_CloseFar[0])")
          .Define("pfcand_close_btagSip3dVal","buildConstVarJet(Jpsi_muon1_Sip3dVal[0],Jpsi_muon2_Sip3dVal[0],JetPFCands_btagSip3dVal,JetPFCands_jetIdx,index_CloseFar[0])")
          .Define("pfcand_close_btagSip2dSig","buildConstVarJet(Jpsi_muon1_Sip2dSig[0],Jpsi_muon2_Sip2dSig[0],JetPFCands_btagSip2dSig,JetPFCands_jetIdx,index_CloseFar[0])")
          .Define("pfcand_close_btagSip2dVal","buildConstVarJet(Jpsi_muon1_Sip2dVal[0],Jpsi_muon2_Sip2dVal[0],JetPFCands_btagSip2dVal,JetPFCands_jetIdx,index_CloseFar[0])")
          .Define("pfcand_close_btagJetDistVal","buildConstVarJet(Jpsi_muon1_DistVal[0],Jpsi_muon2_DistVal[0],JetPFCands_btagJetDistVal,JetPFCands_jetIdx,index_CloseFar[0])")
          ##
          .Define("pfcand_close_ptrel_log","buildConstVar(Jpsi_muon1_pt[0],Jpsi_muon2_pt[0],PFCand_pt,JetPFCands_candIdx,JetPFCands_jetIdx,index_CloseFar[0],jetClose_Pt,0)")
          .Define("pfcand_close_phirel","buildConstVar(Jpsi_muon1_phi[0],Jpsi_muon2_phi[0],PFCand_phi,JetPFCands_candIdx,JetPFCands_jetIdx,index_CloseFar[0],jetClose_Phi,1)")
          .Define("pfcand_close_etarel","buildConstVar(Jpsi_muon1_eta[0],Jpsi_muon2_eta[0],PFCand_eta,JetPFCands_candIdx,JetPFCands_jetIdx,index_CloseFar[0],jetClose_Eta,2)")
          .Define("pfcand_close_charge","buildConstVar(Jpsi_muon1_charge[0],Jpsi_muon2_charge[0],PFCand_pdgId,JetPFCands_candIdx,JetPFCands_jetIdx,index_CloseFar[0],3)")
          .Define("pfcand_close_mass","buildConstVar(Jpsi_muon1_charge[0],Jpsi_muon2_charge[0],PFCand_pdgId,JetPFCands_candIdx,JetPFCands_jetIdx,index_CloseFar[0],4)")
          #####
          .Define("recojet_close_isSig", "mc==1000 ? 1 : 0 ")
#          .Define("recojet_close_isSig", "mc==1000 ? 0 : 1 ")
          .Define("recojet_close_isB", "mc==11 ? 1 : 0 ") #BToJpsi_JPsiToMuMu_BMuonFilter
#          .Define("recojet_close_isB", "1") #BToJpsi_JPsiToMuMu_BMuonFilter
          )
    return df

def addJetProperties(df,mc):

    df = (df.Define("jet1Pt","{0}_pt[goodJets][0]".format(jetDef))
          .Define("jet2Pt","{0}_pt[goodJets][1]".format(jetDef))
          .Define("jet1Eta","{0}_eta[goodJets][0]".format(jetDef))
          .Define("jet2Eta","{0}_eta[goodJets][1]".format(jetDef))
          .Define("jet1CvL","{0}_btag{1}CvL[goodJets][0]".format(jetDef,tagger))
          .Define("jet2CvL","{0}_btag{1}CvL[goodJets][1]".format(jetDef,tagger))
          .Define("jet1CvB","{0}_btag{1}CvB[goodJets][0]".format(jetDef,tagger))
          .Define("jet2CvB","{0}_btag{1}CvB[goodJets][1]".format(jetDef,tagger))
          #
          # for the FAR we apply JEC and c-Regression
          # for the CLOSE no JEC and no c-Regression
          .Define("jetFarPtRaw","{0}_pt[goodJets][index_CloseFar[1]]".format(jetDef))
          #
          .Define("jetClosecRegCorr","{}".format(regressClose))
          .Define("jetFarcRegCorr","{}".format(regress))
          #
          .Define("jetCloseJPsiRatio","Jpsi_kin_pt[0]/jetClose_Pt")
          .Define("jetCloseTrk1Ratio","Jpsi_leadingChargedMaxPt1[0]/jetClose_Pt")
          .Define("jetCloseTrk2Ratio","Jpsi_leadingChargedMaxPt2[0]/jetClose_Pt")
          .Define("jetClose_CvL","{0}_btag{1}CvL[goodJets][index_CloseFar[0]]".format(jetDef,tagger))
          .Define("jetFar_CvL","{0}_btag{1}CvL[goodJets][index_CloseFar[1]]".format(jetDef,tagger))
          .Define("jetClose_CvB","{0}_btag{1}CvB[goodJets][index_CloseFar[0]]".format(jetDef,tagger))
          .Define("jetFar_CvB","{0}_btag{1}CvB[goodJets][index_CloseFar[1]]".format(jetDef,tagger))
          .Define("jetClose_nConst","{0}_nConstituents[goodJets][index_CloseFar[0]]-2.0".format(jetDef))
          .Define("jetFar_nConst","{0}_nConstituents[goodJets][index_CloseFar[1]]".format(jetDef))
          .Define("jetClose_nMuons","{0}_nMuons[goodJets][index_CloseFar[0]]".format(jetDef))
          .Define("jetFar_nMuons","{0}_nMuons[goodJets][index_CloseFar[1]]".format(jetDef))
          .Define("jetClose_nElectrons","{0}_nElectrons[goodJets][index_CloseFar[0]]".format(jetDef))
          .Define("jetFar_nElectrons","{0}_nElectrons[goodJets][index_CloseFar[1]]".format(jetDef))
          .Define("jetFar_chEmEF","{0}_chEmEF[goodJets][index_CloseFar[1]]".format(jetDef))
          .Define("jetFar_chHEF","{0}_chHEF[goodJets][index_CloseFar[1]]".format(jetDef))
          .Define("jetFar_neEmEF","{0}_neEmEF[goodJets][index_CloseFar[1]]".format(jetDef))
          .Define("jetFar_neHEF","{0}_neHEF[goodJets][index_CloseFar[1]]".format(jetDef))
          .Define("jetFar_muEF","{0}_muEF[goodJets][index_CloseFar[1]]".format(jetDef))
          .Define("jetClose_chEmEF","{0}_chEmEF[goodJets][index_CloseFar[0]]".format(jetDef))
          .Define("jetClose_chHEF","{0}_chHEF[goodJets][index_CloseFar[0]]".format(jetDef))
          .Define("jetClose_neEmEF","{0}_neEmEF[goodJets][index_CloseFar[0]]".format(jetDef))
          .Define("jetClose_neHEF","{0}_neHEF[goodJets][index_CloseFar[0]]".format(jetDef))
          .Define("jetClose_muEF","{0}_muEF[goodJets][index_CloseFar[0]]".format(jetDef))
          )

    #chFPV0EF, puIdDisc, puId, jetId, qgl
    if mc>0 and True:
        df= (df.Define("jet1partonFlavour","abs({0}_partonFlavour[goodJets][0])".format(jetDef))
             .Define("jet2partonFlavour","abs({0}_partonFlavour[goodJets][1])".format(jetDef))
             .Define("jetFarPtCharm","(abs({0}_partonFlavour[goodJets][index_CloseFar[1]])==4) ? {1} * {0}_pt[goodJets][index_CloseFar[1]]:-1.".format(jetDef,regress))
             .Define("jetFarPtGluon","(abs({0}_partonFlavour[goodJets][index_CloseFar[1]])==21) ? {1} * {0}_pt[goodJets][index_CloseFar[1]]:-1.".format(jetDef,regress))
             .Define("jetFarCvLCharm","(abs({0}_partonFlavour[goodJets][index_CloseFar[1]])==4) ? {0}_btag{1}CvL[goodJets][index_CloseFar[1]]:-1.".format(jetDef,tagger))
             .Define("jetFarCvLGluon","(abs({0}_partonFlavour[goodJets][index_CloseFar[1]])==21) ? {0}_btag{1}CvL[goodJets][index_CloseFar[1]]:-1.".format(jetDef,tagger))
             .Define("jetClose_partonFlavour","abs({0}_partonFlavour[goodJets][index_CloseFar[0]])".format(jetDef))
             .Define("jetFar_partonFlavour","abs({0}_partonFlavour[goodJets][index_CloseFar[1]])".format(jetDef))
             .Define("jetCloseScale","jetClose_Pt/GenJet_pt[{0}_genJetIdx[goodJets][index_CloseFar[0]]]".format(jetDef))
             .Define("jetFarScale","jetFar_Pt/GenJet_pt[{0}_genJetIdx[goodJets][index_CloseFar[1]]]".format(jetDef))
             .Define("jetClosecRegCorrScale","{0}_pt[goodJets][index_CloseFar[0]]*{1}/{0}_pt[{0}_genJetIdx[goodJets][index_CloseFar[0]]]".format(jetDef,regressClose))
             .Define("jetFarcRegCorrScale","{0}_pt[goodJets][index_CloseFar[1]]*{1}/{0}_pt[{0}_genJetIdx[goodJets][index_CloseFar[1]]]".format(jetDef,regress))
             #
             )

    return df

def analysis(files,year,mc,sumW):
    dfINI = ROOT.RDataFrame("Events", files)

    isData = "false"
    if mc<0: isData = "true"

    lumiIntegrated=0.
    if (isData == "false"):
        lumiIntegrated = lumisJpsiCC[year]
        print('lumiIntegrated = ',lumiIntegrated)

    weight = "{0}".format(1.)
    if mc>=0: weight = "{0}*genWeight*{1}".format(lumiIntegrated,sumW)


    TRIGGER=''
    if year=='12022' or year=='22022' or year=='12023' or year=='22023':
        TRIGGER="HLT_DoubleMu4_3_LowMass or (HLT_Dimuon25_Jpsi or HLT_Dimuon25_Jpsi_noCorrL1) or (HLT_Dimuon0_Jpsi_L1_NoOS or HLT_Dimuon0_Jpsi_NoVertexing_NoOS or HLT_Dimuon0_Jpsi or HLT_Dimuon0_Jpsi_NoVertexing)"
        GOODJPSI = "Jpsi_kin_pt>20"
    if year=='2018':
        TRIGGER="(HLT_Dimuon25_Jpsi or HLT_Dimuon25_Jpsi_noCorrL1)"
        GOODJPSI = "Jpsi_kin_pt>25"
    if year=='2017':
        TRIGGER="HLT_Dimuon25_Jpsi"
        GOODJPSI = "Jpsi_kin_pt>25"
    if year=='12016' or year=='22016':
        TRIGGER="(HLT_Dimuon20_Jpsi or HLT_Dimuon16_Jpsi)"
        GOODJPSI = "Jpsi_kin_pt>20"
    print(TRIGGER)

    dfComm = (dfINI
              .Define("mc","{}".format(mc))
              .Define("isData","{}".format(isData))
#              .Define("applyJson","{}".format(JSON)).Filter("applyJson","pass JSON") # doesn't work for signal ??
              .Define("w","{}".format(weight))
              .Define("wraw","{}".format(weight))
              .Define("lumiIntegrated","{}".format(lumiIntegrated))
              .Filter("PV_npvsGood>0","one good PV")
              .Define("triggerAna","{}".format(TRIGGER))
              .Filter("triggerAna>0","passing trigger")
              )

    df= (dfComm
         .Filter("nJpsi>0","at least one Jpsi;")
         .Define("goodMeson","{}".format(GOODJPSI))
         .Filter("Sum(goodMeson)>0", "one good Jpsi")
         .Define("goodJets","{}".format(GOODJETS))
         .Define("nGoodJets","Sum(goodJets)*1.0f").Filter("Sum(goodJets)>1","two jets for cc")
         .Define("mJJ","Minv({0}_pt[goodJets], {0}_eta[goodJets], {0}_phi[goodJets], {0}_mass[goodJets])".format(jetDef))
         .Define("massHiggs","Minv3massive({0}_pt[goodJets], {0}_eta[goodJets], {0}_phi[goodJets], {0}_mass[goodJets], Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass)".format(jetDef))
         #
         .Define("index_CloseFar","jetCloseFar({0}_pt[goodJets],{0}_eta[goodJets],{0}_phi[goodJets],{0}_mass[goodJets],Jpsi_kin_pt,Jpsi_kin_eta, Jpsi_kin_phi,Jpsi_kin_mass)".format(jetDef))
         .Filter("index_CloseFar[0]!= -1", "at least one close jet")
         # for the FAR we apply JEC and c-Regression
         # for the CLOSE no JEC and no c-Regression
         .Define("jetFar_Pt","{0}_pt[goodJets][index_CloseFar[1]]*{1}".format(jetDef,regress))
         .Filter("jetFar_Pt>20","20 GeV threshould for jetFar")
         .Define("jetFar_Eta","{0}_eta[goodJets][index_CloseFar[1]]".format(jetDef))
         .Define("jetFar_Phi","{0}_phi[goodJets][index_CloseFar[1]]".format(jetDef))
         .Define("jetFar_Mass","{0}_mass[goodJets][index_CloseFar[1]]".format(jetDef))
          #
         .Define("jetClose_Pt","{0}_pt[goodJets][index_CloseFar[0]]".format(jetDef))
         .Define("jetClose_Eta","{0}_eta[goodJets][index_CloseFar[0]]".format(jetDef))
         .Define("jetClose_Phi","{0}_phi[goodJets][index_CloseFar[0]]".format(jetDef))
         .Define("jetClose_Mass","{0}_mass[goodJets][index_CloseFar[0]]".format(jetDef))
#         .Define("jetClosePt","buildCloseJets(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets], Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass, index_CloseFar[0],0)")
         # this is needed if we do the subtraction a posteriori
#         .Define("jetClose_Pt","buildCloseJets({0}_pt[goodJets]*(1-{0}_rawFactor[goodJets]), {0}_eta[goodJets], {0}_phi[goodJets], {0}_mass[goodJets], Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass, index_CloseFar[0],0)".format(jetDef))
#         .Define("jetClose_Eta","buildCloseJets({0}_pt[goodJets], {0}_eta[goodJets], {0}_phi[goodJets], {0}_mass[goodJets], Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass, index_CloseFar[0],1)".format(jetDef))
#         .Define("jetClose_Phi","buildCloseJets({0}_pt[goodJets], {0}_eta[goodJets], {0}_phi[goodJets], {0}_mass[goodJets], Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass, index_CloseFar[0],2)".format(jetDef))
#         .Define("jetClose_Mass","buildCloseJets({0}_pt[goodJets], {0}_eta[goodJets], {0}_phi[goodJets], {0}_mass[goodJets], Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass, index_CloseFar[0],3)".format(jetDef))
         #
         .Define("massHiggsCorr","Minv3body(jetClose_Pt,jetClose_Eta,jetClose_Phi,jetClose_Mass,jetFar_Pt,jetFar_Eta,jetFar_Phi,jetFar_Mass,Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass,0 )")
         .Define("HiggsPt","Minv3body(jetClose_Pt,jetClose_Eta,jetClose_Phi,jetClose_Mass,jetFar_Pt,jetFar_Eta,jetFar_Phi,jetFar_Mass,Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass, 1)")
         #
         .Define("jetFar_PTMHRatio","jetFar_Pt/massHiggsCorr")
         .Define("jpsi_Higgs_Ratio","Jpsi_kin_pt[0]/HiggsPt")
         #
         .Define("minDRjpsi","std::min(deltaR({0}_eta[goodJets][0], {0}_phi[goodJets][0], Jpsi_kin_eta[0], Jpsi_kin_phi[0]),deltaR({0}_eta[goodJets][1], {0}_phi[goodJets][1], Jpsi_kin_eta[0], Jpsi_kin_phi[0]))".format(jetDef))
#         .Filter("minDRjpsi<0.3","jPsi very close to 1 charm-jet")
         #     .Filter("abs(Jpsi_kin_eta[0])>1.4","barrel JPsi")
         #     .Filter("abs(Jpsi_leadingChargedDxy[0])>0.5","bad track")
         .Define("angle_phi1","get_phi1(Jpsi_muon1_pt[0], Jpsi_muon1_eta[0], Jpsi_muon1_phi[0], Jpsi_muon2_pt[0], Jpsi_muon2_eta[0], Jpsi_muon2_phi[0], jetClose_Pt,jetClose_Eta,jetClose_Phi,jetClose_Mass,jetFar_Pt,jetFar_Eta,jetFar_Phi,jetFar_Mass,Jpsi_kin_pt[0],Jpsi_kin_eta[0],Jpsi_kin_phi[0],Jpsi_kin_mass[0])")
         .Define("angle_phi2","get_phi2(Jpsi_muon1_pt[0], Jpsi_muon1_eta[0], Jpsi_muon1_phi[0], Jpsi_muon2_pt[0], Jpsi_muon2_eta[0], Jpsi_muon2_phi[0], jetClose_Pt,jetClose_Eta,jetClose_Phi,jetClose_Mass,jetFar_Pt,jetFar_Eta,jetFar_Phi,jetFar_Mass,Jpsi_kin_pt[0],Jpsi_kin_eta[0],Jpsi_kin_phi[0],Jpsi_kin_mass[0])")
         .Define("angle_phi3","get_phi3(jetClose_Pt,jetClose_Eta,jetClose_Phi,jetClose_Mass,jetFar_Pt,jetFar_Eta,jetFar_Phi,jetFar_Mass,Jpsi_kin_pt[0],Jpsi_kin_eta[0],Jpsi_kin_phi[0],Jpsi_kin_mass[0])")
         .Define("angle_phi4","get_phi4(Jpsi_muon1_pt[0], Jpsi_muon1_eta[0], Jpsi_muon1_phi[0], Jpsi_muon2_pt[0], Jpsi_muon2_eta[0], Jpsi_muon2_phi[0], jetClose_Pt,jetClose_Eta,jetClose_Phi,jetClose_Mass,jetFar_Pt,jetFar_Eta,jetFar_Phi,jetFar_Mass,Jpsi_kin_pt[0],Jpsi_kin_eta[0],Jpsi_kin_phi[0],Jpsi_kin_mass[0])")
         .Define("angle_cos_delta_m_min","get_cos_delta_m_min(Jpsi_muon1_pt[0], Jpsi_muon1_eta[0], Jpsi_muon1_phi[0], Jpsi_muon2_pt[0], Jpsi_muon2_eta[0], Jpsi_muon2_phi[0], jetClose_Pt,jetClose_Eta,jetClose_Phi,jetClose_Mass,jetFar_Pt,jetFar_Eta,jetFar_Phi,jetFar_Mass,Jpsi_kin_pt[0],Jpsi_kin_eta[0],Jpsi_kin_phi[0],Jpsi_kin_mass[0])")
         .Define("angle_cos_theta_jpsi","get_cos_theta_jpsi(Jpsi_muon1_pt[0], Jpsi_muon1_eta[0], Jpsi_muon1_phi[0], Jpsi_muon2_pt[0], Jpsi_muon2_eta[0], Jpsi_muon2_phi[0], jetClose_Pt,jetClose_Eta,jetClose_Phi,jetClose_Mass,jetFar_Pt,jetFar_Eta,jetFar_Phi,jetFar_Mass,Jpsi_kin_pt[0],Jpsi_kin_eta[0],Jpsi_kin_phi[0],Jpsi_kin_mass[0])")
         .Define("angle_cos_theta_dicharm","get_cos_theta_jpsi(Jpsi_muon1_pt[0], Jpsi_muon1_eta[0], Jpsi_muon1_phi[0], Jpsi_muon2_pt[0], Jpsi_muon2_eta[0], Jpsi_muon2_phi[0], jetClose_Pt,jetClose_Eta,jetClose_Phi,jetClose_Mass,jetFar_Pt,jetFar_Eta,jetFar_Phi,jetFar_Mass,Jpsi_kin_pt[0],Jpsi_kin_eta[0],Jpsi_kin_phi[0],Jpsi_kin_mass[0])")
         )

    df = addJetProperties(df,mc)
    df = addPfCands(df,mc)

    if mc==1000 or mc==1001:
        df= (df.Define("idxs", "genMatch_indices(GenPart_pdgId, GenPart_genPartIdxMother)")
             .Define("muon1_idx", 'idxs.at("muon1_idx")')
             .Define("muon2_idx", 'idxs.at("muon2_idx")')
             .Define("neg_muon_idx", 'idxs.at("neg_muon_idx")')
             .Define("pos_muon_idx", 'idxs.at("pos_muon_idx")')
             .Define("c1_idx", 'idxs.at("c1_idx")')
             .Define("c2_idx", 'idxs.at("c2_idx")')
             .Define("higgs_idx", 'idxs.at("higgs_idx")')
             .Define("charm1_gen_Pt",'GenPart_pt[c1_idx]')
             .Define("charm2_gen_Pt",'GenPart_pt[c2_idx]')
             .Define("muon1_gen_Pt",'GenPart_pt[muon1_idx]')
             .Define("muon2_gen_Pt",'GenPart_pt[muon2_idx]')
             .Define("jpsi_gen_Pt",'PT12(GenPart_pt,GenPart_eta,GenPart_phi,GenPart_mass, c1_idx,c2_idx)')
             .Define("minDRjpsiCharm","std::min(deltaR(GenPart_eta[c1_idx], GenPart_phi[c1_idx], Jpsi_kin_eta[0], Jpsi_kin_phi[0]),deltaR(GenPart_eta[c2_idx], GenPart_phi[c2_idx], Jpsi_kin_eta[0], Jpsi_kin_phi[0]))")
             .Define("index_charmCloseFar","genPartCloseFar(GenPart_pt,GenPart_eta,GenPart_phi,GenPart_mass,c1_idx, c2_idx, Jpsi_kin_pt,Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass)")
             .Filter("index_charmCloseFar[0]!= -1", "at least one match")
             .Define("charmClose","get_charmClose(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, index_charmCloseFar[0], muon1_idx, muon2_idx)")
             .Define("jetClose_charm_ptratio","index_charmCloseFar[0]!= -1 ? jetClose_Pt/charmClose.Pt(): -1")
             .Define("jetFar_charm_ptratio","index_charmCloseFar[1]!= -1 ? jetFar_Pt/GenPart_pt[index_charmCloseFar[1]]: -1".format(jetDef))
             #
        )

    df = callMVAclassification(df)
#    df = callMVARegression(df,"close")
#    df = callMVARegression(df,"far")

    if True:
        print("writing plots")
        hists = {
            # take these direction from the nano, no need to Define anything
            "Jpsi_kin_mass":  {"name":"Jpsi_kin_mass","title":"Jpsi mass;m_{#mu^{+}#mu^{-}} (GeV);N_{Events}","bin":100,"xmin":2.9,"xmax":3.4},
            "Jpsi_kin_massErr":  {"name":"Jpsi_kin_massErr","title":"Jpsi mass error;m_{#mu^{+}#mu^{-}} (GeV);N_{Events}","bin":90,"xmin":0.,"xmax":5.},
            "Jpsi_kin_pt":  {"name":"Jpsi_kin_pt","title":"Jpsi p_{T};p_{T}^{#mu^{+}#mu^{-}} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "Jpsi_kin_vtx_prob":  {"name":"Jpsi_kin_vtx_prob","title":"Jpsi vertex probability;N_{Events}","bin":100,"xmin":0.,"xmax":0.1},
            "Jpsi_muon1_pt":  {"name":"Jpsi_muon1_pt","title":"p_{T}(#mu);p_{T}(#mu) (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "Jpsi_muon2_pt":  {"name":"Jpsi_muon2_pt","title":"p_{T}(#mu);p_{T}(#mu) (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "Jpsi_kin_eta":  {"name":"Jpsi_kin_eta","title":"Jpsi #eta;#eta_{#mu^{+}#mu^{-}} (GeV);N_{Events}","bin":100,"xmin":-5.,"xmax":5.},
            "Jpsi_iso":  {"name":"Jpsi_iso","title":"Jpsi isolation;N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "Jpsi_kin_valid":  {"name":"Jpsi_kin_valid","title":"Jpsi kin valid idx;N_{Events}","bin":10,"xmin":0.,"xmax":10.},
            "Jpsi_isoPho":  {"name":"Jpsi_isoPho","title":"Jpsi isolation;N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            #        "Jpsi_neuHad":  {"name":"Jpsi_neuHad","title":"Jpsi isolation;N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "Jpsi_kin_lxy":  {"name":"Jpsi_kin_lxy","title":"vertex displacement in XY plane wrt Beam Spot; Jpsi lxy;N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "Jpsi_kin_sipPV":  {"name":"Jpsi_kin_sipPV","title":"impact parameter significance of the candidate trajectory in 3D wrt PV; sipPV^{#mu^{+}#mu^{-}} ;N_{Events}","bin":100,"xmin":0.,"xmax":5.},
            "Jpsi_muon1_isTighMuon":  {"name":"Jpsi_muon1_isTighMuon","title":"Jpsi muon1_isTighMuon;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            "Jpsi_muon1_isMediumMuon":  {"name":"Jpsi_muon1_isMediumMuon","title":"Jpsi muon1_isMediumMuon;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            "Jpsi_muon2_isTighMuon":  {"name":"Jpsi_muon2_isTighMuon","title":"Jpsi muon2_isTighMuon;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            "Jpsi_muon2_isMediumMuon":  {"name":"Jpsi_muon2_isMediumMuon","title":"Jpsi muon2_isMediumMuon;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            #
            "Jpsi_leadingChargedMaxPt1":  {"name":"Jpsi_leadingChargedMaxPt1","title":"Jpsi Pt leadingCharged within ConedR04; p_{T}","bin":50,"xmin":0.,"xmax":50.},
            "Jpsi_leadingCharged2dSignMaxPt1":  {"name":"Jpsi_leadingCharged2dSignMaxPt1","title":"Jpsi 2dSignMaxPt1 leadingCharged within ConedR04; p_{T}","bin":200,"xmin":-10,"xmax":10.},
            "Jpsi_leadingCharged3dSignMaxPt1":  {"name":"Jpsi_leadingCharged3dSignMaxPt1","title":"Jpsi 3dSignMaxPt1 leadingCharged within ConedR04; p_{T}","bin":200,"xmin":-10,"xmax":10.},
            "Jpsi_leadingChargedDxyMaxPt1":  {"name":"Jpsi_leadingChargedDxyMaxPt1","title":"Jpsi Dxy leadingCharged within ConedR04; d_{xy}","bin":200,"xmin":-0.1,"xmax":0.1},
            "Jpsi_leadingChargedDzMaxPt1":  {"name":"Jpsi_leadingChargedDzMaxPt1","title":"Jpsi Dz leadingCharged within ConedR04; d_{z}","bin":200,"xmin":-0.1,"xmax":0.1},
            "Jpsi_leadingChargedMaxPt2":  {"name":"Jpsi_leadingChargedMaxPt2","title":"Jpsi Pt subLeadingCharged within ConedR04; p_{T}","bin":50,"xmin":0.,"xmax":50.},
            "Jpsi_leadingCharged2dSignMaxPt2":  {"name":"Jpsi_leadingCharged2dSignMaxPt2","title":"Jpsi 2dSignMaxPt2 subLeadingCharged within ConedR04; p_{T}","bin":200,"xmin":-10,"xmax":10.},
            "Jpsi_leadingCharged3dSignMaxPt2":  {"name":"Jpsi_leadingCharged3dSignMaxPt2","title":"Jpsi 3dSignMaxPt2 subLeadingCharged within ConedR04; p_{T}","bin":200,"xmin":-10,"xmax":10.},
            "Jpsi_leadingChargedDxyMaxPt2":  {"name":"Jpsi_leadingChargedDxyMaxPt2","title":"Jpsi Dxy subLeadingCharged within ConedR04; d_{xy}","bin":200,"xmin":-0.1,"xmax":0.1},
            "Jpsi_leadingChargedDzMaxPt2":  {"name":"Jpsi_leadingChargedDzMaxPt2","title":"Jpsi Dz subLeadingCharged within ConedR04; d_{z}","bin":200,"xmin":-0.1,"xmax":0.1},
            ##
            "jet1Pt":  {"name":"jet1Pt","title":"Leading Jet p_{T}; AK4 jet p_{T} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "jet2Pt":  {"name":"jet2Pt","title":"Subleading Jet p_{T}; AK4 jet p_{T} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "jet1Eta":  {"name":"jet1Eta","title":"Leading Jet eta; AK4 jet #eta;N_{Events}","bin":60,"xmin":-3.,"xmax":3.},
            "jet2Eta":  {"name":"jet2Eta","title":"Subleading Jet eta; AK4 jet #eta;N_{Events}","bin":60,"xmin":-3.,"xmax":3.},
            "jet1CvL":  {"name":"jet1CvL","title":"Leading Jet CvL; btag CvL;N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "jet2CvL":  {"name":"jet2CvL","title":"Subleading Jet CvL; btag CvL;N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            #        "jet1nMuons":  {"name":"jet1nMuons","title":"Leading Jet nMuons; nMuons (GeV); N_{Events}","bin":10,"xmin":0.,"xmax":10.},
            #        "jet2nMuons":  {"name":"jet2nMuons","title":"Subleading Jet nMuons; nMuons (GeV); N_{Events}","bin":10,"xmin":0.,"xmax":10.},
#            "jetClose_nMuons":  {"name":"jetClosenMuons","title":"Close nMuons; nMuons; N_{Events}","bin":10,"xmin":0.,"xmax":10.},
#            "jetFar_nMuons":  {"name":"jetFarnMuons","title":"Far nMuons; nMuons; N_{Events}","bin":10,"xmin":0.,"xmax":10.},
#            "jetClose_nElectrons":  {"name":"jetClosenElectrons","title":"Close Electrons; nElectron; N_{Events}","bin":10,"xmin":0.,"xmax":10.},
#            "jetFar_nElectrons":  {"name":"jetFarnElectrons","title":"Far Electrons; nElectrons; N_{Events}","bin":10,"xmin":0.,"xmax":10.},
            #
            "jetFar_Pt":  {"name":"jetFar_Pt","title":"Jet Far Pt; p_{T}^{Far Jet} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "jetFar_Eta":  {"name":"jetFar_Eta","title":"Jet Far eta; #eta^{Far Jet};N_{Events}","bin":60,"xmin":-3.,"xmax":3.},
            "jetClose_Pt":  {"name":"jetClose_Pt","title":"Jet Close Pt; p_{T}^{Close Jet} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "jetClose_Eta":  {"name":"jetClose_Eta","title":"Jet Close eta; #eta^{Close Jet};N_{Events}","bin":60,"xmin":-3.,"xmax":3.},
            "jetFarcRegCorr":  {"name":"jetFarcRegCorr","title":"Jet Far value of the regression; cRegCorr;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            "jetClosecRegCorr":  {"name":"jetClosecRegCorr","title":"Jet Close value of the regression; cRegCorr;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            "jetFar_CvL":  {"name":"jetFar_CvL","title":"Jet Far; btagCvL^{Far Jet};N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "jetClose_CvL":  {"name":"jetClose_CvL","title":"Jet Close; btagCvL^{Close Jet};N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "jetFar_CvB":  {"name":"jetFar_CvB","title":"Jet Far; btagCvB^{Far Jet} ;N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "jetClose_CvB":  {"name":"jetClose_CvB","title":"Jet Close; btagCvB^{Close Jet};N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "jetFar_nConst":  {"name":"jetFar_nConst","title":"Jet Far; n constituent Far;N_{Events}","bin":20,"xmin":0.,"xmax":40.},
            "jetClose_nConst":  {"name":"jetClose_nConst","title":"Jet Close; n constituent Close;N_{Events}","bin":20,"xmin":0.,"xmax":40.},
            "jetCloseJPsiRatio":  {"name":"jetCloseJPsiRatio","title":"Jet Close; p_{T}(JPsi)/p_{T}^{Close Jet};N_{Events}","bin":100,"xmin":0.,"xmax":5.},
            "jetFar_PTMHRatio":  {"name":"jetFar_PTMHRatio","title":"Jet Far; p_{T}^{Jet Close}/MH ;N_{Events}","bin":100,"xmin":0.,"xmax":5.},
            #        "Jet_btagDeepCvL":  {"name":"Jet_btagDeepCvL","title":"btagDeepCvL; FlavCvL} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            #        "Jet_btagDeepFlavCvL":  {"name":"Jet_btagDeepFlavCvL","title":"btagDeepFlavCvL; DeepFlavCvL} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            ##
            "mJJ": {"name":"mJJ","title":"M(jet,jet); M(jet,jet) (GeV);N_{Events}","bin":400,"xmin":50.,"xmax":450},
            "massHiggs": {"name":"massHiggs","title":"M(jet,jet,jpsi); M(jet,jet,jpsi) (GeV);N_{Events}","bin":200,"xmin":50.,"xmax":450.},
            "massHiggsCorr": {"name":"massHiggsCorr","title":"M(jet,jet) [far corrected for cReg and close for 0.88 flat correction]; M(jet,jet) (GeV);N_{Events}","bin":200,"xmin":50.,"xmax":450.},
            "minDRjpsi": {"name":"minDRjpsi","title":"minDR(jet,jpsi); minDR(jet,jpsi) rad;N_{Events}","bin":60,"xmin":.0,"xmax":6},
            ##
#            "HLT_Dimuon20_Jpsi_Barrel_Seagulls":  {"name":"HLT_Dimuon20_Jpsi_Barrel_Seagulls","title":"HLT_Dimuon20_Jpsi_Barrel_Seagulls;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
#            "HLT_Dimuon25_Jpsi":  {"name":"HLT_Dimuon25_Jpsi","title":"HLT_Dimuon25_Jpsi;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
#            "triggerAna":  {"name":"triggerAna","title":"HLT_Dimuon20_Jpsi_Barrel_Seagulls or HLT_Dimuon25_Jpsi;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            "angle_phi1":  {"name":"angle_phi1","title":"angle_phi1; angle_phi1;N_{Events}","bin":100,"xmin":0.,"xmax":3.15},
            "angle_phi2":  {"name":"angle_phi2","title":"angle_phi2; angle_phi2;N_{Events}","bin":100,"xmin":0.,"xmax":3.15},
            "angle_phi3":  {"name":"angle_phi3","title":"angle_phi3; angle_phi3;N_{Events}","bin":100,"xmin":0.,"xmax":3.15},
            "angle_phi4":  {"name":"angle_phi4","title":"angle_phi4; angle_phi4;N_{Events}","bin":100,"xmin":0.,"xmax":3.15},
            "angle_cos_theta_jpsi":  {"name":"angle_cos_theta_jpsi","title":"angle_cos_theta_jpsi; angle_cos_theta_jpsi; N_{Events}","bin":100,"xmin":-1.,"xmax":1.},
            "angle_cos_theta_dicharm":  {"name":"angle_cos_theta_dicharm","title":"angle_cos_theta_dicharm; angle_cos_theta_dicharm; N_{Events}","bin":100,"xmin":-1.,"xmax":1.},
            "angle_cos_delta_m_min":  {"name":"angle_cos_delta_m_min","title":"angle_cos_delta_m_min; angle_cos_delta_m_min; N_{Events}","bin":100,"xmin":-1.,"xmax":1.},
            "discrMVA":  {"name":"discrMVA","title":"discrMVA; discrMVA; N_{Events}","bin":100,"xmin":0.,"xmax":1.},
        }

        hists2 = {
            "jet1partonFlavour":  {"name":"jet1partonFlavour","title":"Leading Jet partonFlavour; partonFlavour; N_{Events}","bin":25,"xmin":0.,"xmax":25.},
            "jet2partonFlavour":  {"name":"jet2partonFlavour","title":"Subleading Jet partonFlavour; partonFlavour; N_{Events}","bin":25,"xmin":0.,"xmax":25.},
            "jetFar_partonFlavour":  {"name":"jetFar_partonFlavour","title":"Jet Far; AK4 partonFlavour Far;N_{Events}","bin":25,"xmin":0.,"xmax":25.},
            "jetClose_partonFlavour":  {"name":"jetClose_partonFlavour","title":"Jet Close; AK4 partonFlavour Close;N_{Events}","bin":25,"xmin":0.,"xmax":25.},
            "jetFarScale":  {"name":"jetFarScale","title":"Jet Far; p_{T}^{RECO} / p_{T}^{GEN} Far;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            "jetCloseScale":  {"name":"jetCloseScale","title":"Jet Close; p_{T}^{RECO} / p_{T}^{GEN} Close;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            "jetFarPtCharm":  {"name":"jetFarPtCharm","title":"Jet Far; p_{T} Far (charm) (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "jetFarPtGluon":  {"name":"jetFarPtGluon","title":"Jet Far; p_{T} Far (gluon) (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "jetFarCvLCharm":  {"name":"jetFarCvLCharm","title":"Jet Far; btag CvL Far -- charm (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "jetFarCvLGluon":  {"name":"jetFarCvLGluon","title":"Jet Far; btag CvL Far -- gluon (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "jetFarcRegCorrScale":  {"name":"jetFarcRegCorrScale","title":"Jet Far; p_{T}^{RECO} * cRegCorr / p_{T}^{GEN} Far;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            "jetClosecRegCorrScale":  {"name":"jetClosecRegCorrScale","title":"Jet Close; p_{T}^{RECO} * cRegCorr / p_{T}^{GEN} Close;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
        }

        if mc>0: hists.update(hists2)

        hists3 = {
            "jetClose_charm_ptratio":  {"name":"jetClose_charm_ptratio","title":"Jet Close; p_{T}^{RECO} / p_{T}^{charm} Close;N_{Events}","bin":100,"xmin":0.,"xmax":10.},
            "jetFar_charm_ptratio":  {"name":"jetFar_charm_ptratio","title":"Jet Far; p_{T}^{RECO} / p_{T}^{charm} Far;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            "minDRjpsiCharm": {"name":"minDRjpsiCharm","title":"minDR(charm{1,2},jpsi); minDR(charm{1,2},jpsi) rad;N_{Events}","bin":60,"xmin":.0,"xmax":6},
            "charm1_gen_Pt":  {"name":"charm1_gen_Pt","title":"charm1 p_{T}; charm1 jet p_{T} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "charm2_gen_Pt":  {"name":"charm2_gen_Pt","title":"chamr2 p_{T}; charm2 p_{T} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "muon1_gen_Pt":  {"name":"muon1_gen_Pt","title":"muon1 p_{T}; muon1 jet p_{T} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "muon2_gen_Pt":  {"name":"muon2_gen_Pt","title":"muon2 p_{T}; muon2 p_{T} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "jpsi_gen_Pt":  {"name":"jpsi_gen_Pt","title":"jpsi p_{T}; jpsi p_{T} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
        }
        if mc==1000 or mc==1001: hists.update(hists3)

        histos = []

        for h in hists:
            model1d = (hists[h]["name"]+"_"+str(year), hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"])
            h1d = df.Histo1D(model1d, hists[h]["name"],"w")
            histos.append(h1d)
#            print("h1d append",h1d.GetName(),' integral',h1d.Integral())

        if False and (mc==1000 or mc==1001):
            for h in histos:
                canv = ROOT.TCanvas("stackcanvas","Stack canvas",800,800)
                #        canv.SetLogy(1)
                h.Draw("hist")

                canv.Draw()
                #        canv.SaveAs("~/public_html/Hrare_Jpsi/"+h.GetName()+".png")
                if mc==1000: canv.SaveAs("~/public_html/Hrare_JpsiSept25/"+h.GetName()+"_H.png")
                if mc==1001: canv.SaveAs("~/public_html/Hrare_JpsiSept25/"+h.GetName()+"_Z.png")

        if True:
            myDir='/work/submit/mariadlf/Hrare_JPsiCC/Sept25/'
            outputFileHisto = myDir+"histoOUTname_"+str(mc)+"_"+str(year)+".root"
            print(outputFileHisto )
            myfile = ROOT.TFile(outputFileHisto,"RECREATE")

            myfile.ls()
            for h in histos:
                h.Write()

            myfile.Close()
            myfile.Write()

        if False:
            branchList = ROOT.vector('string')()
            for branchName in [
                    "mc",
                    "w",
                    "PV_npvsGood",
                    "triggerAna",
                    "run",
                    "event",
                    "luminosityBlock",
                    #
                    "jetFar_Pt",
                    "jetFar_Eta",
                    "jetFar_Phi",
                    "jetFar_Mass",
                    "jetClose_Pt",
                    "jetClose_Eta",
                    "jetClose_Phi",
                    "jetClose_Mass",
                    "Jpsi_kin_pt",
                    "Jpsi_kin_eta",
                    "Jpsi_kin_phi",
                    "Jpsi_kin_mass",
                    "Jpsi_muon1_pt",
                    "Jpsi_muon2_pt",
                    "Jpsi_muon1_eta",
                    "Jpsi_muon2_eta",
                    "Jpsi_muon1_phi",
                    "Jpsi_muon2_phi",
                    # for the classification
                    "nGoodJets",
                    "jetClose_CvL",
                    "jetFar_CvL",
                    "jetClose_CvB",
                    "jetFar_CvB",
                    "jetClose_nConst",
                    "jetFar_nConst",
                    "jetFar_chEmEF",
                    "jetFar_chHEF",
                    "jetFar_neEmEF",
                    "jetFar_neHEF",
                    "jetFar_muEF",
                    "jetClose_chEmEF",
                    "jetClose_chHEF",
                    "jetClose_neEmEF",
                    "jetClose_neHEF",
                    "jetClose_muEF",
                    "jetCloseJPsiRatio",
                    "jetFar_PTMHRatio",
                    "jpsi_Higgs_Ratio",
                    #
                    "Jpsi_kin_sipPV",
                    "Jpsi_iso",
                    "Jpsi_kin_lxy",
                    "Jpsi_leadingCharged3dSignMaxPt1",
                    "Jpsi_leadingChargedMaxPt1",
                    "Jpsi_leadingCharged3dSignMaxPt2",
                    "Jpsi_leadingChargedMaxPt2",
                    # event like quantities
                    "massHiggsCorr",
                    # angles
                    "angle_phi1",
                    "angle_phi2",
                    "angle_phi3",
                    "angle_phi4",
                    "angle_cos_theta_jpsi",
                    "angle_cos_theta_dicharm",
                    "angle_cos_delta_m_min",
                    #
                    "MET_pt",
                    "MET_phi",
                    "discrMVA",
#                    "regr_far",
#                    "regr_close"
            ]:
                branchList.push_back(branchName)

            if mc>0:
                # for the regression
                for branchName in [
                        "GenPart_pt",
                        "GenPart_eta",
                        "GenPart_phi",
                        "GenPart_mass",
                        "GenPart_pdgId",
                        "GenPart_genPartIdxMother",
                        "jetClose_partonFlavour",
                        "jetFar_partonFlavour",
                        "jetClose_charm_ptratio",
                        "jetFar_charm_ptratio",
                        "charm1_gen_Pt",
                        "charm2_gen_Pt",
                        "muon1_gen_Pt",
                        "muon2_gen_Pt",
                        "jpsi_gen_Pt",
                ]:
                    branchList.push_back(branchName)

            myDir='/work/submit/mariadlf/Hrare_JPsiCC/Sept25/'
            outputFile = myDir+"snapshotJpsiCC_"+str(mc)+"_"+str(year)+".root"
            print(outputFile)
            snapshotOptions = ROOT.RDF.RSnapshotOptions()
            snapshotOptions.fCompressionAlgorithm = ROOT.kLZ4
            snapshot_tdf = df.Snapshot("events", outputFile, branchList, snapshotOptions)
            print("snapshot_tdf DONE")
            print(outputFile)

            now = datetime.now()
            print('==> ends: ',now)


        if False:
            branchListTraining = ROOT.vector('string')()

            for branchNameTraining in [
                    "mc",
                    "w",
                    "PV_npvsGood",
                    "triggerAna",
                    "run",
                    "event",
                    "luminosityBlock",
                    #
                    "pfcand_close_pt",
                    "pfcand_close_eta",
                    "pfcand_close_phi",
                    "pfcand_close_dxy",
                    "pfcand_close_dz",
                    "jet_npfcand_close",
                    #
                    "pfcand_close_isMu",
                    "pfcand_close_isEl",
                    "pfcand_close_isChargedHad",
                    "pfcand_close_isGamma",
                    "pfcand_close_isNeutralHad",
                    #
                    "pfcand_close_btagEtaRel",
                    "pfcand_close_btagPtRatio",
                    "pfcand_close_btagPParRatio",
                    "pfcand_close_btagSip3dSig",
                    "pfcand_close_btagSip3dVal",
                    "pfcand_close_btagSip2dSig",
                    "pfcand_close_btagSip2dVal",
                    "pfcand_close_btagJetDistVal",
                    "pfcand_close_ptrel_log",
                    "pfcand_close_phirel",
                    "pfcand_close_etarel",
                    "pfcand_close_charge",
                    "pfcand_close_mass",
                    #
                    "recojet_close_isSig",
                    "recojet_close_isB"
                ]:
                    branchListTraining.push_back(branchNameTraining)

            myDir='/work/submit/mariadlf/Hrare_JPsiCC/JAN2025/'
            outputFile = myDir+"snapshotJpsiCC_"+str(mc)+"_"+str(year)+"_TRAINING.root"
            print(outputFile)
            snapshotOptions = ROOT.RDF.RSnapshotOptions()
            snapshotOptions.fCompressionAlgorithm = ROOT.kLZ4
            snapshot_tdf = df.Snapshot("events", outputFile, branchListTraining, snapshotOptions)
            print("snapshot_tdf DONE")
            print(outputFile)

            now = datetime.now()
            print('==> ends: ',now)

def loopOnDataset(year):

    thisdict = BuildDictJpsiCC(year)

    mc = []
    if year=="2018": mc = [1000]
#    if year=="2017" or year=="12016" or year=="22016": mc = [10,11]
#    if year=="2018": mc = [10,11,1000,1001]
#    if year=="12022": mc = [12]

    for sampleNOW in mc:
        files = SwitchSample(thisdict,sampleNOW)[0]
        print('outside the function: ', len(files))
        rdf = ROOT.RDataFrame("Runs", files) # make sure this is not the distributed
        sumW = computeWeigths(rdf,SwitchSample(thisdict,sampleNOW)[1])
        analysis(files,year,sampleNOW,sumW)

    data = []
    '''
    if year=="12016": data = [-1,-2,-3,-4,-5,-6]
    if year=="22016": data = [-6,-7,-8]
    if year=="2017": data = [-2,-3,-4,-5,-6]
    if year=="2018": data = [-1,-2,-3,-4]
    if year=="12022": data = [-13,-14]
    if year=="22022": data = [-15,-16,-17]
    if year=="12023": data = [-21,-22,-23,-24]
    if year=="22023": data = [-31,-32]
    '''

    readDataQuality(year)

    for sampleNOW in data:
        files = SwitchSample(thisdict,sampleNOW)[0]
        print('outside the function: ', len(files))
        analysis(files,year,sampleNOW,1.)

if __name__ == "__main__":

    now = datetime.now()
    print('==> very beginning: ',now)

    print('year=',sys.argv[1])
    year=sys.argv[1]

    loopOnDataset(year)
