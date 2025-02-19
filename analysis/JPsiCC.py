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

tagger='Deep'
tagger='DeepFlav'
#if year=='12022': tagger='PNet'

#regress='cRegCorr'
#regress='PNetRegPtRawCorrNeutrino'
regress='Jet_cRegCorr[goodJets][index_CloseFar[1]]'
regressClose='Jet_cRegCorr[goodJets][index_CloseFar[0]]'

regress='1.0'
regressClose='1.0'

GOODJETS = "(Jet_pt>30 and abs(Jet_eta)<2.5 and Jet_btag{}CvL>-1)".format(tagger)


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

    dfComm = (dfINI
              .Define("mc","{}".format(mc))
              .Define("isData","{}".format(isData))
#              .Define("applyJson","{}".format(JSON)).Filter("applyJson","pass JSON") # doesn't work for signal ??
              .Define("w","{}".format(weight))
              .Define("wraw","{}".format(weight))
              .Define("lumiIntegrated","{}".format(lumiIntegrated))
              .Filter("PV_npvsGood>0","one good PV")

              )

    GOODJPSI = "Jpsi_kin_pt>25"
    TRIGGER=''
    if year=='2018' or year=='12022' or year=='22022' or year=='12023' or year=='22023': TRIGGER="(HLT_Dimuon25_Jpsi or HLT_Dimuon25_Jpsi_noCorrL1)"
    if year=='2017': TRIGGER="HLT_Dimuon25_Jpsi"
    if year=='12016' or year=='22016':
        TRIGGER="(HLT_Dimuon20_Jpsi or HLT_Dimuon16_Jpsi)"
        GOODJPSI = "Jpsi_kin_pt>20"
    print(TRIGGER)


    df= (dfComm
         .Filter("nJpsi>0","at least one Jpsi;")
         .Define("goodMeson","{}".format(GOODJPSI))
         .Filter("Sum(goodMeson)>0", "one good Jpsi")
         .Define("triggerAna","{}".format(TRIGGER))
         .Filter("triggerAna>0","passing trigger")
         .Define("goodJets","{}".format(GOODJETS))
         .Define("nGoodJets","Sum(goodJets)*1.0f").Filter("Sum(goodJets)>1","two jets for cc")
         .Define("mJJ","Minv(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets])")
         .Define("massHiggs","Minv3massive(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets], Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass)")
         #
         .Define("jet1Pt","Jet_pt[goodJets][0]")
         .Define("jet2Pt","Jet_pt[goodJets][1]")
         .Define("jet1Eta","Jet_eta[goodJets][0]")
         .Define("jet2Eta","Jet_eta[goodJets][1]")
         .Define("jet1CvL","Jet_btag{}CvL[goodJets][0]".format(tagger))
         .Define("jet2CvL","Jet_btag{}CvL[goodJets][1]".format(tagger))
         #
         .Define("index_CloseFar","jetCloseFar(Jet_pt[goodJets],Jet_eta[goodJets],Jet_phi[goodJets],Jet_mass[goodJets],Jpsi_kin_pt,Jpsi_kin_eta, Jpsi_kin_phi,Jpsi_kin_mass)")
         .Filter("index_CloseFar[0]!= -1", "at least one close jet")
         # for the FAR we apply JEC and c-Regression
         # for the CLOSE no JEC and no c-Regression
         .Define("jetFarPtRaw","Jet_pt[goodJets][index_CloseFar[1]]")
         .Define("jetFar_Pt","Jet_pt[goodJets][index_CloseFar[1]]*{}".format(regress))
         .Filter("jetFar_Pt>30","30 GeV threshould for jetFar")
         .Define("jetFar_Eta","Jet_eta[goodJets][index_CloseFar[1]]")
         .Define("jetFar_Phi","Jet_phi[goodJets][index_CloseFar[1]]")
         .Define("jetFar_Mass","Jet_mass[goodJets][index_CloseFar[1]]")
#         .Define("jetClosePt","buildCloseJets(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets], Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass, index_CloseFar[0],0)")
         .Define("jetClose_Pt","buildCloseJets(Jet_pt[goodJets]*(1-Jet_rawFactor[goodJets]), Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets], Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass, index_CloseFar[0],0)")
         .Define("jetClose_Eta","buildCloseJets(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets], Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass, index_CloseFar[0],1)")
         .Define("jetClose_Phi","buildCloseJets(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets], Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass, index_CloseFar[0],2)")
         .Define("jetClose_Mass","buildCloseJets(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets], Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass, index_CloseFar[0],3)")
         #
         .Define("jetClosecRegCorr","{}".format(regressClose))
         .Define("jetFarcRegCorr","{}".format(regress))
         #
         .Define("jetCloseJPsiRatio","Jpsi_kin_pt[0]/jetClose_Pt")
         .Define("jetClose_CvL","Jet_btag{}CvL[goodJets][index_CloseFar[0]]".format(tagger))
         .Define("jetFar_CvL","Jet_btag{}CvL[goodJets][index_CloseFar[1]]".format(tagger))
         .Define("jetClose_nConst","Jet_nConstituents[goodJets][index_CloseFar[0]]-2.0")
         .Define("jetFar_nConst","Jet_nConstituents[goodJets][index_CloseFar[1]]")
         .Define("jetClose_nMuons","Jet_nMuons[goodJets][index_CloseFar[0]]")
         .Define("jetFar_nMuons","Jet_nMuons[goodJets][index_CloseFar[1]]")
         .Define("jetClose_nElectrons","Jet_nElectrons[goodJets][index_CloseFar[0]]")
         .Define("jetFar_nElectrons","Jet_nElectrons[goodJets][index_CloseFar[1]]")
         .Define("jetFar_chEmEF","Jet_chEmEF[goodJets][index_CloseFar[1]]")
         .Define("jetFar_chHEF","Jet_chHEF[goodJets][index_CloseFar[1]]")
         .Define("jetFar_neEmEF","Jet_neEmEF[goodJets][index_CloseFar[1]]")
         .Define("jetFar_neHEF","Jet_neHEF[goodJets][index_CloseFar[1]]")
         .Define("jetFar_muEF","Jet_muEF[goodJets][index_CloseFar[1]]")
         .Define("jetClose_chEmEF","Jet_chEmEF[goodJets][index_CloseFar[0]]")
         .Define("jetClose_chHEF","Jet_chHEF[goodJets][index_CloseFar[0]]")
         .Define("jetClose_neEmEF","Jet_neEmEF[goodJets][index_CloseFar[0]]")
         .Define("jetClose_neHEF","Jet_neHEF[goodJets][index_CloseFar[0]]")
         .Define("jetClose_muEF","Jet_muEF[goodJets][index_CloseFar[0]]")
         #
         .Define("massHiggsCorr","Minv3body(jetClose_Pt,jetClose_Eta,jetClose_Phi,jetClose_Mass,jetFar_Pt,jetFar_Eta,jetFar_Phi,jetFar_Mass,Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass,0 )")
         .Define("HiggsPt","Minv3body(jetClose_Pt,jetClose_Eta,jetClose_Phi,jetClose_Mass,jetFar_Pt,jetFar_Eta,jetFar_Phi,jetFar_Mass,Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass, 1)")
         #
         .Define("jetFar_PTMHRatio","jetFar_Pt/massHiggsCorr")
         .Define("jpsi_Higgs_Ratio","Jpsi_kin_pt[0]/HiggsPt")
         #
         .Define("minDRjpsi","std::min(deltaR(Jet_eta[goodJets][0], Jet_phi[goodJets][0], Jpsi_kin_eta[0], Jpsi_kin_phi[0]),deltaR(Jet_eta[goodJets][1], Jet_phi[goodJets][1], Jpsi_kin_eta[0], Jpsi_kin_phi[0]))")
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

    if mc>0:
        df= (df.Define("jet1partonFlavour","abs(Jet_partonFlavour[goodJets][0])")
         .Define("jet2partonFlavour","abs(Jet_partonFlavour[goodJets][1])")
         .Define("jetFarPtCharm","(abs(Jet_partonFlavour[goodJets][index_CloseFar[1]])==4) ? {} * Jet_pt[goodJets][index_CloseFar[1]]:-1.".format(regress))
         .Define("jetFarPtGluon","(abs(Jet_partonFlavour[goodJets][index_CloseFar[1]])==21) ? {} * Jet_pt[goodJets][index_CloseFar[1]]:-1.".format(regress))
         .Define("jetFarCvLCharm","(abs(Jet_partonFlavour[goodJets][index_CloseFar[1]])==4) ? Jet_btag{}CvL[goodJets][index_CloseFar[1]]:-1.".format(tagger))
         .Define("jetFarCvLGluon","(abs(Jet_partonFlavour[goodJets][index_CloseFar[1]])==21) ? Jet_btag{}CvL[goodJets][index_CloseFar[1]]:-1.".format(tagger))
         .Define("jetClose_partonFlavour","abs(Jet_partonFlavour[goodJets][index_CloseFar[0]])")
         .Define("jetFar_partonFlavour","abs(Jet_partonFlavour[goodJets][index_CloseFar[1]])")
         .Define("jetCloseScale","Jet_pt[goodJets][index_CloseFar[0]]/GenJet_pt[Jet_genJetIdx[goodJets][index_CloseFar[0]]]")
         .Define("jetFarScale","Jet_pt[goodJets][index_CloseFar[1]]/GenJet_pt[Jet_genJetIdx[goodJets][index_CloseFar[1]]]")
         .Define("jetClosecRegCorrScale","Jet_pt[goodJets][index_CloseFar[0]]*{}/GenJet_pt[Jet_genJetIdx[goodJets][index_CloseFar[0]]]".format(regressClose))
         .Define("jetFarcRegCorrScale","Jet_pt[goodJets][index_CloseFar[1]]*{}/GenJet_pt[Jet_genJetIdx[goodJets][index_CloseFar[1]]]".format(regress))
         #
         )

    if mc==1000 or mc==1001:
        df= (df.Define("idxs", "genMatch_indices(GenPart_pdgId, GenPart_genPartIdxMother)")
             .Define("muon1_idx", 'idxs.at("muon1_idx")')
             .Define("muon2_idx", 'idxs.at("muon2_idx")')
             .Define("neg_muon_idx", 'idxs.at("neg_muon_idx")')
             .Define("pos_muon_idx", 'idxs.at("pos_muon_idx")')
             .Define("c1_idx", 'idxs.at("c1_idx")')
             .Define("c2_idx", 'idxs.at("c2_idx")')
             .Define("higgs_idx", 'idxs.at("higgs_idx")')
             .Define("minDRjpsiCharm","std::min(deltaR(GenPart_eta[c1_idx], GenPart_phi[c1_idx], Jpsi_kin_eta[0], Jpsi_kin_phi[0]),deltaR(GenPart_eta[c2_idx], GenPart_phi[c2_idx], Jpsi_kin_eta[0], Jpsi_kin_phi[0]))")
             .Define("index_charmCloseFar","genPartCloseFar(GenPart_pt,GenPart_eta,GenPart_phi,GenPart_mass,c1_idx, c2_idx, Jpsi_kin_pt,Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass)")
             .Filter("index_charmCloseFar[0]!= -1", "at least one match")
             .Define("charmClose","get_charmClose(GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, index_charmCloseFar[0], muon1_idx, muon2_idx)")
             .Define("jetClose_charm_ptratio","index_charmCloseFar[0]!= -1 ? Jet_pt[goodJets][index_CloseFar[0]]/charmClose.Pt(): -1")
             .Define("jetFar_charm_ptratio","index_charmCloseFar[1]!= -1 ?Jet_pt[goodJets][index_CloseFar[1]]/GenPart_pt[index_charmCloseFar[1]]: -1")
             #
        )

    df = callMVAclassification(df)
#    df = callMVARegression(df,"close")
#    df = callMVARegression(df,"far")

    if True:
        print("writing plots")
        hists = {
            # take these direction from the nano, no need to Define anything
            "Jpsi_kin_mass":  {"name":"Jpsi_kin_mass","title":"Jpsi mass;m_{#mu^{+}#mu^{-}} (GeV);N_{Events}","bin":120,"xmin":2.5,"xmax":3.7},
            "Jpsi_kin_massErr":  {"name":"Jpsi_kin_massErr","title":"Jpsi mass error;m_{#mu^{+}#mu^{-}} (GeV);N_{Events}","bin":90,"xmin":0.,"xmax":5.},
            "Jpsi_kin_pt":  {"name":"Jpsi_kin_pt","title":"Jpsi p_{T};p_{T} {#mu^{+}#mu^{-}} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "Jpsi_kin_vtx_prob":  {"name":"Jpsi_kin_vtx_prob","title":"Jpsi vertex probability;N_{Events}","bin":100,"xmin":0.,"xmax":0.1},
            "Jpsi_muon1_pt":  {"name":"Jpsi_muon1_pt","title":"p_{T}(#mu);p_{T}(#mu) (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "Jpsi_muon2_pt":  {"name":"Jpsi_muon2_pt","title":"p_{T}(#mu);p_{T}(#mu)  (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "Jpsi_kin_eta":  {"name":"Jpsi_kin_eta","title":"Jpsi #eta;#eta_{#mu^{+}#mu^{-}} (GeV);N_{Events}","bin":100,"xmin":-5.,"xmax":5.},
            "Jpsi_iso":  {"name":"Jpsi_iso","title":"Jpsi isolation;N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "Jpsi_kin_valid":  {"name":"Jpsi_kin_valid","title":"Jpsi kin valid idx;N_{Events}","bin":10,"xmin":0.,"xmax":10.},
            "Jpsi_isoPho":  {"name":"Jpsi_isoPho","title":"Jpsi isolation;N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            #        "Jpsi_neuHad":  {"name":"Jpsi_neuHad","title":"Jpsi isolation;N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "Jpsi_kin_lxy":  {"name":"Jpsi_kin_lxy","title":"Jpsi kin_lxy;N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "Jpsi_kin_sipPV":  {"name":"Jpsi_kin_sipPV","title":"Jpsi kin_sipPV;N_{Events}","bin":100,"xmin":0.,"xmax":5.},
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
            "Jpsi_leadingCharged2dSignMaxPt2":  {"name":"Jpsi_leadingCharged2dSignMaxPt2","title":"Jpsi 2dSignMaxPt1 subLeadingCharged within ConedR04; p_{T}","bin":200,"xmin":-10,"xmax":10.},
            "Jpsi_leadingCharged3dSignMaxPt2":  {"name":"Jpsi_leadingCharged3dSignMaxPt2","title":"Jpsi 3dSignMaxPt1 subLeadingCharged within ConedR04; p_{T}","bin":200,"xmin":-10,"xmax":10.},
            "Jpsi_leadingChargedDxyMaxPt2":  {"name":"Jpsi_leadingChargedDxyMaxPt2","title":"Jpsi Dxy subLeadingCharged within ConedR04; d_{xy}","bin":200,"xmin":-0.1,"xmax":0.1},
            "Jpsi_leadingChargedDzMaxPt2":  {"name":"Jpsi_leadingChargedDzMaxPt2","title":"Jpsi Dz subLeadingCharged within ConedR04; d_{z}","bin":200,"xmin":-0.1,"xmax":0.1},
            ##
            "jet1Pt":  {"name":"jet1Pt","title":"Leading Jet p_{T}; AK4 jet p_{T}} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "jet2Pt":  {"name":"jet2Pt","title":"Subleading Jet p_{T}; AK4 jet p_{T}} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "jet1Eta":  {"name":"jet1Eta","title":"Leading Jet eta; AK4 jet #eta} (GeV);N_{Events}","bin":100,"xmin":-5.,"xmax":5.},
            "jet2Eta":  {"name":"jet2Eta","title":"Subleading Jet eta; AK4 jet #eta (GeV);N_{Events}","bin":100,"xmin":-5.,"xmax":5.},
            "jet1CvL":  {"name":"jet1CvL","title":"Leading Jet CvL; AK4 btag CvL (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "jet2CvL":  {"name":"jet2CvL","title":"Subleading Jet CvL; AK4 btag CvL (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            #        "jet1nMuons":  {"name":"jet1nMuons","title":"Leading Jet nMuons; nMuons (GeV); N_{Events}","bin":10,"xmin":0.,"xmax":10.},
            #        "jet2nMuons":  {"name":"jet2nMuons","title":"Subleading Jet nMuons; nMuons (GeV); N_{Events}","bin":10,"xmin":0.,"xmax":10.},
#            "jetClose_nMuons":  {"name":"jetClosenMuons","title":"Close nMuons; nMuons; N_{Events}","bin":10,"xmin":0.,"xmax":10.},
#            "jetFar_nMuons":  {"name":"jetFarnMuons","title":"Far nMuons; nMuons; N_{Events}","bin":10,"xmin":0.,"xmax":10.},
#            "jetClose_nElectrons":  {"name":"jetClosenElectrons","title":"Close Electrons; nElectron; N_{Events}","bin":10,"xmin":0.,"xmax":10.},
#            "jetFar_nElectrons":  {"name":"jetFarnElectrons","title":"Far Electrons; nElectrons; N_{Events}","bin":10,"xmin":0.,"xmax":10.},
            #
            "jetFar_Pt":  {"name":"jetFar_Pt","title":"Jet Far Pt; AK4 pt Far (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "jetFar_Eta":  {"name":"jetFar_Eta","title":"Jet Far eta; AK4 jet #eta (GeV);N_{Events}","bin":100,"xmin":-5.,"xmax":5.},
            "jetClose_Pt":  {"name":"jetClose_Pt","title":"Jet Close Pt; AK4 pt Close (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "jetClose_Eta":  {"name":"jetClose_Eta","title":"Jet Close eta; AK4 jet #eta} (GeV);N_{Events}","bin":100,"xmin":-5.,"xmax":5.},
            "jetFarcRegCorr":  {"name":"jetFarcRegCorr","title":"Jet Far value of the regression; cRegCorr ;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            "jetClosecRegCorr":  {"name":"jetClosecRegCorr","title":"Jet Close value of the regression; cRegCorr;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            "jetFar_CvL":  {"name":"jetFar_CvL","title":"Jet Far; AK4 btag CvL Far (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "jetClose_CvL":  {"name":"jetClose_CvL","title":"Jet Close; AK4 btag CvL Close (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "jetFar_nConst":  {"name":"jetFar_nConst","title":"Jet Far; AK4 n constituent Far;N_{Events}","bin":20,"xmin":0.,"xmax":40.},
            "jetClose_nConst":  {"name":"jetClose_nConst","title":"Jet Close; AK4 n constituent Close;N_{Events}","bin":20,"xmin":0.,"xmax":40.},
            "jetCloseJPsiRatio":  {"name":"jetCloseJPsiRatio","title":"Jet Close; pt(JPsi) / AK4 pt(Jet);N_{Events}","bin":100,"xmin":0.,"xmax":5.},
            "jetFar_PTMHRatio":  {"name":"jetFar_PTMHRatio","title":"Jet Far; M-higgs / AK4 pt(Jet);N_{Events}","bin":100,"xmin":0.,"xmax":5.},
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
            "jetFarScale":  {"name":"jetFarScale","title":"Jet Far; AK4 pt / AK4 Gen pt Far;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            "jetCloseScale":  {"name":"jetCloseScale","title":"Jet Close; AK4 pt / AK4 Gen pt Close;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            "jetFarPtCharm":  {"name":"jetFarPtCharm","title":"Jet Far; AK4 pt Far (charm)(GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "jetFarPtGluon":  {"name":"jetFarPtGluon","title":"Jet Far; AK4 pt Far (gluon) (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
            "jetFarCvLCharm":  {"name":"jetFarCvLCharm","title":"Jet Far; AK4 btag CvL Far -- charm (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "jetFarCvLGluon":  {"name":"jetFarCvLGluon","title":"Jet Far; AK4 btag CvL Far -- gluon (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
            "jetFarcRegCorrScale":  {"name":"jetFarcRegCorrScale","title":"Jet Far; AK4 pt * cRegCorr / AK4 Gen pt Far;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            "jetClosecRegCorrScale":  {"name":"jetClosecRegCorrScale","title":"Jet Close; AK4 pt * cRegCorr / AK4 Gen pt Close;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
        }

        if mc>0: hists.update(hists2)

        hists3 = {
            "jetClose_charm_ptratio":  {"name":"jetClose_charm_ptratio","title":"Jet Close; AK4 pt / charm Pt Close;N_{Events}","bin":100,"xmin":0.,"xmax":10.},
            "jetFar_charm_ptratio":  {"name":"jetFar_charm_ptratio","title":"Jet Far; AK4 pt / charm Pt Far;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            "minDRjpsiCharm": {"name":"minDRjpsiCharm","title":"minDR(charm{1,2},jpsi); minDR(charm{1,2},jpsi) rad;N_{Events}","bin":60,"xmin":.0,"xmax":6},
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

def loopOnDataset(year):

    thisdict = BuildDictJpsiCC(year)

    mc = []
    if year=="2017" or year=="12016" or year=="22016": mc = [10,11]
    if year=="2018": mc = [10,11,1000,1001]
    if year=="12022": mc = [12]

    for sampleNOW in mc:
        files = SwitchSample(thisdict,sampleNOW)[0]
        print('outside the function: ', len(files))
        rdf = ROOT.RDataFrame("Runs", files) # make sure this is not the distributed
        sumW = computeWeigths(rdf,SwitchSample(thisdict,sampleNOW)[1])
        analysis(files,year,sampleNOW,sumW)

    data = []
    if year=="12016": data = [-1,-2,-3,-4,-5,-6]
    if year=="22016": data = [-6,-7,-8]
    if year=="2017": data = [-2,-3,-4,-5,-6]
    if year=="2018": data = [-1,-2,-3,-4]
    if year=="12022": data = [-13,-14]
    if year=="22022": data = [-15,-16,-17]
    if year=="12023": data = [-21,-22,-23,-24]
    if year=="22023": data = [-31,-32]

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
