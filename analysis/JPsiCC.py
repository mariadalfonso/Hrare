import ROOT
from utilsAna import loadUserCode
loadUserCode()


def computeWeigths(chain2):

    rdf = ROOT.RDataFrame(chain2)

    genEventSumWeight = rdf.Sum("genEventSumw").GetValue()
    genEventSumNoWeight = rdf.Sum("genEventCount").GetValue()

    print('genEventSumWeight',genEventSumWeight)
    print('genEventSumNoWeight',genEventSumNoWeight)

    # for now neglecting the lumi
    weight = (1./genEventSumWeight)
    print('weight',weight)

ROOT.gROOT.SetBatch()
ROOT.ROOT.EnableImplicitMT()
chain = ROOT.TChain("Events")
chain.Add("/data/submit/cms/store/user/mariadlf/nano/GluGluH_HJPsiCC/NANOAOD_test3/*.root")
chainR = ROOT.TChain("Runs")
chainR.Add("/data/submit/cms/store/user/mariadlf/nano/GluGluH_HJPsiCC/NANOAOD_test3/*.root")

sumw = computeWeigths(chainR)
print('sumw',sumw)
lumi = 1.
weight = "{0}*genWeight*{1}".format(lumi,10000000)
#weight = "{0}*genWeight".format(lumi)

year = 2018
dfINI = ROOT.RDataFrame(chain)
TRIGGER = "HLT_Dimuon20_Jpsi_Barrel_Seagulls or HLT_Dimuon25_Jpsi"
GOODJETS = "(Jet_pt>20 and abs(Jet_eta)<2.5 and Jet_btagDeepCvL>-1)"

df= (dfINI.Filter("nJpsi>0","at least one Jpsi")
     .Define("w","{}".format(weight))
     .Define("triggerAna","{}".format(TRIGGER))
     .Define("goodJets","{}".format(GOODJETS))
     .Define("nGoodJets","Sum(goodJets)*1.0f").Filter("Sum(goodJets)>1","two jets for cc")
     .Define("mJJ","Minv(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets])")
     .Define("massHiggs","Minv3massive(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets], Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass)")
     .Define("ptHiggs","Ptinv3massive(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets], Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass)")     
     #
     .Define("jet1Pt","Jet_pt[goodJets][0]")
     .Define("jet2Pt","Jet_pt[goodJets][1]")
     .Define("jet1Eta","Jet_eta[goodJets][0]")
     .Define("jet2Eta","Jet_eta[goodJets][1]")
     .Define("jet1CvL","Jet_btagDeepCvL[goodJets][0]")
     .Define("jet2CvL","Jet_btagDeepCvL[goodJets][1]")
     .Define("jet1partonFlavour","abs(Jet_partonFlavour[goodJets][0])")
     .Define("jet2partonFlavour","abs(Jet_partonFlavour[goodJets][1])")
     #
     .Define("index_CloseFar","jetCloseFar(Jet_pt[goodJets],Jet_eta[goodJets],Jet_phi[goodJets],Jet_mass[goodJets],Jpsi_kin_pt,Jpsi_kin_eta, Jpsi_kin_phi,Jpsi_kin_mass)")
     .Filter("index_CloseFar[0]!= -1", "at least one close jet")
     .Define("jetCloseCvL","Jet_btagDeepCvL[goodJets][index_CloseFar[0]]")
     .Define("jetFarCvL","Jet_btagDeepCvL[goodJets][index_CloseFar[1]]")
     .Define("jetClosePt","Jet_pt[goodJets][index_CloseFar[0]]")
     .Define("jetFarPt","Jet_pt[goodJets][index_CloseFar[1]]")
     .Define("jetFarPtCharm","(abs(Jet_partonFlavour[goodJets][index_CloseFar[1]])==4) ? Jet_cRegCorr[goodJets][index_CloseFar[1]] * Jet_pt[goodJets][index_CloseFar[1]]:-1.")
     .Define("jetFarPtGluon","(abs(Jet_partonFlavour[goodJets][index_CloseFar[1]])==21) ? Jet_cRegCorr[goodJets][index_CloseFar[1]] * Jet_pt[goodJets][index_CloseFar[1]]:-1.")
     .Define("jetFarCvLCharm","(abs(Jet_partonFlavour[goodJets][index_CloseFar[1]])==4) ? Jet_btagDeepCvL[goodJets][index_CloseFar[1]]:-1.")
     .Define("jetFarCvLGluon","(abs(Jet_partonFlavour[goodJets][index_CloseFar[1]])==21) ? Jet_btagDeepCvL[goodJets][index_CloseFar[1]]:-1.")
     .Define("jetCloseScale","Jet_pt[goodJets][index_CloseFar[0]]/GenJet_pt[Jet_genJetIdx[goodJets][index_CloseFar[0]]]")
     .Define("jetFarScale","Jet_pt[goodJets][index_CloseFar[1]]/GenJet_pt[Jet_genJetIdx[goodJets][index_CloseFar[1]]]")
     .Define("jetClosecRegCorrScale","Jet_pt[goodJets][index_CloseFar[0]]*Jet_cRegCorr[goodJets][index_CloseFar[0]]/GenJet_pt[Jet_genJetIdx[goodJets][index_CloseFar[0]]]")
     .Define("jetFarcRegCorrScale","Jet_pt[goodJets][index_CloseFar[1]]*Jet_cRegCorr[goodJets][index_CloseFar[1]]/GenJet_pt[Jet_genJetIdx[goodJets][index_CloseFar[1]]]")     
     .Define("jetClosecRegCorr","Jet_cRegCorr[goodJets][index_CloseFar[0]]")
     .Define("jetFarcRegCorr","Jet_cRegCorr[goodJets][index_CloseFar[1]]")
     .Define("jetClosepartonFlavour","abs(Jet_partonFlavour[goodJets][index_CloseFar[0]])")
     .Define("jetFarpartonFlavour","abs(Jet_partonFlavour[goodJets][index_CloseFar[1]])")
     .Define("jetClosenConst","Jet_nConstituents[goodJets][index_CloseFar[0]]")
     .Define("jetFarnConst","Jet_nConstituents[goodJets][index_CloseFar[1]]")
     .Define("jetCloseJPsiRatio","Jpsi_kin_pt[0]/Jet_pt[goodJets][index_CloseFar[0]]")
     .Define("jetClosenMuons","Jet_nMuons[goodJets][index_CloseFar[0]]")
     .Define("jetFarnMuons","Jet_nMuons[goodJets][index_CloseFar[1]]")
     .Define("jetClosenElectrons","Jet_nElectrons[goodJets][index_CloseFar[0]]")
     .Define("jetFarnElectrons","Jet_nElectrons[goodJets][index_CloseFar[1]]")
     #
     .Define("massHiggsCorr","Minv3massiveCorr(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets], Jet_cRegCorr[goodJets], Jet_muonIdx1[goodJets][index_CloseFar[0]], Jet_muonIdx2[goodJets][index_CloseFar[0]], Muon_pt, Muon_eta, Muon_phi, Jpsi_muon1_pt, Jpsi_muon1_eta, Jpsi_muon1_phi, Jpsi_muon2_pt, Jpsi_muon2_eta, Jpsi_muon2_phi, Jpsi_kin_pt, Jpsi_kin_eta, Jpsi_kin_phi, Jpsi_kin_mass, index_CloseFar[0])")
     #     .Define("jetCloseCvL","deltaR(Jet_eta[goodJets][0], Jet_phi[goodJets][0], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) < deltaR(Jet_eta[goodJets][1], Jet_phi[goodJets][1], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) ? Jet_btagDeepCvL[goodJets][0] : Jet_btagDeepCvL[goodJets][1]")
#     .Define("jetFarCvL","deltaR(Jet_eta[goodJets][0], Jet_phi[goodJets][0], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) > deltaR(Jet_eta[goodJets][1], Jet_phi[goodJets][1], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) ? Jet_btagDeepCvL[goodJets][0] : Jet_btagDeepCvL[goodJets][1]")
#     .Define("jetClosePt","deltaR(Jet_eta[goodJets][0], Jet_phi[goodJets][0], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) < deltaR(Jet_eta[goodJets][1], Jet_phi[goodJets][1], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) ? Jet_pt[goodJets][0] : Jet_pt[goodJets][1]")
#     .Define("jetFarPt","deltaR(Jet_eta[goodJets][0], Jet_phi[goodJets][0], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) > deltaR(Jet_eta[goodJets][1], Jet_phi[goodJets][1], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) ? Jet_pt[goodJets][0] : Jet_pt[goodJets][1]")
#     .Define("jetClosepartonFlavour","deltaR(Jet_eta[goodJets][0], Jet_phi[goodJets][0], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) < deltaR(Jet_eta[goodJets][1], Jet_phi[goodJets][1], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) ? Jet_partonFlavour[goodJets][0] : Jet_partonFlavour[goodJets][1]")
#     .Define("jetFarpartonFlavour","deltaR(Jet_eta[goodJets][0], Jet_phi[goodJets][0], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) > deltaR(Jet_eta[goodJets][1], Jet_phi[goodJets][1], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) ? Jet_partonFlavour[goodJets][0] : Jet_partonFlavour[goodJets][1]")
#     .Define("jetClosenConst","deltaR(Jet_eta[goodJets][0], Jet_phi[goodJets][0], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) < deltaR(Jet_eta[goodJets][1], Jet_phi[goodJets][1], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) ? Jet_nConstituents[goodJets][0] : Jet_nConstituents[goodJets][1]")
#     .Define("jetFarnConst","deltaR(Jet_eta[goodJets][0], Jet_phi[goodJets][0], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) > deltaR(Jet_eta[goodJets][1], Jet_phi[goodJets][1], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) ? Jet_nConstituents[goodJets][0] : Jet_nConstituents[goodJets][1]")
#     .Define("jetCloseJPsiRatio","deltaR(Jet_eta[goodJets][0], Jet_phi[goodJets][0], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) < deltaR(Jet_eta[goodJets][1], Jet_phi[goodJets][1], Jpsi_kin_eta[0], Jpsi_kin_phi[0]) ? Jpsi_kin_pt[0]/Jet_pt[goodJets][0] : Jpsi_kin_pt[0]/Jet_pt[goodJets][1]")
     #
     .Define("minDRjpsi","std::min(deltaR(Jet_eta[goodJets][0], Jet_phi[goodJets][0], Jpsi_kin_eta[0], Jpsi_kin_phi[0]),deltaR(Jet_eta[goodJets][1], Jet_phi[goodJets][1], Jpsi_kin_eta[0], Jpsi_kin_phi[0]))")
     .Filter("minDRjpsi<0.3","jPsi very close to 1 charm-jet")
#     .Filter("abs(Jpsi_kin_eta[0])>1.4","barrel JPsi")
#     .Filter("abs(Jpsi_leadingChargedDxy[0])>0.5","bad track")
     )

count = df.Count().GetValue()
print('nevents=',count)

if True:
    print("writing plots")
    hists = {
        # take these direction from the nano, no need to Define anything
        "Jpsi_kin_mass":  {"name":"Jpsi_kin_mass","title":"Jpsi mass;m_{#mu^{+}#mu^{-}} (GeV);N_{Events}","bin":90,"xmin":1.5,"xmax":4.5},
        "Jpsi_kin_pt":  {"name":"Jpsi_kin_pt","title":"Jpsi p_{T};p_{T} {#mu^{+}#mu^{-}} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},        
        "Jpsi_muon1_pt":  {"name":"Jpsi_muon1_pt","title":"p_{T}(#mu);p_{T}(#mu) (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
        "Jpsi_muon2_pt":  {"name":"Jpsi_muon2_pt","title":"p_{T}(#mu);p_{T}(#mu)  (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
        "Jpsi_kin_eta":  {"name":"Jpsi_kin_eta","title":"Jpsi #eta;#eta_{#mu^{+}#mu^{-}} (GeV);N_{Events}","bin":100,"xmin":-5.,"xmax":5.},
        "Jpsi_iso":  {"name":"Jpsi_iso","title":"Jpsi isolation;N_{Events}","bin":100,"xmin":0.,"xmax":1.},
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
        "jet2Eta":  {"name":"jet2Eta","title":"Subleading Jet eta; AK4 jet #eta } (GeV);N_{Events}","bin":100,"xmin":-5.,"xmax":5.},        
        "jet1CvL":  {"name":"jet1CvL","title":"Leading Jet CvL; AK4 btagDeepCvL (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
        "jet2CvL":  {"name":"jet2CvL","title":"Subleading Jet CvL; AK4 btagDeepCvL (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
        "jet1partonFlavour":  {"name":"jet1partonFlavour","title":"Leading Jet partonFlavour; partonFlavour; N_{Events}","bin":25,"xmin":0.,"xmax":25.},
        "jet2partonFlavour":  {"name":"jet2partonFlavour","title":"Subleading Jet partonFlavour; partonFlavour; N_{Events}","bin":25,"xmin":0.,"xmax":25.},
#        "jet1nMuons":  {"name":"jet1nMuons","title":"Leading Jet nMuons; nMuons (GeV); N_{Events}","bin":10,"xmin":0.,"xmax":10.},
#        "jet2nMuons":  {"name":"jet2nMuons","title":"Subleading Jet nMuons; nMuons (GeV); N_{Events}","bin":10,"xmin":0.,"xmax":10.},
        "jetClosenMuons":  {"name":"jetClosenMuons","title":"Close nMuons; nMuons; N_{Events}","bin":10,"xmin":0.,"xmax":10.},
        "jetFarnMuons":  {"name":"jetFarnMuons","title":"Far nMuons; nMuons; N_{Events}","bin":10,"xmin":0.,"xmax":10.},
        "jetClosenElectrons":  {"name":"jetClosenElectrons","title":"Close Electrons; nElectron; N_{Events}","bin":10,"xmin":0.,"xmax":10.},
        "jetFarnElectrons":  {"name":"jetFarnElectrons","title":"Far Electrons; nElectrons; N_{Events}","bin":10,"xmin":0.,"xmax":10.},
        #
        "jetFarPtCharm":  {"name":"jetFarPtCharm","title":"Jet Far; AK4 pt Far (charm)(GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
        "jetFarPtGluon":  {"name":"jetFarPtGluon","title":"Jet Far; AK4 pt Far (gluon) (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
        "jetFarPt":  {"name":"jetFarPt","title":"Jet Far; AK4 pt Far (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
        "jetClosePt":  {"name":"jetClosePt","title":"Jet Close; AK4 pt Close (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":100.},
        "jetFarcRegCorrScale":  {"name":"jetFarcRegCorrScale","title":"Jet Far; AK4 pt * cRegCorr / AK4 Gen pt Far;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
        "jetClosecRegCorrScale":  {"name":"jetClosecRegCorrScale","title":"Jet Close; AK4 pt * cRegCorr / AK4 Gen pt Close;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
        "jetFarScale":  {"name":"jetFarScale","title":"Jet Far; AK4 pt / AK4 Gen pt Far;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
        "jetCloseScale":  {"name":"jetCloseScale","title":"Jet Close; AK4 pt / AK4 Gen pt Close;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
        "jetFarcRegCorr":  {"name":"jetFarcRegCorr","title":"Jet Far; cRegCorr ;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
        "jetClosecRegCorr":  {"name":"jetClosecRegCorr","title":"cRegCorr;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
        "jetFarCvL":  {"name":"jetFarCvL","title":"Jet Far; AK4 btagDeepCvL Far (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
        "jetCloseCvL":  {"name":"jetCloseCvL","title":"Jet Close; AK4 btagDeepCvL Close (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
        "jetFarCvLCharm":  {"name":"jetFarCvLCharm","title":"Jet Far; AK4 btagDeepCvL Far -- charm (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
        "jetFarCvLGluon":  {"name":"jetFarCvLGluon","title":"Jet Far; AK4 btagDeepCvL Far -- gluon (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
        "jetFarpartonFlavour":  {"name":"jetFarpartonFlavour","title":"Jet Far; AK4 partonFlavour Far;N_{Events}","bin":25,"xmin":0.,"xmax":25.},
        "jetClosepartonFlavour":  {"name":"jetClosepartonFlavour","title":"Jet Close; AK4 partonFlavour Close;N_{Events}","bin":25,"xmin":0.,"xmax":25.},
        "jetFarnConst":  {"name":"jetFarnConst","title":"Jet Far; AK4 n constituent Far;N_{Events}","bin":20,"xmin":0.,"xmax":40.},
        "jetClosenConst":  {"name":"jetClosenConst","title":"Jet Close; AK4 n constituent Close;N_{Events}","bin":20,"xmin":0.,"xmax":40.},
        "jetCloseJPsiRatio":  {"name":"jetCloseJPsiRatio","title":"Jet Close; pt(JPsi) / AK4 pt(Jet);N_{Events}","bin":40,"xmin":0.,"xmax":2.},
#        "Jet_btagDeepCvL":  {"name":"Jet_btagDeepCvL","title":"btagDeepCvL; FlavCvL} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
#        "Jet_btagDeepFlavCvL":  {"name":"Jet_btagDeepFlavCvL","title":"btagDeepFlavCvL; DeepFlavCvL} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":1.},
        ##
        "mJJ": {"name":"mJJ","title":"M(jet,jet); M(jet,jet) (GeV);N_{Events}","bin":400,"xmin":50.,"xmax":450},
        "massHiggs": {"name":"massHiggs","title":"M(jet,jet,jpsi); M(jet,jet,jpsi) (GeV);N_{Events}","bin":400,"xmin":50.,"xmax":450.},
        "massHiggsCorr": {"name":"massHiggsCorr","title":"M(jet,jet) [far corrected for cReg and close for 0.88 flat correction]; M(jet,jet) (GeV);N_{Events}","bin":400,"xmin":50.,"xmax":450.},
        "minDRjpsi": {"name":"minDRjpsi","title":"minDR(jet,jpsi); minDR(jet,jpsi) (GeV);N_{Events}","bin":100,"xmin":.0,"xmax":10},        
        "ptHiggs": {"name":"ptHiggs","title":"pt(jet,jet,jpsi); pt(jet,jet,jpsi) (GeV);N_{Events}","bin":200,"xmin":0.,"xmax":200},        

        ##        
        "HLT_Dimuon20_Jpsi_Barrel_Seagulls":  {"name":"HLT_Dimuon20_Jpsi_Barrel_Seagulls","title":"HLT_Dimuon20_Jpsi_Barrel_Seagulls;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
        "HLT_Dimuon25_Jpsi":  {"name":"HLT_Dimuon25_Jpsi","title":"HLT_Dimuon25_Jpsi;N_{Events}","bin":100,"xmin":0.,"xmax":2.},                        
        "triggerAna":  {"name":"triggerAna","title":"HLT_Dimuon20_Jpsi_Barrel_Seagulls or HLT_Dimuon25_Jpsi;N_{Events}","bin":100,"xmin":0.,"xmax":2.},                
    }

    histos = []

    for h in hists:
        model1d = (hists[h]["name"]+"_"+str(year), hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"])
        h1d = df.Histo1D(model1d, hists[h]["name"],"w")
        histos.append(h1d)
        print("h1d append")

    for h in histos:

        canv = ROOT.TCanvas("stackcanvas","Stack canvas",800,800)
#        canv.SetLogy(1)
        h.Draw("hist")

        canv.Draw()
#        canv.SaveAs("~/public_html/Hrare_Jpsi/"+h.GetName()+".png")
        canv.SaveAs("~/public_html/Hrare_Jpsi_w_noVTX/"+h.GetName()+".png")

        
if False:

    branchList = ROOT.vector('string')()
    for branchName in [
            "Muon_pt",
    ]:
        branchList.push_back(branchName)
        
        outputFile = "DASKlogs/snapshotOUTname.root"
        print(outputFile)
    snapshotOptions = ROOT.RDF.RSnapshotOptions()
    snapshotOptions.fCompressionAlgorithm = ROOT.kLZ4
    snapshot_tdf = df.Snapshot("events", outputFile, branchList, snapshotOptions)
    print("snapshot_tdf DONE")
    print(outputFile)

    now = datetime.now()
    print('==> ends: ',now)
