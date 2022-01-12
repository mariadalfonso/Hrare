import ROOT
import os
import sys

ROOT.ROOT.EnableImplicitMT()
from utilsHrare import getMClist, getDATAlist
from utilsHrare import plot
from utilsHrare import SwitchSample

doPlot = False

isZ = False
isW = False
isVBF = False

isPhiCat = "false"
isRhoCat = "false"

if sys.argv[1]=='isZtag': isZ = True
if sys.argv[1]=='isWtag': isW = True
if sys.argv[1]=='isVBFtag': isVBF = True

if sys.argv[2]=='isPhiCat': isPhiCat = "true"
if sys.argv[2]=='isRhoCat': isRhoCat = "true"

if isVBF:
    BARRELphotons = "(Photon_pt>75 and Photon_isScEtaEB and (Photon_cutBased & 2) and Photon_electronVeto)"
    ENDCAPphotons = ""
    GOODPHI = "(abs(phi_kin_mass-1.02)<0.01 && phi_kin_pt>5 && phi_trk1_pt>5 && phi_trk2_pt>5 && phi_kin_valid && phi_iso > 0.8)"
    GOODRHO = "(rho_kin_mass>0.6 && rho_kin_mass<0.95 && rho_kin_pt>20 && rho_trk1_pt>5 && rho_trk2_pt>5 && rho_kin_valid && rho_iso > 0.9)"
else:
    BARRELphotons = "Photon_pt>20 and Photon_isScEtaEB and (Photon_cutBased & 2) and Photon_electronVeto"
    ENDCAPphotons = "Photon_pt>20 and Photon_isScEtaEE and (Photon_cutBased & 2) and Photon_electronVeto"
    GOODPHI = "(abs(phi_kin_mass-1.02)<0.01 && phi_kin_pt>10 && phi_trk1_pt>10 && phi_trk2_pt>10 && phi_kin_valid && phi_iso > 0.8)"
    GOODRHO = "(rho_kin_mass>0.6 && rho_kin_mass<0.95 && rho_kin_pt>10 && rho_kin_valid && rho_iso > 0.8)"

GOODJETS = "(Jet_pt>30 && abs(Jet_eta)<5 and (Jet_jetId & 2) and jet_mask and jet_mask2)" ## to add PUid

GOODMUON = "(Muon_pt>20 and abs(Muon_eta)<2.4 and Muon_isGlobal and Muon_isTracker and abs(Muon_dz)<0.10 and abs(Muon_dxy) < 0.05)" # ID and Iso will be asked later depending on W and Z
VETOelectrons = "(Electron_pt>10 and abs(Electron_eta) < 2.5 and Electron_pfRelIso03_all < 0.25 and Electron_mvaFall17V2noIso_WPL)"
LOOSEmuons = "(Muon_pt>10 and abs(Muon_eta)<2.4 and Muon_isGlobal and Muon_isTracker and Muon_pfRelIso04_all < 0.25 and abs(Muon_dz)<0.10 and abs(Muon_dxy) < 0.05 and Muon_looseId)"

TRIGGERincl = "HLT_Photon35_TwoProngs35"
TRIGGERvbf = "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3 or HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3 or HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50"

TRIGGERsMUO = "HLT_IsoMu24"
JSON = "isGoodRunLS(isData, run, luminosityBlock)"

def selectionTAG(df):

    if isZ:
        dftag = (df.Define("goodMuons","{}".format(GOODMUON)+" and Muon_mediumId and Muon_pfRelIso04_all < 0.25")
                 .Define("vetoMuons","{}".format(LOOSEmuons))
                 .Define("vetoElectrons","{}".format(VETOelectrons))
                 .Filter("Sum(goodMuons) >= 2 and Sum(Muon_charge[goodMuons])==0 and Sum(vetoElectrons)==0 and Sum(vetoMuons)==2", "At least two good OS muon, and no extra loose leptons")
                 .Define("V_mass", "Minv(Muon_pt[goodMuons], Muon_eta[goodMuons], Muon_phi[goodMuons], Muon_mass[goodMuons])")
                 .Filter("V_mass>(91-10) and V_mass<(91+15)","At least one good Z")
                 .Define("trigger","{}".format(TRIGGERsMUO))
                 .Define("Mu1_hasTriggerMatch", "hasTriggerMatch(Muon_eta[goodMuons][0], Muon_phi[goodMuons][0], TrigObj_eta, TrigObj_phi)")
                 .Define("Mu2_hasTriggerMatch", "hasTriggerMatch(Muon_eta[goodMuons][1], Muon_phi[goodMuons][1], TrigObj_eta, TrigObj_phi)")
                 .Filter("trigger>0 and ((Mu1_hasTriggerMatch and Muon_pt[goodMuons][0]>26) or (Mu2_hasTriggerMatch and Muon_pt[goodMuons][1]>26))","pass trigger")
        )
        return dftag

    if isW:
        dftag = (df.Define("goodMuons","{}".format(GOODMUON)+" and Muon_tightId and Muon_pfRelIso04_all < 0.15")
                 .Define("vetoElectrons","{}".format(VETOelectrons))
                 .Define("vetoMuons","{}".format(LOOSEmuons))
                 .Filter("Sum(goodMuons) == 1 and MET_pt>20 and Sum(vetoElectrons)==0 and Sum(vetoMuons)==1", "Exactly 1 good muon and MET>20 and no loose leptons")
                 .Define("V_mass","mt(Muon_pt[goodMuons][0], Muon_phi[goodMuons][0], MET_pt, MET_phi)")
                 .Filter("V_mass>20","MT>20")
                 .Define("trigger","{}".format(TRIGGERsMUO))
                 .Define("Mu1_hasTriggerMatch", "hasTriggerMatch(Muon_eta[goodMuons][0], Muon_phi[goodMuons][0], TrigObj_eta, TrigObj_phi)")
                 .Filter("trigger>0 and Mu1_hasTriggerMatch and Muon_pt[goodMuons][0]>26","pass trigger")
        )
        return dftag

    if isVBF:
        dftag = (df.Define("goodJets","{}".format(GOODJETS))
                 .Define("nGoodJets","Sum(goodJets)").Filter("Sum(goodJets)>1","two jets")
                 .Define("mJJ","Minv(Jet_pt[goodJets], Jet_eta[goodJets], Jet_phi[goodJets], Jet_mass[goodJets])")
                 .Define("dEtaJJ","abs(Jet_eta[goodJets][0] - Jet_eta[goodJets][1])")
                 .Define("Y1Y2","Jet_eta[goodJets][0]*Jet_eta[goodJets][1]")
                 .Filter("mJJ>300 and dEtaJJ>3","Filter on MJJ>300 , Deta>3")
                 .Define("vetoElectrons","{}".format(VETOelectrons))
                 .Define("vetoMuons","{}".format(LOOSEmuons))
                 .Filter("Sum(vetoElectrons)==0 and Sum(vetoMuons)==0", "no leptons")
                 .Define("trigger","{}".format(TRIGGERvbf))
                 .Filter("trigger>0", "pass triggers")
        )
    return dftag


def dfGammaMeson(df):

    #    dfbase = (df.Define("goodPhotons", "{}".format(BARRELphotons)+" or {}".format(ENDCAPphotons) )
    dfa= (df.Define("goodPhotons", "{}".format(BARRELphotons))
          .Define("goodPhotons_pt", "Photon_pt[goodPhotons]")
          .Define("goodPhotons_eta", "Photon_eta[goodPhotons]")
          .Define("goodPhotons_phi", "Photon_phi[goodPhotons]")
          .Define("goodPhotons_pfRelIso03_all", "Photon_pfRelIso03_all[goodPhotons]")
          .Define("goodPhotons_hoe", "Photon_hoe[goodPhotons]")
          .Define("goodPhotons_r9", "Photon_r9[goodPhotons]")
          .Define("goodPhotons_sieie", "Photon_sieie[goodPhotons]")
          #
          .Define("jet_mask", "cleaningMask(Photon_jetIdx[goodPhotons],nJet)")
          )
    dfb= (dfa.Define("goodPhi","({}".format(GOODPHI)+" && {}".format(isPhiCat)+")")
          .Define("goodPhi_pt", "phi_kin_pt[goodPhi]")
          .Define("goodPhi_eta", "phi_kin_eta[goodPhi]")
          .Define("goodPhi_phi", "phi_kin_phi[goodPhi]")
          .Define("goodPhi_mass", "phi_kin_mass[goodPhi]")
          .Define("goodPhi_iso", "phi_iso[goodPhi]")
          .Define("goodPhi_vtx_chi2dof", "phi_kin_vtx_chi2dof[goodPhi]")
          .Define("goodPhi_trk1_pt", "phi_trk1_pt[goodPhi]")
          .Define("goodPhi_trk2_pt", "phi_trk2_pt[goodPhi]")
          .Define("goodPhi_trk1_eta", "phi_trk1_eta[goodPhi]")
          .Define("goodPhi_trk2_eta", "phi_trk2_eta[goodPhi]")
          .Define("goodPhiDR","DeltaR(phi_trk1_eta[goodPhi],phi_trk2_eta[goodPhi],phi_trk1_phi[goodPhi],phi_trk2_phi[goodPhi])")
          )

    dfOBJ= (dfb.Define("goodRho","({}".format(GOODRHO)+" && {}".format(isRhoCat)+")")
            .Define("goodRho_pt", "rho_kin_pt[goodRho]")
            .Define("goodRho_eta", "rho_kin_eta[goodRho]")
            .Define("goodRho_phi", "rho_kin_phi[goodRho]")
            .Define("goodRho_iso", "rho_iso[goodRho]")
            .Define("goodRho_mass", "rho_kin_mass[goodRho]")
            .Define("goodRho_vtx_chi2dof", "rho_kin_vtx_chi2dof[goodRho]")
            .Define("goodRho_trk1_pt", "rho_trk1_pt[goodRho]")
            .Define("goodRho_trk2_pt", "rho_trk2_pt[goodRho]")
            .Define("goodRho_trk1_eta", "rho_trk1_eta[goodRho]")
            .Define("goodRho_trk2_eta", "rho_trk2_eta[goodRho]")
            .Define("goodRhoDR","DeltaR(rho_trk1_eta[goodRho],rho_trk2_eta[goodRho],rho_trk1_phi[goodRho],rho_trk2_phi[goodRho])")
        )

    return dfOBJ

def dfHiggsCand(df):

    if(isPhiCat=="true"):

        dfbase = (df.Filter("Sum(goodPhotons)>0", "At least one good Photon")
                  .Filter("Sum(goodPhi)>0", "one good Phi (ptPhi, validfit, ptTracks)")
                  .Define("index_pair","HiggsCandFromRECO(goodPhi_pt,goodPhi_eta,goodPhi_phi,goodPhi_mass,goodPhotons_pt,goodPhotons_eta,goodPhotons_phi)").Filter("index_pair[0]!= -1", "at least a good meson candidate")
                  .Define("jet_mask2", "cleaningJetFromMeson(Jet_eta, Jet_phi, goodPhi_eta[index_pair[0]], goodPhi_phi[index_pair[0]])")
                  .Define("HCandMass", "Minv2(goodPhi_pt[index_pair[0]],goodPhi_eta[index_pair[0]],goodPhi_phi[index_pair[0]],goodPhi_mass[index_pair[0]],goodPhotons_pt[index_pair[1]],goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]]).first")
                  .Define("HCandPT",   "Minv2(goodPhi_pt[index_pair[0]],goodPhi_eta[index_pair[0]],goodPhi_phi[index_pair[0]],goodPhi_mass[index_pair[0]],goodPhotons_pt[index_pair[1]],goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]]).second")
              )
        return dfbase


    if(isRhoCat=="true"):

        dfbase = (df.Filter("Sum(goodPhotons)>0", "At least one good Photon")
                  .Filter("Sum(goodRho)>0", "one good Rho (ptPhi, validfit, ptTracks)")
                  .Define("index_pair","HiggsCandFromRECO(goodRho_pt,goodRho_eta,goodRho_phi,goodRho_mass,goodPhotons_pt,goodPhotons_eta,goodPhotons_phi)").Filter("index_pair[0]!= -1", "at least a good meson candidate")
                  .Define("jet_mask2", "cleaningJetFromMeson(Jet_eta, Jet_phi, goodRho_eta[index_pair[0]], goodRho_phi[index_pair[0]])")
                  .Define("HCandMass", "Minv2(goodRho_pt[index_pair[0]],goodRho_eta[index_pair[0]],goodRho_phi[index_pair[0]],goodRho_mass[index_pair[0]],goodPhotons_pt[index_pair[1]],goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]]).first")
                  .Define("HCandPT",   "Minv2(goodRho_pt[index_pair[0]],goodRho_eta[index_pair[0]],goodRho_phi[index_pair[0]],goodRho_mass[index_pair[0]],goodPhotons_pt[index_pair[1]],goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]]).second")
                  )

        return dfbase


def analysis(df,mc,w,isData):

    dfOBJ= dfGammaMeson(df)
    dfbase = dfHiggsCand(dfOBJ)
    dfcandtag = selectionTAG(dfbase)
    dfFINAL = (dfcandtag.Define("w","{}".format(w))
               .Define("mc","{}".format(mc))
               .Define("isData","{}".format(isData))
               .Define("applyJson","{}".format(JSON)).Filter("applyJson","pass JSON")
    )


    branchList = ROOT.vector('string')()
    for branchName in [
            "HCandMass",
            "HCandPT",
            "index_pair",
            #
            "goodPhotons_pt",
            "goodPhotons_eta",
            "goodPhotons_pfRelIso03_all",
            "goodPhotons_hoe",
            "goodPhotons_r9",
            "goodPhotons_sieie",
            #
            "trigger",
            "SoftActivityJetNjets5",
            "MET_pt",
            #
            "w",
            "mc",
    ]:
        branchList.push_back(branchName)

    if(isPhiCat=="true"):
        for branchName in [
            "goodPhi",
            "goodPhiDR",
            "goodPhi_mass",
            "goodPhi_pt",
            "goodPhi_iso",
            "goodPhi_trk1_pt",
            "goodPhi_trk2_pt",
            "goodPhi_trk1_eta",
            "goodPhi_trk2_eta",
            "goodPhi_vtx_chi2dof",
        ]:
            branchList.push_back(branchName)

    if(isRhoCat=="true"):
        for branchName in [
            "goodRho",
            "goodRhoDR",
            "goodRho_mass",
            "goodRho_pt",
            "goodRho_iso",
            "goodRho_trk1_pt",
            "goodRho_trk2_pt",
            "goodRho_trk1_eta",
            "goodRho_trk2_eta",
            "goodRho_vtx_chi2dof",
        ]:
            branchList.push_back(branchName)

    if isZ or isW:
        for branchName in [
                "V_mass",
        ]:
            branchList.push_back(branchName)

    if isVBF:
        for branchName in [
                "mJJ",
                "nGoodJets",
                "dEtaJJ",
                "Y1Y2",
        ]:
            branchList.push_back(branchName)

    outputFile = "outname_mc%d"%mc+".root"
    cat = ""
    if(isPhiCat=="true"): cat = "_PhiCat"
    if(isRhoCat=="true"): cat = "_RhoCat"
    if isZ : outputFile = "outname_mc%d"%mc+"_Zcat"+cat+".root"
    if isW : outputFile = "outname_mc%d"%mc+"_Wcat"+cat+".root"
    if isVBF: outputFile = "outname_mc%d"%mc+"_VBFcat"+cat+".root"

    snapshot_tdf = dfFINAL.Snapshot("events", outputFile, branchList)
    print("snapshot_tdf DONE")
    print("*** SUMMARY :")
    print(outputFile)

    print("---------------- SUMMARY -------------")
    ## this doens't work with the negative weights
    report = dfcandtag.Report()
    report.Print()

    print("--------------------------------------")

    if doPlot:
        hists = {
            #        "Z_mass":     {"name":"Z_mass","title":"Di Muon mass; m_{#mu^{+}#mu^{-}} (GeV);N_{Events}","bin":500,"xmin":70,"xmax":120},
            "V_mass":     {"name":"V_mass","title":"transverse mass; m_{T}(#mu^{+} MET} (GeV);N_{Events}","bin":80,"xmin":40,"xmax":120},
            "HCandMass":  {"name":"HCandMass","title":"H mass;m_{k^{+}k^{-}#gamma} (GeV);N_{Events}","bin":500,"xmin":100,"xmax":150},
            "phi_num":    {"name":"nphi","title":"Phi N;N {k^{+}k^{-}} (GeV);N_{Events}","bin":10,"xmin":0.,"xmax":10.},
            "Phi_mass":   {"name":"phi_kin_mass","title":"Phi mass;m_{k^{+}k^{-}} (GeV);N_{Events}","bin":200,"xmin":0.95,"xmax":1.15},
            "Phi_pt":     {"name":"phi_kin_pt","title":"Phi pt ;p^{T}_{k^{+}k^{-}} (GeV);N_{Events}","bin":1000,"xmin":0.25,"xmax":50.25},
            "Phi_gen_mass":   {"name":"phi_gen_mass","title":"Phi gen mass;m_{k^{+}k^{-}} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":10.},        
            "Phi_mass_err":   {"name":"phi_kin_massErr","title":"Phi mass error;m_{k^{+}k^{-}} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":0.5},
            "Phi_kin_vtx_chi2dof":   {"name":"phi_kin_vtx_chi2dof","title":"Phi vtx_chi2dof;m_{k^{+}k^{-}} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":5.0},
        }
        
        for h in hists:
            model = (hists[h]["name"], hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"])
            h = snapshot_tdf.Histo1D(model, hists[h]["name"])        
            plot(h,"PLOTS/plot_"+h.GetName()+".png",False,2)

def readMCSample(sampleNOW):

    files = getMClist(sampleNOW)
    print(len(files))
    df = ROOT.RDataFrame("Events", files)

    nevents = df.Count().GetValue()  ## later with negative weights
    w = (SwitchSample(sampleNOW)[1] / nevents)

    lumiEq = (nevents / SwitchSample(sampleNOW)[1])
    print("%s entries in the dataset" %nevents)
    print("lumi equivalent fb %s" %lumiEq)
    analysis(df,sampleNOW,w,"false")


def readDataSample(year,type):

    files = getDATAlist(year,type)

    df = ROOT.RDataFrame("Events", files)

    w=1.
    nevents = df.Count().GetValue()
    print("%s entries in the dataset" %nevents)

    analysis(df,type,w,"true")


def runTest():

    df = ROOT.RDataFrame("Events", "root://eoscms.cern.ch//eos/cms//store/group/phys_higgs/HiggsExo/dalfonso/Hrare/D01/vbf-hphigamma-powheg/NANOAOD_01/step7_VBS_Phigamma_8.root")

    w=1.
    nevents = df.Count().GetValue()
    print("%s entries in the dataset" %nevents)

    sampleNOW=-1
    analysis(df,-1,w,"false")

   
if __name__ == "__main__":

#    runTest()
#    to run: python3 -i VGammaMeson_cat.py isVBFtag isPhiCat 12
#    print(int(sys.argv[3]))

    if(int(sys.argv[3]) < 0): readDataSample(2018,int(sys.argv[3]) )  # SingleMuon
    else: readMCSample(int(sys.argv[3])) # to switch sample
