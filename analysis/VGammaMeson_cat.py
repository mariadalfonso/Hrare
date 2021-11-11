import ROOT
import os

ROOT.ROOT.EnableImplicitMT()
from utilsHrare import getMClist, getDATAlist
from utilsHrare import plot
from utilsHrare import SwitchSample

doPlot = False
lumi = 137.
isZtag = False;
isWtag = True;


BARRELphotons = "Photon_pt>20 and Photon_isScEtaEB and (Photon_cutBased & 2) and Photon_electronVeto"
ENDCAPphotons = "Photon_pt>20 and Photon_isScEtaEE and (Photon_cutBased & 2) and Photon_electronVeto"

GOODPHI = "(phi_kin_pt>10 && phi_trk1_pt>10 && phi_trk2_pt>10 && phi_kin_massErr<0.025 && phi_kin_valid)"

GOODMUON = "(Muon_pt>20 and abs(Muon_eta)<2.4 and Muon_isGlobal and Muon_isTracker and abs(Muon_dz)<0.10 and abs(Muon_dxy) < 0.05)" # ID and Iso will be asked later depending on W and Z
VETOelectrons = "(Electron_pt>10 and abs(Electron_eta) < 2.5 and Electron_pfRelIso03_all < 0.25 and Electron_mvaFall17V2noIso_WPL)"
LOOSEmuons = "(Muon_pt>10 and abs(Muon_eta)<2.4 and Muon_isGlobal and Muon_isTracker and Muon_pfRelIso04_all < 0.25 and abs(Muon_dz)<0.10 and abs(Muon_dxy) < 0.05 and Muon_looseId)"

VETOLEP = "(Muon_pt>10 and abs(Muon_eta)<2.4 and Muon_pfRelIso04_all < 0.25 and Muon_looseId) or (Electron_pt>10 and abs(Electron_eta) < 2.5 and Electron_pfRelIso03_all < 0.25 and Electron_mvaFall17V2noIso_WPL)"

TRIGGERsMUO = "HLT_IsoMu24"

def selectionTAG(df):

    if isZtag:
        dftag = (df.Define("goodMuons","{}".format(GOODMUON)+" and Muon_mediumId and Muon_pfRelIso04_all < 0.25")
                 .Define("vetoMuons","{}".format(LOOSEmuons))
                 .Define("vetoElectrons","{}".format(VETOelectrons))
                 .Filter("Sum(goodMuons) >= 2 and Sum(Muon_charge[goodMuons])==0 and Sum(vetoElectrons)==0 and Sum(vetoMuons)==2", "At least two good OS muon")
                 .Define("V_mass", "Minv(Muon_pt[goodMuons], Muon_eta[goodMuons], Muon_phi[goodMuons], Muon_mass[goodMuons])")
                 .Filter("V_mass>(91-10) and V_mass<(91+15)","At least one good Z")
                 .Define("trigger","{}".format(TRIGGERsMUO))
                 .Define("Mu1_hasTriggerMatch", "hasTriggerMatch(Muon_eta[goodMuons][0], Muon_phi[goodMuons][0], TrigObj_eta, TrigObj_phi)")
                 .Define("Mu2_hasTriggerMatch", "hasTriggerMatch(Muon_eta[goodMuons][1], Muon_phi[goodMuons][1], TrigObj_eta, TrigObj_phi)")
#                 .Filter("trigger>0 and ((Mu1_hasTriggerMatch and Muon_pt[goodMuons][0]>26) or (Mu2_hasTriggerMatch and Muon_pt[goodMuons][1]>26))","pass trigger")
             )
        return dftag

    if isWtag:
        dftag = (df.Define("goodMuons","{}".format(GOODMUON)+" and Muon_tightId and Muon_pfRelIso04_all < 0.15")
                 .Define("vetoElectrons","{}".format(VETOelectrons))
                 .Define("vetoMuons","{}".format(LOOSEmuons))
                 .Filter("Sum(goodMuons) == 1 and MET_pt>20 and Sum(vetoElectrons)==0 and Sum(vetoMuons)==1", "Exactly 1 good muon and MET>20")
                 .Define("V_mass","mt(Muon_pt[goodMuons][0], Muon_phi[goodMuons][0], MET_pt, MET_phi)")
                 .Filter("V_mass>20","MT>20")
                 .Define("trigger","{}".format(TRIGGERsMUO))
                 .Define("Mu1_hasTriggerMatch", "hasTriggerMatch(Muon_eta[goodMuons][0], Muon_phi[goodMuons][0], TrigObj_eta, TrigObj_phi)")
#                 .Filter("trigger>0 and Mu1_hasTriggerMatch and Muon_pt[goodMuons][0]>26","pass trigger")
             )
        return dftag



def analysis(df,mc,w):

    dftag = selectionTAG(df)

    dfbase = (dftag.Define("goodPhotons", "{}".format(BARRELphotons)+" or {}".format(ENDCAPphotons) )              
              .Filter("Sum(goodPhotons)>0", "At least one good Photon")
              .Define("goodPhotons_pt", "Photon_pt[goodPhotons]")
              .Define("goodPhotons_eta", "Photon_eta[goodPhotons]")
              .Define("goodPhotons_phi", "Photon_phi[goodPhotons]")

#              .Filter("nGenIsolatedPhoton>0", "At least one good Photon")  
#              .Define("HCandMass","MesonCandFromRECO(phi_kin_pt,phi_kin_eta,phi_kin_phi,phi_kin_mass,GenIsolatedPhoton_pt,GenIsolatedPhoton_eta,GenIsolatedPhoton_phi)")
              .Define("goodPhi","{}".format(GOODPHI))
              .Filter("Sum(goodPhi)>0", "At least one Phi with pt > 10 GeV")
              .Define("goodPhi_pt", "phi_kin_pt[goodPhi]")
              .Define("goodPhi_eta", "phi_kin_eta[goodPhi]")
              .Define("goodPhi_phi", "phi_kin_phi[goodPhi]")
              .Define("goodPhi_mass", "phi_kin_mass[goodPhi]")
              .Define("goodPhi_trk1_pt", "phi_trk1_pt[goodPhi]")
              .Define("goodPhi_trk2_pt", "phi_trk2_pt[goodPhi]")
              .Define("goodPhiDR","DeltaR(phi_trk1_eta[goodPhi],phi_trk2_eta[goodPhi],phi_trk1_phi[goodPhi],phi_trk2_phi[goodPhi])")
              .Define("HCandMass","MesonCandFromRECO(goodPhi_pt,goodPhi_eta,goodPhi_phi,goodPhi_mass,goodPhotons_pt,goodPhotons_eta,goodPhotons_phi)")
#              .Filter("abs(HCandMass-125)<30","At least one good Higgs candidate")
              .Define("w","{}".format(w))
              .Define("mc","{}".format(mc))
          )
    
    branchList = ROOT.vector('string')()
    for branchName in [
            "V_mass",
            "HCandMass",
            "goodPhotons_pt",
            "goodPhotons_eta",
#            "index_pair",
#            "nPhoton",
#            "Photon_pt",
#            "Photon_eta",
#            "nGenIsolatedPhoton",
#            "GenIsolatedPhoton_pt",
#            "nphi",
            "goodPhiDR",
            "goodPhi_mass",
            "goodPhi_pt",
            "goodPhi_trk1_pt",
            "goodPhi_trk2_pt",
#            "phi_kin_massErr",
#            "phi_kin_mass",
#            "phi_kin_pt",
#            "phi_kin_vtx_chi2dof",
#            "phi_gen_mass",
#            "phi_gen_pdgId",
            "w",
            "mc",
    ]:
        branchList.push_back(branchName)
        
    outputFile = "outname_mc%d"%mc+".root"
    if isZtag : outputFile = "outname_mc%d"%mc+"_Zcat.root"
    if isWtag : outputFile = "outname_mc%d"%mc+"_Wcat.root"
    print(outputFile)
    snapshot_tdf = dfbase.Snapshot("events", outputFile, branchList)

    print("snapshot_tdf DONE")

    print("---------------- SUMMARY -------------")
    report = dfbase.Report()
    report.Print()

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
    analysis(df,sampleNOW,w)

def readDataSample(sampleNOW):

    files = getDATAlist()
    print(len(files))

    df = ROOT.RDataFrame("Events", files)

    w=1.
    nevents = df.Count().GetValue()
    print("%s entries in the dataset" %nevents)

    analysis(df,sampleNOW,w)

def runTest():

    df = ROOT.RDataFrame("Events", "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/HiggsExo/dalfonso/Hrare/vbf-hphigamma-powheg/NANOAOD_00/step7_VBS_Phigamma_480.root")

    w=1.
    nevents = df.Count().GetValue()
    print("%s entries in the dataset" %nevents)

    sampleNOW=-1
    analysis(df,sampleNOW,w)


   
if __name__ == "__main__":

#    runTest()

    readDataSample(100)  # SingleMuon

    readMCSample(12) # signal VBF
    readMCSample(10) # signal Z
    readMCSample(11) # signal W
    readMCSample(1)  # Zgamma
    readMCSample(0)  # DY
    readMCSample(2)  # Wgamma
    readMCSample(3)  # W
    readMCSample(4)  # ttbar 2L
    readMCSample(5)  # ttbar 1L
    readMCSample(6)  # gJets
    readMCSample(7)  # gJets
    readMCSample(8)  # gJets
    readMCSample(9)  # gJets


