import ROOT
import os

ROOT.ROOT.EnableImplicitMT()
from utilsHrare import findDataset, findDIR
from utilsHrare import plot
from utilsHrare import SwitchSample

doPlot = False
lumi = 137.
isZtag = False;
isWtag = True;



BARRELphotons = "(Photon_pt>20 and abs(Photon_eta)<1.4442 and Photon_sieie < 0.01015 and Photon_hoe < 0.02197 and Photon_electronVeto and (Photon_cutBased & 2) and Photon_pfRelIso03_chg < 0.65)"
ENDCAPphotons = "(Photon_pt>20 and abs(Photon_eta)>1.566 and abs(Photon_eta)<2.5 and Photon_sieie < 0.0272 and Photon_hoe < 0.0326 and Photon_electronVeto and (Photon_cutBased & 2) and Photon_pfRelIso03_chg < 0.65)"

GOODPHI = "(phi_kin_pt>5 && phi_kin_massErr<0.025)"

GOODMUON = "(Muon_pt>20 and abs(Muon_eta)<2.4 and Muon_pfRelIso04_all < 0.15 and Muon_tightId)"
VETOLEP = "(Muon_pt>10 and abs(Muon_eta)<2.4 and Muon_pfRelIso04_all < 0.25 and Muon_looseId) or (Electron_pt>10 and abs(Electron_eta) < 2.5 and Electron_pfRelIso03_all < 0.25 and Electron_mvaFall17V2noIso_WPL)"


def selectionTAG(df):

    if isZtag:
        dftag = (df.Define("goodMuons","{}".format(GOODMUON))
                 .Define("vetoLeptons","{}".format(VETOLEP))
                 .Filter("Sum(goodMuons) >= 2 and Sum(Muon_charge[goodMuons])==0 ", "At least two good OS muon")
                 .Define("V_mass", "Minv(Muon_pt[goodMuons], Muon_eta[goodMuons], Muon_phi[goodMuons], Muon_mass[goodMuons])")
                 .Filter("abs(V_mass-91)<15","At least one good Z")
             )
        return dftag

    if isWtag:
        dftag = (df.Define("goodMuons","{}".format(GOODMUON))
                 .Define("vetoLeptons","{}".format(VETOLEP))
                 .Filter("Sum(goodMuons) == 1 and MET_pt>20", "Exactly 1 good muon and MET>20")
                 .Define("V_mass","mt(Muon_pt[goodMuons][0], Muon_phi[goodMuons][0], MET_pt, MET_phi)")
                 .Filter("V_mass>20","MT>20")
             )
        return dftag



def analysis(df,mc,w):

    dftag = selectionTAG(df)

    dfbase = (dftag.Define("goodPhotons", "{}".format(BARRELphotons)+" or {}".format(ENDCAPphotons) )              
              .Filter("Sum(goodPhotons)>0", "At least one good Photon")
#              .Filter("nGenIsolatedPhoton>0", "At least one good Photon")  
#              .Define("HCandMass","MesonCandFromRECO(phi_kin_pt,phi_kin_eta,phi_kin_phi,phi_kin_mass,GenIsolatedPhoton_pt,GenIsolatedPhoton_eta,GenIsolatedPhoton_phi)")
              .Define("goodPhi","{}".format(GOODPHI))
              .Filter("Sum(goodPhi)>0", "At least one Phi with pt > 5 GeV and massErr < 0.025")
              .Define("HCandMass","MesonCandFromRECO(phi_kin_pt[goodPhi],phi_kin_eta[goodPhi],phi_kin_phi[goodPhi],phi_kin_mass[goodPhi],Photon_pt[goodPhotons],Photon_eta[goodPhotons],Photon_phi[goodPhotons])")
#              .Filter("abs(HCandMass-125)<30","At least one good Higgs candidate")
              .Define("w","{}".format(w))
              .Define("mc","{}".format(mc))
          )
    
    branchList = ROOT.vector('string')()
    for branchName in [
            "V_mass",
            "HCandMass",
            "nPhoton",
            "Photon_pt",
            "Photon_eta",
            "nGenIsolatedPhoton",
            "GenIsolatedPhoton_pt",
            "nphi",
            "phi_kin_massErr",
            "phi_kin_mass",
            "phi_kin_pt",
            "phi_kin_vtx_chi2dof",
            "phi_gen_mass",
            "phi_gen_pdgId",
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

def readSample(sampleNOW):

    files = findDIR("/work/submit/mariadlf/Hrare/OCT14/{}".format(SwitchSample(sampleNOW)[0]))

    print(len(files)) 
    df = ROOT.RDataFrame("Events", files)   

    nevents = df.Count().GetValue()  ## later with negative weights
    w = (SwitchSample(sampleNOW)[1] / nevents)

    lumiEq = (nevents / SwitchSample(sampleNOW)[1])
    print("%s entries in the dataset" %nevents)
    print("lumi equivalent fb %s" %lumiEq)
    analysis(df,sampleNOW,w)

   
if __name__ == "__main__":

    readSample(10) # signal Z
    readSample(11) # signal W
    readSample(1)  # Zgamma
    readSample(2)  # Wgamma
    readSample(3)  # W
    readSample(0)  # DY
    readSample(4)  # ttbar 2L



