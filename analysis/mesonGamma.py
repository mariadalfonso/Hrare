import ROOT
import os

ROOT.ROOT.EnableImplicitMT()

if "/functions.so" not in ROOT.gSystem.GetLibraries():
    ROOT.gSystem.CompileMacro("functions.cc","k")

def plot(h,filename,doLogX):

   ROOT.gStyle.SetOptStat(1); 
   ROOT.gStyle.SetTextFont(42)
   c = ROOT.TCanvas("c", "", 800, 700)
   if doLogX: c.SetLogx(); 
#  c.SetLogy()

   h.SetTitle("")
   h.GetXaxis().SetTitleSize(0.04)
   h.GetYaxis().SetTitleSize(0.04)
   h.SetLineColor(4)
#   h.SetLineColor(2)
#   h.SetLineColor(46)
#   h.SetLineColor(30) # green-ish
   h.Draw()

   label = ROOT.TLatex(); label.SetNDC(True)
#   label.DrawLatex(0.175, 0.740, "#eta")
#   label.DrawLatex(0.205, 0.775, "#rho,#omega")
#   label.DrawLatex(0.270, 0.740, "#phi")
#   label.SetTextSize(0.040); label.DrawLatex(0.100, 0.920, "#bf{CMS Open Data}")
#   label.SetTextSize(0.030); label.DrawLatex(0.630, 0.920, "#sqrt{s} = 8 TeV, L_{int} = 11.6 fb^{-1}")

   print("saving file: {} ".format(filename))

   c.SaveAs(filename)
 

PHIid=333
RHOid=113
PIid=211
Kid=321

def selection_1GenPhi(df,pdgId):

    df_1M = df.Define("goodMeson","abs(GenPart_pdgId)=={}".format(pdgId));
    return df_1M;


def selection_2Kaons(df,pdgId):

    df_2k = df.Define("goodkaons", "abs(GenPart_pdgId)=={} and abs(GenPart_eta)<2.5".format(pdgId)).Filter("Sum(goodkaons) >= 2 ")
    return df_2k;


def getPhiCand(df, pdgTest, signalTest):

    # 1 genPhi + 1Photon
    df_genPhiCand = selection_1GenPhi(df.Filter("nGenIsolatedPhoton>=1", "At least one GenIsolatedPhoton"),pdgTest)

    # 2 recoLike + 1Photon
    if(pdgTest==333): df_2k = selection_2Kaons(df,321)
    if(pdgTest==113): df_2k = selection_2Kaons(df,2111)
    df_2k_1g = df_2k.Filter("nGenIsolatedPhoton>=1", "At least one GenIsolatedPhoton")
    
    if(pdgTest==333): phiHyp="true";
    if(pdgTest==113): phiHyp="false";
    mesonCand_ = df_2k_1g.Define("index_pair","mesonCand(GenPart_pt[goodkaons], GenPart_eta[goodkaons], GenPart_phi[goodkaons], GenPart_mass[goodkaons], GenPart_pdgId[goodkaons],GenIsolatedPhoton_pt,GenIsolatedPhoton_eta,GenIsolatedPhoton_phi,{})".format(phiHyp)).Filter("index_pair[0]!= -1 && index_pair[1]!= -1")

    ##$$$$
    # make variables
    ##$$$$

    df_HCand_Mass = mesonCand_.Define("H_mass", "Minv3(GenPart_pt[goodkaons], GenPart_eta[goodkaons], GenPart_phi[goodkaons], GenPart_mass[goodkaons], index_pair, GenIsolatedPhoton_pt,GenIsolatedPhoton_eta,GenIsolatedPhoton_phi)")

    df_phiCand_Mass = mesonCand_.Define("kk_mass", "Minv(GenPart_pt[goodkaons], GenPart_eta[goodkaons], GenPart_phi[goodkaons], GenPart_mass[goodkaons], index_pair)")
    df_phiCand_Pt = mesonCand_.Define("kk_pt", "PT(GenPart_pt[goodkaons], GenPart_eta[goodkaons], GenPart_phi[goodkaons], GenPart_mass[goodkaons], index_pair)")
    df_dPhi_MvsPh = mesonCand_.Define("dPhi_MvsPh", "dPhi_MvsPh(GenPart_pt[goodkaons], GenPart_eta[goodkaons], GenPart_phi[goodkaons], GenPart_mass[goodkaons], index_pair,GenIsolatedPhoton_pt,GenIsolatedPhoton_eta,GenIsolatedPhoton_phi)")
    df_dEta_MvsPh = mesonCand_.Define("dEta_MvsPh", "dEta_MvsPh(GenPart_pt[goodkaons], GenPart_eta[goodkaons], GenPart_phi[goodkaons], GenPart_mass[goodkaons], index_pair,GenIsolatedPhoton_eta)")
    df_dR_MConst = mesonCand_.Define("dR_MConst", "dR_Constituents(GenPart_pt[goodkaons], GenPart_eta[goodkaons], GenPart_phi[goodkaons], GenPart_mass[goodkaons], index_pair)")

    hgPt = df.Define("goodPhotons_pt","GenIsolatedPhoton_pt").Filter("nGenIsolatedPhoton>=1", "At least one GenIsolatedPhoton").Histo1D(("Ph_pt", "Ph_pt (GeV); pt_{gamma}; N_{Events}", 1000, 0.25, 300), "goodPhotons_pt")


    ##$$$$
    #  Histograms
    ##$$$$

    name="_BKGlike"

    if pdgTest == 333:
        if signalTest:
            name="_fromPhi"

        hMass = df_phiCand_Mass.Histo1D(("Dikaon_mass", "Dikaon mass;m_{k^{+}k^{-}} (GeV);N_{Events}", 3000, 0.25, 300), "kk_mass")
        hPt = df_phiCand_Pt.Histo1D(("Dikaon_pt", "Dikaon pt;pT_{k^{+}k^{-}} (GeV);N_{Events}", 3000, 0.25, 300), "kk_pt")
        hdPhi = df_dPhi_MvsPh.Histo1D(("dPhi_MvsPh", "Dikaon gamma ;DeltaPhi(k^{+}k^{-},#gamma);N_{Events}", 1000, 0., 3.14), "dPhi_MvsPh")
        hdEta = df_dEta_MvsPh.Histo1D(("dEta_MvsPh", "Dikaon gamma; DeltaEta(k^{+}k^{-},#gamma);N_{Events}", 1000, -10, 10), "dEta_MvsPh")
        hdRConst = df_dR_MConst.Histo1D(("dR_MConst", "DR kaons ; DeltaR(k^{+},k^{-});N_{Events}", 1000, -10, 10), "dR_MConst")
        hHiggsM = df_HCand_Mass.Histo1D(("DikaonGamma_mass", "DikaonGamma mass;m_{k^{+}k^{-}#gamma} (GeV);N_{Events}", 3000, 0.25, 300), "H_mass")
        plot(hdPhi,"~/www/Hrare/plot_kk_deltaPhiVPh"+name+".png",False)        
        plot(hdEta,"~/www/Hrare/plot_kk_deltaEtaVPh"+name+".png",False)        
        plot(hPt,"~/www/Hrare/plot_kkPt"+name+".png",True)
        plot(hdRConst,"~/www/Hrare/plot_kkhdRConst"+name+".png",False)
        plot(hMass,"~/www/Hrare/plot_kkMass"+name+".png",True)
        plot(hHiggsM,"~/www/Hrare/plot_kkHiggsMass"+name+".png",False)
        plot(hgPt,"~/www/Hrare/plot_gammaPt"+name+".png",True)
    
    if pdgTest == 113:
        if signalTest:
            name="_fromRho"
        hMass = df_phiCand_Mass.Histo1D(("Dipion_mass", "Dipion mass;m_{#pi^{+}#pi^{-}};N_{Events}", 3000, 0.25, 300), "kk_mass")
        hPt = df_phiCand_Pt.Histo1D(("Dipion_pt", "Dipion pt;pT_{#pi^{+}#pi^{-}};N_{Events}", 3000, 0.25, 300), "kk_pt")
        hdPhi = df_dPhi_MvsPh.Histo1D(("dPhi_MvsPh", "Dipion mass;DeltaPhi_(#pi^{+}#pi^{-},#gamma);N_{Events}", 1000, 0., 3.14), "dPhi_MvsPh")
        hdEta = df_dEta_MvsPh.Histo1D(("dEta_MvsPh", "Dipion gamma; DeltaEta(#pi^{+}#pi^{-},#gamma);N_{Events}", 1000, -10, 10), "dEta_MvsPh")
        hdRConst = df_dR_MConst.Histo1D(("dR_MConst", "DR kaons ; DeltaR(#pi^{+},#pi^{-});N_{Events}", 1000, -10, 10), "dR_MConst")
        hHiggsM = df_HCand_Mass.Histo1D(("DikaonGamma_mass", "DikaonGamma mass;m_{#pi^{+}#pi^{-}#gamma} (GeV);N_{Events}", 3000, 0.25, 300), "H_mass")
        plot(hdPhi,"~/www/Hrare/plot_pipi_deltaPhiVPh"+name+".png",False)        
        plot(hdEta,"~/www/Hrare/plot_pipi_deltaEtaVPh"+name+".png",False)        
        plot(hPt,"~/www/Hrare/plot_pipiPt"+name+".png",True)
        plot(hdRConst,"~/www/Hrare/plot_pipihdRConst"+name+".png",False)
        plot(hMass,"~/www/Hrare/plot_pipiMass"+name+".png",True)
        plot(hHiggsM,"~/www/Hrare/plot_pipiHiggsMass"+name+".png",False)
        plot(hgPt,"~/www/Hrare/plot_gammaPt"+name+".png",True)
        

    print("%s entries in the dataset" %df.Count().GetValue())
    print("%s entries with 1#gamma + 2Kaons cand" %mesonCand_.Count().GetValue())
    print("%s entries with 1#gamma + 2Kaons -- H cand" %df_HCand_Mass.Count().GetValue())

def mesonGamma(pdgTest,dosignal):
    
    df = ROOT.RDataFrame("Events", "$CMSSW_BASE/src/genproduction/phiGammaToKpKm-NANOGEN.root")    
    
    if not dosignal: 
        df = ROOT.RDataFrame("Events", "$CMSSW_BASE/src/genproduction/rhoGammaToPipPim-NANOGEN.root")

    if pdgTest==113 and dosignal: 
        df = ROOT.RDataFrame("Events", "$CMSSW_BASE/src/genproduction/rhoGammaToPipPim-NANOGEN.root")
    if pdgTest==113 and not dosignal: 
        df = ROOT.RDataFrame("Events", "$CMSSW_BASE/src/genproduction/phiGammaToKpKm-NANOGEN.root")    
        
    getPhiCand(df,pdgTest,dosignal)
        
   
if __name__ == "__main__":
#    run_fast = True
    mesonGamma(333,False)
