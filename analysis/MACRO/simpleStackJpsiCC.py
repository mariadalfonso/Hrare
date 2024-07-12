import ROOT
import os
import math

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

year = '2018'
redDark = (191, 34, 41)
redLight = (255,82,82)
orange = (255, 204, 153)
azure = (100, 192, 232)
azureDark = (96, 147, 172)

dirname='/work/submit/mariadlf/Hrare_JPsiCC/JUL4/'

def getHisto(sample,histoName):

    fileName=dirname+"histoOUTname_"+str(sample)+"_2018.root"
    inFile = ROOT.TFile.Open(fileName,"READ")
#    inFile.ls()
    histo=histoName+year
    h = inFile.Get(histo)
    h.SetDirectory(0)
    return h

def plots(histoName):

    h1 = getHisto(10,histoName)
    h2 = getHisto(11,histoName)
    h2.Scale(0.75) # this is an ad-Hoc -- B
    
    print('int=',h1.Integral(),' Charm name=',h1.GetName())
    print('int=',h2.Integral(),' B name=',h2.GetName())    
        
    h3 = getHisto(-1,histoName)
    h4 = getHisto(-2,histoName)
    h5 = getHisto(-3,histoName)
    h6 = getHisto(-4,histoName)
    
    for h, color in zip([h1,h2],[azure,orange]):
        if h:
            h.SetLineWidth(3)
            h.SetLineColor(ROOT.TColor.GetColor(*color))
            h.SetFillColor(ROOT.TColor.GetColor(*color))
            
    stack = ROOT.THStack("stack","")
    for h in [h2,h1]:
        stack.Add(h)

    hData = h3.Clone()
    if hData: hData.Add(h4)
    if hData: hData.Add(h5)
    if hData: hData.Add(h6)
    if hData: hData.SetMarkerStyle(21)
    if hData: hData.SetMarkerSize(1.2)
    if hData: hData.SetMarkerColor(ROOT.kBlack)
    if hData: hData.SetLineColor(ROOT.kBlack)
    if hData: hData.SetLineWidth(2)

    hSig = getHisto(1000,histoName)
    hSig.SetLineColor(ROOT.kRed)
    hSig.Scale(hData.Integral()/hSig.Integral())
    
    canv = ROOT.TCanvas("stackcanvas","Stack canvas",800,800)
    if histoName=='Jpsi_kin_massErr_' or histoName=='minDRjpsi_': canv.SetLogy(1)
    hData.Draw("p e")
    stack.Draw("HIST same")
    stack.GetYaxis().SetTitle("Events")
    hData.Draw("p e same")
    hSig.Draw("hist same")
    
    legend = ROOT.TLegend(0.43, 0.6, 0.90, 0.88)
    legend.SetTextFont(42)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)
    legend.SetTextAlign(32)
    legend.AddEntry(hData, "Data-2018" ,"lep")
    legend.AddEntry(h1, "Charmonium" , "f")
    legend.AddEntry(h2, "BToJpsi (x 0.75 ad hoc)" , "f")
    legend.Draw()
    
    canv.Draw()
    canv.SaveAs("~/public_html/Hrare_JpsiJul4/Stack_"+histoName+".png")

if __name__ == "__main__":

    plots('Jpsi_iso_')
    
    plots('Jpsi_kin_pt_')
    plots('Jpsi_kin_eta_')
    plots('Jpsi_kin_valid_')    
    plots('Jpsi_kin_mass_')
    plots('Jpsi_kin_massErr_')
    plots('Jpsi_kin_sipPV_')
    #
    plots('Jpsi_iso_')        
    plots('Jpsi_leadingCharged3dSignMaxPt1_')    
    plots('Jpsi_leadingChargedMaxPt1_')
    plots('Jpsi_muon1_pt_')
    plots('Jpsi_muon2_pt_')    
    #
    plots('jetClosePt_')
    plots('jetFarPt_')        
    plots('jetCloseCvL_')
    plots('jetFarCvL_')
    plots('minDRjpsi_')
    plots('jetCloseJPsiRatio_')
    #
    plots('massHiggsCorr_')    

    
