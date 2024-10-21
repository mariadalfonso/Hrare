import ROOT
import os
import math

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

years = ['12016','22016','2017','2018','12022','12023']
redDark = (191, 34, 41)
redLight = (255,82,82)
orange = (255, 204, 153)
azure = (100, 192, 232)
azureDark = (96, 147, 172)
green = (144, 238, 144)
greenDark = (98, 172, 141)

dirname='/work/submit/mariadlf/Hrare_JPsiCC/Sept25/'

lumisJpsiCC={
    '12016': 19.52, #APV #(B-F for 2016 pre)
    '22016': 16.80, #postVFP
    '2017': 41.5,
    '2018': 59.70,
    '12022':7.98, # C-D
    '22022':26.67, # E, F, G
    '12023':17.794, #C
    '22023':9.451, #D
}


def getHisto(sample,histoName,year):

    fileName=dirname+"histoOUTname_"+str(sample)+"_"+str(year)+".root"
    inFile = ROOT.TFile.Open(fileName,"READ")
#    inFile.ls()
    histo=histoName+year
    h = inFile.Get(histo)
    print(h)
    h.SetDirectory(0)
    return h

def plots(histoName,year):

    ### MC
    
    if year=="12022" or year=="12023":
        h1 = getHisto(12,histoName,"12022")
        stack = ROOT.THStack("stack","")
        h1.SetLineWidth(3)
        h1.SetLineColor(ROOT.TColor.GetColor(*green))
#        h1.SetFillColor(ROOT.TColor.GetColor(*green))
        for h in [h1]:
            print('int=',h.Integral(),' name=',h.GetName())
            stack.Add(h)
    else:
        h1 = getHisto(10,histoName,year)
        h2 = getHisto(11,histoName,year)
        for h, color in zip([h1,h2],[azure,orange]):
            if h:
                h.SetLineWidth(3)
                h.SetLineColor(ROOT.TColor.GetColor(*color))
                h.SetFillColor(ROOT.TColor.GetColor(*color))

        stack = ROOT.THStack("stack","")
        for h in [h2,h1]:
            print('HELLO int=',h.Integral(),' name=',h.GetName())
            stack.Add(h)

    ### DATA

    data = []
    if year=="2018": data = [-1,-2,-3,-4]
    if year=="2017": data = [-2,-3,-4,-5,-6]
    if year=="12016": data = [-1,-2,-3,-4,-5,-6]
    if year=="22016": data = [-6,-7,-8]
    if year=="12022": data = [-13,-14]
    if year=="12023": data = [-21,-22,-23,-24]
    
    first = 'false'
    for sample in data:
        h = getHisto(sample,histoName,year)
        if first == 'false': hData = h.Clone(); first = 'true'
        if hData and h: hData.Add(h)

    if hData: hData.SetMarkerStyle(21)
    if hData: hData.SetMarkerSize(1.2)
    if hData: hData.SetMarkerColor(ROOT.kBlack)
    if hData: hData.SetLineColor(ROOT.kBlack)
    if hData: hData.SetLineWidth(2)

    ### SIGNAL

    hSig = getHisto(1000,histoName,'2018')
    hSig.SetLineColor(ROOT.kRed)
    hSig.Scale(hData.Integral()/hSig.Integral())

    hSigZ = getHisto(1001,histoName,'2018')
    hSigZ.SetLineColor(ROOT.kBlue)
    hSigZ.Scale(hData.Integral()/hSigZ.Integral())

    ### MAKE PLOT
    
    canv = ROOT.TCanvas("stackcanvas","Stack canvas",800,800)
    if histoName=='Jpsi_kin_massErr_' or histoName=='minDRjpsi_': canv.SetLogy(1)
    hData.SetMaximum(1.5*max(hData.GetMaximum(),hSig.GetMaximum()))
    hData.SetTitle("")
    hData.Draw("p e")
    stack.Draw("HIST same")
    stack.GetYaxis().SetTitle("Events")
    hData.Draw("p e same")
    hSig.Draw("hist same")
    hSigZ.Draw("hist same")

    legend = ROOT.TLegend(0.43, 0.6, 0.90, 0.88)
    legend.SetTextFont(42)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)
    legend.SetTextAlign(32)
    legend.AddEntry(hData, "Data"+year ,"lep")
    if year == '12022' or year == '12023':
        legend.AddEntry(h1, "Inclusive Mu" , "f")
    else:
        legend.AddEntry(h1, "Charmonium" , "f")
        legend.AddEntry(h2, "BToJpsi" , "f")

    legend.Draw()

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(72)
    text.SetTextSize(0.04)
    if year == '12022' or year == '12023': text.DrawLatex(0.5, 0.92, "#sqrt{s} = 13.6 TeV,%0.2f fb^{-1}"% (lumisJpsiCC[year]))
    else: text.DrawLatex(0.55, 0.92, "#sqrt{s} = 13 TeV,%0.2f fb^{-1}"% (lumisJpsiCC[year]))
    
    canv.Draw()
    canv.SaveAs("~/public_html/Hrare_JpsiSept25/Stack_"+histoName+"_"+year+".png")

if __name__ == "__main__":

    for year in years:
    
        plots('Jpsi_iso_',year)

        plots('Jpsi_kin_pt_',year)
        plots('Jpsi_kin_eta_',year)
        plots('Jpsi_kin_valid_',year)
        plots('Jpsi_kin_mass_',year)
        plots('Jpsi_kin_massErr_',year)
        plots('Jpsi_kin_sipPV_',year)
        #
        plots('Jpsi_iso_',year)
        plots('Jpsi_leadingCharged3dSignMaxPt1_',year)
        plots('Jpsi_leadingChargedMaxPt1_',year)
        plots('Jpsi_muon1_pt_',year)
        plots('Jpsi_muon2_pt_',year)
        #
        plots('jetClose_Pt_',year)
        plots('jetFar_Pt_',year)
        plots('jetClose_Eta_',year)
        plots('jetFar_Eta_',year)
        plots('jetClose_CvL_',year)
        plots('jetFar_CvL_',year)
        plots('minDRjpsi_',year)
        plots('jetCloseJPsiRatio_',year)
        #
        plots('massHiggsCorr_',year)
        #
        plots('angle_phi1_',year)
        plots('angle_phi2_',year)
        plots('angle_phi3_',year)
        plots('angle_phi4_',year)
        plots('angle_cos_delta_m_min_',year)
        plots('angle_cos_theta_jpsi_',year)
        plots('angle_cos_theta_dicharm_',year)
        #
        plots('discrMVA_',year)
