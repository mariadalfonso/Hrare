import ROOT
import os
import math

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

redDark = (191, 34, 41)
redLight = (255,82,82)
orange = (255, 204, 153)
orangeDark = (255,150,79)
azure = (100, 192, 232)
azureDark = (96, 147, 172)
green = (144, 238, 144)
greenDark = (98, 172, 141)

colors = [azure,azureDark,redDark,redLight,orange,orangeDark,green, greenDark]
years = ['12016','22016','2017','2018','12022','22022','12023','22023']
#years = ['12016','22016']
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

labelsJpsiCC={
    '12016': '2016 preVFP', #APV #(B-F for 2016 pre)
    '22016': '2016', #postVFP
    '2017': '2017',
    '2018': '2018',
    '12022': '2022', # C-D
    '22022': '2022 EE', # E, F, G
    '12023': '2023', #C
    '22023': '2023 BPix', #D
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


def getHistoYear(histoName, year,color):
    
    data = []
    if year=="2018": data = [-1,-2,-3,-4]
    if year=="2017": data = [-2,-3,-4,-5,-6]
    if year=="12016": data = [-1,-2,-3,-4,-5,-6]
    if year=="22016": data = [-6,-7,-8]
    if year=="12022": data = [-13,-14]
    if year=="22022": data = [-15,-16,-17]
    if year=="12023": data = [-21,-22,-23,-24]    
    if year=="22023": data = [-31,-32]    
    
    first = False
    for sample in data:
        h = getHisto(sample,histoName,year)
        if not first: hData = h.Clone(); first = True
        elif hData and h: hData.Add(h)

    if hData: hData.SetMarkerStyle(21)
    if hData: hData.SetMarkerSize(1.2)
    if hData: hData.SetMarkerColor(ROOT.TColor.GetColor(*color))
    if hData: hData.SetLineColor(ROOT.TColor.GetColor(*color))    
#    if hData: hData.SetLineWidth(2)

    return hData


def plots(histoName):

    canv = ROOT.TCanvas("stackcanvas","Stack canvas",800,800)
    canv.cd()
    legend = ROOT.TLegend(0.5, 0.6, 0.9, 0.8)
    legend.SetTextFont(42)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)
    legend.SetTextAlign(32)

    first = False
    for year,color in zip(years,colors):
        
        h = getHistoYear(histoName, year, color)
        print(year, " " , color)
        if not first:
            print('hello 1')            
            if('jetFar_Eta' in histoName): h.SetMaximum(2.*h.GetMaximum())
            else: h.SetMaximum(1.5*h.GetMaximum())
            h.DrawNormalized("p e ")
            first = True
            print(first)
            legend.AddEntry(h, str(labelsJpsiCC[year])+"  ("+str(lumisJpsiCC[year])+"  fb^{-1})", "lep")        
            legend.Draw("same")
            canv.Update()
        else:
            print('hello')
            h.DrawNormalized("p e same")
            legend.AddEntry(h, str(labelsJpsiCC[year])+"  ("+str(lumisJpsiCC[year])+"  fb^{-1})", "lep")        
            legend.Draw("same")
            canv.Update()            

    canv.Draw()
    canv.SaveAs("~/public_html/Hrare_JpsiSept25/CompareData_"+histoName+".png")

if __name__ == "__main__":
    
    plots('Jpsi_kin_mass_')
    plots('Jpsi_kin_pt_')
    plots('Jpsi_kin_sipPV_')
    plots('Jpsi_leadingCharged3dSignMaxPt1_')


    plots('jetFar_Eta_')
    plots('jetFar_Pt_')    
    plots('jetFar_CvL_')

    plots('massHiggsCorr_')
