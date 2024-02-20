import ROOT
import os
import math

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

year = '2018'
histoName='HCandMass_'
redDark = (191, 34, 41)
redLight = (255,82,82)
orange = (255, 204, 153)

inFileName = "../DASKlogs/histoOUTname_test_interactive.root"
inFile = ROOT.TFile.Open(inFileName ,"READ")
h1 = inFile.Get(histoName+year+"_1017")
h2 = inFile.Get(histoName+year+"_1010")
h3 = inFile.Get(histoName+year+"_10")
h4 = inFile.Get(histoName+year+"_11")
h5 = inFile.Get(histoName+year+"_12")
h6 = inFile.Get(histoName+year+"_13")
h7 = inFile.Get(histoName+year+"_14")

h8 = inFile.Get(histoName+year+"_-62")
h9 = inFile.Get(histoName+year+"_-63")
h10 = inFile.Get(histoName+year+"_-64")

hGjets = h3.Clone()
hGjets.Add(h4)
hGjets.Add(h5)
hGjets.Add(h6)
hGjets.Add(h7)

hData = h8.Clone()
if hData: hData.Add(h9)
if hData: hData.Add(h10)
if hData: hData.SetMarkerStyle(20)
if hData: hData.SetMarkerSize(1.2)
if hData: hData.SetMarkerColor(ROOT.kBlack)
if hData: hData.SetLineColor(ROOT.kBlack)
if hData: hData.SetLineWidth(2)

for h, color in zip([h1,h2,hGjets],[redDark,redLight,orange]):
    if h:
        h.Scale(38)
        h.SetLineWidth(3)
        h.SetLineColor(ROOT.TColor.GetColor(*color))
        h.SetFillColor(ROOT.TColor.GetColor(*color))

stack = ROOT.THStack("stack","")
for h in [hGjets,h2,h1]:
    stack.Add(h)

canv = ROOT.TCanvas("stackcanvas","Stack canvas",800,800)
canv.SetLogy(1)
stack.Draw("HIST")
hData.Draw("p e same")

canv.Draw()
canv.SaveAs("~/public_html/Stack.png")
