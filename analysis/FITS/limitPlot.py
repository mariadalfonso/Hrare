import numpy as np
#import pandas as pd
import ROOT
ROOT.gROOT.SetBatch(1)
from ROOT import TCanvas, TH1F, TF1, TLegend, gPad, THStack, TColor
import tdrstyle
import os.path
import os

from ROOT import TCanvas, TGraph
from ROOT import gROOT
from math import sin
from array import array


doChekSync=True
#mesonCat='Phi'
mesonCat='Rho'

dirTO_ = '/work/submit/mariadlf/cards_JUNE22/'
dirMIT_ = '/work/submit/mariadlf/cards_JUNE26/'

if doChekSync:
     ax = range(1,5) #
     labels = [ 'ggH(TO) noMVA'
                , 'ggH(MIT) noMVA'
                , 'ggH(TO) w/ MVA'
                , 'ggH(MIT) w/ MVA'
#                , 'ggH(TO) w/ MVA Phi'
#                , 'ggH(MIT) w/ MVA Phi'
               ]

     if mesonCat=='Rho': limitfiles = [ dirTO_+"230622_Rho/higgsCombineRhoGFpreselection.AsymptoticLimits.mH125.root"
                                        ,dirMIT_+"workspaces_JUNE22_40_40/higgsCombineRhoGFcat.AsymptoticLimits.mH125.root"
                                        ,dirTO_+"230622_Rho/higgsCombineRhoGFcat_Comb.AsymptoticLimits.mH125.root"
                                        ,dirMIT_+"workspaces_MVA_JUNE22_40_40/higgsCombineRhoGFcat.AsymptoticLimits.mH125.root"
                                        ##                    ,"limitsGiulio/higgsCombinePhiGFcat.AsymptoticLimits.mH125.root"
                                        ##                    ,"limitsMVASynch/higgsCombinePhiGFcat.AsymptoticLimits.mH125.root"
                                       ]

     if mesonCat=='Phi': limitfiles = [ dirTO_+"230622_Phi/higgsCombinePhiGFpreselection.AsymptoticLimits.mH125.root"
                                        ,dirMIT_+"workspaces_JUNE22_40_40/higgsCombinePhiGFcat.AsymptoticLimits.mH125.root"
                                        ,dirTO_+"230622_Phi/higgsCombinePhiGFcat_Comb.AsymptoticLimits.mH125.root"
                                        ,dirMIT_+"workspaces_MVA_JUNE22_40_40/higgsCombinePhiGFcat.AsymptoticLimits.mH125.root"
                                        ##                    ,"limitsGiulio/higgsCombinePhiGFcat.AsymptoticLimits.mH125.root"
                                        ##                    ,"limitsMVASynch/higgsCombinePhiGFcat.AsymptoticLimits.mH125.root"
                                       ]


else:
     ax = range(1,5) #
#     ax = range(1,6) #
     labels = [ 'ggH'
                 , 'VBF'
                 , 'VBFlow'
#                 , 'Zinv'
                 , '1l-2l'
               ]
     limitfiles = [ "DATACARDSmva_MARCH11/higgsCombine"+mesonCat+"GFcat.AsymptoticLimits.mH125.root"
                    ,"DATACARDSmva_MARCH11/higgsCombine"+mesonCat+"VBFcat.AsymptoticLimits.mH125.root"
                    ,"DATACARDSmva_MARCH11/higgsCombine"+mesonCat+"VBFcatlow.AsymptoticLimits.mH125.root"
#                    ,"limitsMVADEC28/higgsCombine"+mesonCat+"Zinvcat.AsymptoticLimits.mH125.root"
                    ,"DATACARDS_MARCH11/higgsCombine"+mesonCat+"Vcat.AsymptoticLimits.mH125.root"
                   ]

c44 = ROOT.TCanvas("c44","c44",1200,800)
aexl = [0.5 for a in ax]
aexh = [0.5 for a in ax]
ay      = [  ]
aeyl    = [  ]
aeyh    = [  ]
aeyl_y  = [  ]
aeyh_y  = [  ]
obs = [ ]

for f in limitfiles:
     fo= ROOT.TFile.Open(f)
     tr = fo.Get("limit")
     tr.GetEntry(0)
     aeyl_y.append(tr.limit)
     tr.GetEntry(1)
     aeyl.append(tr.limit)
     tr.GetEntry(2)
     ay.append(tr.limit)
     tr.GetEntry(3)
     aeyh.append(tr.limit)
     tr.GetEntry(4)
     aeyh_y.append(tr.limit)
     tr.GetEntry(5)
     obs.append(tr.limit)
    
n = len(ax)

col = ROOT.TColor()
green = col.GetColor('#19C405')
yellow = col.GetColor('#FEC30A')

legend_args = (0.64, 0.7, 0.92, 0.85, '', 'NDC')

legend = TLegend(*legend_args)
legend.SetFillStyle(0)
legend.SetTextFont(42)

x, y, exl, exh, eyl, eyh, eyl_y, eyh_y, eyl_n, eyh_n = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ) , array( 'd' ), array( 'd' )
#x, y, exl, exh, eyl, eyh = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )
print(ay[0])

yo = array( 'd' )
yow = array( 'd' )
for i in range( n ): 
     
    print('i=',i,' ;limit=',ay[i])
    x.append( ax[i])
    y.append( ay[i])
    
    exh.append( aexh[i] )
    exl.append( aexl[i] )
    
    eyh.append( (aeyh[i]-ay[i]) )
    eyl.append( (ay[i]-aeyl[i]) )
    
    eyh_y.append( (aeyh_y[i]-ay[i]) )
    eyl_y.append( (ay[i]-aeyl_y[i]) )
    
    eyh_n.append( 0 )
    eyl_n.append( 0 )
    
    yow.append(0.5)
    yo.append(obs[i])
    
gae = ROOT.TGraphAsymmErrors(n, x, y, exl, exh, eyl, eyh)
yae = ROOT.TGraphAsymmErrors(n, x, y, exl, exh, eyl_y, eyh_y)
nae = ROOT.TGraphAsymmErrors(n, x, y, exl, exh, eyl_n, eyh_n)
gro = ROOT.TGraphErrors(n, x, yo,yow,eyh_n)
gro.SetMarkerStyle(20)
gro.SetMarkerColor(1)
gro.SetLineColor(1)
gro.SetLineWidth(2)
gro.SetMarkerSize(1.5)

gae.SetFillColor(green)
#gae.SetFillStyle(3001)
gae.SetMarkerStyle(0)

nae.SetLineStyle(2)
nae.SetLineWidth(2)
nae.SetMarkerSize(0)

yae.SetFillColor(yellow)
#yae.SetFillStyle(3001)
yae.SetMarkerStyle(0)


#legend.AddEntry(gro, "Observed", 'lp')
legend.AddEntry(nae, "Median expected", 'lp')
legend.AddEntry(gae, "68% expected", 'f')
legend.AddEntry(yae, "95% expected", 'f')

dummyHistogram = ROOT.TH1F("dummy","",n,0.5,n+.5)
dummyHistogram.GetXaxis().SetTitleOffset(1.8)
dummyHistogram.GetXaxis().SetLabelSize(0.045)
dummyHistogram.GetXaxis().SetRangeUser(0.5, ax[-1]+0.5)
dummyHistogram.GetYaxis().SetRangeUser(0, 5.) # all
if doChekSync: dummyHistogram.GetYaxis().SetRangeUser(0, 0.5) # GF and VBF
dummyHistogram.GetYaxis().SetTitle("95% CL upper limit on #sigma #times B(H#rightarrow "+mesonCat+" #gamma) x 10^-2")
if doChekSync: dummyHistogram.GetYaxis().SetTitle("95% CL upper limit on #sigma #times B(H#rightarrow Meson #gamma) x 10^-2")
dummyHistogram.GetYaxis().SetTitleSize(0.04)
dummyHistogram.GetYaxis().SetTitleOffset(1.8)
if doChekSync and mesonCat=='Rho' : dummyHistogram.SetTitle("Rho")
if doChekSync and mesonCat=='Phi' : dummyHistogram.SetTitle("Phi")


xax = dummyHistogram.GetXaxis()
for i in range(1, n+1):
#for i in range(0, n):     
    binIndex = xax.FindBin(i)
    print(binIndex, " ", str(labels[i-1]))
    xax.SetBinLabel(binIndex, str(labels[i-1]))
#xax.LabelsOption("v")
xax.SetTickLength(0)
dummyHistogram.Draw("AXIS")

#.LabelsOption("v", "X")
gae.Draw("p e2")
yae.Draw("e2 same")
gae.Draw("e2 same")
nae.Draw("p same")
#gro.Draw("p same")  # comment observed

legend.Draw()

gPad.RedrawAxis()
if doChekSync: tdrstyle.cmsPrel(39540, energy= 13,simOnly=False)
else: tdrstyle.cmsPrel(138000, energy= 13,simOnly=False)
c44.Draw()

latex = ROOT.TLatex()
latex.SetTextSize(0.04)
if mesonCat=='Rho': latex.DrawLatex(2, 0.85*dummyHistogram.GetMaximum(), "Rho")
if mesonCat=='Phi': latex.DrawLatex(2, 0.85*dummyHistogram.GetMaximum(), "Phi")

mySTR=''
if doChekSync: mySTR="_synchJune22"
else: mySTR='_final'

c44.SaveAs("limit"+mesonCat+mySTR+".pdf")
c44.SaveAs("limit"+mesonCat+mySTR+".png")
