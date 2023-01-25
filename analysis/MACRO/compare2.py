import ROOT
import os
from array import array
import math

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
#ROOT.SetBatch()

#mesonCat = '_PhiCat'
mesonCat = '_RhoCat'
year = '_2018'


directory = '/home/submit/mariadlf/Hrare/analysis/DEC28/2018/'
mytreeGF = ROOT.TChain('events')
mytreeGFfix = ROOT.TChain('events')

#RHO
#if mesonCat == "_RhoCat":
#   mytreeGF.Add(directory+'outname_mc1027_GFcat'+mesonCat+year+'.root')
#   mytreeGFfix.Add(directory+'outname_mc1051_GFcat'+mesonCat+year+'.root')

#PHI
#if mesonCat == "_PhiCat":
#   mytreeGF.Add(directory+'outname_mc1017_GFcat'+mesonCat+year+'.root')
#   mytreeGFfix.Add(directory+'outname_mc1052_GFcat'+mesonCat+year+'.root')

#mytreeGF.Add(directory+'outname_mc-65_GFcat'+mesonCat+year+'.root')
#mytreeGF.Add(directory+'outname_mc-66_GFcat'+mesonCat+year+'.root')
#stringCat = 'GFcat'

#mytreeGF.Add(directory+'outname_mc-65_Zinvcat'+mesonCat+year+'.root')
#mytreeGF.Add(directory+'outname_mc-66_Zinvcat'+mesonCat+year+'.root')
#stringCat = 'Zinvcat'

#mytreeGF.Add(directory+'outname_mc-65_VBFcatlow'+mesonCat+year+'.root')
#mytreeGF.Add(directory+'outname_mc-66_VBFcatlow'+mesonCat+year+'.root')
#stringCat = 'VBFcatlow'

mytreeGF.Add(directory+'outname_mc-65_VBFcat'+mesonCat+year+'.root')
mytreeGF.Add(directory+'outname_mc-66_VBFcat'+mesonCat+year+'.root')
stringCat = 'VBFcat'

def projection(doTestMVA, mytree, xmin, xmax, rangeLow, rangeHigh, color ):

   nBins = 10
   
   h = ROOT.TH1F('h', ' ', nBins, xmin, xmax)
   if mesonCat == "_RhoCat": h.GetXaxis().SetTitle("m_{#gamma#rho}^{H} [GeV]")
   if mesonCat == "_PhiCat": h.GetXaxis().SetTitle("m_{#gamma#phi}^{H} [GeV]")   
   
   for ev in mytree:
      if ev.MVAdisc[0]>rangeLow and ev.MVAdisc[0]<rangeHigh:
         h.Fill( ev.HCandMass , ev.w )
      
   h.SetLineColor(color)
   h.SetLineWidth(3)   
   return h


if __name__ == "__main__":

   hMH_bin1 = projection(1 , mytreeGF, 100., 170., -1., -0.5, 1)
   hMH_bin2 = projection(1 , mytreeGF, 100., 170., -0.5, 0., 2)
   hMH_bin3 = projection(1 , mytreeGF, 100., 170., 0., 0.5, 3 )
   hMH_bin4 = projection(1 , mytreeGF, 100., 170., 0.5, 1., 4 )
#   hMH_bin5 = projection(1 , mytreeGF, 100., 170., 0.75, 1., 6 )

   c = ROOT.TCanvas("c", "", 600, 600)
   
   hMH_bin1.Scale(1./hMH_bin1.Integral())
   hMH_bin1.SetMaximum(2*hMH_bin1.GetMaximum())
   hMH_bin1.Draw("hist")
   hMH_bin2.DrawNormalized("hist sames")
   hMH_bin3.DrawNormalized("hist sames")
   hMH_bin4.DrawNormalized("hist sames")
#   hMH_bin5.DrawNormalized("hist sames")

   text = ROOT.TText()
   text.SetTextFont(1)
   text.SetTextColor(1)
   text.SetTextSize(0.03)
   text.SetTextAlign(22)
   text.SetTextAngle(0)
   
   text.DrawText(140, 0.18, "MVA[-1,-0.5]");
   text.SetTextColor(2)
   text.DrawText(140, 0.19, "MVA[-0.5,0.]");
   text.SetTextColor(3)
   text.DrawText(140, 0.20, "MVA[0.,0.5]");
   text.SetTextColor(4)
   text.DrawText(140, 0.21, "MVA[0.5,1.]");
#   text.SetTextColor(6)
#   text.DrawText(140, 0.22, "MVA[0.75,1.]");

   text.SetTextColor(1)   
   text.DrawText(110, 0.24, stringCat);


   c.SaveAs("MassHiggs_inBinsOfPT_"+stringCat+"_Rho.png")
