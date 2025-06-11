import ROOT
import os
from array import array
import math

from LoadTree import loadTree
from LoadTree import resetTree
from LoadTree import checkNan


ROOT.ROOT.EnableImplicitMT()
RDataFrame = ROOT.RDataFrame
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

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

def getFWHM(h):

      max_bin = h.GetMaximumBin()
      maxP = h.GetBinContent(max_bin)
      x_maximum = h.GetBinCenter(max_bin)

      print('x_maximum = ',x_maximum)
      
      half_max = maxP / 2

      # Search left for half max
      fwhm_left_bin = max_bin
      while fwhm_left_bin > 1 and h.GetBinContent(fwhm_left_bin) > half_max:
         fwhm_left_bin -= 1
         
      # Search right for half max
      fwhm_right_bin = max_bin
      while fwhm_right_bin < h.GetNbinsX() and h.GetBinContent(fwhm_right_bin) > half_max:
         fwhm_right_bin += 1
            
      # Convert bin to x value
      fwhm_left_x = h.GetBinCenter(fwhm_left_bin)
      fwhm_right_x = h.GetBinCenter(fwhm_right_bin)
      
      # Calculate width and integral
      width = fwhm_right_x - fwhm_left_x
      myRange = h.Integral(fwhm_left_bin, fwhm_right_bin)

      print("FWHM width =", width, "Integral in range =", myRange)

      return fwhm_left_x,fwhm_right_x

def getHisto(item, nbin, low, high, doLog,category, doSignal):

   print("getHisto getting called")

   year = '_2018'
   dirLOCAL_='/work/submit/mariadlf/Hrare_JPsiCC/MARCH2025/'

   mytree = ROOT.TChain('events')   
   mytree = loadTree(mytree, dirLOCAL_ , '_GFcat', year, doSignal) # all sigm, BKG , data,
   mytree.Add(dirLOCAL_+'snapshotJpsiCC_1000'+year+'.root')
   mytree.Add(dirLOCAL_+'snapshotJpsiCC_10'+year+'.root')   
   mytree.Add(dirLOCAL_+'snapshotJpsiCC_11'+year+'.root')
   
   mytree.SetBranchStatus("*", 0)
   mytree.SetBranchStatus("massHiggsCorr", 1)
   mytree.SetBranchStatus("mc", 1)
   mytree.SetBranchStatus("w", 1)
   mytree.SetBranchStatus("lumiIntegrated", 1)

   ROOT.gROOT.cd() 
   h = ROOT.TH1F("hSig", "hSig", nbin, low, high)
   h.SetDirectory(0)
   h1 = ROOT.TH1F("hBkg", "hBkg", nbin, low, high)
   h1.SetDirectory(0)   

   #-------------- BELOW RDF

   '''
 
   print("building RDF")  
   df = RDataFrame(mytree)
   print("building RDF DONE")     
   if item == 4: var = 'massHiggsCorr'

#   print('going check NAN')
#   checkNan(var,mytree)
#   print('end check NAN')

   varFix = "(TMath::IsNaN(massHiggsCorr) ? 0: massHiggsCorr)"

   if doSignal:
      print('INSIDE doSignal')
      h = df.Define("var","{}".format(var)).Define("weight","w*lumiIntegrated").Filter("(mc==1000)").Histo1D(("hSig ","hSig",nbin, low, high),"var","weight")
      print(' integral Signal == ',h.GetIntegral())
      getFWHM(h)
      
   else:

      h = df.Define("var","{}".format(varFix)).Define("weight","w*lumiIntegrated").Filter("(mc==10 or mc==11)").Histo1D(("hBkg ","hBkg",nbin, low, high),"var","weight")
      print(' integral Bkg == ',h.GetIntegral())

   return h

   '''
   nEntries = mytree.GetEntries()
   print('nEntries=',nEntries)

   for i in range(nEntries):
         mytree.GetEntry(i)
         var = mytree.massHiggsCorr
         wei = mytree.w * mytree.lumiIntegrated
         if doSignal:
               if mytree.mc==1000 and not math.isnan(var):
                     h.Fill( var, wei )
         else:
               if (mytree.mc==10 or mytree.mc==11) and not math.isnan(var):
                     h.Fill( var, wei )               

   return h

def plot(item, nbin, low, high, doLog,category):

   h = getHisto(item, nbin, low, high, doLog, category, False)

   # Create canvas with pad
   c = ROOT.TCanvas("c", "", 600, 600)
   pad = ROOT.TPad("upper_pad", "", 0, 0, 1, 1)
#   gStyle.SetOptStat(0)
   if doLog : pad.SetLogy(1)
   pad.SetTickx(False)
   pad.SetTicky(False)
   pad.Draw()
   pad.cd()   

   if h:
      h.GetXaxis().SetTitle("m_{c,c,jpsi}^{H} [GeV]")
      h.GetYaxis().SetTitle("Events")
      h.GetYaxis().SetTitleOffset(1.1)   
      h.GetYaxis().SetLabelSize(0.04)
      h.GetYaxis().SetTitleSize(0.045)
      h.GetYaxis().ChangeLabel(1, -1, 0)
      h.Draw("hist")

      fwhm_left_x,fwhm_right_x = getFWHM(h)
      
      l = ROOT.TLine(125, 0, 125, 1.)
      l.SetLineColor(1)
      l.Draw("SAME") #forcing lines to be drawn on the same canvas as the histogram
      ROOT.gPad.Update()
      l2 = ROOT.TLine(fwhm_left_x, 0, fwhm_left_x, 1.)
      l2.SetLineColor(2)
      l2.Draw("SAME")
      l3 = ROOT.TLine(fwhm_right_x, 0, fwhm_right_x, 1.)
      l3.SetLineColor(2)      
      l3.Draw("SAME")
      
   # Add legend
#   if(category =='_Wcat' or category =='_Zcat'): deltay=0.25
#   if(category =='_VBFcat'): deltay=0.15
#   legend = ROOT.TLegend(0.63-0.5, 0.88-deltay, 0.85-0.5, 0.88)
#   legend = ROOT.TLegend(0.63, 0.88-deltay, 0.85, 0.88)
#   legend = ROOT.TLegend(0.13, 0.88-deltay, 0.35, 0.88)   
#   legend.SetTextFont(42)
#   legend.SetFillStyle(0)
#   legend.SetBorderSize(0)
#   legend.SetTextSize(0.04)
#   legend.SetTextAlign(32)
 
#   if h and category =='_Zcat' and mesonCat == "_PhiCat" and h.Integral()>0: legend.AddEntry(h, "ZH(#gamma#phi)", "f")
#   if h and category =='_Zcat' and mesonCat == "_RhoCat" and h.Integral()>0: legend.AddEntry(h, "ZH(#gamma#phi)", "f")
#   if h and category =='_Wcat' and mesonCat == "_PhiCat" and h.Integral()>0: legend.AddEntry(h, "WH(#gamma#phi)", "f")
#   if h and category =='_Wcat' and mesonCat == "_RhoCat" and h.Integral()>0: legend.AddEntry(h, "WH(#rho#phi)", "f")   
#   if h and category =='_VBFcat' and mesonCat == "_PhiCat" and h.Integral()>0: legend.AddEntry(h, "qqH(#gamma#phi)", "f")
#   if h and category =='_VBFcat' and mesonCat == "_RhoCat" and h.Integral()>0: legend.AddEntry(h, "qqH(#rho#phi)", "f")
#   legend.Draw("SAME")
 
#   # Add label
#   text = ROOT.TLatex()
#   text.SetNDC()
#   text.SetTextFont(72)
#   text.SetTextSize(0.045)
#   text.DrawLatex(0.15, 0.93, "CMS")
#   text.SetTextFont(42)
#   text.DrawLatex(0.15 + 0.10, 0.93, "Simulation")
#   text.SetTextSize(0.04)
#   text.DrawLatex(0.65, 0.93, "#sqrt{s} = 13 TeV,%0.2f fb^{-1}"% (lumis[year]))

#   string = category+mesonCat

   if item==4:
      c.SaveAs("~/public_html/test.png")
      print("HCandmass.png")

if __name__ == "__main__":

   plot(4, 100, 0. , 250.,False,'_GFcat') # HCandMass
#   plot(4, 200, 0. , 200.,True,'_Wcat','_PhiCat') # HCandMass
#   plot(4, 200, 0. , 200.,True,'_VBFcat','_PhiCat') # HCandMass
#   plot(4, 200, 0. , 200.,True,'_VBFcat','_RhoCat') # HCandMass
#   plot(4, 200, 0. , 200.,True,'_VBFcat','_PhiCat') # HCandMass
#   plot(4, 200, 0. , 200.,True,'_Wcat','_RhoCat') # HCandMass
