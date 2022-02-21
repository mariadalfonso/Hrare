import ROOT
import os
from array import array

ROOT.gStyle.SetOptStat(0)

year = '_2018'
directory = '../FEB9sig/'

#lumi=137
lumi=59.70
#lumi=45.97

# Create the plot

def getHisto(item, nbin, low, high, doLog,category,mesonCat, doSignal):

   mytree = ROOT.TChain('events')
   if doSignal:
      if(category =='_VBFcat' and mesonCat == '_PhiCat'): mytree.Add(directory+'outname_mc1010'+category+mesonCat+year+'.root') # VBF
      if(category =='_Wcat' and mesonCat == '_PhiCat'): mytree.Add(directory+'outname_mc1011'+category+mesonCat+year+'.root') # Wp
      if(category =='_Wcat' and mesonCat == '_PhiCat'): mytree.Add(directory+'outname_mc1012'+category+mesonCat+year+'.root') # Wm
      if(category =='_VBFcat' and mesonCat == '_RhoCat'): mytree.Add(directory+'outname_mc1020'+category+mesonCat+year+'.root') # VBF
   else:
#      if(category =='_VBFcat'): mytree.Add(directory+'outname_mc6'+category+mesonCat+year+'.root') # VBF
#      if(category =='_VBFcat'): mytree.Add(directory+'outname_mc7'+category+mesonCat+year+'.root') # VBF
#      if(category =='_VBFcat'): mytree.Add(directory+'outname_mc8'+category+mesonCat+year+'.root') # VBF
#      if(category =='_VBFcat'): mytree.Add(directory+'outname_mc9'+category+mesonCat+year+'.root') # VBF
      if(category =='_VBFcat'): mytree.Add(directory+'outname_mc-31'+category+mesonCat+year+'.root') # VBF
      if(category =='_VBFcat'): mytree.Add(directory+'outname_mc-32'+category+mesonCat+year+'.root') # VBF
      if(category =='_VBFcat'): mytree.Add(directory+'outname_mc-33'+category+mesonCat+year+'.root') # VBF
      if(category =='_Wcat'): mytree.Add(directory+'outname_mc-1'+category+mesonCat+year+'.root') # VH
      if(category =='_Wcat'): mytree.Add(directory+'outname_mc-2'+category+mesonCat+year+'.root') # VH
      if(category =='_Wcat'): mytree.Add(directory+'outname_mc-3'+category+mesonCat+year+'.root') # VH
      if(category =='_Wcat'): mytree.Add(directory+'outname_mc-11'+category+mesonCat+year+'.root') # VH
      if(category =='_Wcat'): mytree.Add(directory+'outname_mc-12'+category+mesonCat+year+'.root') # VH
      if(category =='_Wcat'): mytree.Add(directory+'outname_mc-13'+category+mesonCat+year+'.root') # VH
      if(category =='_Wcat'): mytree.Add(directory+'outname_mc-21'+category+mesonCat+year+'.root') # VH
      if(category =='_Wcat'): mytree.Add(directory+'outname_mc-22'+category+mesonCat+year+'.root') # VH
      if(category =='_Wcat'): mytree.Add(directory+'outname_mc-23'+category+mesonCat+year+'.root') # VH
   
   h = ROOT.TH1F( 'Higgs', '', nbin, low, high )

   print(item)

   for ev in mytree:
      
      if item == 4 :
         var = ev.HCandMass

      # Fill histograms.
      # data
      if ((category =='_VBFcat' and (ev.mc==-31 or ev.mc==-32 or ev.mc==-33)) or
          (category =='_Wcat' and (ev.mc==-1 or ev.mc==-2 or ev.mc==-3 or ev.mc==-11 or ev.mc==-12 or ev.mc==-13 or ev.mc==-21 or ev.mc==-22 or ev.mc==-23)) ):
      if ev.mc==1011 or ev.mc==1012:
         h.Fill( var, ev.w)
      if ev.mc==1010 or ev.mc==1020:
         h.Fill( var, ev.w)

   return h

def plot(item, nbin, low, high, doLog,category,mesonCat):

   h = getHisto(item, nbin, low, high, doLog,category,mesonCat,False)
   
   # Create canvas with pad
   c = ROOT.TCanvas("c", "", 600, 600)
   pad = ROOT.TPad("upper_pad", "", 0, 0, 1, 1)
#   gStyle.SetOptStat(0)
   if doLog : pad.SetLogy(1)
   pad.SetTickx(False)
   pad.SetTicky(False)
   pad.Draw()
   pad.cd()   
   
   h.GetXaxis().SetTitle("m_{#gamma#phi}^{H} [GeV]")
   h.GetYaxis().SetTitle("Events")
   h.GetYaxis().SetTitleOffset(1.1)   
   h.GetYaxis().SetLabelSize(0.04)
   h.GetYaxis().SetTitleSize(0.045)
   h.GetYaxis().ChangeLabel(1, -1, 0)
   h.Draw("hist")
   
   # Add legend
   if(category =='_Wcat' or category =='_Zcat'): deltay=0.25
   if(category =='_VBFcat'): deltay=0.15
#   legend = ROOT.TLegend(0.63-0.5, 0.88-deltay, 0.85-0.5, 0.88)
#   legend = ROOT.TLegend(0.63, 0.88-deltay, 0.85, 0.88)
   legend = ROOT.TLegend(0.13, 0.88-deltay, 0.35, 0.88)   
   legend.SetTextFont(42)
   legend.SetFillStyle(0)
   legend.SetBorderSize(0)
   legend.SetTextSize(0.04)
   legend.SetTextAlign(32)
 
   if h and category =='_Zcat' and mesonCat == "_PhiCat" and h.Integral()>0: legend.AddEntry(h, "ZH(#gamma#phi)", "f")
   if h and category =='_Zcat' and mesonCat == "_RhoCat" and h.Integral()>0: legend.AddEntry(h, "ZH(#gamma#phi)", "f")
   if h and category =='_Wcat' and mesonCat == "_PhiCat" and h.Integral()>0: legend.AddEntry(h, "WH(#gamma#phi)", "f")
   if h and category =='_Wcat' and mesonCat == "_RhoCat" and h.Integral()>0: legend.AddEntry(h, "WH(#rho#phi)", "f")   
   if h and category =='_VBFcat' and mesonCat == "_PhiCat" and h.Integral()>0: legend.AddEntry(h, "qqH(#gamma#phi)", "f")
   if h and category =='_VBFcat' and mesonCat == "_RhoCat" and h.Integral()>0: legend.AddEntry(h, "qqH(#rho#phi)", "f")
   legend.Draw("SAME")
 
   # Add label
   text = ROOT.TLatex()
   text.SetNDC()
   text.SetTextFont(72)
   text.SetTextSize(0.045)
   text.DrawLatex(0.15, 0.93, "CMS")
   text.SetTextFont(42)
   text.DrawLatex(0.15 + 0.10, 0.93, "Simulation")
   text.SetTextSize(0.04)
   text.DrawLatex(0.65, 0.93, "#sqrt{s} = 13 TeV,%0.2f fb^{-1}"% (lumi))   


   string = category+mesonCat+"_isEle"

   if item==4:
      c.SaveAs("HCandMass"+string+".png")
      print("HCandmass.png")
      
if __name__ == "__main__":

   plot(4, 200, 0. , 200.,True,'_Wcat','_PhiCat') # HCandMass
   plot(4, 200, 0. , 200.,True,'_VBFcat','_PhiCat') # HCandMass
