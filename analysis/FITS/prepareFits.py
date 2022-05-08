import ROOT
import os
from array import array

ROOT.gStyle.SetOptStat(0)

#year = '_2018'
directory = '/home/submit/mariadlf/Hrare/analysis/MAY5/2018/'

from LoadTree import loadTree

ROOT.gStyle.SetOptStat(0)

lumis={
    '_12016': 19.52, #APV #(B-F for 2016 pre)
    '_22016': 16.80, #postVFP
    '_2016': 35.9,
    '_2017': 36.4, #41.5, #(C,D,E,F for 2017)
    '_2018': 59.70,
    '_12018': 45.70,
}

# Create the plot

def getHisto(item, nbin, low, high, doLog,category,mesonCat, doSignal, nameSig):

   year = '_2018'
   mytree = ROOT.TChain('events')
   if doSignal:
      if(category =='_GFcat' and mesonCat == '_PhiCat' and nameSig=="ggH"): mytree.Add(directory+'outname_mc1017'+category+mesonCat+year+'.root') # VBF
      if(category =='_GFcat' and mesonCat == '_PhiCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1010'+category+mesonCat+year+'.root') # VBF
      if(category =='_VBFcat' and mesonCat == '_PhiCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1010'+category+mesonCat+year+'.root') # VBF
      if(category =='_VBFcatlow' and mesonCat == '_PhiCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1010'+category+mesonCat+year+'.root') # VBF
      if(category =='_Wcat' and mesonCat == '_PhiCat' and nameSig=="WH"): mytree.Add(directory+'outname_mc1011'+category+mesonCat+year+'.root') # Wp
      if(category =='_Wcat' and mesonCat == '_PhiCat' and nameSig=="WH"): mytree.Add(directory+'outname_mc1012'+category+mesonCat+year+'.root') # Wm
      if(category =='_Wcat' and mesonCat == '_PhiCat' and nameSig=="ZH"): mytree.Add(directory+'outname_mc1013'+category+mesonCat+year+'.root') # Z
      if(category =='_Zcat' and mesonCat == '_PhiCat' and nameSig=="ZH"): mytree.Add(directory+'outname_mc1013'+category+mesonCat+year+'.root') # Z
      if(category =='_Zinvcat' and mesonCat == '_PhiCat' and nameSig=="ZinvH"): mytree.Add(directory+'outname_mc1015'+category+mesonCat+year+'.root') # Zinv
      if(category =='_Zinvcat' and mesonCat == '_PhiCat' and nameSig=="WH"): mytree.Add(directory+'outname_mc1011'+category+mesonCat+year+'.root') # Zinv Wp
      if(category =='_Zinvcat' and mesonCat == '_PhiCat' and nameSig=="WH"): mytree.Add(directory+'outname_mc1012'+category+mesonCat+year+'.root') # Zinv Wp
      #
      if(category =='_GFcat' and mesonCat == '_RhoCat' and nameSig=="ggH"): mytree.Add(directory+'outname_mc1027'+category+mesonCat+year+'.root') # VBF
      if(category =='_GFcat' and mesonCat == '_RhoCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1020'+category+mesonCat+year+'.root') # VBF
      if(category =='_VBFcat' and mesonCat == '_RhoCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1020'+category+mesonCat+year+'.root') # VBF
      if(category =='_VBFcatlow' and mesonCat == '_RhoCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1020'+category+mesonCat+year+'.root') # VBF
      if(category =='_Wcat' and mesonCat == '_RhoCat' and nameSig=="WH"): mytree.Add(directory+'outname_mc1021'+category+mesonCat+year+'.root') # Wp
      if(category =='_Wcat' and mesonCat == '_RhoCat' and nameSig=="WH"): mytree.Add(directory+'outname_mc1022'+category+mesonCat+year+'.root') # Wm
      if(category =='_Wcat' and mesonCat == '_RhoCat' and nameSig=="ZH"): mytree.Add(directory+'outname_mc1023'+category+mesonCat+year+'.root') # Z
      if(category =='_Zcat' and mesonCat == '_RhoCat' and nameSig=="ZH"): mytree.Add(directory+'outname_mc1023'+category+mesonCat+year+'.root') # Z
      if(category =='_Zinvcat' and mesonCat == '_RhoCat' and nameSig=="ZinvH"): mytree.Add(directory+'outname_mc1025'+category+mesonCat+year+'.root') # Zinv
      if(category =='_Zinvcat' and mesonCat == '_RhoCat' and nameSig=="WH"): mytree.Add(directory+'outname_mc1021'+category+mesonCat+year+'.root') # Zinv
      if(category =='_Zinvcat' and mesonCat == '_RhoCat' and nameSig=="WH"): mytree.Add(directory+'outname_mc1022'+category+mesonCat+year+'.root') # Zinv
   else:
      mytree = loadTree(mytree, directory , category, mesonCat, year ) # all sigm, BKG , data, 

   h = ROOT.TH1F( 'Higgs', '', nbin, low, high )

   print(item)

   for ev in mytree:

      if item == 4 :
         var = ev.HCandMass

      idxMeson = ev.index_pair[0]
      idxPh = ev.index_pair[1]

      if ev.goodMeson_iso[idxMeson] < 0.9 :  continue
      # checked for the VBF as well
      leadTrkPt = max(ev.goodMeson_trk1_pt[idxMeson],ev.goodMeson_trk2_pt[idxMeson])
      if leadTrkPt < 20: continue

      ## OPTIMIZED PHASE SPACE
      if(category =='_Zcat' or category =='_Wcat'):
         if ev.goodPhotons_pt[idxPh]<40 :  continue
         if ev.goodMeson_pt[idxMeson]<40 : continue

      ## OPTIMIZED PHASE SPACE
      if(category =='_Zinvcat'):
         if ev.goodPhotons_pt[idxPh]<40 :  continue
         if ev.goodMeson_pt[idxMeson]<40 : continue
#         if min(ev.goodMeson_trk1_pt[idxMeson],ev.goodMeson_trk2_pt[idxMeson]) < 10: continue
         if ev.DeepMETResolutionTune_pt<75 :  continue
         if abs(ev.dPhiGammaMET)<1. : continue
         if abs(ev.dPhiMesonMET)<1. : continue
         if ev.ptRatioMEThiggs>0.8: continue

      ## OPTIMIZED PHASE SPACE
      if(category =='_VBFcatlow'):
         if ev.goodPhotons_pt[idxPh]<40 :  continue # already applied
         if ev.goodPhotons_pt[idxPh]>75 :  continue # already applied
         if ev.goodMeson_pt[idxMeson]<40 : continue

      ## OPTIMIZED PHASE SPACE
      if(category =='_GFcat'):
         if ev.goodPhotons_pt[idxPh]<40 :  continue
         if ev.goodMeson_pt[idxMeson]<40 : continue
#         if ev.nGoodJets < 2: continue

      if mesonCat == '_RhoCat' and abs(ev.goodMeson_mass[idxMeson]-0.77) > 2.*0.07: continue # 2 sigma
      if mesonCat == '_PhiCat' and abs(ev.goodMeson_mass[idxMeson]-1.02) > 3.*0.004: continue # 3. sigma

      # Fill histograms
      if (doSignal) :
         if (nameSig=='WH' and (ev.mc==1011 or ev.mc==1012 or ev.mc==1021 or ev.mc==1022)): #W
            h.Fill( var, ev.w )
         if (nameSig=='ZH' and (ev.mc==1013 or ev.mc==1023)): #Z
            h.Fill( var, ev.w )
         if ((nameSig=='VBFH' or nameSig=='VBFHlow') and (ev.mc==1010 or ev.mc==1020)): #VBF
            h.Fill( var, ev.w )
         if (nameSig=='ZinvH' and (ev.mc==1015 or ev.mc==1025)): #VBF
            h.Fill( var, ev.w )
         if (nameSig=='ggH' and (ev.mc==1017 or ev.mc==1027)): #VBF
            h.Fill( var, ev.w )
      else :
         if ev.mc<0:   # only DATA
            h.Fill( var, ev.w )

   #Loop over tree done
   if ev.mc>=0:
      if(category == '_Zinvcat' or category == '_GFcat' or category == '_VBFcatlow'): year = '_12018'
      print("lumi = ",lumis[year])
      h.Scale(lumis[year])

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
   text.DrawLatex(0.65, 0.93, "#sqrt{s} = 13 TeV,%0.2f fb^{-1}"% (lumis[year]))

   string = category+mesonCat

   if item==4:
      c.SaveAs("HCandMass"+string+".png")
      print("HCandmass.png")
      
if __name__ == "__main__":

   plot(4, 200, 0. , 200.,True,'_Wcat','_PhiCat') # HCandMass
   plot(4, 200, 0. , 200.,True,'_VBFcat','_PhiCat') # HCandMass
#   plot(4, 200, 0. , 200.,True,'_VBFcat','_RhoCat') # HCandMass
#   plot(4, 200, 0. , 200.,True,'_VBFcat','_PhiCat') # HCandMass
#   plot(4, 200, 0. , 200.,True,'_Wcat','_RhoCat') # HCandMass
