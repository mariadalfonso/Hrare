import ROOT
import os
from array import array
import math
import sys

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gROOT.SetBatch()

category = sys.argv[1]
mesonCat = sys.argv[2]
#category = '_GFcat'
#mesonCat = '_RhoCat'
directory = '/work/submit/mariadlf/genStudies/'
year = '_2018'

doPull=False

mytree = ROOT.TChain('events')

if category=='_GFcat':
   if mesonCat == '_PhiCat': mytree.Add(directory+'miniTreeForGNN_mc1017.root') #ggH phi
   if mesonCat == '_RhoCat': mytree.Add(directory+'miniTreeForGNN_mc1027.root') #ggH rho   
elif category=='_VBFcat':
   if mesonCat == '_PhiCat': mytree.Add(directory+'miniTreeForGNN_mc1010.root') #ggH phi
   if mesonCat == '_RhoCat': mytree.Add(directory+'miniTreeForGNN_mc1020.root') #ggH rho
elif category=='_Wcat':
   if mesonCat == '_PhiCat': mytree.Add(directory+'miniTreeForGNN_mc1011.root') #ggH phi
   if mesonCat == '_RhoCat': mytree.Add(directory+'miniTreeForGNN_mc1021.root') #ggH rho
elif category=='_Zcat':
   if mesonCat == '_PhiCat': mytree.Add(directory+'miniTreeForGNN_mc1013.root') #ggH phi
   if mesonCat == '_RhoCat': mytree.Add(directory+'miniTreeForGNN_mc1023.root') #ggH rho            

def plotPull(item, nbin, low, high, doLog):

   h = ROOT.TH1F( 'h', ' ', nbin, low, high)

   print(item)
   nEntries = mytree.GetEntries()
   print('nEntries=',nEntries)
   
   for ev in mytree:

      if item == 1 and mesonCat=='_PhiCat':
         if ev.index_Phi[0]!= -1:
            if doPull: 
               var = (ev.goodPhiFromHiggs_mass-ev.matchedPhiGenPart_mass)/ev.goodPhiFromHiggs_massErr
               h.GetXaxis().SetTitle("m^{RECO}_{#phi}-m^{GEN}_{#phi}/#sigmam^{RECO}_{#phi}")
            else:
               var = (ev.goodPhiFromHiggs_mass-ev.matchedPhiGenPart_mass)
               h.GetXaxis().SetTitle("m^{RECO}_{#phi}-m^{GEN}_{#phi}")
            h.Fill( var, ev.w)
            
      if item == 1 and mesonCat=='_RhoCat':
         if ev.index_Rho[0]!= -1:
            if doPull: 
               var = (ev.goodRhoFromHiggs_mass-ev.matchedRhoGenPart_mass)/ev.goodRhoFromHiggs_massErr            
               h.GetXaxis().SetTitle("m^{RECO}_{#rho}-m^{GEN}_{#rho}/#sigmam^{RECO}_{#rho}")
            else:
               var = (ev.goodRhoFromHiggs_mass-ev.matchedRhoGenPart_mass) 
               h.GetXaxis().SetTitle("m^{RECO}_{#rho}-m^{GEN}_{#rho}")            
            h.Fill( var, ev.w)

      if item == 1 and mesonCat=='_KsCat':
         if ev.index_Ks[0]!= -1:
            var = (ev.goodKsFromHiggs_mass-ev.matchedKsGenPart_mass)/ev.goodKsFromHiggs_massErr            
            h.GetXaxis().SetTitle("m^{RECO}_{K_{S}}-m^{GEN}_{K_{S}}/#sigmam^{RECO}_{K_{S}}")
            h.Fill( var, ev.w)            

   # Create canvas with pad
   c = ROOT.TCanvas("c", "", 600, 600)
   pad = ROOT.TPad("upper_pad", "", 0, 0, 1, 1)
#   gStyle.SetOptStat(0)
   if doLog : pad.SetLogy(1)
   pad.SetTickx(False)
   pad.SetTicky(False)
   pad.Draw()
   pad.cd()

   h.Draw("HIST")
   if doPull:
      if mesonCat == '_RhoCat': h.Fit("gaus","E","",-1.5,1.5)
      if mesonCat == '_PhiCat': h.Fit("gaus","E","",-2.5,2.5)
   else:
      if mesonCat == '_RhoCat': h.Fit("gaus","E","",-0.01,0.01)
      if mesonCat == '_PhiCat': h.Fit("gaus","E","",-0.005,0.005)   
      h.GetFunction("gaus").SetLineColor(4);
      
   h.Draw("pe")
      
   # Add label
   text = ROOT.TLatex()
   text.SetNDC()
   text.SetTextFont(72)
   text.SetTextSize(0.045)
   text.DrawLatex(0.15, 0.93, "CMS")
   text.SetTextFont(42)
   text.DrawLatex(0.15 + 0.10, 0.93, "Simulation")
   text.SetTextSize(0.04)
#   text.DrawLatex(0.65, 0.93, "#sqrt{s} = 13 TeV,%0.2f fb^{-1}"% (lumi))   

   if category=='_GFcat': text.DrawLatex(0.25, 0.7, "ggH")
   if category=='_VBFcat': text.DrawLatex(0.25, 0.7, "VBF")
   if category=='_Wcat': text.DrawLatex(0.25, 0.7, "WH")
   if category=='_Zcat': text.DrawLatex(0.25, 0.7, "ZH")

   if item==1:
      if doPull: c.SaveAs("PullMass_"+mesonCat+"_"+category+".png")
      if not doPull: c.SaveAs("MassRes_"+mesonCat+"_"+category+".png")      

if __name__ == "__main__":

#   plotPull(1, 100, -5. , 5., False) #pull   
#   plotPull(1, 100, -0.05 , 0.05, False) #Mass Rho
   plotPull(1, 100, -0.02 , 0.02, False) #Mass Phi
   

   exit()

'''   
python PullPlots.py '_GFcat' '_RhoCat'
python PullPlots.py '_VBFcat' '_RhoCat'
python PullPlots.py '_Wcat' '_RhoCat'
python PullPlots.py '_Zcat' '_RhoCat'

python PullPlots.py '_GFcat' '_PhiCat'
python PullPlots.py '_VBFcat' '_PhiCat'
python PullPlots.py '_Wcat' '_PhiCat'
python PullPlots.py '_Zcat' '_PhiCat'

'''
