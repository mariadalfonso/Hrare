import ROOT
import os
from array import array
import math

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gROOT.SetBatch()

#mesonCat = '_PhiCat'
mesonCat = '_RhoCat'
#anaCat = '_GFcat'
#anaCat = '_VBFcat'
anaCat = '_VBFcatlow'
#anaCat = '_Wcat'
#anaCat = '_Zcat'

year = '2018'

mytree = ROOT.TChain('events')
#mytree.Add('../TRGstudy/2018/outname_mc-4_Wcat'+mesonCat+'_2018.root')
if mesonCat == "_RhoCat" and anaCat == '_GFcat': mytree.Add('../TRGstudy/'+year+'/outname_mc1027'+anaCat+mesonCat+'_'+year+'.root')
if mesonCat == "_PhiCat" and anaCat == '_GFcat': mytree.Add('../TRGstudy/'+year+'/outname_mc1017'+anaCat+mesonCat+'_'+year+'.root')

if mesonCat == "_RhoCat" and (anaCat == '_VBFcat' or anaCat == '_VBFcatlow'): mytree.Add('../TRGstudy/'+year+'/outname_mc1020'+anaCat+mesonCat+'_'+year+'.root')
if mesonCat == "_PhiCat" and (anaCat == '_VBFcat' or anaCat == '_VBFcatlow'): mytree.Add('../TRGstudy/'+year+'/outname_mc1010'+anaCat+mesonCat+'_'+year+'.root')

if mesonCat == "_RhoCat" and anaCat == '_Wcat': mytree.Add('../TRGstudy/'+year+'/outname_mc1021'+anaCat+mesonCat+'_'+year+'.root')
if mesonCat == "_PhiCat" and anaCat == '_Wcat': mytree.Add('../TRGstudy/'+year+'/outname_mc1011'+anaCat+mesonCat+'_'+year+'.root')

if mesonCat == "_RhoCat" and anaCat == '_Zcat': mytree.Add('../TRGstudy/'+year+'/outname_mc1023'+anaCat+mesonCat+'_'+year+'.root')
if mesonCat == "_PhiCat" and anaCat == '_Zcat': mytree.Add('../TRGstudy/'+year+'/outname_mc1013'+anaCat+mesonCat+'_'+year+'.root')

def plot(itemTRG, item, mytree):

   print(item)
   nEntries = mytree.GetEntries()
   print('nEntries=',nEntries)

   xmin = 0.
   xmax = 1000.
   nBins = int(xmax-xmin)

   xmin2 = 0.
   xmax2 = 400.
   if anaCat == '_Wcat': xmax2 = 200.
   if anaCat == '_Zcat': xmax2 = 100.
   if anaCat == '_Wcat' or anaCat == '_Zcat': xmin2 = 20.
   if anaCat == '_GFcat' or anaCat == '_VBFcatlow': xmin2 = 40.
   if anaCat == '_GFcat': xmax2 = 350.
   if anaCat == '_VBFcatlow': xmax2 = 75.
   if anaCat == '_VBFcat': xmax2 = 500.
   if anaCat == '_VBFcat' and item == 5: xmin2 = 75.
   if anaCat == '_VBFcat' and item == 6: xmin2 = 300.   
   nBins2 = int(xmax2-xmin2)

   h = ROOT.TH1F('h', ' ', nBins, xmin, xmax)
   h1 = ROOT.TH1F('h', ' ', nBins2, xmin2, xmax2)   

   for ev in mytree:

      offline=False
      if item == 1:
         var = ev.photon_pt
         offline = (abs(ev.HCandMass-125)<25)
      elif item == 2:
         var = ev.meson_pt
         offline = (abs(ev.HCandMass-125)<25)
      elif item == 3:
         #muon
         var = ev.LeadingLepton
#         if anaCat == '_Zcat': var = ev.SubLeadingLepton
         offline = (ev.isMuorEle==1 and abs(ev.HCandMass-125)<25)
      elif item == 4:
         #ele
         var = ev.LeadingLepton
#         if anaCat == '_Zcat': var = ev.SubLeadingLepton         
         offline = (ev.isMuorEle==2 and abs(ev.HCandMass-125)<25)         
      elif item == 5:
         #VBFtrigger
         var = ev.photon_pt
         offline = ev.mJJ>300 and ev.dEtaJJ>3 and ev.Y1Y2<0 and abs(ev.HCandMass-125)<25
      elif item == 6:
         #VBFtrigger
         var = ev.mJJ
         offline = ev.mJJ>300 and ev.dEtaJJ>3 and ev.Y1Y2<0 and abs(ev.HCandMass-125)<25         
         
      if itemTRG == 0 and offline: h.Fill( var , ev.w )
      if itemTRG == 1 and ev.triggerAna>0 and offline: h.Fill( var , ev.w )
      if itemTRG == 2 and ev.triggerPhoton>0 and offline: h.Fill( var , ev.w )      
      if itemTRG == 3 and ev.triggerPhotonIso>0 and offline: h.Fill( var , ev.w )      
#      if item == 3 and ev.triggerPhoton>0 and offline: h.Fill( var , ev.w )
#      if item == 2 and ev.triggerVBF>0 and offline: h.Fill( var , ev.w )

   for myBin in range(0,h1.GetNbinsX()+1):  
      h1.SetBinContent(myBin,h.Integral(h.FindBin(h1.GetBinLowEdge(myBin)),h.GetNbinsX()+1))
   
   if (anaCat == '_GFcat' or anaCat == '_VBFcatlow') and mesonCat == '_RhoCat': h1.GetXaxis().SetTitle("p_{T} ( #gamma or #rho ) ")
   if (anaCat == '_GFcat' or anaCat == '_VBFcatlow') and mesonCat == '_PhiCat': h1.GetXaxis().SetTitle("p_{T} ( #gamma or #phi ) ")   
   h1.GetYaxis().SetTitle("efficiency")   
   if anaCat == '_VBFcat': h1.GetXaxis().SetTitle("p_{T}^{#gamma} or M_{jj}")
   if anaCat == '_Wcat': h1.GetXaxis().SetTitle("p_{T} ( leading lepton ) ")
#   if anaCat == '_Zcat': h1.GetXaxis().SetTitle("p_{T} ( subleading lepton ) ")
   if anaCat == '_Zcat': h1.GetXaxis().SetTitle("p_{T} ( leading lepton ) ")   
   h1.SetLineWidth(3)
   
   if anaCat == '_Zcat':
      if mesonCat == "_PhiCat": h1.SetTitle("ZH, H #rightarrow #phi #gamma")
      if mesonCat == "_RhoCat": h1.SetTitle("ZH, H #rightarrow #rho #gamma")
   if anaCat == '_Wcat':
      if mesonCat == "_PhiCat": h1.SetTitle("WH, H #rightarrow #phi #gamma")
      if mesonCat == "_RhoCat": h1.SetTitle("WH, H #rightarrow #rho #gamma")
   if anaCat == '_GFcat':
      if mesonCat == "_PhiCat": h1.SetTitle("H #rightarrow #phi #gamma")
      if mesonCat == "_RhoCat": h1.SetTitle("H #rightarrow #rho #gamma")
   if anaCat == '_VBFcatlow' or anaCat == '_VBFcat':
      if mesonCat == "_PhiCat": h1.SetTitle("qqH, H #rightarrow #phi #gamma")
      if mesonCat == "_RhoCat": h1.SetTitle("qqH, H #rightarrow #rho #gamma")      


   return h1


if __name__ == "__main__":

   if anaCat == '_Wcat' or anaCat == '_Zcat':

      #Muon nominal trigger soup
      h_den1 = plot(0, 3, mytree) # Hmass
      h_num1 = plot(1, 3, mytree) # Hmass

      #Electron nominal trigger soup
      h_den2 = plot(0, 4, mytree) # Hmass
      h_num2 = plot(1, 4, mytree) # Hmass      
      
   elif anaCat == '_GFcat' or anaCat == '_VBFcatlow':

      #triggerAna Photon
      h_den1 = plot(0, 1, mytree) # Hmass
      h_num1 = plot(1, 1, mytree) # Hmass

      #triggerAna Meson
      h_den2 = plot(0, 2, mytree) # Hmass
      h_num2 = plot(1, 2, mytree) # Hmass

   elif  anaCat == '_VBFcat':

      #triggerAna Photon
      h_den1 = plot(0, 5, mytree) # Hmass
      h_num1 = plot(1, 5, mytree) # Hmass

      #triggerAna Meson
      h_den2 = plot(0, 6, mytree) # Hmass
      h_num2 = plot(1, 6, mytree) # Hmass      

   else:
            
      #trigger Iso photon
      h_den3 = plot(0, 1, mytree) # Hmass
      h_num3 = plot(2, 1, mytree) # Hmass
      
      #trigger 200 Photon
      h_den4 = plot(0, 1, mytree) # Hmass
      h_num4 = plot(3, 1, mytree) # Hmass         

   
   c = ROOT.TCanvas("c", "", 600, 600)
#   c.SetLogy(1)
   pad = ROOT.TPad("upper_pad", "", 0, 0, 1, 1)

   print('integral den = ', h_den1.Integral() )
   print('integral num (TRG Ana) = ', h_num1.Integral() )
   print('integral den2 = ', h_den2.Integral() )
   print('integral num2 = ', h_num2.Integral() )      

   h_num1.Divide(h_den1)
   h_num2.Divide(h_den2)
   h_num1.SetLineColor(ROOT.kBlue)
   h_num2.SetLineColor(ROOT.kRed)

   h_num1.GetYaxis().SetRangeUser(0.,1.)
   h_num1.Draw("hist sames")
   if anaCat == '_GFcat' or anaCat == '_VBFcatlow' or anaCat == '_VBFcat': h_num2.Draw("hist sames")
   if anaCat == '_Wcat' or anaCat == '_Zcat': h_num2.Draw("hist sames")
   if anaCat == '_Wcat' or anaCat == '_Zcat':
      testELE = ROOT.TLine(30, 0, 30, 1.)
      testELE.SetLineColor(ROOT.kRed)
      testELE.SetLineStyle(3)   
      testELE.Draw()
      testMU = ROOT.TLine(20, 0, 20, 1.)
      testMU.SetLineColor(ROOT.kBlue)
      testMU.SetLineStyle(3)   
      testMU.Draw()   

   latex = ROOT.TLatex()
   latex.SetTextColor(ROOT.kBlue)
   latex.SetTextSize(0.04)
   if anaCat == '_GFcat': latex.DrawLatex(130, 0.2, " Photon35_TwoProng35 (vs #gamma)")
   if anaCat == '_VBFcatlow': latex.DrawLatex(45, 0.2, " Photon35_TwoProng35 (vs #gamma)")   
   if anaCat == '_Wcat' or anaCat == '_Zcat': latex.DrawLatex(70, 0.2, " muon OR")
   latex.SetTextColor(ROOT.kRed)   
   if anaCat == '_GFcat' and mesonCat == '_RhoCat': latex.DrawLatex(130, 0.3, " Photon35_TwoProng35 (vs #rho)")
   if anaCat == '_GFcat' and mesonCat == '_PhiCat': latex.DrawLatex(130, 0.3, " Photon35_TwoProng35 (vs #phi)")
   if anaCat == '_VBFcatlow' and mesonCat == '_RhoCat': latex.DrawLatex(45, 0.3, " Photon35_TwoProng35 (vs #rho)")
   if anaCat == '_VBFcatlow' and mesonCat == '_PhiCat': latex.DrawLatex(45, 0.3, " Photon35_TwoProng35 (vs #phi)")      
   if anaCat == '_Wcat' or anaCat == '_Zcat': latex.DrawLatex(70, 0.3, " electron OR")

   latex.SetTextColor(ROOT.kBlack)      
   latex.SetTextSize(0.03)
   if anaCat == '_VBFcat': latex.DrawLatex(100, 0.2, "Photon75_*_EBOnly_*MJJ*")
   if anaCat == '_VBFcat': latex.DrawLatex(100, 0.15, "|| Photon35_TwoProng35")
   
   if anaCat == '_GFcat' or anaCat == '_VBFcatlow':
      c.SaveAs('~/public_html/Hrare/TRG/PhotonMesonPt'+anaCat+mesonCat+year+'.png')
   if anaCat == '_VBFcat':
      c.SaveAs('~/public_html/Hrare/TRG/PhotonPtMjj'+anaCat+mesonCat+year+'.png')               
   elif anaCat == '_Zcat':
      c.SaveAs('~/public_html/Hrare/TRG/LeptonPt'+anaCat+mesonCat+year+'.png')      
#      c.SaveAs('~/public_html/Hrare/TRG/SubLeptonPt'+anaCat+mesonCat+year+'.png')
   elif anaCat == '_Wcat':
      c.SaveAs('~/public_html/Hrare/TRG/LeptonPt'+anaCat+mesonCat+year+'.png')      
      
   exit()

   h_num3.Divide(h_den3)
   h_num4.Divide(h_den4)      

   h_num3.SetLineColor(ROOT.kGreen+1)
   h_num4.SetLineColor(ROOT.kMagenta+1)
#   h_den.Draw("hist")
   h_num3.Draw("hist sames")
   h_num4.Draw("hist sames")         

   latex = ROOT.TLatex()
   latex.SetTextColor(ROOT.kBlue)
   latex.SetTextSize(0.04)
   if anaCat == '_GFcat': latex.DrawLatex(130, 0.2, " Photon35_TwoProng35 (vs ph)")
   if anaCat == '_VBFcat': latex.DrawLatex(130, 0.2, " soup (PH75+VBF & PH35+TwoProng35) (vs ph)")
   latex.SetTextColor(ROOT.kRed)
   if anaCat == '_GFcat': latex.DrawLatex(130, 0.3,  " Photon35_TwoProng35 (vs meson)")
   latex.SetTextColor(ROOT.kGreen+1)
   latex.DrawLatex(130, 0.4,  " Photon200")
   latex.SetTextColor(ROOT.kMagenta+1)
   latex.DrawLatex(130, 0.5,  " Photon50_R9Id90_HE10_IsoM")   

#   latex.SetTextColor(ROOT.kGreen+1)
#   latex.DrawLatex(130, 0.2*h_den.GetMaximum(),  "mu TRG + Hcand + Photon50_R9Id90_HE10_IsoM")

#   bmin = hGF_org.FindBin(115)
#   bmax = hGF_org.FindBin(135)
#   integralSig = hGF_org.Integral(hGF_org.FindBin(115),hGF_org.FindBin(135))

#   bmin = hGF_fix.FindBin(115)
#   bmax = hGF_fix.FindBin(135)
#   integralSig2 = hGF_fix.Integral(hGF_fix.FindBin(115),hGF_fix.FindBin(135))

#   print(" ORG: ",integralSig," FIX: ",integralSig2)

   c.SaveAs('~/public_html/Hrare/TRG/PhotonPt'+anaCat+mesonCat+year+'.png')
