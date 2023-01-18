import ROOT
import os
from array import array
import math
import sys

from LoadTree import loadTree

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

category = sys.argv[1]
mesonCat = sys.argv[2]
#year = '_12016'
#year = '_22016'
#year = '_2017'
year = '_2018'
#year = '_all'

if (category=='_GFcat'): directory4 = '/work/submit/mariadlf/Hrare/DEC28/2018/'
if (category=='_Zinvcat'): directory4 = '/work/submit/mariadlf/Hrare/DEC28/2018/'
if (category=='_VBFcatlow'): directory4 = '/work/submit/mariadlf/Hrare/DEC28/2018/'
if (category=='_VBFcat'):
   directory4 = '/work/submit/mariadlf/Hrare/DEC28/2018/'
   directory3 = '/work/submit/mariadlf/Hrare/DEC28/2017/'
   directory1 = '/work/submit/mariadlf/Hrare/DEC28/12016/'
   year = '_all'

if (category=='_Zcat' or category=='_Wcat'):
   directory4 = '/work/submit/mariadlf/Hrare/JAN11/2018/'
   directory3 = '/work/submit/mariadlf/Hrare/JAN11/2017/'
   directory2 = '/work/submit/mariadlf/Hrare/JAN11/22016/'
   directory1 = '/work/submit/mariadlf/Hrare/JAN11/12016/'
   year = '_Run2'

mytree = ROOT.TChain('events')

BTheo=1.
#if (category=='_GFcat'):
#   if mesonCat == '_PhiCat': BTheo = 0.00231  #2.31 * 10e-6  (GF is 10^-3)
#   if mesonCat == '_RhoCat': BTheo = 0.0168  #1.68 * 10e-5
#if (category=='_VBFcat'):
#   if mesonCat == '_PhiCat': BTheo = 5*0.00231  #2.31 * 10e-6  (GF is 5 x 0^-3)
#   if mesonCat == '_RhoCat': BTheo = 0.0168  #1.68 * 10e-5
#if (category=='_VBFcatlow'):
#   if mesonCat == '_PhiCat': BTheo = 0.0231  #2.31 * 10e-6 (VBFlow is 10^-4)
#   if mesonCat == '_RhoCat': BTheo = 5*0.0168  #1.68 * 10e-5
#if (category=='_VBFcat'):
#   if mesonCat == '_PhiCat': BTheo = 5*0.00231  #2.31 * 10e-6  (GF is 10^-3)
#   if mesonCat == '_RhoCat': BTheo = 0.0168  #1.68 * 10e-5
#if (category=='_Zinvcat'):
#   if mesonCat == '_PhiCat': BTheo = 0.0231  #2.31 * 10e-6  (GF is 10^-4)
#   if mesonCat == '_RhoCat': BTheo = 0.168  #1.68 * 10e-5
#print("BTheo=",BTheo)

if( year == '_all' or year == '_Run2' or year == '_2018'): mytree = loadTree(mytree, directory4, category, mesonCat, year ) # add 2018
if( year == '_all' or year == '_Run2' or year == '_2017'): mytree = loadTree(mytree, directory3, category, mesonCat, year ) # add 2017
if( year == '_Run2' or year == '_22016'): mytree = loadTree(mytree, directory2, category, mesonCat, year ) # add 2016
if( year == '_all' or year == '_Run2' or year == '_12016'): mytree = loadTree(mytree, directory1, category, mesonCat, year ) # add 12016

if (category=='_Zinvcat' or category=='_GFcat' or category=='_VBFcatlow'): year = '_12018'
if (category=='_VBFcat' and year == '_2017'): year = '_12017'

lumis={
    '_12016': 19.52, #APV #(B-F for 2016 pre)
    '_22016': 16.80, #postVFP
    '_2016': 35.9,
    '_2017': 41.5,
    '_12017': 7.7, #(F for 2017) for VBF
    '_2018': 59.70,
    '_12018': 39.54,
    '_all': 86.92,      #19.52 + 7.7 + 59.70
    '_Run2': 138.,      #19.52 + 7.7 + 59.70
}

MVAbinRho={
#   '_GFcat': 0.0,
   '_GFcat': 0.2,
   '_Zinvcat': 0.6,
   '_VBFcat': 0.4,
   '_VBFcatlow': 0.5,
}

MVAbinPhi={
#   '_GFcat': 0.0,
   '_GFcat': 0.2,
   '_Zinvcat': 0.7,
   '_VBFcat': 0.5,
   '_VBFcatlow': 0.6,
}

if mesonCat=='_PhiCat': MVAbin = MVAbinPhi
if mesonCat=='_RhoCat': MVAbin = MVAbinRho

def deltaPhi(phi1,phi2):
   result = phi1 - phi2
   if result > float(M_PI): result -= float(2*M_PI)
   if result <= -float(M_PI): result += float(2*M_PI)
   return result


# Create the plot
def createCanvasPads(doLog):

   # Create canvas with pad
   c = ROOT.TCanvas("c", "", 600, 600)
   pad = ROOT.TPad("upper_pad", "", 0, 0, 1, 1)
#   gStyle.SetOptStat(0)
   if doLog : pad.SetLogy(1)
   pad.SetTickx(False)
   pad.SetTicky(False)
   pad.Draw()
   pad.cd()

   ydiv = 0.2
   pad1 = ROOT.TPad("upper_pad", "", 0.,ydiv,1.,1.)
   if doLog : pad1.SetLogy(1)
   pad2 = ROOT.TPad("lower_pad", "", 0.,0.,1.,ydiv)

   pad1.Draw()
   pad2.Draw()

   return c, pad1, pad2

def plot(item, nbin, low, high, doLog):

   hZinv = ROOT.TH1F( 'hZinv', 'This is the distribution', nbin, low, high)
   hDY = ROOT.TH1F( 'hDY', 'This is the distribution', nbin, low, high)
   hZG = ROOT.TH1F( 'hZG', 'This is the distribution', nbin, low, high)
   hW = ROOT.TH1F( 'hW', 'This is the distribution', nbin, low, high)
   hWG = ROOT.TH1F( 'hWG', 'This is the distribution', nbin, low, high)
   hTTG = ROOT.TH1F( 'hTTG', 'This is the distribution', nbin, low, high)
   hTT12L = ROOT.TH1F( 'hTT12L', 'This is the distribution', nbin, low, high)
   hGJet = ROOT.TH1F( 'hGJet', 'This is the distribution', nbin, low, high)
   hVBFGJet = ROOT.TH1F( 'hVBFGJet', 'This is the distribution', nbin, low, high)
   hJet = ROOT.TH1F( 'hJet', 'This is the distribution', nbin, low, high)
   hVBFH = ROOT.TH1F( 'hVBFH', 'This is the distribution', nbin, low, high)
   hggH = ROOT.TH1F( 'hggH', 'This is the distribution', nbin, low, high)
   hZH = ROOT.TH1F( 'hZH', 'This is the distribution', nbin, low, high)
   hZinvH = ROOT.TH1F( 'hZinvH', 'This is the distribution', nbin, low, high)
   hWH = ROOT.TH1F( 'hWH', 'This is the distribution', nbin, low, high)
   hData = ROOT.TH1F( 'hData', '', nbin, low, high)

   print(item)

   for ev in mytree:

      inMass = False
      idxMeson = 0
#      inMass = (abs(ev.HCandMass - 125)<30)
      idxMeson = ev.index_pair[0]
      idxPh = ev.index_pair[1]
#      idxPh = 1
#      idxPh = 0

      ## CHECK EE, EB, BB

#      photonEta = abs(ev.goodPhotons_eta[idxPh])
#      mesonEta = max(abs(ev.goodMeson_trk1_eta[idxMeson]),abs(ev.goodMeson_trk2_eta[idxMeson]))
#      if mesonEta < 1.4: mesString = 'barrel'
#      if mesonEta >= 1.4: mesString = 'endcap'

#      if abs(ev.goodMeson_trk1_eta[idxMeson]) >= 1.4: continue
#      mesString = 'barrel'

#      if(photonEta < 1.4 and mesonEta<1.4): string = 'BB'
#      elif(photonEta > 1.4 and mesonEta>1.4): string = 'EE'
#      elif(photonEta < 1.4 or mesonEta<1.4): string = 'EB'

      ## OPTIMIZED PHASE SPACE
      if(category =='_Zcat' or category =='_Wcat'):
#         if ev.goodPhotons_pt[idxPh]<40 :  continue
#         if ev.goodMeson_pt[idxMeson]<40 : continue
         if ev.photon_pt<40 :  continue
         if ev.meson_pt<40 : continue

      if(category =='_Wcat'):
         if ev.DeepMETResolutionTune_pt<15 :  continue
         if ev.V_mass<15 :  continue

#      if(category =='_Zcat' or category =='_Wcat'):
#         if ev.isMuorEle!=1: continue; # 1 is diMu and 2 is diEle
#         if ev.isMuorEle!=2: continue; # 1 is diMu and 2 is diEle

      ## OPTIMIZED PHASE SPACE
      if(category =='_Zinvcat'):
         if ev.photon_pt<40 :  continue
         if ev.meson_pt<40 : continue
#         if min(ev.goodMeson_trk1_pt[idxMeson],ev.goodMeson_trk2_pt[idxMeson]) < 10: continue
         if item !=13:
            if ev.DeepMETResolutionTune_pt<75 :  continue
#            if ev.DeepMETResolutionTune_pt<50 :  continue
#         if abs(ev.dPhiGammaMET)<1. : continue # already applied
#         if abs(ev.dPhiMesonMET)<1. : continue # already applied
         if item !=22:
            if ev.ptRatioMEThiggs>0.5: continue
         if item !=91:
            if abs(ev.dEtaGammaMesonCand)>2: continue
         if item !=113:
            if ev.nbtag>0: continue

#         isHEM = ev.DeepMETResolutionTune_phi < -0.9 and ev.DeepMETResolutionTune_phi > -1.6
#         if not isHEM: continue

# HEM
#         isHEM = (ev.goodPhotons_eta[idxPh]< -1.39 and ev.goodPhotons_phi[idxPh]>= -1.6 and ev.goodPhotons_phi[idxPh]< -0.9)
#         if not isHEM: continue
#         if ev.run > 319077: continue

      ## OPTIMIZED PHASE SPACE
#      if(category =='_VBFcatlow'):
#         if ev.goodPhotons_pt[idxPh]<40 :  continue # already applied
#         if ev.goodPhotons_pt[idxPh]>75 :  continue # already applied
#         if ev.goodMeson_pt[idxMeson]<40 : continue

      ## OPTIMIZED PHASE SPACE
      if(category =='_VBFcatlow' or category =='_VBFcat'):
         if ev.jet1Pt < 30 : continue
         if ev.jet2Pt < 20 : continue
#         if ev.deltaJetMeson < 1.2 : continue
#         if ev.deltaJetPhoton < 1. : continue
         if category =='_VBFcat':
            if ev.mJJ < 400. : continue
         else:
            if ev.mJJ < 300. : continue

#      ## OPTIMIZED PHASE SPACE
      if(category =='_GFcat'):
         if ev.photon_pt<40 :  continue
         if ev.meson_pt<40 : continue

      if item == 1 :
         var = ev.V_mass
      if item == 114 :
         var = ev.Visr_mass
      if item == 115 :
         var = ev.LeadingLepton
      if item == 116 :
         var = ev.SubLeadingLepton
      if item == 42 :
         var = ev.MVAdisc[0]
      if item == 43:
         var = -1
         if ev.MVAdisc[0]>MVAbin[category]:  var = ev.HCandMass

      if item == 4 :
         var = ev.HCandMass
         #         var = ev.HGenMas
      if item == 41 :
         var = ev.HCandPT
## Photon plots
      if item == 5 :
         var = ev.goodPhotons_pt[idxPh]
      if item == 6 :
         var = ev.goodPhotons_eta[idxPh]
      if item == 205 :
         var = ev.nGoodPhotons
      if item == 206 :
         var = ev.goodPhotons_hoe[idxPh]
      if item == 207 :
         var = ev.goodPhotons_r9[idxPh]
      if item == 208 :
         #         var = ev.goodPhotons_pfRelIso03_all[idxPh]
         var = ev.goodPhotons_pfRelIso03_chg[idxPh]
      if item == 209 :
         var = ev.goodPhotons_sieie[idxPh]
      if item == 101 :
         var = ev.goodPhotons_pixelSeed[idxPh]
## KS plots
      if item == 102 :
         var = ev.goodKs_mass[idxMeson]
      if item == 103 :
         var = ev.goodKs_pt[idxMeson]
      if item == 107 :
         var = ev.goodKs_lxy[idxMeson]
## Meson plots
      if item == 2 :
         var = ev.goodMeson_mass[idxMeson]
      if item == 3 :
         var = ev.goodMeson_pt[idxMeson]
      if item == 30 :
         #leading
         var = max(ev.goodMeson_trk1_pt[idxMeson],ev.goodMeson_trk2_pt[idxMeson])/ev.goodMeson_pt[idxMeson]
      if item == 29 :
         #subleading
         var = min(ev.goodMeson_trk1_pt[idxMeson],ev.goodMeson_trk2_pt[idxMeson])/ev.goodMeson_pt[idxMeson]
      if item == 31 :
         #leading
         var = max(ev.goodMeson_trk1_pt[idxMeson],ev.goodMeson_trk2_pt[idxMeson])
      if item == 32 :
         #subleading
         var = min(ev.goodMeson_trk1_pt[idxMeson],ev.goodMeson_trk2_pt[idxMeson])
      if item == 33 :
         var = ev.goodMeson_trk1_eta[idxMeson]
      if item == 34 :
         var = ev.goodMeson_trk2_eta[idxMeson]
      if item == 35 :
         var = ev.goodMeson_iso[idxMeson]
      if item == 16 :
         var = abs(ev.goodMeson_DR[idxMeson])
      if item == 10 :
         var = (abs(ev.goodMeson_DR[idxMeson])*ev.goodMeson_pt[idxMeson])/(2*ev.goodMeson_mass[idxMeson])
      if item == 7 :
         var = ev.goodMeson_sipPV[idxMeson]
      if item == 8 :
         var = ev.goodMeson_doca_prob[idxMeson]
      if item == 9 :
         var = ev.goodMeson_vtx_prob[idxMeson]/ev.goodMeson_vtx_chi2dof[idxMeson]
      if item == 18 :
         var = ev.goodMeson_vtx_chi2dof[idxMeson]
      if item == 21 :
#         var = ev.goodMeson_massErr[idxMeson]/ev.goodMeson_mass[idxMeson]
         var = ev.goodMeson_massErr[idxMeson]
## event like
      if item == 11 :
         var = ev.SoftActivityJetNjets5
      if item == 12 :
         if(category =='_VBFcat' or category =='_VBFcatlow' or category =='_GFcat' or category =='_Zinvcat'): var = ev.nGoodJets
         else: var = ev.nJet
      if item == 113 :
         if(category =='_Zinvcat'): var = ev.nbtag
      if item == 50 :
         var = ev.deltaLepMeson
      if item == 51 :
         var = ev.deltaJetMeson
      if item == 52 :
         var = ev.deltaJetPhoton
      if item == 53 :
         var = ev.zepVar
      if item == 54 :
         var = ev.detaHigJet1
      if item == 55 :
         var = ev.detaHigJet2
      if item == 56 :
         var = ev.jet1hfsigmaPhiPhi
      if item == 57 :
         var = ev.jet1hfsigmaEtaEta
      if item == 58 :
         var = ev.jet2hfsigmaPhiPhi
      if item == 59 :
         var = ev.jet2hfsigmaEtaEta
      if item == 13 :
#         var = ev.MET_pt
         var = ev.DeepMETResolutionTune_pt
      if item == 36 :
#         var = ev.MET_pt
         var = ev.DeepMETResolutionTune_phi
      if item == 37 :
#         var = ev.MET_pt
         var = ev.HCandPHI
      if item == 14 :
         var = abs(ev.dPhiGammaMET)
      if item == 15 :
         var = abs(ev.dPhiMesonMET)
      if item == 17 :
         var = ev.mJJ
      if item == 19 :
         var = ev.dEtaJJ
      if item == 20 :
         var = ev.Y1Y2
      if item == 22:
         var = ev.ptRatioMEThiggs
      if item == 23:
         var = ev.jet1Pt
      if item == 24:
         var = ev.jet2Pt
      if item == 231:
         var = ev.jet1Eta
      if item == 241:
         var = ev.jet2Eta
      if item == 90:
         var = ev.dPhiGammaMesonCand
      if item == 91:
         var = ev.dEtaGammaMesonCand
      if item == 95:
         var = ev.PV_npvsGood

# Fill histograms.

      if ev.mc > 0 : wei = ev.w * ev.lumiIntegrated
#      if(category=='_VBFcat' or category=='_Wcat' or category=='_Zcat'): wei = ev.w * ev.lumiIntegrated
#      if ev.mc==0:
      if (ev.mc>=37 and ev.mc<=44):
         hZinv.Fill( var, wei )
      if (ev.mc==34 or ev.mc==35 or ev.mc==36 or ev.mc==0):
         hDY.Fill( var, wei )
      if ev.mc==1 or ev.mc==46 or ev.mc==48: #Zll,Zhad,Znn
         hZG.Fill( var, wei )
      if ev.mc==2 or ev.mc==45: #Wln,Whad
         hWG.Fill( var, wei )
      if ev.mc==47:
         hTTG.Fill( var, wei )
      if (ev.mc==31 or ev.mc==32 or ev.mc==33):
         hW.Fill( var, wei)
      if (ev.mc==4 or ev.mc==5):
         hTT12L.Fill( var, wei)
      if (ev.mc==6 or ev.mc==7 or ev.mc==8 or ev.mc==9 or ev.mc==10 or ev.mc==11 or ev.mc==12 or ev.mc==13 or ev.mc==14):
         hGJet.Fill( var, wei)
      if (ev.mc==15):
         hVBFGJet.Fill( var, wei)
      if (ev.mc==20 or ev.mc==21 or ev.mc==22 or ev.mc==23 or ev.mc==24 or ev.mc==25):
         hJet.Fill( var, wei)
###
      if ev.mc==1013 or ev.mc==1023 or ev.mc==1014 or ev.mc==1024:
         hZH.Fill( var, wei)
      if ev.mc==1015 or ev.mc==1025 or ev.mc==1016 or ev.mc==1026:
         hZinvH.Fill( var, wei)
      if ev.mc==1011 or ev.mc==1012 or ev.mc==1021 or ev.mc==1022:
         hWH.Fill( var, wei)
      if ev.mc==1010 or ev.mc==1020:
         hVBFH.Fill( var, wei)
      if ev.mc==1017 or ev.mc==1027:
         hggH.Fill( var, wei)
      if ev.mc<0:
#      if (ev.mc>=37 and ev.mc<=44):
         if ((item==4 or item==43) and var>115 and var<135) : continue
         if item==42:
            if var>MVAbin[category]: continue
         hData.Fill( var )

   hZinvH.Scale(BTheo)
   hZH.Scale(BTheo)
   hVBFH.Scale(BTheo)
   hggH.Scale(BTheo)
   hWH.Scale(BTheo)

   redDark = (191, 34, 41)
#   redMed = (255, 40, 0)
   redMed = (237, 41, 57)
#   redMed = (220,0,5)
#   redLight = (255,61,65)
   redLight = (255,82,82)
   orange = (255, 204, 153)
   orangeDark = (255,150,79)
   gray = (136,139,141)
   azure = (100, 192, 232)
   azureDark = (96, 147, 172)
   green = (144, 238, 144)
   greenDark = (98, 172, 141)
   gold = (212 ,175, 55)
   violet = (181 ,100, 227)

   c, pad1, pad2 = createCanvasPads(doLog)

   # Draw stackXS with MC contributions
   stack = ROOT.THStack()
   BKGstack = ROOT.THStack()
   for h in [hZG, hDY, hTT12L, hTTG, hWG, hW, hGJet]:
      BKGstack.Add(h)

   for h, color in zip([hZG, hDY, hTT12L, hWG, hW, hTTG, hGJet, hVBFGJet, hJet, hZinv, hZH, hWH, hVBFH, hZinvH, hggH], [azure, azureDark, gray, green, greenDark, gray, orange, orangeDark, gold, violet, redMed, redDark, redLight, redMed, redDark]):
#   for h, color in zip([hZH, hWH, hVBFH, hZinvH, hggH], [redDark, redMed, redLight, redDark, redDark]):
      h.SetLineWidth(1)
      h.SetLineColor(1)
      h.SetFillColor(ROOT.TColor.GetColor(*color))
#      print('lumi',lumis[year])
#      if not category == '_VBFcat': h.Scale(lumis[year])
#      if not category == '_VBFcat': h.Scale(lumis[year])
      stack.Add(h)
   if item==4 or item==43:
      if BTheo==1. and doLog:
         stack.SetMaximum(10*max(hWH.GetMaximum(),hZH.GetMaximum(), hVBFH.GetMaximum(), hZinvH.GetMaximum(), hggH.GetMaximum()))
         if category == '_GFcat' and mesonCat == '_PhiCat': stack.SetMinimum(1)
         if category == '_VBFcat' and mesonCat == '_PhiCat': stack.SetMinimum(1)
         if category == '_VBFcatlow' and mesonCat == '_PhiCat': stack.SetMinimum(1)
      else:
         stack.SetMaximum(2*max(hData.GetMaximum(), hWH.GetMaximum(),hZH.GetMaximum(), hVBFH.GetMaximum(), hZinvH.GetMaximum(), hggH.GetMaximum()))

   if item==42 or item==14 or item==15 or item==22: stack.SetMaximum(2*max(hData.GetMaximum(), hWH.GetMaximum(),hZH.GetMaximum(), hVBFH.GetMaximum(), hZinvH.GetMaximum()))
   if item==102: stack.SetMaximum(10*max(hWG.GetMaximum(),hZG.GetMaximum(), hGJet.GetMaximum()))
   if item==1 or item==101 or item==51 or item==52 or item==23 or item==24 or item==12 or item==41: stack.SetMaximum(2*max(hData.GetMaximum(), hWH.GetMaximum(),hZH.GetMaximum(), hVBFH.GetMaximum(), hZinvH.GetMaximum(), hggH.GetMaximum()))
   if item==2 and category == '_GFcat': stack.SetMaximum(30000)     #Meson Mass
   if item==17 and mesonCat == '_RhoCat': stack.SetMaximum(4000)    # MJJ
   if item==17 and mesonCat == '_PhiCat': stack.SetMaximum(600)    # MJJ
   if item==2 and (category == '_VBFcat' or category == '_VBFcatlow') : stack.SetMaximum(2000)
   if item==9 and (category == '_VBFcat' or category == '_VBFcatlow' and category == '_GFcat') : stack.SetMaximum(2*max(hData.GetMaximum(),hGJet.GetMaximum()))
   if (item==36 or item==37) and category == '_Zinvcat': stack.SetMaximum(100)    # Zinv
   if item==6 and category == '_Zinvcat': stack.SetMaximum(300)    # Zinv
   if item==42 and category == '_GFcat' and mesonCat == '_PhiCat': stack.SetMaximum(8000)    # GF
   if item==1 and category == '_Zcat' and mesonCat == '_RhoCat': stack.SetMaximum(300)    # MZ

   print("ALL ev: sig(ggH) = ", hggH.Integral(), ' sig(hVBFH) = ', hVBFH.Integral())

   pad1.cd()
   # Draw data first
   if hData: hData.SetMarkerStyle(20)
   if hData: hData.SetMarkerSize(1.2)
   if hData: hData.SetLineWidth(2)
   if hData: hData.SetLineColor(ROOT.kBlack)
   if item==35: hData.GetXaxis().SetTitle("meson Isolation")
   if item==9: hData.GetXaxis().SetTitle("vertex compatibility")
   if item==21 and mesonCat == "_PhiCat": hData.GetXaxis().SetTitle("#sigmaM_{#phi}/M_{#phi}")
   if item==21 and mesonCat == "_RhoCat": hData.GetXaxis().SetTitle("#sigmaM_{#rho}/M_{#rho}")
   if item==41: hData.GetXaxis().SetTitle("p_{T}^{#phi,#gamma}")
   if item==42: hData.GetXaxis().SetTitle("MVA discriminator")
   #   if item==35 or item==9 or item==21:
   if item==35 or item==21:
      if hData: hData.GetYaxis().SetTitle("normalized events")
      if hData: hData.Scale(1./hData.Integral())
      if hData and item==35: hData.SetMaximum(100*hData.GetMaximum())
      if hData: hData.Draw("E")
   else:
      if hData: hData.Draw("E ")

#   print("ALL ev: sig(ggH) = ", hggH.Integral(), ' sig(hVBFH) = ', hVBFH.Integral())
   if hData: print("ALL Data (All except the SR) : ", hData.Integral())
#   std::cout << "INTEGRAL data" << hData.GetIntegral() << std::endl;

#   if (not doLog) and hData: hData.SetMinimum(0)
#   if hData: hData.SetMaximum(2*hData.GetMaximum())
#   if hData: hData.Draw("")
#   if not hData: stack.Draw("HIST")
#   if item==35 or item==9 or item==21:
   if item==35 or item==21:
      h = stack.GetStack().Last()
      h.DrawNormalized("HIST same")
   else:

      pad1.cd()
      stack.Draw("HIST")

      stack.GetXaxis().SetLabelSize(0.04)
      stack.GetXaxis().SetTitleSize(0.045)
      #stack.GetXaxis().SetTitleOffset(1.3)
      if item==1:
         if(category =='_Wcat'): stack.GetXaxis().SetTitle("mT_{1l,MET}^{W#rightarrow l#nu} [GeV]")
         if(category =='_Zcat'): stack.GetXaxis().SetTitle("m_{2l}^{Z#rightarrow ll} [GeV]")
      if item==2:
         if mesonCat == "_PhiCat": stack.GetXaxis().SetTitle("m_{2trk}^{ #phi #rightarrow k k } [GeV]")
         if mesonCat == "_RhoCat": stack.GetXaxis().SetTitle("m_{2trk}^{ #rho #rightarrow #pi #pi } [GeV]")
      if item==3:
         if mesonCat == "_PhiCat": stack.GetXaxis().SetTitle("pT_{2trk}^{ #phi #rightarrow k k } [GeV]")
         if mesonCat == "_RhoCat": stack.GetXaxis().SetTitle("pT_{2trk}^{ #rho #rightarrow #pi #pi } [GeV]")
      if item==16:
         if mesonCat == "_PhiCat": stack.GetXaxis().SetTitle("#DeltaR ( k k )")
         if mesonCat == "_RhoCat": stack.GetXaxis().SetTitle("#DeltaR ( #pi #pi )")
      if (item==4 or item==43) and mesonCat == "_PhiCat": stack.GetXaxis().SetTitle("m_{#gamma#phi}^{H} [GeV]")
      if (item==4 or item==43) and mesonCat == "_RhoCat": stack.GetXaxis().SetTitle("m_{#gamma#rho}^{H} [GeV]")

      if item==11: stack.GetXaxis().SetTitle("SoftActivityJetNjets5")
      if item==12: stack.GetXaxis().SetTitle("nJet")
      if item==5: stack.GetXaxis().SetTitle("P_{T}^{#gamma}")
      if item==6: stack.GetXaxis().SetTitle("#eta^{#gamma}")
      if item==206: stack.GetXaxis().SetTitle("H/E")
      if item==207: stack.GetXaxis().SetTitle("r9")
      if item==208: stack.GetXaxis().SetTitle("pfRelIso")
      if item==209: stack.GetXaxis().SetTitle("#sigma_{i#eta i#eta}")
      if item==22: stack.GetXaxis().SetTitle("|p_{T}^{miss}-p_{T}^{H_{#rho#gamma}}|/p_{T}^{H_{#rho#gamma}}")
      if item==14: stack.GetXaxis().SetTitle("#Delta #phi (p_{T}^{miss},#gamma)")
      if item==15: stack.GetXaxis().SetTitle("Delta #phi (p_{T}^{miss},meson)")
      if item==42: stack.GetXaxis().SetTitle("MVA discriminator")
      if item==13: stack.GetXaxis().SetTitle("E_{T}^{miss}")
      if item==17: stack.GetXaxis().SetTitle("M_{JJ}")
      if item==36: stack.GetXaxis().SetTitle("#phi E_{T}^{miss}")
      if item==37: stack.GetXaxis().SetTitle("#phi Higgs")

      if item==31: stack.GetXaxis().SetTitle("p_{T}^{leadTrk}")
      if item==32: stack.GetXaxis().SetTitle("p_{T}^{subLeadTrk}")

      stack.GetYaxis().SetTitle("Events")
      if item == 4: stack.GetYaxis().SetTitle("Events/ 1 [GeV]")
      stack.GetYaxis().SetTitleOffset(1.1)
      stack.GetYaxis().SetLabelSize(0.04)
      stack.GetYaxis().SetTitleSize(0.045)
      stack.GetYaxis().ChangeLabel(1, -1, 0)

      # Draw data again
#      if hData: hData.DrawNormalized("E SAME")
      if hData: hData.Draw("E SAME")

      pad2.cd()
      ratio = hData.Clone("dataratio")
      mcTOT = BKGstack.GetStack().Last()
      ratio.Divide(mcTOT)
      ratio.GetYaxis().SetTitle("data/MC")
#      ratio.GetYaxis().SetRangeUser(0.,2.)
      ratio.GetYaxis().SetRangeUser(0.5,1.5)
      ratio.GetXaxis().SetTitleOffset(4.)
      ratio.GetXaxis().SetTitleSize(0.15)
      ratio.GetXaxis().SetLabelSize(0.12)

      ratio.GetYaxis().SetTitleOffset(0.3)
      ratio.GetYaxis().SetTitleSize(0.15)
      ratio.GetYaxis().SetLabelSize(0.12)

      ratio.Draw("pe")
      lineZero = ROOT.TLine(mcTOT.GetXaxis().GetXmin(), 1.,  mcTOT.GetXaxis().GetXmax(), 1.)
      lineZero.SetLineColor(11)
      lineZero.Draw("same")

   # Add legend
   deltay=0
   if(category =='_Wcat' or category =='_Zcat'): deltay=0.25
   if(category =='_VBFcat' or category =='_VBFcatlow' or category =='_GFcat'): deltay=0.15
   if(category =='_VBFcat' or category =='_VBFcatlow' or category =='_GFcat'): deltay=0.25
   if(category =='_Zinvcat'): deltay=0.25
   if(category =='_Zinvcat'): lenght = 0.35
   deltax=0
   if(item==35 or item==4 or item==42): deltax=0.4
   lenght = 0.22
   if BTheo!=1.: lenght = 0.35
   if BTheo!=1.: deltax = 0.5
   legend = ROOT.TLegend(0.63-deltax, 0.88-deltay, 0.66-deltax+lenght, 0.88)
   legend.SetTextFont(42)
   legend.SetFillStyle(0)
   legend.SetBorderSize(0)
   legend.SetTextSize(0.04)
   legend.SetTextAlign(32)
#   if hData and hData.Integral()>0: legend.AddEntry(hData, "Data("+str(math.ceil(hData.Integral()))+")" ,"lep")
   if hData and hData.Integral()>0: legend.AddEntry(hData, "Data" ,"lep")

   if item!=43:
      print("hZinv=",hZinv.Integral())
      if hZinv and hZinv.Integral()>0: legend.AddEntry(hZinv, "Z #nu #nu", "f")
      if hZG and hZG.Integral()>0: legend.AddEntry(hZG, "ZG", "f")
      if hDY and hDY.Integral()>0: legend.AddEntry(hDY, "DY + jets", "f")
      if hWG and hWG.Integral()>0: legend.AddEntry(hWG, "WG", "f")
      if hW and hW.Integral()>0: legend.AddEntry(hW, "W + jets", "f")
      if hTTG and hTTG.Integral()>0: legend.AddEntry(hTTG, "TTG", "f")
      if hTT12L and hTT12L.Integral()>0: legend.AddEntry(hTT12L, "TT12l", "f")
      if hGJet and hGJet.Integral()>0: legend.AddEntry(hGJet, "#gamma + jets", "f")
      if hVBFGJet and hVBFGJet.Integral()>0: legend.AddEntry(hVBFGJet, "#gamma + jets (EWK)", "f")
      ######
   if hZinvH and hZinvH.Integral()>0 and mesonCat == "_PhiCat":
      if BTheo!=1.:
         legend.AddEntry(hZinvH, "Z#nu#nuH(#gamma#phi) SM x 10^{4}", "f")
      else:
         legend.AddEntry(hZinvH, "Z#nu#nuH(#gamma#phi)", "f")
   if hZH and hZH.Integral()>0 and mesonCat == "_PhiCat":
      if BTheo!=1.:
         legend.AddEntry(hZH, "ZH(#gamma#phi) SM x 10^{3}", "f")
      else:
         legend.AddEntry(hZH, "ZH(#gamma#phi)", "f")
   if hWH and hWH.Integral()>0 and mesonCat == "_PhiCat":
      if BTheo!=1.:
         legend.AddEntry(hWH, "WH(#gamma#phi) SM x 10^{4}", "f")
      else:
         legend.AddEntry(hWH, "WH(#gamma#phi)", "f")
   ######
   if hZinvH and hZinvH.Integral()>0 and mesonCat == "_RhoCat":
      if BTheo!=1.:
         legend.AddEntry(hZinvH, "Z#nu#nuH(#gamma#rho) SM x 10^{4}", "f")
      else:
         legend.AddEntry(hZinvH, "Z#nu#nuH(#gamma#rho)", "f")
   if hZH and hZH.Integral()>0 and mesonCat == "_RhoCat":
      if BTheo!=1.:
         legend.AddEntry(hZH, "ZH(#gamma#rho) SM x 10^{4}", "f")
      else:
         legend.AddEntry(hZH, "ZH(#gamma#rho)", "f")
   if hWH and hWH.Integral()>0 and mesonCat == "_RhoCat":
      if BTheo!=1.:
         legend.AddEntry(hWH, "WH(#gamma#rho) SM x 10^{4}", "f")
      else:
         legend.AddEntry(hWH, "WH(#gamma#rho)", "f")
   ######
   if hVBFH and mesonCat == "_PhiCat" and hVBFH.Integral()>0:
      if BTheo!=1.:
         if category=='_VBFcatlow': legend.AddEntry(hVBFH, "qqH(#gamma#phi) SM x 10^{4}", "f")
         if category=='_VBFcat': legend.AddEntry(hVBFH, "qqH(#gamma#phi) SM x 5 x 10^{3}", "f")
         if category=='_GFcat': legend.AddEntry(hVBFH, "qqH(#gamma#phi) SM x 10^{3}", "f")
      else:
         legend.AddEntry(hVBFH, "qqH(#gamma#phi)", "f")
   if hVBFH and mesonCat == "_RhoCat" and hVBFH.Integral()>0:
      if BTheo!=1.:
         if (category=='_VBFcatlow'): legend.AddEntry(hVBFH, "qqH(#gamma#rho) SM x 5 x 10^{3}", "f")
         if category=='_GFcat' or category=='_VBFcat': legend.AddEntry(hVBFH, "qqH(#gamma#rho) SM x 10^{3}", "f")
      else:
         legend.AddEntry(hVBFH, "qqH(#gamma#rho)", "f")
   #####
   if hggH and mesonCat == "_PhiCat" and hggH.Integral()>0:
      if BTheo!=1.:
         if (category=='_VBFcatlow'): legend.AddEntry(hggH, "ggH(#gamma#phi) SM x 10^{4}", "f")
         if category=='_GFcat' or category=='_VBFcat': legend.AddEntry(hggH, "ggH(#gamma#phi) SM x 10^{3}", "f")
      else:
         legend.AddEntry(hggH, "ggH(#gamma#phi)", "f")
   if hggH and mesonCat == "_RhoCat" and hggH.Integral()>0:
      if BTheo!=1.:
         if category=='_VBFcatlow': legend.AddEntry(hggH, "ggH(#gamma#rho) SM x 10^{4}", "f")
         if category=='_GFcat' or category=='_VBFcat': legend.AddEntry(hggH, "ggH(#gamma#rho) SM x 10^{3}", "f")
      else:
         legend.AddEntry(hggH, "ggH(#gamma#rho)", "f")

   pad1.cd()
   legend.Draw("SAME")

#   bmin = hZinvH.FindBin(110)
#   bmax = hZinvH.FindBin(140)
#   integralSig = hZinvH.Integral(hZinvH.FindBin(110),hZinvH.FindBin(140))
#   bmin = hZH.FindBin(110)
#   bmax = hZH.FindBin(140)
#   integralSig = hZH.Integral(hZH.FindBin(110),hZH.FindBin(140))
#   bmin = hWH.FindBin(110)
#   bmax = hWH.FindBin(140)
#   integralSig = hWH.Integral(hWH.FindBin(110),hWH.FindBin(140))
   bmin = hggH.FindBin(115)
   bmax = hggH.FindBin(135)
   integralSig = hggH.Integral(hggH.FindBin(115),hggH.FindBin(135))
   print("integralSig(115:135)=",integralSig)

   if item == 42:
      print('MVAbin=',MVAbin[category])
      if category=='_GFcat':
         integralSigCut = hggH.Integral(hggH.FindBin(MVAbin[category]),hggH.FindBin(2.))
         integralSig = hggH.Integral(hggH.FindBin(-1.),hggH.FindBin(2.))
         frac = round((integralSigCut/integralSig)*100,2)
         print("MVA integralSig(0.:1.)=",frac)
      if category=='_VBFcatlow' or category=='_VBFcat':
         integralSigCut = hVBFH.Integral(hggH.FindBin(MVAbin[category]),hVBFH.FindBin(2.))
         integralSig = hVBFH.Integral(hggH.FindBin(-1.),hVBFH.FindBin(2.))
         frac = round((integralSigCut/integralSig)*100,2)
         print("MVA integralSig(0.3:1.)=",frac)
      if category=='_Zinvcat':
         integralSigCut = hZinvH.Integral(hZinvH.FindBin(MVAbin[category]),hZinvH.FindBin(2.))
         integralSig = hZinvH.Integral(hZinvH.FindBin(-1.),hZinvH.FindBin(2.))
         frac = round((integralSigCut/integralSig)*100,2)
         print("MVA integralSig(0.6:1.)=",frac)
#      print("MVA integralSig(0.3:2.)=",integralSigCut)
#      print("MVA(denom) integralSig(-1.:2.)=",integralSig)

#   bminDllow = hData.FindBin(100)
#   bmaxDlow = hData.FindBin(115)
   integralDataLow = hData.Integral(hData.FindBin(100),hData.FindBin(115))
   integralDataHigh = hData.Integral(hData.FindBin(135),hData.FindBin(170))
   integralDataLowOUT = hData.Integral(hData.FindBin(0),hData.FindBin(100))
   integralDataHighOUT = hData.Integral(hData.FindBin(171),hData.FindBin(172))
   print("integralData(100:115)=",integralDataLow)
   print("integralData(135:170)=",integralDataHigh)

   print("integralData(0:100)=",integralDataLowOUT)
   print("integralData(170:--> Inf)=",integralDataHighOUT)

   # Add label
   text = ROOT.TLatex()
   text.SetNDC()
   text.SetTextFont(72)
   text.SetTextSize(0.045)
   text.DrawLatex(0.15, 0.93, "CMS")
   text.SetTextFont(42)
#   text.DrawLatex(0.15 + 0.10, 0.93, "Simulation")
   text.DrawLatex(0.15 + 0.10, 0.93, "Internal")
   text.SetTextSize(0.04)
   text.DrawLatex(0.65, 0.93, "#sqrt{s} = 13 TeV,%0.2f fb^{-1}"% (lumis[year]))

#   n = 4
#   x, y = array( 'd' ), array( 'd' )
#   x.append(100)
#   x.append(150)
#   x.append(150)
#   x.append(150)
#   y.append(0)
#   y.append(0)
#   y.append(2000)
#   y.append(2000)

   # Add TLine blind
   line1 = ROOT.TLine( 115, 0., 115., 500.)
   line1.SetLineColor(11);
   if item==4: line1.Draw()
   line2 = ROOT.TLine( 135, 0, 135., 500.)
   line2.SetLineColor(11);
   if item==4: line2.Draw()

#   stack.Draw("HIST")
#   hZinv.Draw()

#   gr = ROOT.TGraph( n, x, y )
#   gr.SetLineColor(5);
#   gr.SetFillStyle(3002);
#   gr.SetFillColor(2);

#   string = category+mesonCat+year+"_isMu"
#   string = category+mesonCat+year+"_isEle"
   string = category+mesonCat+year
#   if category=='_VBFcat' or category=='_VBFcatlow': string = category+mesonCat

   mydir = "PLOTS/"

   plotString = 'default'

   # Save the plot
   if item==41: plotString = "HCandPT"
   if item==42: plotString = "MVAdisc"
   if item==43: plotString = "HCandMassWithCuts"
   if item==1: plotString = "VMass"
   if item==4: plotString = "HCandMass"

   #####
   if item==2: plotString = "MesonCandidate_Mass"
   if item==3: plotString = "MesonCandidate_Pt"
   if item==35: plotString = "MesonCandidate_Iso"

   if item==7: plotString = "MesonSipPV"
   if item==8: plotString = "MesonDOCA"
   if item==9: plotString = "MesonVTXprob"

   if item==16: plotString = "DRtk12"
   if item==10: plotString = "DRtk12norm"
   if item==18: plotString = "chi2dof"
   if item==21: plotString = "massErr"

   if item==29: plotString = "MesonCandidate_Trk1Pt_norm"
   if item==30: plotString = "MesonCandidate_Trk2Pt_norm"
   if item==31: plotString = "MesonCandidate_Trk1Pt"
   if item==32: plotString = "MesonCandidate_Trk2Pt"
   if item==33: plotString = "MesonCandidate_Trk1Eta"
   if item==34: plotString = "MesonCandidate_Trk2Eta"

   if item==5: plotString = "PhotonPt"
   if item==6: plotString = "Photoneta"
   if item==205: plotString = "NPhotons"
   if item==206: plotString = "Photon_hoe"
   if item==207: plotString = "Photon_r9"
   if item==208: plotString = "Photon_relIso"
   if item==209: plotString = "Photons_sieie"
   if item==101: plotString = "Photons_pixelSeed"


   if item==11: plotString = "SoftActivityJetNjets5"
   if item==12: plotString = "nJet"
   if item==51: plotString = "drJet"
   if item==52: plotString = "drPHO"

   ### VBF - like
   if item==17: plotString = "Mjj"
   if item==19: plotString = "Detajj"
   if item==20: plotString = "Y1Y2jj"
   if item==23: plotString = "ptJ1"
   if item==24: plotString = "ptJ2"
   if item==231: plotString = "etaJ1"
   if item==241: plotString = "etaJ2"
   if item==53: plotString = "zepVar"
   if item==54: plotString = "dEtaHIGjet1"
   if item==55: plotString = "dEtaHIGjet2"
   if item==56: plotString = "jet1hfsigmaPhiPhi"
   if item==57: plotString = "jet1hfsigmaEtaEta"
   if item==58: plotString = "jet2hfsigmaPhiPhi"
   if item==59: plotString = "jet2hfsigmaEtaEta"

   ## Zinv like
   if item==13: plotString = "METpt"
   if item==14: plotString = "dPhiMETphoton"
   if item==15: plotString = "dPhiMETmeson"
   if item==36: plotString = "METphi"
   if item==37: plotString = "HiggsPHI"
   if item==22: plotString = "ptRatio"

   ##
   if item==90: plotString = "dPhiPHMes"
   if item==91: plotString = "dEtaPHMes"
   if item==95: plotString = "PV_npvsGood"
   if item==113: plotString = "nBjet"

   if item==114: plotString = "Visr_mass"
   if item==115: plotString = "LeadingLepton"
   if item==116: plotString = "SubLeadingLepton"

#   if not isHEM:
#      c.SaveAs(mydir+plotString+string+".png")
#   else:
   c.SaveAs(mydir+plotString+string+".png")
#   c.SaveAs(mydir+plotString+string+".png")
   print(plotString+".png")

   ##
   if item==102:
      c.SaveAs(mydir+"KSCandidate_Mass"+string+".png")
      print("KSMass.png")
   if item==103:
      c.SaveAs(mydir+"KSCandidate_Pt"+string+".png")
      print("KSPt.png")
   if item==107:
      c.SaveAs(mydir+"KSCandidate_Lxy"+string+".png")
      print("KSCandidate_Lxy.png")

def plotPhoton():

#   plot(101, 10, 0. , 10., False) # Photon pixel seed
#   plot(206, 100, 0. , 0.1, True) # ph h/e
#   plot(207, 60, 0.5 , 1.1, False) # ph
#   plot(208, 100, 0. , 0.1, True) # ph relIso all
#   plot(209, 100, 0. , 0.1, False) # ph sieta sieta

   plot(5, 30, 0. , 150., False) # Photon_pt
   plot(6, 20, -5. , 5., False) # Photon_eta

def plotMeson():

   if mesonCat == '_PhiCat': plot(2, 40, 1. , 1.04, False) # phi_mass
   if mesonCat == '_RhoCat': plot(2, 50, 0.5 , 1.00, False) # rho_mass

   plot(3, 20, 0. , 100., True) # meson_pt

   plot(31, 80, 0. , 80., False) # max pt_trk1
   plot(32, 80, 0. , 80., False) # min pt_trk2

#   plot(29, 100, 0. , 1., True) # pt_trk1_norm min
#   plot(30, 100, 0. , 1., True) # pt_trk2_norm max

#   plot(33, 30, -3. , 3., False) # meson_trk1_eta
#   plot(34, 30, -3. , 3., False) # meson_trk2_eta

#   if mesonCat == '_PhiCat': plot(10, 100, 0. , 1., True) # DRtk12 norm
#   if mesonCat == '_RhoCat': plot(10, 100, 0. , 5., True) # DRtk12 norm

   if mesonCat == '_PhiCat': plot(16, 50, 0. , 0.05, False) # DRtk12
   if mesonCat == '_RhoCat': plot(16, 50, 0. , 0.1, False) # DRtk12

   plot(35, 120, 0. , 1.2, False) # meson_iso

   plot(7, 100, 0. , 10., True) # phi_lxy
   plot(9, 30, 0. , 1.5, False) # phi_VTXprob
#   plot(8, 100, 0. , 0.1, True) # phi_doca
#   plot(18, 100, 0. , 1., False) # vtx_chi2dof

#   if(category =='_Wcat' or category =='_Zcat'):
#      plot(50, 100, 0. , 10., False) # minDR (Lepton,Meson)

#   plot(21, 100, 0. , 0.1, True) # massErr meson

def plotWZlep():

   if(category =='_Wcat' or category =='_Zcat'):

      if(category =='_Zcat'): plot(1, 30, 75. , 105., False) # Z_mass
      if(category =='_Wcat'): plot(1, 21, 15. , 120., False) # W_mt

      if(category =='_Zcat'): plot(114, 30, 75. , 105., False) # Z_mass (llgamma)
      if(category =='_Zcat' or category =='_Wcat'): plot(115, 65, 0. , 65., False) # leadingLep
      if(category =='_Zcat'): plot(116, 65, 0. , 65., False)   #subLeadinLep

      if(category =='_Wcat'): plot(14, 20, 0. , 3.14, False) # dPhiMETmeson
      if(category =='_Wcat'): plot(15, 20, 0. , 3.14, False) # dPhiMETphoton

   if(category =='_Wcat'): plot(13, 40, 0. , 200., False) # MET_pt

def plotZinv():

   if (category=='_Zinvcat'):
      plot(13, 40, 0. , 200., False) # MET_pt

      plot(22, 50, 0. , 2., False) # PTratio
      plot(14, 20, 0. , 3.14, False) # dPhiMETmeson
      plot(15, 20, 0. , 3.14, False) # dPhiMETphoton

      plot(91, 100, 0. , 10, False) # Deta Meson-Gamma
      plot(90, 100, 0. , 6.28, False) # Dphi Meson-Gamma

      plot(37, 50, -6.28 , 6.28, False) # Higgs phi
      plot(36, 100, -6.28 , 6.28, False) # Met phi

def plotVBF():

   if (category=='_VBFcatlow' or category=='_VBFcat'):
      plot(17, 60, 0. , 3000., False) # MJJ_mass

      plot(51, 100, 0. , 10., False) # minDR (jet,Meson)
      plot(52, 100, 0. , 10., False) # minDR (jet,Photon)

      plot(23, 20, 0. , 200., False) # jet1 Pt
      plot(24, 20, 0. , 200., False) # jet2 Pt

      plot(231, 20, -5. , 5., False) # jet1 Eta
      plot(241, 20, -5. , 5., False) # jet2 Eta

      plot(19, 100, 0. , 10, False) # Deta
      plot(20, 50, -5 , 5, False) # Y1Y2

      plot(53, 50, 0. , 5., False) # zepVar
      plot(54, 100, 0. , 10., False) # dEta HIG - jet1
      plot(55, 100, 0. , 10., False) # dEta HIG - jet2
      plot(56, 100, 0. , 0.5, False) # dPhiPhi jet1
      plot(57, 100, 0. , 0.5, False) # dEtaEta jet2
      plot(58, 100, 0. , 0.5, False) # dPhiPhi jet1
      plot(59, 100, 0. , 0.5, False) # dEtaEta jet2

if __name__ == "__main__":

#   plot(4, 171, 0. , 171., True) # HCandMass
#   if (category=='_VBFcatlow' or category=='_VBFcat' or category=='_GFcat' or category=='_Zinvcat'):
#      plot(43, 70, 100. , 170., False) # HCandMassWithCut
#      plot(43, 200, 0. , 200., False) # HCandMassWithCut
#      plot(42, 40, -1. , 1., False) # MVAdiscr

   plotPhoton()
   plotMeson()
   plotWZlep()
   plotZinv()
   plotVBF()

   exit()

   plot(41, 150, 0. , 150.,False) # HCandPT
   plot(22, 50, 0. , 2., False) # PTratio
   plot(12, 15, 0. , 15., False) # nJet
   plot(113, 15, 0. , 15., False) # nbjet
   plot(95, 100, 0. , 100., False) # NVTX
   plot(13, 40, 0. , 200., False) # MET_pt
   plot(11, 15, 0. , 15., False) # SoftActivityJetNjets5

if __name__ == "__plots__":

#    plot(102, 100, 0.45, 0.55, False) # ks_mass
#    plot(103, 80, 0. , 80., True) # ks_pt
#    plot(107, 100, 0. , 10., True) # ks_lxy
