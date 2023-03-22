import ROOT
import os
from array import array
import math
import sys

from LoadTree import loadTree

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()


class MyHisto():

    def __init__(self, name, hOBJ):
        self.name = name
        self.hOBJ = hOBJ

    def __repr__(self):
        return self.name


MVAbinRho={
#   '_GFcat': 0.0,                                                                                                                                                                                            
   '_GFcat': 0.5,
   '_Zinvcat': 0.6,
   '_VBFcat': 0.8,
   '_VBFcatlow': 0.8,
}

MVAbinPhi={
#   '_GFcat': 0.0,                                                                                                                                                                                            
   '_GFcat': 0.4,
   '_Zinvcat': 0.7,
   '_VBFcat': 0.6,
   '_VBFcatlow': 0.6,
}

    
def getHisto(mytree, category, mesonCat, item, nbin, low, high):

   if mesonCat=='_PhiCat': MVAbin = MVAbinPhi
   if mesonCat=='_RhoCat': MVAbin = MVAbinRho    

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
   hTTH = ROOT.TH1F( 'hTTH', 'This is the distribution', nbin, low, high)   
   hData = ROOT.TH1F( 'hData', '', nbin, low, high)

   for h, color in zip([hZG, hDY, hTT12L, hWG, hW, hTTG, hGJet, hVBFGJet, hJet, hZinv, hZH, hWH, hTTH, hVBFH, hZinvH, hggH], [azure, azureDark, gray, green, greenDark, gray, orange, orangeDark, gold, violet, redMed, redDark, redLight, redLight, redMed, redDark]):
#   for h, color in zip([hZH, hWH, hVBFH, hZinvH, hggH], [redDark, redMed, redLight, redDark, redDark]):
      h.SetLineWidth(1)
      h.SetLineColor(1)
      h.SetFillColor(ROOT.TColor.GetColor(*color))

   if hData: hData.SetMarkerStyle(20)
   if hData: hData.SetMarkerSize(1.2)
   if hData: hData.SetLineWidth(2)      

   hZinv_ = MyHisto('hZinv', hZinv)
   hDY_ = MyHisto('hDY', hDY)
   hZG_ = MyHisto('hZG', hZG)
   hW_ = MyHisto('hW', hW)
   hWG_ = MyHisto('hWG', hWG)
   hTTG_ = MyHisto('hTTG', hTTG)
   hTT12L_ = MyHisto('hTT12L',hTT12L)
   hGJet_ = MyHisto('hGJet', hGJet)
   hVBFGJet_ = MyHisto('hVBFGJet', hVBFGJet)
   hJet_ = MyHisto('hJet', hJet)
   hVBFH_ = MyHisto('hVBFH', hVBFH)
   hggH_ = MyHisto('hggH', hggH)
   hZH_ = MyHisto('hZH', hZH)
   hZinvH_ = MyHisto('hZinvH', hZinvH)   
   hWH_ = MyHisto('hWH', hWH)
   hTTH_ = MyHisto('hTTH', hTTH)   
   hData_ = MyHisto('hData', hData)

   
   listHisto = [hZG_, hDY_, hWG_, hW_, hTT12L_, hTTG_, hVBFGJet_, hGJet_, hJet_, hZinv_, hZH_, hWH_, hTTH_, hVBFH_, hZinvH_, hggH_, hData_]
   
   
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

      ## added already in MAY11
      if(category =='_Wcat'):
         if ev.DeepMETResolutionTune_pt<15 :  continue
         if ev.V_mass<15 :  continue
         if ev.deltaLepMeson<0.5 :  continue
         
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
#         if ev.dEtaJJ < 3. : continue # already applied

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
#        temp = HCandMass
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
         var = ev.nPhotonsVeto
      if item == 206 :
         var = ev.goodPhotons_hoe[idxPh]
      if item == 207 :
         var = ev.goodPhotons_r9[idxPh]
      if item == 208 :
         var = ev.goodPhotons_pfRelIso03_all[idxPh]
         #var = ev.goodPhotons_pfRelIso03_chg[idxPh]
      if item == 209 :
         var = ev.goodPhotons_sieie[idxPh]
      if item == 210 :
         var = ev.goodPhotons_mvaID[idxPh]
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
      if item == 39 :
         #subleading
         var = min(ev.goodMeson_trk1_pt[idxMeson],ev.goodMeson_trk2_pt[idxMeson])/max(ev.goodMeson_trk1_pt[idxMeson],ev.goodMeson_trk2_pt[idxMeson])
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

#      if ev.mc >= 0 : wei = ev.w * ev.lumiIntegrated
      if ev.mc >= 0 : wei = ev.w_allSF * ev.lumiIntegrated      
#      if(category=='_VBFcat' or category=='_Wcat' or category=='_Zcat'): wei = ev.w * ev.lumiIntegrated
#      if ev.mc==0:
      if (ev.mc>=37 and ev.mc<=44):
         hZinv.Fill( var, wei )
      if (ev.mc==34 or ev.mc==35 or ev.mc==36 or ev.mc==0):
#      if (ev.mc==34 or ev.mc==35 or ev.mc==36):
         hDY.Fill( var, wei )
      if ev.mc==1 or ev.mc==46 or ev.mc==48: #Zll,Zhad,Znn + Photon
         hZG.Fill( var, wei )
      if ev.mc==2 or ev.mc==45: #Wln,Whad + Photon
         hWG.Fill( var, wei )
      if (ev.mc==31 or ev.mc==32 or ev.mc==33 or ev.mc==3): # W fakes
         hW.Fill( var, wei)
###
      if ev.mc==47 or ev.mc==49 or ev.mc==51:
         hTTG.Fill( var, wei )
      if (ev.mc==4 or ev.mc==5):
         hTT12L.Fill( var, wei)
      if ((ev.mc==15 or ev.mc==16 or ev.mc==17 or ev.mc==18 or ev.mc==19) or (ev.mc==10 or ev.mc==11 or ev.mc==12 or ev.mc==13 or ev.mc==14)):
         hGJet.Fill( var, wei)
      if ev.mc==9: #VBF-EWK
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
      if ev.mc==1018 or ev.mc==1028:
         hTTH.Fill( var, wei)
      if ev.mc==1010 or ev.mc==1020:
         hVBFH.Fill( var, wei)
      if ev.mc==1017 or ev.mc==1027:
         hggH.Fill( var, wei)
      if ev.mc<0:
#      if (ev.mc>=37 and ev.mc<=44):
         if False: # temporary for the CR
#         if 'CR' in dirLOCAL_:
            hData.Fill( var )
         else:
            if ((item==4 or item==43) and var>115 and var<135) : continue
            if item==42:
               if var>MVAbin[category]: continue
            hData.Fill( var )

   return listHisto
