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

ROOT.ROOT.EnableImplicitMT()
RDataFrame = ROOT.RDataFrame
    
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
 
   ####

   df = RDataFrame(mytree)

   inMass = False
   #   inMass = (abs(ev.HCandMass - 125)<30)

   if(category =='_Zcat'):
       offline = "photon_pt>40 and meson_pt>40"

   if(category =='_Wcat'):
       offline = "photon_pt>40 and meson_pt>40 and DeepMETResolutionTune_pt>15 and V_mass>15 and deltaLepMeson>0.5"
#       offline = "photon_pt>40 and meson_pt>40 and DeepMETResolutionTune_pt>15 and deltaLepMeson>0.5"
# VMass can be removed to gain stat

   if(category =='_GFcat'):
       offline = "photon_pt>40 and meson_pt>40"

   if(category =='_VBFcat'):
       offline = "jet1Pt>30 and jet2Pt>20 and mJJ>400"

   if(category =='_VBFcatlow'):
       offline = "jet1Pt>30 and jet2Pt>20 and mJJ>300"

   if(category =='_Zinvcat'):
       offline = "photon_pt>40 and meson_pt>40"

#      if(category =='_Zcat' or category =='_Wcat'):
#         if ev.isMuorEle!=1: continue; # 1 is diMu and 2 is diEle
#         if ev.isMuorEle!=2: continue; # 1 is diMu and 2 is diEle

#      ## OPTIMIZED PHASE SPACE
#      if(category =='_Zinvcat'):
#         if ev.photon_pt<40 :  continue
#         if ev.meson_pt<40 : continue
##         if min(ev.goodMeson_trk1_pt[ev.index_pair[0]],ev.goodMeson_trk2_pt[ev.index_pair[0]]) < 10: continue
#         if item !=13:
#            if ev.DeepMETResolutionTune_pt<75 :  continue
##            if ev.DeepMETResolutionTune_pt<50 :  continue
##         if abs(ev.dPhiGammaMET)<1. : continue # already applied
##         if abs(ev.dPhiMesonMET)<1. : continue # already applied
#         if item !=22:
#            if ev.ptRatioMEThiggs>0.5: continue
#         if item !=91:
#            if abs(ev.dEtaGammaMesonCand)>2: continue
#         if item !=113:
#            if ev.nbtag>0: continue

#         isHEM = ev.DeepMETResolutionTune_phi < -0.9 and ev.DeepMETResolutionTune_phi > -1.6
#         if not isHEM: continue

# HEM
#         isHEM = (ev.goodPhotons_eta[idxPh]< -1.39 and ev.goodPhotons_phi[idxPh]>= -1.6 and ev.goodPhotons_phi[idxPh]< -0.9)
#         if not isHEM: continue
#         if ev.run > 319077: continue

   if item == 1 : var = "V_mass"
   if item == 114 : var = "Visr_mass"
   if item == 115 : var = "LeadingLeptonPt"
   if item == 116 : var = "SubLeadingLeptonPt"
   if item == 42 : var = "MVAdisc[0]"
   if item == 43:
       var = -1
       if ev.MVAdisc[0]>MVAbin[category]:  var = ev.HCandMass

   if item == 4 : var = "HCandMass" #"HGenMass"
   if item == 41 : var = "HCandPT"
   ## Photon plots
   if item == 5 : var = "goodPhotons_pt[index_pair[1]]"
   if item == 6 : var = "goodPhotons_eta[index_pair[1]]"
   if item == 205 : var = "nPhotonsVeto"
   if item == 206 : var = "goodPhotons_hoe[index_pair[1]]"
   if item == 207 : var = "goodPhotons_r9[index_pair[1]]"
   if item == 208 : var = "goodPhotons_pfRelIso03_all[index_pair[1]]"
   if item == 209 : var = "goodPhotons_sieie[index_pair[1]]"
   if item == 210 : var = "goodPhotons_mvaID[index_pair[1]]"
   if item == 101 : var = "goodPhotons_pixelSeed[index_pair[1]]"
   ## KS plots
   if item == 102 : var = "goodKs_mass[index_pair[0]]"
   if item == 103 : var = "goodKs_pt[index_pair[0]]"
   if item == 107 : var = "goodKs_lxy[index_pair[0]] "
   ## Meson plots
   if item == 2 : var = "goodMeson_mass[index_pair[0]]"
   if item == 3 : var = "goodMeson_pt[index_pair[0]]"
   #leading
   if item == 30 : var = "max(goodMeson_trk1_pt[index_pair[0]],goodMeson_trk2_pt[index_pair[0]])/goodMeson_pt[index_pair[0]]"
   #subleading
   if item == 29 : var = "min(goodMeson_trk1_pt[index_pair[0]],goodMeson_trk2_pt[index_pair[0]])/goodMeson_pt[index_pair[0]]"
   #subleading
   if item == 39 : var = "min(goodMeson_trk1_pt[index_pair[0]],goodMeson_trk2_pt[index_pair[0]])/max(goodMeson_trk1_pt[index_pair[0]],goodMeson_trk2_pt[index_pair[0]])"
   #leading
   if item == 31 : var = "max(goodMeson_trk1_pt[index_pair[0]],goodMeson_trk2_pt[index_pair[0]])"
   #subleading
   if item == 32 : var = "min(goodMeson_trk1_pt[index_pair[0]],goodMeson_trk2_pt[index_pair[0]])"
   if item == 33 : var = "var = goodMeson_trk1_eta[index_pair[0]]"
   if item == 34 : var = "goodMeson_trk2_eta[index_pair[0]]"
   if item == 35 : var = "goodMeson_iso[index_pair[0]]"
   if item == 16 : var = "abs(goodMeson_DR[index_pair[0]])"
   if item == 10 : var = "(abs(goodMeson_DR[index_pair[0]])*goodMeson_pt[index_pair[0]])/(2*goodMeson_mass[index_pair[0]])"
   if item == 7 : var = "goodMeson_sipPV[index_pair[0]]"
   if item == 8 : var = "goodMeson_doca_prob[index_pair[0]]"
   if item == 9 : var = "goodMeson_vtx_prob[index_pair[0]]/goodMeson_vtx_chi2dof[index_pair[0]]"
   if item == 18 : var = "goodMeson_vtx_chi2dof[index_pair[0]]"
   if item == 21 : var = "goodMeson_massErr[index_pair[0]]"
   #         var = ev.goodMeson_massErr[idxMeson]/ev.goodMeson_mass[idxMeson]

   ## event like
   if item == 11 : var = "SoftActivityJetNjets5"
   if item == 12 :
       if(category =='_VBFcat' or category =='_VBFcatlow' or category =='_GFcat' or category =='_Zinvcat'): var = "nGoodJets"
       else: var = "nJet"
   if item == 113 :
         if(category =='_Zinvcat'): var = "nbtag"
   if item == 50 : var = "deltaLepMeson"
   if item == 51 : var = "deltaJetMeson"
   if item == 52 : var = "deltaJetPhoton"
   if item == 53 : var = "zepVar"
   if item == 54 : var = "detaHigJet1"
   if item == 55 : var = "detaHigJet2"
   if item == 56 : var = "jet1hfsigmaPhiPhi"
   if item == 57 : var = "jet1hfsigmaEtaEta"
   if item == 58 : var = "jet2hfsigmaPhiPhi"
   if item == 59 : var = "jet2hfsigmaEtaEta"
   if item == 13 : var = "DeepMETResolutionTune_pt"
   #         var = ev.MET_pt
   if item == 36 : var = "DeepMETResolutionTune_phi"
   #         var = ev.MET_pt
   if item == 37 : var = "HCandPHI"
   #         var = ev.MET_pt
   if item == 14 : var = "abs(dPhiGammaMET)"
   if item == 15 : var = "abs(dPhiMesonMET)"
   if item == 17 : var = "mJJ"
   if item == 19 : var = "dEtaJJ"
   if item == 20 : var = "Y1Y2"
   if item == 22: var = "ptRatioMEThiggs"
   if item == 23: var = "jet1Pt"
   if item == 24: var = "jet2Pt"
   if item == 231: var = "jet1Eta"
   if item == 241: var = "jet2Eta"
   if item == 90: var = "dPhiGammaMesonCand"
   if item == 91: var = "dEtaGammaMesonCand"
   if item == 95: var = "PV_npvsGood"

   df_common = df.Define("var","{}".format(var)).Define("offlineSel","{}".format(offline)).Define("weight","mc>=0 ? w_allSF*lumiIntegrated:1").Filter("offlineSel")

   hZinv = df_common.Filter("mc>=37 and mc<=44").Histo1D(("hZinv","h",nbin, low, high),"var","weight")
   hDY = df_common.Filter("mc==34 or mc==35 or mc==36 or mc==0").Histo1D(("hDY","h",nbin, low, high),"var","weight")
   hZG = df_common.Filter("mc==1 or mc==46 or mc==48").Histo1D(("hZG","h",nbin, low, high),"var","weight")
   hW = df_common.Filter("mc==31 or mc==32 or mc==33 or mc==3").Histo1D(("hW","h",nbin, low, high),"var","weight")
   hWG = df_common.Filter("mc==2 or mc==45").Histo1D(("hWG","h",nbin, low, high),"var","weight")

   hTTG = df_common.Filter("mc==47 or mc==49 or mc==51").Histo1D(("hTTG","h",nbin, low, high),"var","weight")
   hTT12L = df_common.Filter("mc==4 or mc==5").Histo1D(("hTT12L","h",nbin, low, high),"var","weight")
   hGJet = df_common.Filter("(mc==15 or mc==16 or mc==17 or mc==18 or mc==19) or (mc==10 or mc==11 or mc==12 or mc==13 or mc==14)").Histo1D(("hGJet","h",nbin, low, high),"var","weight")
   hVBFGJet = df_common.Filter("mc==9").Histo1D(("hVBFGJet","h",nbin, low, high),"var","weight")
   hJet = df_common.Filter("mc==20 or mc==21 or mc==22 or mc==23 or mc==24 or mc==25").Histo1D(("hJet","h",nbin, low, high),"var","weight")

   hZH = df_common.Filter("mc==1013 or mc==1023 or mc==1014 or mc==1024").Histo1D(("hZH","h",nbin, low, high),"var","weight")
   hZinvH = df_common.Filter("mc==1015 or mc==1025 or mc==1016 or mc==1026").Histo1D(("hZinvH","h",nbin, low, high),"var","weight")
   hWH = df_common.Filter("mc==1011 or mc==1012 or mc==1021 or mc==1022").Histo1D(("hWH","h",nbin, low, high),"var","weight")
   hTTH = df_common.Filter("mc==1018 or mc==1028").Histo1D(("hTTH","h",nbin, low, high),"var","weight")
   hVBFH = df_common.Filter("mc==1010 or mc==1020").Histo1D(("hVBFH","h",nbin, low, high),"var","weight")
   hggH = df_common.Filter("mc==1017 or mc==1027").Histo1D(("hggH ","h",nbin, low, high),"var","weight")
   if False: hData = df_common.Filter("mc<0").Histo1D(("hData","h",nbin, low, high),"var","weight")
   else:
       print('HELLO')
       if (item==4 or item==43): hData = df_common.Filter("mc<0 and (var<115 or var>135)").Histo1D(("hData","h",nbin, low, high),"var","weight")
       elif (item==42): hData = df_common.Filter("mc<0 and var<{}".format(MVAbin[category])).Histo1D(("hData","h",nbin, low, high),"var","weight")
       else: hData = df_common.Filter("mc<0").Histo1D(("hData","h",nbin, low, high),"var","weight")

   if hData: hData.SetMarkerStyle(20)
   if hData: hData.SetMarkerSize(1.2)
   if hData: hData.SetLineWidth(2)      

#   for h, color in zip([hZH, hWH, hVBFH, hZinvH, hggH], [redDark, redMed, redLight, redDark, redDark]):
   for h, color in zip([hZG, hDY, hTT12L, hWG, hW, hTTG, hGJet, hVBFGJet, hJet, hZinv, hZH, hWH, hTTH, hVBFH, hZinvH, hggH], [azure, azureDark, gray, green, greenDark, gray, orange, orangeDark, gold, violet, redMed, redDark, redLight, redLight, redMed, redDark]):
       if h:
           h.SetLineWidth(3)
           h.SetLineColor(ROOT.TColor.GetColor(*color))
           h.SetFillColor(ROOT.TColor.GetColor(*color))

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

   return listHisto
