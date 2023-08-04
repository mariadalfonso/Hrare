import ROOT
import os
from array import array
import math
import sys

from LoadTree import loadTree
from prepareHisto import getHisto, MVAbinRho, MVAbinPhi

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

category = sys.argv[1]
mesonCat = sys.argv[2]
#year = '_12016'
#year = '_22016'
#year = '_2017'
year = '_2018'
#year = '_all'

if (category=='_GFcat' or category=='_VBFcatlow' or category=='_VBFcat'):
   dirLOCAL_= '/work/submit/mariadlf/JUL22_1Thread/'
#   dirLOCAL_= '/home/submit/mariadlf/Hrare/CMSSW_10_6_27_new/src/Hrare/analysis/JUL22_1Thread/'
else:
   dirLOCAL_= '/work/submit/mariadlf/JUL22/'

if (category=='_Wcat' and mesonCat == '_K0StarCat'):
   dirLOCAL_= '/home/submit/mariadlf/Hrare/CMSSW_10_6_27_new/src/Hrare/analysis/JUL27/'

if (category=='_Zinvcat'): directory4 = '/work/submit/mariadlf/Hrare/DEC28/2018/'
if (category=='_GFcat'): directory4 = dirLOCAL_+'2018/'
if (category=='_VBFcatlow'): directory4 = dirLOCAL_+'2018/'
if (category=='_VBFcat'):
   directory4 = dirLOCAL_+'2018/'
   directory3 = dirLOCAL_+'2017/'
   directory1 = dirLOCAL_+'12016/'
   year = '_all'
if (category=='_Zcat' or category=='_Wcat'):
   directory4 = dirLOCAL_+'2018/'
   directory3 = dirLOCAL_+'2017/'
   directory2 = dirLOCAL_+'22016/'
   directory1 = dirLOCAL_+'12016/'
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
if 'CR' in dirLOCAL_: BTheo=0.1
#BTheo=0.1

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

if mesonCat=='_PhiCat': MVAbin = MVAbinPhi
if mesonCat=='_RhoCat': MVAbin = MVAbinRho
if mesonCat=='_K0StarCat': MVAbin = MVAbinRho

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

   if BTheo==0.1: doLog=False # this is both for the meson mass SB and the signal region

   listHisto = getHisto(mytree, category, mesonCat, item, nbin, low, high)

   print(listHisto)
   for obj in listHisto:
      if obj.name == 'hData': hData = obj.hOBJ
      #
      if obj.name == 'hZinvH': hZinvH = obj.hOBJ
      if obj.name == 'hZH': hZH = obj.hOBJ
      if obj.name == 'hWH': hWH = obj.hOBJ
      if obj.name == 'hTTH': hTTH = obj.hOBJ
      if obj.name == 'hVBFH': hVBFH = obj.hOBJ
      if obj.name == 'hggH': hggH = obj.hOBJ
      #
      if obj.name == 'hZG': hZG = obj.hOBJ
      if obj.name == 'hDY': hDY = obj.hOBJ
      if obj.name == 'hWG': hWG = obj.hOBJ
      if obj.name == 'hW': hW = obj.hOBJ
      if obj.name == 'hTT12L': hTT12L = obj.hOBJ
      if obj.name == 'hTTG': hTTG = obj.hOBJ
      if obj.name == 'hGJet': hGJet = obj.hOBJ
      if obj.name == 'hVBFGJet': hVBFGJet = obj.hOBJ
      if obj.name == 'hJet': hJet = obj.hOBJ
      if obj.name == 'hZinv': hZinv = obj.hOBJ


   hZinvH.Scale(BTheo)
   hZH.Scale(BTheo)
   hTTH.Scale(BTheo)
   hVBFH.Scale(BTheo)
   hggH.Scale(BTheo)
   hWH.Scale(BTheo)

   c, pad1, pad2 = createCanvasPads(doLog)

   # Draw stackXS with MC contributions
   allstack = ROOT.THStack()
   BKGstack = ROOT.THStack()
   SIGstack = ROOT.THStack()
   for h in [hZG, hDY, hTT12L, hTTG, hWG, hW, hGJet]:
      BKGstack.Add(h.GetValue())

   for h in [hZH, hWH, hTTH, hVBFH, hZinvH, hggH]:
      if item!=4 and item!=43 and item!=42: h.SetFillStyle(0)
      SIGstack.Add(h.GetValue())

   for h in [hZG, hDY, hTT12L, hWG, hW, hTTG, hGJet, hVBFGJet, hJet, hZinv, hZH, hWH, hTTH, hVBFH, hZinvH, hggH]:
#   for h in [hVBFGJet, hJet, hZinv, hZH, hWH, hTTH, hVBFH, hZinvH, hggH]:
#      print('lumi',lumis[year])
#      if not category == '_VBFcat': h.Scale(lumis[year])
#      if not category == '_VBFcat': h.Scale(lumis[year])
      allstack.Add(h.GetValue())

   if item==4 or item==43 or item==42: stack = allstack
   else : stack = BKGstack

   if item==4 or item==43:
      if BTheo==1. and doLog:
         stack.SetMaximum(10*max(hWH.GetValue().GetMaximum(),hZH.GetValue().GetMaximum(), hTTH.GetValue().GetMaximum(), hVBFH.GetValue().GetMaximum(), hZinvH.GetValue().GetMaximum(), hggH.GetValue().GetMaximum()))
         if category == '_GFcat' and mesonCat == '_PhiCat': stack.SetMinimum(1)
         if category == '_VBFcat' and mesonCat == '_PhiCat': stack.SetMinimum(1)
         if category == '_VBFcatlow' and mesonCat == '_PhiCat': stack.SetMinimum(1)
      else:
         stack.SetMaximum(2*max(hData.GetValue().GetMaximum(), hWH.GetValue().GetMaximum(),hZH.GetValue().GetMaximum(), hTTH.GetValue().GetMaximum(), hVBFH.GetValue().GetMaximum(), hZinvH.GetValue().GetMaximum(), hggH.GetValue().GetMaximum()))

   if item==42: stack.SetMaximum(2*max(hData.GetValue().GetMaximum(), hWH.GetValue().GetMaximum(),hZH.GetValue().GetMaximum(), hVBFH.GetValue().GetMaximum(), hZinvH.GetValue().GetMaximum()))
   if item==42 and category == '_GFcat' and mesonCat == '_PhiCat': stack.SetMaximum(15000)    # GF

   '''
   if item==14 or item==15 or item==22: stack.SetMaximum(2*max(hData.GetValue().GetMaximum(), hWH.GetValue().GetMaximum(),hZH.GetValue().GetMaximum(), hVBFH.GetValue().GetMaximum(), hZinvH.GetValue().GetMaximum()))
   if item==102: stack.SetMaximum(10*max(hWG.GetValue().GetMaximum(),hZG.GetValue().GetMaximum(), hGJet.GetValue().GetMaximum()))
   if item==1 or item==101 or item==51 or item==52 or item==23 or item==24 or item==12 or item==41: stack.SetMaximum(2*max(hData.GetValue().GetMaximum(), hWH.GetValue().GetMaximum(),hZH.GetValue().GetMaximum(), hVBFH.GetValue().GetMaximum(), hZinvH.GetValue().GetMaximum(), hggH.GetValue().GetMaximum()))
   if item==2 and category == '_GFcat': stack.SetMaximum(30000)     #Meson Mass
   if item==17 and mesonCat == '_RhoCat': stack.SetMaximum(5000)    # MJJ
   if item==17 and mesonCat == '_PhiCat': stack.SetMaximum(600)    # MJJ
   if item==2 and (category == '_VBFcat' or category == '_VBFcatlow') : stack.SetMaximum(2000)
   if item==9 and (category == '_VBFcat' or category == '_VBFcatlow' and category == '_GFcat') : stack.SetMaximum(2*max(hData.GetMaximum(),hGJet.GetMaximum()))
   if (item==36 or item==37) and category == '_Zinvcat': stack.SetMaximum(100)    # Zinv
   if item==6 and category == '_Zinvcat': stack.SetMaximum(300)    # Zinv
   if item==1 and category == '_Zcat' and mesonCat == '_RhoCat': stack.SetMaximum(300)    # MZ
   '''
   if item!=4 and item!=43 and item!=42: stack.SetMaximum(2*hData.GetValue().GetMaximum())
   print("ALL ev: sig(ggH) = ", hggH.Integral(), ' sig(hVBFH) = ', hVBFH.Integral())

   pad1.cd()
   # Draw data first
#   if hData: hData.SetMarkerStyle(20)
#   if hData: hData.SetMarkerSize(1.2)
#   if hData: hData.SetLineWidth(2)
   if hData: hData.SetLineColor(ROOT.kBlack)
   if item==35: hData.GetXaxis().SetTitle("meson Isolation")
   if item==36: hData.GetXaxis().SetTitle("meson neutral Isolation")
   if item==37: hData.GetXaxis().SetTitle("meson neutral Isolation")
   if item==9: hData.GetXaxis().SetTitle("vertex compatibility")
   if item==21 and mesonCat == "_PhiCat": hData.GetXaxis().SetTitle("#sigmaM_{#phi}/M_{#phi}")
   if item==21 and mesonCat == "_RhoCat": hData.GetXaxis().SetTitle("#sigmaM_{#rho}/M_{#rho}")
   if item==21 and mesonCat == "_K0StarCat": hData.GetXaxis().SetTitle("#sigmaM_{K^{0*}}/M_{K^{0*}}")
   if item==41: hData.GetXaxis().SetTitle("p_{T}^{#phi,#gamma}")
   if item==42: hData.GetXaxis().SetTitle("MVA discriminator")
   #   if item==35 or item==9 or item==21:
   if item==35 or item==21:
      if hData: hData.GetYaxis().SetTitle("normalized events")
      if hData: hData.Scale(1./hData.Integral())
      if hData and item==35: hData.GetValue().SetMaximum(100*hData.GetMaximum())
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
         if mesonCat == "_K0StarCat": stack.GetXaxis().SetTitle("m_{2trk}^{ K^{0*} #rightarrow k #pi } [GeV]")
      if item==3:
         if mesonCat == "_PhiCat": stack.GetXaxis().SetTitle("pT_{2trk}^{ #phi #rightarrow k k } [GeV]")
         if mesonCat == "_RhoCat": stack.GetXaxis().SetTitle("pT_{2trk}^{ #rho #rightarrow #pi #pi } [GeV]")
         if mesonCat == "_K0StarCat": stack.GetXaxis().SetTitle("pT_{2trk}^{ K^{0*} #rightarrow k #pi } [GeV]")
      if item==16:
         if mesonCat == "_PhiCat": stack.GetXaxis().SetTitle("#DeltaR ( k k )")
         if mesonCat == "_RhoCat": stack.GetXaxis().SetTitle("#DeltaR ( #pi #pi )")
         if mesonCat == "_K0StarCat": stack.GetXaxis().SetTitle("#DeltaR ( k #pi )")
      if (item==4 or item==43) and mesonCat == "_PhiCat": stack.GetXaxis().SetTitle("m_{#gamma#phi}^{H} [GeV]")
      if (item==4 or item==43) and mesonCat == "_RhoCat": stack.GetXaxis().SetTitle("m_{#gamma#rho}^{H} [GeV]")
      if (item==4 or item==43) and mesonCat == "_K0StarCat": stack.GetXaxis().SetTitle("m_{#gammaK^{0*}}^{H} [GeV]")

      if item==11: stack.GetXaxis().SetTitle("SoftActivityJetNjets5")
      if item==12: stack.GetXaxis().SetTitle("nJet")
      if item==5: stack.GetXaxis().SetTitle("P_{T}^{#gamma}")
      if item==6: stack.GetXaxis().SetTitle("#eta^{#gamma}")
      if item==206: stack.GetXaxis().SetTitle("H/E")
      if item==207: stack.GetXaxis().SetTitle("r9")
      if item==208: stack.GetXaxis().SetTitle("pfRelIso")
      if item==209: stack.GetXaxis().SetTitle("#sigma_{i#eta i#eta}")
      if item==210: stack.GetXaxis().SetTitle("photon MVA-id")
      if item==22: stack.GetXaxis().SetTitle("|p_{T}^{miss}-p_{T}^{H_{#rho#gamma}}|/p_{T}^{H_{#rho#gamma}}")
      if item==14: stack.GetXaxis().SetTitle("#Delta #phi (p_{T}^{miss},#gamma)")
      if item==15: stack.GetXaxis().SetTitle("#Delta #phi (p_{T}^{miss},meson)")
      if item==42: stack.GetXaxis().SetTitle("MVA discriminator")
      if item==13: stack.GetXaxis().SetTitle("E_{T}^{miss}")
      if item==17: stack.GetXaxis().SetTitle("M_{JJ}")
      if item==23: stack.GetXaxis().SetTitle("p_{T}^{j1}")
      if item==24: stack.GetXaxis().SetTitle("p_{T}^{j2}")
      if item==25: stack.GetXaxis().SetTitle("#phi E_{T}^{miss}")
      if item==26: stack.GetXaxis().SetTitle("#phi Higgs")

      if item==31: stack.GetXaxis().SetTitle("p_{T}^{leadTrk}")
      if item==32: stack.GetXaxis().SetTitle("p_{T}^{subLeadTrk}")

      if item==38: stack.GetXaxis().SetTitle("p_{T}^{#pi}")
      if item==39: stack.GetXaxis().SetTitle("p_{T}^{k}")

      stack.GetYaxis().SetTitle("Events")
      if item == 4: stack.GetYaxis().SetTitle("Events/ 1 [GeV]")
      stack.GetYaxis().SetTitleOffset(1.1)
      stack.GetYaxis().SetLabelSize(0.04)
      stack.GetYaxis().SetTitleSize(0.045)
      stack.GetYaxis().ChangeLabel(1, -1, 0)

      # Draw data again
#      if hData: hData.DrawNormalized("E SAME")
      if hData: hData.Draw("E SAME")

      if item!=4 and item!=43 and item!=42:
         sigTOT = SIGstack.GetStack().Last()
         sigTOT.Scale(hData.Integral()/sigTOT.Integral())
         sigTOT.SetLineColor(ROOT.kRed)
         if sigTOT: sigTOT.Draw("hist SAME")

      pad2.cd()
      ratio = hData.Clone("dataratio")
      mcTOT = BKGstack.GetStack().Last()
      print("ALL mcTOT integral(): ",mcTOT.Integral())
      print("ALL data integral(): ",hData.Integral())

      ratio.Divide(mcTOT)
      ratio.GetYaxis().SetTitle("data/MC")
      ratio.GetYaxis().SetRangeUser(0.5,1.5)
      if (item==205): ratio.GetYaxis().SetRangeUser(0.,2.)
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
   if hData and hData.Integral()>0: legend.AddEntry(hData.GetValue(), "Data" ,"lep")

   if item!=43:
      print("hZinv=",hZinv.Integral())
      if hZinv and hZinv.Integral()>0: legend.AddEntry(hZinv.GetValue(), "Z #nu #nu", "f")
      if hZG and hZG.Integral()>0: legend.AddEntry(hZG.GetValue(), "ZG", "f")
      if hDY and hDY.Integral()>0: legend.AddEntry(hDY.GetValue(), "DY + jets", "f")
      if hWG and hWG.Integral()>0: legend.AddEntry(hWG.GetValue(), "WG", "f")
      if hW and hW.Integral()>0: legend.AddEntry(hW.GetValue(), "W + jets", "f")
      if hTTG and hTTG.Integral()>0: legend.AddEntry(hTTG.GetValue(), "TTG", "f")
      if hTT12L and hTT12L.Integral()>0: legend.AddEntry(hTT12L.GetValue(), "TT12l", "f")
      if hGJet and hGJet.Integral()>0: legend.AddEntry(hGJet.GetValue(), "#gamma + jets", "f")
      if hVBFGJet and hVBFGJet.Integral()>0: legend.AddEntry(hVBFGJet.GetValue(), "#gamma + jets (EWK)", "f")
      ######
   if hZinvH and hZinvH.Integral()>0 and mesonCat == "_PhiCat":
      if BTheo!=1.:
         legend.AddEntry(hZinvH.GetValue(), "Z#nu#nuH(#gamma#phi) SM x 10^{4}", "f")
      else:
         legend.AddEntry(hZinvH.GetValue(), "Z#nu#nuH(#gamma#phi)", "f")
   if hZH and hZH.Integral()>0 and mesonCat == "_PhiCat":
      if BTheo!=1.:
         legend.AddEntry(hZH.GetValue(), "ZH(#gamma#phi) SM x 10^{3}", "f")
      else:
         if item!=4 and item!=43 and item!=42: legend.AddEntry(hZH.GetValue(), "ZH(#gamma#phi)", "l")
         else: legend.AddEntry(hZH.GetValue(), "ZH(#gamma#phi)", "f")
   if hWH and hWH.Integral()>0 and mesonCat == "_PhiCat":
      if BTheo!=1.:
         legend.AddEntry(hWH.GetValue(), "WH(#gamma#phi) SM x 10^{4}", "f")
      else:
         if item!=4 and item!=43 and item!=42: legend.AddEntry(hWH.GetValue(), "ttH(#gamma#phi)", "l")
         else: legend.AddEntry(hWH.GetValue(), "ttH(#gamma#phi)", "f")
   if hTTH and hTTH.Integral()>0 and mesonCat == "_PhiCat":
      if BTheo!=1.:
         legend.AddEntry(hTTH.GetValue(), "ttH(#gamma#phi) SM x 10^{4}", "f")
      else:
         if item!=4 and item!=43 and item!=42: legend.AddEntry(hTTH.GetValue(), "ttH(#gamma#phi)", "l")
         else: legend.AddEntry(hTTH.GetValue(), "ttH(#gamma#phi)", "f")
   ######
   if hZinvH and hZinvH.Integral()>0 and mesonCat == "_RhoCat":
      if BTheo!=1.:
         legend.AddEntry(hZinvH.GetValue(), "Z#nu#nuH(#gamma#rho) SM x 10^{4}", "f")
      else:
         legend.AddEntry(hZinvH.GetValue(), "Z#nu#nuH(#gamma#rho)", "f")
   if hZH and hZH.Integral()>0 and mesonCat == "_RhoCat":
      if BTheo!=1.:
         legend.AddEntry(hZH.GetValue(), "ZH(#gamma#rho) SM x 10^{4}", "f")
      else:
         if item!=4 and item!=43 and item!=42: legend.AddEntry(hZH.GetValue(), "ZH(#gamma#rho)", "l")
         else: legend.AddEntry(hZH.GetValue(), "ZH(#gamma#rho)", "f")
   if hWH and hWH.Integral()>0 and mesonCat == "_RhoCat":
      if BTheo!=1.:
         legend.AddEntry(hWH.GetValue(), "WH(#gamma#rho) SM x 10^{4}", "f")
      else:
         legend.AddEntry(hWH.GetValue(), "WH(#gamma#rho)", "f")
   if hTTH and hTTH.Integral()>0 and mesonCat == "_RhoCat":
      if BTheo!=1.:
         legend.AddEntry(hTTH.GetValue(), "ttH(#gamma#rho) SM x 10^{4}", "f")
      else:
         if item!=4 and item!=43 and item!=42: legend.AddEntry(hTTH.GetValue(), "ttH(#gamma#rho)", "l")
         else: legend.AddEntry(hTTH.GetValue(), "ttH(#gamma#rho)", "f")
   ######
   if hZinvH and hZinvH.Integral()>0 and mesonCat == "_K0StarCat":
      if BTheo!=1.:
         legend.AddEntry(hZinvH.GetValue(), "Z#nu#nuH(#gammak^{0*}) SM x 10^{4}", "f")
      else:
         legend.AddEntry(hZinvH.GetValue(), "Z#nu#nuH(#gammak^{0*})", "f")
   if hZH and hZH.Integral()>0 and mesonCat == "_K0StarCat":
      if BTheo!=1.:
         legend.AddEntry(hZH.GetValue(), "ZH(#gammak^{0*}) SM x 10^{4}", "f")
      else:
         if item!=4 and item!=43 and item!=42: legend.AddEntry(hZH.GetValue(), "ZH(#gammak^{0*})", "l")
         else: legend.AddEntry(hZH.GetValue(), "ZH(#gammak^{0*})", "f")
   if hWH and hWH.Integral()>0 and mesonCat == "_K0StarCat":
      if BTheo!=1.:
         legend.AddEntry(hWH.GetValue(), "WH(#gammak^{0*}) SM x 10^{4}", "f")
      else:
         if item!=4 and item!=43 and item!=42: legend.AddEntry(hWH.GetValue(), "WH(#gammak^{0*})", "l")
         legend.AddEntry(hWH.GetValue(), "WH(#gammak^{0*})", "f")
   if hTTH and hTTH.Integral()>0 and mesonCat == "_K0StarCat":
      if BTheo!=1.:
         legend.AddEntry(hTTH.GetValue(), "ttH(#gammak^{0*}) SM x 10^{4}", "f")
      else:
         if item!=4 and item!=43 and item!=42: legend.AddEntry(hTTH.GetValue(), "ttH(#gammak^{0*})", "l")
         else: legend.AddEntry(hTTH.GetValue(), "ttH(#gammak^{0*})", "f")
   ######
   if hVBFH and mesonCat == "_PhiCat" and hVBFH.Integral()>0:
      if BTheo!=1.:
         if BTheo==0.1:
            if category=='_VBFcatlow': legend.AddEntry(hVBFH, "qqH(#gamma#phi) (BR=0.1)", "f")
            if category=='_VBFcat': legend.AddEntry(hVBFH, "qqH(#gamma#phi) (BR=0.1)", "f")
            if category=='_GFcat': legend.AddEntry(hVBFH, "qqH(#gamma#phi) (BR=0.1)", "f")
         else:
            if category=='_VBFcatlow': legend.AddEntry(hVBFH, "qqH(#gamma#phi) SM x 10^{4}", "f")
            if category=='_VBFcat': legend.AddEntry(hVBFH, "qqH(#gamma#phi) SM x 5 x 10^{3}", "f")
            if category=='_GFcat': legend.AddEntry(hVBFH, "qqH(#gamma#phi) SM x 10^{3}", "f")
      else:
         if item!=4 and item!=43 and item!=42: legend.AddEntry(hVBFH.GetValue(), "qqH(#gamma#phi)", "l")
         else: legend.AddEntry(hVBFH.GetValue(), "qqH(#gamma#phi)", "f")
   if hVBFH and mesonCat == "_RhoCat" and hVBFH.Integral()>0:
      if BTheo!=1.:
         if BTheo==0.1:
            if (category=='_VBFcatlow'): legend.AddEntry(hVBFH, "qqH(#gamma#rho) (BR=0.1)", "f")
            if category=='_GFcat' or category=='_VBFcat': legend.AddEntry(hVBFH, "qqH(#gamma#rho) (BR=0.1)", "f")
         else:
            if (category=='_VBFcatlow'): legend.AddEntry(hVBFH, "qqH(#gamma#rho) SM x 5 x 10^{3}", "f")
            if category=='_GFcat' or category=='_VBFcat': legend.AddEntry(hVBFH.GetValue(), "qqH(#gamma#rho) SM x 10^{3}", "f")
      else:
         if item!=4 and item!=43 and item!=42: legend.AddEntry(hVBFH.GetValue(), "qqH(#gamma#rho)", "l")
         else: legend.AddEntry(hVBFH.GetValue(), "qqH(#gamma#rho)", "f")
   if hVBFH and mesonCat == "_K0StarCat" and hVBFH.Integral()>0:
      if BTheo!=1.:
         if BTheo==0.1:
            if (category=='_VBFcatlow'): legend.AddEntry(hVBFH, "qqH(#gammak^{0*}) (BR=0.1)", "f")
            if category=='_GFcat' or category=='_VBFcat': legend.AddEntry(hVBFH, "qqH(#gammak^{0*}) (BR=0.1)", "f")
         else:
            if (category=='_VBFcatlow'): legend.AddEntry(hVBFH, "qqH(#gammak^{0*}) SM x 5 x 10^{3}", "f")
            if category=='_GFcat' or category=='_VBFcat': legend.AddEntry(hVBFH.GetValue(), "qqH(#gammak^{0*}) SM x 10^{3}", "f")
      else:
         if item!=4 and item!=43 and item!=42: legend.AddEntry(hVBFH.GetValue(), "qqH(#gammak^{0*})", "l")
         else: legend.AddEntry(hVBFH.GetValue(), "qqH(#gammak^{0*})", "f")
   #####
   if hggH and mesonCat == "_PhiCat" and hggH.Integral()>0:
      if BTheo!=1.:
         if BTheo==0.1:
            if (category=='_VBFcatlow'): legend.AddEntry(hggH, "ggH(#gamma#phi) (BR=0.1)", "f")
            if category=='_GFcat' or category=='_VBFcat': legend.AddEntry(hggH.GetValue(), "ggH(#gamma#phi) (BR=0.1)", "f")
         else:
            if (category=='_VBFcatlow'): legend.AddEntry(hggH, "ggH(#gamma#phi) SM x 10^{4}", "f")
            if category=='_GFcat' or category=='_VBFcat': legend.AddEntry(hggH.GetValue(), "ggH(#gamma#phi) SM x 10^{3}", "f")
      else:
         if item!=4 and item!=43 and item!=42: legend.AddEntry(hggH.GetValue(), "ggH(#gamma#phi)", "l")
         else: legend.AddEntry(hggH.GetValue(), "ggH(#gamma#phi)", "f")
   if hggH and mesonCat == "_RhoCat" and hggH.Integral()>0:
      if BTheo!=1.:
         if BTheo==0.1:
            if category=='_VBFcatlow': legend.AddEntry(hggH, "ggH(#gamma#rho) (BR=0.1)", "f")
            if category=='_GFcat' or category=='_VBFcat': legend.AddEntry(hggH, "ggH(#gamma#rho) (BR=0.1)", "f")
         else:
            if category=='_VBFcatlow': legend.AddEntry(hggH, "ggH(#gamma#rho) SM x 10^{4}", "f")
            if category=='_GFcat' or category=='_VBFcat': legend.AddEntry(hggH.GetValue(), "ggH(#gamma#rho) SM x 10^{3}", "f")
      else:
         if item!=4 and item!=43 and item!=42: legend.AddEntry(hggH.GetValue(), "ggH(#gamma#rho)", "l")
         else: legend.AddEntry(hggH.GetValue(), "ggH(#gamma#rho)", "f")
   if hggH and mesonCat == "_K0StarCat" and hggH.Integral()>0:
      if BTheo!=1.:
         if BTheo==0.1:
            if category=='_VBFcatlow': legend.AddEntry(hggH, "ggH(#gammak^{0*}) (BR=0.1)", "f")
            if category=='_GFcat' or category=='_VBFcat': legend.AddEntry(hggH, "ggH(#gammak^{0*}) (BR=0.1)", "f")
         else:
            if category=='_VBFcatlow': legend.AddEntry(hggH, "ggH(#gammak^{0*}) SM x 10^{4}", "f")
            if category=='_GFcat' or category=='_VBFcat': legend.AddEntry(hggH.GetValue(), "ggH(#gammak^{0*}) SM x 10^{3}", "f")
      else:
         if item!=4 and item!=43 and item!=42: legend.AddEntry(hggH.GetValue(), "ggH(#gammak^{0*})", "l")
         else: legend.AddEntry(hggH.GetValue(), "ggH(#gammak^{0*})", "f")

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

#   mydir = "PLOTS/"
   mydir = "/home/submit/mariadlf/public_html/Hrare/PLOTS/"

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
   if item==36: plotString = "MesonCandidate_PhoIso"
   if item==37: plotString = "MesonCandidate_NeuIso"

   if item==7: plotString = "MesonSipPV"
   if item==8: plotString = "MesonDOCA"
   if item==9: plotString = "MesonVTXprob"

   if item==16: plotString = "DRtk12"
   if item==10: plotString = "DRtk12norm"
   if item==18: plotString = "chi2dof"
   if item==21: plotString = "massErr"

   if item==40: plotString = "MesonCandidate_TrkPtRel_"
   if item==29: plotString = "MesonCandidate_Trk1Pt_norm"
   if item==30: plotString = "MesonCandidate_Trk2Pt_norm"
   if item==31: plotString = "MesonCandidate_Trk1Pt"
   if item==32: plotString = "MesonCandidate_Trk2Pt"
   if item==33: plotString = "MesonCandidate_Trk1Eta"
   if item==34: plotString = "MesonCandidate_Trk2Eta"
   if item==38: plotString = "MesonCandidate_TrkPionPt" # pion
   if item==39: plotString = "MesonCandidate_TrkKaonPt" # kaon

   if item==5: plotString = "PhotonPt"
   if item==6: plotString = "Photoneta"
   if item==205: plotString = "NPhotons"
   if item==206: plotString = "Photon_hoe"
   if item==207: plotString = "Photon_r9"
   if item==208: plotString = "Photon_relIso"
   if item==209: plotString = "Photons_sieie"
   if item==210: plotString = "Photons_mvaID"
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
   if item==25: plotString = "METphi"
   if item==26: plotString = "HiggsPHI"
   if item==22: plotString = "ptRatio"
   if item==50: plotString = "dRLepMes"

   ##
   if item==90: plotString = "dPhiPHMes"
   if item==91: plotString = "dEtaPHMes"
   if item==95: plotString = "PV_npvsGood"
   if item==96: plotString = "meson_bestVtx_idx"
   if item==113: plotString = "nBjet"

   if item==114: plotString = "Visr_mass"
   if item==115: plotString = "LeadingLeptonPt"
   if item==116: plotString = "SubLeadingLeptonPt"

#   if not isHEM:
#      c.SaveAs(mydir+plotString+string+".png")
#   else:
   if 'CR' in dirLOCAL_:
      c.SaveAs(mydir+plotString+string+"MesonMassSD.png")
   else:
      c.SaveAs(mydir+plotString+string+".png")
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

##   plot(101, 10, 0. , 10., False) # Photon pixel seed
   plot(206, 100, 0. , 0.1, True) # ph h/e
   plot(207, 60, 0.5 , 1.1, False) # ph
   plot(208, 100, 0. , 0.1, True) # ph relIso all
##   plot(209, 100, 0. , 0.1, False) # ph sieta sieta
   plot(210, 100, 0. , 1., True) # ph mvaID

   plot(5, 30, 0. , 150., False) # Photon_pt
   plot(6, 30, -3. , 3., False) # Photon_eta

   plot(205, 5, 0. , 5., True) # nPhotonLoose

def plotMeson():

   if 'CR' in dirLOCAL_:
      if mesonCat == '_PhiCat': plot(2, 40, 0.98 , 1.08, False) # phi_mass for CR
      #      if mesonCat == '_RhoCat': plot(2, 50, 0.5 , 1.04, False) # rho_mass
      if mesonCat == '_RhoCat': plot(2, 50, 0.44 , 1.1, False) # rho_mass for CR
   else:
      if mesonCat == '_PhiCat': plot(2, 40, 1. , 1.04, False) # phi_mass
      if mesonCat == '_RhoCat': plot(2, 50, 0.5 , 1.00, False) # rho_mass
      if mesonCat == '_K0StarCat': plot(2, 40, 0.8 , 1., False) # K0Star_mass for CR

   plot(96, 50, 0. , 50., True) # NVTX

   plot(3, 20, 0. , 100., False) # meson_pt

   if mesonCat == '_K0StarCat':
      plot(38, 80, 0. , 80., False) # pt_trk1
      plot(39, 80, 0. , 80., False) # pt_trk2
   else:
      if(category =='_Wcat' or category =='_Zcat'):
         plot(31, 16, 0. , 80., False) # max pt_trk1
         plot(32, 16, 0. , 80., False) # min pt_trk2
      else:
         plot(31, 80, 0. , 80., False) # max pt_trk1
         plot(32, 80, 0. , 80., False) # min pt_trk2

   return

#   plot(29, 100, 0. , 1., True) # pt_trk1_norm min
#   plot(30, 100, 0. , 1., True) # pt_trk2_norm max
   plot(40, 100, 0. , 1., True) # pt_trk_rel

   plot(33, 30, -3. , 3., False) # meson_trk1_eta
   plot(34, 30, -3. , 3., False) # meson_trk2_eta

#   if mesonCat == '_PhiCat': plot(10, 100, 0. , 1., True) # DRtk12 norm
#   if mesonCat == '_RhoCat': plot(10, 100, 0. , 5., True) # DRtk12 norm

   if mesonCat == '_PhiCat': plot(16, 50, 0. , 0.05, False) # DRtk12
   if mesonCat == '_RhoCat': plot(16, 50, 0. , 0.1, False) # DRtk12

   plot(35, 60, 0.6 , 1.2, True) # meson_iso
#   plot(36, 60, 0.6 , 1.2, True) # meson_iso neu
#   plot(37, 60, 0.6 , 1.2, True) # meson_iso neu

   plot(7, 100, 0. , 10., True) # phi_lxy
   plot(9, 30, 0. , 1.5, False) # phi_VTXprob
#   plot(8, 100, 0. , 0.1, True) # phi_doca
#   plot(18, 100, 0. , 1., False) # vtx_chi2dof

#   plot(21, 100, 0. , 0.1, True) # massErr meson

def plotWZlep():

   if(category =='_Wcat' or category =='_Zcat'):

      if(category =='_Zcat'): plot(1, 30, 75. , 105., False) # Z_mass
      if(category =='_Wcat'): plot(1, 24, 0. , 120., False) # W_mt

      if(category =='_Zcat'): plot(114, 50, 75. , 125., False) # Z_mass (llgamma)
      if(category =='_Zcat' or category =='_Wcat'): plot(115, 65, 0. , 65., False) # leadingLep
      if(category =='_Zcat'): plot(116, 65, 0. , 65., False)   #subLeadinLep

   if(category =='_Wcat'):
      plot(13, 40, 0. , 200., False) # MET_pt
      plot(14, 20, 0. , 3.14, False) # dPhiMETmeson
      plot(15, 20, 0. , 3.14, False) # dPhiMETphoton
      plot(50, 100, 0. , 10., False) # minDR (Lepton,Meson)

def plotZinv():

   if (category=='_Zinvcat'):
      plot(13, 40, 0. , 200., False) # MET_pt

      plot(22, 50, 0. , 2., False) # PTratio
      plot(14, 20, 0. , 3.14, False) # dPhiMETmeson
      plot(15, 20, 0. , 3.14, False) # dPhiMETphoton

      plot(91, 100, 0. , 10, False) # Deta Meson-Gamma
      plot(90, 100, 0. , 6.28, False) # Dphi Meson-Gamma

      plot(25, 50, -6.28 , 6.28, False) # Higgs phi
      plot(26, 100, -6.28 , 6.28, False) # Met phi

def plotVBF():

   if (category=='_VBFcatlow' or category=='_VBFcat'):
      plot(17, 60, 0. , 3000., False) # MJJ_mass

      plot(23, 20, 0. , 200., False) # jet1 Pt
      plot(24, 20, 0. , 200., False) # jet2 Pt

      plot(231, 20, -5. , 5., False) # jet1 Eta
      plot(241, 20, -5. , 5., False) # jet2 Eta

      plot(19, 100, 0. , 10, False) # Deta

      if False:
         plot(51, 100, 0. , 10., False) # minDR (jet,Meson)
         plot(52, 100, 0. , 10., False) # minDR (jet,Photon)

         plot(20, 50, -5 , 5, False) # Y1Y2

         plot(53, 50, 0. , 5., False) # zepVar
         plot(54, 100, 0. , 10., False) # dEta HIG - jet1
         plot(55, 100, 0. , 10., False) # dEta HIG - jet2
         plot(56, 100, 0. , 0.5, False) # dPhiPhi jet1
         plot(57, 100, 0. , 0.5, False) # dEtaEta jet2
         plot(58, 100, 0. , 0.5, False) # dPhiPhi jet1
         plot(59, 100, 0. , 0.5, False) # dEtaEta jet2

if __name__ == "__main__":

   plot(4, 171, 0. , 171., True)  # normal plot for the AN
#   plot(4, 71, 100. , 171., True) # HCandMass $# this is for the mass SD

   if (category=='_VBFcatlow' or category=='_VBFcat' or category=='_GFcat' or category=='_Zinvcat'):
#      plot(43, 70, 100. , 170., False) # HCandMassWithCut
#      plot(43, 100, 70. , 170., False) # HCandMassWithCut
#      plot(43, 200, 0. , 200., False) # HCandMassWithCut
      plot(42, 40, -1. , 1., False) # MVAdiscr

   plotMeson()
   plotWZlep()
   plotPhoton()
   plotVBF()

   exit()

#   plotZinv()

   exit()

   plot(41, 150, 0. , 150.,False) # HCandPT
   plot(22, 50, 0. , 2., False) # PTratio
   plot(12, 15, 0. , 15., False) # nJet
   plot(113, 15, 0. , 15., False) # nbjet
   plot(95, 100, 0. , 100., False) # NVTX
   plot(13, 40, 0. , 200., False) # MET_pt
   plot(11, 15, 0. , 15., False) # SoftActivityJetNjets5

#if __name__ == "__plots__":

#    plot(102, 100, 0.45, 0.55, False) # ks_mass
#    plot(103, 80, 0. , 80., True) # ks_pt
#    plot(107, 100, 0. , 10., True) # ks_lxy
