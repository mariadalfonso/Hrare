import ROOT
import os
from array import array
import math
import sys

from LoadTree import loadTree
from prepareHisto import getHisto, MVAbinRho, MVAbinPhi

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()


def dumpSF(inFile, histoName, is1d, binY): 

  h=inFile.Get(histoName)
  if(is1d): print(h.GetName(),": ",h.Integral())

  if not is1d:
    myname = h.GetName()+"_"+str(binY)
    if binY==1: myname = h.GetName()+"_up"
    if binY==2: myname = h.GetName()+"_dn"
    h2=h.ProjectionX(myname,1+binY,1+binY)
    hn=h.ProjectionX('nom',1,1)
    print(h2.GetName(),": ",h2.Integral(), '  --- ',round((h2.Integral()/hn.Integral()-1)*100,2),'%')

    return h2
    
if __name__ == "__main__":

  isW=False
  isZ=False

#  inFile = ROOT.TFile.Open("../TEST/2018/histoname_mc1027_GFcat_RhoCat_2018_wSF.root","READ")
#  inFile = ROOT.TFile.Open("../TEST/2018/histoname_mc1020_GFcat_RhoCat_2018_wSF.root","READ")
  inFile = ROOT.TFile.Open("../TEST/2018/histoname_mc1020_VBFcat_RhoCat_2018_wSF.root","READ")

#  inFile = ROOT.TFile.Open("../TEST/2018/histoname_mc1017_GFcat_PhiCat_2018_wSF.root","READ")
#  inFile = ROOT.TFile.Open("../TEST/2018/histoname_mc1027_GFcat_RhoCat_2018_wSF.root","READ")

#  inFile = ROOT.TFile.Open("../TEST/2018/histoname_mc1020_VBFcat_RhoCat_2018_wSF.root","READ")
#  inFile = ROOT.TFile.Open("../TEST/2018/histoname_mc1010_VBFcat_PhiCat_2018_wSF.root","READ")

#  inFile = ROOT.TFile.Open("../TEST/2017/histoname_mc1020_VBFcat_RhoCat_2017_wSF.root","READ")
#  inFile = ROOT.TFile.Open("../TEST/2017/histoname_mc1010_VBFcat_PhiCat_2017_wSF.root","READ")

#  inFile = ROOT.TFile.Open("../TEST/12016/histoname_mc1020_VBFcat_RhoCat_12016_wSF.root","READ")
#  inFile = ROOT.TFile.Open("../TEST/12016/histoname_mc1010_VBFcat_PhiCat_12016_wSF.root","READ")

#  inFile = ROOT.TFile.Open("../TEST/2018/histoname_mc1022_Wcat_RhoCat_2018_wSF.root","READ")
#  inFile = ROOT.TFile.Open("../TEST/2017/histoname_mc1022_Wcat_RhoCat_2017_wSF.root","READ")
#  inFile = ROOT.TFile.Open("../TEST/22016/histoname_mc1022_Wcat_RhoCat_22016_wSF.root","READ")
#  inFile = ROOT.TFile.Open("../TEST/12016/histoname_mc1022_Wcat_RhoCat_12016_wSF.root","READ")

#  inFile = ROOT.TFile.Open("../TEST/2018/histoname_mc1023_Zcat_RhoCat_2018_wSF.root","READ")
#  inFile = ROOT.TFile.Open("../TEST/2017/histoname_mc1023_Zcat_RhoCat_2017_wSF.root","READ")
#  inFile = ROOT.TFile.Open("../TEST/22016/histoname_mc1023_Zcat_RhoCat_22016_wSF.root","READ")
#  inFile = ROOT.TFile.Open("../TEST/12016/histoname_mc1023_Zcat_RhoCat_12016_wSF.root","READ")

#  inFile = ROOT.TFile.Open("../SYSTstudy/2018/histoname_mc1023_Zcat_RhoCat_2018_wSF.root","READ")
#  inFile = ROOT.TFile.Open("../SYSTstudy/2018/histoname_mc1027_GFcat_RhoCat_2018_wSF.root","READ")
  inFile.ls()

  dumpSF(inFile, "HCandMass",True,-1)

  dumpSF(inFile, "HCandMass:JetSYST:dn",True,-1)
  dumpSF(inFile, "HCandMass:JetSYST:up",True,-1)
  print('---------------------------------')
  dumpSF(inFile, "HCandMass:PhotonSYST:dn",True,-1)
  dumpSF(inFile, "HCandMass:PhotonSYST:up",True,-1)    
  print('---------------------------------')
  #
  dumpSF(inFile, "HCandMass:PU",False,0)
  dumpSF(inFile, "HCandMass:PU",False,1)
  dumpSF(inFile, "HCandMass:PU",False,2)
  print('---------------------------------')
  #
  dumpSF(inFile, "HCandMass:L1",False,0)
  dumpSF(inFile, "HCandMass:L1",False,1)
  dumpSF(inFile, "HCandMass:L1",False,2)
  print('---------------------------------')
  #
  dumpSF(inFile, "HCandMass:phoID",False,0)
  dumpSF(inFile, "HCandMass:phoID",False,1)
  dumpSF(inFile, "HCandMass:phoID",False,2)
  print('---------------------------------')
  #
  hNom=dumpSF(inFile, "HCandMass:mesonRECO",False,0)
  hUp=dumpSF(inFile, "HCandMass:mesonRECO",False,1)
  hDn=dumpSF(inFile, "HCandMass:mesonRECO",False,2)
  print('---------------------------------')
  #
  hNom=dumpSF(inFile, "HCandMass:PSisr",False,0)
  hUp=dumpSF(inFile, "HCandMass:PSisr",False,1)
  hDn=dumpSF(inFile, "HCandMass:PSisr",False,2)
  print('---------------------------------')
  #
  hNom=dumpSF(inFile, "HCandMass:PSfsr",False,0)
  hUp=dumpSF(inFile, "HCandMass:PSfsr",False,1)
  hDn=dumpSF(inFile, "HCandMass:PSfsr",False,2)
  print('---------------------------------')

  canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
  hNom.SetMaximum(1.10*max(hNom.GetMaximum(),hUp.GetMaximum(),hDn.GetMaximum()))
  hNom.Draw("hist")
  hUp.SetLineColor(ROOT.kRed)
  hUp.Draw("hist same")
  hDn.SetLineColor(ROOT.kRed)
  hDn.Draw("hist same")
  canvas.SaveAs("MesonRECO.png")


  if isW or isZ:
    #
    dumpSF(inFile, "HCandMass:eleID",False,0)
    dumpSF(inFile, "HCandMass:eleID",False,1)
    dumpSF(inFile, "HCandMass:eleID",False,2)
    print('---------------------------------')
    dumpSF(inFile, "HCandMass:eleRECO",False,0)
    dumpSF(inFile, "HCandMass:eleRECO",False,1)
    dumpSF(inFile, "HCandMass:eleRECO",False,2)
    print('---------------------------------')
    #
    dumpSF(inFile, "HCandMass:muoID",False,0)
    dumpSF(inFile, "HCandMass:muoID",False,1)
    dumpSF(inFile, "HCandMass:muoID",False,2)
    print('---------------------------------')
    #
    dumpSF(inFile, "HCandMass:muoISO",False,0)
    dumpSF(inFile, "HCandMass:muoISO",False,1)
    dumpSF(inFile, "HCandMass:muoISO",False,2)
