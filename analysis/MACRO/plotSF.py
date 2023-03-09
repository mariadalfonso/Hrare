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
    
if __name__ == "__main__":

  inFile = ROOT.TFile.Open("../SYSTstudy/2018/histoname_mc1022_Wcat_RhoCat_2018.root","READ")
#  inFile = ROOT.TFile.Open("../SYSTstudy/2018/histoname_mc1023_Zcat_RhoCat_2018.root","READ")
#  inFile = ROOT.TFile.Open("../SYSTstudy/2018/histoname_mc1027_GFcat_RhoCat_2018.root","READ")
  inFile.ls()

  dumpSF(inFile, "HCandMass",True,-1)
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
