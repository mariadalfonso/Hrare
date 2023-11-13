import ROOT
import os
from array import array
import math
import sys

from LoadTree import loadTree
from prepareHisto import getHisto, MVAbinRho, MVAbinPhi

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

mesonCat = '_K0StarCat'
#mesonCat = '_RhoCat'
#mesonCat = '_PhiCat'
anaCat = '_GFcat'

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

def dumpSF(inFile, histoName, is1d, binY): 

  h=inFile.Get(histoName)
  if(is1d): print(h.GetName(),": ",h.Integral())

  if not is1d:
    myname = h.GetName()+"_"+str(binY)
    if binY==1: myname = h.GetName()+"_up"
    if binY==2: myname = h.GetName()+"_dn"
    h2=h.ProjectionX(myname,1+binY,1+binY)
    hn=h.ProjectionX('nom',1,1)
#    print(h2.GetName(),": ",h2.Integral(), '  --- ',round((h2.Integral()/hn.Integral()-1)*100,2),'%')

    fileToWrite="preselection"+mesonCat+anaCat+"_GGHsig.txt"
    with open(fileToWrite, "a") as f:
       str2 = h2.GetName()+": "+str(h2.Integral())+'  --- '+str(round((h2.Integral()/hn.Integral()-1)*100,2))+'% \n'
       f.write(str2)

    return h2

def plotSF(inFile,stringToPlot):

  print('---------------------------------')

  hNom=dumpSF(inFile, stringToPlot, False,0)
  hUp=dumpSF(inFile, stringToPlot, False,1)
  hDn=dumpSF(inFile, stringToPlot,False,2)

  doLog = False
  c, pad1, pad2 = createCanvasPads(doLog)
  pad1.cd()
  hNom.SetLineColor(ROOT.kRed)
  hNom.SetMaximum(1.10*max(hNom.GetMaximum(),hUp.GetMaximum(),hDn.GetMaximum()))
  hNom.Draw("hist")
  hUp.SetLineColor(ROOT.kGreen)
  hUp.Draw("hist same")
  hDn.SetLineColor(ROOT.kBlue)
  hDn.Draw("hist same")

  pad2.cd()
  ratio = hUp.Clone("dataratio")
  ratio.GetYaxis().SetTitle("Var/Nom")
  ratio.GetYaxis().SetRangeUser(0.5,1.5)
  ratio.Divide(hNom)
  ratio.GetXaxis().SetTitleOffset(4.)
  ratio.GetXaxis().SetTitleSize(0.15)
  ratio.GetXaxis().SetLabelSize(0.12)

  ratio.GetYaxis().SetTitleOffset(0.3)
  ratio.GetYaxis().SetTitleSize(0.15)
  ratio.GetYaxis().SetLabelSize(0.12)

  ratio.Draw("hist same")
  ratio2 = hDn.Clone("dataratio")
  ratio2.Divide(hNom)
  ratio2.Draw("hist same")

  lineZero = ROOT.TLine(hNom.GetXaxis().GetXmin(), 1.,  hNom.GetXaxis().GetXmax(), 1.)
  lineZero.SetLineColor(11)
  lineZero.Draw("same")

  c.SaveAs("~/public_html/Hrare/SYSTstudies/"+stringToPlot+".png")


if __name__ == "__main__":

  isW=False
  isZ=False

  if anaCat == '_GFcat':
     if mesonCat == '_RhoCat': inFile = ROOT.TFile.Open("/work/submit/mariadlf/OCT_SYSTvar/2018/histoname_mc1027"+anaCat+"_RhoCat_2018_wSF.root","READ")
     if mesonCat == '_PhiCat': inFile = ROOT.TFile.Open("/work/submit/mariadlf/OCT_SYSTvar/2018/histoname_mc1017"+anaCat+"_PhiCat_2018_wSF.root","READ")
     if mesonCat == '_K0StarCat': inFile = ROOT.TFile.Open("/work/submit/mariadlf/OCT_SYSTvar/2018/histoname_mc1037"+anaCat+"_K0StarCat_2018_wSF.root","READ")

#  if anaCat == '_VBFcatlow' or anaCat == '_VBFcat' or anaCat == '_GFcat':
#     if mesonCat == '_RhoCat': inFile = ROOT.TFile.Open("/work/submit/mariadlf/OCT_SYSTvar/2018/histoname_mc1020"+anaCat+"_RhoCat_2018_wSF.root","READ")
#     if mesonCat == '_PhiCat': inFile = ROOT.TFile.Open("/work/submit/mariadlf/OCT_SYSTvar/2018/histoname_mc1010"+anaCat+"_PhiCat_2018_wSF.root","READ")
#     if mesonCat == '_K0StarCat': inFile = ROOT.TFile.Open("/work/submit/mariadlf/OCT_SYSTvar/2018/histoname_mc1030"+anaCat+"_K0StarCat_2018_wSF.root","READ")

  inFile.ls()

  dumpSF(inFile, "HCandMass",True,-1)

  dumpSF(inFile, "HCandMass:JetSYST:dn",True,-1)
  dumpSF(inFile, "HCandMass:JetSYST:up",True,-1)
  print('---------------------------------')
  dumpSF(inFile, "HCandMass:PhotonSYST:dn",True,-1)
  dumpSF(inFile, "HCandMass:PhotonSYST:up",True,-1)    
  print('---------------------------------')

  plotSF(inFile,"HCandMass:PU")
  plotSF(inFile,"HCandMass:L1")
  plotSF(inFile,"HCandMass:phoID")
  plotSF(inFile,"HCandMass:PSfsr")
  plotSF(inFile,"HCandMass:PSisr")
  plotSF(inFile,"HCandMass:mesonRECO")
  plotSF(inFile,"HCandMass:mesonChISO")
  plotSF(inFile,"HCandMass:phoTrig")
  plotSF(inFile,"HCandMass:tauTrig")

  '''
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
  #
  hNom=dumpSF(inFile, "HCandMass:mesonChISO",False,0)
  hUp=dumpSF(inFile, "HCandMass:mesonChISO",False,1)
  hDn=dumpSF(inFile, "HCandMass:mesonChISO",False,2)
  print('---------------------------------')
  #
  hNom=dumpSF(inFile, "HCandMass:phoTrig",False,0)
  hUp=dumpSF(inFile, "HCandMass:phoTrig",False,1)
  hDn=dumpSF(inFile, "HCandMass:phoTrig",False,2)
  print('---------------------------------')
  #
  hNom=dumpSF(inFile, "HCandMass:tauTrig",False,0)
  hUp=dumpSF(inFile, "HCandMass:tauTrig",False,1)
  hDn=dumpSF(inFile, "HCandMass:tauTrig",False,2)
  '''

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
