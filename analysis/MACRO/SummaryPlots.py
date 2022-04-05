import ROOT
import os

from LoadTree import loadTree

category = '_Wcat'
mesonCat= '_PhiCat'
year = '_2018'
directory = '/home/submit/mariadlf/Hrare/analysis/2018/APR5/'

mytree = ROOT.TChain('events')
mytree = loadTree(mytree, directory , category, mesonCat, year )

lumis={
    '_12016': 19.52, #APV #(B-F for 2016 pre)
    '_22016': 16.80, #postVFP
    '_2016': 35.9,
    '_2017': 36.4, #41.5, #(C,D,E,F for 2017)
    '_2018': 59.70,
}

# Create the plot

def plot(item, nbin, low, high, doLog):

   hDY = ROOT.TH1F( 'hDY', 'This is the distribution', nbin, low, high)
   hZG = ROOT.TH1F( 'hZG', 'This is the distribution', nbin, low, high )
   hW = ROOT.TH1F( 'hW', 'This is the distribution', nbin, low, high )
   hWG = ROOT.TH1F( 'hWG', 'This is the distribution', nbin, low, high )
   hTT12L = ROOT.TH1F( 'hTT12L', 'This is the distribution', nbin, low, high )
   hZH = ROOT.TH1F( 'hZH', 'This is the distribution', nbin, low, high )
   hWH = ROOT.TH1F( 'hWH', 'This is the distribution', nbin, low, high )
   hGJet = ROOT.TH1F( 'hGJet', 'This is the distribution', nbin, low, high )
   hJet = ROOT.TH1F( 'hJet', 'This is the distribution', nbin, low, high )
   hVBFH = ROOT.TH1F( 'hVBFH', 'This is the distribution', nbin, low, high )
   hData = ROOT.TH1F( 'hData', '', nbin, low, high )


   print(item)

   for ev in mytree:
      if item == 1 :
         var = ev.V_mass
      if item == 2 :
         var = ev.phi_kin_mass[0]
      if item == 3 :
         var = ev.phi_kin_pt[0]
      if item == 4 :
         var = ev.HCandMass
      if item == 5 :
         var = ev.Photon_pt[0]
      if item == 6 :
         var = ev.Photon_eta[0]
   # Fill histograms.
      if ev.mc==0:
         hDY.Fill( var, ev.w)
      if ev.mc==1:
         hZG.Fill( var, ev.w)
      if ev.mc==2:
         hWG.Fill( var, ev.w)
      if (ev.mc==31 or ev.mc==32 or ev.mc==33):
         hW.Fill( var, ev.w)
      if (ev.mc==4 or ev.mc==5):
         hTT12L.Fill( var, ev.w)
      if (ev.mc==6 or ev.mc==7 or ev.mc==8 or ev.mc==9):
         hGJet.Fill( var, ev.w)
      ####
      if ev.mc==1013 or ev.mc==1023:
         hZH.Fill( var, ev.w)
      if ev.mc==1011 or ev.mc==1012 or ev.mc==1021 or ev.mc==1022:
         hWH.Fill( var, ev.w)
      if ev.mc==1010 or ev.mc==1020:
         hVBFH.Fill( var, ev.w)
      ####
      if ev.mc<0:
         if (item==4 and var>115 and var<135) : continue
         hData.Fill( var )

   # Create canvas with pad
   c = ROOT.TCanvas("c", "", 600, 600)
   pad = ROOT.TPad("upper_pad", "", 0, 0, 1, 1)
   if doLog : pad.SetLogy(1)
   pad.SetTickx(False)
   pad.SetTicky(False)
   pad.Draw()
   pad.cd()
 
   redDark = (191, 34, 41)
#   redMed = (255, 40, 0)
   redMed = (237, 41, 57)
#   redMed = (220,0,5)
#   redLight = (255,61,65)
   redLight = (255,82,82)
   orange = (255, 204, 153)
   gray = (136,139,141)
   azure = (100, 192, 232)
   azureDark = (96, 147, 172)
   green = (144, 238, 144)
   greenDark = (98, 172, 141)
   gold = (212 ,175, 55)

   # Draw stackXS with MC contributions
   stack = ROOT.THStack()
   for h, color in zip([hZG, hDY, hTT12L, hWG, hW, hGJet, hJet, hZH, hWH, hVBFH], [azure, azureDark, gray, green, greenDark, orange, gold, redDark, redMed, redLight]):
      h.SetLineWidth(1)
      h.SetLineColor(1)
      h.SetFillColor(ROOT.TColor.GetColor(*color))
      h.Scale(lumis[year])
      stack.Add(h)
      if item==4: stack.SetMaximum(10*max(hWH.GetMaximum(),hZH.GetMaximum()))
      if item==2: stack.SetMaximum(6*max(hWG.GetMaximum(),hZG.GetMaximum()))
   stack.Draw("HIST")
   stack.GetXaxis().SetLabelSize(0.04)
   stack.GetXaxis().SetTitleSize(0.045)
   #stack.GetXaxis().SetTitleOffset(1.3)
   if item==1: 
      if(category =='_Wcat'): stack.GetXaxis().SetTitle("mT_{1l,MET}^{W#rightarrow #mu#nu} [GeV]")
      if(category =='_Zcat'): stack.GetXaxis().SetTitle("m_{2l}^{Z#rightarrow #mu#mu} [GeV]")
   if item==2: 
      stack.GetXaxis().SetTitle("m_{2trk}^{ #phi #rightarrow k k } [GeV]")
   if item==3: 
      stack.GetXaxis().SetTitle("pT_{2trk}^{ #phi #rightarrow k k } [GeV]")
   if item==4: stack.GetXaxis().SetTitle("m_{#gamma#phi}^{H} [GeV]")
   stack.GetYaxis().SetTitle("Events")
   stack.GetYaxis().SetLabelSize(0.04)
   stack.GetYaxis().SetTitleSize(0.045)
   stack.GetYaxis().ChangeLabel(1, -1, 0)
 
   # Draw data
   #data.SetMarkerStyle(20)
   #data.SetMarkerSize(1.2)
   #data.SetLineWidth(2)
   #data.SetLineColor(ROOT.kBlack)
   #data.Draw("E SAME")
 
   # Add legend
   legend = ROOT.TLegend(0.63, 0.65+0.05, 0.85, 0.83+0.05)
   legend.SetTextFont(42)
   legend.SetFillStyle(0)
   legend.SetBorderSize(0)
   legend.SetTextSize(0.04)
   legend.SetTextAlign(32)
   #legend.AddEntry(data, "Data" ,"lep")

   if hZG and hZG.Integral()>0: legend.AddEntry(hZG, "ZG", "f")
   if hDY and hDY.Integral()>0: legend.AddEntry(hDY, "DY + jets", "f")
   if hWG and hWG.Integral()>0: legend.AddEntry(hWG, "WG", "f")
   if hW and hW.Integral()>0: legend.AddEntry(hW, "W + jets", "f")
   if hTT12L and hTT12L.Integral()>0: legend.AddEntry(hTT12L, "TT12l", "f")
   if hGJet and hGJet.Integral()>0: legend.AddEntry(hGJet, "#gamma + jets", "f")
   if hZH and hZH.Integral()>0: legend.AddEntry(hZH, "ZH(#gamma#phi)", "f")
   if hWH and hWH.Integral()>0: legend.AddEntry(hWH, "WH(#gamma#phi)", "f")
   if hVBFH and hVBFH.Integral()>0: legend.AddEntry(hVBFH, "qqH(#gamma#phi)", "f")

   if hData and hData.Integral()>0: legend.AddEntry(hData, "Data" ,"lep")

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
   text.DrawLatex(0.65, 0.93, "#sqrt{s} = 13 TeV, X fb^{-1}")
   
   # Save the plot
   if item==1:
      c.SaveAs("VMass"+category+".png")
      print("Vmass.png")
   if item==4:
      c.SaveAs("HCandMass"+category+".png")
      print("HCandmass.png")
   if item==2:
      c.SaveAs("PhiMass"+category+".png")
      print("PhiMass.png")
   if item==3:
      c.SaveAs("PhiPt"+category+".png")
      print("PhiPt.png")
   if item==5:
      c.SaveAs("PhotonPt"+category+".png")
      print("PhotonPt.png")
   if item==6:
      c.SaveAs("PhotonEta"+category+".png")
      print("PhotonEta.png")

if __name__ == "__main__":

    if(category =='_Zcat'): plot(1, 30, 75. , 105., False) # Z_mass
    if(category =='_Wcat'): plot(1, 100, 20. , 120., False) # W_mass
    plot(4, 150, 0. , 150.,True) # HCandMass
    plot(2, 40, 1. , 1.04, False) # phi_mass
    plot(3, 80, 0. , 80., True) # phi_pt
    plot(5, 20, 0. , 100., True) # Photon_pt
    plot(6, 20, -5. , 5., False) # Photon_eta

