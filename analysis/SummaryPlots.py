import ROOT
import os

category = '_Wcat'
directory = ''
mytree = ROOT.TChain('events')
mytree.Add(directory+'outname_mc1'+category+'.root')  #ZG
mytree.Add(directory+'outname_mc2'+category+'.root')  #WG
mytree.Add(directory+'outname_mc10'+category+'.root') #ZH
mytree.Add(directory+'outname_mc11'+category+'.root') #WH
mytree.Add(directory+'outname_mc3'+category+'.root')  #Wjets
mytree.Add(directory+'outname_mc4'+category+'.root')  #tt2l
mytree.Add(directory+'outname_mc0'+category+'.root')  #DY

lumi=137

# Create the plot

def plot(item, nbin, low, high, doLog):

   hDY = ROOT.TH1F( 'hDY', 'This is the distribution', nbin, low, high)
   hZG = ROOT.TH1F( 'hZG', 'This is the distribution', nbin, low, high )
   hW = ROOT.TH1F( 'hW', 'This is the distribution', nbin, low, high )
   hWG = ROOT.TH1F( 'hWG', 'This is the distribution', nbin, low, high )
   hTT2L = ROOT.TH1F( 'hTT2L', 'This is the distribution', nbin, low, high )
   hZH = ROOT.TH1F( 'hZH', 'This is the distribution', nbin, low, high )
   hWH = ROOT.TH1F( 'hWH', 'This is the distribution', nbin, low, high )

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
      if ev.mc==3:
         hW.Fill( var, ev.w)
      if ev.mc==4:
         hTT2L.Fill( var, ev.w)
      if ev.mc==10:
         hZH.Fill( var, ev.w)
      if ev.mc==11:
         hWH.Fill( var, ev.w)
                  
   # Create canvas with pad
   c = ROOT.TCanvas("c", "", 600, 600)
   pad = ROOT.TPad("upper_pad", "", 0, 0, 1, 1)
   if doLog : pad.SetLogy(1)
   pad.SetTickx(False)
   pad.SetTicky(False)
   pad.Draw()
   pad.cd()
 
   # Draw stackXS with MC contributions
   stack = ROOT.THStack()
   for h, color in zip([hZG, hDY, hTT2L, hWG, hW, hZH, hWH], [(100, 192, 232), (96, 147, 172), (136,139,141), (144, 238, 144), (98,  172, 141), (191, 34, 41), (237, 41, 57)]):
      h.SetLineWidth(1)
      h.SetLineColor(1)
      h.SetFillColor(ROOT.TColor.GetColor(*color))
      h.Scale(lumi)
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

   if hZG: legend.AddEntry(hZG, "ZG", "f")
   if hDY: legend.AddEntry(hDY, "DY + jets", "f")
   if hWG: legend.AddEntry(hWG, "WG", "f")
   if hW: legend.AddEntry(hW, "W + jets", "f")
   if hTT2L: legend.AddEntry(hTT2L, "TT2l", "f")
   if hZH: legend.AddEntry(hZH, "ZH(#gamma#phi)", "f")
   if hWH: legend.AddEntry(hWH, "WH(#gamma#phi)", "f")
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

