import ROOT
import os
from array import array

from LoadTree import loadTree
from LoadTree import resetTree

ROOT.gStyle.SetOptStat(0)

lumis={
    '_12016': 19.52, #APV #(B-F for 2016 pre)
    '_22016': 16.80, #postVFP
    '_2016': 35.9,
    '_2017': 36.4, #41.5, #(C,D,E,F for 2017)
    '_2018': 59.70,
    '_12018': 39.54,
}

# Create the plot

def addSigHisto(mytree, directory, category, mesonCat, nameSig, year ):


   if(category =='_GFcat' and mesonCat == '_K0StarCat' and nameSig=="ggH"): mytree.Add(directory+'outname_mc1037'+category+mesonCat+year+'.root') # GF k0Star
#   if(category =='_GFcat' and mesonCat == '_K0StarCat' and nameSig=="ggH"): mytree.Add(directory+'outname_mc1027'+category+mesonCat+year+'.root') # GF rho
#   if(category =='_GFcat' and mesonCat == '_K0StarCat' and nameSig=="ggH"): mytree.Add(directory+'outname_mc1017'+category+mesonCat+year+'.root') # GF phi
   if(category =='_GFcat' and mesonCat == '_K0StarCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1030'+category+mesonCat+year+'.root') # VBF k0Star
#   if(category =='_GFcat' and mesonCat == '_K0StarCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1020'+category+mesonCat+year+'.root') # VBF rho
#   if(category =='_GFcat' and mesonCat == '_K0StarCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1010'+category+mesonCat+year+'.root') # VBF phi
   if(category =='_VBFcat' and mesonCat == '_K0StarCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1030'+category+mesonCat+year+'.root') # VBF k0Star
#   if(category =='_VBFcat' and mesonCat == '_K0StarCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1020'+category+mesonCat+year+'.root') # VBF rho
#   if(category =='_VBFcat' and mesonCat == '_K0StarCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1010'+category+mesonCat+year+'.root') # VBF phi
   if(category =='_VBFcatlow' and mesonCat == '_K0StarCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1030'+category+mesonCat+year+'.root') # VBF k0Star
#   if(category =='_VBFcatlow' and mesonCat == '_K0StarCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1020'+category+mesonCat+year+'.root') # VBF rho
#   if(category =='_VBFcatlow' and mesonCat == '_K0StarCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1010'+category+mesonCat+year+'.root') # VBF phi
   # for phi should add the 1027 1020 series
   if(category =='_GFcat' and mesonCat == '_PhiCat' and nameSig=="ggH"): mytree.Add(directory+'outname_mc1017'+category+mesonCat+year+'.root') # VBF
   if(category =='_GFcat' and mesonCat == '_PhiCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1010'+category+mesonCat+year+'.root') # VBF
   if(category =='_VBFcat' and mesonCat == '_PhiCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1010'+category+mesonCat+year+'.root') # VBF
   if(category =='_VBFcatlow' and mesonCat == '_PhiCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1010'+category+mesonCat+year+'.root') # VBF
   if(category =='_Wcat' and mesonCat == '_PhiCat' and nameSig=="WH"): mytree.Add(directory+'outname_mc1011'+category+mesonCat+year+'.root') # Wp
   if(category =='_Wcat' and mesonCat == '_PhiCat' and nameSig=="WH"): mytree.Add(directory+'outname_mc1012'+category+mesonCat+year+'.root') # Wm
   if(category =='_Wcat' and mesonCat == '_PhiCat' and nameSig=="ZHl"): mytree.Add(directory+'outname_mc1013'+category+mesonCat+year+'.root') # Z (lostL)
   if(category =='_Wcat' and mesonCat == '_PhiCat' and nameSig=="ZHl"): mytree.Add(directory+'outname_mc1014'+category+mesonCat+year+'.root') # Z (lostL)
   if(category =='_Wcat' and mesonCat == '_PhiCat' and nameSig=="TTH"): mytree.Add(directory+'outname_mc1018'+category+mesonCat+year+'.root') # ttH
   if(category =='_Zcat' and mesonCat == '_PhiCat' and nameSig=="ZH"): mytree.Add(directory+'outname_mc1013'+category+mesonCat+year+'.root') # Z
   if(category =='_Zcat' and mesonCat == '_PhiCat' and nameSig=="ZH"): mytree.Add(directory+'outname_mc1014'+category+mesonCat+year+'.root') # Z
   if(category =='_Zinvcat' and mesonCat == '_PhiCat' and nameSig=="ZinvH"): mytree.Add(directory+'outname_mc1015'+category+mesonCat+year+'.root') # Zinv
   if(category =='_Zinvcat' and mesonCat == '_PhiCat' and nameSig=="ZinvH"): mytree.Add(directory+'outname_mc1016'+category+mesonCat+year+'.root') # Zinv
   if(category =='_Zinvcat' and mesonCat == '_PhiCat' and nameSig=="WHl"): mytree.Add(directory+'outname_mc1011'+category+mesonCat+year+'.root') # Zinv Wp (lostL)
   if(category =='_Zinvcat' and mesonCat == '_PhiCat' and nameSig=="WHl"): mytree.Add(directory+'outname_mc1012'+category+mesonCat+year+'.root') # Zinv Wm (lostL)
   if(category =='_Vcat' and mesonCat == '_PhiCat' and nameSig=="WH"): mytree.Add(directory+'outname_mc1011_Wcat'+mesonCat+year+'.root') # Wp
   if(category =='_Vcat' and mesonCat == '_PhiCat' and nameSig=="WH"): mytree.Add(directory+'outname_mc1012_Wcat'+mesonCat+year+'.root') # Wm
   if(category =='_Vcat' and mesonCat == '_PhiCat' and nameSig=="TTH"): mytree.Add(directory+'outname_mc1018_Wcat'+mesonCat+year+'.root') # ttH
   if(category =='_Vcat' and mesonCat == '_PhiCat' and nameSig=="ZHl"): mytree.Add(directory+'outname_mc1013_Wcat'+mesonCat+year+'.root') # Z (lostL)
   if(category =='_Vcat' and mesonCat == '_PhiCat' and nameSig=="ZHl"): mytree.Add(directory+'outname_mc1014_Wcat'+mesonCat+year+'.root') # Z (lostL)
   if(category =='_Vcat' and mesonCat == '_PhiCat' and nameSig=="ZH"): mytree.Add(directory+'outname_mc1013_Zcat'+mesonCat+year+'.root') # Z
   if(category =='_Vcat' and mesonCat == '_PhiCat' and nameSig=="ZH"): mytree.Add(directory+'outname_mc1014_Zcat'+mesonCat+year+'.root') # Z
   #
   # for rho should add the 1017 1010 series
   if(category =='_GFcat' and mesonCat == '_RhoCat' and nameSig=="ggH"): mytree.Add(directory+'outname_mc1027'+category+mesonCat+year+'.root') # VBF
   if(category =='_GFcat' and mesonCat == '_RhoCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1020'+category+mesonCat+year+'.root') # VBF
   if(category =='_VBFcat' and mesonCat == '_RhoCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1020'+category+mesonCat+year+'.root') # VBF
   if(category =='_VBFcatlow' and mesonCat == '_RhoCat' and nameSig=="VBFH"): mytree.Add(directory+'outname_mc1020'+category+mesonCat+year+'.root') # VBF
   if(category =='_Wcat' and mesonCat == '_RhoCat' and nameSig=="WH"): mytree.Add(directory+'outname_mc1021'+category+mesonCat+year+'.root') # Wp
   if(category =='_Wcat' and mesonCat == '_RhoCat' and nameSig=="WH"): mytree.Add(directory+'outname_mc1022'+category+mesonCat+year+'.root') # Wm
   if(category =='_Wcat' and mesonCat == '_RhoCat' and nameSig=="ZHl"): mytree.Add(directory+'outname_mc1023'+category+mesonCat+year+'.root') # Z (lostL)
   if(category =='_Wcat' and mesonCat == '_RhoCat' and nameSig=="ZHl"): mytree.Add(directory+'outname_mc1024'+category+mesonCat+year+'.root') # Z (lostL)
   if(category =='_Wcat' and mesonCat == '_RhoCat' and nameSig=="TTH"): mytree.Add(directory+'outname_mc1028'+category+mesonCat+year+'.root') # ttH
   if(category =='_Zcat' and mesonCat == '_RhoCat' and nameSig=="ZH"): mytree.Add(directory+'outname_mc1023'+category+mesonCat+year+'.root') # Z
   if(category =='_Zcat' and mesonCat == '_RhoCat' and nameSig=="ZH"): mytree.Add(directory+'outname_mc1024'+category+mesonCat+year+'.root') # Z
   if(category =='_Zinvcat' and mesonCat == '_RhoCat' and nameSig=="ZinvH"): mytree.Add(directory+'outname_mc1025'+category+mesonCat+year+'.root') # Zinv
   if(category =='_Zinvcat' and mesonCat == '_RhoCat' and nameSig=="ZinvH"): mytree.Add(directory+'outname_mc1026'+category+mesonCat+year+'.root') # Zinv
   if(category =='_Zinvcat' and mesonCat == '_RhoCat' and nameSig=="WHl"): mytree.Add(directory+'outname_mc1021'+category+mesonCat+year+'.root') # Zinv Wp (lostL)
   if(category =='_Zinvcat' and mesonCat == '_RhoCat' and nameSig=="WHl"): mytree.Add(directory+'outname_mc1022'+category+mesonCat+year+'.root') # Zinv Wm (lostL)
   if(category =='_Vcat' and mesonCat == '_RhoCat' and nameSig=="WH"): mytree.Add(directory+'outname_mc1021_Wcat'+mesonCat+year+'.root') # Wp
   if(category =='_Vcat' and mesonCat == '_RhoCat' and nameSig=="WH"): mytree.Add(directory+'outname_mc1022_Wcat'+mesonCat+year+'.root') # Wm
   if(category =='_Vcat' and mesonCat == '_RhoCat' and nameSig=="TTH"): mytree.Add(directory+'outname_mc1028_Wcat'+mesonCat+year+'.root') # ttH
   if(category =='_Vcat' and mesonCat == '_RhoCat' and nameSig=="ZHl"): mytree.Add(directory+'outname_mc1023_Wcat'+mesonCat+year+'.root') # Z (lostL)
   if(category =='_Vcat' and mesonCat == '_RhoCat' and nameSig=="ZHl"): mytree.Add(directory+'outname_mc1024_Wcat'+mesonCat+year+'.root') # Z (lostL)
   if(category =='_Vcat' and mesonCat == '_RhoCat' and nameSig=="ZH"): mytree.Add(directory+'outname_mc1023_Zcat'+mesonCat+year+'.root') # Z
   if(category =='_Vcat' and mesonCat == '_RhoCat' and nameSig=="ZH"): mytree.Add(directory+'outname_mc1024_Zcat'+mesonCat+year+'.root') # Z
   # maybe add the Zinv to the GF

   resetTree(mytree, category)

MVAbinRho={
   '_GFcat': 0.4,
   '_VBFcat': 0.6,
   '_VBFcatlow': 0.3,
}

MVAbinPhi={
   '_GFcat': 0.4,
   '_VBFcat': 0.6,
   '_VBFcatlow': 0.3,
}

#MVAbinRho={
#   '_GFcat': 0.5,
#   '_Zinvcat': 0.6,
#   '_VBFcat': 0.8,
#   '_VBFcatlow': 0.8,
#}

#MVAbinPhi={
#   '_GFcat': 0.4,
#   '_Zinvcat': 0.7,
#   '_VBFcat': 0.6,
#   '_VBFcatlow': 0.6,
#}

def getHisto(item, nbin, low, high, doLog,category,mesonCat, doSignal, nameSig):

   if mesonCat == '_K0StarCat': MVAbin = MVAbinRho
   if mesonCat == '_RhoCat': MVAbin = MVAbinRho
   if mesonCat == '_PhiCat': MVAbin = MVAbinPhi

   if category =='_VBFcat' or category =='_Wcat' or category =='_Zcat' :
      year = '_all'
   else:
      year = '_2018'

   mytree = ROOT.TChain('events')

   dirLOCAL_='/work/submit/mariadlf/JUL22/'

   directory1 = dirLOCAL_+'12016/'
   directory2 = dirLOCAL_+'22016/'
   directory3 = dirLOCAL_+'2017/'
   directory4 = dirLOCAL_+'2018/'

   if doSignal:
      if(year == '_2018' or year == '_all'): addSigHisto(mytree, directory4, category, mesonCat, nameSig, '_2018')
      if(year == '_2017' or year == '_all'): addSigHisto(mytree, directory3, category, mesonCat, nameSig, '_2017')
      if(category =='_Wcat' or category =='_Zcat'):
         if(year == '_22016' or year == '_all'): addSigHisto(mytree, directory2, category, mesonCat, nameSig, '_22016')
      if(year == '_12016' or year == '_all'): addSigHisto(mytree, directory1, category, mesonCat, nameSig, '_12016')
   else:
      if(category =='_Vcat'):
         if(year == '_2018' or year == '_all'):
            mytree = loadTree(mytree, directory4 , '_Wcat', mesonCat, '_2018', doSignal) # all sigm, BKG , data,
            mytree = loadTree(mytree, directory4 , '_Zcat',  mesonCat, '_2018', doSignal) # all sigm, BKG , data,
         if(year == '_2017' or year == '_all'):
            mytree = loadTree(mytree, directory3 , '_Wcat', mesonCat, '_2017', doSignal) # all sigm, BKG , data,
            mytree = loadTree(mytree, directory3 , '_Zcat', mesonCat, '_2017', doSignal) # all sigm, BKG , data,
         if(category =='_Wcat' or category =='_Zcat'):
            if(year == '_22016' or year == '_all'):
               mytree = loadTree(mytree, directory2 , '_Wcat', mesonCat, '_22016', doSignal) # all sigm, BKG , data,
               mytree = loadTree(mytree, directory2 , '_Zcat', mesonCat, '_22016', doSignal) # all sigm, BKG , data,
         if(year == '_12016' or year == '_all'):
            mytree = loadTree(mytree, directory1 , '_Wcat', mesonCat, '_12016', doSignal) # all sigm, BKG , data,
            mytree = loadTree(mytree, directory1 , '_Zcat', mesonCat, '_12016', doSignal) # all sigm, BKG , data,
      else:
         if(year == '_2018' or year == '_all'): mytree = loadTree(mytree, directory4 , category, mesonCat, '_2018', doSignal) # all sigm, BKG , data,
         if(year == '_2017' or year == '_all'): mytree = loadTree(mytree, directory3 , category, mesonCat, '_2017', doSignal) # all sigm, BKG , data,
         if(category =='_Wcat' or category =='_Zcat'):
            if(year == '_22016' or year == '_all'): mytree = loadTree(mytree, directory2 , category, mesonCat, '_22016', doSignal) # all sigm, BKG , data,
         if(year == '_12016' or year == '_all'): mytree = loadTree(mytree, directory1 , category, mesonCat, '_12016', doSignal) # all sigm, BKG , data,

   h = ROOT.TH1F( 'Higgs', '', nbin, low, high )

   PH_ptCut = 40
   MES_ptCut = 40

   for ev in mytree:

      if abs(item) == 43:
         var = -1
         if item> 0:
            if ev.MVAdisc[0]>MVAbin[category]:  var = ev.HCandMass
         else:
            if ev.MVAdisc[0]<MVAbin[category]:  var = ev.HCandMass

      if item == 4:
         var = ev.HCandMass

#      idxMeson = ev.index_pair[0]
#      idxPh = ev.index_pair[1]

#      leadTrkPt = max(ev.goodMeson_trk1_pt[idxMeson],ev.goodMeson_trk2_pt[idxMeson]) # already applied
#      if leadTrkPt < 20: continue

#      if ev.goodMeson_iso[idxMeson] < 0.9 :  continue  # already applied 

#      if mesonCat == '_RhoCat' and abs(ev.goodMeson_mass[idxMeson]-0.77) > 2.*0.07: continue # 2 sigma # already applied  
#      if mesonCat == '_PhiCat' and abs(ev.goodMeson_mass[idxMeson]-1.02) > 3.*0.004: continue # 3. sigma # already applied

      ## OPTIMIZED PHASE SPACE
      if(category =='_Zcat' or category =='_Wcat' or category =='_Vcat'):
#         if ev.goodPhotons_pt[idxPh]<40 :  continue
#         if ev.goodMeson_pt[idxMeson]<40 : continue
         if ev.photon_pt<PH_ptCut :  continue
         if ev.meson_pt<MES_ptCut : continue

#      if(category =='_Wcat'):
#         if (abs(ev.Z_veto-91) < 5): continue # default cutting at +- 10

      ## ALREADY applied for Vcat in MARCH11
#      if(category =='_Wcat'):
#         if ev.V_mass < 15 : continue
#         if ev.DeepMETResolutionTune_pt<15 :  continue
#         if ev.deltaLepMeson<0.5 :  continue

      ## OPTIMIZED PHASE SPACE
      if(category =='_Zinvcat'):
#         if ev.goodPhotons_pt[idxPh]<40 :  continue
#         if ev.goodMeson_pt[idxMeson]<40 : continue
         if ev.photon_pt<PH_ptCut :  continue
         if ev.meson_pt<MES_ptCut : continue
#         if min(ev.goodMeson_trk1_pt[idxMeson],ev.goodMeson_trk2_pt[idxMeson]) < 10: continue
         if ev.DeepMETResolutionTune_pt<75 :  continue # at skim level is 50
         if abs(ev.dPhiGammaMET)<1. : continue # already applied
         if abs(ev.dPhiMesonMET)<1. : continue # already applied
         if abs(item) == 4:
            if ev.ptRatioMEThiggs>0.5: continue
            if abs(ev.dEtaGammaMesonCand)>2: continue
         if ev.nbtag>0.: continue

      ## OPTIMIZED PHASE SPACE
      if(category =='_VBFcatlow'):
#         if ev.goodPhotons_pt[idxPh]<40 :  continue # already applied
#         if ev.goodPhotons_pt[idxPh]>75 :  continue # already applied
#         if ev.goodMeson_pt[idxMeson]<40 : continue
         if ev.photon_pt<PH_ptCut :  continue # already applied
         if ev.photon_pt>75 :  continue # already applied
         if ev.meson_pt<MES_ptCut : continue

      ## OPTIMIZED PHASE SPACE
      if(category =='_VBFcat' or category =='_VBFcatlow'):
         if ev.jet1Pt < 30 : continue
         if ev.jet2Pt < 20 : continue
         if ev.deltaJetMeson < 1.2 : continue
         if ev.deltaJetPhoton < 1. : continue
         if ev.Y1Y2 > 0. : continue
         if category =='_VBFcat':
            if ev.mJJ < 400. : continue
         else:
            if ev.mJJ < 300. : continue
         if ev.dEtaJJ < 3. : continue
#         if ev.dEtaGammaMesonCand > 2. : continue

      ## OPTIMIZED PHASE SPACE
      if(category =='_VBFcat'):
         if ev.photon_pt<PH_ptCut :  continue
         if ev.meson_pt<30 : continue # relaxed this

      ## OPTIMIZED PHASE SPACE
      if(category =='_GFcat'):
#         if ev.goodPhotons_pt[idxPh]<40 :  continue
#         if ev.goodMeson_pt[idxMeson]<40 : continue
         if ev.photon_pt<PH_ptCut :  continue
         if ev.meson_pt<MES_ptCut : continue
#         if ev.nGoodJets < 2: continue
#         if ev.nGoodPhotons > 1: continue
#         if ev.nPhotonsVeto>1: continue (COMMENT FOR NOW)

      # Fill histograms
      if (doSignal) :
#         wei = ev.w * ev.lumiIntegrated
         wei = ev.w_allSF * ev.lumiIntegrated
         if ((nameSig=='WH' or nameSig=='WHl') and (ev.mc==1011 or ev.mc==1012 or ev.mc==1021 or ev.mc==1022 or ev.mc==1031 or ev.mc==1032)): #W
            h.Fill( var, wei )
         if ((nameSig=='ZH' or nameSig=='ZHl') and (ev.mc==1013 or ev.mc==1023 or ev.mc==1014 or ev.mc==1024 or ev.mc==1033 or ev.mc==1034)): #Z
            h.Fill( var, wei )
         if ((nameSig=='TTH') and (ev.mc==1018 or ev.mc==1028 or ev.mc==1038)): #TTH
            h.Fill( var, wei )
         if ((nameSig=='VBFH' or nameSig=='VBFHlow') and (ev.mc==1010 or ev.mc==1020 or ev.mc==1030)): #VBF
            h.Fill( var, wei )
         if (nameSig=='ZinvH' and (ev.mc==1015 or ev.mc==1025 or ev.mc==1016 or ev.mc==1026)): #Zinv
            h.Fill( var, wei )
         if (nameSig=='ggH' and (ev.mc==1017 or ev.mc==1027 or ev.mc==1037)): #GF
            h.Fill( var, wei )
      else :
         wei = ev.w
         if ev.mc<0:   # only DATA
            h.Fill( var, wei )
#            if (var>100 and var < 170):
#               mystr=str(ev.run)+":"+str(ev.luminosityBlock)+":"+str(ev.event)+"\n"
#               with open("example_VBFcatlow_RhoCat_2018_mHwindow.txt", "a") as f:
#                  f.write(mystr)

      #Loop over tree done
#   if ev.mc>=0:
#      if(category == '_Zinvcat' or category == '_GFcat' or category == '_VBFcatlow'): year = '_12018'
#      print("lumi = ",lumis[year])
#      h.Scale(lumis[year])
   return h

def plot(item, nbin, low, high, doLog,category,mesonCat):

   h = getHisto(item, nbin, low, high, doLog,category,mesonCat,False)
   
   # Create canvas with pad
   c = ROOT.TCanvas("c", "", 600, 600)
   pad = ROOT.TPad("upper_pad", "", 0, 0, 1, 1)
#   gStyle.SetOptStat(0)
   if doLog : pad.SetLogy(1)
   pad.SetTickx(False)
   pad.SetTicky(False)
   pad.Draw()
   pad.cd()   
   
   h.GetXaxis().SetTitle("m_{#gamma#phi}^{H} [GeV]")
   h.GetYaxis().SetTitle("Events")
   h.GetYaxis().SetTitleOffset(1.1)   
   h.GetYaxis().SetLabelSize(0.04)
   h.GetYaxis().SetTitleSize(0.045)
   h.GetYaxis().ChangeLabel(1, -1, 0)
   h.Draw("hist")

   # Add legend
   if(category =='_Wcat' or category =='_Zcat'): deltay=0.25
   if(category =='_VBFcat'): deltay=0.15
#   legend = ROOT.TLegend(0.63-0.5, 0.88-deltay, 0.85-0.5, 0.88)
#   legend = ROOT.TLegend(0.63, 0.88-deltay, 0.85, 0.88)
   legend = ROOT.TLegend(0.13, 0.88-deltay, 0.35, 0.88)   
   legend.SetTextFont(42)
   legend.SetFillStyle(0)
   legend.SetBorderSize(0)
   legend.SetTextSize(0.04)
   legend.SetTextAlign(32)
 
   if h and category =='_Zcat' and mesonCat == "_PhiCat" and h.Integral()>0: legend.AddEntry(h, "ZH(#gamma#phi)", "f")
   if h and category =='_Zcat' and mesonCat == "_RhoCat" and h.Integral()>0: legend.AddEntry(h, "ZH(#gamma#phi)", "f")
   if h and category =='_Wcat' and mesonCat == "_PhiCat" and h.Integral()>0: legend.AddEntry(h, "WH(#gamma#phi)", "f")
   if h and category =='_Wcat' and mesonCat == "_RhoCat" and h.Integral()>0: legend.AddEntry(h, "WH(#rho#phi)", "f")   
   if h and category =='_VBFcat' and mesonCat == "_PhiCat" and h.Integral()>0: legend.AddEntry(h, "qqH(#gamma#phi)", "f")
   if h and category =='_VBFcat' and mesonCat == "_RhoCat" and h.Integral()>0: legend.AddEntry(h, "qqH(#rho#phi)", "f")
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
   text.DrawLatex(0.65, 0.93, "#sqrt{s} = 13 TeV,%0.2f fb^{-1}"% (lumis[year]))

   string = category+mesonCat

   if item==4:
      c.SaveAs("HCandMass"+string+".png")
      print("HCandmass.png")
      
if __name__ == "__main__":

   plot(4, 200, 0. , 200.,True,'_GFcat','_PhiCat') # HCandMass
#   plot(4, 200, 0. , 200.,True,'_Wcat','_PhiCat') # HCandMass
#   plot(4, 200, 0. , 200.,True,'_VBFcat','_PhiCat') # HCandMass
#   plot(4, 200, 0. , 200.,True,'_VBFcat','_RhoCat') # HCandMass
#   plot(4, 200, 0. , 200.,True,'_VBFcat','_PhiCat') # HCandMass
#   plot(4, 200, 0. , 200.,True,'_Wcat','_RhoCat') # HCandMass
