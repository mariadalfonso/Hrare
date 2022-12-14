import ROOT
import os


def loadTree(mytree, directory , category, mesonCat, year, doSignal ):

   if category=='_VBFcatlow' or category=='_GFcat':
      if doSignal:
         mytree.Add(directory+'outname_mc1020'+category+mesonCat+year+'.root') #VBFH -- rho
         mytree.Add(directory+'outname_mc1010'+category+mesonCat+year+'.root') #VBFH -- phi
         mytree.Add(directory+'outname_mc1027'+category+mesonCat+year+'.root') #GFH -- rho
         mytree.Add(directory+'outname_mc1017'+category+mesonCat+year+'.root') #GFH -- phi
      else:
         mytree.Add(directory+'outname_mc10'+category+mesonCat+year+'.root')  #GJets
         mytree.Add(directory+'outname_mc11'+category+mesonCat+year+'.root')  #GJets
         mytree.Add(directory+'outname_mc12'+category+mesonCat+year+'.root')  #GJets
         mytree.Add(directory+'outname_mc13'+category+mesonCat+year+'.root')  #GJets
         mytree.Add(directory+'outname_mc14'+category+mesonCat+year+'.root')  #GJets
         mytree.Add(directory+'outname_mc15'+category+mesonCat+year+'.root')  #VBF GJets
         mytree.Add(directory+'outname_mc-62'+category+mesonCat+year+'.root')  #DATA
         mytree.Add(directory+'outname_mc-63'+category+mesonCat+year+'.root')  #DATA
         mytree.Add(directory+'outname_mc-64'+category+mesonCat+year+'.root')  #DATA
   elif category=='_VBFcat':
      if directory.find("2018")!=-1:
         if doSignal:
            mytree.Add(directory+'outname_mc1020'+category+mesonCat+'_2018.root') #VBFH -- rho
            mytree.Add(directory+'outname_mc1010'+category+mesonCat+'_2018.root') #VBFH -- phi
         else:
            mytree.Add(directory+'outname_mc10'+category+mesonCat+'_2018.root')  #GJets
            mytree.Add(directory+'outname_mc11'+category+mesonCat+'_2018.root')  #GJets
            mytree.Add(directory+'outname_mc12'+category+mesonCat+'_2018.root')  #GJets
            mytree.Add(directory+'outname_mc13'+category+mesonCat+'_2018.root')  #GJets
            mytree.Add(directory+'outname_mc14'+category+mesonCat+'_2018.root')  #GJets
            mytree.Add(directory+'outname_mc15'+category+mesonCat+'_2018.root')  # VBF GJets
            mytree.Add(directory+'outname_mc45'+category+mesonCat+'_2018.root')  # W + GJets
            mytree.Add(directory+'outname_mc46'+category+mesonCat+'_2018.root')  # Z + GJets
            mytree.Add(directory+'outname_mc47'+category+mesonCat+'_2018.root')  # TT + GJets
            ##
            mytree.Add(directory+'outname_mc-31'+category+mesonCat+'_2018.root')  #DATA EG
            mytree.Add(directory+'outname_mc-32'+category+mesonCat+'_2018.root')  #DATA
            mytree.Add(directory+'outname_mc-33'+category+mesonCat+'_2018.root')  #DATA
            mytree.Add(directory+'outname_mc-34'+category+mesonCat+'_2018.root')  #DATA
            mytree.Add(directory+'outname_mc-62'+category+mesonCat+'_2018.root')  #DATA TAU
            mytree.Add(directory+'outname_mc-63'+category+mesonCat+'_2018.root')  #DATA
            mytree.Add(directory+'outname_mc-64'+category+mesonCat+'_2018.root')  #DATA
      if directory.find("2017")!=-1:
         if doSignal:
            mytree.Add(directory+'outname_mc1020'+category+mesonCat+'_2017.root') #VBFH -- rho
            mytree.Add(directory+'outname_mc1010'+category+mesonCat+'_2017.root') #VBFH -- phi
         else:
            mytree.Add(directory+'outname_mc10'+category+mesonCat+'_2017.root')  #GJets
            mytree.Add(directory+'outname_mc11'+category+mesonCat+'_2017.root')  #GJets
            mytree.Add(directory+'outname_mc12'+category+mesonCat+'_2017.root')  #GJets
            mytree.Add(directory+'outname_mc13'+category+mesonCat+'_2017.root')  #GJets
            mytree.Add(directory+'outname_mc14'+category+mesonCat+'_2017.root')  #GJets
            mytree.Add(directory+'outname_mc15'+category+mesonCat+'_2017.root')  # VBF GJets
            mytree.Add(directory+'outname_mc45'+category+mesonCat+'_2017.root')  # W + GJets
            mytree.Add(directory+'outname_mc46'+category+mesonCat+'_2017.root')  # Z + GJets
            mytree.Add(directory+'outname_mc47'+category+mesonCat+'_2017.root')  # TT + GJets
            ##
            mytree.Add(directory+'outname_mc-76'+category+mesonCat+'_2017.root')  #DATA PH
      if directory.find("_12016")!=-1:
         if doSignal:
            mytree.Add(directory+'outname_mc1020'+category+mesonCat+'_12016.root') #VBFH -- rho
            mytree.Add(directory+'outname_mc1010'+category+mesonCat+'_12016.root') #VBFH -- phi
         else:
            mytree.Add(directory+'outname_mc10'+category+mesonCat+'_12016.root')  #GJets
            mytree.Add(directory+'outname_mc11'+category+mesonCat+'_12016.root')  #GJets
            mytree.Add(directory+'outname_mc12'+category+mesonCat+'_12016.root')  #GJets
            mytree.Add(directory+'outname_mc13'+category+mesonCat+'_12016.root')  #GJets
            mytree.Add(directory+'outname_mc14'+category+mesonCat+'_12016.root')  #GJets
            mytree.Add(directory+'outname_mc15'+category+mesonCat+'_12016.root')  # VBF GJets
            #
            mytree.Add(directory+'outname_mc-81'+category+mesonCat+'_12016.root')  #DATA PH
            mytree.Add(directory+'outname_mc-82'+category+mesonCat+'_12016.root')  #DATA PH
            mytree.Add(directory+'outname_mc-83'+category+mesonCat+'_12016.root')  #DATA PH
            mytree.Add(directory+'outname_mc-84'+category+mesonCat+'_12016.root')  #DATA PH
            mytree.Add(directory+'outname_mc-85'+category+mesonCat+'_12016.root')  #DATA PH
            mytree.Add(directory+'outname_mc-86'+category+mesonCat+'_12016.root')  #DATA PH
            #
            mytree.Add(directory+'outname_mc45'+category+mesonCat+'_12016.root')  # W + GJets
            mytree.Add(directory+'outname_mc46'+category+mesonCat+'_12016.root')  # Z + GJets
            mytree.Add(directory+'outname_mc47'+category+mesonCat+'_12016.root')  # TT + GJets
   elif category=='_Zinvcat':
      if doSignal:
         mytree.Add(directory+'outname_mc1015'+category+mesonCat+year+'.root') #Z-H phi
         mytree.Add(directory+'outname_mc1025'+category+mesonCat+year+'.root') #Z-H rho
         mytree.Add(directory+'outname_mc1016'+category+mesonCat+year+'.root') #ggZ-H phi
         mytree.Add(directory+'outname_mc1026'+category+mesonCat+year+'.root') #ggZ-H rho
         mytree.Add(directory+'outname_mc1011'+category+mesonCat+year+'.root') #W+H phi
         mytree.Add(directory+'outname_mc1012'+category+mesonCat+year+'.root') #W-H
         mytree.Add(directory+'outname_mc1021'+category+mesonCat+year+'.root') #W+H rho
         mytree.Add(directory+'outname_mc1022'+category+mesonCat+year+'.root') #W-H
      else:
         mytree.Add(directory+'outname_mc10'+category+mesonCat+year+'.root')  #GJets
         mytree.Add(directory+'outname_mc11'+category+mesonCat+year+'.root')  #GJets
         mytree.Add(directory+'outname_mc12'+category+mesonCat+year+'.root')  #GJets
         mytree.Add(directory+'outname_mc13'+category+mesonCat+year+'.root')  #GJets
         mytree.Add(directory+'outname_mc14'+category+mesonCat+year+'.root')  #GJets
         mytree.Add(directory+'outname_mc1'+category+mesonCat+year+'.root')   #ZG
         #
         mytree.Add(directory+'outname_mc2'+category+mesonCat+year+'.root')   #WG
         mytree.Add(directory+'outname_mc31'+category+mesonCat+year+'.root')  #W0Jet
         mytree.Add(directory+'outname_mc32'+category+mesonCat+year+'.root')  #W1Jet
         mytree.Add(directory+'outname_mc33'+category+mesonCat+year+'.root')  #W2Jet
         mytree.Add(directory+'outname_mc37'+category+mesonCat+year+'.root')  # Z1nn
         mytree.Add(directory+'outname_mc38'+category+mesonCat+year+'.root')  #Z1nn
         mytree.Add(directory+'outname_mc39'+category+mesonCat+year+'.root')  #Z1nn
         mytree.Add(directory+'outname_mc40'+category+mesonCat+year+'.root')  #Z1nn
         mytree.Add(directory+'outname_mc41'+category+mesonCat+year+'.root')  # Z2nn
         mytree.Add(directory+'outname_mc42'+category+mesonCat+year+'.root')  #Z2nn
         mytree.Add(directory+'outname_mc43'+category+mesonCat+year+'.root')  #Z2nn
         mytree.Add(directory+'outname_mc44'+category+mesonCat+year+'.root')  #Z2nn
         #      mytree.Add(directory+'outname_mc20'+category+mesonCat+year+'.root')  #PtJet too low weight
         #      mytree.Add(directory+'outname_mc21'+category+mesonCat+year+'.root')  #PtJet too low weight
         #      mytree.Add(directory+'outname_mc22'+category+mesonCat+year+'.root')  #PtJet too low weight
         mytree.Add(directory+'outname_mc23'+category+mesonCat+year+'.root')  #PtJet
         mytree.Add(directory+'outname_mc24'+category+mesonCat+year+'.root')  #PtJet
         mytree.Add(directory+'outname_mc25'+category+mesonCat+year+'.root')  #PtJet
         #
         mytree.Add(directory+'outname_mc-62'+category+mesonCat+year+'.root')  #DATA
         mytree.Add(directory+'outname_mc-63'+category+mesonCat+year+'.root')  #DATA
         mytree.Add(directory+'outname_mc-64'+category+mesonCat+year+'.root')  #DATA
   elif category=='_Wcat' or category=='_Zcat':
      loalTime = '_all'
      if directory.find("2018")!=-1:
         localTime = '_2018'
         ##
         mytree.Add(directory+'outname_mc-1'+category+mesonCat+'_2018.root')  #DATA oneMu
         mytree.Add(directory+'outname_mc-2'+category+mesonCat+'_2018.root')  #DATA
         mytree.Add(directory+'outname_mc-3'+category+mesonCat+'_2018.root')  #DATA
         mytree.Add(directory+'outname_mc-4'+category+mesonCat+'_2018.root')  #DATA
         #
         mytree.Add(directory+'outname_mc-11'+category+mesonCat+'_2018.root')  #DATA diMu
         mytree.Add(directory+'outname_mc-12'+category+mesonCat+'_2018.root')  #DATA
         mytree.Add(directory+'outname_mc-13'+category+mesonCat+'_2018.root')  #DATA
         mytree.Add(directory+'outname_mc-14'+category+mesonCat+'_2018.root')  #DATA
         #
         mytree.Add(directory+'outname_mc-21'+category+mesonCat+'_2018.root')  #DATA diMu
         mytree.Add(directory+'outname_mc-22'+category+mesonCat+'_2018.root')  #DATA
         mytree.Add(directory+'outname_mc-23'+category+mesonCat+'_2018.root')  #DATA
         mytree.Add(directory+'outname_mc-24'+category+mesonCat+'_2018.root')  #DATA
         #
         mytree.Add(directory+'outname_mc-31'+category+mesonCat+'_2018.root')  #DATA EG
         mytree.Add(directory+'outname_mc-32'+category+mesonCat+'_2018.root')  #DATA
         mytree.Add(directory+'outname_mc-33'+category+mesonCat+'_2018.root')  #DATA
         mytree.Add(directory+'outname_mc-34'+category+mesonCat+'_2018.root')  #DATA
         #
         mytree.Add(directory+'outname_mc-21'+category+mesonCat+'_2018.root')  #DATA muEG
         mytree.Add(directory+'outname_mc-22'+category+mesonCat+'_2018.root')  #DATA
         mytree.Add(directory+'outname_mc-23'+category+mesonCat+'_2018.root')  #DATA
         mytree.Add(directory+'outname_mc-24'+category+mesonCat+'_2018.root')  #DATA
      if directory.find("2017")!=-1:
         localTime = '_2017'
         ##
         mytree.Add(directory+'outname_mc-2'+category+mesonCat+'_2017.root')  #DATA oneMu
         mytree.Add(directory+'outname_mc-3'+category+mesonCat+'_2017.root')  #DATA
         mytree.Add(directory+'outname_mc-4'+category+mesonCat+'_2017.root')  #DATA
         mytree.Add(directory+'outname_mc-5'+category+mesonCat+'_2017.root')  #DATA
         mytree.Add(directory+'outname_mc-6'+category+mesonCat+'_2017.root')  #DATA
         #
         mytree.Add(directory+'outname_mc-12'+category+mesonCat+'_2017.root')  #DATA diMu
         mytree.Add(directory+'outname_mc-13'+category+mesonCat+'_2017.root')  #DATA
         mytree.Add(directory+'outname_mc-14'+category+mesonCat+'_2017.root')  #DATA
         mytree.Add(directory+'outname_mc-15'+category+mesonCat+'_2017.root')  #DATA
         mytree.Add(directory+'outname_mc-16'+category+mesonCat+'_2017.root')  #DATA
         #
         mytree.Add(directory+'outname_mc-22'+category+mesonCat+'_2017.root')  #DATA diMu
         mytree.Add(directory+'outname_mc-23'+category+mesonCat+'_2017.root')  #DATA
         mytree.Add(directory+'outname_mc-24'+category+mesonCat+'_2017.root')  #DATA
         mytree.Add(directory+'outname_mc-25'+category+mesonCat+'_2017.root')  #DATA
         mytree.Add(directory+'outname_mc-26'+category+mesonCat+'_2017.root')  #DATA
         #
         mytree.Add(directory+'outname_mc-42'+category+mesonCat+'_2017.root')  #DATA diEG
         mytree.Add(directory+'outname_mc-43'+category+mesonCat+'_2017.root')  #DATA
         mytree.Add(directory+'outname_mc-44'+category+mesonCat+'_2017.root')  #DATA
         mytree.Add(directory+'outname_mc-45'+category+mesonCat+'_2017.root')  #DATA
         mytree.Add(directory+'outname_mc-46'+category+mesonCat+'_2017.root')  #DATA
         #
         mytree.Add(directory+'outname_mc-52'+category+mesonCat+'_2017.root')  #DATA oneEL
         mytree.Add(directory+'outname_mc-53'+category+mesonCat+'_2017.root')  #DATA
         mytree.Add(directory+'outname_mc-54'+category+mesonCat+'_2017.root')  #DATA
         mytree.Add(directory+'outname_mc-55'+category+mesonCat+'_2017.root')  #DATA
         mytree.Add(directory+'outname_mc-56'+category+mesonCat+'_2017.root')  #DATA
      if directory.find("12016")!=-1:
         localTime = '_12016'
         mytree.Add(directory+'outname_mc-1'+category+mesonCat+'_12016.root')  #DATA oneMu
         mytree.Add(directory+'outname_mc-2'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-3'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-4'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-5'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-6'+category+mesonCat+'_12016.root')  #DATA
         #
         mytree.Add(directory+'outname_mc-11'+category+mesonCat+'_12016.root')  #DATA diMu
         mytree.Add(directory+'outname_mc-12'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-13'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-14'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-15'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-16'+category+mesonCat+'_12016.root')  #DATA
         #
         mytree.Add(directory+'outname_mc-21'+category+mesonCat+'_12016.root')  #DATA diMu
         mytree.Add(directory+'outname_mc-22'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-23'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-24'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-25'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-26'+category+mesonCat+'_12016.root')  #DATA
         #
         mytree.Add(directory+'outname_mc-41'+category+mesonCat+'_12016.root')  #DATA diEG
         mytree.Add(directory+'outname_mc-42'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-43'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-44'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-45'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-46'+category+mesonCat+'_12016.root')  #DATA
         #
         mytree.Add(directory+'outname_mc-51'+category+mesonCat+'_12016.root')  #DATA oneEL
         mytree.Add(directory+'outname_mc-52'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-53'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-54'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-55'+category+mesonCat+'_12016.root')  #DATA
         mytree.Add(directory+'outname_mc-56'+category+mesonCat+'_12016.root')  #DATA
      if directory.find("22016")!=-1:
         localTime = '_22016'
         mytree.Add(directory+'outname_mc-6'+category+mesonCat+'_22016.root')  #DATA oneMu
         mytree.Add(directory+'outname_mc-7'+category+mesonCat+'_22016.root')  #DATA
         mytree.Add(directory+'outname_mc-8'+category+mesonCat+'_22016.root')  #DATA
         mytree.Add(directory+'outname_mc-16'+category+mesonCat+'_22016.root')  #DATA diMu
         mytree.Add(directory+'outname_mc-17'+category+mesonCat+'_22016.root')  #DATA
         mytree.Add(directory+'outname_mc-18'+category+mesonCat+'_22016.root')  #DATA
         mytree.Add(directory+'outname_mc-26'+category+mesonCat+'_22016.root')  #DATA MuEG
         mytree.Add(directory+'outname_mc-27'+category+mesonCat+'_22016.root')  #DATA
         mytree.Add(directory+'outname_mc-28'+category+mesonCat+'_22016.root')  #DATA
         mytree.Add(directory+'outname_mc-46'+category+mesonCat+'_22016.root')  #DATA diEG
         mytree.Add(directory+'outname_mc-47'+category+mesonCat+'_22016.root')  #DATA
         mytree.Add(directory+'outname_mc-48'+category+mesonCat+'_22016.root')  #DATA
         mytree.Add(directory+'outname_mc-56'+category+mesonCat+'_22016.root')  #DATA oneEL
         mytree.Add(directory+'outname_mc-57'+category+mesonCat+'_22016.root')  #DATA
         mytree.Add(directory+'outname_mc-58'+category+mesonCat+'_22016.root')  #DATA
      mytree.Add(directory+'outname_mc1011'+category+mesonCat+localTime+'.root') #W+H phi
      mytree.Add(directory+'outname_mc1012'+category+mesonCat+localTime+'.root') #W-H
      mytree.Add(directory+'outname_mc1021'+category+mesonCat+localTime+'.root') #W+H rho
      mytree.Add(directory+'outname_mc1022'+category+mesonCat+localTime+'.root') #W-H
      mytree.Add(directory+'outname_mc1013'+category+mesonCat+localTime+'.root') #Z-H
      mytree.Add(directory+'outname_mc1014'+category+mesonCat+localTime+'.root') #ggZ-H
      mytree.Add(directory+'outname_mc1023'+category+mesonCat+localTime+'.root') #Z-H
      mytree.Add(directory+'outname_mc1024'+category+mesonCat+localTime+'.root') #ggZ-H
#      mytree.Add(directory+'outname_mc1'+category+mesonCat+localTime+'.root')  #ZG
#      mytree.Add(directory+'outname_mc2'+category+mesonCat+localTime+'.root')  #WG
#      mytree.Add(directory+'outname_mc31'+category+mesonCat+localTime+'.root')  #W0Jet
#      mytree.Add(directory+'outname_mc32'+category+mesonCat+localTime+'.root')  #W1Jet
#      mytree.Add(directory+'outname_mc33'+category+mesonCat+localTime+'.root')  #W2Jet
#      mytree.Add(directory+'outname_mc34'+category+mesonCat+localTime+'.root')  #DY0Jet
#      mytree.Add(directory+'outname_mc35'+category+mesonCat+localTime+'.root')  #DY1Jet
#      mytree.Add(directory+'outname_mc36'+category+mesonCat+localTime+'.root')  #DY2Jet
      ##   #mytree.Add(directory+'outname_mc0'+category+mesonCat+year+'.root')  #DY
      ##   #mytree.Add(directory+'outname_mc3'+category+mesonCat+year+'.root')  #Wjets
      ##   mytree.Add(directory+'outname_mc4'+category+mesonCat+year+'.root')  #tt2l
      ##   mytree.Add(directory+'outname_mc5'+category+mesonCat+year+'.root')  #tt1l
   
   else:
      print("ERROR: category not specified")

   resetTree(mytree,category)

   return mytree

def resetTree(mytree,category):

   mytree.SetBranchStatus('*',0)
   mytree.SetBranchStatus('w',1)
   mytree.SetBranchStatus('lumiIntegrated',1)
   mytree.SetBranchStatus('mc',1)

   mytree.SetBranchStatus('HCandMass',1)
   mytree.SetBranchStatus('photon_pt',1)
   mytree.SetBranchStatus('meson_pt',1)
   mytree.SetBranchStatus('index_pair',1)

   if(category =='_GFcat' or category =='_Zinvcat' or category =='_VBFcat' or category =='_VBFcatlow'):
      mytree.SetBranchStatus('MVAdisc',1)

   if(category =='_Wcat'):
      mytree.SetBranchStatus('Z_veto',1)

   if(category =='_VBFcat' or category =='_VBFcatlow'):
      mytree.SetBranchStatus('jet1Pt',1)
      mytree.SetBranchStatus('jet2Pt',1)
      mytree.SetBranchStatus('deltaJetMeson',1)
      mytree.SetBranchStatus('deltaJetPhoton',1)
      mytree.SetBranchStatus('Y1Y2',1)
      mytree.SetBranchStatus('mJJ',1)
      mytree.SetBranchStatus('dEtaJJ',1)

   if(category =='_Zinvcat' or category =='_Wcat'):
      mytree.SetBranchStatus('DeepMETResolutionTune_pt',1)

   if(category =='_Wcat'):
      mytree.SetBranchStatus('V_mass',1)

   if(category =='_Zinvcat'):
      mytree.SetBranchStatus('dPhiGammaMET',1)
      mytree.SetBranchStatus('dPhiMesonMET',1)
      mytree.SetBranchStatus('ptRatioMEThiggs',1)
      mytree.SetBranchStatus('nbtag',1)

   return mytree
