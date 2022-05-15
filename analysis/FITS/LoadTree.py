import ROOT
import os

def loadTree(mytree, directory , category, mesonCat, year ):

   if category=='_VBFcatlow' or category=='_GFcat':
      mytree.Add(directory+'outname_mc1020'+category+mesonCat+year+'.root') #VBFH -- rho
      mytree.Add(directory+'outname_mc1010'+category+mesonCat+year+'.root') #VBFH -- phi
      mytree.Add(directory+'outname_mc1027'+category+mesonCat+year+'.root') #GFH -- rho
      mytree.Add(directory+'outname_mc1017'+category+mesonCat+year+'.root') #GFH -- phi
      mytree.Add(directory+'outname_mc6'+category+mesonCat+year+'.root')  #GJets
      mytree.Add(directory+'outname_mc7'+category+mesonCat+year+'.root')  #GJets
      mytree.Add(directory+'outname_mc8'+category+mesonCat+year+'.root')  #GJets
      mytree.Add(directory+'outname_mc9'+category+mesonCat+year+'.root')  #GJets
      mytree.Add(directory+'outname_mc-62'+category+mesonCat+year+'.root')  #DATA
      mytree.Add(directory+'outname_mc-63'+category+mesonCat+year+'.root')  #DATA
      mytree.Add(directory+'outname_mc-64'+category+mesonCat+year+'.root')  #DATA
   elif category=='_VBFcat':
      mytree.Add(directory+'outname_mc1020'+category+mesonCat+year+'.root') #VBFH -- rho
      mytree.Add(directory+'outname_mc1010'+category+mesonCat+year+'.root') #VBFH -- phi
      mytree.Add(directory+'outname_mc6'+category+mesonCat+year+'.root')  #GJets
      mytree.Add(directory+'outname_mc7'+category+mesonCat+year+'.root')  #GJets
      mytree.Add(directory+'outname_mc8'+category+mesonCat+year+'.root')  #GJets
      mytree.Add(directory+'outname_mc9'+category+mesonCat+year+'.root')  #GJets
      ##   mytree.Add(directory5+'outname_mc20'+category+mesonCat+year+'.root')  #PtJet too low weight
      ##   mytree.Add(directory5+'outname_mc21'+category+mesonCat+year+'.root')  #PtJet too low weight
      ##   mytree.Add(directory5+'outname_mc22'+category+mesonCat+year+'.root')  #PtJet too low weight
      #   mytree.Add(directory5+'outname_mc23'+category+mesonCat+year+'.root')  #PtJet
      #   mytree.Add(directory5+'outname_mc24'+category+mesonCat+year+'.root')  #PtJet
      #   mytree.Add(directory5+'outname_mc25'+category+mesonCat+year+'.root')  #PtJet
      mytree.Add(directory+'outname_mc-31'+category+mesonCat+year+'.root')  #DATA EG
      mytree.Add(directory+'outname_mc-32'+category+mesonCat+year+'.root')  #DATA EG
      mytree.Add(directory+'outname_mc-33'+category+mesonCat+year+'.root')  #DATA EG
      mytree.Add(directory+'outname_mc-34'+category+mesonCat+year+'.root')  #DATA EG
   elif category=='_Zinvcat':
      mytree.Add(directory+'outname_mc1015'+category+mesonCat+year+'.root') #Z-H phi
      mytree.Add(directory+'outname_mc1025'+category+mesonCat+year+'.root') #Z-H rho
      mytree.Add(directory+'outname_mc1016'+category+mesonCat+year+'.root') #ggZ-H phi
      mytree.Add(directory+'outname_mc1026'+category+mesonCat+year+'.root') #ggZ-H rho
      mytree.Add(directory+'outname_mc1011'+category+mesonCat+year+'.root') #W+H phi
      mytree.Add(directory+'outname_mc1012'+category+mesonCat+year+'.root') #W-H
      mytree.Add(directory+'outname_mc1021'+category+mesonCat+year+'.root') #W+H rho
      mytree.Add(directory+'outname_mc1022'+category+mesonCat+year+'.root') #W-H
      mytree.Add(directory+'outname_mc1'+category+mesonCat+year+'.root')  #ZG
      mytree.Add(directory+'outname_mc2'+category+mesonCat+year+'.root')  #WG
      mytree.Add(directory+'outname_mc6'+category+mesonCat+year+'.root')  #GJets
      mytree.Add(directory+'outname_mc7'+category+mesonCat+year+'.root')  #GJets
      mytree.Add(directory+'outname_mc8'+category+mesonCat+year+'.root')  #GJets
      mytree.Add(directory+'outname_mc9'+category+mesonCat+year+'.root')  #GJets
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
      mytree.Add(directory+'outname_mc-62'+category+mesonCat+year+'.root')  #DATA
      mytree.Add(directory+'outname_mc-63'+category+mesonCat+year+'.root')  #DATA
      mytree.Add(directory+'outname_mc-64'+category+mesonCat+year+'.root')  #DATA
   elif category=='_Zcat':
      mytree.Add(directory+'outname_mc1013'+category+mesonCat+year+'.root') #Z-H phi
      mytree.Add(directory+'outname_mc1023'+category+mesonCat+year+'.root') #Z-H rho
      mytree.Add(directory+'outname_mc1014'+category+mesonCat+year+'.root') #ggZ-H phi
      mytree.Add(directory+'outname_mc1024'+category+mesonCat+year+'.root') #ggZ-H rho
      mytree.Add(directory+'outname_mc1'+category+mesonCat+year+'.root')  #ZG   
      mytree.Add(directory+'outname_mc34'+category+mesonCat+year+'.root')  #DY0Jet
      mytree.Add(directory+'outname_mc35'+category+mesonCat+year+'.root')  #DY1Jet
      mytree.Add(directory+'outname_mc36'+category+mesonCat+year+'.root')  #DY2Jet
      mytree.Add(directory+'outname_mc-1'+category+mesonCat+year+'.root')  #DATA singleMu
      mytree.Add(directory+'outname_mc-2'+category+mesonCat+year+'.root')  #DATA
      mytree.Add(directory+'outname_mc-3'+category+mesonCat+year+'.root')  #DATA
      mytree.Add(directory+'outname_mc-4'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017'): mytree.Add(directory+'outname_mc-5'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017'): mytree.Add(directory+'outname_mc-6'+category+mesonCat+year+'.root')  #DATA
      
      mytree.Add(directory+'outname_mc-11'+category+mesonCat+year+'.root')  #DATA diMu
      mytree.Add(directory+'outname_mc-12'+category+mesonCat+year+'.root')  #DATA
      mytree.Add(directory+'outname_mc-13'+category+mesonCat+year+'.root')  #DATA
      mytree.Add(directory+'outname_mc-14'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017'): mytree.Add(directory+'outname_mc-15'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017'): mytree.Add(directory+'outname_mc-16'+category+mesonCat+year+'.root')  #DATA
      
      mytree.Add(directory+'outname_mc-21'+category+mesonCat+year+'.root')  #DATA muEG
      mytree.Add(directory+'outname_mc-22'+category+mesonCat+year+'.root')  #DATA
      mytree.Add(directory+'outname_mc-23'+category+mesonCat+year+'.root')  #DATA
      mytree.Add(directory+'outname_mc-24'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017'): mytree.Add(directory+'outname_mc-25'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017'): mytree.Add(directory+'outname_mc-26'+category+mesonCat+year+'.root')  #DATA
      
      mytree.Add(directory+'outname_mc-31'+category+mesonCat+year+'.root')  #DATA EG
      mytree.Add(directory+'outname_mc-32'+category+mesonCat+year+'.root')  #DATA
      mytree.Add(directory+'outname_mc-33'+category+mesonCat+year+'.root')  #DATA
      mytree.Add(directory+'outname_mc-34'+category+mesonCat+year+'.root')  #DATA
      
   elif category=='_Wcat':
      mytree.Add(directory+'outname_mc1011'+category+mesonCat+year+'.root') #W+H phi
      mytree.Add(directory+'outname_mc1012'+category+mesonCat+year+'.root') #W-H
      mytree.Add(directory+'outname_mc1013'+category+mesonCat+year+'.root') #Z-H
      mytree.Add(directory+'outname_mc1014'+category+mesonCat+year+'.root') #ggZ-H
      mytree.Add(directory+'outname_mc1021'+category+mesonCat+year+'.root') #W+H rho
      mytree.Add(directory+'outname_mc1022'+category+mesonCat+year+'.root') #W-H      
      mytree.Add(directory+'outname_mc1023'+category+mesonCat+year+'.root') #Z-H
      mytree.Add(directory+'outname_mc1024'+category+mesonCat+year+'.root') #ggZ-H
      mytree.Add(directory+'outname_mc1'+category+mesonCat+year+'.root')  #ZG
      mytree.Add(directory+'outname_mc2'+category+mesonCat+year+'.root')  #WG
      mytree.Add(directory+'outname_mc31'+category+mesonCat+year+'.root')  #W0Jet
      mytree.Add(directory+'outname_mc32'+category+mesonCat+year+'.root')  #W1Jet
      mytree.Add(directory+'outname_mc33'+category+mesonCat+year+'.root')  #W2Jet
      mytree.Add(directory+'outname_mc34'+category+mesonCat+year+'.root')  #DY0Jet
      mytree.Add(directory+'outname_mc35'+category+mesonCat+year+'.root')  #DY1Jet
      mytree.Add(directory+'outname_mc36'+category+mesonCat+year+'.root')  #DY2Jet
      ##   #mytree.Add(directory+'outname_mc0'+category+mesonCat+year+'.root')  #DY
      ##   #mytree.Add(directory+'outname_mc3'+category+mesonCat+year+'.root')  #Wjets
      ##   mytree.Add(directory+'outname_mc4'+category+mesonCat+year+'.root')  #tt2l
      ##   mytree.Add(directory+'outname_mc5'+category+mesonCat+year+'.root')  #tt1l
      
      if (year == '_2018'): mytree.Add(directory+'outname_mc-1'+category+mesonCat+year+'.root')  #DATA singleMu
      if (year == '_2017' or year == '_2018'): mytree.Add(directory+'outname_mc-2'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017' or year == '_2018'): mytree.Add(directory+'outname_mc-3'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017' or year == '_2018'): mytree.Add(directory+'outname_mc-4'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017'): mytree.Add(directory+'outname_mc-5'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017' or year == '_2016'): mytree.Add(directory+'outname_mc-6'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2016'): mytree.Add(directory+'outname_mc-7'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2016'): mytree.Add(directory+'outname_mc-8'+category+mesonCat+year+'.root')  #DATA
      
      if (year == '_2018'): mytree.Add(directory+'outname_mc-11'+category+mesonCat+year+'.root')  #DATA diMu
      if (year == '_2017' or year == '_2018'): mytree.Add(directory+'outname_mc-12'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017' or year == '_2018'): mytree.Add(directory+'outname_mc-13'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017' or year == '_2018'): mytree.Add(directory+'outname_mc-14'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017'): mytree.Add(directory+'outname_mc-15'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017' or year == '_2016'): mytree.Add(directory+'outname_mc-16'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2016'): mytree.Add(directory+'outname_mc-17'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2016'): mytree.Add(directory+'outname_mc-18'+category+mesonCat+year+'.root')  #DATA
      
      if (year == '_2018'): mytree.Add(directory+'outname_mc-21'+category+mesonCat+year+'.root')  #DATA muEG
      if (year == '_2017' or year == '_2018'): mytree.Add(directory+'outname_mc-22'+category+mesonCat+year+'.root') 
      if (year == '_2017' or year == '_2018'): mytree.Add(directory+'outname_mc-23'+category+mesonCat+year+'.root')
      if (year == '_2017' or year == '_2018'): mytree.Add(directory+'outname_mc-24'+category+mesonCat+year+'.root')
      if (year == '_2017'): mytree.Add(directory+'outname_mc-25'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017' or year == '_2016'): mytree.Add(directory+'outname_mc-26'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2016'): mytree.Add(directory+'outname_mc-27'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2016'): mytree.Add(directory+'outname_mc-28'+category+mesonCat+year+'.root')  #DATA   
      
      if (year == '_2018'): mytree.Add(directory+'outname_mc-31'+category+mesonCat+year+'.root')  #DATA EG
      if (year == '_2018'): mytree.Add(directory+'outname_mc-32'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2018'): mytree.Add(directory+'outname_mc-33'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2018'): mytree.Add(directory+'outname_mc-34'+category+mesonCat+year+'.root')  #DATA
      
      if (year == '_2017'): mytree.Add(directory+'outname_mc-42'+category+mesonCat+year+'.root')  #DATA B
      if (year == '_2017'): mytree.Add(directory+'outname_mc-43'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017'): mytree.Add(directory+'outname_mc-44'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017'): mytree.Add(directory+'outname_mc-45'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017' or year == '_2016'): mytree.Add(directory+'outname_mc-46'+category+mesonCat+year+'.root')  #DATA      
      if (year == '_2016'):  mytree.Add(directory+'outname_mc-47'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2016'):  mytree.Add(directory+'outname_mc-48'+category+mesonCat+year+'.root')  #DATA      
      
      if (year == '_2017'): mytree.Add(directory+'outname_mc-52'+category+mesonCat+year+'.root')  #DATA B 
      if (year == '_2017'): mytree.Add(directory+'outname_mc-53'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017'): mytree.Add(directory+'outname_mc-54'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017'): mytree.Add(directory+'outname_mc-55'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2017' or year == '_2016'): mytree.Add(directory+'outname_mc-56'+category+mesonCat+year+'.root')  #DATA      
      if (year == '_2016'):  mytree.Add(directory+'outname_mc-57'+category+mesonCat+year+'.root')  #DATA
      if (year == '_2016'):  mytree.Add(directory+'outname_mc-58'+category+mesonCat+year+'.root')  #DATA      
   
   else:
      print("ERROR: category not specified")

   return mytree

