import ROOT
import os
import math

'''
-rw-r--r-- 1 mariadlf mariadlf   62892430 May 20 11:06 /work/submit/mariadlf/Hrare_JPsiCC/MARCH2025/snapshotJpsiCC_10_2018.root
-rw-r--r-- 1 mariadlf mariadlf  189176550 May 20 11:22 /work/submit/mariadlf/Hrare_JPsiCC/MARCH2025/snapshotJpsiCC_11_2018.root
-rw-r--r-- 1 mariadlf mariadlf   23416648 May 20 11:29 /work/submit/mariadlf/Hrare_JPsiCC/MARCH2025/snapshotJpsiCC_1000_2018.root
'''

def loadTree(mytree, directory , category, year, doSignal ):
   ROOT.TH1.AddDirectory(False)
   filenames = []
   
   if category in ['_VBFcatlow', '_GFcat']:
      filenames += [
         f"{directory}snapshotJpsiCC_1000{year}.root",
         f"{directory}snapshotJpsiCC_10{year}.root",
         f"{directory}snapshotJpsiCC_11{year}.root"
      ]
      
   for fname in filenames:
      f = ROOT.TFile.Open(fname)
      if f:
         f.Close()  # Just ensure it's not loading any objects
         mytree.Add(fname)

   return mytree            
            
def loadTreeOLD(mytree, directory , category, year, doSignal ):

   print(loadTree)
   
   doLoadMC = True
   
   if category=='_VBFcatlow' or category=='_GFcat':
      mytree.Add(directory+'snapshotJpsiCC_1000'+year+'.root') #signal ggH
      mytree.Add(directory+'snapshotJpsiCC_10'+year+'.root') #signal ggH
      mytree.Add(directory+'snapshotJpsiCC_11'+year+'.root') #signal ggH         

   # note: there is another call for the signal
#   resetTree(mytree,category,doSignal)

   return mytree

def resetTree(mytree,category,doSignal):

   print(resetTree)
   
   mytree.SetBranchStatus('*',0)
   mytree.SetBranchStatus('w',1)
#   mytree.SetBranchStatus('w_allSF',1)
   mytree.SetBranchStatus('lumiIntegrated',1)
   mytree.SetBranchStatus('mc',1)
   mytree.SetBranchStatus('run',1)
   mytree.SetBranchStatus('luminosityBlock',1)
   mytree.SetBranchStatus('event',1)

   mytree.SetBranchStatus('massHiggsCorr',1)

   return mytree

def checkNan(var,mytree):

  print('cheching NaN in ',var)

  # Loop over entries like a normal TTree
  
  for i in range(mytree.GetEntries()):
    mytree.GetEntry(i)
    value = getattr(mytree, var)  # Access the branch dynamically
    
  # If value is iterable (e.g. vector<float>)
  
    if hasattr(value, "__len__") and not isinstance(value, str):
        for j, v in enumerate(value):
            if isinstance(v, float) and math.isnan(v):
                print(f"Found NaN in {var}[{j}] at entry {i}: {v}")
    else:
        if isinstance(value, float) and math.isnan(value):
            print(f"Found NaN in scalar {var} at entry {i}: {value}")

