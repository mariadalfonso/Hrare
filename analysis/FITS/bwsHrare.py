#!env python
import sys, os
import re
#from array import array
import math
from optparse import OptionParser,OptionGroup

parser= OptionParser()

parser.add_option("","--inputFileSig",type='string',help="Input ROOT file. [%default]", default="WS/Signal_Wcat__RhoCat_2018_workspace.root")
parser.add_option("","--inputFileBKG",type='string',help="Input ROOT file bkg model. [%default]", default="WS/Bkg_Wcat__RhoCat_2018_workspace.root")

parser.add_option("-c","--whichCat",type='string',help="Which category (Wcat, Zcat Zinvcatm, VBFcat)", default="Wcat")
parser.add_option("-o","--output",type='string',help="Output ROOT file. [%default]", default="workspace_STAT_Rho_2018.root")
parser.add_option("-d","--datCardName",type='string',help="Output txt file. [%default]", default="datacard_STAT_Rho_2018.txt")

opts, args = parser.parse_args()

sys.argv=[]
import ROOT
ROOT.gROOT.SetBatch()

############ CONFIGURABLES ###########

lumis={
    12016: 19.52, #APV
    22016: 16.80, #postVFP
    2016: 35.9,
    2017: 41.5,
    2018: 60.0,
}

BkgPdf={
    'Wcat': 'exp1',
    'Zcat': 'exp1',
    'VBFcat': 'bxg',
    'Zinvcat': 'bxg',
    'VBFcatlow': 'bxg',
    'GFcat': 'bxg',
}

SigPdf={
    'Wcat': 'crystal_ball',
    'Zcat': 'crystal_ball',
    'VBFcat': 'crystal_ball',
    'Zinvcat': 'crystal_ball',
    'VBFcatlow': 'crystal_ball',
    'GFcat': 'crystal_ball',
}

ENUM={
    'ggH': 0,
    'VBFH': -1,
    'WH': -2,
    'ZH': -3,
    'ZinvH': -4,
}


print(opts.whichCat)

if opts.whichCat=='GFcat':
    sigAll = ['ggH','VBFH']
    mcAll = ['ggH','VBFH','bkg']
    category = ['GFcat']

if opts.whichCat=='Wcat':
    sigAll = ['WH','ZH']
    mcAll = ['WH','ZH','bkg']
    category = ['Wcat']

if opts.whichCat=='Zcat':
    sigAll = ['ZH']
    mcAll = ['ZH','bkg']
    category = ['Zcat']

if opts.whichCat=='VBFcat':
    sigAll = ['VBFH']
    mcAll = ['VBFH','bkg']
    category = ['VBFcat']

if opts.whichCat=='VBFcatlow':
    sigAll = ['VBFH']
    mcAll = ['VBFH','bkg']
    category = ['VBFcatlow']

if opts.whichCat=='Zinvcat':
    sigAll = ['WH','ZinvH']
    mcAll = ['WH','ZinvH','bkg']
    category = ['Zinvcat']

################### OPEN OUTPUT ############
w = ROOT.RooWorkspace("w","w")

############ DATACARD ###########
#datName = "cms_datacard_ws"
#datName += ".txt"
datName = opts.datCardName

numChannel = len(mcAll)-1
print("-> Opening datacard",datName)
datacard=open(datName,"w")
datacard.write("-------------------------------------\n")
datacard.write("imax "+str(len(category))+" number of channels\n")
datacard.write("jmax "+ str(numChannel)+" number of background minus 1\n")
datacard.write("kmax * number of nuisance parameters\n")
datacard.write("-------------------------------------\n")

########################## IMPORT DATA #############
w.factory("mh[100,250]"); # RooRealVar
mh=w.var("mh")
arglist_obs = ROOT.RooArgList(mh)
argset_obs = ROOT.RooArgSet(mh)

################# Import SIGNAL/BKG CONTRIBUTIONS ##############

fSigIn = ROOT.TFile.Open(opts.inputFileSig,"READ")
if fSigIn == None: 
    print("ERROR: file",opts.inputFileSig,"doesn't exist")
    exit(1)    

fBkgIn = ROOT.TFile.Open(opts.inputFileBKG,"READ")
if fBkgIn == None: 
    print("ERROR: file",opts.inputFileBKG,"doesn't exist")
    exit(1)        
    
for cat in category:
    for proc in mcAll:
        
        if proc=='bkg':
            if opts.inputFileBKG != "":
                wInput=fBkgIn.Get("w")
                name = BkgPdf[cat]+"_"+cat+"_"+proc
                nameNorm = name+"_norm"
        else:
            if opts.inputFileSig != "":
                wInput=fSigIn.Get("w")
                name = SigPdf[cat]+"_"+cat+"_"+proc
                nameNorm = name+"_norm"
            
        print("proc=",proc," cat=",cat," name=",name)
        func = wInput.pdf(name)
        if func == None: raise IOError("Unable to get func" + name)
        getattr(w,'import')(func)
        funcNorm = wInput.var(nameNorm)
        if funcNorm == None: raise IOError("Unable to get func normalization " + nameNorm)
        getattr(w,'import')(funcNorm)

        datacard.write("shapes")
        datacard.write("\t"+proc )
        datacard.write("\t"+cat )
        datacard.write("\t"+opts.output )
        datacard.write("\tw:"+name )
        datacard.write("\n")

    wInput=fBkgIn.Get("w")
    wInput.Print()
    hist_data = wInput.data("datahist"+'_'+cat)
    hist_data.SetName("observed_data")
    getattr(w,'import')(hist_data)
    datacard.write("shapes")
    datacard.write("\tdata_obs" )
    datacard.write("\t"+cat )
    datacard.write("\t"+opts.output )
    datacard.write("\tw:"+"observed_data")
    datacard.write("\n")

#### OBSERVATION
datacard.write("-------------------------------------\n")
datacard.write("bin")
for cat in category:
    datacard.write("\t"+cat)
datacard.write("\n")
datacard.write("observation")
for cat in category:
    datacard.write("\t-1")
datacard.write("\n")

#### RATE
datacard.write("-------------------------------------\n")
datacard.write("bin\t")
#for cat in range(0,opts.ncat):
for cat in category:    
    for proc in mcAll:
        datacard.write("\t"+cat)
datacard.write("\n")
datacard.write("process\t")
for cat in category:    
    for proc in mcAll:
        datacard.write("\t"+proc)
datacard.write("\n")
datacard.write("process\t")
for cat in category:
    for idx,proc in enumerate(mcAll): ## TRICK, put the only signal first
        print(idx,proc,(-1)*idx)
        if (proc=='bkg'): datacard.write("\t1")
        else:
#            newIDX=(-1)*idx
            newIDX=ENUM[proc]
            datacard.write("\t%d"%newIDX)
datacard.write("\n")
datacard.write("rate\t")
for cat in category:
    for proc in mcAll:
        datacard.write("\t1")
#	datacard.write("\t1")
#        if savePdf:
#	    datacard.write("\t%.0f"%(opts.lumi) )
#        else:
#	    datacard.write("\t-1")
datacard.write("\n")
############ SYST########
datacard.write("-------------------------------------\n")

w.writeToFile(opts.output)
print("->Done")


