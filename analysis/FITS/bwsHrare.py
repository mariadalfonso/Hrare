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
parser.add_option("-m","--whichMeson",type='string',help="Which meson (_RhoCat or _PhiCat)", default="_RhoCat")
parser.add_option("-b","--whichBIN",type='string',help="Which bin of BDT (_bin1 or empty)", default="")
parser.add_option("-o","--output",type='string',help="Output ROOT file. [%default]", default="workspace_STAT_Rho_2018.root")
parser.add_option("-d","--datCardName",type='string',help="Output txt file. [%default]", default="datacard_STAT_Rho_2018.txt")

opts, args = parser.parse_args()

sys.argv=[]
import ROOT
ROOT.gROOT.SetBatch()

doSyst=True
MultiPdf=True

############ CONFIGURABLES ###########

if MultiPdf:
    BkgPdf={
        'Vcat': 'multipdf'+opts.whichBIN+opts.whichMeson,
        'Wcat': 'multipdf'+opts.whichBIN+opts.whichMeson,
        'Zcat': 'multipdf'+opts.whichBIN+opts.whichMeson,
        'VBFcat': 'multipdf'+opts.whichBIN+opts.whichMeson,
        'Zinvcat': 'multipdf'+opts.whichBIN+opts.whichMeson,
        'VBFcatlow': 'multipdf'+opts.whichBIN+opts.whichMeson,
        'GFcat': 'multipdf'+opts.whichBIN+opts.whichMeson,
    }
else:
    BkgPdf={
        'Vcat': 'exp1'+opts.whichBIN+opts.whichMeson,
        'Wcat': 'exp1'+opts.whichBIN+opts.whichMeson,
        'Zcat': 'exp1'+opts.whichBIN+opts.whichMeson,
        'VBFcat': 'bxg'+opts.whichBIN+opts.whichMeson,
        'Zinvcat': 'exp1'+opts.whichBIN+opts.whichMeson,
        'VBFcatlow': 'bxg'+opts.whichBIN+opts.whichMeson,
        'GFcat': 'bxg'+opts.whichBIN+opts.whichMeson,
    }

SigPdf={
    'Vcat': 'crystal_ball'+opts.whichBIN+opts.whichMeson,
    'Wcat': 'crystal_ball'+opts.whichBIN+opts.whichMeson,
    'Zcat': 'crystal_ball'+opts.whichBIN+opts.whichMeson,
    'VBFcat': 'crystal_ball'+opts.whichBIN+opts.whichMeson,
    'Zinvcat': 'crystal_ball'+opts.whichBIN+opts.whichMeson,
    'VBFcatlow': 'crystal_ball'+opts.whichBIN+opts.whichMeson,
    'GFcat': 'crystal_ball'+opts.whichBIN+opts.whichMeson,
}

ENUM={
    'ggH': 0,
    'VBFH': -1,
    'WH': -2,
    'ZH': -3,
    'ZinvH': -4,
    'WHl': -5,
    'ZHl': -6,
    'TTH': -7,
}

QCDscale={
    'ggH': '0.961/1.0039',
    'VBFH': '0.997/1.004',
    'WH': '0.993/1.006',
    'ZH': '0.995/1.005',
    'ZinvH': '0.995/1.005',
    'WHl': '0.993/1.006',
    'ZHl': '0.995/1.005',
    'TTH': '0.995/1.005', ## TMP
}

pdf_Higgs={
    'ggH': '0.968/1.032',
    'VBFH': '0.979/1.021',
    'WH': '0.98/1.020',
    'ZH': '0.981/1.019', #assume the majority is not ggZH
    'ZinvH': '0.981/1.019', #assume the majority is not ggZH
    'WHl': '0.98/1.020',
    'ZHl': '0.981/1.019', #assume the majority is not ggZH
    'TTH': '0.981/1.019', ## TMP
}

lumi={
    '_2016': '1.007',
    '_2017': '1.008',
    '_2018': '1.011',
}

def addSystematics(systname, systtype, value, whichProc, category, mcAll, datacard):

    datacard.write(systname+" \t"+systtype)
    for cat in category:
        print("cat",cat)
        for proc in mcAll:
            print("proc",proc)
            if (proc==whichProc): datacard.write("\t"+value)
            else:
                datacard.write("\t-")
    datacard.write("\n")

############ cat and processes ###########

print(opts.whichCat)

if opts.whichCat=='GFcat':
    sigAll = ['ggH','VBFH']
    mcAll = ['ggH','VBFH','bkgGF']
    category = ['GFcat']

if opts.whichCat=='Vcat':
    sigAll = ['WH','ZH','ZHl']
    mcAll = ['WH','ZH','ZHl','TTH','bkgV']
    category = ['Vcat']

if opts.whichCat=='Wcat':
    sigAll = ['WH','ZHl']
    mcAll = ['WH','ZHl','bkg']
    category = ['Wcat']

if opts.whichCat=='Zcat':
    sigAll = ['ZH']
    mcAll = ['ZH','bkg']
    category = ['Zcat']

if opts.whichCat=='VBFcat':
    sigAll = ['VBFH']
    mcAll = ['VBFH','bkgVBF']
    category = ['VBFcat']

if opts.whichCat=='VBFcatlow':
    sigAll = ['VBFH']
    mcAll = ['VBFH','bkgVBFlow']
    category = ['VBFcatlow']

if opts.whichCat=='Zinvcat':
    sigAll = ['WHl','ZinvH']
    mcAll = ['WHl','ZinvH','bkg']
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
if opts.whichCat=='Vcat':
    w.factory("mhVcat[100,150]"); # RooRealVar
    mh=w.var("mhVcat")
elif opts.whichCat=='GFcat':
    w.factory("mhGF[110,160]"); # RooRealVar
    mh=w.var("mhGF")
elif opts.whichCat=='VBFcat':
    w.factory("mhVBFcat[100,170]"); # RooRealVar
    mh=w.var("mhVBFcat")
elif opts.whichCat=='VBFcatlow':
    w.factory("mhVBFcatlow[100,170]"); # RooRealVar
    mh=w.var("mhVBFcatlow")
else:
    print('mg not specified')

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
        
        if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow':
            if opts.inputFileBKG != "":
                wInput=fBkgIn.Get("w")
                if opts.whichCat=='Vcat': name = BkgPdf[cat]+"_"+cat+"_"+'bkg'
                else: name = BkgPdf[cat]+"_"+cat+"_"+'bkg'
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
    hist_data = wInput.data("datahist"+opts.whichBIN+opts.whichMeson+'_'+cat)
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
        if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t1")
        else:
#            newIDX=(-1)*idx
            newIDX=ENUM[proc]
            datacard.write("\t%d"%newIDX)
datacard.write("\n")
datacard.write("rate\t")
for cat in category:
    for proc in mcAll:
#        datacard.write("\t1")
        if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t1")
        else:
            datacard.write("\t0.01")
#	datacard.write("\t1")
#        if savePdf:
#	    datacard.write("\t%.0f"%(opts.lumi) )
#        else:
#	    datacard.write("\t-1")
datacard.write("\n")

############ SYST ########
datacard.write("-------------------------------------\n")
if opts.whichCat=='GFcat' or opts.whichCat=='VBFcatlow': datacard.write("lumi_13TeV_2018 \tlnN ")
else: datacard.write("lumi_13TeV \tlnN ")
for cat in category:
    for proc in mcAll:
        if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t-")
        else:
            if opts.whichCat=='GFcat' or opts.whichCat=='VBFcatlow': datacard.write("\t%.3f"%(1.025) )
            else: datacard.write("\t%.3f"%(1.016) )
    datacard.write("\n")

if doSyst:

    if opts.whichCat=='GFcat' or opts.whichCat=='VBFcatlow' or opts.whichCat=='VBFcat':
        datacard.write("CMS_trig_PhotonLeg  \tlnN ")
        for cat in category:
            for proc in mcAll:
                if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t-")
                else:
                    if opts.whichMeson == '_RhoCat': datacard.write("\t%.3f"%(1.02) )
                    if opts.whichMeson == '_PhiCat': datacard.write("\t%.3f"%(1.02) )
                    if opts.whichMeson == '_K0StarCat': datacard.write("\t%.3f"%(1.02) )
            datacard.write("\n")

    if opts.whichCat=='GFcat' or opts.whichCat=='VBFcatlow':
        datacard.write("CMS_trig_TwoProngs  \tlnN ")
        for cat in category:
            for proc in mcAll:
                if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t-")
                else:
                    if opts.whichMeson == '_RhoCat': datacard.write("\t%.3f"%(1.05) )
                    if opts.whichMeson == '_PhiCat': datacard.write("\t%.3f"%(1.05) )
                    if opts.whichMeson == '_K0StarCat': datacard.write("\t%.3f"%(1.05) )
            datacard.write("\n")

    datacard.write("CMS_trackingEff  \tlnN ")
    for cat in category:
        for proc in mcAll:
            if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t-")
            else:
                if opts.whichMeson == '_RhoCat': datacard.write("\t%.3f"%(1.05) )
                if opts.whichMeson == '_PhiCat': datacard.write("\t%.3f"%(1.03) )
                if opts.whichMeson == '_K0StarCat': datacard.write("\t%.3f"%(1.05) )
        datacard.write("\n")

    datacard.write("CMS_IsoEff  \tlnN ")
    for cat in category:
        for proc in mcAll:
            if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t-")
            else:
                if opts.whichMeson == '_RhoCat': datacard.write("\t%.3f"%(1.01) )
                if opts.whichMeson == '_PhiCat': datacard.write("\t%.3f"%(1.01) )
                if opts.whichMeson == '_K0StarCat': datacard.write("\t%.3f"%(1.01) )
        datacard.write("\n")

    datacard.write("CMS_photonID \tlnN ")
    for cat in category:
        for proc in mcAll:
            if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t-")
            else:
                datacard.write("\t%.3f"%(1.01) )
        datacard.write("\n")

    if opts.whichCat=='Vcat' or opts.whichCat=='VBFcat':
        datacard.write("CMS_prefiring \tlnN ")
        for cat in category:
            for proc in mcAll:
                if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t-")
                else:
                    datacard.write("\t%.3f"%(1.005) )
            datacard.write("\n")

    datacard.write("CMS_pileup  \tlnN ")
    for cat in category:
        for proc in mcAll:
            if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t-")
            else:
                datacard.write("\t%.3f"%(1.01) )
        datacard.write("\n")

    for proc in mcAll:
        print("proc",proc)
        if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': continue
        else:
            # theo
            addSystematics("QCDscale_"+proc, 'lnN' ,QCDscale[proc], proc, category, mcAll, datacard)
            addSystematics("pdf_Higgs_"+proc, 'lnN' ,pdf_Higgs[proc], proc, category, mcAll, datacard)
#            if (proc=='ggH'):
#                addSystematics("pdf_Higgs_gg", 'lnN' ,pdf_Higgs[proc], proc, category, mcAll, datacard)
#            else: addSystematics("pdf_Higgs_qq", 'lnN' ,pdf_Higgs[proc], proc, category, mcAll, datacard)

if opts.whichCat=='GFcat' and doSyst:
    datacard.write("CMS_jes  \tlnN ")
    for cat in category:
        for proc in mcAll:
            if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t-")
            else:
                if(proc=='ggH'): datacard.write("\t0.99/1.01")
                if(proc=='VBFH'): datacard.write("\t0.97/1.03")
        datacard.write("\n")
    datacard.write("CMS_isr \tlnN ")
    for cat in category:
        for proc in mcAll:
            if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t-")
            else:
                if(proc=='ggH'): datacard.write("\t0.97/1.03")
                if(proc=='VBFH'): datacard.write("\t0.995/1.005")
        datacard.write("\n")
    datacard.write("CMS_fsr \tlnN ")
    for cat in category:
        for proc in mcAll:
            if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t-")
            else:
                if(proc=='ggH'): datacard.write("\t0.99/1.01")
                if(proc=='VBFH'): datacard.write("\t0.985/1.015")
        datacard.write("\n")

if (opts.whichCat=='VBFcat' or opts.whichCat=='VBFcatlow') and doSyst:
    datacard.write("CMS_jes  \tlnN ")
    for cat in category:
        for proc in mcAll:
            if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t-")
            else:
                if(proc=='VBFH'): datacard.write("\t0.97/1.03")
        datacard.write("\n")
    datacard.write("CMS_isr \tlnN ")
    for cat in category:
        for proc in mcAll:
            if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t-")
            else:
                if(proc=='VBFH'): datacard.write("\t0.995/1.005")
        datacard.write("\n")
    datacard.write("CMS_fsr \tlnN ")
    for cat in category:
        for proc in mcAll:
            if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t-")
            else:
                if(proc=='VBFH'): datacard.write("\t0.99/1.01")
        datacard.write("\n")

if opts.whichCat=='Vcat' and doSyst:
    datacard.write("CMS_eff_e  \tlnN ")
    for cat in category:
        for proc in mcAll:
            if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t-")
            else:
                datacard.write("\t%.3f"%(1.03) )
        datacard.write("\n")
    datacard.write("CMS_eff_m  \tlnN ")
    for cat in category:
        for proc in mcAll:
            if proc=='bkgGF' or proc=='bkgV' or proc=='bkgVBF' or proc=='bkgVBFlow': datacard.write("\t-")
            else:
                datacard.write("\t%.3f"%(1.01) )
        datacard.write("\n")

#    ### Add autoMCStats
#    datacard.write("\n")
#    datacard.write("* autoMCStats 0\n")

datacard.write("-------------------------------------\n")

if MultiPdf:
    pdfindexSTR='pdfindex'+opts.whichBIN+opts.whichMeson+"_"+opts.whichCat
    datacard.write("\n")
    datacard.write("%s discrete\n"%pdfindexSTR)

############ DONE ########
w.writeToFile(opts.output)
print("->Done")
