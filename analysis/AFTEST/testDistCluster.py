import ROOT
import json
import sys
import os
from datetime import datetime


#AF="MIT"
#AF="MIT-dask"
AF="MIT-daskgateway"
#AF="EAF-daskgateway"
#AF="Purdue"
#AF="cern-dask"
#AF="cern-spark"

#MIT
if AF=="MIT-dask" or AF=="MIT-daskgateway" or AF=="MIT": myDir = '/home/submit/mariadlf/Hrare/CMSSW_10_6_27_new/src/Hrare/analysis/AFTEST/'
#FNAL-EAF
if AF=="EAF-daskgateway": myDir = '/home/dalfonso/Hrare/analysis/AFTEST/'
#Purdue
if AF=="Purdue": myDir = '/home/dalfonso-cern/Hrare/analysis/AFTEST/'
if AF=="cern-dask" or AF=="cern-spark": myDir = '/eos/user/d/dalfonso/SWAN_projects/Hrare/JULY_exp/'

sys.path.insert(0, myDir)

from utilsAna import pickTRG,findDIR,getMesonFromJson,getMVAFromJson
from utilsAna import loadUserCode, loadCorrectionSet, loadUserCodeGateway
from utilsAna import BuildDictMgamma, SwitchSample
#from utilsAna import readDataQuality

doINTERACTIVE=False
doLocalCluster=False

doLocalDisk=True
doRemoteAccess=False

###-----------------------------------------
###-----------------------------------------

with open("../config/selection.json") as jsonFile:
    jsonObject = json.load(jsonFile)
    jsonFile.close()

with open("../config/trigger.json") as trgJsonFile:
    trgObject = json.load(trgJsonFile)
    trgJsonFile.close()

with open("../config/qualityData.json") as qualityJsonFile:
    qualObject = json.load(qualityJsonFile)
    qualityJsonFile.close()

JSON = qualObject['JSON']

BARRELphotons = jsonObject['BARRELphotons']
ENDCAPphotons = jsonObject['ENDCAPphotons']
TRIGGERS = trgObject['triggers']
mesons = jsonObject['mesons']
GOODJETS = jsonObject['GOODJETS']

year = 2018

GOODphotons = "({0} or {1}) and Photon_pt>38 and Photon_electronVeto and Photon_mvaID_WP80 and abs(Photon_eta)<2.1".format(BARRELphotons,ENDCAPphotons) #90-80
GOODPHI = "{}".format(getMesonFromJson(mesons, "isZinv", "isPhiCat"))
GOODRHO = "{}".format(getMesonFromJson(mesons, "isZinv", "isRhoCat"))
GOODK0STAR = "{}".format(getMesonFromJson(mesons, "isZinv", "isK0StarCat"))

TRIGGER=pickTRG(TRIGGERS,year,"NULL",True,False,False,False,False)

#print(GOODphotons)
#print(GOODPHI)
#print(TRIGGER)

lumis={
    '12016': 19.52, #APV #(B-F for 2016 pre)
    '22016': 16.80, #postVFP
    '2016': 35.9,
    '2017': 41.5,
    '12017': 7.7, #(F for 2017)
    '2018': 59.70,
    '12018': 39.54,
    'all': 86.92,      #19.52 + 7.7 + 59.70
}

###-----------------------------------------
###-----------------------------------------

def initSpark():
    from pathlib import Path
    from pyspark import SparkFiles
    print('loadUserCode.h')
    localdir = SparkFiles.getRootDirectory()
    lib_path = Path(localdir) / "functions.h"
#    lib_path = Path(localdir) / "myLibrary.h"
    ROOT.gInterpreter.ProcessLine(f'#include "{lib_path}"')
#    ROOT.gInterpreter.Declare(f'#include "{lib_path}"')
    loadCorrectionSet(2018)

def init():

    if doINTERACTIVE or doLocalCluster:
        loadUserCode()
        loadCorrectionSet(2018)
    else:
        if AF=="MIT-daskgateway"  or AF=="Purdue" or AF=="EAF-daskgateway":
            loadUserCodeGateway("functions.h")
#            loadUserCodeGateway("sfCorrLib.h") # comment becuase of #include "correction.h"
        else:
            loadUserCode()
            loadCorrectionSet(2018)
            #    readDataQuality(2018)

###-----------------------------------------
###-----------------------------------------

if doINTERACTIVE:
    ROOT.ROOT.EnableImplicitMT()
    RDataFrame = ROOT.RDataFrame
    RunGraphs = ROOT.RDF.RunGraphs
else:
#    import pkg_resources

    if AF=="cern-spark":
        import pyspark
        RDataFrame = ROOT.RDF.Experimental.Distributed.Spark.RDataFrame # for CERN
    else:
        RDataFrame = ROOT.RDF.Experimental.Distributed.Dask.RDataFrame
    RunGraphs = ROOT.RDF.Experimental.Distributed.RunGraphs
    if doLocalCluster:
        print('Hello doing Local Cluster')
        from utilsDask import create_local_connection
    else:
#
        if AF=="MIT-dask": from utilsDask import create_remote_Dask
        elif AF=="MIT-daskgateway": from utilsDask import create_DaskGateway_MIT
        elif AF=="EAF-daskgateway": from utilsDask import create_EAF_connection
        elif AF=="Purdue": from utilsDask import create_Purdue_connection

def makeRDF(files):

    if doINTERACTIVE:
        df = ROOT.RDataFrame("Events", files)
        init()
    else:
        #files_per_partition = nfiles / npartitions
#        NCLUSTERS=100 # this is the numbers of cluster to open
        if doLocalCluster:
            NPARTITIONS=5 # assuming one launch 7 Graphs (7*15 = 105 cores)
            nCLUSTER=10 #(like scale for remote)
#            connection = create_local_connection(NPARTITIONS)
            connection = create_local_connection(nCLUSTER)
            print(connection)
            # The `npartitions` optional argument tells the RDataFrame how many tasks are desired
            df = RDataFrame("Events", files, daskclient=connection, npartitions=NPARTITIONS)
            ROOT.RDF.Experimental.Distributed.initialize(init)
        else:

            #Create the client
            if AF=="cern-spark":
                sc.addPyFile(myDir+"utilsAna.py")
            else:
                if AF=="MIT-dask": client = create_remote_Dask()
                elif AF=="MIT-daskgateway": client = create_DaskGateway_MIT()
                elif AF=="EAF-daskgateway": client = create_EAF_connection()
                elif AF=="cern-dask":
                    from dask.distributed import Client
                    client = Client("tls://10.100.18.182:31037")
                client.upload_file(myDir+"utilsAna.py")

            print(client)
            NPARTITIONS=5 # assuming one launch 10 Graphs (10*30 = 200 cores)

            #Create the RDF and initialize user lib
            if AF=="cern-spark":
                df = RDataFrame("Events", files, sparkcontext=sc, npartitions=NPARTITIONS) # at CERN
            else:
                df = RDataFrame("Events", files, daskclient=client, npartitions=NPARTITIONS)

            if AF=="cern-spark" or AF=="cern-dask":
                df._headnode.backend.distribute_unique_paths(
                    [
                        myDir+"config/functions.h",
                    ]
                )
            else:
                df._headnode.backend.distribute_unique_paths(
                    [
                        myDir+"../config/functions.h",
#                        myDir+"../config/sfCorrLib.h",
                    ]
                )

            if AF=="cern-spark":
                ROOT.RDF.Experimental.Distributed.initialize(initSpark)
            else:
                ROOT.RDF.Experimental.Distributed.initialize(init)

    return df
    
###-----------------------------------------
###-----------------------------------------

def callSFVariation(df):

    dfVaryPh=(df
              .Define("idx_nom_up_down", "indices(3)")
              .Define("SFpu_A",'corr_sf.eval_puSF(10,"nominal")')
              .Define("SFpu_Nom",'corr_sf.eval_puSF(Pileup_nTrueInt,"nominal")')
              .Define("SFpu_Up",'corr_sf.eval_puSF(Pileup_nTrueInt,"up")')
              .Define("SFpu_Dn",'corr_sf.eval_puSF(Pileup_nTrueInt,"down")')
              #
              .Define("w_allSF", "w*SFpu_Nom")
              .Define("pu_weights", "NomUpDownVar(SFpu_Nom, SFpu_Up, SFpu_Dn, w_allSF)")
              #
              )
    return dfVaryPh

def callVary(df):


#   if df.HasColumn("Photon_dEsigmaUp") and df.HasColumn("Photon_dEsigmaDown"):
    
    dfVaryPh=(df
              .Define("goodPhotons_dEsigmaUp", 'Photon_dEsigmaUp[goodPhotons]')
              .Define("goodPhotons_dEsigmaDown", 'Photon_dEsigmaDown[goodPhotons]')
              .Define("photon_dEsigmaUp",'(1.f+goodPhotons_dEsigmaUp[index_pair[1]])')
              .Define("photon_dEsigmaDown",'(1.f+goodPhotons_dEsigmaDown[index_pair[1]])')
              .Vary("photon_pt", "ROOT::RVecF{photon_pt*photon_dEsigmaDown,photon_pt*photon_dEsigmaUp}", variationTags=["dn","up"], variationName = "PhotonSYST")
#              .Define("goodPhoton_pt", 'goodPhotons_pt[0]')
#              .Vary("photon_pt", "ROOT::RVecF{photon_pt*0.9f,photon_pt*1.1f}", variationTags=["dn","up"], variationName = "PhotonSYST")
	      )
    return dfVaryPh

histos = []
evtcounts = []

def analysis(files,year,mc,sumW):

    now = datetime.now()
    print('==> start: ',now)

    dfINI=makeRDF(files)

    lumi = 1.
    weight = "{0}".format(1.)
    if mc>=0: weight = "{0}*genWeight*{1}".format(lumi,sumW)
    if year==2018: lumiIntegrated = lumis['2018']

    isData = "false"
    if mc<0: isData = "true"

    df = (dfINI
          .Filter("nPhoton>0 and PV_npvsGood>0","photon from nano >0 and PV_npvsGood > 0")
          .Define("w","{}".format(weight))
          .Define("isData","{}".format(isData))
#          .Define("applyJson","{}".format(JSON)).Filter("applyJson","pass JSON") # this doens't wotk to check why
          .Define("lumiIntegrated","{}".format(lumiIntegrated))
          .Define("triggerAna","{}".format(TRIGGER))
          .Filter("triggerAna>0", "pass triggers")
          .Define("goodPhotons", "{}".format(GOODphotons))
          .Define("nGoodPhotons","Sum(goodPhotons)*1.0f")
          .Filter("Sum(goodPhotons)>0", "At least one good Photon")
          .Define("goodPhotons_pt", "Photon_pt[goodPhotons]")
          .Define("goodPhotons_eta", "Photon_eta[goodPhotons]")
          .Define("goodPhotons_phi", "Photon_phi[goodPhotons]")
          .Filter("nphi>0").Define("goodMeson","({}".format(GOODPHI)+")")
          .Filter("Sum(goodMeson)>0", "one good Phi (ptPhi, validfit, ptTracks)")
          .Define("goodMeson_pt", "phi_kin_pt[goodMeson]")
          .Define("goodMeson_eta", "phi_kin_eta[goodMeson]")
          .Define("goodMeson_phi", "phi_kin_phi[goodMeson]")
          .Define("goodMeson_mass", "phi_kin_mass[goodMeson]")
          #
          .Define("goodMeson_iso", "phi_iso[goodMeson]")
          .Define("goodMeson_DR","DeltaR(phi_trk1_eta[goodMeson],phi_trk2_eta[goodMeson],phi_trk1_phi[goodMeson],phi_trk2_phi[goodMeson])")
          .Define("goodMeson_trk1_pt", "phi_trk1_pt[goodMeson]")
          .Define("goodMeson_trk2_pt", "phi_trk2_pt[goodMeson]")
          .Define("wrongMeson","({}".format(GOODRHO)+")")
          .Define("wrongMeson_pt","Sum(wrongMeson) > 0 ? rho_kin_pt[wrongMeson]: ROOT::VecOps::RVec<float>(0.f)")
          .Define("wrongMeson2","({}".format(GOODK0STAR)+")")
          .Define("wrongMeson2_pt","Sum(wrongMeson2) > 0 ? K0Star_kin_pt[wrongMeson2]: ROOT::VecOps::RVec<float>(0.f)")
          #
          .Define("index_pair","HiggsCandFromRECO(goodMeson_pt, goodMeson_eta, goodMeson_phi, goodMeson_mass, goodMeson_trk1_pt, goodMeson_trk2_pt, wrongMeson_pt, wrongMeson2_pt, goodPhotons_pt, goodPhotons_eta, goodPhotons_phi)").Filter("index_pair[0]!= -1", "at least a good meson candidate")
          .Define("meson_pt", "(index_pair[0]!= -1) ? goodMeson_pt[index_pair[0]]: 0.f")
          .Define("photon_pt", "(index_pair[1]!= -1) ? goodPhotons_pt[index_pair[1]]: 0.f")
          .Define("HCandMass", "compute_HiggsVars_var(meson_pt,goodMeson_eta[index_pair[0]],goodMeson_phi[index_pair[0]],goodMeson_mass[index_pair[0]],photon_pt,goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]],0)")
          .Define("HCandPT",   "compute_HiggsVars_var(meson_pt,goodMeson_eta[index_pair[0]],goodMeson_phi[index_pair[0]],goodMeson_mass[index_pair[0]],photon_pt,goodPhotons_eta[index_pair[1]],goodPhotons_phi[index_pair[1]],1)")
          #
          )

#    if mc>0:
#        df = callSFVariation(df)
#        df = callVary(df)
        
    evtcounts.append(df.Count())

    if True:
        print("writing plots")
        hists = {
#            "photon_pt":  {"name":"photon_pt","title":"Photon PT; pt_{#gamma} (GeV);N_{Events}","bin":200,"xmin":0,"xmax":200},
            "HCandMass":  {"name":"HCandMass","title":"H mass;m_{k^{+}k^{-}#gamma} (GeV);N_{Events}","bin":70,"xmin":100,"xmax":170},
#            "MVAdisc":  {"name":"MVAdisc","title":"MVA; MVA discr;N_{Events}","bin":200,"xmin":-1.,"xmax":1.},                        
        }

        for h in hists:

            # 1D is for nom only
            model1d = (hists[h]["name"]+"_"+str(year)+"_"+str(mc), hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"])
            h1d = df.Histo1D(model1d, hists[h]["name"], "w")
            histos.append(h1d)
            print("h1d.GetName()",h1d.GetName())
            print("h1d append")

            if False and mc>1000:
                model2d_pu = (hists[h]["name"]+"_"+str(year)+"_"+str(mc)+":PU", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
                histos.append(df.Histo2D(model2d_pu, hists[h]["name"], "idx_nom_up_down", "pu_weights"))

                ## to use the SYST that change the variable
#                if doINTERACTIVE: hx = ROOT.RDF.Experimental.VariationsFor(h1d)
#                else: hx = ROOT.RDF.Experimental.Distributed.VariationsFor(h1d)
#                for hist_handle in hx:
#                    print(hist_handle[1].GetName())
#                hx["PhotonSYST:dn"].SetName(hists[h]["name"]+":PhotonSYST:dn")
#                histos.append(hx["PhotonSYST:dn"])
#                hx["PhotonSYST:up"].SetName(hists[h]["name"]+":PhotonSYST:up")
#                histos.append(hx["PhotonSYST:up"])
#                histos.append(hx["PhotonSYST:up"])

                ## to use the SYST that change the variable
#                hx = ROOT.RDF.Experimental.VariationsFor(h1d.GetValue())
                #                hx = ROOT.RDF.Experimental.VariationsFor(h1d.GetValue()) 
                # pair<const string,shared_ptr<TH1D> >
#                for hist_handle in hx:
#                    print(hist_handle[1].GetName())                                
#                hx["PhotonSYST:dn"].SetName(hists[h]["name"]+"_"+str(year)+"_"+str(mc)+":PhotonSYST:dn")
#                histos.append(hx["PhotonSYST:dn"])
#                hx["PhotonSYST:up"].SetName(hists[h]["name"]+"_"+str(year)+"_"+str(mc)+":PhotonSYST:up")
#                histos.append(hx["PhotonSYST:up"])


    if False:

        branchList = ROOT.vector('string')()
        for branchName in [
                "goodPhotons_pt",
                "HCandMass",
                "MVAdisc"
        ]:
            branchList.push_back(branchName)

        outputFile = "DASKlogs/snapshotOUTname.root"
        print(outputFile)
        snapshotOptions = ROOT.RDF.RSnapshotOptions()
        snapshotOptions.fCompressionAlgorithm = ROOT.kLZ4
        snapshot_tdf = df.Snapshot("events", outputFile, branchList, snapshotOptions)
        print("snapshot_tdf DONE")
        print(outputFile)

    now = datetime.now()
    print('==> ends: ',now)

def loopRemoteDataset():

    # used also in Purdue
    from utilsAna import computeWeigths
    
    thisdict = BuildDictMgamma()

    year=2018
    mc = [1017,1010,10,11,12,13,14]
#    mc = []

    for sampleNOW in mc:
        files = SwitchSample(thisdict,sampleNOW)[0]
        xsec = SwitchSample(thisdict,sampleNOW)[1]
        print('outside the function: ',len(files))
        rdf = ROOT.RDataFrame("Runs", files) # make sure this is not the distributed
        sumW = computeWeigths(rdf,xsec)
        analysis(files,year,sampleNOW,sumW)

    data = [-62,-63,-64]
#    data = []

    for datasetNumber in data:
        files = SwitchSample(thisdict,datasetNumber)[0]
        print(len(files))
        analysis(files,year,datasetNumber,1.)

def loopOnDatasetLocalMIT():

    print('inside loopOnDatasetLocalMIT')
    # localDataset -- as used in the VGammaAnalysis
    from utilsHrare import getMClist, getSkims
    from utilsHrare import computeWeigths

    useD03=False
    year=2018

    mc = [1017,1010,10,11,12,13,14]
#    mc = []

    for sampleNOW in mc:
        files = getMClist(year,sampleNOW,useD03)
        print(len(files))
        rdf = ROOT.RDataFrame("Runs", files) # make sure this is not the distributed
        sumW = computeWeigths(rdf,sampleNOW, year, True, "False")
        analysis(files,year,sampleNOW,sumW)

    data = [-62,-63,-64]
#    data = []
#    readDataQuality(year)

    for datasetNumber in data:
         pair = getSkims(datasetNumber,year,"Zinv",useD03)
         files = pair[0]
         PDType = pair[1]
         print(len(files))
         print(PDType)
         analysis(files,year,datasetNumber,1.)


# In a Python script the Dask client needs to be initalized in a context
# Jupyter notebooks / Python session don't need this
if __name__ == "__main__":

    now = datetime.now()
    print('==> very beginning: ',now)

    try:
        if doLocalDisk: loopOnDatasetLocalMIT()
        elif doRemoteAccess: loopRemoteDataset()
    except:
        print("exception caught, sleeping for 10 minutes")
        import time
        time.sleep(600)


    if True:

        now = datetime.now()
        print('==> start runGraph: ',now)
        RunGraphs(evtcounts);

        now = datetime.now()
        print('==> done runGraph: ',now)

        if doINTERACTIVE: outputFileHisto = "DASKlogs/histoOUTname_test_interactive.root"
        else:
            if doLocalCluster: outputFileHisto = "DASKlogs/histoOUTname_test_LocalCluster.root"
            else:
                if AF=="MIT-dask": outputFileHisto = "DASKlogs/histoOUTname_test_RemoteCluster_Dask.root"
                elif AF=="MIT-daskgateway": outputFileHisto = "DASKlogs/histoOUTname_test_RemoteCluster_Gateway.root"
                else: outputFileHisto = "DASKlogs/histoOUTname_test_RemoteCluster.root"

        myfile = ROOT.TFile(outputFileHisto,"RECREATE")
        myfile.ls()

        print(histos)
        for h in histos:
            h.Write()
        myfile.Close()
        myfile.Write()

        now = datetime.now()
        print('==> histo done: ',now)
