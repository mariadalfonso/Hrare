import ROOT
import json
#import sys
import os
from datetime import datetime

from utilsAna import pickTRG,findDIR,getMesonFromJson,getMVAFromJson
from utilsAna import loadUserCode, loadCorrectionSet
from utilsHrare import computeWeigths
from utilsHrare import getMClist, getDATAlist, getSkims
from utilsAna import listDir

doINTERACTIVE=True
useXROOTD=False
doLocal=True

###-----------------------------------------
###-----------------------------------------

with open("./config/selection.json") as jsonFile:
    jsonObject = json.load(jsonFile)
    jsonFile.close()

with open("./config/trigger.json") as trgJsonFile:
    trgObject = json.load(trgJsonFile)
    trgJsonFile.close()

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

TRIGGER=pickTRG(TRIGGERS,year,"NULL",True,False,False,False)

MVA = jsonObject['MVAweights']
MVAweights = "{}".format(getMVAFromJson(MVA, "isGF" , "isPhiCat"))
#from urllib import request
#request.urlretrieve("root://submit30.mit.edu//mariadlf/utilFiles/weights_mva_oct/ggH_phi/TMVAClassification_wp80_6vars_mh110-160_BDTG.weights.xml")

print(GOODphotons)
print(GOODPHI)
print(TRIGGER)
print(MVAweights)

###-----------------------------------------
###-----------------------------------------

def init():
    loadUserCode()
    loadCorrectionSet(2018)
#    loadtmva_local()

def loadtmva_local():
    print('loadtmva_local() ')
    s ='''
    TMVA::Experimental::RReader model("weights_mva_oct/ggH_phi/TMVAClassification_wp80_6vars_mh110-160_BDTG.weights.xml");
    auto computeModel = TMVA::Experimental::Compute<6, float>(model);
    auto prediction = model.Compute({0.5, 1.0, -0.2, 1.5, 0.5, 1.0});
    '''
#    s ='''
#    TMVA::Experimental::RReader model("{0}");
#    auto computeModel = TMVA::Experimental::Compute<{1}, float>(model);
#    auto prediction = model.Compute({0.5, 1.0, -0.2, 1.5, 0.5, 1.0});
#    '''

    NVar=6
#    print(s.format(MVAweights,NVar))
    ROOT.gSystem.AccessPathName(MVAweights)
#    ROOT.gInterpreter.ProcessLine(s.format(MVAweights,NVar))
    ROOT.gInterpreter.ProcessLine(s)
    variables = ROOT.model.GetVariableNames()
    print('just inside the load',variables)

###-----------------------------------------
###-----------------------------------------

if doINTERACTIVE:
    ROOT.ROOT.EnableImplicitMT()
    RDataFrame = ROOT.RDataFrame
    RunGraphs = ROOT.RDF.RunGraphs
    init()
else:
    RDataFrame = ROOT.RDF.Experimental.Distributed.Dask.RDataFrame
    RunGraphs = ROOT.RDF.Experimental.Distributed.RunGraphs
    import pkg_resources
    from utilsDask import create_local_connection
#    from utilsDask import client
    ROOT.RDF.Experimental.Distributed.initialize(init)

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

def callMVA(df):

    if False:
        #need to be multicore
        tmva_helper_ = tmva_helper_xml.TMVAHelperXML(MVAweights)
        print(tmva_helper_.variables)

        print('printing tmva_helper_ inside callMVA: ',tmva_helper_)
    else:
        variables = ROOT.model.GetVariableNames()
        print('just inside call MVA: ',variables)
        print(list(variables))

    dfWithMVA = (df.Define("goodJets","{}".format(GOODJETS))
                 .Define("nGoodJets","Sum(goodJets)*1.0f")
#                 .Define("meson_DR", "goodMeson_DR[0]")
#                 .Define("meson_iso","goodMeson_iso[0]")
                 .Define("meson_DR", "(index_pair[0]!= -1) ? goodMeson_DR[index_pair[0]]: 0.f")
                 .Define("meson_iso","(index_pair[0]!= -1) ? goodMeson_iso[index_pair[0]]: 0.f")
                 .Define("HCandPT__div_HCandMass", "(HCandMass>0) ? HCandPT/HCandMass: 0.f")
#                .Define("photon_pt__daiv_HCandMass", "(HCandMass>0) ? photon_pt/HCandMass: 0.f")
#                 .Define("meson_pt__div_HCandMass", "(HCandMass>0) ? goodMeson_pt[index_pair[0]]/HCandMass: 0.f")
                 .Define("photon_pt__div_HCandMass", "(index_pair[1]!= -1 && HCandMass>0) ? photon_pt/HCandMass: 0.f")
                 .Define("meson_pt__div_HCandMass", "(index_pair[0]!= -1 && HCandMass>0) ? goodMeson_pt[index_pair[0]]/HCandMass: 0.f")
                 )

    if False:
        dfWithMVA = tmva_helper_.run_inference(dfWithMVA,"MVAdisc")
    else:
        dfWithMVA = dfWithMVA.Define("MVAdisc", ROOT.computeModel, list(ROOT.model.GetVariableNames()))

    return dfWithMVA

def callVary(df):

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

def callVariation(df):

    dfVaryPh=(df
              .Define("goodPhotons_dEsigmaUp", 'Photon_dEsigmaUp[goodPhotons]')
              .Define("goodPhotons_dEsigmaDown", 'Photon_dEsigmaDown[goodPhotons]')
              .Define("photon_dEsigmaUp",'(1.f+goodPhotons_dEsigmaUp[index_pair[1]])')
              .Define("photon_dEsigmaDown",'(1.f+goodPhotons_dEsigmaDown[index_pair[1]])')
              .Vary("photon_pt", "ROOT::RVecF{photon_pt*photon_dEsigmaDown,photon_pt*photon_dEsigmaUp}", variationTags=["dn","up"], variationName = "PhotonSYST")
#              .Define("goodPhoton_pt", 'goodPhotons_pt[0]')
#              .Vary("goodPhoton_pt", "ROOT::RVecF{goodPhoton_pt*0.9f,goodPhoton_pt*1.1f}", variationTags=["dn","up"], variationName = "PhotonSYST")
	      )
    return dfVaryPh

histos = []
evtcounts = []

def analysis(files,year,mc,sumW):

    now = datetime.now()
    print('==> start: ',now)

    if doINTERACTIVE:
        dfINI = RDataFrame("Events", files)
        #        loadtmva_helper()
    else:
        if doLocal:
            connection = create_local_connection(1)
            print(connection)
            # The `npartitions` optional argument tells the RDataFrame how many tasks are desired
            NPARTITIONS=300
            dfINI = RDataFrame("Events", files, daskclient=connection, npartitions=NPARTITIONS)
        else:
            dfINI = RDataFrame("Events", files, daskclient=client)

    lumi = 1.
    weight = "{0}".format(1.)
    if mc>=0: weight = "{0}*genWeight*{1}".format(lumi,sumW)
    if year==2018: lumiIntegrated = lumis['2018']

    df = (dfINI
          .Filter("nPhoton>0 and PV_npvsGood>0","photon from nano >0 and PV_npvsGood > 0")
          .Define("w","{}".format(weight))
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

#    df = callMVA(df)
#    df = callVary(df)
#    df = callVariation(df) # only works for photon_pt (both interactive and distr)

    evtcounts.append(df.Count())

    if True:
        print("writing plots")
        hists = {
#            "photon_pt":  {"name":"photon_pt","title":"Photon PT; pt_{#gamma} (GeV);N_{Events}","bin":200,"xmin":0,"xmax":200},
            "HCandMass":  {"name":"HCandMass","title":"H mass;m_{k^{+}k^{-}#gamma} (GeV);N_{Events}","bin":70,"xmin":100,"xmax":170},
        }

        for h in hists:

            # 1D is for nom only
            model1d = (hists[h]["name"]+"_"+str(year)+"_"+str(mc), hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"])
            h1d = df.Histo1D(model1d, hists[h]["name"], "w")
            histos.append(h1d)
            print("h1d append")

            ## SYST with a weight
#            model2d_pu = (hists[h]["name"]+":PU", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
#            histos.append(df.Histo2D(model2d_pu, hists[h]["name"], "idx_nom_up_down", "pu_weights"))
#            print("h2d append")

#            all_hs = ROOT.RDF.Experimental.VariationsFor(h1d);
#            all_hs.GetKeys();
#            print(all_hs.GetKeys())

            ## to use the SYST that change the variable  (for some reason this is associated only to the photon+pt and not the HCandMass)
#            if doINTERACTIVE: hx = ROOT.RDF.Experimental.VariationsFor(h1d)
#            else: hx = ROOT.RDF.Experimental.Distributed.VariationsFor(h1d)
#            hx["PhotonSYST:dn"].SetName(hists[h]["name"]+":PhotonSYST:dn")
#            histos.append(hx["PhotonSYST:dn"])
#            hx["PhotonSYST:up"].SetName(hists[h]["name"]+":PhotonSYST:up")
#            histos.append(hx["PhotonSYST:up"])


#    if not doINTERACTIVE:
##    if False:
#        outputFileHisto = "DASKlogs/histoOUTname_mc{0}_{1}.root".format(mc,year)
#        myfile = ROOT.TFile(outputFileHisto,"RECREATE")
#        myfile.ls()

#        print(histos)
#        for h in histos:
#            h.Write()
#        myfile.Close()
#        myfile.Write()

    if False:

        branchList = ROOT.vector('string')()
        for branchName in [
                "goodPhotons_pt",
                "HCandMass",
#                "MVAdisc",
        ]:
            branchList.push_back(branchName)

        outputFile = "DASKlogs/snapshotOUTname.root"
        print(outputFile)
        snapshotOptions = ROOT.RDF.RSnapshotOptions()
        snapshotOptions.fCompressionAlgorithm = ROOT.kLZ4
        snapshot_tdf = df.Snapshot("events", outputFile, branchList, snapshotOptions)
        print("snapshot_tdf DONE")
        print(outputFile)

#    return histos
    now = datetime.now()
    print('==> ends: ',now)

def BuildDict():

    isT2=False
    if isT2:
        dirName="/store/user/paus/nanohr/D02/"
    else:
        dirName="/store/user/mariadlf/nano/D02/"

    campaignv3 = "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3+MINIAODSIM/"
    campaignv1 = "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM/"
    campaignv2 = "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM/"
    campaignFIX = "RunIISummer20UL18MiniAODv2-4cores5k_106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM/"

    thisdict = {
        1010: (listDir(isT2,dirName+"VBF_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaignv1),3781.7*0.49),
        1017: (listDir(isT2,dirName+"GluGlu_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+"+campaignv3),48580*0.49),
        10: (listDir(isT2,dirName+"GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaignv2),18540.0*1000*1.26), #LO *1.26
        11: (listDir(isT2,dirName+"GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaignFIX),8644.0*1000*1.26), #LO *1.26
        12: (listDir(isT2,dirName+"GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaignv2),2183.0*1000*1.26), #LO *1.26
        13: (listDir(isT2,dirName+"GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaignv2),260.2*1000*1.26), #LO *1.26
        14: (listDir(isT2,dirName+"GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8+"+campaignv2),86.58*1000*1.26) #LO *1.26
    }

    return thisdict

def SwitchSample(thisdict,argument):

    return thisdict.get(argument, "BKGdefault, xsecDefault")

def loopOnDataset():

    thisdict = BuildDict()

    year=2018
    mc = [1017,1010,10,11,12,13,14]

    for sampleNOW in mc:
        files = SwitchSample(thisdict,sampleNOW)[0]
        print('outside the function: ',len(files))
        rdf = ROOT.RDataFrame("Runs", files) # make sure this is not the distributed
        sumW = computeWeigths(rdf,mc, year, True, "False")
        analysis(files,year,sampleNOW,sumW)

#    print('DONE with MC going for Data')
#    data = [-62,-63,-64]

#    for datasetNumber in data:
#         pair = getSkims(datasetNumber,year,"Zinv",useD03)
#         files = pair[0]
#         PDType = pair[1]
#         print(len(files))
#         print(PDType)
#         analysis(files,year,datasetNumber,1.)

def loopOnDatasetLocal():

    useD03=False
    year=2018

    mc = [1017,1010,10,11,12,13,14]

    for sampleNOW in mc:
        files = getMClist(year,sampleNOW,useD03)
        print(len(files))
        rdf = ROOT.RDataFrame("Runs", files) # make sure this is not the distributed
        sumW = computeWeigths(rdf,sampleNOW, year, True, "False")
        analysis(files,year,sampleNOW,sumW)

    data = [-62,-63,-64]
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

    loopOnDatasetLocal()

    if True:

        now = datetime.now()
        print('==> start runGraph: ',now)
        RunGraphs(evtcounts);

        now = datetime.now()
        print('==> done runGraph: ',now)

        if doINTERACTIVE: outputFileHisto = "DASKlogs/histoOUTname_test_interactive.root"
        else: outputFileHisto = "DASKlogs/histoOUTname_test_LocalCluster.root"
        myfile = ROOT.TFile(outputFileHisto,"RECREATE")
        myfile.ls()

        print(histos)
        for h in histos:
            h.Write()
        myfile.Close()
        myfile.Write()

        now = datetime.now()
        print('==> histo done: ',now)
