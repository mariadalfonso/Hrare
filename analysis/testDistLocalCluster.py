import ROOT
import json
#import sys
import os
from datetime import datetime

from utilsAna import pickTRG,findDIR,getMesonFromJson,getMVAFromJson
from utilsAna import loadUserCode, loadCorrectionSet
from utilsHrare import computeWeigths

doINTERACTIVE=False
useXROOTD=True
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

year = 2018
GOODphotons = "({0} or {1}) and Photon_pt>38 and Photon_electronVeto and Photon_mvaID_WP80 and abs(Photon_eta)<2.1".format(BARRELphotons,ENDCAPphotons) #90-80
GOODPHI = "{}".format(getMesonFromJson(mesons, "isZinv", "isPhiCat"))
GOODRHO = "{}".format(getMesonFromJson(mesons, "isZinv", "isRhoCat"))
GOODK0STAR = "{}".format(getMesonFromJson(mesons, "isZinv", "isK0StarCat"))

TRIGGER=pickTRG(TRIGGERS,year,"NULL",True,False,False,False)

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

def init():
    loadUserCode()
    loadCorrectionSet(2018)

if doINTERACTIVE:
    ROOT.ROOT.EnableImplicitMT()
    RDataFrame = ROOT.RDataFrame

else:
    RDataFrame = ROOT.RDF.Experimental.Distributed.Dask.RDataFrame
    import pkg_resources
    from utilsDask import create_local_connection
    from utilsDask import client_

def callVary(df):

    ROOT.gInterpreter.ProcessLine('corr_sf.eval_puSF(10,"nominal")')

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

def analysis():

    now = datetime.now()
    print('==> start: ',now)

    files = findDIR("/data/submit/cms/store/user/mariadlf/nano/D02/GluGlu_HToPhiGamma_M125_TuneCP5_PSWeights_13TeV_powheg_pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3+MINIAODSIM/",useXROOTD)

    if doINTERACTIVE:
        dfINI = RDataFrame("Events", files)
        rdf = RDataFrame("Runs", files)
        sumW = computeWeigths(dfINI, rdf, 1017, year, True, "False")
        init()
    else:
        if doLocal:
            connection = create_local_connection(10)
            dfINI = RDataFrame("Events", files, daskclient=connection)
            rdf = RDataFrame("Runs", files, daskclient=connection)
            sumW = computeWeigths(dfINI, rdf, 1017, year, True, "False")
        else:
            dfINI = RDataFrame("Events", files, daskclient=client_)
#            rdf = RDataFrame("Runs", files, daskclient=client_)
            sumW = 1. # fixme the computeWeigths( doens't work for remote slurm

        ROOT.RDF.Experimental.Distributed.initialize(init)

    print(GOODphotons)
    print(GOODPHI)
    print(TRIGGER)

    lumi = 1.
    weight = "{0}*genWeight*{1}".format(lumi,sumW)
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

    MassSum = df.Sum("goodPhotons_pt")
    res = MassSum.GetValue()
    print('sum of HCandMass  = ',res)

    df = callVary(df)
    df = callVariation(df) # only works for photon_pt (both interactive and distr)

    if True:
        print("writing plots")
        hists = {
#            "photon_pt":  {"name":"photon_pt","title":"Photon PT; pt_{#gamma} (GeV);N_{Events}","bin":200,"xmin":0,"xmax":200},
            "HCandMass":  {"name":"HCandMass","title":"H mass;m_{k^{+}k^{-}#gamma} (GeV);N_{Events}","bin":70,"xmin":100,"xmax":170},
        }

        histos = []
        for h in hists:

            # 1D is for nom only
            model1d = (hists[h]["name"], hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"])
            h1d = df.Histo1D(model1d, hists[h]["name"], "w")
            histos.append(h1d)
            print("h1d append")

            ## SYST with a weight
            model2d_pu = (hists[h]["name"]+":PU", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
            histos.append(df.Histo2D(model2d_pu, hists[h]["name"], "idx_nom_up_down", "pu_weights"))
            print("h2d append")

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

            outputFileHisto = "DASKlogs/histoOUTname.root"
            myfile = ROOT.TFile(outputFileHisto,"RECREATE")
            myfile.ls()

            print(histos)
            for h in histos:
                h.Write()
            myfile.Close()
            myfile.Write()

    if True:

        branchList = ROOT.vector('string')()
        for branchName in [
                "goodPhotons_pt",
                "HCandMass",
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

# In a Python script the Dask client needs to be initalized in a context
# Jupyter notebooks / Python session don't need this
if __name__ == "__main__":

    analysis()