import ROOT
import json
import os
from datetime import datetime

from utilsHrare import pickTRG,findDIR

doINTERACTIVE=False
useXROOTD=False

if doINTERACTIVE:
    ROOT.ROOT.EnableImplicitMT()
    from utilsHrare import loadCorrectionSet
    RDataFrame = ROOT.RDataFrame
    loadCorrectionSet(2018)
    
else:
    RDataFrame = ROOT.RDF.Experimental.Distributed.Dask.RDataFrame
    import pkg_resources
    from utilsDask import client_

    import correctionlib
    # do I really need this pyroot_binding since I'm importing the libraries with pkg_resources
    correctionlib.register_pyroot_binding()
    def sf():
        ROOT.gInterpreter.Declare('''
        #ifndef MYFUN
        #define MYFUN
        auto mySFfun(){
        auto corr_sf = MyCorrections(2018);
        return corr_sf;
        }
        #endif
        ''')

with open("config/selection.json") as jsonFile:
    jsonObject = json.load(jsonFile)
    jsonFile.close()

with open("config/trigger.json") as trgJsonFile:
    trgObject = json.load(trgJsonFile)
    trgJsonFile.close()

BARRELphotons = jsonObject['BARRELphotons']
ENDCAPphotons = jsonObject['ENDCAPphotons']
TRIGGERS = trgObject['triggers']

year = 2018
GOODphotons = "({0} or {1}) and Photon_pt>38 and Photon_electronVeto and abs(Photon_eta)<2.1".format(BARRELphotons,ENDCAPphotons) #90-80

TRIGGER=pickTRG(TRIGGERS,year,"NULL",True,False,False,False)

def analysis():

    now = datetime.now()
    print('==> start: ',now)
    
    files = findDIR("/data/submit/cms/store/user/mariadlf/nano/D01/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8+RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1+MINIAODSIM/",useXROOTD)

    if doINTERACTIVE:
        dfINI = RDataFrame("Events", files)
    else:
        dfINI = RDataFrame("Events", files, daskclient=client_)
        #user code (not sure I need the .json files since those strings are already looked up)
        dfINI._headnode.backend.distribute_shared_libraries("functions_cc.so")
        dfINI._headnode.backend.distribute_shared_libraries("functions_cc_ACLiC_dict_rdict.pcm")
        dfINI._headnode.backend.distribute_headers("config/selection.json")
        dfINI._headnode.backend.distribute_headers("config/trigger.json")

        #correctionlib for SF
        dfINI._headnode.backend.distribute_headers(pkg_resources.resource_filename("correctionlib", "include/correction.h"))
        dfINI._headnode.backend.distribute_headers(pkg_resources.resource_filename("correctionlib", "include/correctionlib_version.h"))
        dfINI._headnode.backend.distribute_headers("config/mysf.h")
        dfINI._headnode.backend.distribute_shared_libraries("mysf.so")

        ROOT.RDF.Experimental.Distributed.initialize(sf)

    print(GOODphotons)

    df = (dfINI
          .Filter("nPhoton>0 and PV_npvsGood>0","photon from nano >0 and PV_npvsGood > 0")
          .Define("w","{}".format(1))
#          .Define("triggerAna","{}".format(TRIGGER))
#          .Filter("triggerAna>0", "pass triggers")  ## comment while doing trigger studies
          .Define("goodPhotons", "{}".format(GOODphotons))
          .Define("nGoodPhotons","Sum(goodPhotons)*1.0f")
          .Filter("Sum(goodPhotons)>0", "At least one good Photon")
          .Define("goodPhotons_pt", "Photon_pt[goodPhotons][0]")
          .Define("goodPhotons_eta", "Photon_eta[goodPhotons][0]")
#          ##
          .Define("SFphoton_ID_Nom",'corr_sf.eval_photonSF("{0}", "sf", "{1}", goodPhotons_eta, goodPhotons_pt)'.format(year,"wp90"))
          .Define("SFphoton_ID_Up",'corr_sf.eval_photonSF("{0}", "sfup", "{1}" , goodPhotons_eta, goodPhotons_pt)'.format(year,"wp90"))
          .Define("SFphoton_ID_Dn",'corr_sf.eval_photonSF("{0}", "sfdown", "{1}" , goodPhotons_eta, goodPhotons_pt)'.format(year,"wp90"))
#          .Define("SFphoton_ID_Nom",'mySFfun().eval_photonSF("{0}", "sf", "{1}", goodPhotons_eta[0], goodPhotons_pt[0])'.format(year,"wp90"))
#          .Define("SFphoton_ID_Up",'mySFfun().eval_photonSF("{0}", "sfup", "{1}" , goodPhotons_eta[0], goodPhotons_pt[0])'.format(year,"wp90"))
#          .Define("SFphoton_ID_Dn",'mySFfun().eval_photonSF("{0}", "sfdown", "{1}" , goodPhotons_eta[0], goodPhotons_pt[0])'.format(year,"wp90"))
          .Define("idx_nom_up_down", "indices(3)")
          .Define("phoID_weights", "NomUpDownVar(SFphoton_ID_Nom, SFphoton_ID_Up, SFphoton_ID_Dn, w)")
          )

    ptsum = df.Sum("Muon_pt")
    res = ptsum.GetValue()
    print('sum of muons = ',res)

    if False:
        print("writing plots")
        hists = {
            "goodPhotons_pt":  {"name":"goodPhotons_pt","title":"Photon PT; pt_{#gamma} (GeV);N_{Events}","bin":200,"xmin":0,"xmax":200}
        }

        histos = []
        for h in hists:

            # 1D is for nom only
            model = (hists[h]["name"], hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"])
            h1d = df.Histo1D(model, hists[h]["name"], "w")
            histos.append(h1d)

            ## those that change the weights only
            # 2D is for nom, up, down
            model2d_phoID = (hists[h]["name"]+":phoID", hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"], 3, 0, 3)
            histos.append(df.Histo2D(model2d_phoID, hists[h]["name"], "idx_nom_up_down", "phoID_weights"))

            outputFileHisto = "DASKlogs/histoOUTname.root".format(year)
            myfile = ROOT.TFile(outputFileHisto,"RECREATE")

            for h in histos:
                h.Write()
            myfile.Close()
            myfile.Write()

    if False:
        branchList = ROOT.vector('string')()

        for branchName in [
                "goodPhotons_pt",
                "Muon_pt",
                "PV_npvsGood",
                "run",
                "luminosityBlock",
                "event",
                "SFphoton_ID_Nom",
        ]:
            branchList.push_back(branchName)

            snapshotOptions = ROOT.RDF.RSnapshotOptions()
            snapshotOptions.fCompressionAlgorithm = ROOT.kLZ4
            snapshot_tdf = df.Snapshot("events", "DASKlogs/output.root", branchList, snapshotOptions)

    now = datetime.now()
    print('==> ends: ',now)

# In a Python script the Dask client needs to be initalized in a context
# Jupyter notebooks / Python session don't need this
if __name__ == "__main__":

    analysis()
