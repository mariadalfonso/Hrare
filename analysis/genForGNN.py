import ROOT
import os
import sys

ROOT.ROOT.EnableImplicitMT()
ROOT.ROOT.EnableThreadSafety()
from utilsHrare import getMClist
from utilsHrare import plot, computeWeigths

isVBF = True
isW = False
isZ = False

BARRELphotons = "Photon_pt>20 and Photon_isScEtaEB and (Photon_cutBased & 2) and Photon_electronVeto"
ENDCAPphotons = "Photon_pt>20 and Photon_isScEtaEE and (Photon_cutBased & 2) and Photon_electronVeto"
GOODphotons = ""
if(isVBF): GOODphotons = "{} && Photon_pt>75".format(BARRELphotons)
if(isW or isZ): GOODphotons = "{0} or {1}".format(BARRELphotons,ENDCAPphotons)

if isVBF:
    GOODPHI = "(abs(phi_kin_mass-1.02)<0.01 && phi_kin_pt>5 && phi_trk1_pt>5 && phi_trk2_pt>5 && phi_kin_valid && phi_iso > 0.8)"
    GOODRHO = "(rho_kin_mass>0.6 && rho_kin_mass<0.95 && rho_kin_pt>20 && rho_trk1_pt>5 && rho_trk2_pt>5 && rho_kin_valid && rho_iso > 0.9)"
    GOODKS = "(abs(ks_kin_mass-0.498)<0.05 && ks_kin_pt>5 && ks_trk1_pt>5 && ks_trk2_pt>5 && ks_kin_valid)"
else:
    GOODPHI = "(abs(phi_kin_mass-1.02)<0.01 && phi_kin_pt>10 && phi_trk1_pt>10 && phi_trk2_pt>10 && phi_kin_valid && phi_iso > 0.8)"
    GOODRHO = "(rho_kin_mass>0.6 && rho_kin_mass<0.95 && rho_kin_pt>10 && rho_trk1_pt>5 && rho_trk2_pt>5 && rho_kin_valid && rho_iso > 0.8)"

def dfGammaCand(df):

    dfOBJ = (df.Filter("nPhoton>0")
             .Define("goodPhotons", "{}".format(GOODphotons))
             .Filter("Sum(goodPhotons)>0", "At least one good Photon")
             .Define("goodPhotons_pt", "Photon_pt[goodPhotons]")
             .Define("goodPhotons_eta", "Photon_eta[goodPhotons]")
             .Define("goodPhotons_phi", "Photon_phi[goodPhotons]")
             .Define("goodPhotons_mass", "Photon_mass[goodPhotons]")             
             .Define("index_PH","genMatchRECO(goodPhotons_pt,goodPhotons_eta,goodPhotons_phi,goodPhotons_mass,GenPart_eta,GenPart_phi,GenPart_pdgId,GenPart_genPartIdxMother,22,25)")
             .Define("goodPhotonsFromHiggs_pt", "(index_PH[0]!= -1) ? goodPhotons_pt[index_PH[0]] : 0 ")
             .Define("goodPhotonsFromHiggs_eta", "(index_PH[0]!= -1) ? goodPhotons_eta[index_PH[0]] : 0 ")
             .Define("goodPhotonsFromHiggs_phi", "(index_PH[0]!= -1) ? goodPhotons_phi[index_PH[0]] : 0")
             )
    return dfOBJ
    

def dfMesonCand(df):

    dfbase = (df.Filter("nphi>0 or nrho>0")
              .Define("goodPhi","{}".format(GOODPHI))
              .Define("goodPhi_pt", "phi_kin_pt[goodPhi]")
              .Define("goodPhi_eta", "phi_kin_eta[goodPhi]")
              .Define("goodPhi_phi", "phi_kin_phi[goodPhi]")
              .Define("goodPhi_mass", "phi_kin_mass[goodPhi]")
              .Define("index_Phi","genMatchRECO(goodPhi_pt,goodPhi_eta,goodPhi_phi,goodPhi_mass,GenPart_eta,GenPart_phi,GenPart_pdgId,GenPart_genPartIdxMother,333,25)")
              .Define("goodPhiFromHiggs_pt", "(index_Phi[0]!= -1) ? goodPhi_pt[index_Phi[0]] : 0")
              .Define("goodPhiFromHiggs_eta", "(index_Phi[0]!= -1) ? goodPhi_eta[index_Phi[0]] : 0")
              .Define("goodPhiFromHiggs_phi", "(index_Phi[0]!= -1) ? goodPhi_phi[index_Phi[0]] : 0")
              .Define("goodPhiFromHiggs_mass", "(index_Phi[0]!= -1) ? goodPhi_mass[index_Phi[0]] : 0")                            
              ##
              .Define("goodRho","{}".format(GOODRHO))
              .Define("goodRho_pt", "rho_kin_pt[goodRho]")
              .Define("goodRho_eta", "rho_kin_eta[goodRho]")
              .Define("goodRho_phi", "rho_kin_phi[goodRho]")
              .Define("goodRho_mass", "rho_kin_mass[goodRho]")
              .Define("index_Rho","genMatchRECO(goodRho_pt,goodRho_eta,goodRho_phi,goodRho_mass,GenPart_eta,GenPart_phi,GenPart_pdgId,GenPart_genPartIdxMother,113,25)")
              .Define("goodRhoFromHiggs_pt", "(index_Rho[0]!= -1) ? goodRho_pt[index_Rho[0]] : 0")
              .Define("goodRhoFromHiggs_eta", "(index_Rho[0]!= -1) ? goodRho_eta[index_Rho[0]] : 0")
              .Define("goodRhoFromHiggs_phi", "(index_Rho[0]!= -1) ? goodRho_phi[index_Rho[0]] : 0")
              .Define("goodRhoFromHiggs_mass", "(index_Rho[0]!= -1) ? goodRho_mass[index_Rho[0]] : 0")                                          
              ##
              .Define("goodKs","{}".format(GOODKS))
              .Define("goodKs_pt", "ks_kin_pt[goodKs]")
              .Define("goodKs_eta", "ks_kin_eta[goodKs]")
              .Define("goodKs_phi", "ks_kin_phi[goodKs]")
              .Define("goodKs_mass", "ks_kin_mass[goodKs]")
              .Define("index_Ks","genMatchRECO(goodKs_pt,goodKs_eta,goodKs_phi,goodKs_mass,GenPart_eta,GenPart_phi,GenPart_pdgId,GenPart_genPartIdxMother,310,333)")
              .Define("goodKsFromHiggs_pt", "(index_Ks[0]!= -1) ? goodKs_pt[index_Ks[0]] : 0")
              .Define("goodKsFromHiggs_eta", "(index_Ks[0]!= -1) ? goodKs_eta[index_Ks[0]] : 0")
              .Define("goodKsFromHiggs_phi", "(index_Ks[0]!= -1) ? goodKs_phi[index_Ks[0]] : 0")
              .Define("goodKsFromHiggs_mass", "(index_Ks[0]!= -1) ? goodKs_mass[index_Ks[0]] : 0")
              )
    return dfbase


def analysis(df,year,mc,sumw,isData):

    lumi=1.
    weight = "{0}".format(1.)
    if mc>0: weight = "{0}*genWeight*{1}".format(lumi,sumw)

    dfOBJ= dfGammaCand(df)
    dfbase = dfMesonCand(dfOBJ)
    dfFINAL = (dfbase.Define("w","{}".format(weight))
               .Define("mc","{}".format(mc))
    )

    branchList = ROOT.vector('string')()
    for branchName in [
            "goodPhotons_pt",
            "goodPhotons_eta",
            "goodPhotons_phi",
            "index_PH",
            "goodPhotonsFromHiggs_pt",
            "goodPhotonsFromHiggs_eta",
            "goodPhotonsFromHiggs_phi",
            #
            "goodPhi",
            "goodPhi_mass",
            "goodPhi_pt",
            "goodPhi_eta",
            "goodPhi_phi",
            "index_Phi",
            "goodPhiFromHiggs_pt",
            "goodPhiFromHiggs_eta",
            "goodPhiFromHiggs_phi",
            "goodPhiFromHiggs_mass",
            #
            "goodRho",
            "goodRho_pt",
            "goodRho_eta",
            "goodRho_phi",
            "goodRho_mass",
            "index_Rho",
            "goodRhoFromHiggs_pt",
            "goodRhoFromHiggs_eta",
            "goodRhoFromHiggs_phi",
            "goodRhoFromHiggs_mass",
            #
            "goodKs",
            "goodKs_pt",
            "goodKs_eta",
            "goodKs_phi",
            "goodKs_mass",
            "index_Ks",
            "goodKsFromHiggs_pt",
            "goodKsFromHiggs_eta",
            "goodKsFromHiggs_phi",
            "goodKsFromHiggs_mass",
            #
            "w",
            "mc",
    ]:
        branchList.push_back(branchName)


    outputFile = "miniTreeForGNN_mc%d"%mc+".root"
    print(outputFile)
                    
    if True:            
        snapshot_tdf = dfFINAL.Snapshot("events", outputFile, branchList)
        print("snapshot_tdf DONE")
        
def readMCSample(year,sampleNOW):

    files = getMClist(year,sampleNOW)
    print(len(files))
    df = ROOT.RDataFrame("Events", files)

    w = computeWeigths(df, files, sampleNOW, year, True)

    analysis(df,year,sampleNOW,w,"false")


if __name__ == "__main__":

#   to run: python3 -i genForGNN.py
    readMCSample(2018,1010)    #vbs-phi
    readMCSample(2018,1020)    #vbs-rho
    readMCSample(2018,1013)    #vbs-kskl
    
