import ROOT
import os
import sys
import time
from datetime import datetime
import json

ROOT.ROOT.EnableImplicitMT()
ROOT.ROOT.EnableThreadSafety()
from utilsHrare import getMClist, getDATAlist
from utilsHrare import computeWeigths, getMesonFromJson

isGF = True
isVBF = False
isW = False
isZ = False
isZinv = False
isVBFlow = False

BARRELphotons = "Photon_pt>10 and Photon_isScEtaEB and Photon_mvaID_WP90 and Photon_electronVeto"
ENDCAPphotons = "Photon_pt>10 and Photon_isScEtaEE and Photon_mvaID_WP90 and Photon_electronVeto"
#GOODphotons = ""
if(isVBF): GOODphotons = "{} && Photon_pt>20".format(BARRELphotons)
#if(isW or isZ or isZinv or isGF or isVBF): GOODphotons = "{0} or {1}".format(BARRELphotons,ENDCAPphotons)
GOODphotons = "{0} or {1}".format(BARRELphotons,ENDCAPphotons)

with open("/home/submit/mariadlf/Hrare/CMSSW_10_6_27/src/Hrare/analysis/config/selection.json") as jsonFile:
    jsonObject = json.load(jsonFile)
    jsonFile.close()
mesons = jsonObject['mesons']

if(isVBF): GOODRHO = "{}".format(getMesonFromJson(mesons, "isVBF", "isRhoCat"))
if(isVBFlow): GOODRHO = "{}".format(getMesonFromJson(mesons, "isVBFlow" , "isRhoCat"))
if(isZinv or isGF): GOODRHO = "{}".format(getMesonFromJson(mesons, "isZinv", "isRhoCat"))
if(isW or isZ): GOODRHO = "{}".format(getMesonFromJson(mesons, "VH", "isRhoCat"))

if(isVBF): GOODPHI = "{}".format(getMesonFromJson(mesons, "isVBF", "isPhiCat"))
if(isVBFlow): GOODPHI = "{}".format(getMesonFromJson(mesons, "isVBFlow" , "isPhiCat"))
if(isZinv or isGF): GOODPHI = "{}".format(getMesonFromJson(mesons, "isZinv", "isPhiCat"))
if(isW or isZ): GOODPHI = "{}".format(getMesonFromJson(mesons, "VH", "isPhiCat"))

if(isZinv or isGF): GOODK0Star = "{}".format(getMesonFromJson(mesons, "isZinv", "isK0StarCat"))
if(isZinv or isGF): GOODD0Star = "{}".format(getMesonFromJson(mesons, "isZinv", "isD0StarCat"))
if(isZinv or isGF): GOODOMEGA = "{}".format(getMesonFromJson(mesons, "isZinv", "isOmegaCat"))

if False:
#    GOODPHI = "(abs(phi_kin_mass-1.02)<0.01 && phi_kin_pt>10 && phi_trk1_pt>5 && phi_trk2_pt>5 && phi_kin_valid && phi_iso > 0.9 && phi_kin_vtx_prob > 0.05)"
#    GOODRHO = "(rho_kin_mass>0.6 && rho_kin_mass<0.95 && rho_kin_pt>20 && rho_trk1_pt>5 && rho_trk2_pt>5 && rho_kin_valid && rho_iso > 0.9 && rho_kin_vtx_prob > 0.05)"
    GOODKS = "(abs(ks_kin_mass-0.498)<0.05 && ks_kin_pt>5 && ks_kin_valid && ks_iso > 0.8)"
    GOODJPSI = "(abs(Jpsi_kin_mass-3.0969)<0.5 && Jpsi_kin_pt>5  && Jpsi_kin_valid && (Jpsi_muon1_isMediumMuon & 1) && (Jpsi_muon2_isMediumMuon & 1))"
#else:
#    GOODPHI = "(abs(phi_kin_mass-1.02)<0.01 && phi_kin_pt>20 && phi_trk1_pt>5 && phi_trk2_pt>5 && phi_kin_valid && phi_iso > 0.8 && phi_kin_massErr<0.2)"
#    GOODRHO = "(rho_kin_mass>0.6 && rho_kin_mass<0.95 && rho_kin_pt>20 && rho_trk1_pt>5 && rho_trk2_pt>5 && rho_kin_valid && rho_iso > 0.9 && rho_kin_vtx_prob > 0.05)"

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

    if True:
        dfbase = (df.Define("goodPhi","{}".format(GOODPHI))
              .Define("goodPhi_pt", "phi_kin_pt[goodPhi]")
              .Define("goodPhi_eta", "phi_kin_eta[goodPhi]")
              .Define("goodPhi_phi", "phi_kin_phi[goodPhi]")
              .Define("goodPhi_mass", "phi_kin_mass[goodPhi]")
              .Define("goodPhi_massErr", "phi_kin_massErr[goodPhi]")
              .Define("index_Phi","genMatchRECO(goodPhi_pt,goodPhi_eta,goodPhi_phi,goodPhi_mass,GenPart_eta,GenPart_phi,GenPart_pdgId,GenPart_genPartIdxMother,333,25)")
              .Define("goodPhiFromHiggs_pt", "(index_Phi[0]!= -1) ? goodPhi_pt[index_Phi[0]] : 0")
              .Define("goodPhiFromHiggs_eta", "(index_Phi[0]!= -1) ? goodPhi_eta[index_Phi[0]] : 0")
              .Define("goodPhiFromHiggs_phi", "(index_Phi[0]!= -1) ? goodPhi_phi[index_Phi[0]] : 0")
              .Define("goodPhiFromHiggs_mass", "(index_Phi[0]!= -1) ? goodPhi_mass[index_Phi[0]] : 0")
              .Define("goodPhiFromHiggs_massErr", "(index_Phi[0]!= -1) ? goodPhi_massErr[index_Phi[0]] : 0")
              .Define("allPhiGenPart_mass", "(index_Phi[1]!= -1) ? GenPart_mass[index_Phi[1]] : 0")
              .Define("allPhiGenPart_pt", "(index_Phi[1]!= -1) ? GenPart_pt[index_Phi[1]] : 0")
              .Define("allPhiGenPart_eta", "(index_Phi[1]!= -1) ? GenPart_eta[index_Phi[1]] : 0")
              .Define("matchedPhiGenPart_mass", "(index_Phi[1]!= -1 and index_Phi[0]!= -1) ? GenPart_mass[index_Phi[1]] : 0")
              .Define("matchedPhiGenPart_pt", "(index_Phi[1]!= -1 and index_Phi[0]!= -1) ? GenPart_pt[index_Phi[1]] : 0")
              .Define("matchedPhiGenPart_eta", "(index_Phi[1]!= -1 and index_Phi[0]!= -1) ? GenPart_eta[index_Phi[1]] : 0")
              ##
              .Define("goodRho","{}".format(GOODRHO))
              .Define("goodRho_pt", "rho_kin_pt[goodRho]")
              .Define("goodRho_eta", "rho_kin_eta[goodRho]")
              .Define("goodRho_phi", "rho_kin_phi[goodRho]")
              .Define("goodRho_mass", "rho_kin_mass[goodRho]")
              .Define("goodRho_massErr", "rho_kin_massErr[goodRho]")
              .Define("index_Rho","genMatchRECO(goodRho_pt,goodRho_eta,goodRho_phi,goodRho_mass,GenPart_eta,GenPart_phi,GenPart_pdgId,GenPart_genPartIdxMother,113,25)")
              .Define("goodRhoFromHiggs_pt", "(index_Rho[0]!= -1) ? goodRho_pt[index_Rho[0]] : 0")
              .Define("goodRhoFromHiggs_eta", "(index_Rho[0]!= -1) ? goodRho_eta[index_Rho[0]] : 0")
              .Define("goodRhoFromHiggs_phi", "(index_Rho[0]!= -1) ? goodRho_phi[index_Rho[0]] : 0")
              .Define("goodRhoFromHiggs_mass", "(index_Rho[0]!= -1) ? goodRho_mass[index_Rho[0]] : 0")
              .Define("goodRhoFromHiggs_massErr", "(index_Rho[0]!= -1) ? goodRho_massErr[index_Rho[0]] : 0")
              .Define("matchedRhoGenPart_mass", "(index_Rho[1]!= -1) ? GenPart_mass[index_Rho[1]] : 0")
              .Define("matchedRhoGenPart_pt", "(index_Rho[1]!= -1) ? GenPart_pt[index_Rho[1]] : 0")
              .Define("matchedRhoGenPart_eta", "(index_Rho[1]!= -1) ? GenPart_eta[index_Rho[1]] : 0")
              .Define("allRhoGenPart_mass", "(index_Rho[1]!= -1 and index_Rho[0]!= -1) ? GenPart_mass[index_Rho[1]] : 0")
              .Define("allRhoGenPart_pt", "(index_Rho[1]!= -1 and index_Rho[0]!= -1) ? GenPart_pt[index_Rho[1]] : 0")
              .Define("allRhoGenPart_eta", "(index_Rho[1]!= -1 and index_Rho[0]!= -1) ? GenPart_eta[index_Rho[1]] : 0")
#              )
              ##
#    if False:
#        dfADD_ = (dfbase.Define("goodK0Star","{}".format(GOODK0Star))
              .Define("goodK0Star","{}".format(GOODK0Star))
              .Define("goodK0Star_pt", "K0Star_kin_pt[goodK0Star]")
              .Define("goodK0Star_eta", "K0Star_kin_eta[goodK0Star]")
              .Define("goodK0Star_phi", "K0Star_kin_phi[goodK0Star]")
              .Define("goodK0Star_mass", "K0Star_kin_mass[goodK0Star]")
              .Define("goodK0Star_massErr", "K0Star_kin_massErr[goodK0Star]")
              .Define("index_K0Star","genMatchRECO(goodK0Star_pt,goodK0Star_eta,goodK0Star_phi,goodK0Star_mass,GenPart_eta,GenPart_phi,GenPart_pdgId,GenPart_genPartIdxMother,313,25)")
              .Define("goodK0StarFromHiggs_pt", "(index_K0Star[0]!= -1) ? goodK0Star_pt[index_K0Star[0]] : 0")
              .Define("goodK0StarFromHiggs_eta", "(index_K0Star[0]!= -1) ? goodK0Star_eta[index_K0Star[0]] : 0")
              .Define("goodK0StarFromHiggs_phi", "(index_K0Star[0]!= -1) ? goodK0Star_phi[index_K0Star[0]] : 0")
              .Define("goodK0StarFromHiggs_mass", "(index_K0Star[0]!= -1) ? goodK0Star_mass[index_K0Star[0]] : 0")
              .Define("goodK0StarFromHiggs_massErr", "(index_K0Star[0]!= -1) ? goodK0Star_massErr[index_K0Star[0]] : 0")
              .Define("matchedK0StarGenPart_mass", "(index_K0Star[1]!= -1) ? GenPart_mass[index_K0Star[1]] : 0")
              .Define("matchedK0StarGenPart_pt", "(index_K0Star[1]!= -1) ? GenPart_pt[index_K0Star[1]] : 0")
              .Define("matchedK0StarGenPart_eta", "(index_K0Star[1]!= -1) ? GenPart_eta[index_K0Star[1]] : 0")
              .Define("allK0StarGenPart_mass", "(index_K0Star[1]!= -1 and index_K0Star[0]!= -1) ? GenPart_mass[index_K0Star[1]] : 0")
              .Define("allK0StarGenPart_pt", "(index_K0Star[1]!= -1 and index_K0Star[0]!= -1) ? GenPart_pt[index_K0Star[1]] : 0")
              .Define("allK0StarGenPart_eta", "(index_K0Star[1]!= -1 and index_K0Star[0]!= -1) ? GenPart_eta[index_K0Star[1]] : 0")
              ##
              .Define("goodD0Star","{}".format(GOODD0Star))
              .Define("goodD0Star_mass", "d0_d0Star_3body_mass[goodD0Star]")
              .Define("goodD0Star_Nphotons", "d0_d0Star_Nphotons[goodD0Star]")
              .Define("goodD0_pt", "d0_kin_pt[goodD0Star]")
              .Define("goodD0_eta", "d0_kin_eta[goodD0Star]")
              .Define("goodD0_phi", "d0_kin_phi[goodD0Star]")
              .Define("goodD0_mass", "d0_kin_mass[goodD0Star]")
              .Define("goodD0_massErr", "d0_kin_massErr[goodD0Star]")
              .Define("index_D0Star","genMatchRECO(goodD0_pt,goodD0_eta,goodD0_phi,goodD0_mass,GenPart_eta,GenPart_phi,GenPart_pdgId,GenPart_genPartIdxMother,421,423)")
              .Define("goodD0FromHiggs_pt", "(index_D0Star[0]!= -1) ? goodD0_pt[index_D0Star[0]] : 0")
              .Define("goodD0FromHiggs_eta", "(index_D0Star[0]!= -1) ? goodD0_eta[index_D0Star[0]] : 0")
              .Define("goodD0FromHiggs_phi", "(index_D0Star[0]!= -1) ? goodD0_phi[index_D0Star[0]] : 0")
              .Define("goodD0FromHiggs_mass", "(index_D0Star[0]!= -1) ? goodD0_mass[index_D0Star[0]] : 0")
              .Define("goodD0FromHiggs_massErr", "(index_D0Star[0]!= -1) ? goodD0_massErr[index_D0Star[0]] : 0")
              .Define("matchedD0GenPart_mass", "(index_D0Star[1]!= -1) ? GenPart_mass[index_D0Star[1]] : 0")
              .Define("matchedD0GenPart_pt", "(index_D0Star[1]!= -1) ? GenPart_pt[index_D0Star[1]] : 0")
              .Define("matchedD0GenPart_eta", "(index_D0Star[1]!= -1) ? GenPart_eta[index_D0Star[1]] : 0")
              .Define("allD0GenPart_mass", "(index_D0Star[1]!= -1 and index_D0Star[0]!= -1) ? GenPart_mass[index_D0Star[1]] : 0")
              .Define("allD0GenPart_pt", "(index_D0Star[1]!= -1 and index_D0Star[0]!= -1) ? GenPart_pt[index_D0Star[1]] : 0")
              .Define("allD0GenPart_eta", "(index_D0Star[1]!= -1 and index_D0Star[0]!= -1) ? GenPart_eta[index_D0Star[1]] : 0")
              ##
              .Define("goodOmega","{}".format(GOODOMEGA))
              .Define("goodOmega_pt", "omega_kin_pt[goodOmega]")
              .Define("goodOmega_eta", "omega_kin_eta[goodOmega]")
              .Define("goodOmega_phi", "omega_kin_phi[goodOmega]")
              .Define("goodOmega_mass", "omega_kin_mass[goodOmega]")
              .Define("goodOmega_threemass", "omega_threemass[goodOmega]")
              .Define("goodOmega_massErr", "omega_kin_massErr[goodOmega]")
              .Define("index_Omega","genMatchRECO(goodOmega_pt,goodOmega_eta,goodOmega_phi,goodOmega_mass,GenPart_eta,GenPart_phi,GenPart_pdgId,GenPart_genPartIdxMother,223,25)")
              .Define("goodOmegaFromHiggs_pt", "(index_Omega[0]!= -1) ? goodOmega_pt[index_Omega[0]] : 0")
              .Define("goodOmegaFromHiggs_eta", "(index_Omega[0]!= -1) ? goodOmega_eta[index_Omega[0]] : 0")
              .Define("goodOmegaFromHiggs_phi", "(index_Omega[0]!= -1) ? goodOmega_phi[index_Omega[0]] : 0")
              .Define("goodOmegaFromHiggs_mass", "(index_Omega[0]!= -1) ? goodOmega_mass[index_Omega[0]] : 0")
              .Define("goodOmegaFromHiggs_massErr", "(index_Omega[0]!= -1) ? goodOmega_massErr[index_Omega[0]] : 0")
              .Define("goodOmegaFromHiggs_threemass", "(index_Omega[0]!= -1) ? goodOmega_threemass[index_Omega[0]] : 0")
              .Define("matchedOmegaGenPart_mass", "(index_Omega[1]!= -1) ? GenPart_mass[index_Omega[1]] : 0")
              .Define("matchedOmegaGenPart_pt", "(index_Omega[1]!= -1) ? GenPart_pt[index_Omega[1]] : 0")
              .Define("matchedOmegaGenPart_eta", "(index_Omega[1]!= -1) ? GenPart_eta[index_Omega[1]] : 0")
              .Define("allOmegaGenPart_mass", "(index_Omega[1]!= -1 and index_Omega[0]!= -1) ? GenPart_mass[index_Omega[1]] : 0")
              .Define("allOmegaGenPart_pt", "(index_Omega[1]!= -1 and index_Omega[0]!= -1) ? GenPart_pt[index_Omega[1]] : 0")
              .Define("allOmegaGenPart_eta", "(index_Omega[1]!= -1 and index_Omega[0]!= -1) ? GenPart_eta[index_Omega[1]] : 0")
            )
              ##
    if False:
        dfADD_ = (dfbase.Define("goodKs","{}".format(GOODKS))
                  .Define("goodKs_pt", "ks_kin_pt[goodKs]")
                  .Define("goodKs_eta", "ks_kin_eta[goodKs]")
                  .Define("goodKs_phi", "ks_kin_phi[goodKs]")
                  .Define("goodKs_mass", "ks_kin_mass[goodKs]")
                  .Define("goodKs_massErr", "ks_kin_massErr[goodKs]")
                  .Define("goodKs_trk1_pt", "ks_trk1_pt[goodKs]")
                  .Define("goodKs_trk2_pt", "ks_trk2_pt[goodKs]")
                  .Define("index_Ks","genMatchRECO(goodKs_pt,goodKs_eta,goodKs_phi,goodKs_mass,GenPart_eta,GenPart_phi,GenPart_pdgId,GenPart_genPartIdxMother,310,333)")
                  .Define("goodKsFromHiggs_pt", "(index_Ks[0]!= -1) ? goodKs_pt[index_Ks[0]] : 0")
                  .Define("goodKsFromHiggs_eta", "(index_Ks[0]!= -1) ? goodKs_eta[index_Ks[0]] : 0")
                  .Define("goodKsFromHiggs_phi", "(index_Ks[0]!= -1) ? goodKs_phi[index_Ks[0]] : 0")
                  .Define("goodKsFromHiggs_mass", "(index_Ks[0]!= -1) ? goodKs_mass[index_Ks[0]] : 0")
                  .Define("goodKsFromHiggs_massErr", "(index_Ks[0]!= -1) ? goodKs_massErr[index_Ks[0]] : 0")
                  .Define("matchedKsGenPart_mass", "(index_Ks[1]!= -1) ? GenPart_mass[index_Ks[1]] : 0")
                  ##
                  .Define("goodJPsi","{}".format(GOODJPSI))
                  .Define("goodJPsi_pt", "Jpsi_kin_pt[goodJPsi]")
                  .Define("goodJPsi_eta", "Jpsi_kin_eta[goodJPsi]")
                  .Define("goodJPsi_phi", "Jpsi_kin_phi[goodJPsi]")
                  .Define("goodJPsi_mass", "Jpsi_kin_mass[goodJPsi]")
                  .Define("index_JPsi","genMatchRECO(goodJPsi_pt,goodJPsi_eta,goodJPsi_phi,goodJPsi_mass,GenPart_eta,GenPart_phi,GenPart_pdgId,GenPart_genPartIdxMother,443,25)")
                  .Define("goodJPsiFromHiggs_mass", "(index_JPsi[0]!= -1) ? goodJPsi_mass[index_JPsi[0]] : 0")
                  #
                  .Define("dPhiGammaPhimatched","(index_PH[0]!= -1 and index_Phi[0]!= -1) ? deltaPhi(goodPhotons_phi[index_PH[0]], goodPhi_phi[index_Phi[0]]):-999")
                  .Define("dEtaGammaPhimatched","(index_PH[0]!= -1 and index_Phi[0]!= -1) ? abs(goodPhotons_eta[index_PH[0]] - goodPhi_eta[index_Phi[0]]):-999")
                  .Define("dPhiMETPhimatched","(index_Phi[0]!= -1) ? deltaPhi(MET_phi, goodPhi_phi[index_Phi[0]]):-999")
                  .Define("dPhiMETRhomatched","(index_Rho[0]!= -1) ? deltaPhi(MET_phi, goodRho_phi[index_Rho[0]]):-999")
                  .Define("dPhiMETKsmatched","(index_Ks[0]!= -1) ? deltaPhi(MET_phi, goodKs_phi[index_Ks[0]]):-999")
                  .Define("genKS","abs(GenPart_pdgId)=={0} && checkMother(GenPart_pdgId,GenPart_genPartIdxMother,{0},{1})".format(310,333))
                  .Define("genKL","abs(GenPart_pdgId)=={0} && checkMother(GenPart_pdgId,GenPart_genPartIdxMother,{0},{0})".format(130,333))
                  .Define("genKS_pt","GenPart_pt[genKS]")
                  .Define("genKL_pt","GenPart_pt[genKL]")
                  .Define("genKS_eta","GenPart_eta[genKS]")
                  .Define("genKL_eta","GenPart_eta[genKL]")
                  .Define("genKS_phi","GenPart_phi[genKS]")
                  .Define("genKL_phi","GenPart_phi[genKL]")
                  #              .Define("dRkSkL","DeltaR(genKS_eta,genKL_eta,genKS_phi,genKL_phi)")
                  #              .Define("dRkSkL","deltaR(genKS_eta[0],genKS_phi[0],genKL_eta[0],genKL_phi[0])")
                  #              float deltaR(float eta1, float phi1, float eta2, float phi2)
                  .Define("genJet_KSmask","cleaningJetFromMeson(GenJet_eta,GenJet_phi,genKS_eta[0],genKS_phi[0])") # when false, we have a good match
                  #              .Define("genGetMATCH_pt", "(GenJet_pt>100 and (genJet_KSmask == false)) ? GenJet_pt : 0")
                  #              .Define("genGetMATCH_eta", "(GenJet_pt>100 and (genJet_KSmask == false)) ? GenJet_eta : 0")
                  #              .Define("genGetMATCH_phi", "(GenJet_pt>100 and (genJet_KSmask == false)) ? GenJet_phi : 0")
                  #              .Define("HcandMassKenJetKS", "Minv2(genGetMATCH_pt[0],genGetMATCH_eta[0],genGetMATCH_phi[0],genGetMATCH_mass[0],goodPhotons_pt[index_PH[0]],goodPhotons_eta[index_PH[0]],goodPhotons_phi[index_PH[0]]).first")
                  #
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
            "nGenPart",
            "GenPart_mass",
            "w",
            "mc",
            #
            "goodPhotons_pt",
            "goodPhotons_eta",
            "goodPhotons_phi",
            "index_PH",
            "goodPhotonsFromHiggs_pt",
            "goodPhotonsFromHiggs_eta",
            "goodPhotonsFromHiggs_phi",
    ]:
            branchList.push_back(branchName)
            #

    if True:
        for branchName in [
            #
            "goodPhi",
            "goodPhi_mass",
            "goodPhi_massErr",
            "goodPhi_pt",
            "goodPhi_eta",
            "goodPhi_phi",
            "index_Phi",
            "goodPhiFromHiggs_pt",
            "goodPhiFromHiggs_eta",
            "goodPhiFromHiggs_phi",
            "goodPhiFromHiggs_mass",
            "goodPhiFromHiggs_massErr",
            "matchedPhiGenPart_mass",
            "matchedPhiGenPart_pt",
            "matchedPhiGenPart_eta",
            "allPhiGenPart_mass",
            "allPhiGenPart_pt",
            "allPhiGenPart_eta",
            #
            "goodRho",
            "goodRho_pt",
            "goodRho_eta",
            "goodRho_phi",
            "goodRho_mass",
            "goodRho_massErr",
            "index_Rho",
            "goodRhoFromHiggs_pt",
            "goodRhoFromHiggs_eta",
            "goodRhoFromHiggs_phi",
            "goodRhoFromHiggs_mass",
            "goodRhoFromHiggs_massErr",
            "matchedRhoGenPart_mass",
            "matchedRhoGenPart_pt",
            "matchedRhoGenPart_eta",
            "allRhoGenPart_mass",
            "allRhoGenPart_pt",
            "allRhoGenPart_eta",
    ]:
            branchList.push_back(branchName)
            #

    if True:
        for branchName in [
                "goodK0Star",
                "goodK0Star_pt",
                "goodK0Star_eta",
                "goodK0Star_phi",
                "goodK0Star_mass",
                "goodK0Star_massErr",
                "index_K0Star",
                "goodK0StarFromHiggs_pt",
                "goodK0StarFromHiggs_eta",
                "goodK0StarFromHiggs_phi",
                "goodK0StarFromHiggs_mass",
                "goodK0StarFromHiggs_massErr",
                "matchedK0StarGenPart_mass",
                "matchedK0StarGenPart_pt",
                "matchedK0StarGenPart_eta",
                "allK0StarGenPart_mass",
                "allK0StarGenPart_pt",
                "allK0StarGenPart_eta",
                ##
                "goodD0Star",
                "goodD0Star_mass",
                "goodD0_pt",
                "goodD0_eta",
                "goodD0_phi",
                "goodD0_mass",
                "goodD0_massErr",
                "goodD0Star_Nphotons",
                "index_D0Star",
                "goodD0FromHiggs_pt",
                "goodD0FromHiggs_eta",
                "goodD0FromHiggs_phi",
                "goodD0FromHiggs_mass",
                "goodD0FromHiggs_massErr",
                "matchedD0GenPart_mass",
                "matchedD0GenPart_pt",
                "matchedD0GenPart_eta",
                "allD0GenPart_mass",
                "allD0GenPart_pt",
                "allD0GenPart_eta",
                ##
                "goodOmega",
                "goodOmega_pt",
                "goodOmega_eta",
                "goodOmega_phi",
                "goodOmega_mass",
                "goodOmega_massErr",
                "goodOmega_threemass",
                "index_Omega",
                "goodOmegaFromHiggs_pt",
                "goodOmegaFromHiggs_eta",
                "goodOmegaFromHiggs_phi",
                "goodOmegaFromHiggs_mass",
                "goodOmegaFromHiggs_massErr",
                "goodOmegaFromHiggs_threemass",
                "matchedOmegaGenPart_mass",
                "matchedOmegaGenPart_pt",
                "matchedOmegaGenPart_eta",
                "allOmegaGenPart_mass",
                "allOmegaGenPart_pt",
                "allOmegaGenPart_eta",
        ]:
            branchList.push_back(branchName)

    if False:
        for branchName in [
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
                "goodKsFromHiggs_massErr",
                "matchedKsGenPart_mass",
                #
                "goodJPsi",
                "goodJPsi_pt",
                "goodJPsi_eta",
                "goodJPsi_phi",
                "goodJPsi_mass",
                "index_JPsi",
                "goodJPsiFromHiggs_mass",
                #
                "MET_pt",
                "MET_phi",
                "dPhiMETPhimatched",
                "dPhiMETRhomatched",
                "dPhiMETKsmatched",
                #            "dRkSkL",
                "genKS_pt",
                "genKL_pt",
                #            "HcandMassKenJetKS"
                #
                "dPhiGammaPhimatched",
                "dEtaGammaPhimatched",
        ]:
            branchList.push_back(branchName)


    outputFile = "/work/submit/mariadlf/genStudies/miniTreeForGNN_mc%d"%mc+".root"
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
#    readMCSample(2018,12)
#    readMCSample(2018,13)

    readMCSample(2018,1037)    #ggH-Omega
    readMCSample(2018,1039)    #ggH-D0Star
    readMCSample(2018,1038)    #ggH-K0Star

    exit()

    readMCSample(2018,1010)    #vbs-phi
    readMCSample(2018,1020)    #vbs-rho
    
    readMCSample(2018,1013)    #zh-phi
    readMCSample(2018,1023)    #zh-rho

    readMCSample(2018,1011)    #wh-phi
    readMCSample(2018,1021)    #wh-rho

    readMCSample(2018,1017)    #ggH-phi
    readMCSample(2018,1027)    #ggH-rho

    #
#    readMCSample(2018,1012)    #wh-phi
#    readMCSample(2018,1022)    #wh-rho
    #
#    readMCSample(2018,1018)    #vbs-kskl
#    readMCSample(2018,1030)    #zh Jpsi
