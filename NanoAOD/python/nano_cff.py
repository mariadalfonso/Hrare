from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_cff import *


def nanoAOD_customizeMesons_Run3(process):
    process.load('Hrare.NanoAOD.MesonsReco_cff')
    process.load('Hrare.NanoAOD.DiMuonReco_cff')


    if hasattr(process.triggerObjectTable, 'selections') and hasattr(process.triggerObjectTable.selections, 'id') and process.triggerObjectTable.selections.id.value() == 22:
        process.triggerObjectTable.selections.name = cms.string("Photon")
        process.triggerObjectTable.selections.sel = cms.string("type(92) && pt > 20 && coll('hltEgammaCandidates')")
        process.triggerObjectTable.selections.qualityBits = cms.string(
            "filter('hltEG33L1EG26HEFilter') + " \
            "2*filter('hltEG50HEFilter') + " \
            "4*filter('hltEG75HEFilter') + " \
            "8*filter('hltEG90HEFilter') + " \
            "16*filter('hltEG120HEFilter') + " \
            "32*filter('hltEG150HEFilter') + " \
            "64*filter('hltEG175HEFilter') + " \
            "128*filter('hltEG200HEFilter') + " \
            "256*filter('hltHtEcal800') + " \
            "512*filter('hltEG110EBTightIDTightIsoTrackIsoFilter') + " \
            "1024*filter('hltEG120EBTightIDTightIsoTrackIsoFilter')+ " \
            "2048*filter('hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter')")
        process.triggerObjectTable.selections.qualityBitsDoc = cms.string("Single Photon filters: 1 = hltEG33L1EG26HEFilter, 2 = hltEG50HEFilter, 4 = hltEG75HEFilter, 8 = hltEG90HEFilter, 16 = hltEG120HEFilter, 32 = hltEG150HEFilter, 64 = hltEG175HEFilter, 128 = hltEG200HEFilter, 256 = hltHtEcal800, 512 = hltEG110EBTightIDTightIsoTrackIsoFilter, 1024 = hltEG120EBTightIDTightIsoTrackIsoFilter, 2048 = 1mu-1photon")


    process.photonTable.variables.x_calo = Var("superCluster().seed().position().x()",float,doc="photon supercluster position on calorimeter, x coordinate (cm)",precision=10),
    process.photonTable.variables.y_calo = Var("superCluster().seed().position().y()",float,doc="photon supercluster position on calorimeter, y coordinate (cm)",precision=10),
    process.photonTable.variables.z_calo = Var("superCluster().seed().position().z()",float,doc="photon supercluster position on calorimeter, z coordinate (cm)",precision=10),

    process.genParticleTable.src = "prunedGenParticles"
    process.genParticleTable.variables.pt=Var("pt",  float, precision=-1)
    process.genParticleTable.variables.phi=Var("phi",  float, precision=-1)
    process.genParticleTable.variables.eta=Var("eta",  float, precision=-1)
    process.genParticleTable.variables.mass=Var("mass",  float, precision=-1)

    finalGenParticles.select +=[
        "keep (4 <= abs(pdgId) <= 5) && statusFlags().isLastCopy()", # BTV: keep b/c quarks in their last copy
    ]

    # Data
#    process.nanoSequence   = cms.Sequence(process.slimmedMuons + process.nanoSequence + process.V0Sequence + process.V0ForMuonFakeTables)
    process.nanoSequence   = cms.Sequence(process.nanoSequence + process.V0Sequence + process.V0Tables + process.DiMuProdSequence + process.DiMuTables )
    # MC
#    process.nanoSequenceMC = cms.Sequence(process.slimmedMuons + process.nanoSequenceMC + process.V0ForMuonFakeMcSequence + process.V0ForMuonFakeMcTables)
    process.nanoSequenceMC = cms.Sequence(process.nanoSequenceMC + process.V0McSequence + process.V0McTables + process.DiMuProdMcSequence + process.DiMuMcTables)
    process.muonTable.variables.softMva = Var("softMvaValue()",float,doc="soft MVA ID score",precision=6)
    return process

def nanoAOD_customizeMesons(process):
    process.load('Hrare.NanoAOD.MesonsReco_cff')
    process.load('Hrare.NanoAOD.DiMuonReco_cff')

    # remove: slow and unused
    (run2_nanoAOD_106Xv1 | run2_nanoAOD_106Xv2).toReplaceWith(process.nanoSequenceMC, process.nanoSequenceMC.copyAndExclude([btagWeightTable]))

    (run2_nanoAOD_106Xv1 | run2_nanoAOD_106Xv2).toModify(photonTable.variables,
                                                         x_calo = Var("superCluster().seed().position().x()",float,doc="photon supercluster position on calorimeter, x coordinate (cm)",precision=10),
                                                         y_calo = Var("superCluster().seed().position().y()",float,doc="photon supercluster position on calorimeter, y coordinate (cm)",precision=10),
                                                         z_calo = Var("superCluster().seed().position().z()",float,doc="photon supercluster position on calorimeter, z coordinate (cm)",precision=10),
                                                         )

    process.genParticleTable.src = "prunedGenParticles"
    (run2_nanoAOD_106Xv1 | run2_nanoAOD_106Xv2).toModify(genParticleTable.variables,
                                                         pt  = Var("pt",  float, precision=-1),
                                                         phi  = Var("phi",  float, precision=-1),
                                                         eta  = Var("eta",  float, precision=-1),
                                                         px  = Var("px",  float, precision=-1),
                                                         py  = Var("py",  float, precision=-1),
                                                         pz  = Var("pz",  float, precision=-1),
                                                         energy  = Var("energy",  float, precision=-1),
                                                         mass = Var("mass",  float, precision=-1),
#                                                         mass = Var("?!((abs(pdgId)>=1 && abs(pdgId)<=5) || (abs(pdgId)>=11 && abs(pdgId)<=16) || pdgId==21 || pdgId==111 || abs(pdgId)==211 || abs(pdgId)==421 || abs(pdgId)==411 || (pdgId==22 && mass<1))?mass:0", float,precision="?((abs(pdgId)==6 || abs(pdgId)>1000000) && statusFlags().isLastCopy())?20:-1",doc="Mass stored for all particles with the exception of quarks (except top), leptons/neutrinos, photons with mass < 1 GeV, gluons, pi0(111), pi+(211), D0(421), and D+(411). For these particles, you can lookup the value from PDG."),
    )

    if triggerObjectTable.selections[1].id.value()==22:
        process.triggerObjectTable.selections[1].name=cms.string("Photon")
        process.triggerObjectTable.selections[1].sel=cms.string("type(92) && pt > 20 && coll('hltEgammaCandidates')")
        process.triggerObjectTable.selections[1].qualityBits=cms.string(
            "filter('hltEG33L1EG26HEFilter') + " \
            "2*filter('hltEG50HEFilter') + " \
            "4*filter('hltEG75HEFilter') + " \
            "8*filter('hltEG90HEFilter') + " \
            "16*filter('hltEG120HEFilter') + " \
            "32*filter('hltEG150HEFilter') + " \
            "64*filter('hltEG175HEFilter') + " \
            "128*filter('hltEG200HEFilter') + " \
            "256*filter('hltHtEcal800') + " \
            "512*filter('hltEG110EBTightIDTightIsoTrackIsoFilter') + " \
            "1024*filter('hltEG120EBTightIDTightIsoTrackIsoFilter')+ " \
            "2048*filter('hltMu17Photon30IsoCaloIdPhotonlegTrackIsoFilter')")
        process.triggerObjectTable.selections[1].qualityBitsDoc = cms.string("Single Photon filters: 1 = hltEG33L1EG26HEFilter, 2 = hltEG50HEFilter, 4 = hltEG75HEFilter, 8 = hltEG90HEFilter, 16 = hltEG120HEFilter, 32 = hltEG150HEFilter, 64 = hltEG175HEFilter, 128 = hltEG200HEFilter, 256 = hltHtEcal800, 512 = hltEG110EBTightIDTightIsoTrackIsoFilter, 1024 = hltEG120EBTightIDTightIsoTrackIsoFilter, 2048 = 1mu-1photon")

    # Data
#    process.nanoSequence   = cms.Sequence(process.slimmedMuons + process.nanoSequence + process.V0Sequence + process.V0ForMuonFakeTables)
    process.nanoSequence   = cms.Sequence(process.nanoSequence + process.V0Sequence + process.V0Tables + process.DiMuProdSequence + process.DiMuTables )
    # MC
#    process.nanoSequenceMC = cms.Sequence(process.slimmedMuons + process.nanoSequenceMC + process.V0ForMuonFakeMcSequence + process.V0ForMuonFakeMcTables)
    process.nanoSequenceMC = cms.Sequence(process.nanoSequenceMC + process.V0McSequence + process.V0McTables + process.DiMuProdMcSequence + process.DiMuMcTables)
    process.muonTable.variables.softMva = Var("softMvaValue()",float,doc="soft MVA ID score",precision=6)
    return process
