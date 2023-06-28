from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_cff import *


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
