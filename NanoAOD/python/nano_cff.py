from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_cff import *


def nanoAOD_customizeMesons(process):
    process.load('Hrare.NanoAOD.MesonsReco_cff')

    # remove: slow and unused
    (run2_nanoAOD_106Xv1 | run2_nanoAOD_106Xv2).toReplaceWith(process.nanoSequenceMC, process.nanoSequenceMC.copyAndExclude([btagWeightTable]))

    process.genParticleTable.src = "prunedGenParticles"

    # Data 
#    process.nanoSequence   = cms.Sequence(process.slimmedMuons + process.nanoSequence + process.V0Sequence + process.V0ForMuonFakeTables)
    process.nanoSequence   = cms.Sequence(process.nanoSequence + process.V0Sequence + process.V0Tables)
    # MC
#    process.nanoSequenceMC = cms.Sequence(process.slimmedMuons + process.nanoSequenceMC + process.V0ForMuonFakeMcSequence + process.V0ForMuonFakeMcTables)
    process.nanoSequenceMC = cms.Sequence(process.nanoSequenceMC + process.V0McSequence + process.V0McTables)
    process.muonTable.variables.softMva = Var("softMvaValue()",float,doc="soft MVA ID score",precision=6)
    return process
