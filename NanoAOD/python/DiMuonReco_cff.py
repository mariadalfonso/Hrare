from PhysicsTools.NanoAOD.common_cff import *
import FWCore.ParameterSet.Config as cms

# NOTE: 
#    All instances of FlatTableProducers must end with Table in their
#    names so that their product match keep patterns in the default
#    event content. Otherwise you need to modify outputCommands in
#    NanoAODEDMEventContent or provide a custom event content to the
#    output module

def merge_psets(*argv):
    result = cms.PSet()
    for pset in argv:
        if isinstance(pset, cms._Parameterizable):
            for name in pset.parameters_().keys():
                value = getattr(pset,name)
                type = value.pythonTypeName()
                setattr(result,name,value)
    return result

DiMuProd = cms.EDProducer(
    "DiMuonProducer",
    beamSpot=cms.InputTag("offlineBeamSpot"),
    vertexCollection=cms.InputTag("offlineSlimmedPrimaryVertices"),
    muonCollection = cms.InputTag("linkedObjects","muons"),
    PFCandCollection = cms.InputTag("packedPFCandidates"),
    packedGenParticleCollection = cms.InputTag("packedGenParticles"),
    minMuonPt  = cms.double(3.5),
    maxMuonEta = cms.double(2.4),
    minPionPt  = cms.double(1.0),
    maxPionEta = cms.double(2.4),
    minJpsisMass  = cms.double(2.9), # jpsi mass 3.096
    maxJpsisMass  = cms.double(3.4),
    minJpsisPreselectMass = cms.double(2.5),
    maxJpsisPreselectMass = cms.double(3.8),
    minUpsilonsMass  = cms.double(9.00), # Y(1S) mass 9.460
    maxUpsilonsMass  = cms.double(10.00),
    minUpsilonsPreselectMass = cms.double(8.),
    maxUpsilonsPreselectMass = cms.double(11.),
    maxTwoTrackDOCA = cms.double(0.1),
    maxLxy = cms.double(999),
    minSigLxy = cms.double(5),
    minVtxProb = cms.double(0.001),
    minCosAlpha = cms.double(0.9),
    minDisplaceTrackSignificance = cms.double(1),
    isMC = cms.bool(False)
)

DiMuProdMC = DiMuProd.clone( isMC = cms.bool(True) )

# JpsisToMuMu

JpsisVariables = cms.PSet(
    mass         = Var("mass",                         float, doc = "Unfit invariant mass"),
    doca         = Var("userFloat('doca')",            float, doc = "Distance of closest approach of tracks"),
    iso          = Var("userFloat('iso')",             float, doc = "tracks isolation (pt/(pt+sum))"),
    muon1_pt     = Var("userFloat('muon1_pt')",        float, doc = "Muon 1 pt"),
    muon1_eta    = Var("userFloat('muon1_eta')",       float, doc = "Muon 1 eta"),
    muon1_phi    = Var("userFloat('muon1_phi')",       float, doc = "Muon 1 phi"),
#    trk1_mu_index = Var("userInt('trk1_mu_index')",      int, doc = "Matched muon index for track 1"),
    muon2_pt     = Var("userFloat('muon2_pt')",        float, doc = "Muon 2 pt"),
    muon2_eta    = Var("userFloat('muon2_eta')",       float, doc = "Muon 2 eta"),
    muon2_phi    = Var("userFloat('muon2_phi')",       float, doc = "Muon 2 phi"),
#    trk2_mu_index = Var("userInt('trk2_mu_index')",      int, doc = "Matched muon index for track 2"),
#    trk1_sip     = Var("userFloat('trk1_sip')",        float, doc = "Track 1 2D impact parameter significance wrt Beam Spot"),
#    trk2_sip     = Var("userFloat('trk2_sip')",        float, doc = "Track 2 2D impact parameter significance wrt Beam Spot"),
    kin_valid    = Var("userInt('kin_valid')",         int,   doc = "Kinematic fit: vertex validity"),
    kin_vtx_prob = Var("userFloat('kin_vtx_prob')",    float, doc = "Kinematic fit: vertex probability"),
    kin_vtx_chi2dof = Var("userFloat('kin_vtx_chi2dof')", float, doc = "Kinematic fit: vertex normalized Chi^2"),
    kin_mass     = Var("userFloat('kin_mass')",        float, doc = "Kinematic fit: vertex refitted mass"),
    kin_pt       = Var("userFloat('kin_pt')",          float, doc = "Kinematic fit: vertex refitted pt"),
    kin_eta      = Var("userFloat('kin_eta')",         float, doc = "Kinematic fit: vertex refitted eta"),
    kin_phi      = Var("userFloat('kin_phi')",         float, doc = "Kinematic fit: vertex refitted phi"),
    kin_massErr  = Var("userFloat('kin_massErr')",     float, doc = "Kinematic fit: vertex refitted mass error"),
    kin_lxy      = Var("userFloat('kin_lxy')",         float, doc = "Kinematic fit: vertex displacement in XY plane wrt Beam Spot"),
    kin_slxy     = Var("userFloat('kin_sigLxy')",      float, doc = "Kinematic fit: vertex displacement significance in XY plane wrt Beam Spot"),
    kin_cosAlphaXY = Var("userFloat('kin_cosAlphaXY')",    float, doc = "Kinematic fit: cosine of pointing angle in XY wrt BS"),
    kin_sipBS    = Var("userFloat('kin_sipBS')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in XY wrt BS"),
    kin_sipPV    = Var("userFloat('kin_sipPV')",       float, doc = "Kinematic fit: impact parameter significance of the candidate trajectory in 3D wrt PV"),
)

JpsisVariablesMC = merge_psets(
    JpsisVariables,
    cms.PSet(
        gen_trk1_pdgId  = Var("userInt(  'gen_muon1_pdgId')",    int,   doc = "Gen match: first track pdg Id"),
        gen_trk1_mpdgId = Var("userInt(  'gen_muon1_mpdgId')",   int,   doc = "Gen match: first track mother pdg Id"),
        gen_trk1_pt     = Var("userFloat('gen_muon1_pt')",     float,   doc = "Gen match: first track pt"),
        gen_trk2_pdgId  = Var("userInt(  'gen_muon2_pdgId')",    int,   doc = "Gen match: second track pdg Id"),
        gen_trk2_mpdgId = Var("userInt(  'gen_muon2_mpdgId')",   int,   doc = "Gen match: second track mother pdg Id"),
        gen_trk2_pt     = Var("userFloat('gen_muon2_pt')",     float,   doc = "Gen match: second track pt"),
        gen_pdgId       = Var("userInt(  'gen_pdgId')",         int,   doc = "Gen match: ditrack pdg Id"),
        gen_mass        = Var("userFloat('gen_mass')",        float,   doc = "Gen match: ditrack mass"),
        gen_pt          = Var("userFloat('gen_pt')",          float,   doc = "Gen match: ditrack pt"),
        ),
)


JpsisTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DiMuProd","Jpsi"),
    cut=cms.string(""),
    name=cms.string("Jpsi"),
    doc=cms.string("Jpsi Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = JpsisVariables
)

JpsisMcTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DiMuProdMC","Jpsi"),
    cut=cms.string(""),
    name=cms.string("Jpsi"),
    doc=cms.string("Jpsi Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = JpsisVariablesMC
)

UpsilonsTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DiMuProd","Upsilon"),
    cut=cms.string(""),
    name=cms.string("Upsilon"),
    doc=cms.string("Upsilons Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = JpsisVariables
)

UpsilonsMcTable=cms.EDProducer("SimpleCompositeCandidateFlatTableProducer", 
    src=cms.InputTag("DiMuProdMC","Upsilon"),
    cut=cms.string(""),
    name=cms.string("Upsilon"),
    doc=cms.string("Upsilons Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables = JpsisVariablesMC
)



DiMuProdSequence   = cms.Sequence(DiMuProd)
DiMuProdMcSequence = cms.Sequence(DiMuProdMC)
DiMuTables     = cms.Sequence(JpsisTable + UpsilonsTable)
DiMuMcTables   = cms.Sequence(JpsisMcTable + UpsilonsMcTable)
