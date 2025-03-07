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
    iso          = Var("userFloat('iso')",             float, doc = "tracks isolation (pt/(pt+sum)) around dimuon of R=03"),
    isoPho       = Var("userFloat('isoPho')",          float, doc = "photon isolation (pt/(pt+sum)) around dimuon of R=03"),
    isoNeuHad    = Var("userFloat('isoNeuHad')",       float, doc = "neutral Hadron isolation (pt/(pt+sum)) around dimuon of R=03"),
    phoPtSumOutsideSignalConedR03       = Var("userFloat('phoPtSumOutsideSignalConedR03')",       float, doc = "phoPtSumOutsideSignalConedR03"),
    ##
    leadingChargedMaxPt1 = Var("userFloat('leadingChargedMaxPt1')", float, doc = "leadingCharged Pt in a signal cone around dimuon of R=04"),
    leadingChargedDxyMaxPt1 = Var("userFloat('leadingChargedDxyMaxPt1')", float, doc = "leadingCharged Dxy in a signal cone around dimuon of R=04"),
    leadingChargedDzMaxPt1 = Var("userFloat('leadingChargedDzMaxPt1')", float, doc = "leadingCharged Dxy in a signal cone around dimuon of R=04"),
    leadingCharged2dSignMaxPt1 = Var("userFloat('leadingCharged2dSignMaxPt1')", float, doc = "leadingCharged 2-sign in a signal cone around dimuon of R=04"),
    leadingCharged3dSignMaxPt1 = Var("userFloat('leadingCharged3dSignMaxPt1')", float, doc = "leadingCharged 3-sign in a signal cone around dimuon of R=04"),
    leadingChargedEtaMaxPt1 = Var("userFloat('leadingChargedEtaMaxPt1')", float, doc = "leadingCharged eta in a signal cone around dimuon of R=04"),
    leadingChargedPhiMaxPt1 = Var("userFloat('leadingChargedPhiMaxPt1')", float, doc = "leadingCharged phi in a signal cone around dimuon of R=04"),
    ##
    leadingChargedMaxPt2 = Var("userFloat('leadingChargedMaxPt2')", float, doc = "SubLeadingCharged Pt in a signal cone around dimuon of R=04"),
    leadingChargedDxyMaxPt2 = Var("userFloat('leadingChargedDxyMaxPt2')", float, doc = "SubLeadingCharged Dxy in a signal cone around dimuon of R=04"),
    leadingChargedDzMaxPt2 = Var("userFloat('leadingChargedDzMaxPt2')", float, doc = "SubLeadingCharged Dxy in a signal cone around dimuon of R=04"),
    leadingCharged2dSignMaxPt2 = Var("userFloat('leadingCharged2dSignMaxPt2')", float, doc = "SubLeadingCharged 2-sign in a signal cone around dimuon of R=04"),
    leadingCharged3dSignMaxPt2 = Var("userFloat('leadingCharged3dSignMaxPt2')", float, doc = "SubLeadingCharged 3-sign in a signal cone around dimuon of R=04"),
    leadingChargedEtaMaxPt2 = Var("userFloat('leadingChargedEtaMaxPt2')", float, doc = "SubLeadingCharged eta in a signal cone around dimuon of R=04"),
    leadingChargedPhiMaxPt2 = Var("userFloat('leadingChargedPhiMaxPt2')", float, doc = "SubLeadingCharged phi in a signal cone around dimuon of R=04"),
    ##
    muon1_pt       = Var("userFloat('muon1_pt')",       float, doc = "Muon 1 pt"),
    muon1_eta      = Var("userFloat('muon1_eta')",      float, doc = "Muon 1 eta"),
    muon1_phi      = Var("userFloat('muon1_phi')",      float, doc = "Muon 1 phi"),
    muon1_charge   = Var("userInt('muon1_charge')",     int, doc = "Muon 1 charge"),
    muon1_dz       = Var("userFloat('muon1_dz')",       float, doc = "Muon 1 dz"),
    muon1_dxy      = Var("userFloat('muon1_dxy')",      float, doc = "Muon 1 dxy"),
    muon1_Sip3dSig = Var("userFloat('muon1_Sip3dSig')", float, doc = "Muon 1 Sip3dSig"),
    muon1_Sip3dVal = Var("userFloat('muon1_Sip3dVal')", float, doc = "Muon 1 Sip3dVal"),
    muon1_Sip2dSig = Var("userFloat('muon1_Sip2dSig')", float, doc = "Muon 1 Sip2dSig"),
    muon1_Sip2dVal = Var("userFloat('muon1_Sip2dVal')", float, doc = "Muon 1 Sip2dVal"),
    muon1_DistVal  = Var("userFloat('muon1_DistVal')",  float, doc = "Muon 1 DistVal"),
    #    trk1_mu_index = Var("userInt('trk1_mu_index')",      int, doc = "Matched muon index for track 1"),
    muon2_pt       = Var("userFloat('muon2_pt')",       float, doc = "Muon 2 pt"),
    muon2_eta      = Var("userFloat('muon2_eta')",      float, doc = "Muon 2 eta"),
    muon2_phi      = Var("userFloat('muon2_phi')",      float, doc = "Muon 2 phi"),
    muon2_charge   = Var("userInt('muon2_charge')",     int, doc = "Muon 2 charge"),
    muon2_dz       = Var("userFloat('muon2_dz')",       float, doc = "Muon 2 dz"),
    muon2_dxy      = Var("userFloat('muon2_dxy')",      float, doc = "Muon 2 dxy"),
    muon2_Sip3dSig = Var("userFloat('muon2_Sip3dSig')", float, doc = "Muon 2 Sip3dSig"),
    muon2_Sip3dVal = Var("userFloat('muon2_Sip3dVal')", float, doc = "Muon 2 Sip3dVal"),
    muon2_Sip2dSig = Var("userFloat('muon2_Sip2dSig')", float, doc = "Muon 2 Sip2dSig"),
    muon2_Sip2dVal = Var("userFloat('muon2_Sip2dVal')", float, doc = "Muon 2 Sip2dVal"),
    muon2_DistVal  = Var("userFloat('muon2_DistVal')",  float, doc = "Muon 2 DistVal"),
    #    trk2_mu_index = Var("userInt('trk2_mu_index')",      int, doc = "Matched muon index for track 2"),
    muon1_isTighMuon     = Var("userInt('muon1_isTightMuon')",        int, doc = "isTightMuon"),
    muon2_isTighMuon     = Var("userInt('muon2_isTightMuon')",        int, doc = "isTightMuon"),
    muon1_isMediumMuon   = Var("userInt('muon1_isMediumMuon')",      int, doc = "isMediumMuon"),
    muon2_isMediumMuon   = Var("userInt('muon2_isMediumMuon')",      int, doc = "isMediumMuon"),
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
    kin_l3d      = Var("userFloat('kin_l3d')",         float, doc = "Kinematic fit: decay lenght in 3D wrt PV"),
    kin_sl3d     = Var("userFloat('kin_sl3d')",        float, doc = "Kinematic fit: decay lenght significance in 3D wrt PV"),
#
    kin_bestVtx_x  = Var("userFloat('kin_bestVtx_x')",       float, doc = "best vtx SV x"),
    kin_bestVtx_y  = Var("userFloat('kin_bestVtx_y')",       float, doc = "best vtx SV y"),
    kin_bestVtx_z  = Var("userFloat('kin_bestVtx_z')",       float, doc = "best vtx SV z"),
    kin_pvIndex    = Var("userInt('kin_pvIndex')",           int,   doc = "Kinematic fit: index of the PV"),

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
