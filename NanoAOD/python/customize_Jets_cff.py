import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_cff import *

def myRun2Jets(process,runOnMC):

    #
    # Check DiMuProdMC
    #
    from PhysicsTools.PatAlgos.tools.helpers  import getPatAlgosToolsTask, addToProcessAndTask
    task = getPatAlgosToolsTask(process)

    #
    # Select CHS particleFlow and create packedCandidates
    #

    NumOfPUVtxsForCharged=2
    DzCutForChargedFromPUVtxs=0.2

    process.mypfCHS = cms.EDFilter("CandPtrSelector",
                         src = cms.InputTag("DiMuProdMC","PFCandNoMeson"),
#                         src = cms.InputTag("packedPFCandidates"),
                         cut = cms.string("fromPV(0)>0 || (vertexRef().key<={} && abs(dz(0))<{})".format(NumOfPUVtxsForCharged,DzCutForChargedFromPUVtxs))
                         )

    #
    # Recluster jets
    #
    process.load("RecoJets.JetProducers.ak4PFJets_cfi")
    process.ak4PFJetsCHS.src = cms.InputTag("mypfCHS")

    process.jetTask = cms.Task(process.mypfCHS,process.ak4PFJetsCHS)
    task.add(process.jetTask)

    pvLabel = "offlineSlimmedPrimaryVertices"
    svLabel = "slimmedSecondaryVertices"
    muLabel = "slimmedMuons"
    elLabel = "slimmedElectrons"
    gpLabel = "prunedGenParticles"
    genJetsCollection = "slimmedGenJets"

    ## and add them to the event content
    from RecoBTag.ONNXRuntime.pfParticleNetAK4_cff import _pfParticleNetAK4JetTagsAll as pfParticleNetAK4JetTagsAll
    from RecoBTag.ONNXRuntime.pfDeepFlavourJetTags_cfi import pfDeepFlavourJetTags

    btagDiscriminatorsDeepFlavor = [
            'pfDeepFlavourJetTags:probb',
            'pfDeepFlavourJetTags:probbb',
            'pfDeepFlavourJetTags:problepb',
            'pfDeepFlavourJetTags:probc',
            'pfDeepFlavourJetTags:probuds',
            'pfDeepFlavourJetTags:probg',
    ]

    bTagCSVV2    = ['pfCombinedInclusiveSecondaryVertexV2BJetTags']
    bTagDeepCSV  = ['pfDeepCSVJetTags:probb','pfDeepCSVJetTags:probbb','pfDeepCSVJetTags:probc','pfDeepCSVJetTags:probudsg']

    #
    # PATify jets
    #
    tmppostfix = "NoMuFromJpsi"


#    this create patJetsCHSNoMuFromJpsi,selectedPatJetsCHSNoMuFromJpsi
    from PhysicsTools.PatAlgos.tools.jetTools import supportedJetAlgos, addJetCollection, updateJetCollection
    addJetCollection(
        process,
        postfix            = tmppostfix,
        labelName          = "CHS",
        jetSource          = cms.InputTag("ak4PFJetsCHS"),
        algo               = "ak",
        rParam             = 0.4,
        pfCandidates       = cms.InputTag("mypfCHS"),
        pvSource           = cms.InputTag(pvLabel),
        svSource           = cms.InputTag(svLabel),
        muSource           = cms.InputTag(muLabel),
        elSource           = cms.InputTag(elLabel),
        genJetCollection   = cms.InputTag(genJetsCollection),
        genParticles       = cms.InputTag(gpLabel),
        jetCorrections     = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
        getJetMCFlavour    = runOnMC
    )


    addToProcessAndTask('selectedPatJetsCHS'+tmppostfix+'NoFlavour', process.selectedPatJetsCHSNoMuFromJpsi.clone(), process, task)
    

#    this create process.patJetsCHSNoMuFromJpsi
    updateJetCollection(
        process,
        jetSource = cms.InputTag('selectedPatJetsCHS'+tmppostfix+'NoFlavour'),
        pvSource = cms.InputTag(pvLabel),
        svSource = cms.InputTag(svLabel),
        muSource = cms.InputTag(muLabel),
        elSource = cms.InputTag(elLabel),
        jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
        btagDiscriminators = btagDiscriminatorsDeepFlavor + pfParticleNetAK4JetTagsAll + bTagCSVV2 + bTagDeepCSV,
        postfix='CHS'+tmppostfix+'WithDeepInfo',
    )

    # slimmedJets with DeepFlavour (remove DeepFlavour-less)
    delattr(process, 'selectedPatJetsCHS'+tmppostfix)
    addToProcessAndTask('selectedPatJetsCHS'+tmppostfix, getattr(process,'selectedUpdatedPatJetsCHS'+tmppostfix+'WithDeepInfo').clone(), process, task)
    # delete module not used anymore (slimmedJets substitutes)
    delattr(process, 'selectedUpdatedPatJetsCHS'+tmppostfix+'WithDeepInfo')

    #TagInfos
    process.pfParticleNetAK4TagInfosCHSNoMuFromJpsiWithDeepInfo.pf_candidates = cms.InputTag("DiMuProdMC","PFCandNoMeson")
    process.pfDeepFlavourTagInfosCHSNoMuFromJpsiWithDeepInfo.candidates = cms.InputTag("DiMuProdMC","PFCandNoMeson")
    process.pfDeepCSVTagInfosCHSNoMuFromJpsiWithDeepInfo.pf_candidates = cms.InputTag("DiMuProdMC","PFCandNoMeson")
    process.pfImpactParameterTagInfosCHSNoMuFromJpsiWithDeepInfo.candidates = cms.InputTag("DiMuProdMC","PFCandNoMeson")

    ### $$$$
    ### $$$$
    # below is to add the userfloats puIdDisc,puId,jetId,qgl
    # https://github.com/cms-sw/cmssw/blob/b8776f576fdc5aa8e36d0396fd3aa567ea05877f/PhysicsTools/NanoAOD/python/jetsAK4_CHS_cff.py#L39

    # Run2
    from PhysicsTools.NanoAOD.jets_cff import jetTable, updatedJetsWithUserData, jercVars
    # Run3
    #    from PhysicsTools.NanoAOD.jetsAK4_CHS_cff import jetTable, updatedJetsWithUserData, jercVars

    process.updatedJetsWithUserDataNoMuFromJpsi = updatedJetsWithUserData.clone(
        src = cms.InputTag('selectedPatJetsCHS'+tmppostfix),
        userFloats = cms.PSet(),
        userInts = cms.PSet()
    )

    task.add(process.updatedJetsWithUserDataNoMuFromJpsi)

    ##############################################################
    # Build jet and MC table
    ###############################################################

    process.jetFinalTable = jetTable.clone(
        src = cms.InputTag("updatedJetsWithUserData"+ tmppostfix),
        doc = cms.string('ak4 PFJets CHS with JECs applied, after basic selection (pt > 15), and after the removal of the two muons'),
        name = cms.string("Jet"+ tmppostfix),
        externalVariables = cms.PSet(),
        variables = cms.PSet(P4Vars,
                             area = Var("jetArea()", float, doc="jet catchment area, for JECs",precision=10),
                             nMuons = Var("?hasOverlaps('muons')?overlaps('muons').size():0", int, doc="number of muons in the jet"),
                             muonIdx1 = Var("?overlaps('muons').size()>0?overlaps('muons')[0].key():-1", int, doc="index of first matching muon"),
                             muonIdx2 = Var("?overlaps('muons').size()>1?overlaps('muons')[1].key():-1", int, doc="index of second matching muon"),
                             electronIdx1 = Var("?overlaps('electrons').size()>0?overlaps('electrons')[0].key():-1", int, doc="index of first matching electron"),
                             electronIdx2 = Var("?overlaps('electrons').size()>1?overlaps('electrons')[1].key():-1", int, doc="index of second matching electron"),
                             nElectrons = Var("?hasOverlaps('electrons')?overlaps('electrons').size():0", int, doc="number of electrons in the jet"),
                             ####
                             btagPNetB = Var("?bDiscriminator('pfParticleNetAK4DiscriminatorsJetTags:BvsAll')>0?bDiscriminator('pfParticleNetAK4DiscriminatorsJetTags:BvsAll'):-1",float,precision=10,doc="ParticleNet b vs. udscg"),
                             btagPNetCvL = Var("?bDiscriminator('pfParticleNetAK4DiscriminatorsJetTags:CvsL')>0?bDiscriminator('pfParticleNetAK4DiscriminatorsJetTags:CvsL'):-1",float,precision=10,doc="ParticleNet c vs. udsg"),
                             btagPNetCvB = Var("?bDiscriminator('pfParticleNetAK4DiscriminatorsJetTags:CvsB')>0?bDiscriminator('pfParticleNetAK4DiscriminatorsJetTags:CvsB'):-1",float,precision=10,doc="ParticleNet c vs. b"),
                             btagPNetQvsG = Var("?bDiscriminator('pfParticleNetAK4DiscriminatorsJetTags:QvsG')>0?bDiscriminator('pfParticleNetAK4DiscriminatorsJetTags:QvsG'):-1",float,precision=10,doc="ParticleNet uds vs. g"),
                             #### those below are in the standard nano V9 UL
                             btagDeepB = Var("?(bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb'))>=0?bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb'):-1",float,doc="DeepCSV b+bb tag discriminator",precision=10),
                             btagDeepFlavB = Var("bDiscriminator('pfDeepFlavourJetTags:probb')+bDiscriminator('pfDeepFlavourJetTags:probbb')+bDiscriminator('pfDeepFlavourJetTags:problepb')",float,doc="DeepJet b+bb+lepb tag discriminator",precision=10),
                             btagCSVV2 = Var("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')",float,doc=" pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)",precision=10),
                             btagDeepCvL = Var("?bDiscriminator('pfDeepCSVJetTags:probc')>=0?bDiscriminator('pfDeepCSVJetTags:probc')/(bDiscriminator('pfDeepCSVJetTags:probc')+bDiscriminator('pfDeepCSVJetTags:probudsg')):-1", float,doc="DeepCSV c vs udsg discriminator",precision=10),
                             btagDeepCvB = Var("?bDiscriminator('pfDeepCSVJetTags:probc')>=0?bDiscriminator('pfDeepCSVJetTags:probc')/(bDiscriminator('pfDeepCSVJetTags:probc')+bDiscriminator('pfDeepCSVJetTags:probb')+bDiscriminator('pfDeepCSVJetTags:probbb')):-1",float,doc="DeepCSV c vs b+bb discriminator",precision=10),
                             btagDeepFlavCvL = Var("?(bDiscriminator('pfDeepFlavourJetTags:probc')+bDiscriminator('pfDeepFlavourJetTags:probuds')+bDiscriminator('pfDeepFlavourJetTags:probg'))>0?bDiscriminator('pfDeepFlavourJetTags:probc')/(bDiscriminator('pfDeepFlavourJetTags:probc')+bDiscriminator('pfDeepFlavourJetTags:probuds')+bDiscriminator('pfDeepFlavourJetTags:probg')):-1",float,doc="DeepJet c vs uds+g discriminator",precision=10),
                             btagDeepFlavCvB = Var("?(bDiscriminator('pfDeepFlavourJetTags:probc')+bDiscriminator('pfDeepFlavourJetTags:probb')+bDiscriminator('pfDeepFlavourJetTags:probbb')+bDiscriminator('pfDeepFlavourJetTags:problepb'))>0?bDiscriminator('pfDeepFlavourJetTags:probc')/(bDiscriminator('pfDeepFlavourJetTags:probc')+bDiscriminator('pfDeepFlavourJetTags:probb')+bDiscriminator('pfDeepFlavourJetTags:probbb')+bDiscriminator('pfDeepFlavourJetTags:problepb')):-1",float,doc="DeepJet c vs b+bb+lepb discriminator",precision=10),
                             btagDeepFlavQG = Var("?(bDiscriminator('pfDeepFlavourJetTags:probg')+bDiscriminator('pfDeepFlavourJetTags:probuds'))>0?bDiscriminator('pfDeepFlavourJetTags:probg')/(bDiscriminator('pfDeepFlavourJetTags:probg')+bDiscriminator('pfDeepFlavourJetTags:probuds')):-1",float,doc="DeepJet g vs uds discriminator",precision=10),
                             ####
#                             puIdDisc = Var("userFloat('puId106XUL18Disc')", float,doc="Pileup ID discriminant with 106X (2018) training",precision=10),
#                             puId = Var("userInt('puId106XUL18Id')", int,doc="Pileup ID flags with 106X (2018) training"),
#                             jetId = Var("userInt('tightId')*2+4*userInt('tightIdLepVeto')",int,doc="Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto"),
#                             qgl = Var("?userFloat('qgl')>0?userFloat('qgl'):-1",float,doc="Quark vs Gluon likelihood discriminator",precision=10),
                             nConstituents = Var("numberOfDaughters()","uint8",doc="Number of particles in the jet"),
                             rawFactor = Var("1.-jecFactor('Uncorrected')",float,doc="1 - Factor to get back to raw pT",precision=6),
                             chHEF = Var("chargedHadronEnergyFraction()", float, doc="charged Hadron Energy Fraction", precision= 6),
                             neHEF = Var("neutralHadronEnergyFraction()", float, doc="neutral Hadron Energy Fraction", precision= 6),
                             chEmEF = Var("chargedEmEnergyFraction()", float, doc="charged Electromagnetic Energy Fraction", precision= 6),
                             neEmEF = Var("neutralEmEnergyFraction()", float, doc="neutral Electromagnetic Energy Fraction", precision= 6),
                             muEF = Var("muonEnergyFraction()", float, doc="muon Energy Fraction", precision= 6),
#                             chFPV0EF = Var("userFloat('chFPV0EF')", float, doc="charged fromPV==0 Energy Fraction (energy excluded from CHS jets). Previously called betastar.", precision= 6),
                             )
    )

    from PhysicsTools.NanoAOD.jets_cff import jetMCTable
#    from PhysicsTools.NanoAOD.jetMC_cff import jetMCTable # in late releases not in the 10_6_27
    
    process.jetFinalMCTable = jetMCTable.clone(
        src = cms.InputTag("updatedJetsWithUserData"+ tmppostfix),
        name = cms.string("Jet"+ tmppostfix),
        extension = cms.bool(True), # this is an extension  table for the jets
        variables = cms.PSet(
            partonFlavour = Var("partonFlavour()", int, doc="flavour from parton matching"),
            hadronFlavour = Var("hadronFlavour()", int, doc="flavour from hadron ghost clustering"),
            genJetIdx = Var("?genJetFwdRef().backRef().isNonnull()?genJetFwdRef().backRef().key():-1", int, doc="index of matched gen jet"),
        )
    )

    process.pfTestMC = cms.EDAnalyzer("PFCandsAnalyzer",
                                      inputCandidates = cms.InputTag("packedPFCandidates"),
                                      slimmedCandidates = cms.InputTag("DiMuProdMC","PFCandNoMeson"),
                                      orgJet = cms.InputTag("slimmedJets"),
                                      modJet = cms.InputTag("patJetsCHS"+tmppostfix),
                                      )

    process.jetTablesTask = cms.Task(process.jetFinalTable,process.jetFinalMCTable)
    task.add(process.jetTablesTask)

    ##############################################################
    # Take AK4 jets and collect their PF constituents
    ###############################################################

    process.finalJetsForConst = cms.EDProducer("PatJetConstituentPtrSelector",
                                               src = cms.InputTag("updatedJetsWithUserData"+ tmppostfix),
                                               cut = cms.string("abs(eta) <= 2.5 && pt>20")
                                               )

    ##############################################################
    # Setup PF candidates table
    ##############################################################

    process.finalJetsConstituents = cms.EDProducer("PackedCandidatePtrMerger",
                                                   src = cms.VInputTag(cms.InputTag("finalJetsForConst", "constituents")), ### candList,
                                                   skipNulls = cms.bool(True),
                                                   warnOnSkip = cms.bool(True)
                                                   )
    process.customConstituentsExtTable = cms.EDProducer("SimplePATCandidateFlatTableProducer",
                                                        src = cms.InputTag("finalJetsConstituents"), ### candInput,
                                                        cut = cms.string(""),
                                                        name = cms.string("PFCand"),
                                                        doc = cms.string("PF candidate constituents of AK4 jets with |eta| <= 2.5 and pt>20"),
                                                        singleton = cms.bool(False),
                                                        extension = cms.bool(False),
                                                        variables = cms.PSet(
                                                            pt = Var("pt", float, doc="pt", precision=10),
                                                            mass = Var("mass", float, doc="mass", precision=10),
                                                            eta = Var("eta", float, precision=12),
                                                            phi = Var("phi", float, precision=12),
                                                            pdgId  = Var("pdgId", int, doc="PF candidate type (+/-211 = ChgHad, 130 = NeuHad, 22 = Photon, +/-11 = Electron, +/-13 = Muon, 1 = HFHad, 2 = HFEM)"),
                                                            dz = Var("?hasTrackDetails()?dz():-1", float, doc="pf dz", precision=10),
                                                            dzErr = Var("?hasTrackDetails()?dzError():-1", float, doc="pf dz err", precision=10),
                                                            dxy = Var("?hasTrackDetails()?dxy():-1", float, doc="pf d0", precision=10),
                                                            dxyErr = Var("?hasTrackDetails()?dxyError():-1", float, doc="pf d0 err", precision=10),
                                                        )
                                                        )

    process.customAK4ConstituentsTable = cms.EDProducer("Run2PatJetConstituentTableProducer",
                                                        candidates = cms.InputTag("finalJetsConstituents"), ### candInput,
                                                        jets = cms.InputTag("finalJetsForConst"),
                                                        jet_radius = cms.double(0.4),
                                                        name = cms.string("JetPFCands"),
#                                                        idx_name = cms.string("pFCandsIdx"),
                                                        nameSV = cms.string("JetSVs"),
#                                                        idx_nameSV = cms.string("sVIdx"),
                                                        )

    process.jetConstituentsTask = cms.Task(process.finalJetsForConst,process.finalJetsConstituents)
    process.jetConstituentsTablesTask = cms.Task(process.finalJetsConstituents,process.customConstituentsExtTable,process.customAK4ConstituentsTable)
    task.add(process.jetConstituentsTask)
    task.add(process.jetConstituentsTablesTask)

#    process.myJetSequence= cms.Sequence(cms.Sequence(task)*process.pfTestMC)
    process.myJetSequence= cms.Sequence(cms.Sequence(task))    
    
    return process

