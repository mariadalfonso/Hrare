#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"

#include "DataFormats/Candidate/interface/Candidate.h"
typedef SimpleFlatTableProducer<reco::Candidate> SimpleCandidateFlatTableProducer;

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
typedef EventSingletonSimpleFlatTableProducer<GenEventInfoProduct> SimpleGenEventFlatTableProducer;

#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
typedef EventSingletonSimpleFlatTableProducer<HTXS::HiggsClassification> SimpleHTXSFlatTableProducer;

// already reo-->pat in Run3
//#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
//typedef SimpleFlatTableProducer<pat::CompositeCandidate> SimpleCompositeCandidateFlatTableProducer;

// this we need only in the 10_6_27 is already there for Run3
//#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
//typedef SimpleFlatTableProducer<pat::PackedCandidate> SimplePATCandidateFlatTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(SimpleCompositeCandidateFlatTableProducer);
//DEFINE_FWK_MODULE(SimplePATCandidateFlatTableProducer);
