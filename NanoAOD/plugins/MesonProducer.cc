#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/PointingKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TMVA/Reader.h"

#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include <TLorentzVector.h>
#include <TVector.h>
#include <TMatrix.h>
#include "TMatrixDSym.h"
#include <algorithm>

// 
// MesonProducer is designed for Bs/d->mumu analysis
//

typedef reco::Candidate::LorentzVector LorentzVector;

namespace {
  const float muon_mass_     = 0.10565837;
  const float kaon_mass_     = 0.493677;
  const float mass_err_      = 1.6e-5;
  const float pion_mass_     = 0.139570;
  const float pion0_mass_    = 0.134977;
  const float jpsi_mass_     = 3.0969;
  const float proton_mass_   = 0.9382;
};

struct KinematicFitResult{
  bool treeIsValid;
  bool vertexIsValid;
  RefCountedKinematicVertex      refitVertex;
  RefCountedKinematicParticle    refitMother;
  RefCountedKinematicTree        refitTree;
  std::vector<RefCountedKinematicParticle> refitDaughters;
  float lxy, lxyErr, sigLxy, cosAlpha, ipSigBS, ipSigPV;
  KinematicFitResult():treeIsValid(false),vertexIsValid(false),
		       lxy(-1.0), lxyErr(-1.0), sigLxy(-1.0), cosAlpha(-999.),
		       ipSigBS(999.), ipSigPV(999.)
  {}

  bool valid() const {
    return treeIsValid and vertexIsValid;
  }

  void postprocess(const reco::BeamSpot& beamSpot,
		   const reco::VertexCollection& vertices)
  {
    if ( not valid() ) return;
    // displacement information
    TVector v(2);
    v[0] = refitVertex->position().x()-beamSpot.position().x();
    v[1] = refitVertex->position().y()-beamSpot.position().y();

    TMatrix errVtx(2,2);
    errVtx(0,0) = refitVertex->error().cxx();
    errVtx(0,1) = refitVertex->error().matrix()(0,1);
    errVtx(1,0) = errVtx(0,1);
    errVtx(1,1) = refitVertex->error().cyy();

    TMatrix errBS(2,2);
    errBS(0,0) = beamSpot.covariance()(0,0);
    errBS(0,1) = beamSpot.covariance()(0,1);
    errBS(1,0) = beamSpot.covariance()(1,0);
    errBS(1,1) = beamSpot.covariance()(1,1);
    
    lxy = sqrt(v.Norm2Sqr());
    lxyErr = sqrt( v*(errVtx*v) + v*(errBS*v) ) / lxy;
    if (lxyErr > 0) sigLxy = lxy/lxyErr;
    
    // compute cosAlpha 2D wrt BeamSpot
    v[0] = refitVertex->position().x()-beamSpot.position().x();
    v[1] = refitVertex->position().y()-beamSpot.position().y();
    TVector w(2);
    w[0] = refitMother->currentState().globalMomentum().x();
    w[1] = refitMother->currentState().globalMomentum().y();
    cosAlpha = v*w/sqrt(v.Norm2Sqr()*w.Norm2Sqr());

    // comput impact parameter
    auto transientTrack = refitMother->refittedTransientTrack();
    ipSigBS = transientTrack.stateAtBeamLine().transverseImpactParameter().significance();

    ipSigPV = 999.;
    for (const auto& vertex: vertices){
      auto impactParameter3D = IPTools::absoluteImpactParameter3D(transientTrack, vertex);
      if (impactParameter3D.first)
	if (ipSigPV > impactParameter3D.second.significance())
	  ipSigPV = impactParameter3D.second.significance();
    }
  }
  
  float mass() const
  {
    if ( not valid() ) return -1.0;
    return refitMother->currentState().mass();
  }

  float refit_mass(unsigned int i, unsigned int j) const
  {
    if ( not valid() ) return -1.0;
    if (i >= refitDaughters.size()) return -2.0;
    if (j >= refitDaughters.size()) return -3.0;
    if (refitDaughters.at(i)->currentState().globalMomentum().mag2()<0) return -4.0;
    if (refitDaughters.at(j)->currentState().globalMomentum().mag2()<0) return -5.0;
    auto momentum = refitDaughters.at(i)->currentState().globalMomentum() + 
      refitDaughters.at(j)->currentState().globalMomentum();
    auto energy1 = sqrt(refitDaughters.at(i)->currentState().globalMomentum().mag2() + 
			pow(refitDaughters.at(i)->currentState().mass(),2));
    auto energy2 = sqrt(refitDaughters.at(j)->currentState().globalMomentum().mag2() + 
			pow(refitDaughters.at(j)->currentState().mass(),2));
    return sqrt(pow(energy1+energy2,2)-momentum.mag2());
  }

  GlobalVector p3() const
  {
    if ( not valid() ) return GlobalVector();
    return refitMother->currentState().globalMomentum();
  }

  GlobalVector dau_p3(unsigned int i) const
  {
    if ( not valid() or i>=refitDaughters.size() ) return GlobalVector();
    return refitDaughters.at(i)->currentState().globalMomentum();
  }

  float massErr() const
  {
    if ( not valid() ) return -1.0;
    return sqrt(refitMother->currentState().kinematicParametersError().matrix()(6,6));
  }

  float chi2() const
  {
    if ( not valid() ) return -1.0;
    return refitVertex->chiSquared();
  }

  float ndof() const
  {
    return refitVertex->degreesOfFreedom();
  }

  float vtxProb() const
  {
    if ( not valid() ) return -1.0;
    return TMath::Prob((double)refitVertex->chiSquared(), int(rint(refitVertex->degreesOfFreedom())));
  }
  
};

struct DisplacementInformationIn3D{
  double decayLength, decayLengthErr, decayLength2, decayLength2Err,
    distaceOfClosestApproach, distaceOfClosestApproachErr, distaceOfClosestApproachSig,
    distaceOfClosestApproach2, distaceOfClosestApproach2Err, distaceOfClosestApproach2Sig,
    longitudinalImpactParameter, longitudinalImpactParameterErr, longitudinalImpactParameterSig,
    longitudinalImpactParameter2, longitudinalImpactParameter2Err,longitudinalImpactParameter2Sig,
    decayTime, decayTimeError, decayTimeXY, decayTimeXYError,
    alpha, alphaErr, alphaXY, alphaXYErr;
  const reco::Vertex *pv,*pv2;
  int pvIndex,pv2Index;
  DisplacementInformationIn3D():decayLength(-1.0), decayLengthErr(0.), decayLength2(-1.0), decayLength2Err(0.),
				distaceOfClosestApproach(-1.0), distaceOfClosestApproachErr(0.0), distaceOfClosestApproachSig(0.0),
				distaceOfClosestApproach2(-1.0), distaceOfClosestApproach2Err(0.0), distaceOfClosestApproach2Sig(0.0),
				longitudinalImpactParameter(0.0), longitudinalImpactParameterErr(0.), longitudinalImpactParameterSig(0.),
				longitudinalImpactParameter2(0.0), longitudinalImpactParameter2Err(0.), longitudinalImpactParameter2Sig(0.),
    decayTime(-999.), decayTimeError(-999.),
    decayTimeXY(-999.), decayTimeXYError(-999.),
    alpha(-999.), alphaErr(-999.),
    alphaXY(-999.), alphaXYErr(-999.),
    pv(0), pv2(0),
    pvIndex(-1), pv2Index(-1)
  {};
};


struct GenMatchInfo{
  const pat::PackedGenParticle* mc_trk1;
  const pat::PackedGenParticle* mc_trk2;
  const reco::Candidate*   match;
  GenMatchInfo():mc_trk1(0), mc_trk2(0), match(0)
  {}
};

using namespace std;

///////////////////////////////////////////////////////////////////////////
///                             P L U G I N
///////////////////////////////////////////////////////////////////////////

class MesonProducer : public edm::stream::EDProducer<> {
    
public:
    
  explicit MesonProducer(const edm::ParameterSet &iConfig);
    
  ~MesonProducer() override {};
    
    
private:
    
  virtual void produce(edm::Event&, const edm::EventSetup&);

  //  bool isGoodMuon(const pat::Muon& muon);
  bool isGoodTrack(const pat::PackedCandidate& track);
  bool displacedTrack(const pat::PackedCandidate& track);
  //  bool isGoodMuonProbe(const pat::PackedCandidate& track);
  bool isGoodPion(const pat::PackedCandidate& track);
  bool isGoodPair(const pat::PackedCandidate& track1,
		  const pat::PackedCandidate& track2);
    
  float
  trackImpactParameterSignificance( const pat::PackedCandidate& track);
  
  pat::CompositeCandidate
  getKsToPiPi(const edm::Event& iEvent,
	      const pat::PackedCandidate& pfCand1,
	      const pat::PackedCandidate& pfCand2);
  pat::CompositeCandidate
  getRhosToPiPi(const edm::Event& iEvent,
		const pat::PackedCandidate& pfCand1,
		const pat::PackedCandidate& pfCand2);
  pat::CompositeCandidate
  getOmegasToPiPiPi0(const edm::Event& iEvent,
		     const pat::PackedCandidate& pfCand1,
		     const pat::PackedCandidate& pfCand2);
  pat::CompositeCandidate
  getK0StarToKPi(const edm::Event& iEvent,
		 const pat::PackedCandidate& kaonCand,
		 const pat::PackedCandidate& pion);
  pat::CompositeCandidate
  getD0ToKPi(const edm::Event& iEvent,
	     const pat::PackedCandidate& kaonCand,
	     const pat::PackedCandidate& pion);
  pat::CompositeCandidate
  getD0ToKPiPi0(const edm::Event& iEvent,
		const pat::PackedCandidate& kaonCand,
		const pat::PackedCandidate& pion);
  pat::CompositeCandidate
  getPhiToKK(const edm::Event& iEvent,
	     const pat::PackedCandidate& pfCand1,
	     const pat::PackedCandidate& pfCand2);
  pat::CompositeCandidate
  getLambdaToPPi(const edm::Event& iEvent,
		 const pat::PackedCandidate& pfCand1,
		 const pat::PackedCandidate& pfCand2);

  KinematicFitResult 
  vertexWithKinematicFitter(std::vector<const reco::Track*> trks,
			    std::vector<float> masses);

  //  pair<double,double> computeDCA(const pat::PackedCandidate &kaon,
  //   				 reco::BeamSpot beamSpot);
  GenMatchInfo getGenMatchInfo( const pat::PackedCandidate& track1,
				    const pat::PackedCandidate& track2);
  // Two track DOCA
  float 
  distanceOfClosestApproach( const reco::Track* track1,
			     const reco::Track* track2 );
  float 
  distanceOfClosestApproach( const pat::PackedGenParticle* track1,
			     const pat::PackedGenParticle* track2);
  // Track to vertex DOCA
  Measurement1D
  distanceOfClosestApproach( const reco::Track* track,
			     RefCountedKinematicVertex vertex);
  Measurement1D 
  distanceOfClosestApproach( const reco::Track* track,
			     const reco::Vertex& vertex);

  DisplacementInformationIn3D
  compute3dDisplacement(const KinematicFitResult& fit,
			const reco::VertexCollection& vertices,
			bool closestIn3D = true);

  float
  computeCandIsolation(	const pat::PackedCandidate& pfCand1,
			const pat::PackedCandidate& pfCand2,
			unsigned int primaryVertexIndex,
			float minPt=0.9, float dR=0.7, int particleType=0,
			bool returnSum = false,
			std::vector<const pat::PackedCandidate*> ignoreTracks =
			std::vector<const pat::PackedCandidate*>());

  std::pair<float, float>
  getAlpha(const GlobalPoint& vtx_position, const GlobalError& vtx_error,
	   const GlobalPoint& ip_position,  const GlobalError& ip_error,
	   const GlobalVector &momentum,
	   bool transverse = true);

  KinematicFitResult 
  fillInfo(pat::CompositeCandidate& ksCand,
	   const edm::Event& iEvent,
	   const pat::PackedCandidate & cand1,
	   const pat::PackedCandidate & cand2, 
	   string cand1_name = "trk1",
	   string cand2_name = "trk2");

  reco::Vertex refit_vertex(edm::Event& iEvent, const edm::EventSetup& iSetup, size_t ipv, bool beamSpotContrant);
  reco::Track fix_track(const reco::Track *tk, double delta=1e-8);
  reco::Track fix_track(const reco::TrackRef& tk);
  //  reco::Track fix_track(const reco::Track *tk, double delta);

  // ----------member data ---------------------------
    

  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBuilderToken_;
  const TransientTrackBuilder* theTTBuilder_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> idealMagneticFieldRecordToken_;
  const MagneticField* bField_;

  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  const reco::BeamSpot* beamSpot_;

  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const reco::VertexCollection* primaryVertices_;
  vector<reco::Vertex> possibleVtxs_ = {};

  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandToken_;
  const pat::PackedCandidateCollection* pfCands_;
  edm::EDGetTokenT<std::vector<pat::PackedGenParticle> >   packedGenToken_;
  const std::vector<pat::PackedGenParticle>* packedGenParticles_;

  const AnalyticalImpactPointExtrapolator* impactPointExtrapolator_;

  bool isMC_;
  bool doFitxTrackCov = false;

  //  double minMuonPt_;
  //  double maxMuonEta_;
  double minPionPt_;
  double maxPionEta_;

  double minKsPreselectMass_;
  double maxKsPreselectMass_;
  double minKsMass_;
  double maxKsMass_;
  double minPhiPreselectMass_;
  double maxPhiPreselectMass_;
  double minPhiMass_;
  double maxPhiMass_;
  double minRhosPreselectMass_;
  double maxRhosPreselectMass_;
  double minRhosMass_;
  double maxRhosMass_;
  double minOmegasPreselectMass_;
  double maxOmegasPreselectMass_;
  double minDsPreselectMass_;
  double maxDsPreselectMass_;
  double minDsMass_;
  double maxDsMass_;
  double minD0pi0PreselectMass_;
  double maxD0pi0PreselectMass_;
  double minD0PreselectMass_;
  double maxD0PreselectMass_;
  double minD0Mass_;
  double maxD0Mass_;
  double minD0StarMass_;
  double maxD0StarMass_;
  double minK0StarPreselectMass_;
  double maxK0StarPreselectMass_;
  double minK0StarMass_;
  double maxK0StarMass_;
  double minLambdaPreselectMass_;
  double maxLambdaPreselectMass_;
  double minLambdaMass_;
  double maxLambdaMass_;
  double maxTwoTrackDOCA_;
  double minDisplaceTrackSignificance_;
  double maxLxy_;
  double minSigLxy_;
  double minCosAlpha_;
  double minVtxProb_;
};

MesonProducer::MesonProducer(const edm::ParameterSet &iConfig):
theTTBuilderToken_(esConsumes<>(edm::ESInputTag("", "TransientTrackBuilder"))),
theTTBuilder_(nullptr),
idealMagneticFieldRecordToken_(esConsumes<>()),
bField_(nullptr),
beamSpotToken_( consumes<reco::BeamSpot> ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
beamSpot_(nullptr),
vertexToken_( consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) ),
primaryVertices_(nullptr),
muonToken_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
pfCandToken_( consumes<pat::PackedCandidateCollection> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
pfCands_(nullptr),
packedGenToken_( consumes<std::vector<pat::PackedGenParticle>> ( iConfig.getParameter<edm::InputTag>( "packedGenParticleCollection" ) ) ),
packedGenParticles_(nullptr),
impactPointExtrapolator_(0),
isMC_(             iConfig.getParameter<bool>( "isMC" ) ),
//minMuonPt_(         iConfig.getParameter<double>( "minMuonPt" ) ),
//maxMuonEta_(        iConfig.getParameter<double>( "maxMuonEta" ) ),
minPionPt_(         iConfig.getParameter<double>( "minPionPt" ) ),
maxPionEta_(        iConfig.getParameter<double>( "maxPionEta" ) ),
minKsPreselectMass_(     iConfig.getParameter<double>( "minKsPreselectMass" ) ),
maxKsPreselectMass_(     iConfig.getParameter<double>( "maxKsPreselectMass" ) ),
minKsMass_(     iConfig.getParameter<double>( "minKsMass" ) ),
maxKsMass_(     iConfig.getParameter<double>( "maxKsMass" ) ),
minPhiPreselectMass_(     iConfig.getParameter<double>( "minPhiPreselectMass" ) ),
maxPhiPreselectMass_(     iConfig.getParameter<double>( "maxPhiPreselectMass" ) ),
minPhiMass_(     iConfig.getParameter<double>( "minPhiMass" ) ),
maxPhiMass_(     iConfig.getParameter<double>( "maxPhiMass" ) ),
minRhosPreselectMass_(     iConfig.getParameter<double>( "minRhosPreselectMass" ) ),
maxRhosPreselectMass_(     iConfig.getParameter<double>( "maxRhosPreselectMass" ) ),
minRhosMass_(     iConfig.getParameter<double>( "minRhosMass" ) ),
maxRhosMass_(     iConfig.getParameter<double>( "maxRhosMass" ) ),
minOmegasPreselectMass_(     iConfig.getParameter<double>( "minOmegasPreselectMass" ) ),
maxOmegasPreselectMass_(     iConfig.getParameter<double>( "maxOmegasPreselectMass" ) ),
minDsPreselectMass_(     iConfig.getParameter<double>( "minDsPreselectMass" ) ),
maxDsPreselectMass_(     iConfig.getParameter<double>( "maxDsPreselectMass" ) ),
minDsMass_(     iConfig.getParameter<double>( "minDsMass" ) ),
maxDsMass_(     iConfig.getParameter<double>( "maxDsMass" ) ),
minD0pi0PreselectMass_(     iConfig.getParameter<double>( "minD0pi0PreselectMass" ) ),
maxD0pi0PreselectMass_(     iConfig.getParameter<double>( "maxD0pi0PreselectMass" ) ),
minD0PreselectMass_(     iConfig.getParameter<double>( "minD0PreselectMass" ) ),
maxD0PreselectMass_(     iConfig.getParameter<double>( "maxD0PreselectMass" ) ),
minD0Mass_(     iConfig.getParameter<double>( "minD0Mass" ) ),
maxD0Mass_(     iConfig.getParameter<double>( "maxD0Mass" ) ),
minD0StarMass_(     iConfig.getParameter<double>( "minD0StarMass" ) ),
maxD0StarMass_(     iConfig.getParameter<double>( "maxD0StarMass" ) ),
minK0StarPreselectMass_(     iConfig.getParameter<double>( "minK0StarPreselectMass" ) ),
maxK0StarPreselectMass_(     iConfig.getParameter<double>( "maxK0StarPreselectMass" ) ),
minK0StarMass_(     iConfig.getParameter<double>( "minK0StarMass" ) ),
maxK0StarMass_(     iConfig.getParameter<double>( "maxK0StarMass" ) ),
minLambdaPreselectMass_(     iConfig.getParameter<double>( "minLambdaPreselectMass" ) ),
maxLambdaPreselectMass_(     iConfig.getParameter<double>( "maxLambdaPreselectMass" ) ),
minLambdaMass_(     iConfig.getParameter<double>( "minLambdaMass" ) ),
maxLambdaMass_(     iConfig.getParameter<double>( "maxLambdaMass" ) ),
maxTwoTrackDOCA_( iConfig.getParameter<double>( "maxTwoTrackDOCA" ) ),
minDisplaceTrackSignificance_( iConfig.getParameter<double>( "minDisplaceTrackSignificance" ) ),
maxLxy_( iConfig.getParameter<double>( "maxLxy" ) ),
minSigLxy_( iConfig.getParameter<double>( "minSigLxy" ) ),
minCosAlpha_( iConfig.getParameter<double>( "minCosAlpha" ) ),
minVtxProb_( iConfig.getParameter<double>( "minVtxProb" ) )
{
    produces<pat::CompositeCandidateCollection>("Ks");
    produces<pat::CompositeCandidateCollection>("D0");
    produces<pat::CompositeCandidateCollection>("D0pi0");
    produces<pat::CompositeCandidateCollection>("K0Star");
    produces<pat::CompositeCandidateCollection>("Phi");
    produces<pat::CompositeCandidateCollection>("Lambda");
    produces<pat::CompositeCandidateCollection>("Rho");
    produces<pat::CompositeCandidateCollection>("Omega");
}

bool MesonProducer::isGoodTrack(const pat::PackedCandidate& track){
  if ( track.charge() == 0 ) return false;
  if ( not track.hasTrackDetails() ) return false;
  if ( not track.bestTrack()->quality(reco::Track::highPurity) ) return false; 
  return true;
}

float MesonProducer::trackImpactParameterSignificance(const pat::PackedCandidate& track){
  return track.bestTrack()->dxyError()>0 ? fabs(track.bestTrack()->dxy(*beamSpot_))/track.bestTrack()->dxyError():0.0;
}

bool MesonProducer::displacedTrack(const pat::PackedCandidate& track){
  return trackImpactParameterSignificance(track) > minDisplaceTrackSignificance_;
}

/*
bool MesonProducer::isGoodMuonProbe(const pat::PackedCandidate& track){
  return fabs(track.eta()) < maxMuonEta_ and track.pt()>minMuonPt_;
}
*/

bool MesonProducer::isGoodPion(const pat::PackedCandidate& track){
  return fabs(track.eta()) < maxPionEta_ and track.pt()>minPionPt_;
}

bool MesonProducer::isGoodPair(const pat::PackedCandidate& track1,
			    const pat::PackedCandidate& track2){

  /*
  // MARIA: this will not work for me 
  // try some small DR , opposite to the photon
  return (isGoodMuonProbe(track1) and isGoodPion(track2)) or 
    (isGoodMuonProbe(track2) and isGoodPion(track1));
  */
  return (isGoodPion(track2) and isGoodPion(track1));

}

namespace {
  void addFitInfo( pat::CompositeCandidate& cand, const KinematicFitResult& fit, std::string name ){
    cand.addUserInt(   name+"_valid",       fit.valid() );
    cand.addUserFloat( name+"_vtx_prob",    fit.vtxProb() );
    cand.addUserFloat( name+"_vtx_chi2dof", fit.chi2()>0?fit.chi2()/fit.ndof():-1);
    cand.addUserFloat( name+"_mass",        fit.mass() );
    cand.addUserFloat( name+"_massErr",     fit.massErr() );
    cand.addUserFloat( name+"_lxy",         fit.lxy );
    cand.addUserFloat( name+"_sigLxy",      fit.sigLxy );
    cand.addUserFloat( name+"_cosAlphaXY",  fit.cosAlpha );
    cand.addUserFloat( name+"_sipBS",       fit.ipSigBS );
    cand.addUserFloat( name+"_sipPV",       fit.ipSigPV );
    cand.addUserFloat( name+"_pt",          fit.p3().perp() );
    cand.addUserFloat( name+"_eta",         fit.p3().eta() );
    cand.addUserFloat( name+"_phi",         fit.p3().phi() );
  }    
}

reco::Vertex MesonProducer::refit_vertex(edm::Event& iEvent, const edm::EventSetup& iSetup, size_t ipv, bool beamSpotContrant)
{
  std::vector<reco::TransientTrack> mytracks;

  //  edm::ESHandle<TransientTrackBuilder> TTBuilder;
  //  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);

  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", recoBeamSpotHandle);
  reco::BeamSpot vertexBeamSpot = *recoBeamSpotHandle;

  // edm::Handle<vector<reco::Vertex>> vtxHandle;
  // iEvent.getByLabel("offlineSlimmedPrimaryVertices", vtxHandle);

  // cout << "========> Vtx: " << ipv << endl;
  for (const auto& ptk : *pfCands_){

    if (!ptk.hasTrackDetails()) continue;

    if (ptk.fromPV(ipv) < 3) continue;

    // Equivalent to the below
    // if (ptk.pvAssociationQuality() != pat::PackedCandidate::PVAssociationQuality::UsedInFitTight) continue;
    // auto vtxTK = ptk.vertexRef();
    // auto vtxMINI = (*vtxHandle)[ipv];
    // if (vtxTK->x() != vtxMINI.x() || vtxTK->y() != vtxMINI.y() || vtxTK->z() != vtxMINI.z()) continue;

    auto tk = ptk.bestTrack();
    reco::TransientTrack transientTrack = theTTBuilder_->build(fix_track(tk));
    transientTrack.setBeamSpot(vertexBeamSpot);
    mytracks.push_back(transientTrack);
  }

  if (mytracks.size() < 2) {
    // cout << "[WARNING]: Less than 2 tracks for vertex "<<ipv<<" fit" << endl;
    return reco::Vertex(); // Invalid vertex
  }

  AdaptiveVertexFitter theFitter;
  if (beamSpotContrant) {
    TransientVertex tmp = theFitter.vertex(mytracks, vertexBeamSpot);
    return tmp;
  }
  else {
    TransientVertex tmp = theFitter.vertex(mytracks);
    return tmp;
  }

}

reco::Track MesonProducer::fix_track(const reco::TrackRef& tk)
{
  reco::Track t = reco::Track(*tk);
  return fix_track(&t);
}

/* Check for a not positive definite covariance matrix. If the covariance
 * matrix is not positive definite, we force it to be positive definite by
 * adding the minimum eigenvalue to the diagonal of the covariance matrix plus
 * `delta`.
 *
 * See https://nhigham.com/2020/12/22/what-is-a-modified-cholesky-factorization/.
 *
 * Note: There may be better ways of doing this, but since most of these tracks
 * don't end up in the final sample, this is probably good enough. */
reco::Track MesonProducer::fix_track(const reco::Track *tk, double delta)
{
  unsigned int i, j;
  double min_eig = 1;

  /* Get the original covariance matrix. */
  reco::TrackBase::CovarianceMatrix cov = tk->covariance();

  /* Convert it from an SMatrix to a TMatrixD so we can get the eigenvalues. */
  TMatrixDSym new_cov(cov.kRows);
  for (i = 0; i < cov.kRows; i++) {
    for (j = 0; j < cov.kRows; j++) {
      /* Need to check for nan or inf, because for some reason these
       * cause a segfault when calling Eigenvectors().
       *
       * No idea what to do here or why this happens. */
      if (std::isnan(cov(i,j)) || std::isinf(cov(i,j)))
	cov(i,j) = 1e-6;
      new_cov(i,j) = cov(i,j);
    }
  }

  /* Get the eigenvalues. */
  TVectorD eig(cov.kRows);
  new_cov.EigenVectors(eig);
  for (i = 0; i < cov.kRows; i++)
    if (eig(i) < min_eig)
      min_eig = eig(i);

  /* If the minimum eigenvalue is less than zero, then subtract it from the
   * diagonal and add `delta`. */
  if (min_eig < 0) {
    for (i = 0; i < cov.kRows; i++)
      cov(i,i) -= min_eig - delta;
  }

  /* Correct for the fact that the impact parameter error is significantly
   * different in data and MC. Technically we should be *decreasing* the data
   * uncertainty since according to Marco Musich:
   *
   * > we suspect that the track uncertainty in data is over-estimated, while
   * > it should be well calibrated in MC.  To corroborate this, here is a
   * > plot from the recent alignment paper which shows the pull of the pixel
   * > residuals for various reconstruction passes and MC:
   * > for an ideally calibrated detector you expect a peak around 1, but for
   * > the "alignment during data-taking" which is what was used for the
   * > BParking Prompt reconstruction you can see it's pretty shifted to
   * > lower values, meaning the error is overestimated.
   * > So I would recommend to scale down the data uncertainties.
   *
   * However, we can't do this because this would cause the impact parameter
   * significance to increase which means that some events which didn't
   * trigger should have triggered. Therefore, we *increase* the error on the
   * MC.
   *
   * This number comes from looking at the impact parameter error for the
   * trigger muon for the B0 -> J/psi K* events. I tested both a scaling and
   * an offset and the offset seemed to match really well, whereas the
   * scaling didn't work. */
  if (isMC_)
    cov(reco::TrackBase::i_dxy, reco::TrackBase::i_dxy) = pow(sqrt(cov(reco::TrackBase::i_dxy, reco::TrackBase::i_dxy)) + 0.0003,2);

  return reco::Track(tk->chi2(), tk->ndof(), tk->referencePoint(), tk->momentum(), tk->charge(), cov, tk->algo(), (reco::TrackBase::TrackQuality) tk->qualityMask());
}


KinematicFitResult 
MesonProducer::fillInfo(pat::CompositeCandidate& v0Cand,
			const edm::Event& iEvent,
			const pat::PackedCandidate & cand1,
			const pat::PackedCandidate & cand2,
			string cand1_name,
			string cand2_name)
{
  std::vector<const reco::Track*> trks;
  std::vector<float> masses;
  trks.push_back(cand1.bestTrack());
  trks.push_back(cand2.bestTrack());
  masses.push_back(cand1.mass());
  masses.push_back(cand2.mass());
  auto vtxFit = vertexWithKinematicFitter(trks, masses);
  if(doFitxTrackCov) vtxFit.postprocess(*beamSpot_, possibleVtxs_);
  else  vtxFit.postprocess(*beamSpot_, *primaryVertices_);

  // printf("vtxFit (x,y,z): (%7.3f,%7.3f,%7.3f)\n", 
  // 	 vtxFit.refitVertex->position().x(),
  // 	 vtxFit.refitVertex->position().y(),
  // 	 vtxFit.refitVertex->position().z());
  addFitInfo(v0Cand, vtxFit, "kin");

  int pvIndex = -1;
  if(doFitxTrackCov) {
    auto displacement3D = compute3dDisplacement(vtxFit, possibleVtxs_ ,true);
    pvIndex = displacement3D.pvIndex;
  } else {
    auto displacement3D = compute3dDisplacement(vtxFit, *primaryVertices_ ,true);
    pvIndex = displacement3D.pvIndex;
  }

  // 0 is charged Isolation; 22 is photon, -1 is neuHad
  v0Cand.addUserFloat( "iso", computeCandIsolation(cand1,cand2,pvIndex,       0.9,0.3,0)); //minPt and DR=0.3 as for muons
  v0Cand.addUserFloat( "isoPho", computeCandIsolation(cand1,cand2,pvIndex,    0.9,0.3,22)); //minPt and DR=0.3 as for muons
  v0Cand.addUserFloat( "isoNeuHad", computeCandIsolation(cand1,cand2,pvIndex, 0.9,0.3,-1)); //minPt and DR=0.3 as for muons
  v0Cand.addUserFloat( "phoPtSumOutsideSignalConedR03", computeCandIsolation(cand1,cand2,pvIndex, 0.9,0.3,22, true)); //minPt and DR=0.3 as for muons

  if( vtxFit.valid() ) {

    const auto & vtx_point = vtxFit.refitVertex->vertexState().position();
    v0Cand.addUserFloat( "bestVtx_X", vtx_point.x());
    v0Cand.addUserFloat( "bestVtx_Y", vtx_point.y());
    v0Cand.addUserFloat( "bestVtx_Z", vtx_point.z());
    v0Cand.addUserInt( "bestVtx_idx", pvIndex);

  }

  if (isMC_){
    auto gen_tt = getGenMatchInfo( cand1, cand2 );
    if (gen_tt.mc_trk1){
      v0Cand.addUserInt(  "gen_" + cand1_name + "_pdgId",   gen_tt.mc_trk1->pdgId() );
      v0Cand.addUserInt(  "gen_" + cand1_name + "_mpdgId",  gen_tt.mc_trk1->numberOfMothers()>0?gen_tt.mc_trk1->mother(0)->pdgId():0 );
      v0Cand.addUserFloat("gen_" + cand1_name + "_pt",      gen_tt.mc_trk1->pt() );
    } else {
      v0Cand.addUserInt(  "gen_" + cand1_name + "_pdgId",   0 );
      v0Cand.addUserInt(  "gen_" + cand1_name + "_mpdgId",  0 );
      v0Cand.addUserFloat("gen_" + cand1_name + "_pt",      0 );
    }
    if (gen_tt.mc_trk2){
      v0Cand.addUserInt(  "gen_" + cand2_name + "_pdgId",   gen_tt.mc_trk2->pdgId() );
      v0Cand.addUserInt(  "gen_" + cand2_name + "_mpdgId",  gen_tt.mc_trk2->numberOfMothers()>0?gen_tt.mc_trk2->mother(0)->pdgId():0 );
      v0Cand.addUserFloat("gen_" + cand2_name + "_pt",      gen_tt.mc_trk2->pt() );
    } else {
      v0Cand.addUserInt(  "gen_" + cand2_name + "_pdgId",   0 );
      v0Cand.addUserInt(  "gen_" + cand2_name + "_mpdgId",  0 );
      v0Cand.addUserFloat("gen_" + cand2_name + "_pt",      0 );
    }
    if (gen_tt.match){
      v0Cand.addUserFloat("gen_mass",         gen_tt.match->mass() );
      v0Cand.addUserFloat("gen_pt",           gen_tt.match->pt() );
      v0Cand.addUserInt(  "gen_pdgId",        gen_tt.match->pdgId() );
    } else {
      v0Cand.addUserFloat("gen_mass",         0 );
      v0Cand.addUserFloat("gen_pt",           0 );
      v0Cand.addUserInt(  "gen_pdgId",        0 );
    }
  }
  return vtxFit;
}

namespace {
  typedef ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3> > cov33_t;
  typedef ROOT::Math::SMatrix<double,6,6,ROOT::Math::MatRepSym<double,6> > cov66_t;
  typedef ROOT::Math::SMatrix<double,7,7,ROOT::Math::MatRepSym<double,7> > cov77_t;
  typedef ROOT::Math::SMatrix<double,9,9,ROOT::Math::MatRepSym<double,9> > cov99_t;
  typedef ROOT::Math::SVector<double,9> jac9_t;

  cov33_t GlobalError2SMatrix_33(GlobalError m_in)
  {
    cov33_t m_out;
    for (int i=0; i<3; i++) {
      for (int j=i; j<3; j++)  {
	m_out(i,j) = m_in.matrix()(i,j);
      }
    }
    return m_out;
  }

  cov99_t makeCovarianceMatrix(const cov33_t cov_vtx1,
			       const cov77_t cov_vtx2)
  {
    cov99_t cov;
    cov.Place_at(cov_vtx1,0,0);
    cov.Place_at(cov_vtx2.Sub<cov66_t>(0,0),3,3);
    return cov;
  }

  jac9_t makeJacobianVector2d(const AlgebraicVector3 &vtx1, const AlgebraicVector3 &vtx2,
			      const AlgebraicVector3 &momentum) {
    jac9_t jac;
    const double momentumMag = ROOT::Math::Mag(momentum);
    const AlgebraicVector3 dist = vtx2 - vtx1;
    const double distMag = ROOT::Math::Mag(dist);
    const double factorPositionComponent = 1./(distMag*momentumMag);
    const double factorMomentumComponent = 1./pow(momentumMag,3);
    jac(0)=-dist(0)*factorPositionComponent;
    jac(1)=-dist(1)*factorPositionComponent;
    jac(3)= dist(0)*factorPositionComponent;
    jac(4)= dist(1)*factorPositionComponent;
    jac(6)= momentum(0)*factorMomentumComponent;
    jac(7)= momentum(1)*factorMomentumComponent;
    return jac;
  }

  jac9_t makeJacobianVector2d(const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,
			      ROOT::Math::DefaultCoordinateSystemTag> &vtx1,
			      const GlobalPoint &vtx2, const TVector3 &tv3momentum) {
    return makeJacobianVector2d(AlgebraicVector3(vtx1.X(),vtx1.Y(),vtx1.Z()),
				AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
				AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
  }

  jac9_t makeJacobianVector3d(const AlgebraicVector3 &vtx1,
			      const AlgebraicVector3 &vtx2,
			      const AlgebraicVector3 &momentum)
  {
    jac9_t jac;
    const AlgebraicVector3 dist = vtx2 - vtx1;
    const double factor2 = 1. / ROOT::Math::Mag2(momentum);
    const double lifetime = ROOT::Math::Dot(dist, momentum) * factor2;
    jac.Place_at(-momentum*factor2,0);
    jac.Place_at( momentum*factor2,3);
    jac.Place_at( factor2*(dist-2*lifetime*momentum*factor2),6);
    return jac;
  }

  jac9_t makeJacobianVector3d(const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,
			      ROOT::Math::DefaultCoordinateSystemTag> &vtx1,
			      const GlobalPoint &vtx2, const TVector3 &tv3momentum)
  {
    return makeJacobianVector3d(AlgebraicVector3(vtx1.X(),vtx1.Y(),vtx1.Z()),
				AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
				AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
  }



}

/*
namespace{
  int match_to_muon(const pat::PackedCandidate& pfCand,
		    const std::vector<pat::Muon>& muons){
    if ( abs(pfCand.pdgId())!=13 ) return -1;
    for ( unsigned int i=0; i<muons.size(); ++i ){
      if ( deltaR(muons.at(i), pfCand) < 0.01 ) return i;
    }
    return -1;
  }
}
*/

void MesonProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    auto const& bFieldHandle_ = iSetup.getHandle(idealMagneticFieldRecordToken_);
    bField_ = bFieldHandle_.product();
    auto const& theTTBuilderHandle_ = iSetup.getHandle(theTTBuilderToken_);
    theTTBuilder_ = theTTBuilderHandle_.product();

    AnalyticalImpactPointExtrapolator extrapolator(bField_);
    impactPointExtrapolator_ = &extrapolator;

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(beamSpotToken_, beamSpotHandle);
    if ( ! beamSpotHandle.isValid() ) {
        edm::LogError("MesonProducer") << "No beam spot available from EventSetup" ;
    }
    beamSpot_ = beamSpotHandle.product();

    edm::Handle<reco::VertexCollection> pvHandle;
    iEvent.getByToken(vertexToken_, pvHandle);
    primaryVertices_ = pvHandle.product();

    edm::Handle<pat::PackedCandidateCollection > pfCandHandle_;
    //    iEvent.getByToken(muonToken_, muonHandle_);
    iEvent.getByToken(pfCandToken_, pfCandHandle_);
    pfCands_ = pfCandHandle_.product();

    edm::Handle<std::vector<pat::PackedGenParticle> > packedGenParticleHandle;
    if ( isMC_ ) {
      iEvent.getByToken(packedGenToken_,packedGenParticleHandle);
      packedGenParticles_ = packedGenParticleHandle.product();
    } else {
      packedGenParticles_ = nullptr;
    }

    // auto nMuons = muonHandle_->size();
    auto nPFCands = pfCands_->size();
    
    // REFIT VERTICES as in
    // https://github.com/ocerri/BPH_RDntuplizer/blob/1fef8bb05ebf8eb99a3d75343197eb72da52ab83/plugins/B2DstMuDecayTreeProducer.cc#L142
    for(uint i_vtx = 0; i_vtx<pvHandle->size(); i_vtx++) {
      auto vtx = (*pvHandle)[i_vtx];
      /* Cut on the number of degrees of freedom. The number of degrees of
       * freedom is calculated as:
       *
       *     ndof = -3 + 2*(w_1 + w_2 + ...)
       *
       *  where w_1, w_2, etc. are weights for each track associated with the
       *  vertex. Quoting from https://arxiv.org/pdf/1405.6569.pdf:
       *
       *  > In the adaptive vertex fit, each track in the vertex is assigned
       *  > a weight between 0 and 1, which reflects the likelihood that it
       *  > genuinely > belongs to the vertex.
       *
       *  Therefore, the number of degrees of freedom is strongly correlated
       *  with the number of tracks. We cut here requiring more than 4 which
       *  translates roughly to at least 4 "good" tracks. This is because for
       *  ndof <= 4 the data and MC do not seem to agree, and there is an
       *  excess of vertices with low ndof in data. */
      if (vtx.ndof() <= 4) continue;
      reco::Vertex tmp = refit_vertex(iEvent, iSetup, i_vtx, 1);
      if (tmp.isValid()) possibleVtxs_.push_back(tmp);
      // else {
      // cout << "[WARNING] Invalid vertex refit for " << i_vtx << endl;
      // possibleVtxs.push_back(vtx);
      // }
    }

    // Output collection
    auto kss  = std::make_unique<pat::CompositeCandidateCollection>();
    auto d0s  = std::make_unique<pat::CompositeCandidateCollection>();
    auto d0pi0s  = std::make_unique<pat::CompositeCandidateCollection>();
    auto k0Stars = std::make_unique<pat::CompositeCandidateCollection>();
    auto phis = std::make_unique<pat::CompositeCandidateCollection>();
    auto lambdas = std::make_unique<pat::CompositeCandidateCollection>();
    auto rhos = std::make_unique<pat::CompositeCandidateCollection>();
    auto omegas = std::make_unique<pat::CompositeCandidateCollection>();

    // Build V0 candidates first

    if ( nPFCands > 1 ){
      for ( unsigned int i = 0; i < nPFCands-1; ++i ) {
	const pat::PackedCandidate& pfCand1( (*pfCandHandle_)[i] );
	if ( not isGoodTrack(pfCand1) ) continue;
	for ( unsigned int j = i+1; j < nPFCands; ++j ) {
	  const pat::PackedCandidate& pfCand2( (*pfCandHandle_)[j] );
	  if ( pfCand1.charge()*pfCand2.charge() >= 0 ) continue;
	  if ( not isGoodTrack(pfCand2) ) continue;
	  if ( not isGoodPair(pfCand1,pfCand2) ) continue;

	  auto tt_doca = distanceOfClosestApproach(pfCand1.bestTrack(),
						   pfCand2.bestTrack());
	  if ( maxTwoTrackDOCA_>0 and tt_doca > maxTwoTrackDOCA_ )
	    continue;
	  
	  // $$$$$$$$$$
	  // $$$$$$ BELOW diTrack

	  // PhiToKK (1020) [ also Look for DsToPhiPi Ds (1968./27)]
	  auto phiCand = getPhiToKK(iEvent, pfCand1, pfCand2);
	  if (phiCand.numberOfDaughters() > 0 and phiCand.pt() > 5 ){
	    phiCand.addUserFloat( "doca", tt_doca);
	    phis->push_back(phiCand);
	  }

	  // RhoToPiPi (750)
	  auto rhoCand = getRhosToPiPi(iEvent, pfCand1, pfCand2);
	  if (rhoCand.numberOfDaughters() > 0 and rhoCand.pt() > 5 ){
	    rhoCand.addUserFloat( "doca", tt_doca);
	    rhos->push_back(rhoCand);
	  }

	  // K0starToKPi (892)
	  auto k0StarCand1 = getK0StarToKPi(iEvent, pfCand1, pfCand2);
	  if (k0StarCand1.numberOfDaughters() > 0 and k0StarCand1.pt() > 5){
	    k0StarCand1.addUserFloat( "doca", tt_doca);
	    k0Stars->push_back(k0StarCand1);
	  }
	  auto k0StarCand2 = getK0StarToKPi(iEvent, pfCand2, pfCand1);
	  if (k0StarCand2.numberOfDaughters() > 0 and k0StarCand2.pt() > 5){
	    k0StarCand2.addUserFloat( "doca", tt_doca);
	    k0Stars->push_back(k0StarCand2);
	  }

	  // $$$$$$$$$$
	  // $$$$$$ BELOW diTrack + pi0/gammas

	  // OmegaToPiPiPi0 (782)
	  auto omegaCand = getOmegasToPiPiPi0(iEvent, pfCand1, pfCand2);
	  if (omegaCand.numberOfDaughters() > 0 and omegaCand.pt() > 5 ){
	    omegaCand.addUserFloat( "doca", tt_doca);
	    omegas->push_back(omegaCand);
	  }

	  //D0->K-pi+: Br=3.95% <-- this is the target (doing both for the anti-particle, to pick up one per event)
	  //D0->K+pi-: Br=0.015%

	  // D0ToKPi (1864) D0*ToD0Gammas (2007)
	  auto d0Cand1 = getD0ToKPi(iEvent, pfCand1, pfCand2);
	  if (d0Cand1.numberOfDaughters() > 0 and d0Cand1.pt() > 5){
	    d0Cand1.addUserFloat( "doca", tt_doca);
	    d0s->push_back(d0Cand1);
	  }
	  auto d0Cand2 = getD0ToKPi(iEvent, pfCand2, pfCand1);
	  if (d0Cand2.numberOfDaughters() > 0 and d0Cand2.pt() > 5){
	    d0Cand2.addUserFloat( "doca", tt_doca);
	    d0s->push_back(d0Cand2);
	  }

	  //D0->K-pi+pi0: Br=14% <-- this is the target
	  // D0ToKPi (1864) D0*ToD0Gammas (2007)
	  auto d0pi0Cand1 = getD0ToKPiPi0(iEvent, pfCand1, pfCand2);
	  if (d0pi0Cand1.numberOfDaughters() > 0 and d0pi0Cand1.pt() > 5){
	    d0pi0Cand1.addUserFloat( "doca", tt_doca);
	    d0pi0s->push_back(d0pi0Cand1);
	  }
	  auto d0pi0Cand2 = getD0ToKPiPi0(iEvent, pfCand2, pfCand1);
	  if (d0pi0Cand2.numberOfDaughters() > 0 and d0pi0Cand2.pt() > 5){
	    d0pi0Cand2.addUserFloat( "doca", tt_doca);
	    d0pi0s->push_back(d0pi0Cand2);
	  }

	  // Look for V0s built from displaced tracks

	  if ( not displacedTrack(pfCand1) or 
	       not displacedTrack(pfCand2) ) continue;

	  // KsToPiPi (493)
	  auto ksCand = getKsToPiPi(iEvent, pfCand1, pfCand2);
	  if (ksCand.numberOfDaughters() > 0){
	    ksCand.addUserFloat( "doca", tt_doca);
	    kss->push_back(ksCand);
	  }

	  // LambdaToPPi (1520)
	  auto lambdaCand1 = getLambdaToPPi(iEvent, pfCand1, pfCand2);
	  if (lambdaCand1.numberOfDaughters() > 0){
	    lambdaCand1.addUserFloat( "doca", tt_doca);
	    lambdas->push_back(lambdaCand1);
	  }
	  auto lambdaCand2 = getLambdaToPPi(iEvent, pfCand2, pfCand1);
	  if (lambdaCand2.numberOfDaughters() > 0){
	    lambdaCand2.addUserFloat( "doca", tt_doca);
	    lambdas->push_back(lambdaCand2);
	  }

	}
      }
    }
    
    iEvent.put(std::move(kss), "Ks");
    iEvent.put(std::move(d0s), "D0");
    iEvent.put(std::move(d0pi0s), "D0pi0");
    iEvent.put(std::move(k0Stars), "K0Star");
    iEvent.put(std::move(phis), "Phi");
    iEvent.put(std::move(lambdas), "Lambda");
    iEvent.put(std::move(rhos), "Rho");
    iEvent.put(std::move(omegas), "Omega");
}

pat::CompositeCandidate
MesonProducer::getKsToPiPi(const edm::Event& iEvent,
			   const pat::PackedCandidate& ipfCand1,
			   const pat::PackedCandidate& ipfCand2)
{
  pat::CompositeCandidate ksCand;
  pat::PackedCandidate pfCand1(ipfCand1);
  pfCand1.setMass(pion_mass_);
  pat::PackedCandidate pfCand2(ipfCand2);
  pfCand2.setMass(pion_mass_);
  ksCand.addDaughter( pfCand1 , "trk1" );
  ksCand.addDaughter( pfCand2 , "trk2" );
  AddFourMomenta addP4;
  addP4.set( ksCand );

  if ( ksCand.mass() < minKsPreselectMass_ or ksCand.mass() > maxKsPreselectMass_ )
    return pat::CompositeCandidate();

  ksCand.addUserFloat( "trk1_pt",  pfCand1.pt() );
  ksCand.addUserFloat( "trk1_eta", pfCand1.eta() );
  ksCand.addUserFloat( "trk1_phi", pfCand1.phi() );
  ksCand.addUserFloat( "trk2_pt",  pfCand2.pt() );
  ksCand.addUserFloat( "trk2_eta", pfCand2.eta() );
  ksCand.addUserFloat( "trk2_phi", pfCand2.phi() );
  ksCand.addUserFloat( "trk1_sip", trackImpactParameterSignificance(pfCand1) );
  ksCand.addUserFloat( "trk2_sip", trackImpactParameterSignificance(pfCand2) );
  //  ksCand.addUserInt( "trk1_mu_index", match_to_muon(pfCand1,*muonHandle_));
  //  ksCand.addUserInt( "trk2_mu_index", match_to_muon(pfCand2,*muonHandle_));

  auto ksVtxFit = fillInfo(ksCand, iEvent, pfCand1, pfCand2);

  if ( not ksVtxFit.valid()  or 
       ksVtxFit.vtxProb() < minVtxProb_ or
       ksVtxFit.mass() < minKsMass_ or  
       ksVtxFit.mass() > maxKsMass_ or
       ksVtxFit.lxy > maxLxy_ or 
       ksVtxFit.sigLxy < minSigLxy_ ) return pat::CompositeCandidate();

  return ksCand;
}


pat::CompositeCandidate
MesonProducer::getRhosToPiPi(const edm::Event& iEvent,
			     const pat::PackedCandidate& ipfCand1,
			     const pat::PackedCandidate& ipfCand2)
{
  pat::CompositeCandidate rhosCand;
  pat::PackedCandidate pfCand1(ipfCand1);
  pfCand1.setMass(pion_mass_);
  pat::PackedCandidate pfCand2(ipfCand2);
  pfCand2.setMass(pion_mass_);
  rhosCand.addDaughter( pfCand1 , "trk1" );
  rhosCand.addDaughter( pfCand2 , "trk2" );
  AddFourMomenta addP4;
  addP4.set( rhosCand );

  if ( rhosCand.mass() < minRhosPreselectMass_ or rhosCand.mass() > maxRhosPreselectMass_ )
    return pat::CompositeCandidate();

  rhosCand.addUserFloat( "trk1_pt",  pfCand1.pt() );
  rhosCand.addUserFloat( "trk1_eta", pfCand1.eta() );
  rhosCand.addUserFloat( "trk1_phi", pfCand1.phi() );
  rhosCand.addUserFloat( "trk2_pt",  pfCand2.pt() );
  rhosCand.addUserFloat( "trk2_eta", pfCand2.eta() );
  rhosCand.addUserFloat( "trk2_phi", pfCand2.phi() );
  rhosCand.addUserFloat( "trk1_sip", trackImpactParameterSignificance(pfCand1) );
  rhosCand.addUserFloat( "trk2_sip", trackImpactParameterSignificance(pfCand2) );
  //  ksCand.addUserInt( "trk1_mu_index", match_to_muon(pfCand1,*muonHandle_));
  //  ksCand.addUserInt( "trk2_mu_index", match_to_muon(pfCand2,*muonHandle_));

  auto rhosVtxFit = fillInfo(rhosCand, iEvent, pfCand1, pfCand2);

  if ( not rhosVtxFit.valid()  or
       rhosVtxFit.vtxProb() < minVtxProb_ or
       rhosVtxFit.mass() < minRhosMass_ or
       rhosVtxFit.mass() > maxRhosMass_ ) return pat::CompositeCandidate();

  return rhosCand;
}


pat::CompositeCandidate
MesonProducer::getOmegasToPiPiPi0(const edm::Event& iEvent,
				  const pat::PackedCandidate& ipfCand1,
				  const pat::PackedCandidate& ipfCand2)
{
  float signalCone_ = 0.05; // to collect the pi0
  float ptPhotonSignal_ = 5.;

  pat::CompositeCandidate omegasCand;
  pat::PackedCandidate pfCand1(ipfCand1);
  pfCand1.setMass(pion_mass_);
  pat::PackedCandidate pfCand2(ipfCand2);
  pfCand2.setMass(pion_mass_);
  omegasCand.addDaughter( pfCand1 , "trk1" );
  omegasCand.addDaughter( pfCand2 , "trk2" );
  AddFourMomenta addP4;
  addP4.set( omegasCand );

  // this is done to preselect the two charged Pis
  if ( omegasCand.mass() < minOmegasPreselectMass_ or omegasCand.mass() > maxOmegasPreselectMass_ )
    return pat::CompositeCandidate();

  omegasCand.addUserFloat( "trk1_pt",  pfCand1.pt() );
  omegasCand.addUserFloat( "trk1_eta", pfCand1.eta() );
  omegasCand.addUserFloat( "trk1_phi", pfCand1.phi() );
  omegasCand.addUserFloat( "trk2_pt",  pfCand2.pt() );
  omegasCand.addUserFloat( "trk2_eta", pfCand2.eta() );
  omegasCand.addUserFloat( "trk2_phi", pfCand2.phi() );
  omegasCand.addUserFloat( "trk1_sip", trackImpactParameterSignificance(pfCand1) );
  omegasCand.addUserFloat( "trk2_sip", trackImpactParameterSignificance(pfCand2) );

  auto omegasVtxFit = fillInfo(omegasCand, iEvent, pfCand1, pfCand2);

  if ( not omegasVtxFit.valid() or
       omegasVtxFit.vtxProb() < minVtxProb_ or
       omegasVtxFit.mass() < minOmegasPreselectMass_ or
       omegasVtxFit.mass() > maxOmegasPreselectMass_
       ) return pat::CompositeCandidate();

  // look for a Photon
  int nPhotons = 0;
  float photon_pt_ = 0.f;
  float photon_eta_ = 0.f;
  float photon_phi_ = 0.f;
  float photon2_pt_ = 0.f;
  float photon2_eta_ = 0.f;
  float photon2_phi_ = 0.f;

  const auto & vtx_point = omegasVtxFit.refitVertex->vertexState().position();
  pat::CompositeCandidate omegaFullCand;
  omegaFullCand.addDaughter( pfCand1 , "trk1" );
  omegaFullCand.addDaughter( pfCand2 , "trk2" );

  for (const auto& iphoton: *pfCands_){

    if (iphoton.charge() != 0 ) continue;
    if (iphoton.mass() > 1 ) continue; // some uninitialized mass
    if (abs(iphoton.pdgId()) !=22) continue; // otherwise some "KL" enters
    if (iphoton.pt() < ptPhotonSignal_ ) continue;
    if (deltaR((pfCand1.p4()+pfCand2.p4()),iphoton.p4()) > signalCone_) continue; // photon should be collimated

    pat::PackedCandidate photon(iphoton);
    photon.setMass(0.);
    photon.setVertex(math::XYZPoint(vtx_point.x(), vtx_point.y(), vtx_point.z()));

    // Omega: aim to get two photons from Pi0
    nPhotons++;
    if (nPhotons==1) omegaFullCand.addDaughter( photon , "photon" );
    if (nPhotons==2) omegaFullCand.addDaughter( photon , "photon2" );

  } // if Photon

  addP4.set(omegaFullCand);

  if ( omegaFullCand.mass() < minRhosPreselectMass_ or omegaFullCand.mass() > maxRhosPreselectMass_ )
    return pat::CompositeCandidate();

  TLorentzVector tot_4;
  tot_4.SetPtEtaPhiM(omegasVtxFit.p3().perp(),omegasVtxFit.p3().eta(),omegasVtxFit.p3().phi(),omegasVtxFit.mass());

  unsigned int ndau = omegaFullCand.numberOfDaughters();

  if (ndau>2) {
    const reco::Candidate* ph1 = omegaFullCand.daughter("photon");
    photon_pt_ = ph1->pt();
    photon_eta_ = ph1->eta();
    photon_phi_ = ph1->phi();
    TLorentzVector ph1_;
    ph1_.SetPtEtaPhiM(photon_pt_, photon_eta_, photon_phi_, 0.);
    tot_4 = tot_4 + ph1_;
  }
  if (ndau>3) {
    const reco::Candidate*  ph2 = omegaFullCand.daughter("photon2");
    photon2_pt_ = ph2->pt();
    photon2_eta_ = ph2->eta();
    photon2_phi_ = ph2->phi();
    TLorentzVector ph2_;
    ph2_.SetPtEtaPhiM(photon2_pt_, photon2_eta_, photon2_phi_, 0.);
    tot_4 = tot_4 + ph2_;
  }

  omegasCand.addUserInt( "Nphotons", nPhotons );
  omegasCand.addUserFloat( "Nbody_mass", tot_4.M() );
  omegasCand.addUserFloat( "Nbody_pt", tot_4.Pt() );
  omegasCand.addUserFloat( "Nbody_eta", tot_4.Eta() );
  omegasCand.addUserFloat( "Nbody_phi", tot_4.Phi() );

  omegasCand.addUserFloat( "photon_pt", photon_pt_ );
  omegasCand.addUserFloat( "photon_eta", photon_eta_);
  omegasCand.addUserFloat( "photon_phi", photon_phi_);
  omegasCand.addUserFloat( "photon2_pt", photon2_pt_ );
  omegasCand.addUserFloat( "photon2_eta", photon2_eta_);
  omegasCand.addUserFloat( "photon2_phi", photon2_phi_);

  return omegasCand;

}

pat::CompositeCandidate
MesonProducer::getK0StarToKPi(const edm::Event& iEvent,
			      const pat::PackedCandidate& ikaon,
			      const pat::PackedCandidate& ipion)
{
  pat::CompositeCandidate Cand;
  pat::PackedCandidate kaon(ikaon);
  kaon.setMass(kaon_mass_);
  pat::PackedCandidate pion(ipion);
  pion.setMass(pion_mass_);
  Cand.addDaughter( kaon , "kaon" );
  Cand.addDaughter( pion , "pion" );
  AddFourMomenta addP4;
  addP4.set( Cand );

  if ( Cand.mass() < minK0StarPreselectMass_ or Cand.mass() > maxK0StarPreselectMass_ )
    return pat::CompositeCandidate();

  Cand.addUserFloat( "kaon_pt",  kaon.pt() );
  Cand.addUserFloat( "kaon_eta", kaon.eta() );
  Cand.addUserFloat( "kaon_phi", kaon.phi() );
  Cand.addUserFloat( "kaon_charge", kaon.charge() );
  Cand.addUserFloat( "pion_pt",  pion.pt() );
  Cand.addUserFloat( "pion_eta", pion.eta() );
  Cand.addUserFloat( "pion_phi", pion.phi() );
  Cand.addUserFloat( "pion_charge", pion.charge() );
  Cand.addUserFloat( "kaon_sip", trackImpactParameterSignificance(kaon) );
  Cand.addUserFloat( "pion_sip", trackImpactParameterSignificance(pion) );
  //  d0Cand.addUserInt( "kaon_mu_index", match_to_muon(kaon, *muonHandle_));
  //  d0Cand.addUserInt( "pion_mu_index", match_to_muon(pion, *muonHandle_));

  auto candVtxFit = fillInfo(Cand, iEvent, kaon, pion, "kaon", "pion");

  if ( not candVtxFit.valid()  or
       candVtxFit.vtxProb() < minVtxProb_ or
       candVtxFit.mass() < minK0StarMass_ or
       candVtxFit.mass() > maxK0StarMass_ ) return pat::CompositeCandidate();
  //       candVtxFit.lxy > maxLxy_ or
  //       candVtxFit.sigLxy < minSigLxy_ ) return pat::CompositeCandidate();

  return Cand;
}

pat::CompositeCandidate
MesonProducer::getD0ToKPi(const edm::Event& iEvent,
		       const pat::PackedCandidate& ikaon,
		       const pat::PackedCandidate& ipion)
{
  float ptPhotonSignal_ = 5.;
  float signalCone_ = 0.05; // to collect the pi0

  pat::CompositeCandidate d0Cand;
  pat::PackedCandidate kaon(ikaon);
  kaon.setMass(kaon_mass_);
  pat::PackedCandidate pion(ipion);
  pion.setMass(pion_mass_);
  d0Cand.addDaughter( kaon , "kaon" );
  d0Cand.addDaughter( pion , "pion" );
  AddFourMomenta addP4;
  addP4.set( d0Cand );
  
  if ( d0Cand.mass() < minD0PreselectMass_ or d0Cand.mass() > maxD0PreselectMass_ )
    return pat::CompositeCandidate();

  d0Cand.addUserFloat( "kaon_pt",  kaon.pt() );
  d0Cand.addUserFloat( "kaon_eta", kaon.eta() );
  d0Cand.addUserFloat( "kaon_phi", kaon.phi() );
  d0Cand.addUserFloat( "kaon_charge", kaon.charge() );
  d0Cand.addUserFloat( "pion_pt",  pion.pt() );
  d0Cand.addUserFloat( "pion_eta", pion.eta() );
  d0Cand.addUserFloat( "pion_phi", pion.phi() );
  d0Cand.addUserFloat( "pion_charge", pion.charge() );
  d0Cand.addUserFloat( "kaon_sip", trackImpactParameterSignificance(kaon) );
  d0Cand.addUserFloat( "pion_sip", trackImpactParameterSignificance(pion) );
  //  d0Cand.addUserInt( "kaon_mu_index", match_to_muon(kaon, *muonHandle_));
  //  d0Cand.addUserInt( "pion_mu_index", match_to_muon(pion, *muonHandle_));

  auto d0VtxFit = fillInfo(d0Cand, iEvent, kaon, pion, "kaon", "pion");
	  
  if ( not d0VtxFit.valid()  or
       d0VtxFit.vtxProb() < minVtxProb_ or
       d0VtxFit.mass() < minD0Mass_ or
       d0VtxFit.mass() > maxD0Mass_ ) return pat::CompositeCandidate();
  //       d0VtxFit.lxy > maxLxy_ or
  //       d0VtxFit.sigLxy < minSigLxy_ ) return pat::CompositeCandidate();

  // look for a Photon
  //  const pat::PackedCandidate* d0Star_photon(nullptr);
  int nPhotons = 0;
  float photon_pt_ = 0.f;
  float photon_eta_ = 0.f;
  float photon_phi_ = 0.f;
  float photon2_pt_ = 0.f;
  float photon2_eta_ = 0.f;
  float photon2_phi_ = 0.f;

  const auto & vtx_point = d0VtxFit.refitVertex->vertexState().position();
  pat::CompositeCandidate d0StarCand;
  d0StarCand.addDaughter( kaon , "kaon" );
  d0StarCand.addDaughter( pion , "pion" );

  for (const auto& iphoton: *pfCands_){

    if (iphoton.charge() != 0 ) continue;
    if (iphoton.mass() > 1 ) continue; // some uninitialized mass
    if (abs(iphoton.pdgId()) !=22) continue; // otherwise some "KL" enters
    if (iphoton.pt() < ptPhotonSignal_ ) continue;
    if (deltaR((kaon.p4()+pion.p4()),iphoton.p4()) > signalCone_) continue; // photon should be collimated

    pat::PackedCandidate photon(iphoton);
    photon.setMass(0.);
    photon.setVertex(math::XYZPoint(vtx_point.x(), vtx_point.y(), vtx_point.z()));

    // D0->k/pi; there the two photons are aiming to get one or two photons from Pi0 from the D0*->D0
    nPhotons++;
    if (nPhotons==1) d0StarCand.addDaughter( photon , "photon" );
    if (nPhotons==2) d0StarCand.addDaughter( photon , "photon2" );
  }

  addP4.set( d0StarCand);

  // keep the lower bound to catch the cases with 0 photons
  if ( d0StarCand.mass() < minD0Mass_ or d0StarCand.mass() > maxD0StarMass_ )
    return pat::CompositeCandidate();

  TLorentzVector tot_4;
  tot_4.SetPtEtaPhiM(d0VtxFit.p3().perp(),d0VtxFit.p3().eta(),d0VtxFit.p3().phi(),d0VtxFit.mass());

  unsigned int ndau = d0StarCand.numberOfDaughters();

  if (ndau>2) {
    const reco::Candidate* ph1 = d0StarCand.daughter("photon");
    photon_pt_ = ph1->pt();
    photon_eta_ = ph1->eta();
    photon_phi_ = ph1->phi();
    TLorentzVector ph1_;
    ph1_.SetPtEtaPhiM(photon_pt_, photon_eta_, photon_phi_, 0.);
    tot_4 = tot_4 + ph1_;
  }
  if (ndau>3) {
    const reco::Candidate*  ph2 = d0StarCand.daughter("photon2");
    photon2_pt_ = ph2->pt();
    photon2_eta_ = ph2->eta();
    photon2_phi_ = ph2->phi();
    TLorentzVector ph2_;
    ph2_.SetPtEtaPhiM(photon2_pt_, photon2_eta_, photon2_phi_, 0.);
    tot_4 = tot_4 + ph2_;
  }

  d0Cand.addUserFloat( "d0Star_Nbody_mass", tot_4.M() );
  d0Cand.addUserFloat( "d0Star_Nbody_pt", tot_4.Pt() );
  d0Cand.addUserFloat( "d0Star_Nbody_eta", tot_4.Eta() );
  d0Cand.addUserFloat( "d0Star_Nbody_phi", tot_4.Phi() );

  d0Cand.addUserInt( "d0Star_Nphotons", nPhotons );
  d0Cand.addUserFloat( "d0Star_photon_pt", photon_pt_ );
  d0Cand.addUserFloat( "d0Star_photon_eta", photon_eta_ );
  d0Cand.addUserFloat( "d0Star_photon_phi", photon_phi_ );
  d0Cand.addUserFloat( "d0Star_photon2_pt", photon2_pt_ );
  d0Cand.addUserFloat( "d0Star_photon2_eta", photon2_eta_ );
  d0Cand.addUserFloat( "d0Star_photon2_phi", photon2_phi_ );

  return d0Cand;
}

pat::CompositeCandidate
MesonProducer::getD0ToKPiPi0(const edm::Event& iEvent,
			     const pat::PackedCandidate& ikaon,
			     const pat::PackedCandidate& ipion)
{

  float ptPhotonSignal_ = 5.;
  float signalCone_ = 0.1; // to collect the pi0
  pat::CompositeCandidate d0Cand;
  pat::PackedCandidate kaon(ikaon);
  kaon.setMass(kaon_mass_);
  pat::PackedCandidate pion(ipion);
  pion.setMass(pion_mass_);
  d0Cand.addDaughter( kaon , "kaon" );
  d0Cand.addDaughter( pion , "pion" );
  AddFourMomenta addP4;
  addP4.set( d0Cand );

  if ( d0Cand.mass() < minD0pi0PreselectMass_ or d0Cand.mass() > maxD0pi0PreselectMass_ )
    return pat::CompositeCandidate();

  d0Cand.addUserFloat( "kaon_pt",  kaon.pt() );
  d0Cand.addUserFloat( "kaon_eta", kaon.eta() );
  d0Cand.addUserFloat( "kaon_phi", kaon.phi() );
  d0Cand.addUserFloat( "kaon_charge", kaon.charge() );
  d0Cand.addUserFloat( "pion_pt",  pion.pt() );
  d0Cand.addUserFloat( "pion_eta", pion.eta() );
  d0Cand.addUserFloat( "pion_phi", pion.phi() );
  d0Cand.addUserFloat( "pion_charge", pion.charge() );
  d0Cand.addUserFloat( "kaon_sip", trackImpactParameterSignificance(kaon) );
  d0Cand.addUserFloat( "pion_sip", trackImpactParameterSignificance(pion) );
  //  d0Cand.addUserInt( "kaon_mu_index", match_to_muon(kaon, *muonHandle_));
  //  d0Cand.addUserInt( "pion_mu_index", match_to_muon(pion, *muonHandle_));

  auto d0VtxFit = fillInfo(d0Cand, iEvent, kaon, pion, "kaon", "pion");

  if ( not d0VtxFit.valid()  or
       d0VtxFit.vtxProb() < minVtxProb_ or
       d0VtxFit.mass() < minD0pi0PreselectMass_ or
       d0VtxFit.mass() > maxD0pi0PreselectMass_ ) return pat::CompositeCandidate();
  //       d0VtxFit.lxy > maxLxy_ or
  //       d0VtxFit.sigLxy < minSigLxy_ ) return pat::CompositeCandidate();

  // look for a Photon
  //  const pat::PackedCandidate* d0Star_photon(nullptr);
  int nPhotons = 0;
  float photon_pt_ = 0.f;
  float photon_eta_ = 0.f;
  float photon_phi_ = 0.f;
  float photon2_pt_ = 0.f;
  float photon2_eta_ = 0.f;
  float photon2_phi_ = 0.f;

  const auto & vtx_point = d0VtxFit.refitVertex->vertexState().position();
  pat::CompositeCandidate d0StarCand;
  d0StarCand.addDaughter( kaon , "kaon" );
  d0StarCand.addDaughter( pion , "pion" );

  for (const auto& iphoton: *pfCands_){

    if (iphoton.charge() != 0 ) continue;
    if (iphoton.mass() > 1 ) continue; // some uninitialized mass
    if (abs(iphoton.pdgId()) !=22) continue; // otherwise some "KL" enters
    if (iphoton.pt() < ptPhotonSignal_ ) continue;
    if (deltaR((kaon.p4()+pion.p4()),iphoton.p4()) > signalCone_) continue; // photon should be collimated

    pat::PackedCandidate photon(iphoton);
    photon.setMass(0.);
    photon.setVertex(math::XYZPoint(vtx_point.x(), vtx_point.y(), vtx_point.z()));

    nPhotons++;
    if (nPhotons==1) d0StarCand.addDaughter( photon , "photon" );
    if (nPhotons==2) d0StarCand.addDaughter( photon , "photon2" );
  }

  addP4.set( d0StarCand);

  TLorentzVector tot_4;
  tot_4.SetPtEtaPhiM(d0VtxFit.p3().perp(),d0VtxFit.p3().eta(),d0VtxFit.p3().phi(),d0VtxFit.mass());

  if ( d0StarCand.mass() < minD0StarMass_ and d0StarCand.mass() > maxD0StarMass_ )
    return pat::CompositeCandidate();

  unsigned int ndau = d0StarCand.numberOfDaughters();

  if (ndau>2) {
    const reco::Candidate* ph1 = d0StarCand.daughter("photon");
    photon_pt_ = ph1->pt();
    photon_eta_ = ph1->eta();
    photon_phi_ = ph1->phi();
    TLorentzVector ph1_;
    ph1_.SetPtEtaPhiM(photon_pt_, photon_eta_, photon_phi_, 0.);
    tot_4 = tot_4 + ph1_;
  }
  if (ndau>3) {
    const reco::Candidate*  ph2 = d0StarCand.daughter("photon2");
    photon2_pt_ = ph2->pt();
    photon2_eta_ = ph2->eta();
    photon2_phi_ = ph2->phi();
    TLorentzVector ph2_;
    ph2_.SetPtEtaPhiM(photon2_pt_, photon2_eta_, photon2_phi_, 0.);
    tot_4 = tot_4 + ph2_;
  }

  d0Cand.addUserFloat( "d0Star_Nbody_mass", tot_4.M() );
  d0Cand.addUserFloat( "d0Star_Nbody_pt", tot_4.Pt() );
  d0Cand.addUserFloat( "d0Star_Nbody_eta", tot_4.Eta() );
  d0Cand.addUserFloat( "d0Star_Nbody_phi", tot_4.Phi() );

  d0Cand.addUserInt( "d0Star_Nphotons", nPhotons );
  d0Cand.addUserFloat( "d0Star_photon_pt", photon_pt_ );
  d0Cand.addUserFloat( "d0Star_photon_eta", photon_eta_ );
  d0Cand.addUserFloat( "d0Star_photon_phi", photon_phi_ );
  d0Cand.addUserFloat( "d0Star_photon2_pt", photon2_pt_ );
  d0Cand.addUserFloat( "d0Star_photon2_eta", photon2_eta_ );
  d0Cand.addUserFloat( "d0Star_photon2_phi", photon2_phi_ );

  return d0Cand;

}


pat::CompositeCandidate
MesonProducer::getPhiToKK(const edm::Event& iEvent,
			  const pat::PackedCandidate& ipfCand1,
			  const pat::PackedCandidate& ipfCand2)
{
  pat::CompositeCandidate phiCand;
  pat::PackedCandidate pfCand1(ipfCand1);
  pfCand1.setMass(kaon_mass_);
  pat::PackedCandidate pfCand2(ipfCand2);
  pfCand2.setMass(kaon_mass_);
  phiCand.addDaughter( pfCand1 , "trk1" );
  phiCand.addDaughter( pfCand2 , "trk2" );
  AddFourMomenta addP4;
  addP4.set( phiCand );

  if ( phiCand.mass() < minPhiPreselectMass_ or phiCand.mass() > maxPhiPreselectMass_ )
    return pat::CompositeCandidate();

  phiCand.addUserFloat( "trk1_pt",  pfCand1.pt() );
  phiCand.addUserFloat( "trk1_eta", pfCand1.eta() );
  phiCand.addUserFloat( "trk1_phi", pfCand1.phi() );
  phiCand.addUserFloat( "trk2_pt",  pfCand2.pt() );
  phiCand.addUserFloat( "trk2_eta", pfCand2.eta() );
  phiCand.addUserFloat( "trk2_phi", pfCand2.phi() );
  phiCand.addUserFloat( "trk1_sip", trackImpactParameterSignificance(pfCand1) );
  phiCand.addUserFloat( "trk2_sip", trackImpactParameterSignificance(pfCand2) );
  //  phiCand.addUserInt( "trk1_mu_index", match_to_muon(pfCand1,*muonHandle_));
  //  phiCand.addUserInt( "trk2_mu_index", match_to_muon(pfCand2,*muonHandle_));

  auto ksVtxFit = fillInfo(phiCand, iEvent, pfCand1, pfCand2);

  if ( not ksVtxFit.valid()  or 
       ksVtxFit.vtxProb() < minVtxProb_ or
       ksVtxFit.mass() < minPhiMass_ or  
       ksVtxFit.mass() > maxPhiMass_ ) return pat::CompositeCandidate();

  // Look for DsToPhiPi
  // Keep the best candidate by vertex probability if there are multiple
  KinematicFitResult dsVtx;
  const pat::PackedCandidate* ds_pion(nullptr);
  for (const pat::PackedCandidate& ipion: *pfCands_){
    if (&ipion == &ipfCand1 or &ipion == &ipfCand2) continue;
    if (not isGoodTrack(ipion) ) continue;
    if (not isGoodPion(ipion) ) continue;
    pat::CompositeCandidate dsCand;
    pat::PackedCandidate pion(ipion);
    pion.setMass(pion_mass_);
    dsCand.addDaughter( pfCand1 , "trk1" );
    dsCand.addDaughter( pfCand2 , "trk2" );
    dsCand.addDaughter( pion , "pion" );
    addP4.set( dsCand );
    if ( dsCand.mass() < minDsPreselectMass_ or
	 dsCand.mass() > maxDsPreselectMass_ ) continue;

    // fit 3-body vertex
    std::vector<const reco::Track*> trks;
    std::vector<float> masses;
    trks.push_back(pfCand1.bestTrack());
    trks.push_back(pfCand2.bestTrack());
    trks.push_back(pion.bestTrack());
    masses.push_back(pfCand1.mass());
    masses.push_back(pfCand2.mass());
    masses.push_back(pion.mass());

    auto vtxFit = vertexWithKinematicFitter(trks, masses);
    if(doFitxTrackCov) vtxFit.postprocess(*beamSpot_, possibleVtxs_);
    else  vtxFit.postprocess(*beamSpot_, *primaryVertices_);

    if ( vtxFit.valid()  and vtxFit.vtxProb() > minVtxProb_ and
	 vtxFit.mass() > minDsMass_ and vtxFit.mass() < maxDsMass_ )
      {
	if (not dsVtx.valid() or vtxFit.vtxProb() > dsVtx.vtxProb()){
	  dsVtx = vtxFit;
	  ds_pion = &ipion;
	}
      }
  }
  if (ds_pion){
    phiCand.addUserFloat( "ds_pion_pt", ds_pion->pt() );
    phiCand.addUserFloat( "ds_pion_eta", ds_pion->eta() );
    phiCand.addUserFloat( "ds_pion_phi", ds_pion->phi() );
    //    phiCand.addUserInt(   "ds_pion_mu_index", match_to_muon(*ds_pion, *muonHandle_));
  } else {
    phiCand.addUserFloat( "ds_pion_pt", 0.0 );
    phiCand.addUserFloat( "ds_pion_eta", 0.0 );
    phiCand.addUserFloat( "ds_pion_phi", 0.0 );
    phiCand.addUserInt(   "ds_pion_mu_index", -1);
  }
  addFitInfo(phiCand, dsVtx, "ds");
    
  return phiCand;
}

pat::CompositeCandidate
MesonProducer::getLambdaToPPi(const edm::Event& iEvent,
			   const pat::PackedCandidate& iproton,
			   const pat::PackedCandidate& ipion)
{
  pat::CompositeCandidate lambdaCand;
  pat::PackedCandidate proton( iproton );
  proton.setMass( proton_mass_ );
  pat::PackedCandidate pion( ipion );
  pion.setMass( pion_mass_ );
  lambdaCand.addDaughter( proton, "proton" );
  lambdaCand.addDaughter( pion, "pion" );
  AddFourMomenta addP4;
  addP4.set( lambdaCand );

  if ( lambdaCand.mass() < minLambdaPreselectMass_ or lambdaCand.mass() > maxLambdaPreselectMass_ )
    return pat::CompositeCandidate();

  lambdaCand.addUserFloat( "proton_pt",  proton.pt() );
  lambdaCand.addUserFloat( "proton_eta", proton.eta() );
  lambdaCand.addUserFloat( "proton_phi", proton.phi() );
  lambdaCand.addUserFloat( "pion_pt",  pion.pt() );
  lambdaCand.addUserFloat( "pion_eta", pion.eta() );
  lambdaCand.addUserFloat( "pion_phi", pion.phi() );
  lambdaCand.addUserFloat( "proton_sip", trackImpactParameterSignificance(proton) );
  lambdaCand.addUserFloat( "pion_sip", trackImpactParameterSignificance(pion) );
  //  lambdaCand.addUserInt( "proton_mu_index", match_to_muon(proton,*muonHandle_));
  //  lambdaCand.addUserInt( "pion_mu_index", match_to_muon(pion,*muonHandle_));

  auto ksVtxFit = fillInfo(lambdaCand, iEvent, proton, pion, "proton", "pion");
	  
  if ( not ksVtxFit.valid()  or 
       ksVtxFit.vtxProb() < minVtxProb_ or
       ksVtxFit.mass() < minLambdaMass_ or  
       ksVtxFit.mass() > maxLambdaMass_ or
       ksVtxFit.lxy > maxLxy_ or 
       ksVtxFit.sigLxy < minSigLxy_
       ) return pat::CompositeCandidate();

  return lambdaCand;
}

DisplacementInformationIn3D
  MesonProducer::compute3dDisplacement(const KinematicFitResult& fit,
				       const reco::VertexCollection& vertices,
				       bool closestIn3D)
{
  DisplacementInformationIn3D result;
  if (not fit.valid()) return result;

  // Potential issue: tracks used to build the candidate could
  // also be used in the primary vertex fit. One can refit the vertices
  // excluding tracks from the cadndidate. It's not done at the moment
  // due to non-trivial linkig between primary vertex and its tracks
  // in MiniAOD. Also not all muons associated to a vertex are really
  // used in the fit, so the potential bias most likely small.

  auto candTransientTrack = fit.refitMother->refittedTransientTrack();

  const reco::Vertex* bestVertex(0);
  int bestVertexIndex(-1);
  double minDistance(999.);

  for ( unsigned int i = 0; i<vertices.size(); ++i ){
    const auto & vertex = vertices.at(i);
    if (closestIn3D){
      auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, vertex);
      if (impactParameter3D.first and impactParameter3D.second.value() < minDistance){
	minDistance = impactParameter3D.second.value();
	bestVertex = &vertex;
	bestVertexIndex = i;
      }
    } else{
      auto impactParameterZ  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), vertex);
      double distance = fabs(impactParameterZ.second.value());
      if (impactParameterZ.first and distance < minDistance){
	minDistance = distance;
	bestVertex = &vertex;
	bestVertexIndex = i;
      }
    }
  }

  // find second best vertex
  const reco::Vertex* bestVertex2(0);
  int bestVertexIndex2(-1);
  double minDistance2(999.);
  for ( unsigned int i = 0; i<vertices.size(); ++i ){
    const auto & vertex = vertices.at(i);
    if (closestIn3D){
      auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, vertex);
      if (impactParameter3D.first and impactParameter3D.second.value() < minDistance2 and impactParameter3D.second.value() > minDistance){
	minDistance2 = impactParameter3D.second.value();
	bestVertex2 = &vertex;
	bestVertexIndex2 = i;
      }
    } else{
      auto impactParameterZ  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), vertex);
      double distance = fabs(impactParameterZ.second.value());
      if (impactParameterZ.first and distance < minDistance2 and distance > minDistance){
	minDistance2 = distance;
	bestVertex2 = &vertex;
	bestVertexIndex2 = i;
      }
    }
  }

  if (! bestVertex) return result;

  auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, *bestVertex);
  auto impactParameterZ  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), *bestVertex);
  result.pv = bestVertex;
  result.pvIndex = bestVertexIndex;
  if (impactParameterZ.first) {
    result.longitudinalImpactParameter    = impactParameterZ.second.value();
    result.longitudinalImpactParameterSig = impactParameterZ.second.significance();
    result.longitudinalImpactParameterErr = impactParameterZ.second.error();
  }
  if (impactParameter3D.first and not isnan(impactParameter3D.second.error())) {
    result.distaceOfClosestApproach       = impactParameter3D.second.value();
    result.distaceOfClosestApproachSig    = impactParameter3D.second.significance();
    result.distaceOfClosestApproachErr    = impactParameter3D.second.error();
  }

  // compute decay length
  VertexDistance3D distance3D;
  auto dist = distance3D.distance(*bestVertex, fit.refitVertex->vertexState() );
  result.decayLength    = dist.value();
  result.decayLengthErr = dist.error();

  VertexDistanceXY distanceXY;
  auto distXY = distanceXY.distance(*bestVertex, fit.refitVertex->vertexState() );

  if (bestVertex2){
    auto impactParameter3D2 = IPTools::absoluteImpactParameter3D(candTransientTrack, *bestVertex2);
    auto impactParameterZ2  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), *bestVertex2);
    result.pv2 = bestVertex2;
    result.pv2Index = bestVertexIndex2;
    if (impactParameterZ2.first) {
      result.longitudinalImpactParameter2    = impactParameterZ2.second.value();
      result.longitudinalImpactParameter2Sig = impactParameterZ2.second.significance();
      result.longitudinalImpactParameter2Err = impactParameterZ2.second.error();
    }
    if (impactParameter3D2.first) {
      result.distaceOfClosestApproach2       = impactParameter3D2.second.value();
      result.distaceOfClosestApproach2Sig    = impactParameter3D2.second.value();
      result.distaceOfClosestApproach2Err    = impactParameter3D2.second.error();
    }

    // compute decay length
    VertexDistance3D distance3D;
    auto dist = distance3D.distance(*bestVertex2, fit.refitVertex->vertexState() );
    result.decayLength2    = dist.value();
    result.decayLength2Err = dist.error();

  }

  //
  // Pointing angle
  //
  auto alpha = getAlpha(fit.refitVertex->vertexState().position(),
			fit.refitVertex->vertexState().error(),
			GlobalPoint(Basic3DVector<float>(bestVertex->position())),
			GlobalError(bestVertex->covariance()),
			fit.refitMother->currentState().globalMomentum());

  auto alphaXY = getAlpha(fit.refitVertex->vertexState().position(),
			  fit.refitVertex->vertexState().error(),
			  GlobalPoint(Basic3DVector<float>(bestVertex->position())),
			  GlobalError(bestVertex->covariance()),
			  fit.refitMother->currentState().globalMomentum(),
			  true);

  result.alpha    = alpha.first;
  result.alphaErr = alpha.second;

  result.alphaXY    = alphaXY.first;
  result.alphaXYErr = alphaXY.second;

  //
  // Decay time information
  //
  TVector3 plab(fit.refitMother->currentState().globalMomentum().x(),
		fit.refitMother->currentState().globalMomentum().y(),
		fit.refitMother->currentState().globalMomentum().z());
  const double massOverC = fit.mass()/TMath::Ccgs();

  // get covariance matrix for error propagation in decayTime calculation
  auto vtxDistanceCov = makeCovarianceMatrix(GlobalError2SMatrix_33(bestVertex->error()),
					     fit.refitMother->currentState().kinematicParametersError().matrix());
  auto vtxDistanceJac3d = makeJacobianVector3d(bestVertex->position(), fit.refitVertex->vertexState().position(), plab);
  auto vtxDistanceJac2d = makeJacobianVector2d(bestVertex->position(), fit.refitVertex->vertexState().position(), plab);

  result.decayTime = dist.value() / plab.Mag() * cos(result.alpha) * massOverC;
  result.decayTimeError = TMath::Sqrt(ROOT::Math::Similarity(vtxDistanceCov, vtxDistanceJac3d)) * massOverC;

  result.decayTimeXY = distXY.value() / plab.Perp() * cos(result.alphaXY) * massOverC;
  result.decayTimeXYError = TMath::Sqrt(ROOT::Math::Similarity(vtxDistanceCov, vtxDistanceJac2d)) * massOverC;

  return result;

}

float
MesonProducer::computeCandIsolation(const pat::PackedCandidate& pfCand1, const pat::PackedCandidate& pfCand2, 
				    unsigned int primaryVertexIndex,
				    float minPt, float dR, int particleType,
				    bool returnSum,
				    std::vector<const pat::PackedCandidate*> ignoreTracks)
{
    float sumPt(0);
    auto b_p4 = pfCand1.p4()+pfCand2.p4();
    for (const auto& pfCandIso: *pfCands_){
      bool ignore_track = false;
      for (auto trk: ignoreTracks){
	if (trk==&pfCandIso){
	  ignore_track = true;
	  break;
	}
      }
      if (ignore_track) continue;

      if (deltaR(pfCand1, pfCandIso) < 0.001 || deltaR(pfCand2, pfCandIso) < 0.001) continue;
      if (particleType==0) {
	if (pfCandIso.charge() == 0 ) continue; // for chargedIsolation
	if (!pfCandIso.hasTrackDetails()) continue;
	if (pfCandIso.vertexRef().key()!=primaryVertexIndex) continue;
      } else if (particleType==22 and pfCandIso.pdgId()==22) {
	if (pfCandIso.charge() != 0 ) continue; // only photons for neutralIsolation
      } else {
	if (pfCandIso.pdgId() == 22) continue; // only neutral hadrons
	if (pfCandIso.charge() != 0 ) continue; // for neutralIsolation
      }
      if (pfCandIso.pt()<minPt) continue;
      if (deltaR(b_p4, pfCandIso) > dR) continue;
      sumPt += pfCandIso.pt();
    }

    if(returnSum) return sumPt;
    return b_p4.pt()/(b_p4.pt()+sumPt);
  }


KinematicFitResult 
MesonProducer::vertexWithKinematicFitter(std::vector<const reco::Track*> trks,
					    std::vector<float> masses)
{
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideKinematicVertexFit
  if ( trks.size() != masses.size() ) 
    throw cms::Exception("Error") << "number of tracks and number of masses should match";

  std::vector<reco::TransientTrack> transTrks;

  KinematicParticleFactoryFromTransientTrack factory;
  KinematicParticleVertexFitter fitter;
    
  std::vector<RefCountedKinematicParticle> particles;

  double chi = 0.;
  double ndf = 0.;
  float mass_err(mass_err_);
  for (unsigned int i=0; i<trks.size(); ++i){
    if(doFitxTrackCov) transTrks.push_back((*theTTBuilder_).build(fix_track(trks[i])));
    else transTrks.push_back((*theTTBuilder_).build(trks[i]));
    particles.push_back(factory.particle(transTrks.back(),masses[i],chi,ndf,mass_err));
  }

  RefCountedKinematicTree vertexFitTree = fitter.fit(particles);
  KinematicFitResult result;
    
  if ( !vertexFitTree->isValid()) return result;
  
  result.treeIsValid = true;

  vertexFitTree->movePointerToTheTop();
  result.refitVertex = vertexFitTree->currentDecayVertex();
  result.refitMother = vertexFitTree->currentParticle();
  result.refitTree   = vertexFitTree;
  if ( !result.refitVertex->vertexIsValid()) return result;

  result.vertexIsValid = true;

  // extract the re-fitted tracks
  vertexFitTree->movePointerToTheTop();
  
  if ( vertexFitTree->movePointerToTheFirstChild() ){
    do {
      result.refitDaughters.push_back(vertexFitTree->currentParticle());
    } while (vertexFitTree->movePointerToTheNextChild());
  }
  return result;
}


/*
pair<double,double> MesonProducer::computeDCA(const pat::PackedCandidate &kaon,
                                                 reco::BeamSpot beamSpot){

  const reco::TransientTrack trackTT((*(kaon.bestTrack())), &(*bFieldHandle_));

  TrajectoryStateClosestToPoint theDCAXBS = trackTT.trajectoryStateClosestToPoint( GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()) );  
  
  double DCABS = theDCAXBS.perigeeParameters().transverseImpactParameter();
  double DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
  pair<double,double> DCA = make_pair(DCABS,DCABSErr);
    
  return DCA;
}
*/
namespace{

  const pat::PackedGenParticle*
    gen_match(const std::vector<pat::PackedGenParticle>& packedGenParticles,
	      const LorentzVector& reco)
  {
    const pat::PackedGenParticle* best_match(nullptr);
    float min_delta = 1e9;
    for (auto const & gen: packedGenParticles){
      if (deltaR(reco,gen) > 0.15) continue;
      float delta = fabs(reco.pt()-gen.pt())/gen.pt();
      if (delta > min_delta) continue;
      min_delta = delta;
      best_match = &gen;
    }
    return best_match;
  }

  std::vector<unsigned int> 
    get_depth_from_permutation(const std::vector<unsigned int>& elements){
    std::vector<unsigned int> result;
    unsigned int counter(0);
    for (auto element: elements){
      if (element==0){
	counter++;
      } else {
	result.push_back(counter);
	counter = 0;
      }
    }
    result.push_back(counter);
    return result;
  }

  bool is_acceptable(const reco::Candidate* cand){
    if ( not cand) return false; 
    // skip quarks
    if ( abs(cand->pdgId())<10 ) return false;
    // skip protons
    if ( abs(cand->pdgId())==2212 ) return false;
    // skip gluons
    if ( abs(cand->pdgId())==21 ) return false;
    return true;
  }

  // depth 0 - first mother

  const reco::Candidate* get_mother(const reco::Candidate* cand, unsigned int depth){
    if (not cand) return 0;
    const reco::Candidate* mother = cand->mother();
    unsigned int i = 0;
    while ( is_acceptable(mother) and i<depth ){
      i++;
      mother = mother->mother();
    }
    if (is_acceptable(mother))
      return mother;
    else
      return 0;
  }

  const reco::Candidate* 
    find_common_ancestor(const std::vector<const reco::Candidate*>& particles, 
			 unsigned int max_depth=10){
    auto n = particles.size();
    for (unsigned int depth=0; depth<max_depth; ++depth){
      // make a list of elements (0) and separators (1) and
      // find all possible permutations of the elements
      std::vector<unsigned int> elements;
      for (unsigned int i=0; i<depth; ++i)
	elements.push_back(0);
      for (unsigned int i=0; i<n-1; ++i)
	elements.push_back(1);
      do {
	auto depth_vector = get_depth_from_permutation(elements);
	const reco::Candidate* common_mother(0);
	for (unsigned int i=0; i<n; ++i){
	  auto mother = get_mother(particles[i],depth_vector[i]);
	  if (not mother) {
	    common_mother = 0;
	    break;
	  }
	  if (not common_mother) common_mother = mother;
	  if (common_mother != mother) {
	    common_mother = 0;
	    break;
	  }	  
	}
	if (common_mother) return common_mother;
      } while(std::next_permutation(elements.begin(), elements.end()));
    }
    return 0;
  }

}

GenMatchInfo MesonProducer::getGenMatchInfo( const pat::PackedCandidate& track1,
					     const pat::PackedCandidate& track2 )
{
  GenMatchInfo result;
  std::vector<const reco::Candidate*> daughters;
  result.mc_trk1 = gen_match(*packedGenParticles_, track1.p4());
  if (result.mc_trk1)
    daughters.push_back(result.mc_trk1);
  result.mc_trk2 = gen_match(*packedGenParticles_, track2.p4());
  if (result.mc_trk2)
    daughters.push_back(result.mc_trk2);
  if (daughters.size()==2){
    const auto* mother = find_common_ancestor(daughters);
    if (mother) result.match        = mother;
  }
  return result;
}

float MesonProducer::distanceOfClosestApproach( const reco::Track* track1,
					     const reco::Track* track2)
{
  TwoTrackMinimumDistance md;
  reco::TransientTrack tt1;
  reco::TransientTrack tt2;

  if(doFitxTrackCov) { 
    tt1 = theTTBuilder_->build(fix_track(track1));
    tt2 = theTTBuilder_->build(fix_track(track2));
  } else {
    tt1 = theTTBuilder_->build(track1);
    tt2 = theTTBuilder_->build(track2);
  }
  if ( not md.calculate( tt1.initialFreeState(), tt2.initialFreeState() ) ) return -1.0;
  return md.distance();
}

Measurement1D 
MesonProducer::distanceOfClosestApproach( const reco::Track* track,
					     RefCountedKinematicVertex vertex)
{
  if (not vertex->vertexIsValid()) return Measurement1D(-1.0,-1.0);
  VertexDistance3D distance3D;
  reco::TransientTrack tt;
  if(doFitxTrackCov) tt = theTTBuilder_->build(fix_track(track)); 
  else tt = theTTBuilder_->build(track);

  assert(impactPointExtrapolator_);
  auto tsos = impactPointExtrapolator_->extrapolate(tt.initialFreeState(), vertex->position());
  if ( not tsos.isValid()) return Measurement1D(-1.0,-1.0);
  Measurement1D doca = distance3D.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), vertex->vertexState());
  return doca;
}

Measurement1D 
MesonProducer::distanceOfClosestApproach( const reco::Track* track,
					     const reco::Vertex& vertex)
{
  VertexDistance3D distance3D;
  reco::TransientTrack tt;
  if(doFitxTrackCov) tt = theTTBuilder_->build(fix_track(track));
  else tt = theTTBuilder_->build(track);

  assert(impactPointExtrapolator_);
  auto tsos = impactPointExtrapolator_->extrapolate(tt.initialFreeState(), GlobalPoint(Basic3DVector<float>(vertex.position())));
  if ( not tsos.isValid()) return Measurement1D(-1.0,-1.0);
  Measurement1D doca = distance3D.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), vertex);
  return doca;
}

std::pair<float, float>
MesonProducer::getAlpha(const GlobalPoint& vtx_position, const GlobalError& vtx_error,
			const GlobalPoint& ip_position,  const GlobalError& ip_error,
			const GlobalVector &momentum,
			bool transverse)
{
  AlgebraicSymMatrix33 error_matrix(vtx_error.matrix() + ip_error.matrix());
  GlobalVector dir(vtx_position - ip_position);
  if (dir.mag() == 0)
    return std::pair<float, float>(999., 999.);

  GlobalVector p(momentum);
  if (transverse){
    dir = GlobalVector(dir.x(), dir.y(), 0);
    p = GlobalVector(p.x(), p.y(), 0);
  }

  double dot_product = dir.dot(p);
  double cosAlpha = dot_product / p.mag() / dir.mag();
  if (cosAlpha > 1) cosAlpha = 1;
  if (cosAlpha < -1) cosAlpha = -1;

  // Error propagation

  double c1 = 1 / dir.mag() / p.mag();
  double c2 = dot_product / pow(dir.mag(), 3) / p.mag();

  double dfdx = p.x() * c1 - dir.x() * c2;
  double dfdy = p.y() * c1 - dir.y() * c2;
  double dfdz = p.z() * c1 - dir.z() * c2;

  double err2_cosAlpha =
    pow(dfdx, 2) * error_matrix(0, 0) +
    pow(dfdy, 2) * error_matrix(1, 1) +
    pow(dfdz, 2) * error_matrix(2, 2) +
    2 * dfdx * dfdy * error_matrix(0, 1) +
    2 * dfdx * dfdz * error_matrix(0, 2) +
    2 * dfdy * dfdz * error_matrix(1, 2);

  float err_alpha = fabs(cosAlpha) <= 1 and err2_cosAlpha >=0 ? sqrt(err2_cosAlpha) / sqrt(1-pow(cosAlpha, 2)) : 999;
  float alpha = acos(cosAlpha);
  if (isnan(alpha) or isnan(err_alpha))
    return std::pair<float, float>(999., 999.);
  else
    return std::pair<float, float>(alpha, err_alpha);
}




DEFINE_FWK_MODULE(MesonProducer);

//  LocalWords:  vertices
