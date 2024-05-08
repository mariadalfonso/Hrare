#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
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
#include "TMVA/Reader.h"

#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include <TLorentzVector.h>
#include <TVector.h>
#include <TMatrix.h>
#include <algorithm>

// 
// DiMuon Producer is designed for J/psi, Y ->mumu analysis
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

class DiMuonProducer : public edm::EDProducer {
    
public:
    
  explicit DiMuonProducer(const edm::ParameterSet &iConfig);
    
  ~DiMuonProducer() override {};
    
    
private:
    
  virtual void produce(edm::Event&, const edm::EventSetup&);

  bool isGoodMuon(const pat::Muon& muon);
  bool isGoodTrack(const pat::PackedCandidate& track);
  bool displacedTrack(const pat::PackedCandidate& track);
  bool isGoodMuonProbe(const pat::PackedCandidate& track);
  bool isGoodPion(const pat::PackedCandidate& track);
  bool isGoodPair(const pat::PackedCandidate& track1,
		  const pat::PackedCandidate& track2);
    
  float
  trackImpactParameterSignificance( const pat::PackedCandidate& track);
  
  pat::CompositeCandidate
  getJpsiToMuMu(const edm::Event& iEvent,
		const pat::Muon& pfCand1,
		const pat::Muon& pfCand2);
  pat::CompositeCandidate
  getUpsilonToMuMu(const edm::Event& iEvent,
		   const pat::Muon& pfCand1,
		   const pat::Muon& pfCand2);
  pat::CompositeCandidate
  getJpsiToTkTk(const edm::Event& iEvent,
		const pat::PackedCandidate& pfCand1,
		const pat::PackedCandidate& pfCand2);


  KinematicFitResult 
  vertexWithKinematicFitter(std::vector<const reco::Track*> trks,
			    std::vector<float> masses);

  pair<double,double> computeDCA(const pat::PackedCandidate &kaon,
   				 reco::BeamSpot beamSpot);
  GenMatchInfo getGenMatchInfo( const pat::Muon& track1,
				const pat::Muon& track2);
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
  computeCandIsolation(	const pat::Muon& pfCand1,
			const pat::Muon& pfCand2,
			unsigned int primaryVertexIndex,
			float minPt=0.9, float dR=0.7,
			int particleType=0,
			unsigned int returnType = 0,
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
	   const pat::Muon & cand1,
	   const pat::Muon & cand2, 
	   string cand1_name = "muon1",
	   string cand2_name = "muon2");

  // ----------member data ---------------------------
    
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  const reco::BeamSpot* beamSpot_;

  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const reco::VertexCollection* primaryVertices_;

  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> pfCandToken_;
  edm::EDGetTokenT<std::vector<pat::PackedGenParticle> >   packedGenToken_;
  const std::vector<pat::PackedGenParticle>* packedGenParticles_;

  edm::ESHandle<TransientTrackBuilder> theTTBuilder_;
  edm::ESHandle<MagneticField> bFieldHandle_;
  edm::Handle<std::vector<pat::PackedCandidate> > pfCandHandle_;
  edm::Handle<std::vector<pat::Muon>> muonHandle_;

  const AnalyticalImpactPointExtrapolator* impactPointExtrapolator_;

  bool isMC_;

  // for J/psi and Upsilon
  double minMuonPt_;
  double maxMuonEta_;
  // for isolation 
  double minPionPt_;
  double maxPionEta_;

  double minJpsisPreselectMass_;
  double maxJpsisPreselectMass_;
  double minJpsisMass_;
  double maxJpsisMass_;

  double minUpsilonsPreselectMass_;
  double maxUpsilonsPreselectMass_;
  double minUpsilonsMass_;
  double maxUpsilonsMass_;

  double maxTwoTrackDOCA_;
  double minDisplaceTrackSignificance_;
  double maxLxy_;
  double minSigLxy_;
  double minCosAlpha_;
  double minVtxProb_;
};

DiMuonProducer::DiMuonProducer(const edm::ParameterSet &iConfig):
beamSpotToken_( consumes<reco::BeamSpot> ( iConfig.getParameter<edm::InputTag>( "beamSpot" ) ) ),
beamSpot_(nullptr),
vertexToken_( consumes<reco::VertexCollection> ( iConfig.getParameter<edm::InputTag>( "vertexCollection" ) ) ),
primaryVertices_(nullptr),
muonToken_( consumes<std::vector<pat::Muon>> ( iConfig.getParameter<edm::InputTag>( "muonCollection" ) ) ),
pfCandToken_( consumes<std::vector<pat::PackedCandidate>> ( iConfig.getParameter<edm::InputTag>( "PFCandCollection" ) ) ),
packedGenToken_( consumes<std::vector<pat::PackedGenParticle>> ( iConfig.getParameter<edm::InputTag>( "packedGenParticleCollection" ) ) ),
packedGenParticles_(nullptr),
impactPointExtrapolator_(0),
isMC_(             iConfig.getParameter<bool>( "isMC" ) ),
minMuonPt_(         iConfig.getParameter<double>( "minMuonPt" ) ),
maxMuonEta_(        iConfig.getParameter<double>( "maxMuonEta" ) ),

minJpsisPreselectMass_(     iConfig.getParameter<double>( "minJpsisPreselectMass" ) ),
maxJpsisPreselectMass_(     iConfig.getParameter<double>( "maxJpsisPreselectMass" ) ),
minJpsisMass_(     iConfig.getParameter<double>( "minJpsisMass" ) ),
maxJpsisMass_(     iConfig.getParameter<double>( "maxJpsisMass" ) ),

minUpsilonsPreselectMass_(     iConfig.getParameter<double>( "minUpsilonsPreselectMass" ) ),
maxUpsilonsPreselectMass_(     iConfig.getParameter<double>( "maxUpsilonsPreselectMass" ) ),
minUpsilonsMass_(     iConfig.getParameter<double>( "minUpsilonsMass" ) ),
maxUpsilonsMass_(     iConfig.getParameter<double>( "maxUpsilonsMass" ) ),

maxTwoTrackDOCA_( iConfig.getParameter<double>( "maxTwoTrackDOCA" ) ),
minDisplaceTrackSignificance_( iConfig.getParameter<double>( "minDisplaceTrackSignificance" ) ),
maxLxy_( iConfig.getParameter<double>( "maxLxy" ) ),
minSigLxy_( iConfig.getParameter<double>( "minSigLxy" ) ),
minCosAlpha_( iConfig.getParameter<double>( "minCosAlpha" ) ),
minVtxProb_( iConfig.getParameter<double>( "minVtxProb" ) )
{
    produces<pat::CompositeCandidateCollection>("Jpsi");
    produces<pat::CompositeCandidateCollection>("Upsilon");
}

bool DiMuonProducer::isGoodTrack(const pat::PackedCandidate& track){
  if ( track.charge() == 0 ) return false;
  if ( not track.hasTrackDetails() ) return false;
  if ( not track.bestTrack()->quality(reco::Track::highPurity) ) return false; 
  return true;
}

float DiMuonProducer::trackImpactParameterSignificance(const pat::PackedCandidate& track){
  return track.bestTrack()->dxyError()>0 ? fabs(track.bestTrack()->dxy(*beamSpot_))/track.bestTrack()->dxyError():0.0;
}

bool DiMuonProducer::displacedTrack(const pat::PackedCandidate& track){
  return trackImpactParameterSignificance(track) > minDisplaceTrackSignificance_;
}

bool DiMuonProducer::isGoodMuon(const pat::Muon& muon){
  if ( not muon.innerTrack().isNonnull()) return false;
  if ( not muon.isLooseMuon() ) return false;
  if ( not muon.isTrackerMuon() ) return false;
  if ( not muon.innerTrack()->quality(reco::Track::highPurity) ) return false; 
  if ( muon.pt() < minMuonPt_ || fabs(muon.eta()) > maxMuonEta_ ) return false;
  return true;
}

bool DiMuonProducer::isGoodPion(const pat::PackedCandidate& track){
  return fabs(track.eta()) < maxPionEta_ and track.pt()>minPionPt_;
}

bool DiMuonProducer::isGoodPair(const pat::PackedCandidate& track1,
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

KinematicFitResult 
DiMuonProducer::fillInfo(pat::CompositeCandidate& v0Cand,
			const edm::Event& iEvent,
			const pat::Muon & cand1,
			const pat::Muon & cand2,
			string cand1_name,
			string cand2_name)
{
  std::vector<const reco::Track*> trks;
  std::vector<float> masses;
  trks.push_back(&(*cand1.innerTrack()));
  trks.push_back(&(*cand2.innerTrack()));
  masses.push_back(cand1.mass());
  masses.push_back(cand2.mass());
  auto vtxFit = vertexWithKinematicFitter(trks, masses);
  vtxFit.postprocess(*beamSpot_, *primaryVertices_);

  // printf("vtxFit (x,y,z): (%7.3f,%7.3f,%7.3f)\n", 
  // 	 vtxFit.refitVertex->position().x(),
  // 	 vtxFit.refitVertex->position().y(),
  // 	 vtxFit.refitVertex->position().z());
  addFitInfo(v0Cand, vtxFit, "kin");

  auto displacement3D = compute3dDisplacement(vtxFit, *primaryVertices_ ,true);
  int pvIndex = displacement3D.pvIndex;

  // particleType: 0 is charged Isolation; 22 is photon, -1 is neuHad
  // returnType: 0 isolation, 1 scalar sum, 2 leading Track Pt, 2 leading Track Dxy
  v0Cand.addUserFloat( "iso", computeCandIsolation(cand1,cand2,pvIndex,      0.9,0.3,0)); //minPt and DR=0.3 as for muons
  v0Cand.addUserFloat( "isoPho", computeCandIsolation(cand1,cand2,pvIndex,   0.9,0.3,22)); //minPt and DR=0.3 as for muons
  v0Cand.addUserFloat( "isoNeuHad", computeCandIsolation(cand1,cand2,pvIndex,0.9,0.3,-1)); //minPt and DR=0.3 as for muons
  v0Cand.addUserFloat( "phoPtSumOutsideSignalConedR03", computeCandIsolation(cand1,cand2,pvIndex, 0.9,0.3,22,  1)); //minPt and DR=0.3 as for muons
  //
  v0Cand.addUserFloat( "leadingChargedMaxPt1", computeCandIsolation(cand1,cand2,pvIndex,      0.9,0.4,0,  2)); //minPt and DR=0.4 as for muons
  v0Cand.addUserFloat( "leadingCharged2dSignMaxPt1", computeCandIsolation(cand1,cand2,pvIndex,0.9,0.4,0,  3)); //minPt and DR=0.4 as for muons
  v0Cand.addUserFloat( "leadingCharged3dSignMaxPt1", computeCandIsolation(cand1,cand2,pvIndex,0.9,0.4,0,  4)); //minPt and DR=0.4 as for muons
  v0Cand.addUserFloat( "leadingChargedDxyMaxPt1", computeCandIsolation(cand1,cand2,pvIndex,   0.9,0.4,0,  5)); //minPt and DR=0.4 as for muons
  v0Cand.addUserFloat( "leadingChargedDzMaxPt1", computeCandIsolation(cand1,cand2,pvIndex,    0.9,0.4,0,  6)); //minPt and DR=0.4 as for muons
  //
  v0Cand.addUserFloat( "leadingChargedMaxPt2", computeCandIsolation(cand1,cand2,pvIndex,      0.9,0.4,0,  7)); //minPt and DR=0.4 as for muons
  v0Cand.addUserFloat( "leadingCharged2dSignMaxPt2", computeCandIsolation(cand1,cand2,pvIndex,0.9,0.4,0,  8)); //minPt and DR=0.4 as for muons
  v0Cand.addUserFloat( "leadingCharged3dSignMaxPt2", computeCandIsolation(cand1,cand2,pvIndex,0.9,0.4,0,  9)); //minPt and DR=0.4 as for muons
  v0Cand.addUserFloat( "leadingChargedDxyMaxPt2", computeCandIsolation(cand1,cand2,pvIndex,   0.9,0.4,0,  10)); //minPt and DR=0.4 as for muons
  v0Cand.addUserFloat( "leadingChargedDzMaxPt2", computeCandIsolation(cand1,cand2,pvIndex,    0.9,0.4,0,  11)); //minPt and DR=0.4 as for muons

  v0Cand.addUserInt( "muon1_isTightMuon", (pvIndex!=-1) ? cand1.isTightMuon((*primaryVertices_).at(pvIndex)): false );
  v0Cand.addUserInt( "muon2_isTightMuon", (pvIndex!=-1) ? cand2.isTightMuon((*primaryVertices_).at(pvIndex)): false );

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

void DiMuonProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle_);
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder_);

    AnalyticalImpactPointExtrapolator extrapolator(bFieldHandle_.product());
    impactPointExtrapolator_ = &extrapolator;

    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(beamSpotToken_, beamSpotHandle);
    if ( ! beamSpotHandle.isValid() ) {
        edm::LogError("DiMuonProducer") << "No beam spot available from EventSetup" ;
    }
    beamSpot_ = beamSpotHandle.product();

    edm::Handle<reco::VertexCollection> pvHandle;
    iEvent.getByToken(vertexToken_, pvHandle);
    primaryVertices_ = pvHandle.product();
    
    iEvent.getByToken(muonToken_, muonHandle_);
    iEvent.getByToken(pfCandToken_, pfCandHandle_);
    
    edm::Handle<std::vector<pat::PackedGenParticle> > packedGenParticleHandle;
    if ( isMC_ ) {
      iEvent.getByToken(packedGenToken_,packedGenParticleHandle);
      packedGenParticles_ = packedGenParticleHandle.product();
    } else {
      packedGenParticles_ = nullptr;
    }

    auto nMuons = muonHandle_->size();
    //    auto nPFCands = pfCandHandle_->size();
    
    // Output collection
    auto jpsis = std::make_unique<pat::CompositeCandidateCollection>();
    auto upsilons = std::make_unique<pat::CompositeCandidateCollection>();

    // Build V0 candidates first from Muon Collections

    if ( nMuons > 1 ){
      for ( unsigned int i = 0; i < nMuons-1; ++i ) {
	const pat::Muon& patMuon1( (*muonHandle_)[i] );
	if ( not isGoodMuon(patMuon1) ) continue;

	for ( unsigned int j = i+1; j < nMuons; ++j ) {
	  const pat::Muon& patMuon2( (*muonHandle_)[j] );
	  if ( not isGoodMuon(patMuon2) ) continue;
	  if ( patMuon1.charge()*patMuon2.charge() >= 0 ) continue;

	  //	  if ( not isGoodTrack(pfCand2) ) continue;
	  //	  if ( not isGoodPair(pfCand1,pfCand2) ) continue;

	  /*
	  if (patMuon1->innerTrack().isNonnull()) {
            trackref = muon->innerTrack();
	  }
	  const reco::Track *trk = &(*trackref);
	  */


	  // two track Doca
	  auto tt_doca = distanceOfClosestApproach(&(*patMuon1.innerTrack()),
						   &(*patMuon2.innerTrack()));
	  if ( maxTwoTrackDOCA_>0 and tt_doca > maxTwoTrackDOCA_ )
	    continue;
	  
	  auto jpsiCand = getJpsiToMuMu(iEvent, patMuon1, patMuon2);
	  if (jpsiCand.numberOfDaughters() > 0){
	    jpsiCand.addUserFloat( "doca", tt_doca);
	    jpsis->push_back(jpsiCand);
	  }

	  // Upsilon
	  auto upsilonCand = getUpsilonToMuMu(iEvent, patMuon1, patMuon2);
	  if (upsilonCand.numberOfDaughters() > 0){
	    upsilonCand.addUserFloat( "doca", tt_doca);
	    upsilons->push_back(upsilonCand);
	  }

	}
      }
    }

    /*
    // Build V0 candidates first from PFCand Collections

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
	  
	  // J/psi
	  auto jpsiCand = getJpsiToTkTk(iEvent, pfCand1, pfCand2);
	  if (jpsiCand.numberOfDaughters() > 0){
	    jpsiCand.addUserFloat( "doca", tt_doca);
	    jpsis->push_back(jpsiCand);
	  }

	}
      }
    }
    */
    
    iEvent.put(std::move(jpsis), "Jpsi");
    iEvent.put(std::move(upsilons), "Upsilon");

}

pat::CompositeCandidate
DiMuonProducer::getJpsiToMuMu(const edm::Event& iEvent,
			     const pat::Muon& iMuon1,
			     const pat::Muon& iMuon2)
{
  pat::CompositeCandidate diMusCand;
  pat::Muon pfCand1(iMuon1);
  pfCand1.setMass(muon_mass_);
  pat::Muon pfCand2(iMuon2);
  pfCand2.setMass(muon_mass_);
  diMusCand.addDaughter( pfCand1 , "muon1" );
  diMusCand.addDaughter( pfCand2 , "muon2" );
  AddFourMomenta addP4;
  addP4.set( diMusCand );

  if ( diMusCand.mass() < minJpsisPreselectMass_ or diMusCand.mass() > maxJpsisPreselectMass_ )
    return pat::CompositeCandidate();

  diMusCand.addUserFloat( "muon1_pt",  pfCand1.pt() );
  diMusCand.addUserFloat( "muon1_eta", pfCand1.eta() );
  diMusCand.addUserFloat( "muon1_phi", pfCand1.phi() );
  diMusCand.addUserFloat( "muon2_pt",  pfCand2.pt() );
  diMusCand.addUserFloat( "muon2_eta", pfCand2.eta() );
  diMusCand.addUserFloat( "muon2_phi", pfCand2.phi() );
  diMusCand.addUserInt( "muon1_isMediumMuon", pfCand1.isMediumMuon() );
  diMusCand.addUserInt( "muon2_isMediumMuon", pfCand2.isMediumMuon() );

  //  diMusCand.addUserFloat( "muon1_sip", trackImpactParameterSignificance(pfCand1) );
  //  diMusCand.addUserFloat( "muon2_sip", trackImpactParameterSignificance(pfCand2) );
  //  ksCand.addUserInt( "trk1_mu_index", match_to_muon(pfCand1,*muonHandle_));
  //  ksCand.addUserInt( "trk2_mu_index", match_to_muon(pfCand2,*muonHandle_));

  auto diMusVtxFit = fillInfo(diMusCand, iEvent, pfCand1, pfCand2);

  if ( not diMusVtxFit.valid()  or
       diMusVtxFit.vtxProb() < minVtxProb_ or
       diMusVtxFit.mass() < minJpsisMass_ or
       diMusVtxFit.mass() > maxJpsisMass_ ) return pat::CompositeCandidate();

  return diMusCand;
}


pat::CompositeCandidate
DiMuonProducer::getUpsilonToMuMu(const edm::Event& iEvent,
				 const pat::Muon& iMuon1,
				 const pat::Muon& iMuon2)
{
  pat::CompositeCandidate diMusCand;
  pat::Muon pfCand1(iMuon1);
  pfCand1.setMass(muon_mass_);
  pat::Muon pfCand2(iMuon2);
  pfCand2.setMass(muon_mass_);
  diMusCand.addDaughter( pfCand1 , "muon1" );
  diMusCand.addDaughter( pfCand2 , "muon2" );
  AddFourMomenta addP4;
  addP4.set( diMusCand );

  if ( diMusCand.mass() < minUpsilonsPreselectMass_ or diMusCand.mass() > maxUpsilonsPreselectMass_ )
    return pat::CompositeCandidate();

  diMusCand.addUserFloat( "muon1_pt",  pfCand1.pt() );
  diMusCand.addUserFloat( "muon1_eta", pfCand1.eta() );
  diMusCand.addUserFloat( "muon1_phi", pfCand1.phi() );
  diMusCand.addUserFloat( "muon2_pt",  pfCand2.pt() );
  diMusCand.addUserFloat( "muon2_eta", pfCand2.eta() );
  diMusCand.addUserFloat( "muon2_phi", pfCand2.phi() );
  diMusCand.addUserInt( "muon1_isMediumMuon", pfCand1.isMediumMuon() );
  diMusCand.addUserInt( "muon2_isMediumMuon", pfCand2.isMediumMuon() );
  //  diMusCand.addUserFloat( "muon1_sip", trackImpactParameterSignificance(pfCand1) );
  //  diMusCand.addUserFloat( "muon2_sip", trackImpactParameterSignificance(pfCand2) );
  //  ksCand.addUserInt( "trk1_mu_index", match_to_muon(pfCand1,*muonHandle_));
  //  ksCand.addUserInt( "trk2_mu_index", match_to_muon(pfCand2,*muonHandle_));

  auto diMusVtxFit = fillInfo(diMusCand, iEvent, pfCand1, pfCand2);

  if ( not diMusVtxFit.valid()  or
       diMusVtxFit.vtxProb() < minVtxProb_ or
       diMusVtxFit.mass() < minUpsilonsMass_ or
       diMusVtxFit.mass() > maxUpsilonsMass_ ) return pat::CompositeCandidate();

  return diMusCand;
}

//
//
//

DisplacementInformationIn3D
  DiMuonProducer::compute3dDisplacement(const KinematicFitResult& fit,
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
DiMuonProducer::computeCandIsolation(const pat::Muon& muon1, const pat::Muon& muon2, 
				     unsigned int primaryVertexIndex,
				     float minPt, float dR, int particleType,
				     unsigned int returnType,
				     std::vector<const pat::PackedCandidate*> ignoreTracks)
{
    float sumPt(0);

    float maxPt1(0);
    float dxyForMaxPt1(-999.);
    float dzForMaxPt1(-999.);
    float signedIP2dSignMaxPt1(-999.);
    float signedIP3dSignMaxPt1(-999.);

    float maxPt2(0);
    float dxyForMaxPt2(-999.);
    float dzForMaxPt2(-999.);
    float signedIP2dSignMaxPt2(-999.);
    float signedIP3dSignMaxPt2(-999.);

    auto b_p4 = muon1.p4()+muon2.p4();
    auto b_dir = GlobalVector(b_p4.x(),b_p4.y(),b_p4.z());
    for (const auto& pfCandIso: *pfCandHandle_.product()){
      bool ignore_track = false;
      for (auto trk: ignoreTracks){
	if (trk==&pfCandIso){
	  ignore_track = true;
	  break;
	}
      }

      if (ignore_track) continue;
      if (deltaR(muon1, pfCandIso) < 0.001 || deltaR(muon2, pfCandIso) < 0.001) continue;
      if (pfCandIso.pt()<minPt) continue;
      if (deltaR(b_p4, pfCandIso) > dR) continue;
      if (particleType==0) {
	if (pfCandIso.charge() == 0 ) continue; // for charged quantities
	if (!pfCandIso.hasTrackDetails()) continue;
	//	if (pfCandIso.vertexRef().key()!=primaryVertexIndex) continue; // check what happens for displaced
	if (pfCandIso.pt() > maxPt1 ) {
	  maxPt1 = pfCandIso.pt();
	  if (returnType==4) {
	    const reco::TransientTrack trackCand((*(pfCandIso.bestTrack())), &(*bFieldHandle_));
	    auto signedIP3d = IPTools::signedImpactParameter3D(trackCand, b_dir, (*primaryVertices_).at(primaryVertexIndex));
	    signedIP3dSignMaxPt1 = signedIP3d.second.significance();
	  } else if (returnType==3){
	    const reco::TransientTrack trackCand((*(pfCandIso.bestTrack())), &(*bFieldHandle_));
	    auto signedIP2d = IPTools::signedTransverseImpactParameter(trackCand, b_dir, (*primaryVertices_).at(primaryVertexIndex));
	    signedIP2dSignMaxPt1 = signedIP2d.second.significance();
	  } else if (returnType==5 or returnType==6){
	    auto vtx_point = (*primaryVertices_).at(primaryVertexIndex).position();
	    dxyForMaxPt1 = pfCandIso.dxy(vtx_point);
	    dzForMaxPt1 = pfCandIso.dz(vtx_point);
	  }
	} else if (pfCandIso.pt() > maxPt2) {
	  maxPt2 = pfCandIso.pt();
	  if (returnType==9) {
	    const reco::TransientTrack trackCand((*(pfCandIso.bestTrack())), &(*bFieldHandle_));
	    auto signedIP3d = IPTools::signedImpactParameter3D(trackCand, b_dir, (*primaryVertices_).at(primaryVertexIndex));
	    signedIP3dSignMaxPt2 = signedIP3d.second.significance();
	  } else if (returnType==8){
	    const reco::TransientTrack trackCand((*(pfCandIso.bestTrack())), &(*bFieldHandle_));
	    auto signedIP2d = IPTools::signedTransverseImpactParameter(trackCand, b_dir, (*primaryVertices_).at(primaryVertexIndex));
	    signedIP2dSignMaxPt2 = signedIP2d.second.significance();
	  } else if (returnType==10 or returnType==11){
	    auto vtx_point = (*primaryVertices_).at(primaryVertexIndex).position();
	    dxyForMaxPt2 = pfCandIso.dxy(vtx_point);
	    dzForMaxPt2 = pfCandIso.dz(vtx_point);
	  }
	}
      } else if (particleType==22 and pfCandIso.pdgId()==22) {
	if (pfCandIso.charge() != 0 ) continue; // only photons for neutralIsolation
      } else {
	if (pfCandIso.pdgId() == 22) continue; // only neutral hadrons
	if (pfCandIso.charge() != 0 ) continue; // for neutralIsolation
      }
      sumPt += pfCandIso.pt();
    }

    if(returnType==1) return sumPt;
    //
    else if (returnType==2) return maxPt1;
    else if (returnType==3) return signedIP2dSignMaxPt1;
    else if (returnType==4) return signedIP3dSignMaxPt1;
    else if (returnType==5) return dxyForMaxPt1;
    else if (returnType==6) return dzForMaxPt1;
    //
    else if (returnType==7) return maxPt2;
    else if (returnType==8) return signedIP2dSignMaxPt2;
    else if (returnType==9) return signedIP3dSignMaxPt2;
    else if (returnType==10) return dxyForMaxPt2;
    else if (returnType==11) return dzForMaxPt2;
    return b_p4.pt()/(b_p4.pt()+sumPt);
}


KinematicFitResult 
DiMuonProducer::vertexWithKinematicFitter(std::vector<const reco::Track*> trks,
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
    transTrks.push_back((*theTTBuilder_).build(trks[i]));
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


pair<double,double> DiMuonProducer::computeDCA(const pat::PackedCandidate &kaon,
                                                 reco::BeamSpot beamSpot){

  const reco::TransientTrack trackTT((*(kaon.bestTrack())), &(*bFieldHandle_));

  TrajectoryStateClosestToPoint theDCAXBS = trackTT.trajectoryStateClosestToPoint( GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()) );  
  
  double DCABS = theDCAXBS.perigeeParameters().transverseImpactParameter();
  double DCABSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
  pair<double,double> DCA = make_pair(DCABS,DCABSErr);
    
  return DCA;
}
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

GenMatchInfo DiMuonProducer::getGenMatchInfo( const pat::Muon& track1,
					     const pat::Muon& track2 )
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

float DiMuonProducer::distanceOfClosestApproach( const reco::Track* track1,
					     const reco::Track* track2)
{
  TwoTrackMinimumDistance md;
  const reco::TransientTrack tt1 = theTTBuilder_->build(track1);
  const reco::TransientTrack tt2 = theTTBuilder_->build(track2);
  if ( not md.calculate( tt1.initialFreeState(), tt2.initialFreeState() ) ) return -1.0;
  return md.distance();
}

Measurement1D 
DiMuonProducer::distanceOfClosestApproach( const reco::Track* track,
					     RefCountedKinematicVertex vertex)
{
  if (not vertex->vertexIsValid()) return Measurement1D(-1.0,-1.0);
  VertexDistance3D distance3D;
  const reco::TransientTrack tt = theTTBuilder_->build(track);
  assert(impactPointExtrapolator_);
  auto tsos = impactPointExtrapolator_->extrapolate(tt.initialFreeState(), vertex->position());
  if ( not tsos.isValid()) return Measurement1D(-1.0,-1.0);
  Measurement1D doca = distance3D.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), vertex->vertexState());
  return doca;
}

Measurement1D 
DiMuonProducer::distanceOfClosestApproach( const reco::Track* track,
					     const reco::Vertex& vertex)
{
  VertexDistance3D distance3D;
  const reco::TransientTrack tt = theTTBuilder_->build(track);
  assert(impactPointExtrapolator_);
  auto tsos = impactPointExtrapolator_->extrapolate(tt.initialFreeState(), GlobalPoint(Basic3DVector<float>(vertex.position())));
  if ( not tsos.isValid()) return Measurement1D(-1.0,-1.0);
  Measurement1D doca = distance3D.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), vertex);
  return doca;
}

std::pair<float, float>
DiMuonProducer::getAlpha(const GlobalPoint& vtx_position, const GlobalError& vtx_error,
			 const GlobalPoint& ip_position,  const GlobalError& ip_error,
			 const GlobalVector &momentum,
			 bool transverse){
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

DEFINE_FWK_MODULE(DiMuonProducer);

//  LocalWords:  vertices
