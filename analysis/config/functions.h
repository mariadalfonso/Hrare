#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include <TTree.h>
#include <TDirectory.h>
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TH2Poly.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"
#include "TVector2.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib> //as stdlib.h      
#include <cstdio>
#include <cmath>
#include <array>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <algorithm>
//#include <boost/algorithm/string/join.hpp>
//#include <boost/algorithm/string.hpp>
//#include <boost/functional/hash.hpp>
#include <limits>
#include <map>

//#include <ROOT/RVec.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
//#include <ROOT/RDF/RInterface.hxx>

using Vec_b = ROOT::VecOps::RVec<bool>;
using Vec_d = ROOT::VecOps::RVec<double>;
using Vec_f = ROOT::VecOps::RVec<float>;
using Vec_i = ROOT::VecOps::RVec<int>;
using Vec_ui = ROOT::VecOps::RVec<unsigned int>;

using stdVec_i = std::vector<int>;
using stdVec_b = std::vector<bool>;
using stdVec_f = std::vector<float>;

typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<float> > XYZVectorF;
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > PtEtaPhiMVector;
std::unordered_map< UInt_t, std::vector< std::pair<UInt_t,UInt_t> > > jsonMap;

float massTheo = 0.892;

Vec_b cleaningPair(const ULong64_t event, const UInt_t nK0Star,
		   Vec_f K0Star_kaon_phi, Vec_f K0Star_kaon_eta , Vec_f K0Star_kaon_pt,
		   Vec_f K0Star_pion_phi, Vec_f K0Star_pion_eta , Vec_f K0Star_pion_pt,
		   Vec_f K0Star_kin_mass) {
  std::vector<std::pair<unsigned int, unsigned int>> idxPair;

  // at least one combination
  if (nK0Star > 1) {
    // find the pairs
    for (unsigned int i = 0; i < nK0Star ; i++){
      for (unsigned int j = i+1; j < nK0Star ; j++){
	if (K0Star_kaon_eta[j]==K0Star_pion_eta[i] && K0Star_kaon_eta[i]==K0Star_pion_eta[j]) {
	  // found i and j that are the same
	  idxPair.push_back({i,j});
	}
      }
    }
  }

  // Result vector to store selected elements, pick one element from each pair
  std::vector<unsigned int> selectedElements;
  for (const auto& pair : idxPair) {
    if (abs(K0Star_kin_mass[pair.first] - massTheo) < abs(K0Star_kin_mass[pair.second]-massTheo)) {
      selectedElements.push_back(pair.second);
    } else {
      selectedElements.push_back(pair.first);
    }
  }

  Vec_b mask(nK0Star, true);
  for (unsigned int idx = 0; idx < nK0Star; ++idx) {
    for (unsigned int iSel = 0; iSel < selectedElements.size(); ++iSel) {
      if(idx == selectedElements[iSel]) mask[idx] = false;
    }
  }

  return mask;

}

bool firedTrigger(Vec_i id, int v1, Vec_f pt, float v2){

  bool fired = false;
  for (unsigned int i = 0; i < id.size(); i++){
    if (id[i]==v1 && pt[i]>v2) fired = true;
  }
  return fired;
}


Vec_f getMinimum(Vec_f v1, Vec_f v2){
	Vec_f output = {};
	if (!v1.empty() && !v2.empty() && v1.size() == v2.size()){
		for (unsigned int i = 0; i < v1.size(); i++){
			output.push_back(min(v1[i], v2[i]));
		}
	}
	return output;
}

Vec_f getMaximum(Vec_f v1, Vec_f v2){
	Vec_f output = {};
	if (!v1.empty() && !v2.empty() && v1.size() == v2.size()){
		for (unsigned int i = 0; i < v1.size(); i++){
			output.push_back(max(v1[i], v2[i]));
		}
	}
	return output;
}

Vec_i indices(const int& size, const int& start = 0) {
    Vec_i res(size, 0);
    std::iota(std::begin(res), std::end(res), start);
    return res;
}

Vec_f NomUpDownVar(const float nom, const float up, const float down, float weight) {

  Vec_f res(3, 1);
  res[0] = weight;  // nom - already mutliplied for the Nom
  res[1] = (nom!=0) ? weight*up/nom : weight;  // up
  res[2] = (nom!=0) ? weight*down/nom : weight;  // down
  return res;
}



bool isGoodRunLS(const bool isData, const UInt_t run, const UInt_t lumi) {

  if(not isData) return true;

  if(jsonMap.find(run) == jsonMap.end()) return false; // run not found

  auto& validlumis = jsonMap.at(run);
  auto match = std::lower_bound(std::begin(validlumis), std::end(validlumis), lumi,
				[](std::pair<unsigned int, unsigned int>& range, unsigned int val) { return range.second < val; });
  return match->first <= lumi && match->second >= lumi;
}

float deltaPhi(float phi1, float phi2) {
  float result = phi1 - phi2;
  while (result > float(M_PI)) result -= float(2*M_PI);
  while (result <= -float(M_PI)) result += float(2*M_PI);
  return result;
}

float deltaR2(float eta1, float phi1, float eta2, float phi2) {
  float deta = eta1-eta2;
  float dphi = deltaPhi(phi1,phi2);
  return deta*deta + dphi*dphi;
}

float deltaR(float eta1, float phi1, float eta2, float phi2) {
  return std::sqrt(deltaR2(eta1,phi1,eta2,phi2));
}

Vec_b cleaningMask(Vec_i indices, int size) {

  Vec_b mask(size, true);
  for (int idx : indices) {
    if(idx < 0) continue;
    mask[idx] = false;
  }
  return mask;
}

Vec_b cleaningJetFromOBJ(Vec_f & Jeta, Vec_f & Jphi, float & eta, float & phi) {

  Vec_b mask(Jeta.size(), true);
  for (unsigned int idx = 0; idx < Jeta.size(); ++idx) {
    if(deltaR(Jeta[idx], Jphi[idx], eta, phi)<0.5) mask[idx] = false;
  }
  return mask;
}

//####

bool checkMother(const Vec_i& genPart_pdgId, Vec_i& genPart_genPartIdxMother,
		 int pdgToMatch,
		 int pdgMotherToMatch
		 ) {

  int idxMother = -1;
  // loop over all the genPartCand
  for (unsigned int i=0; i<genPart_pdgId.size(); i++) {
    if(genPart_pdgId[i]==25) idxMother = i;
  }

  int idxIntermediate = -1;
  if(pdgMotherToMatch==333 and (pdgToMatch==310 or pdgToMatch==130)) {
    for (unsigned int i=0; i<genPart_pdgId.size(); i++) {
      if(genPart_pdgId[i]==333 and genPart_genPartIdxMother[i]==idxMother) { idxIntermediate = i; }
    }
  }

  if(pdgMotherToMatch==423 and pdgToMatch==421) {
    for (unsigned int i=0; i<genPart_pdgId.size(); i++) {
      if(genPart_pdgId[i]==423 and genPart_genPartIdxMother[i]==idxMother) { idxIntermediate = i; }
    }
  }

  int idxMatch = -1.;
  // loop over all the genPartCand
  for (unsigned int i=0; i<genPart_pdgId.size(); i++) {
    if(abs(genPart_pdgId[i])==pdgToMatch) {
      if(genPart_pdgId[i] == pdgToMatch && genPart_genPartIdxMother[i] == idxIntermediate ) {
	idxMatch=i;
      }
    }
  }

  bool foundMother = (idxMatch!=-1);
  return foundMother;

}

Vec_i genMatchRECO(const Vec_f& reco_pt, const Vec_f& reco_eta, const Vec_f& reco_phi, const Vec_f& reco_mass,
		   const Vec_f& genPart_eta, const Vec_f& genPart_phi,
		   const Vec_i& genPart_pdgId, Vec_i& genPart_genPartIdxMother,
		   int pdgToMatch,
		   int pdgMotherToMatch
		   ) {

  int idxMother = -1;
  // loop over all the genPartCand
  for (unsigned int i=0; i<genPart_pdgId.size(); i++) {
    if(genPart_pdgId[i]==25) idxMother = i;
  }

  int idxIntermediate = -1;
  if(pdgMotherToMatch==333 and pdgToMatch==310) {
    for (unsigned int i=0; i<genPart_pdgId.size(); i++) {
      if(genPart_pdgId[i]==333 and genPart_genPartIdxMother[i]==idxMother) { idxIntermediate = i; }
    }
  }

  if(pdgMotherToMatch==423 and pdgToMatch==421) {
    for (unsigned int i=0; i<genPart_pdgId.size(); i++) {
      if(genPart_pdgId[i]==423 and genPart_genPartIdxMother[i]==idxMother) { idxIntermediate = i; }
    }
  }

  int idxMatch = -1;
  float etaGen = 999.;
  float phiGen = 999.;
  // loop over all the genPartCand
  for (unsigned int i=0; i<genPart_pdgId.size(); i++) {
    // 310 is K0s (fromPhi) 421 is D0 (from D0Star)
    if((pdgToMatch==310 or pdgToMatch==421) and abs(genPart_pdgId[i])==pdgToMatch) {
      if(abs(genPart_pdgId[i]) == pdgToMatch && genPart_genPartIdxMother[i] == idxIntermediate ) {
	idxMatch=i;  etaGen=genPart_eta[i]; phiGen=genPart_phi[i];
      }
    } else if(genPart_pdgId[i]==22 or abs(genPart_pdgId[i])==333 or abs(genPart_pdgId[i])==113 or abs(genPart_pdgId[i])==223 or abs(genPart_pdgId[i])==443 or abs(genPart_pdgId[i])==313) {
      if(abs(genPart_pdgId[i]) == pdgToMatch && genPart_genPartIdxMother[i] == idxMother ) {
	idxMatch=i;  etaGen=genPart_eta[i]; phiGen=genPart_phi[i];
      }
    }
  }

  Vec_i idx(2, -1); // initialize with -1 a vector of size 2
  // idx[0] for the RECO index 
  // idx[1] for the GEN index
  
  if (pdgToMatch==421) idx[1] = idxIntermediate; // save D0Star instead of D0(no real mass)
  else idx[1] = idxMatch; // index of the genPart (-1 if not found gen cand)

  // loop over all the recoCand
  for (unsigned int i=0; i<reco_pt.size(); i++) {
    PtEtaPhiMVector p_reco(reco_pt[i], reco_eta[i], reco_phi[i], reco_mass[i]);
    if(idxMother<0) continue; // no good Mother (i.e. Higgs)
    if(idxMatch<0) continue; // no good Match
    if(deltaR(etaGen,phiGen,p_reco.eta(),p_reco.phi())>0.01) continue;
    idx[0] = i;
  }

  return idx;

}

Vec_i genMatch(const float& reco_pt, const float& reco_eta, const float& reco_phi, const float& reco_mass,
	       const Vec_f& genPart_eta, const Vec_f& genPart_phi,
	       const Vec_i& genPart_pdgId, Vec_i& genPart_genPartIdxMother,
	       int pdgToMatch
	       //		   int pdgMotherToMatch
	       ) {

  // examine the recoCand
  PtEtaPhiMVector p_reco(reco_pt, reco_eta, reco_phi, reco_mass);
  Vec_i idx(1, -1); // initialize with -1 a vector of size 1
  // idx[1] for the GEN index
  
  int idxMatch = -1;
  float etaGen = 999.;
  float phiGen = 999.;
  // loop over all the genPartCand
  for (unsigned int i=0; i<genPart_pdgId.size(); i++) {
    // 310 is K0s (fromPhi) 421 is D0 (from D0Star)
    //    if(genPart_pdgId[i]==22 or abs(genPart_pdgId[i])==333 or abs(genPart_pdgId[i])==113 or abs(genPart_pdgId[i])==223 or abs(genPart_pdgId[i])==443 or abs(genPart_pdgId[i])==313) {
    if(abs(genPart_pdgId[i])==333 or abs(genPart_pdgId[i])==113 or abs(genPart_pdgId[i])==223 or abs(genPart_pdgId[i])==443 or abs(genPart_pdgId[i])==313) {      
      if(abs(genPart_pdgId[i]) == pdgToMatch ) {
	idxMatch=i;  etaGen=genPart_eta[i]; phiGen=genPart_phi[i];
	if(deltaR(etaGen,phiGen,p_reco.eta(),p_reco.phi())<0.05) {    
	  idx[0] = idxMatch; // index of the genPart (-1 if not found gen cand)
	}
      }
    }
  }

  return idx;

}

stdVec_i HiggsCandFromRECO(const Vec_f& meson_pt, const Vec_f& meson_eta, const Vec_f& meson_phi, const Vec_f& meson_mass,
			   const Vec_f& meson_trk1_pt, const Vec_f& meson_trk2_pt,
			   const Vec_f& wrong_meson_pt, const Vec_f& wrong_meson2_pt,
			   const Vec_f& ph_pt, const Vec_f& ph_eta, const Vec_f& ph_phi) {

  float Minv = -1;
  float ptHiggs = -1;
  float ptCandMax=0;
  PtEtaPhiMVector p_ph(ph_pt[0], ph_eta[0], ph_phi[0], 0);
  unsigned int indexPhoton = 0;
  stdVec_i idx(2, -1); // initialize with -1 a vector of size 2

  if(ph_pt.size()> 1) {
    if(ph_pt[1] > ph_pt[0]) p_ph.SetPt(ph_pt[1]);
    if(ph_pt[1] > ph_pt[0]) p_ph.SetEta(ph_eta[1]);
    if(ph_pt[1] > ph_pt[0]) p_ph.SetPhi(ph_phi[1]);
    indexPhoton = 1;
  }

  float ptWrongMax=0;
  for (unsigned int j=0; j<wrong_meson_pt.size(); j++) {
    if(wrong_meson_pt[j] <  ptWrongMax) continue;
    ptWrongMax = wrong_meson_pt[j];
  }

  float ptWrong2Max=0;
  for (unsigned int j=0; j<wrong_meson2_pt.size(); j++) {
    if(wrong_meson2_pt[j] <  ptWrong2Max) continue;
    ptWrong2Max = wrong_meson2_pt[j];
  }

  // loop over all the phi/rho Cand
  for (unsigned int i=0; i<meson_pt.size(); i++) {

    // if(max(meson_trk1_pt[i], meson_trk2_pt[i]) < 20) continue;

    PtEtaPhiMVector p_meson(meson_pt[i], meson_eta[i], meson_phi[i], meson_mass[i]);
    if((p_meson + p_ph).M()<5.) continue; // object disambiguation, also remove the omega/tau resonances

    // save the leading Pt
    float ptCand = p_meson.pt();

    if( ptCand < ptCandMax ) continue;
    ptCandMax=ptCand;
    if(ptCandMax < ptWrongMax) continue; // we want the leading meson to the of the right flavor
    if(ptCandMax < ptWrong2Max) continue; // we want the leading meson to the of the right flavor
    Minv = (p_meson + p_ph).M();
    ptHiggs = (p_meson + p_ph).pt();
    idx[0] = i;
    idx[1] = indexPhoton;
  }

  return idx;
}

int topology(float eta1, float eta2) {

  int topology = 0;
  if (abs(eta1)<1.4 and abs(eta2)<1.4) topology = 1; // this is BB
  if (abs(eta1)<1.4 and abs(eta2)>1.4) topology = 2; // this is BE
  if (abs(eta1)>1.4 and abs(eta2)<1.4) topology = 3; // this is EB
  if (abs(eta1)>1.4 and abs(eta2)>1.4) topology = 4; // this is EE

  return topology;
}

Vec_i mesonCand(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m, const Vec_f& ch,
		const Vec_f& ph_pt, const Vec_f& ph_eta, const Vec_f& ph_phi,
		bool phiHyp
		) {
  
  Vec_i idx(2, -1); // initialize with -1 a vector of size 2

  float minPtTracks_= phiHyp ? 10: 0.f;
  float minPtMeson_= phiHyp ? 10: 20.f;

  float phiLeadingPhoton = ph_phi[0];
  if(ph_pt[1] > ph_pt[0]) phiLeadingPhoton = ph_phi[1];

  float mass=-100;
  float ptCandMax=0;
  for (unsigned int i=0; i<pt.size(); i++) {

    for (unsigned int j=0; j<i; j++) {
      
      if(i==j) continue;
      if(ch[i]*ch[j]>0) continue; // opposite sign kaons

      if(pt[i]<minPtTracks_ and pt[j]<minPtTracks_) continue; //at least a kaons of X GeV

      PtEtaPhiMVector p1(pt[i], eta[i], phi[i], m[i]);
      PtEtaPhiMVector p2(pt[j], eta[j], phi[j], m[j]);

      if(abs(deltaPhi(phiLeadingPhoton, (p1 + p2).phi()))<float(M_PI/2)) continue; // M,gamma opposite hemishpere

      if(deltaR(eta[i], phi[i], eta[j], phi[j])>0.5) continue; // meson's decay product inside a narrow cone

      float ptCand = (p1 + p2).pt();
      if( ptCand < minPtMeson_ ) continue; // pt of the Cand
      if( ptCand < ptCandMax ) continue; // leading in pt candidate in case of multiple combination

      ptCandMax=ptCand;

      if (idx[0] == -1) { idx[0] = i; idx[1] = j; }

    }
  }

  return idx;

}

bool hasTriggerMatch(const float& eta, const float& phi, const Vec_f& TrigObj_eta, const Vec_f& TrigObj_phi) {

  for (unsigned int jtrig = 0; jtrig < TrigObj_eta.size(); ++jtrig) {
    if (deltaR(eta, phi, TrigObj_eta[jtrig], TrigObj_phi[jtrig]) < 0.3) return true;
  }
  return false;
}

float mt(float pt1, float phi1, float pt2, float phi2) {
  return std::sqrt(2*pt1*pt2*(1-std::cos(phi1-phi2)));
}

Vec_f sum2BodyMass(const Vec_f& pt1, const Vec_f& eta1, const Vec_f& phi1, const Vec_f& m1,
		   const Vec_f& pt2, const Vec_f& eta2, const Vec_f& phi2) {

  const float pion_mass_ = 0.139570;
  Vec_f mass(pt1.size());

  for (unsigned int idx = 0; idx < pt1.size(); ++idx) {
    PtEtaPhiMVector p_1(pt1[idx], eta1[idx], phi1[idx], m1[idx]);
    PtEtaPhiMVector p_2(pt2[idx], eta2[idx], phi2[idx], pion_mass_);
    mass[idx] = (p_1 + p_2).M();
  }
  return mass;

}

float Minv(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m) {
  PtEtaPhiMVector p1(pt[0], eta[0], phi[0], m[0]);
  PtEtaPhiMVector p2(pt[1], eta[1], phi[1], m[1]);
  return (p1 + p2).M();
}

Vec_f Pair12PT(const Vec_f& pt1, const Vec_f& eta1, const Vec_f& phi1,
	       const Vec_f& pt2, const Vec_f& eta2, const Vec_f& phi2) {

  const float pion_mass_     = 0.139570;
  Vec_f pair(pt1.size(), true);

  for (unsigned int idx = 0; idx < pt1.size(); ++idx) {
    PtEtaPhiMVector p_1(pt1[idx], eta1[idx], phi1[idx], pion_mass_);
    PtEtaPhiMVector p_2(pt2[idx], eta2[idx], phi2[idx], pion_mass_);
    pair[idx] = (p_1 + p_2).pt();
  }
  return pair;
}

Vec_f Pair12ETA(const Vec_f& pt1, const Vec_f& eta1, const Vec_f& phi1,
		const Vec_f& pt2, const Vec_f& eta2, const Vec_f& phi2) {

  const float pion_mass_     = 0.139570;
  Vec_f pair(pt1.size(), true);

  for (unsigned int idx = 0; idx < pt1.size(); ++idx) {
    PtEtaPhiMVector p_1(pt1[idx], eta1[idx], phi1[idx], pion_mass_);
    PtEtaPhiMVector p_2(pt2[idx], eta2[idx], phi2[idx], pion_mass_);
    pair[idx] = (p_1 + p_2).eta();
  }
  return pair;
}

Vec_f Pair12PHI(const Vec_f& pt1, const Vec_f& eta1, const Vec_f& phi1,
		const Vec_f& pt2, const Vec_f& eta2, const Vec_f& phi2) {

  const float pion_mass_     = 0.139570;
  Vec_f pair(pt1.size(), true);

  for (unsigned int idx = 0; idx < pt1.size(); ++idx) {
    PtEtaPhiMVector p_1(pt1[idx], eta1[idx], phi1[idx], pion_mass_);
    PtEtaPhiMVector p_2(pt2[idx], eta2[idx], phi2[idx], pion_mass_);
    pair[idx] = (p_1 + p_2).phi();
  }
  return pair;
}

float Minv3(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m,
	    const float& ph_pt, const float& ph_eta, const float& ph_phi) {
 
  PtEtaPhiMVector p1(pt[0], eta[0], phi[0], m[0]);
  PtEtaPhiMVector p2(pt[1], eta[1], phi[1], m[1]);
  PtEtaPhiMVector p_ph(ph_pt, ph_eta, ph_phi, 0);
  return (p1 + p2 + p_ph).M();
  
}

float PT(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m, Vec_i idx) {
  PtEtaPhiMVector p1(pt[idx[0]], eta[idx[0]], phi[idx[0]], m[idx[0]]);
  PtEtaPhiMVector p2(pt[idx[1]], eta[idx[1]], phi[idx[1]], m[idx[1]]);
  return (p1 + p2).pt();
}

float dPhi_MvsPh(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m, Vec_i idx,
		 const Vec_f& ph_pt, const Vec_f& ph_eta, const Vec_f& ph_phi) {

  PtEtaPhiMVector p1(pt[idx[0]], eta[idx[0]], phi[idx[0]], m[idx[0]]);
  PtEtaPhiMVector p2(pt[idx[1]], eta[idx[1]], phi[idx[1]], m[idx[1]]);
  
  return  abs(deltaPhi(ph_phi[0], (p1 + p2).phi()));
}

float dEta_MvsPh(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m, Vec_i idx,
		 const Vec_f& ph_eta) {

  PtEtaPhiMVector p1(pt[idx[0]], eta[idx[0]], phi[idx[0]], m[idx[0]]);
  PtEtaPhiMVector p2(pt[idx[1]], eta[idx[1]], phi[idx[1]], m[idx[1]]);
  
  return  ph_eta[0]-(p1 + p2).eta();
}

std::pair<float, float>  Minv2(const float& pt, const float& eta, const float& phi, const float& m,
                               const float& ph_pt, const float& ph_eta, const float& ph_phi) {

  PtEtaPhiMVector p_M(pt, eta, phi, m);
  PtEtaPhiMVector p_ph(ph_pt, ph_eta, ph_phi, 0);

  float Minv = (p_M + p_ph).M();
  float ptPair = (p_M + p_ph).pt();

  std::pair<float, float> pairRECO = std::make_pair(Minv , ptPair);
  return pairRECO;

}

float dR_Constituents(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m, Vec_i idx) {
  return deltaR(eta[idx[0]], phi[idx[0]], eta[idx[1]], phi[idx[1]]);
}

float compute_HiggsVars_var_VtxCorr(const float mes_pt, const float mes_eta, const float mes_phi, const float mes_mass,
				    const float mes_vtx_X, const float mes_vtx_Y, const float mes_vtx_Z,
				    const float ph_pt, const float ph_eta, const float ph_phi,
				    const float ph_calo_X, const float ph_calo_Y, const float ph_calo_Z,
				    unsigned int var)
{

  // passing only the one that make the Higgs candidate
  // Here I have to do the correction to the photon 4 momentum
  PtEtaPhiMVector p4_ph(ph_pt, ph_eta, ph_phi, 0);

  XYZVectorF calorimiterPos(ph_calo_X, ph_calo_Y, ph_calo_Z);
  XYZVectorF mesonVertex(mes_vtx_X, mes_vtx_Y, mes_vtx_Z);

  XYZVectorF p_ph_origin = calorimiterPos.unit() * p4_ph.E();
  XYZVectorF p_ph_vtx = (calorimiterPos - mesonVertex).unit() * p4_ph.E();

  PtEtaPhiMVector p4_ph_origin;
  p4_ph_origin.SetPxPyPzE(p_ph_origin.X(), p_ph_origin.Y(), p_ph_origin.Z(), p4_ph.E());

  PtEtaPhiMVector p4_ph_vtx;
  p4_ph_vtx.SetPxPyPzE(p_ph_vtx.X(), p_ph_vtx.Y(), p_ph_vtx.Z(), p4_ph.E());

  //Correct them so they have  m=0
  PtEtaPhiMVector p4_ph_origin_M(p4_ph_origin.Pt(), p4_ph_origin.Eta(), p4_ph_origin.Phi(), 0);
  PtEtaPhiMVector p4_ph_vtx_M(p4_ph_vtx.Pt(), p4_ph_vtx.Eta(), p4_ph_vtx.Phi(), 0);

  PtEtaPhiMVector p_mes(mes_pt, mes_eta, mes_phi, mes_mass);
  PtEtaPhiMVector p_Hig_ORG = (p4_ph + p_mes);

  //  cout << "--------------------------------------------------------------------------------------------" << endl;
  //  cout << "(XYZ) Calo coords [mm]:                     " << calorimiterPos.X()  << " " << calorimiterPos.Y()    << " " << calorimiterPos.Z()    << endl;
  //  cout << "(XYZ) Meson vtx coords [mm]:                " << mesonVertex.X()     << " " << mesonVertex.Y()       << " " << mesonVertex.Z()       << endl;
  //  cout << "--------------------------------------------" << endl;
  //  cout << "(PxPyPzE|p|) (0-calo) photon [GeV]:         " << p4_ph_origin.X()    << " " << p4_ph_origin.Y()      << " " << p4_ph_origin.Z()      << " " << p4_ph_origin.E()    << " " << p4_ph_origin.P()    << endl;
  //  cout << "(PxPyPzE|p|) (meson_vtx-calo) photon [GeV]: " << p4_ph_vtx.X()       << " " << p4_ph_vtx.Y()         << " " << p4_ph_vtx.Z()         << " " << p4_ph_vtx.E()       << " " << p4_ph_vtx.P()       << endl;
  //  cout << "(PxPyPzE|p|) Original photon [GeV]:         " << p4_ph.X()           << " " << p4_ph.Y()             << " " << p4_ph.Z()             << " " << p4_ph.E()           << " " << p4_ph.P()           << endl;
  //  cout << "--------------------------------------------" << endl;
  //  cout << "(PtEtaPhi) (0-calo) photon manually comp:   " << std::sqrt(p_ph_origin.X()*p_ph_origin.X()+p_ph_origin.Y()*p_ph_origin.Y())     << " " << std::atanh(p_ph_origin.Z()/p_ph_origin.P())    << " " << std::atan2(p_ph_origin.Y(), p_ph_origin.X()) << endl;
  //  cout << "(PtEtaPhiME) (0-calo) photon:               " << p4_ph_origin.Pt()   << " " << p4_ph_origin.Eta()    << " " << p4_ph_origin.Phi()    << " " << p4_ph_origin.M()    << " " << p4_ph_origin.E()    << endl;
  //  cout << "(PtEtaPhiME) (0-calo) photon M=0:           " << p4_ph_origin_M.Pt() << " " << p4_ph_origin_M.Eta()  << " " << p4_ph_origin_M.Phi()  << " " << p4_ph_origin_M.M()  << " " << p4_ph_origin_M.E()  << endl;
  //  cout << "(PtEtaPhiME) (meson_vtx-calo) photon:       " << p4_ph_vtx.Pt()      << " " << p4_ph_vtx.Eta()       << " " << p4_ph_vtx.Phi()       << " " << p4_ph_vtx.M()       << " " << p4_ph_vtx.E()       << endl;
  //  cout << "(PtEtaPhiME) (meson_vtx-calo) photon M=0:   " << p4_ph_vtx_M.Pt()    << " " << p4_ph_vtx_M.Eta()     << " " << p4_ph_vtx_M.Phi()     << " " << p4_ph_vtx_M.M()     << " " << p4_ph_vtx_M.E()     << endl;
  //  cout << "(PtEtaPhiME) Original photon [GeV]:             " << p4_ph.Pt()          << " " << p4_ph.Eta()           << " " << p4_ph.Phi()           << " " << p4_ph.M()           << " " << p4_ph.E()           << endl;
  //  cout << "--------------------------------------------------------------------------------------------" << endl;

  PtEtaPhiMVector p_Hig = (p4_ph_vtx_M + p_mes);
  float theVar = 0;
  if     (var == 0) theVar = p_Hig.M();
  else if(var == 1) theVar = p_Hig.Pt();
  else if(var == 2) theVar = p_Hig.Phi();

  return theVar;

}


float compute_HiggsVars_var(const float mes_pt, const float mes_eta, const float mes_phi, const float mes_mass,
			    const float ph_pt, const float ph_eta, const float ph_phi,
			    unsigned int var)
{

  // passing only the one that make the Higgs candidate
  PtEtaPhiMVector p_ph(ph_pt, ph_eta, ph_phi, 0);
  PtEtaPhiMVector p_mes(mes_pt, mes_eta, mes_phi, mes_mass);

  PtEtaPhiMVector p_Hig = (p_ph + p_mes);
  float theVar = 0;
  if     (var == 0) theVar = p_Hig.M();
  else if(var == 1) theVar = p_Hig.Pt();
  else if(var == 2) theVar = p_Hig.Phi();

  return theVar;

}

float compute_jet_HiggsVars_var(const Vec_f& jet_pt, const Vec_f& jet_eta, const Vec_f& jet_phi, const Vec_f& jet_mass, 
				const float ph_pt, const float ph_eta, const float ph_phi,
				const float mes_pt, const float mes_eta, const float mes_phi, const float mes_mass,
				unsigned int var)
{

  if(jet_pt.size() < 2) return -1;
  // passing only the one that make the Higgs candidate
  PtEtaPhiMVector p_ph(ph_pt, ph_eta, ph_phi, 0);
  PtEtaPhiMVector p_mes(mes_pt, mes_eta, mes_phi, mes_mass);

  PtEtaPhiMVector p_Hig = (p_ph + p_mes);
  
  PtEtaPhiMVector p_j1(jet_pt[0], jet_eta[0], jet_phi[0], jet_mass[0]);
  PtEtaPhiMVector p_j2(jet_pt[1], jet_eta[1], jet_phi[1], jet_mass[1]);

  if(p_j1.Pt() < p_j2.Pt()) printf("Pt jet reversed!\n");

  float deltaEtaJJ = fabs(p_j1.Eta()-p_j2.Eta());

  float theVar = 0;
  if     (var == 0) theVar = fabs(p_Hig.Eta()-(p_j1.Eta()+p_j2.Eta())/2.)/deltaEtaJJ; //zeppenfeld variable
  else if(var == 1) theVar = fabs(p_Hig.Eta()-p_j1.Eta());
  else if(var == 2) theVar = fabs(p_Hig.Eta()-p_j2.Eta());
  else if(var == 3) theVar = (p_Hig + p_j1 + p_j2).Pt();

  return theVar;

}

// functions for polarization reweighting

// Get four vectors
Vec_d getHPhiKaonFourVectors(const Vec_i& genPart_pdgId, Vec_i& genPart_genPartIdxMother,
        const Vec_d& genPart_pt,
        const Vec_d& genPart_eta,
        const Vec_d& genPart_phi,
        const Vec_d& genPart_mass) {

  int idxHiggsfromKPlus = -1;
  int idxHiggsfromKMinus = -1;
  int idxPhifromKPlus = -1;
  int idxPhifromKMinus = -1;
  int idxKPlus = -1;
  int idxKMinus = -1;
  for (int i=genPart_pdgId.size(); i>=0; i--) {
    if (genPart_pdgId[i]==321 &&
        genPart_pdgId[genPart_genPartIdxMother[i]]==333 &&
        genPart_pdgId[genPart_genPartIdxMother[genPart_genPartIdxMother[i]]] == 25) {
          idxKPlus = i;
          idxPhifromKPlus = genPart_genPartIdxMother[i];
          idxHiggsfromKPlus = genPart_genPartIdxMother[genPart_genPartIdxMother[i]];
          break;
        }
  }
  for (int i=genPart_pdgId.size(); i>=0; i--) {
    if (genPart_pdgId[i]==-321 &&
        genPart_pdgId[genPart_genPartIdxMother[i]]==333 &&
        genPart_pdgId[genPart_genPartIdxMother[genPart_genPartIdxMother[i]]] == 25) {
          idxKMinus = i;
          idxPhifromKMinus = genPart_genPartIdxMother[i];
          idxHiggsfromKMinus = genPart_genPartIdxMother[genPart_genPartIdxMother[i]];
          break;
        }
  }

  bool foundMother = (idxKPlus!=-1 && idxKMinus!=-1 &&
                      idxPhifromKPlus!=-1 && idxPhifromKPlus==idxPhifromKMinus &&
                      idxHiggsfromKPlus!=-1 && idxHiggsfromKPlus==idxHiggsfromKMinus);

  Vec_f res(8, 0.);

  if (foundMother) {
    int idxPhi = idxPhifromKPlus;
    PtEtaPhiMVector p_kplus(genPart_pt[idxKPlus], genPart_eta[idxKPlus],
                            genPart_phi[idxKPlus], genPart_mass[idxKPlus]);
    PtEtaPhiMVector p_kminus(genPart_pt[idxKMinus], genPart_eta[idxKMinus],
                             genPart_phi[idxKMinus], genPart_mass[idxKMinus]);
    PtEtaPhiMVector p_ditrack = p_kplus + p_kminus;

    res[0] = p_ditrack.Pt();
    res[1] = p_ditrack.Eta();
    res[2] = p_ditrack.Phi();
    res[3] = p_ditrack.M();

    PtEtaPhiMVector p_phi(genPart_pt[idxPhi], genPart_eta[idxPhi],
                          genPart_phi[idxPhi], genPart_mass[idxPhi]);

    res[4] = p_phi.Pt();
    res[5] = p_phi.Eta();
    res[6] = p_phi.Phi();
    res[7] = p_phi.M();
  }
  return res;
}

// Checks if H->Phi->Kaon
bool isHPhiKaon(
  const Vec_i& genPart_pdgId,
  const Vec_i& genPart_genPartIdxMother
) {
  int idxHiggs = -1;
  int idxPhi = -1;
  int idxKPlus = -1;
  for (int i=genPart_pdgId.size(); i>=0; i--) {
    if (genPart_pdgId[i]==321 &&
        genPart_pdgId[genPart_genPartIdxMother[i]]==333 &&
        genPart_pdgId[genPart_genPartIdxMother[genPart_genPartIdxMother[i]]]==25) {
      idxKPlus = i;
      idxPhi = genPart_genPartIdxMother[i];
      idxHiggs = genPart_genPartIdxMother[genPart_genPartIdxMother[i]];
      break;
    }
  }

  bool foundMother = (idxHiggs!=-1 && idxPhi!=-1 && idxKPlus!=-1);

  return foundMother;
}

// Checks if H->Phi->Phi->Kaon
bool isHPhiPhiKaon(
  const Vec_i& genPart_pdgId,
  const Vec_i& genPart_genPartIdxMother
) {
  int idxHiggs = -1;
  int idxPhi = -1;
  int idxKPlus = -1;
  for (int i=genPart_pdgId.size(); i>=0; i--) {
    if (genPart_pdgId[i]==321 &&
        genPart_pdgId[genPart_genPartIdxMother[i]]==333 &&
        genPart_pdgId[genPart_genPartIdxMother[genPart_genPartIdxMother[i]]]==333 &&
        genPart_pdgId[genPart_genPartIdxMother[genPart_genPartIdxMother[genPart_genPartIdxMother[i]]]]==25) {
      idxKPlus = i;
      idxPhi = genPart_genPartIdxMother[i];
      idxHiggs = genPart_genPartIdxMother[genPart_genPartIdxMother[genPart_genPartIdxMother[i]]];
      break;
    }
  }

  bool foundMother = (idxHiggs!=-1 && idxPhi!=-1 && idxKPlus!=-1);

  return foundMother;
}

float getPhiPolarizationAngle(
  const Vec_i& genPart_pdgId,
  const Vec_i& genPart_genPartIdxMother,
  const Vec_d& genPart_pt,
  const Vec_d& genPart_eta,
  const Vec_d& genPart_phi,
  const Vec_d& genPart_mass) {

  ///// Part 1. Find the K+K- tracks from the phi
  int idxHiggs = -1;
  int idxPhi = -1;
  int idxKPlus = -1;
  for (int i=genPart_pdgId.size(); i>=0; i--) {
    if (genPart_pdgId[i]==321 &&
        genPart_pdgId[genPart_genPartIdxMother[i]]==333 &&
        genPart_pdgId[genPart_genPartIdxMother[genPart_genPartIdxMother[i]]] == 25) {
      idxKPlus = i;
      idxPhi = genPart_genPartIdxMother[i];
      idxHiggs = genPart_genPartIdxMother[genPart_genPartIdxMother[i]];
      break;
    } /*else if (genPart_pdgId[i]==321 &&
               genPart_pdgId[genPart_genPartIdxMother[i]]==333 &&
               genPart_pdgId[genPart_genPartIdxMother[genPart_genPartIdxMother[i]]]==333 &&
               genPart_pdgId[genPart_genPartIdxMother[genPart_genPartIdxMother[genPart_genPartIdxMother[i]]]]==25) {
      idxKPlus = i;
      idxPhi = genPart_genPartIdxMother[i];
      idxHiggs = genPart_genPartIdxMother[genPart_genPartIdxMother[genPart_genPartIdxMother[i]]];
      break;
    }*/
  }

  bool foundMother = (idxHiggs!=-1 && idxPhi!=-1 && idxKPlus!=-1);

  ///// Part 2. Get polarization angle
  if (foundMother) {
    TLorentzVector p_phi;
    TLorentzVector p_kplus;
    p_phi.SetPtEtaPhiM(genPart_pt[idxPhi], genPart_eta[idxPhi],
                       genPart_phi[idxPhi], genPart_mass[idxPhi]);
    p_kplus.SetPtEtaPhiM(genPart_pt[idxKPlus], genPart_eta[idxKPlus],
                         genPart_phi[idxKPlus], genPart_mass[idxKPlus]);
    TVector3 phiBoost = p_phi.BoostVector();
    p_kplus.Boost(-phiBoost);
    float theta = p_phi.Vect().Angle(p_kplus.Vect());
    return theta;
  } else {
    return std::numeric_limits<float>::quiet_NaN();
  }
}

TH1F * myTH1_photon;
TH1F * myTH1_twoProng;

float SF_HLT_leg(float pt, int idParticle, int nomUpDowm) {

  TH1* h_sf;

  if(idParticle == 22) h_sf = myTH1_photon;
  if(idParticle == 15) h_sf = myTH1_twoProng;

  float trgSF = 1.0;
  int ptbin = std::max(1, std::min(h_sf->GetNbinsX(), h_sf->GetXaxis()->FindBin(pt)));

  if(nomUpDowm==0) trgSF = h_sf->GetBinContent(ptbin);
  if(nomUpDowm==1) trgSF = h_sf->GetBinContent(ptbin) + std::min(0.10, h_sf->GetBinError(ptbin));
  if(nomUpDowm==-1) trgSF = h_sf->GetBinContent(ptbin) - std::min(0.10, h_sf->GetBinError(ptbin));

  return trgSF;

}

void initTrigSF() {

  const char* filename_ph = "/work/submit/mariadlf/Hrare/utilFiles/sfFiles/Photon35TriggerSF_FSR.root";
  //  const char* filename_ph = "/work/submit/mariadlf/Hrare/utilFiles/sfFiles/Photon35TriggerSF.root";
  const char* filename_had = "/work/submit/mariadlf/Hrare/utilFiles/sfFiles/TwoProngsTriggerSF.root";

  TFile* myFile_photon = TFile::Open(filename_ph, "READ");
  myTH1_photon = static_cast<TH1F*>(myFile_photon->Get("h_triggerEff_eT"));

  TFile* myFile_twoProng = TFile::Open(filename_had, "READ");
  myTH1_twoProng = static_cast<TH1F*>(myFile_twoProng->Get("h_efficiency_Data"));

}

TH1F * myTH1_chIso_barrel;
TH1F * myTH1_chIso_endcap;

float SF_chIso(float pt, float eta, int nomUpDowm) {

  TH1* h_sf;

  if(abs(eta)<1.4) h_sf = myTH1_chIso_barrel;
  if(abs(eta)>=1.4) h_sf = myTH1_chIso_endcap;

  float isoSF = 1.0;
  int ptbin = std::max(1, std::min(h_sf->GetNbinsX(), h_sf->GetXaxis()->FindBin(pt)));

  if(nomUpDowm==0) isoSF = h_sf->GetBinContent(ptbin);
  if(nomUpDowm==1) isoSF = h_sf->GetBinContent(ptbin) + std::min(0.10, h_sf->GetBinError(ptbin));
  if(nomUpDowm==-1) isoSF = h_sf->GetBinContent(ptbin) - std::min(0.10, h_sf->GetBinError(ptbin));

  return isoSF;

}

void initIsoSF() {

  const char* filename_barrel = "/work/submit/mariadlf/Hrare/utilFiles/sfFiles/IsoChEfficiencySF_barrel.root";
  const char* filename_endcap = "/work/submit/mariadlf/Hrare/utilFiles/sfFiles/IsoChEfficiencySF_endcap.root";

  TFile* myFile_chIso_barrel = TFile::Open(filename_barrel, "READ");
  myTH1_chIso_barrel = static_cast<TH1F*>(myFile_chIso_barrel->Get("h_efficiency_Data"));

  TFile* myFile_chIso_endcap = TFile::Open(filename_endcap, "READ");
  myTH1_chIso_endcap = static_cast<TH1F*>(myFile_chIso_endcap->Get("h_efficiency_Data"));

}

std::vector<std::vector<float>> theta_pol_dict_VEC;

void initPol(int mc, int year, int nSlot) {

  TString prod = "";
  if(mc==1027 or mc==1017 or mc==1037) prod = "ggH";
  if(mc==1020 or mc==1010 or mc==1030) prod = "VBF";
  if(mc==1022 or mc==1012 or mc==1032) prod = "WminusH_WToLNu";
  if(mc==1021 or mc==1011 or mc==1031) prod = "WplusH_WToLNu";
  if(mc==1023 or mc==1013 or mc==1033) prod = "ZH_ZToLL";
  if(mc==1024 or mc==1014 or mc==1034) prod = "ggZH_ZToLL";
  if(mc==1028 or mc==1018 or mc==1038) prod = "TTH";

  TString meson = "";
  if(mc>=1010 and mc<=1018) meson = "Phi";
  if(mc>=1020 and mc<=1028) meson = "Rho";
  if(mc>=1030 and mc<=1038) meson = "K0s";

  const char* filename = Form("/work/submit/mariadlf/Hrare/utilFiles/polFiles/H%sGammaAnalysis_Signal_%s_%s_%d_M125_AOD.root",meson.Data(),meson.Data(),prod.Data(),year);
  const char* dirname = Form("H%sGammaAOD",meson.Data());
  std::cout << "reading polarization file: " << filename << std::endl;

  TFile* myFile = TFile::Open(filename, "READ");
  TDirectory* myDir = static_cast<TDirectory*>(myFile->Get(dirname));
  TTree* myTree = static_cast<TTree*>(myDir->Get("mytree"));
  myTree->SetBranchStatus("*",0);
  myTree->SetBranchStatus("event_number",1);
  myTree->SetBranchStatus("theta_pol",1);
  myTree->SetBranchStatus("genPhoton_eT",1);
  myTree->SetBranchStatus("genMeson_pT",1);

  std::vector<float> theta_pol_dict(3000000, 0.); // initialize with so that sin of 0 is 1.

  Int_t event_var;
  float theta_pol_var;
  float genPhoton_eT;
  float genMeson_pT;

  myTree->SetBranchAddress("event_number", &event_var);
  myTree->SetBranchAddress("theta_pol", &theta_pol_var);

  myTree->SetBranchAddress("genPhoton_eT", &genPhoton_eT);
  myTree->SetBranchAddress("genMeson_pT", &genMeson_pT);

  int countNan=0;
  Long64_t nentries = myTree->GetEntries();
  for (Long64_t i=0;i<nentries;++i) {
    myTree->GetEntry(i);
    theta_pol_dict.at(int(event_var)) = theta_pol_var;
    if (genMeson_pT < -9 and genPhoton_eT < -9) countNan++; 
  }

  std::cout << " countNan " << countNan << " out of tot events" << nentries << std::endl;

  for ( int i=0; i<nSlot;i++) {
    theta_pol_dict_VEC.push_back(theta_pol_dict);
  }
  std::cout << " theta_pol_dict_VEC.size() " << theta_pol_dict_VEC.size() << std::endl;

  myTree->Delete();
  myFile->Close();
  delete myFile;

}

float getPolAngle(const ULong64_t event_, int nSlot) {

  auto theta_pol_dict = theta_pol_dict_VEC[nSlot];
  return theta_pol_dict[event_];

}

stdVec_i jetCloseFar(const Vec_f& jet_pt, const Vec_f& jet_eta, const Vec_f& jet_phi, const Vec_f& jet_m,
                     const Vec_f& jpsi_pt, const Vec_f& jpsi_eta, const Vec_f& jpsi_phi, const Vec_f& jpsi_m) {

  unsigned int indexClose = -1;

  for (unsigned int idx = 0; idx < jet_pt.size(); ++idx) {
    float drJetPhi = deltaR(jet_eta[idx], jet_phi[idx], jpsi_eta[0], jpsi_phi[0]);
    if (drJetPhi < 0.4 and indexClose==-1) indexClose=idx;
  }

  stdVec_i idx(2, -1); // initialize with -1 a vector of size 2

  idx[0] = indexClose;
  idx[1] = (indexClose==1) ? 0 : 1 ;

  return idx;

}

float Minv3massiveCorr(const Vec_f& jet_pt, const Vec_f& jet_eta, const Vec_f& jet_phi, const Vec_f& jet_m,
                       const Vec_f& jet_cRegCorr,
                       int& jet_muonIdx1, int& jet_muonIdx2,
                       const Vec_f& muon_pt, const Vec_f& muon_eta, const Vec_f& muon_phi,
                       const Vec_f& jpsi_muon1_pt, const Vec_f& jpsi_muon1_eta, const Vec_f& jpsi_muon1_phi,
                       const Vec_f& jpsi_muon2_pt, const Vec_f& jpsi_muon2_eta, const Vec_f& jpsi_muon2_phi,
                       const Vec_f& jpsi_pt, const Vec_f& jpsi_eta, const Vec_f& jpsi_phi, const Vec_f& jpsi_m,
                       int indexCloseScale) {

  const float muon_mass_ = 0.1057;
  const float jetCloseScale_ = 0.90;

  float invMass = 0.;
  PtEtaPhiMVector mu1(muon_pt[jet_muonIdx1], muon_eta[jet_muonIdx1], muon_phi[jet_muonIdx1], muon_mass_);
  PtEtaPhiMVector mu2(muon_pt[jet_muonIdx2], muon_eta[jet_muonIdx2], muon_phi[jet_muonIdx2], muon_mass_);
  PtEtaPhiMVector p_jpsi(jpsi_pt[0], jpsi_eta[0], jpsi_phi[0], jpsi_m[0]);

  // to add is a DR between the two muons
  if (indexCloseScale == 0) {
    PtEtaPhiMVector p1(jetCloseScale_*jet_pt[0], jet_eta[0], jet_phi[0], jet_m[0]);
    PtEtaPhiMVector p2(jet_pt[1]*jet_cRegCorr[1], jet_eta[1], jet_phi[1], jet_m[1]);
    invMass = (p1 + p2 + p_jpsi - mu1 - mu2).M();
  } else {
    PtEtaPhiMVector p1(jet_pt[0]*jet_cRegCorr[0], jet_eta[0], jet_phi[0], jet_m[0]);
    PtEtaPhiMVector p2(jetCloseScale_*jet_pt[1], jet_eta[1], jet_phi[1], jet_m[1]);
    invMass = (p1 + p2 + p_jpsi - mu1 - mu2).M();
  }
  return invMass;
}


#endif
