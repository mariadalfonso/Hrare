#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "TROOT.h"
#include "TFile.h"
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


#include <ROOT/RVec.hxx>
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

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > PtEtaPhiMVector;
std::unordered_map< UInt_t, std::vector< std::pair<UInt_t,UInt_t> > > jsonMap;

Vec_i indices(const int& size, const int& start = 0) {
    Vec_i res(size, 0);
    std::iota(std::begin(res), std::end(res), start);
    return res;
}


Vec_f NomUpDownVar(const float up, const float down, const float nom, float weight) {

  Vec_f res(3, 1);
  res[0] = weight;  // nom
  res[1] = weight*up/nom;  // up
  res[2] = weight*down/nom;  // down
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

  int idxMatch = -1;
  float etaGen = 999.;
  float phiGen = 999.;
  // loop over all the genPartCand
  for (unsigned int i=0; i<genPart_pdgId.size(); i++) {
    if(pdgToMatch==310  and abs(genPart_pdgId[i])==pdgToMatch) {
      if(genPart_pdgId[i] == pdgToMatch && genPart_genPartIdxMother[i] == idxIntermediate ) {
	idxMatch=i;  etaGen=genPart_eta[i]; phiGen=genPart_phi[i];
      }
    } else if(genPart_pdgId[i]==22 or abs(genPart_pdgId[i])==333 or abs(genPart_pdgId[i])==113  or abs(genPart_pdgId[i])==443) {
      if(genPart_pdgId[i] == pdgToMatch && genPart_genPartIdxMother[i] == idxMother ) {
	idxMatch=i;  etaGen=genPart_eta[i]; phiGen=genPart_phi[i];
      }
    }
  }

  Vec_i idx(2, -1); // initialize with -1 a vector of size 2

  idx[1] = idxMatch; // index of the genPart (-1 if not found gen cand)

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

stdVec_i HiggsCandFromRECO(const Vec_f& meson_pt, const Vec_f& meson_eta, const Vec_f& meson_phi, const Vec_f& meson_mass,
			const Vec_f& meson_trk1_pt, const Vec_f& meson_trk2_pt,
			const Vec_f& wrong_meson_pt,
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

  // loop over all the phiCand
  for (unsigned int i=0; i<meson_pt.size(); i++) {

    if(max(meson_trk1_pt[i], meson_trk2_pt[i]) < 20) continue;

    PtEtaPhiMVector p_meson(meson_pt[i], meson_eta[i], meson_phi[i], meson_mass[i]);
    if((p_meson + p_ph).M()<5.) continue; // object disambiguiation, also remove the omega/tau resonances

    // save the leading Pt
    float ptCand = p_meson.pt();

    if(ptCandMax < ptWrongMax) continue; // we want the leading meson to the of the right flavor
    if( ptCand < ptCandMax ) continue;
    ptCandMax=ptCand;
    Minv = (p_meson + p_ph).M();
    ptHiggs = (p_meson + p_ph).pt();
    idx[0] = i;
    idx[1] = indexPhoton;
  }

  return idx;
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

#endif
