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


typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;


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


float MesonCandFromRECO(const Vec_f& meson_pt, const Vec_f& meson_eta, const Vec_f& meson_phi, const Vec_f& meson_mass,
                        const Vec_f& ph_pt, const Vec_f& ph_eta, const Vec_f& ph_phi) {

  float Minv = -1;
  int index = -1;
  float ptCandMax=0;
  PtEtaPhiMVector p_ph(ph_pt[0], ph_eta[0], ph_phi[0], 0);

  if(ph_pt.size()> 0) {
    if(ph_pt[1] > ph_pt[0]) p_ph.SetPt(ph_pt[1]);
    if(ph_pt[1] > ph_pt[0]) p_ph.SetEta(ph_eta[1]);
    if(ph_pt[1] > ph_pt[0]) p_ph.SetPhi(ph_phi[1]);
  }

  // loop over all the phiCand
  for (unsigned int i=0; i<meson_pt.size(); i++) {

    PtEtaPhiMVector p_meson(meson_pt[i], meson_eta[i], meson_phi[i], meson_mass[i]);

    if(abs(deltaPhi(p_ph.phi(), p_meson.phi()))<float(M_PI/2)) continue; // M,gamma opposite hemishpere
    if(p_meson.pt()<5) continue;  // Higgs candidate is a two body decay

    // save the leading Pt
    float ptCand = p_meson.pt();
    if( ptCand < ptCandMax ) continue;
    ptCandMax=ptCand;
    Minv = (p_meson + p_ph).mass();
    index = i;
  }

  return Minv;

}

Vec_i mesonCand(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m, const Vec_f& ch,
		const Vec_f& ph_pt, const Vec_f& ph_eta, const Vec_f& ph_phi,
		bool phiHyp
		) {
  
  Vec_i idx(2, -1); // initialize with -1 a vector of size 2

  float minPtTracks_= phiHyp ? 10: 0.f;
  float minPtMeson_= phiHyp ? 10: 20.f;

  PtEtaPhiMVector pPhoton(ph_pt[0], ph_eta[0], ph_phi[0], 0.);

  //  PtEtaPhiMVector MesonCand(0,0,0,0);
  float mass=-100;
  float ptCandMax=0;
  for (unsigned int i=0; i<pt.size(); i++) {

    for (unsigned int j=0; j<i; j++) {
      
      if(i==j) continue;
      if(ch[i]*ch[j]>0) continue; // opposite sign kaons

      if(pt[i]<minPtTracks_ and pt[j]<minPtTracks_) continue; //at least a kaons of X GeV

      PtEtaPhiMVector p1(pt[i], eta[i], phi[i], m[i]);
      PtEtaPhiMVector p2(pt[j], eta[j], phi[j], m[j]);

      if(abs(deltaPhi(ph_phi[0], (p1 + p2).phi()))<float(M_PI/2)) continue; // M,gamma opposite hemishpere

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

float mt(float pt1, float phi1, float pt2, float phi2) {
  return std::sqrt(2*pt1*pt2*(1-std::cos(phi1-phi2)));
}


float Minv(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m) {
    PtEtaPhiMVector p1(pt[0], eta[0], phi[0], m[0]);
    PtEtaPhiMVector p2(pt[1], eta[1], phi[1], m[1]);
    return (p1 + p2).mass();
}

float Minv3(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m, Vec_i idx,
	    const Vec_f& ph_pt, const Vec_f& ph_eta, const Vec_f& ph_phi) {
 
  PtEtaPhiMVector p1(pt[idx[0]], eta[idx[0]], phi[idx[0]], m[idx[0]]);
  PtEtaPhiMVector p2(pt[idx[1]], eta[idx[1]], phi[idx[1]], m[idx[1]]);
  PtEtaPhiMVector p_ph(ph_pt[0], ph_eta[0], ph_phi[0], 0);
  return (p1 + p2 + p_ph).mass();
  
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


float dR_Constituents(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m, Vec_i idx) {
  return deltaR(eta[idx[0]], phi[idx[0]], eta[idx[1]], phi[idx[1]]);
}

#endif
