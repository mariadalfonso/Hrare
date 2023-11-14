#ifndef mysf_h
#define mysf_h

#include "correction.h"
#include <stdio.h>
#include <string.h>
#include <ostream>
#include <iostream>

//g++ $(correction config --cflags --ldflags) mysf.cc -shared -fPIC -o mysf.so
//g++ $(correction config --cflags --ldflags --rpath) mysf.cc -shared -fPIC -o mysf.so
//g++ $(correction config --cflags --ldflags --rpath) demo.cc -o demo

class MyCorrections {
public:
  MyCorrections(int year);

  double eval_jetCORR   (double area, double eta, double pt, double rho);
  double eval_jesUnc    (double eta, double pt, int type);
  double eval_jer       (double pt, double eta, double rho, double area);
  double eval_jetVeto   (std::string str1, double pt, double eta);
  double eval_puSF      (double NumTrueInteractions, std::string weights);
  double eval_photonSF  (std::string year, std::string valType, std::string workingPoint, double eta, double pt);
  double eval_photonPixVetoSF  (std::string year, std::string valType, std::string workingPoint, double eta, double pt);
  double eval_electronSF(std::string year, std::string valType, std::string workingPoint, double eta, double pt);
  double eval_muonTRKSF (std::string year, std::string valType, double eta, double pt);
  double eval_muonIDSF  (std::string year, std::string valType, double eta, double pt, std::string workingPoint);
  double eval_muonISOSF (std::string year, std::string valType, double eta, double pt, std::string workingPoint);

private:
  correction::Correction::Ref puSF_;
  correction::Correction::Ref photonSF_;
  correction::Correction::Ref photonPixVetoSF_;
  correction::Correction::Ref electronSF_;
  correction::Correction::Ref muonTRKSF_;
  correction::Correction::Ref muonIDMSF_;
  correction::Correction::Ref muonISOLSF_;
  correction::Correction::Ref muonIDTSF_;
  correction::Correction::Ref muonISOTSF_;
  correction::Correction::Ref JER_;
  correction::Correction::Ref JERsf_;
  correction::CompoundCorrection::Ref JEC_;
  correction::Correction::Ref jesUnc_;
  correction::Correction::Ref vetoMaps_;

};

#endif
