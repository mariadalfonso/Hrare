//#include "mysf.h"
//#include <ostream>
//#include <iostream>

#ifndef sfDistr_h
#define sfDistr_h

#include "correction.h"
#include <stdio.h>
#include <string.h>
#include <ostream>
#include <iostream>

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


MyCorrections::MyCorrections(int year) {

  std::string dirName = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/";

  std::string subDirName = "";
  std::string dataName = "";

  if(year == 2018)  { subDirName += "2018_UL/"; dataName = "UL18"; }
  if(year == 2017)  { subDirName += "2017_UL/"; dataName = "UL17"; }
  if(year == 22016)  { subDirName += "2016postVFP_UL/"; dataName = "UL16"; }
  if(year == 12016)  { subDirName += "2016preVFP_UL/"; dataName = "UL16APV"; }
  if(year == 12022)  { subDirName = "2022_Summer22/"; }
  if(year == 22022)  { subDirName = "2022_Summer22EE/"; }

  const std::string fileNameLUM = dirName+"LUM/"+subDirName+"puWeights.json.gz";

  std::string corrNameLUM = "";  
  if(year == 2018) corrNameLUM = "Collisions18_UltraLegacy_goldenJSON";
  if(year == 2017) corrNameLUM = "Collisions17_UltraLegacy_goldenJSON";
  if(year == 22016 or year == 12016) corrNameLUM = "Collisions16_UltraLegacy_goldenJSON";
  // PU missing for 2022/2023

  auto csetPU = correction::CorrectionSet::from_file(fileNameLUM);
  puSF_ = csetPU->at(corrNameLUM);

  const std::string fileNamePH = dirName+"EGM/"+subDirName+"photon.json.gz";
  auto csetPH = correction::CorrectionSet::from_file(fileNamePH);
  photonSF_ = csetPH->at("UL-Photon-ID-SF");
  photonPixVetoSF_ = csetPH->at("UL-Photon-PixVeto-SF");
  // EGM missing forr 2022/2023

  const std::string fileNameELE = dirName+"EGM/"+subDirName+"electron.json.gz";
  auto csetELE = correction::CorrectionSet::from_file(fileNameELE);
  electronSF_ = csetELE->at("UL-Electron-ID-SF");
  // ELE missing for 2023/2023

  std::string fileNameMU = dirName+"MUO/"+subDirName+"muon_Z_v2.json.gz";
  if(year == 12022) fileNameMU = dirName+"MUO/"+"2022_27Jun2023/"+"muon_Z.json.gz";
  if(year == 22022) fileNameMU = dirName+"MUO/"+"2022EE_27Jun2023/"+"muon_Z.json.gz";

  auto csetMu = correction::CorrectionSet::from_file(fileNameMU);
  muonTRKSF_ = csetMu->at("NUM_TrackerMuons_DEN_genTracks");
  muonIDTSF_ = csetMu->at("NUM_TightID_DEN_genTracks");
  muonIDMSF_ = csetMu->at("NUM_MediumID_DEN_genTracks");
  muonISOTSF_ = csetMu->at("NUM_TightRelIso_DEN_TightIDandIPCut");
  muonISOLSF_ = csetMu->at("NUM_LooseRelIso_DEN_MediumID");

  const std::string fileNameJEC = dirName+"JME/"+subDirName+"jet_jerc.json.gz";
  auto csetJEC = correction::CorrectionSet::from_file(fileNameJEC);

  std::string tagName = "Summer19"+dataName+"_V5_MC_L1L2L3Res_AK4PFchs";
  if(year == 22016 or year == 12016) tagName = "Summer19"+dataName+"_V7_MC_L1L2L3Res_AK4PFchs";
  if(year == 12022)  tagName = "Summer22_22Sep2023_V2_MC_L1L2L3Res_AK4PFPuppi";
  if(year == 22022)  tagName = "Summer22EE_22Sep2023_V2_MC_L1L2L3Res_AK4PFPuppi";
  // likely also the data are needed for 2022
  JEC_ = csetJEC->compound().at(tagName);

  std::string tagNameUnc = "Summer19"+dataName+"_V5_MC_Total_AK4PFchs";
  if(year == 22016 or year == 12016) tagNameUnc = "Summer19"+dataName+"_V7_MC_Total_AK4PFchs";
  if(year == 12022)  tagNameUnc = "Summer22_22Sep2023_V2_MC_L1L2L3Res_AK4PFPuppi";
  if(year == 22022)  tagNameUnc = "Summer22EE_22Sep2023_V2_MC_Total_AK4PFPuppi";
  jesUnc_ = csetJEC->at(tagNameUnc);

  std::string tagNameR = "Summer19"+dataName+"_JRV2_MC_PtResolution_AK4PFchs";
  if(year == 22016 or year == 12016) tagNameR = "Summer20"+dataName+"_JRV3_MC_PtResolution_AK4PFchs";
  JER_ = csetJEC->at(tagNameR);

  std::string tagNameRsf = "Summer19"+dataName+"_JRV2_MC_ScaleFactor_AK4PFchs";
  if(year == 22016 or year == 12016) tagNameRsf = "Summer20"+dataName+"_JRV3_MC_ScaleFactor_AK4PFchs";
  JERsf_ = csetJEC->at(tagNameRsf);

  /*
  // veto the jet
  std::string fileNameJetVeto = dirName+"JME/"+subDirName+"jetvetomaps.json.gz";
  auto csetVeto = correction::CorrectionSet::from_file(fileNameJetVeto);
  std::string tagNameVeto = "Summer19"+dataName+"_V1";
  if(year == 12016) tagNameVeto = "Summer19UL16_V1";   // the same ?
  vetoMaps_ = csetVeto->at(tagNameVeto);

  // puJetID
  std::string fileNamePUJetID = dirName+"JME/"+subDirName+"jmar.json.gz";
  auto csetPUJetID = correction::CorrectionSet::from_file(fileNamePUJetID);
  puJetIDSF_ = csetPUJetID->at("PUJetID_eff");
  */

};

double MyCorrections::eval_jetCORR(double area, double eta, double pt, double rho) {
  return JEC_->evaluate({area, eta, pt, rho});
};

double MyCorrections::eval_jesUnc(double eta, double pt, int type) {
  if(type == 0) return jesUnc_->evaluate({eta, pt});
  return 0.0;
};

double MyCorrections::eval_jer(double double1, double double2, double double3, double double4) {
  return JER_->evaluate({double1, double2, double3, double4});
};

double MyCorrections::eval_jetVeto(std::string str1, double double1, double double2) {
  return vetoMaps_->evaluate({str1,double1, double2});
};

double MyCorrections::eval_electronSF(std::string year, std::string valType,  std::string workingPoint, double eta, double pt) {
  pt = std::max(pt,10.001);
  return electronSF_->evaluate({year, valType, workingPoint, eta, pt});
};

double MyCorrections::eval_photonSF(std::string year, std::string valType,  std::string workingPoint, double eta, double pt) {
  pt = std::max(pt,20.001);
  return photonSF_->evaluate({year, valType, workingPoint, eta, pt});
};

double MyCorrections::eval_photonPixVetoSF(std::string year, std::string valType,  std::string workingPoint, double eta, double pt) {
  pt = std::max(pt,20.001);
  return photonPixVetoSF_->evaluate({year, valType, workingPoint, eta, pt});
};

double MyCorrections::eval_muonTRKSF(std::string year, std::string valType, double eta, double pt) {
  eta = std::min(std::abs(eta),2.399);
  pt = std::max(pt,20.001);
  return muonTRKSF_->evaluate({year, eta, pt, valType});
};

double MyCorrections::eval_muonIDSF(std::string year, std::string valType, double eta, double pt, std::string workingPoint) {
  eta = std::min(std::abs(eta),2.399);
  pt = std::max(pt,15.001);
  
  if (workingPoint=="T") {
    return muonIDTSF_->evaluate({year, eta, pt, valType});
  } else if (workingPoint=="M") {
    return muonIDMSF_->evaluate({year, eta, pt, valType});
  }
  return 1.;
};

double MyCorrections::eval_muonISOSF(std::string year, std::string valType, double eta, double pt, std::string workingPoint) {
  eta = std::min(std::abs(eta),2.399);
  pt = std::max(pt,15.001);
  
  if (workingPoint=="T") {
    return muonISOTSF_->evaluate({year, eta, pt, valType});
  } else if (workingPoint=="L") {
    return muonISOLSF_->evaluate({year, eta, pt, valType});
  }
  return 1.;
};

double MyCorrections::eval_puSF(double int1, std::string str1) {
  return puSF_->evaluate({int1, str1});
};

#endif


