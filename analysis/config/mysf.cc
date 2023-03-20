#include "mysf.h"
#include <ostream>
#include <iostream>

MyCorrections::MyCorrections(int year) {

  std::string dirName = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/";

  std::string subDirName = "";
  std::string dataName = "";

  if(year == 2018)  { subDirName += "2018_UL/"; dataName = "UL18"; }
  if(year == 2017)  { subDirName += "2017_UL/"; dataName = "UL17"; }
  if(year == 22016)  { subDirName += "2016postVFP_UL/"; dataName = "UL16"; }
  if(year == 12016)  { subDirName += "2016preVFP_UL/"; dataName = "UL16APV"; }
  
  std::string fileNameLUM = dirName+"LUM/"+subDirName+"puWeights.json.gz";

  std::string corrNameLUM = "";  
  if(year == 2018) corrNameLUM = "Collisions18_UltraLegacy_goldenJSON";
  if(year == 2017) corrNameLUM = "Collisions17_UltraLegacy_goldenJSON";
  if(year == 22016 or year == 12016) corrNameLUM = "Collisions16_UltraLegacy_goldenJSON";
  
  auto csetPU = correction::CorrectionSet::from_file(fileNameLUM);
  puSF_ = csetPU->at(corrNameLUM);
  
  std::string fileNamePH = dirName+"EGM/"+subDirName+"photon.json.gz";
  auto csetPH = correction::CorrectionSet::from_file(fileNamePH);
  photonSF_ = csetPH->at("UL-Photon-ID-SF");
  photonPixVetoSF_ = csetPH->at("UL-Photon-PixVeto-SF");

  std::string fileNameELE = dirName+"EGM/"+subDirName+"electron.json.gz";
  auto csetELE = correction::CorrectionSet::from_file(fileNameELE);
  electronSF_ = csetELE->at("UL-Electron-ID-SF");

  std::string fileNameMU = dirName+"MUO/"+subDirName+"muon_Z.json.gz";
  auto csetMu = correction::CorrectionSet::from_file(fileNameMU);
  muonTRKSF_ = csetMu->at("NUM_TrackerMuons_DEN_genTracks");
  muonIDTSF_ = csetMu->at("NUM_TightID_DEN_genTracks");
  muonIDMSF_ = csetMu->at("NUM_MediumID_DEN_genTracks");
  muonISOTSF_ = csetMu->at("NUM_TightRelIso_DEN_TightIDandIPCut");
  muonISOLSF_ = csetMu->at("NUM_LooseRelIso_DEN_MediumID");

  std::string fileNameJEC = dirName+"JME/"+subDirName+"jet_jerc.json.gz";
  auto csetJEC = correction::CorrectionSet::from_file(fileNameJEC);

  std::string tagName = "Summer19"+dataName+"_V5_MC_L1L2L3Res_AK4PFchs";
  if(year == 22016 or year == 12016) tagName = "Summer19"+dataName+"_V7_MC_L1L2L3Res_AK4PFchs";
  JEC_ = csetJEC->compound().at(tagName);

  std::string tagNameUnc = "Summer19"+dataName+"_V5_MC_Total_AK4PFchs";
  if(year == 22016 or year == 12016) tagNameUnc = "Summer19"+dataName+"_V7_MC_Total_AK4PFchs";
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

double MyCorrections::eval_jetCORR(
    double area, double eta, double pt, double rho) {
  return JEC_->evaluate({area, eta, pt, rho});
};

double MyCorrections::eval_jesUnc(
    double eta, double pt, int type) {
  if(type == 0) return jesUnc_->evaluate({eta, pt});
  return 0.0;
};

double MyCorrections::eval_jer(
    double double1, double double2, double double3, double double4) {
  return JER_->evaluate({double1, double2, double3, double4});
};

double MyCorrections::eval_jetVeto(
    std::string str1, double double1, double double2) {
  return vetoMaps_->evaluate({str1,double1, double2});
};

double MyCorrections::eval_electronSF(
    std::string year, std::string valType,  std::string workingPoint, double eta, double pt) {
  pt = std::max(pt,10.001);
  return electronSF_->evaluate({year, valType, workingPoint, eta, pt});
};

double MyCorrections::eval_photonSF(
    std::string year, std::string valType,  std::string workingPoint, double eta, double pt) {
  pt = std::max(pt,20.001);
  return photonSF_->evaluate({year, valType, workingPoint, eta, pt});
};

double MyCorrections::eval_photonPixVetoSF(
    std::string year, std::string valType,  std::string workingPoint, double eta, double pt) {
  pt = std::max(pt,20.001);
  return photonPixVetoSF_->evaluate({year, valType, workingPoint, eta, pt});
};

double MyCorrections::eval_muonTRKSF(std::string year, std::string valType, double eta, double pt) {
  eta = std::min(std::abs(eta),2.399);
  pt = std::max(pt,15.001);
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

double MyCorrections::eval_puSF(
    double int1, std::string str1) {
  return puSF_->evaluate({int1, str1});
};
