#include "mysf.h"
#include <ostream>
#include <iostream>

MyCorrections::MyCorrections(int year) {

  std::string dirName = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/";

  std::string subDirName = "";
  if(year == 2018)  subDirName += "2018_UL/";
  if(year == 2017)  subDirName += "2017_UL/";
  if(year == 22016)  subDirName += "2016postVFP_UL/";
  if(year == 12016)  subDirName += "2016preVFP_UL/";  
  
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

  std::string fileNameELE = dirName+"EGM/"+subDirName+"electron.json.gz";
  auto csetELE = correction::CorrectionSet::from_file(fileNameELE);
  electronSF_ = csetELE->at("UL-Electron-ID-SF");

  std::string fileNameJEC = dirName+"JME/"+subDirName+"jet_jerc.json.gz";
  auto csetJEC = correction::CorrectionSet::from_file(fileNameJEC);

  std::string tagName = "";
  if(year == 2018) tagName = "Summer19UL18_V5_MC_L1L2L3Res_AK4PFchs";
  if(year == 2017) tagName = "Summer19UL17_V5_MC_L1L2L3Res_AK4PFchs";
  if(year == 22016) tagName = "Summer19UL16_V5_MC_L1L2L3Res_AK4PFchs";
  if(year == 12016) tagName = "Summer19UL16APV_V5_MC_L1L2L3Res_AK4PFchs";

  JEC_ = csetJEC->compound().at(tagName);

};

double MyCorrections::eval_jetCORR(
    double double1, double double2, double double3, double double4) {
  return JEC_->evaluate({double1, double2, double3, double4});
};

double MyCorrections::eval_electronSF(
    std::string str1, std::string str2,  std::string str3, double double1, double double2) {
  return electronSF_->evaluate({str1, str2, str3, double1, double2});
};

double MyCorrections::eval_photonSF(
    std::string str1, std::string str2,  std::string str3, double double1, double double2) {
  return photonSF_->evaluate({str1, str2, str3, double1, double2});
};

double MyCorrections::eval_puSF(
    double int1, std::string str1) {
  return puSF_->evaluate({int1, str1});
};
