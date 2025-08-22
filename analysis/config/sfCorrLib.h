//#include "mysf.h"
//#include <ostream>
//#include <iostream>

#ifndef sfDistr_h
#define sfDistr_h

//#include "correction.h"
#include <stdio.h>
#include <string.h>
#include <ostream>
#include <iostream>

class MyCorrections {

public:
  MyCorrections(int year);

  double eval_jetCORR   (double area, double eta, double phi, double pt, double rho, bool isData, int run, std::string year, std::string mc);
  double eval_jesUnc    (double eta, double pt, int type);
  double eval_jer       (double pt, double eta, double rho, double area);
  double eval_jetID     (float eta, float chHEF, float neHEF, float chEmEF, float neEmEF, float muEF, int chMultiplicity, int neMultiplicity);
  double eval_jetVeto   (double eta, double phi);
  double eval_puSF      (double NumTrueInteractions, std::string weights);
  double eval_photonSF  (std::string year, std::string valType, std::string workingPoint, double eta, double pt);
  double eval_photonPixVetoSF  (std::string year, std::string valType, std::string workingPoint, double eta, double pt);
  double eval_electronSF(std::string year, std::string valType, std::string workingPoint, double eta, double pt, double minVal);
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
  correction::CompoundCorrection::Ref JECdata_;
  correction::Correction::Ref jesUnc_;
  correction::Correction::Ref vetoMaps_;
  correction::Correction::Ref jetTightID_;

};


MyCorrections::MyCorrections(int year) {

  std::string dirName = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/";

  std::string subDirName = "";
  std::string dataName = "";

  if(year == 2018)  { subDirName += "2018_UL/"; dataName = "UL18"; }
  if(year == 2017)  { subDirName += "2017_UL/"; dataName = "UL17"; }
  if(year == 22016)  { subDirName += "2016postVFP_UL/"; dataName = "UL16"; }
  if(year == 12016)  { subDirName += "2016preVFP_UL/"; dataName = "UL16APV"; }
  if(year == 12022)  { subDirName += "2022_Summer22/"; }
  if(year == 22022)  { subDirName += "2022_Summer22EE/"; }
  if(year == 12023)  { subDirName += "2023_Summer23/"; }
  if(year == 22023)  { subDirName += "2023_Summer23BPix/"; }
  if(year == 2024)   { subDirName += "2024_Summer24/"; }

  const std::string fileNameLUM = dirName+"LUM/"+subDirName+"puWeights.json.gz";

  std::string corrNameLUM = "";  
  if(year == 2018) corrNameLUM = "Collisions18_UltraLegacy_goldenJSON";
  if(year == 2017) corrNameLUM = "Collisions17_UltraLegacy_goldenJSON";
  if(year == 22016 or year == 12016) corrNameLUM = "Collisions16_UltraLegacy_goldenJSON";
  if(year == 12022) corrNameLUM = "Collisions2022_355100_357900_eraBCD_GoldenJson";
  if(year == 22022) corrNameLUM = "Collisions2022_359022_362760_eraEFG_GoldenJson";
  if(year == 12023) corrNameLUM = "Collisions2023_366403_369802_eraBC_GoldenJson";
  if(year == 22023) corrNameLUM = "Collisions2023_369803_370790_eraD_GoldenJson";

  auto csetPU = correction::CorrectionSet::from_file(fileNameLUM);
  puSF_ = csetPU->at(corrNameLUM);

  const std::string fileNamePH = dirName+"EGM/"+subDirName+"photon.json.gz";

  auto csetPH = correction::CorrectionSet::from_file(fileNamePH);
  std::string corrNamePH = "Photon-ID-SF";
  std::string corrNamePH2 = "Photon-PixVeto-SF";
  if(year == 2018 or year == 2017 or year == 12016 or year == 22016) corrNamePH = "UL-Photon-ID-SF";
  if(year == 2018 or year == 2017 or year == 12016 or year == 22016) corrNamePH2 = "UL-Photon-PixVeto-SF";
  photonSF_ = csetPH->at(corrNamePH);
  photonPixVetoSF_ = csetPH->at(corrNamePH2);
  // note scale and smearing need to be applied in Run3 from the JSON file

  const std::string fileNameELE = dirName+"EGM/"+subDirName+"electron.json.gz";
  auto csetELE = correction::CorrectionSet::from_file(fileNameELE);
  std::string corrNameEGM = "Electron-ID-SF";
  if(year == 2018 or year == 2017 or year == 12016 or year == 22016) corrNameEGM = "UL-Electron-ID-SF";
  electronSF_ = csetELE->at(corrNameEGM);

  std::cout << " ELE done " << std::endl;

  std::string fileNameMU = dirName+"MUO/"+subDirName+"muon_Z.json.gz";
  if(year == 2018 or year == 2017 or year == 12016 or year == 22016) fileNameMU = dirName+"MUO/"+subDirName+"muon_Z_v2.json.gz";
  if(year == 2018 or year == 2017 or year == 12016 or year == 22016) fileNameMU = "/home/submit/mariadlf/Hrare/CMSSW_10_6_27_new/src/Hrare/analysis/config/POG/MUO/"+subDirName+"muon_Z.json.gz";

  auto csetMu = correction::CorrectionSet::from_file(fileNameMU);
  if(year == 2018 or year == 2017 or year == 12016 or year == 22016) {
    muonTRKSF_ = csetMu->at("NUM_TrackerMuons_DEN_genTracks");
    muonIDTSF_ = csetMu->at("NUM_TightID_DEN_genTracks");
    muonIDMSF_ = csetMu->at("NUM_MediumID_DEN_genTracks");
    muonISOTSF_ = csetMu->at("NUM_TightRelIso_DEN_TightIDandIPCut");
    muonISOLSF_ = csetMu->at("NUM_LooseRelIso_DEN_MediumID");
  }

  std::cout << " MUO to do for Run3 " << std::endl;

  const std::string fileNameJEC = dirName+"JME/"+subDirName+"jet_jerc.json.gz";
  auto csetJEC = correction::CorrectionSet::from_file(fileNameJEC);

  if(year == 12022 or year == 22022 or year == 12023 or year == 22023) {

    const std::string jetType="AK4PFPuppi";

    std::string tagName = "";
    if(year == 12022)  tagName = "Summer22_22Sep2023_RunCD_V2";
    if(year == 22022)  tagName = "Summer22EE_22Sep2023_RunE_V2";
    if(year == 12023)  tagName = "Summer23Prompt23_RunCv123_V1";
    if(year == 22023)  tagName = "Summer23BPixPrompt23_RunD_V1";
    JECdata_ = csetJEC->compound().at(tagName+"_DATA_L1L2L3Res_"+jetType);

    std::string tagNameMC = "";
    if(year == 12022)  tagNameMC = "Summer22_22Sep2023_V2";
    if(year == 22022)  tagNameMC = "Summer22EE_22Sep2023_V2";
    if(year == 12023)  tagNameMC = "Summer23Prompt23_V1";
    if(year == 22023)  tagNameMC = "Summer23BPixPrompt23_V1";
    JEC_ = csetJEC->compound().at(tagNameMC+"_MC_L1L2L3Res_"+jetType);
    jesUnc_ = csetJEC->at(tagNameMC+"_MC_Total_"+jetType);

  }

  if(year == 2018 or year == 2017 or year == 12016 or year == 22016) {

    const std::string jetType="AK4PFchs";
    const std::string prefix="Summer19";

    std::string tagName = prefix+dataName+"_RunA_V5";
    //    if(year == 22016 or year == 12016) tagName = prefix+dataName+"_V7";
    if(year == 12016) tagName = prefix+dataName+"_RunBCD_V7";
    if(year == 22016) tagName = prefix+dataName+"_RunFGH_V7";
    if(year == 2017) tagName = prefix+dataName+"_RunB_V5";
    // likely also the data are needed for 2022
    JECdata_ = csetJEC->compound().at(tagName+"_DATA_L1L2L3Res_"+jetType);

    std::string tagNameMC = prefix+dataName+"_V5";
    if(year == 22016 or year == 12016) tagNameMC = prefix+dataName+"_V7";
    // likely also the data are needed for 2022
    JEC_ = csetJEC->compound().at(tagNameMC+"_MC_L1L2L3Res_"+jetType);

    std::string tagNameUnc = "Summer19"+dataName+"_V5";
    if(year == 22016 or year == 12016) tagNameUnc = prefix+dataName+"_V7";
    jesUnc_ = csetJEC->at(tagNameUnc+"_MC_Total_"+jetType);

    std::string tagNameR = prefix+dataName+"_JRV2";
    if(year == 22016 or year == 12016) tagNameR = "Summer20"+dataName+"_JRV3";
    JER_ = csetJEC->at(tagNameR+"_MC_PtResolution_"+jetType);

    std::string tagNameRsf = prefix+dataName+"_JRV2";
    if(year == 22016 or year == 12016) tagNameRsf = "Summer20"+dataName+"_JRV3";
    JERsf_ = csetJEC->at(tagNameRsf+"_MC_ScaleFactor_"+jetType);

  }

  // veto Map the jet
  std::string fileNameJetVeto = dirName+"JME/"+subDirName+"jetvetomaps.json.gz";
  auto csetVeto = correction::CorrectionSet::from_file(fileNameJetVeto);
  std::string tagNameVeto = "";
  if(year == 22016 or year == 2017 or year == 2018) tagNameVeto = "Summer19"+dataName+"_V1";
  if(year == 12016) tagNameVeto = "Summer19UL16_V1";
  if(year == 12022) tagNameVeto = "Summer22_23Sep2023_RunCD_V1";
  if(year == 22022) tagNameVeto = "Summer22EE_23Sep2023_RunEFG_V1";
  if(year == 12023) tagNameVeto = "Summer23Prompt23_RunC_V1";
  if(year == 22023) tagNameVeto = "Summer23BPixPrompt23_RunD_V1";
  if(year == 2024) tagNameVeto = "Summer24Prompt24_RunBCDEFGHI_V1";
  vetoMaps_ = csetVeto->at(tagNameVeto);

  // jetID
  std::string fileNameJetID = dirName+"JME/"+subDirName+"jetid.json.gz";
  auto csetJetID = correction::CorrectionSet::from_file(fileNameJetID);
  const std::string tagNameJetID = "AK4PUPPI_Tight";
  //  std::string tagNameJetID = "AK4PUPPI_TightLeptonVeto";
  jetTightID_           = csetJetID->at(tagNameJetID);

  /*
  // puJetID
  std::string fileNamePUJetID = dirName+"JME/"+subDirName+"jmar.json.gz";
  auto csetPUJetID = correction::CorrectionSet::from_file(fileNamePUJetID);
  puJetIDSF_ = csetPUJetID->at("PUJetID_eff");
  */

};

double MyCorrections::eval_jetCORR(double area, double eta, double phi, double pt, double rho, bool isData, int run, std::string year, std::string mc) {

  if (year == "2024" or year == "22023") {
    if(isData) return JECdata_->evaluate({area, eta,  pt, rho, phi, (float) run});
    else JEC_->evaluate({area, eta, pt, rho, phi});
  } else if (year == "12023") {
    if(isData) return JECdata_->evaluate({area, eta, pt, rho, (float) run});
    else return JEC_->evaluate({area, eta, pt, rho});
    /*
  } else if (year == "22022") {
    if(isData and mc == "-15") return JECdata22E_->evaluate({area, eta, pt, rho});
    if(isData and mc == "-16") return JECdata22F_->evaluate({area, eta, pt, rho});
    if(isData and mc == "-17") return JECdata22G_->evaluate({area, eta, pt, rho});
    else return JEC_->evaluate({area, eta, pt, rho});
    */
  } else {
    if(isData) return JECdata_->evaluate({area, eta, pt, rho});
    else return JEC_->evaluate({area, eta, pt, rho});
  }

  return 1.0;

};

double MyCorrections::eval_jesUnc(double eta, double pt, int type) {
  if(type == 0) return jesUnc_->evaluate({eta, pt});
  return 0.0;
};

double MyCorrections::eval_jer(double double1, double double2, double double3, double double4) {
  return JER_->evaluate({double1, double2, double3, double4});
};

double MyCorrections::eval_jetVeto(double eta, double phi) {
  std::string typeMaps = "jetvetomap";
  return vetoMaps_->evaluate({typeMaps,eta, phi});
};

double MyCorrections::eval_jetID(float eta, float chHEF, float neHEF, float chEmEF, float neEmEF, float muEF, int chMultiplicity, int neMultiplicity) {
  eta = fabs(eta);
  int multiplicity = chMultiplicity + neMultiplicity;
  //"chEmEF" and "muEF" unused in 2024 but still needed
  return jetTightID_->evaluate({eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity});
};

double MyCorrections::eval_electronSF(std::string year, std::string valType,  std::string workingPoint, double eta, double pt, double minVal) {
  //  pt = std::max(pt,10.001); for wp80 is 10 while for the above10 is 20
  pt = std::max(pt,minVal);
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


