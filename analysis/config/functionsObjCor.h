#ifndef MYFUN                                                                                                                                                                       
#define MYFUN                                                                                                                                                                       

using Vec_b = ROOT::VecOps::RVec<bool>;
using Vec_f = ROOT::VecOps::RVec<float>;

Vec_f computeJECcorrection(MyCorrections corrSFs, Vec_f jet_pt, Vec_f jet_rawFactor, Vec_f jet_eta, Vec_f jet_phi, Vec_f jet_area, float rho, float run, bool isData, string year, \
string mc){                                                                                                                                                                                 
  Vec_f new_jet; new_jet.resize(jet_pt.size());
  Vec_f raw_jet; raw_jet.resize(jet_pt.size());
  for (unsigned int idx = 0; idx < jet_pt.size(); ++idx) {                                                                                                                            
    raw_jet[idx] = jet_pt[idx] * (1.0 - jet_rawFactor[idx]);                                                                                                                       
    new_jet[idx] = raw_jet[idx] * corrSFs.eval_jetCORR(jet_area[idx], jet_eta[idx], jet_phi[idx], raw_jet[idx], rho, isData, run, year, mc );                                      
  }                                                                                                                                                                                   
  return new_jet;                                                                                                                                                                     
}                                                                                                                                                                                   

Vec_f computeJECuncertainties(MyCorrections corrSFs, Vec_f jet_pt, Vec_f jet_eta){                                                                                                  
  Vec_f new_jet_delta; new_jet_delta.resize(jet_pt.size());
  int type = 0;                                                                                                                                                                       
  for (unsigned int idx = 0; idx < jet_pt.size(); ++idx) new_jet_delta[idx] = corrSFs.eval_jesUnc(jet_eta[idx], jet_pt[idx], type );                                                  
  return new_jet_delta;                                                                                                                                                               
}                                                                                                                                                                                   

Vec_b cleaningJetVetoMapMask(const Vec_f& jet_eta, const Vec_f& jet_phi, const string year) {                                                                                       
  Vec_b jet_vetoMap_mask(jet_eta.size(), true);                                                                                                                                       
  for (unsigned int idx = 0; idx < jet_eta.size(); ++idx) {                                                                                                                           
    double jetVetoMap = corr_sf.eval_jetVeto(jet_eta[idx], jet_phi[idx]);                                                                                                               
    if(jetVetoMap > 0) jet_vetoMap_mask[idx] = false;                                                                                                                                   
  }                                                                                                                                                                                   
  return jet_vetoMap_mask;                                                                                                                                                            
}                                                                                                                                                                                   

Vec_b cleaningJetIDMask(Vec_f jet_eta, Vec_f jet_chHEF, Vec_f jet_neHEF, Vec_f jet_chEmEF, Vec_f jet_neEmEF, Vec_f jet_muEF, Vec_f jet_chMultiplicity, Vec_f jet_neMultiplicity, string year) {                                                                                                                                                                                
  Vec_b jetID_mask(jet_eta.size(), true);                                                                                                                                             
  for (unsigned int idx = 0; idx < jet_eta.size(); ++idx) {                                                                                                                           
    double jetID = corr_sf.eval_jetID(jet_eta[idx], jet_chHEF[idx], jet_neHEF[idx], jet_chEmEF[idx], jet_neEmEF[idx], jet_muEF[idx], jet_chMultiplicity[idx], jet_neMultiplicity[idx]); 
    if (jetID < 0) jetID_mask[idx] = false;                                                                                                                                             
  }                                                                                                                                                                                   
  return jetID_mask;                                                                                                                                                                  
}                                                                                                                                                                                   
#endif                 
