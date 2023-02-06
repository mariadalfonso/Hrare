#include "correction.h"
#include <stdio.h>
#include <string.h>

//g++ $(correction config --cflags --ldflags) mysf.cc -shared -fPIC -o mysf.so
//g++ $(correction config --cflags --ldflags --rpath) mysf.cc -shared -fPIC -o mysf.so
//g++ $(correction config --cflags --ldflags --rpath) demo.cc -o demo

class MyCorrections {
public:
  MyCorrections(int year);

  double eval_jetCORR(double pt, double eta, double rho, double area);
  double eval_puSF(double NumTrueInteractions, std::string weights);
  double eval_photonSF(std::string year, std::string valType, std::string workingPoint, double eta, double pt);
  double eval_electronSF(std::string year, std::string valType, std::string workingPoint, double eta, double pt);

private:
  correction::Correction::Ref puSF_;
  correction::Correction::Ref photonSF_;
  correction::Correction::Ref electronSF_;
  correction::CompoundCorrection::Ref JEC_;

};
