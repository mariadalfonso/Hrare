#include "correction.h"
#include <stdio.h>
#include <string.h>

//g++ $(correction config --cflags --ldflags) mysf.cc -shared -fPIC -o mysf.so 

class MyCorrections {
public:
  MyCorrections(int year);
  
  double eval_puSF(double NumTrueInteractions, std::string weights);
  double eval_photonSF(std::string year, std::string valType, std::string workingPoint, double eta, double pt);
  
private:
  correction::Correction::Ref puSF_;
  correction::Correction::Ref photonSF_;  

};
