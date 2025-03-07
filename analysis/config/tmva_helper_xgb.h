#ifndef HELPER_TMVA_H
#define HELPER_TMVA_H

// ROOT                                                                                                                                                                                                                                                                                                                  
#include "ROOT/RVec.hxx"

using Vec_f = ROOT::VecOps::RVec<float>;


class tmva_helper_xgb {
    public:
        tmva_helper_xgb(const std::string &filename, const std::string &name, const unsigned &nvars, const unsigned int nslots = 1) :
            model_(name, filename), nvars_(nvars) {

            const unsigned int nslots_actual = std::max(nslots, 1U);

            interpreters_.reserve(nslots_actual);
            for (unsigned int islot = 0; islot < nslots_actual; ++islot) {
                interpreters_.emplace_back(model_);
            }
        }

        Vec_f operator()(unsigned int slot, const Vec_f vars) {

	    return interpreters_[slot].Compute(vars);
        }

    private:
        unsigned int nvars_;
        TMVA::Experimental::RBDT model_;
        std::vector<TMVA::Experimental::RBDT> interpreters_;

};

#endif
