class tmva_xml {

    public:
        tmva_xml(const std::string &filename) {

            auto c = TMVA::Experimental::Internal::ParseXMLConfig(filename);
			fVariables = c.variables;
			fExpressions = c.expressions;
			fAnalysisType = c.analysisType;
			fNumClasses = c.numClasses;

            fReader = std::make_unique<TMVA::Reader>("Silent");
			const auto numVars = fVariables.size();
            fValues = std::vector<float>(numVars);

            for(std::size_t i = 0; i < numVars; i++) {
                fReader->AddVariable(TString(fExpressions[i]), &fValues[i]);
            }
            fReader->BookMVA(name, filename.c_str());

        }

        std::vector<float> Compute(const Vec_f &x) {

            if (x.size() != fVariables.size())
                throw std::runtime_error("Size of input vector is not equal to number of variables.");

            // Copy over inputs to memory used by TMVA reader
            for (std::size_t i = 0; i < x.size(); i++) {
                fValues[i] = x[i];
            }

            // Evaluate TMVA model
            // Classification
            if (fAnalysisType == TMVA::Experimental::Internal::AnalysisType::Classification) {
                return std::vector<float>({static_cast<float>(fReader->EvaluateMVA(name))});
            }
            // Regression
            else if (fAnalysisType == TMVA::Experimental::Internal::AnalysisType::Regression) {
                return fReader->EvaluateRegression(name);
            }
            // Multiclass
            else if (fAnalysisType == TMVA::Experimental::Internal::AnalysisType::Multiclass) {
                return fReader->EvaluateMulticlass(name);
            }
            // Throw error
            else {
                throw std::runtime_error("RReader has undefined analysis type.");
                return std::vector<float>();
            }
        }

    private:
        std::unique_ptr<TMVA::Reader> fReader;
		std::vector<float> fValues;
		std::vector<std::string> fVariables;
		std::vector<std::string> fExpressions;
		unsigned int fNumClasses;
		const char *name = "RReader";
		TMVA::Experimental::Internal::AnalysisType fAnalysisType;
};

class tmva_helper_xml {
    public:
        tmva_helper_xml(const std::string &filename, const unsigned int nslots = 1) {

            const unsigned int nslots_actual = std::max(nslots, 1U);

			for (unsigned int islot = 0; islot < nslots_actual; ++islot) {
                tmva_xml *tmp = new tmva_xml(filename);
                interpreters_.emplace_back(tmp);
            }
        }

        std::vector<float> operator()(unsigned int slot, const Vec_f vars) {
			return interpreters_[slot]->Compute(vars);
        }

    private:
        std::vector<tmva_xml *> interpreters_;

};
