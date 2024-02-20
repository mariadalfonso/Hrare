import ROOT

class TMVAHelperXML():

    def __init__(self, model_input, model_name=""):
        model_tmp = ROOT.TMVA.Experimental.RReader(model_input)
        self.variables = [str(var) for var in model_tmp.GetVariableNames()]
        self.model_input = model_input
        self.model_name = model_name
        self.nthreads = ROOT.GetThreadPoolSize()

        self.tmva_helper = ROOT.tmva_helper_xml(self.model_input, self.nthreads)
        self.var_col = f"tmva_vars_{self.model_name}"

    def run_inference(self, df, col_name = "mva_score"):

        # check if columns exist in the dataframe
        cols = df.GetColumnNames()
        for var in self.variables:
            if not var in cols:
                raise Exception(f"Variable {var} not defined in dataframe.")

        vars_str = ', '.join(self.variables)
        df = df.Define(self.var_col, f"ROOT::VecOps::RVec<float>{{{vars_str}}}")
        df = df.DefineSlot(col_name, self.tmva_helper, [self.var_col])
        return df
