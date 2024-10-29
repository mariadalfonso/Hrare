import ROOT
import json

import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score

myDir='/work/submit/mariadlf/Hrare_JPsiCC/Sept25/'

ROOT.gROOT.SetBatch(True)

weight_sf = 1e9

# lower learning rate = a little better? wait actually much worse oops
paramsRegress = {
    #'objective': 'reg:pseudohubererror',
    'objective': 'reg:squarederror',
    'eval_metric': 'mae',
    'max_depth': 6,
    'subsample': 0.5,
    'colsample_bytree': 0.5,
    'seed': 42,
    'n_estimators': 900,
    'learning_rate': 0.1,
    'gamma': 3,
    'min_child_weight': 10,
    'max_delta_step': 0,
    'early_stopping_rounds': 40,
}

# Use GridSearchCV for parameter tuning
param_grid = {
    'max_depth': [4, 5, 6],
    'learning_rate': [0.01, 0.1, 0.2],
    'n_estimators': [100, 300, 500],
    'gamma': [0, 1, 3],
    'min_child_weight': [1, 5, 10],
    'subsample': [0.5, 0.7, 1],
    'colsample_bytree': [0.5, 0.7, 1]
}

paramsClass = {
    'objective': 'binary:logistic',
    'eval_metric': 'logloss',
    'eta': 0.1,
    'max_depth': 4,
    'subsample': 0.5,
    'colsample_bytree': 0.5,
    'seed': 42,
    'n_estimators': 4000,
    'early_stopping_rounds': 25,
    'num_rounds': 20,
    'learning_rate': 0.1,
    'gamma': 3,
    'min_child_weight': 10,
    'max_delta_step': 0,
}

variables_noangles = [
    "jetClose_CvL",
    "jetFar_CvL",
    "jetCloseJPsiRatio",
    "jetFar_PTMHRatio",
    "jpsi_Higgs_Ratio"
]

variables_close =  [
    "jetClose_Pt",
    "jetClose_Eta",
    "jetClose_Phi",
    "jetClose_Mass",

    "jetClose_CvL",
    "jetClose_nConst",

    "MET_pt",
    "MET_phi",

]

variables_far =  [
    "jetFar_Pt",
    "jetFar_Eta",
    "jetFar_Phi",
    "jetFar_Mass",

    "jetFar_CvL",
    "jetFar_nConst",

    "MET_pt",
    "MET_phi",
]


def load_process_class(fIn, variables, target=0, weight=1.):

    # Heavy-lifting in C++ and remote access of data
    df = ROOT.RDataFrame("events", fIn)

    cols = df.AsNumpy(variables)
    # Push data to scipy ecosystem

    pdf = pd.DataFrame(cols)
    pdf['target'] = target # add a target column to indicate signal (1) and background (0)
    pdf['weight'] = weight
    return pdf

def load_process_regr(fIn, variables, target=0, weight=1.):

    # Heavy-lifting in C++ and remote access of data
    df = ROOT.RDataFrame("events", fIn)

    cols = df.AsNumpy(variables+[target])
    # Push data to scipy ecosystem

    pdf = pd.DataFrame(cols)
    pdf['target'] = target # add a target column to indicate signal (1) and background (0)
    pdf['weight'] = weight
    return pdf

def _test_XGB_regr(label, testCase):

    # TO DO: 1) refine the variables

    variables = variables_close
    target_var = "jetClose_charm_ptratio"
    fOutName = "HJpsiCC_regression_model_close"

    if testCase=='far':
        variables = variables_far
        target_var = "jetFar_charm_ptratio"
        fOutName = "HJpsiCC_regression_model_far"

    sig_df = load_process_regr(myDir+"snapshotJpsiCC_1000_2018.root", variables, target=target_var, weight=weight_sf)

    data = pd.concat([sig_df], ignore_index=True)

    train_data, test_data, train_targets, test_targets, train_weights, test_weights  = train_test_split(
        data[variables], data[target_var], data['weight'], test_size = 0.1, random_state = 42)

    train_data = train_data.to_numpy()
    test_data = test_data.to_numpy()
    train_targets = train_targets.to_numpy()
    test_targets = test_targets.to_numpy()
    train_weights = train_weights.to_numpy()
    test_weights = test_weights.to_numpy()
    eval_set = [(train_data, train_targets), (test_data, test_targets)]

    regressor = xgb.XGBRegressor(**paramsRegress)
    regressor.fit(train_data, train_targets, eval_set=eval_set, sample_weight=train_weights)

    model_name = "regress_model"
    print("variables",variables)
    print("Export model ",model_name)

    ROOT.TMVA.Experimental.SaveXGBoost(regressor, model_name, fOutName+".root", num_inputs=len(variables))
    print(f"output written to {fOutName} with name {model_name}")
    regressor.save_model(f'{fOutName}.json')

    variables_ = ROOT.TList()
    for var in variables+[target_var]:
        print(var)
        variables_.Add(ROOT.TObjString(var))
    fOut = ROOT.TFile(fOutName+".root", "UPDATE")
    fOut.WriteObject(variables_, "variables")

def _test_XGB_class(label):

    # TO DO: 1) refine the variable 2) split even and odd

    # ROC with angles: 0.89/0.88
    variables = variables_noangles
    variables_plus_HM = variables + ["massHiggsCorr"]

    #####
    sig_df_hm = load_process_class(myDir+"snapshotJpsiCC_1000_2018.root", variables_plus_HM, weight=weight_sf, target=1)
    bkg1_df_hm = load_process_class(myDir+"snapshotJpsiCC_10_2018.root", variables_plus_HM, weight=weight_sf)
    bkg2_df_hm = load_process_class(myDir+"snapshotJpsiCC_11_2018.root", variables_plus_HM, weight=weight_sf)

    data = pd.concat([sig_df_hm, bkg1_df_hm, bkg2_df_hm], ignore_index=True)

    train_data_hm, test_data_hm, train_labels_hm, test_labels_hm, train_weights_hm, test_weights_hm  = train_test_split(
        data[variables_plus_HM], data['target'], data['weight'], test_size=0.2, random_state=42)

    train_data_hm = train_data_hm.to_numpy()
    #####

    correlation_matrix = np.corrcoef(train_data_hm, rowvar=False)
    correlation_df = pd.DataFrame(correlation_matrix, index=variables_plus_HM, columns=variables_plus_HM)
    correlation_df_rounded = correlation_df.round(3)
    correlation_df_rounded.to_csv("correlation_matrix.txt", sep='\t')
    print("Correlation matrix saved to correlation_matrix.txt")

    ### repeated block
    sig_df = load_process(myDir+"snapshotJpsiCC_1000_2018.root", variables, weight=weight_sf, target=1)
    bkg1_df = load_process(myDir+"snapshotJpsiCC_10_2018.root", variables, weight=weight_sf)
    bkg2_df = load_process(myDir+"snapshotJpsiCC_11_2018.root", variables, weight=weight_sf)
    
    data = pd.concat([sig_df, bkg1_df, bkg2_df], ignore_index=True)
    
    train_data, test_data, train_labels, test_labels, train_weights, test_weights  = train_test_split(
        data[variables], data['target'], data['weight'], test_size=0.2, random_state=42)

    train_data = train_data.to_numpy()
    test_data = test_data.to_numpy()
    train_labels = train_labels.to_numpy()
    test_labels = test_labels.to_numpy()
    train_weights = train_weights.to_numpy()
    test_weights = test_weights.to_numpy()
    eval_set = [(train_data, train_labels), (test_data, test_labels)]
    ###

    print("Start training")
    bdt = xgb.XGBClassifier(**paramsClass)
    bdt.fit(train_data, train_labels, verbose=True, eval_set=eval_set, sample_weight=train_weights)
    
    fOutName = "HJpsiCC_classification_model.root"
    model_name = "bdt_model"
    print("variables",variables)
    print("Export model ",model_name)

    ROOT.TMVA.Experimental.SaveXGBoost(bdt, model_name, fOutName, num_inputs=len(variables))
    print(f"output written to {fOutName} with name {model_name}")
    
    # append the variables
    
    variables_ = ROOT.TList()
    for var in variables:
        print(var)
        variables_.Add(ROOT.TObjString(var))
    fOut = ROOT.TFile(fOutName, "UPDATE")
    fOut.WriteObject(variables_, "variables")

if __name__ == "__main__":

    _test_XGB_class("default")
    _test_XGB_regr("default","close")
    _test_XGB_regr("default","far")
