import ROOT
import json

import matplotlib
matplotlib.use("Agg")     # IMPORTANT for batch mode / no display
import matplotlib.pyplot as plt
import scipy.stats

import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve, auc

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score

from LoadTree import loadTreeRun3

myDir='/work/submit/mariadlf/HrareRun3/JAN26/2024/VBFcat/'

years = ['_2024']
category="_VBFcat"
#category="ggHcat"
mesonCat = '_RhoCat'
mytree = ROOT.TChain('events')
for year in years:
    print('year=',year)
    mytree = loadTreeRun3(mytree, myDir, category, mesonCat, year )

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT()

paramsClassSUGGESTnotworking = {
    "objective": "binary:logistic",
    "eval_metric": "auc",

    # Tree complexity
    "max_depth": 3,              # shallower trees generalize better
    "min_child_weight": 5,      # prevents overfitting small fluctuations
    "gamma": 1.0,               # require meaningful split improvement

    # Learning control
    "eta": 0.05,                # slower learning = better separation
    "n_estimators": 1200,       # compensate lower eta

    # Sampling (VERY important)
    "subsample": 0.8,
    "colsample_bytree": 0.8,

    # Regularization (important)
    "lambda": 2.0,              # stronger L2 improves robustness
    "alpha": 0.5,               # small L1 helps ignore useless vars

    # Performance
    "tree_method": "hist",
    "max_bin": 256,

    "seed": 42,
}


paramsClass = {
    "lambda": 0, # to make the discriminator larger for VL
    "objective": "binary:logistic",
    "eval_metric": "auc",
    "max_depth": 4,
    "eta": 0.1,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "n_estimators": 400,
    "tree_method": "hist",
    "seed": 42,
    "min_child_weight": 1,
}

signal_map = {
    "_RhoCat": {
        "_VBFcat": ["1020"],
        "GF":  ["1027", "1029"],
    },
    "_PhiCat": {
        "_VBFcat": ["1010"],
        "GF":  ["1017", "1019"],
    },
    "_K0StarCat": {
        "_VBFcat": ["1030"],
        "GF":  ["1037"],
    },
    "_D0StarCat": {
        "_VBFcat": ["1040"],
        "GF":  ["1047"],
    },
}

bkg_map = {
    "_VBFcat": [str(i) for i in range(10, 25, 1)],
    "_GFcat": [str(i) for i in range(10, 25, 1)],    
}

dir_map = {
    "_VBFcat": "/home/submit/mariadlf/public_html/HrareRun3_MVA/VBF/",
    "_ggHcat": "/home/submit/mariadlf/public_html/HrareRun3_MVA/ggH/",
}

variables_map = {
    "_VBFcat": ["Rpt","mJJ","dEtaJJ","zepVar","HCandPT","photon_pt","meson_pt","sigmaHCandMass_Rel2"],
    "_ggHcat": ["HiggsCandCorrPt", "HiggsCandMassErr", "Muon1_pt", "Muon2_pt"],
}

# here the no discrimination and very weak
variables_notUseful_map = {
    "_ggHcat": ["HiggsCandCorrRapidity", "cosThetaCS", "phiStarCS", "Muon1_eta","Muon2_eta", "PuppiMET_pt"],
    "_VBFcat": ["HiggsCandCorrRapidity", "HiggsCandMassErr", "cosThetaCS","phiStarCS","Muon1_eta","Muon2_eta","dPhiJJ"],
}

def get_var_range(df, var, pad=0.05):
    """Return (xmin, xmax) for a variable using RDataFrame Min/Max."""
    min_val = df.Min(var).GetValue()
    max_val = df.Max(var).GetValue()

    # Add small padding around edges (5% by default)
    span = max_val - min_val
    xmin = min_val - span * pad
    xmax = max_val + span * pad

    return xmin, xmax

def make_plots(df):

    # Define weight
    df = df.Define("weight", "w_allSF")

    # Define signal & background filters
    sig_ids = signal_map[mesonCat][category]
    bkg_ids = bkg_map[category]
    sig_filter = " || ".join([f"mc == {mid}" for mid in sig_ids])
    bkg_filter = " || ".join([f"mc == {mid}" for mid in bkg_ids])
    df_sig = df.Filter(sig_filter, f"Signal ({sig_filter})")
    df_bkg = df.Filter(bkg_filter, f"Signal ({bkg_filter})")
    df_sig = df_sig.Define("absweight", "fabs(weight)")
    df_bkg = df_bkg.Define("absweight", "fabs(weight)")

    for var in variables_map[category]: #+ variables_notUseful_map[category]:
        print(f"Plotting {var}...")

        # Get automatic min/max from **both** samples together
        xmin_s, xmax_s = get_var_range(df_sig, var)
        xmin_b, xmax_b = get_var_range(df_bkg, var)
        xmin = min(xmin_s, xmin_b)
        xmax = max(xmax_s, xmax_b)

        nbins = 40  # if you want different binning per variable, I can automate this too

        # Create histograms
        hsig = df_sig.Histo1D(
            (f"h_{var}_sig", f"{var} (Signal);{var};Normalized Events", nbins, xmin, xmax),
            var, "weight"
        )
        hbkg = df_bkg.Histo1D(
            (f"h_{var}_bkg", f"{var} (Background);{var};Normalized Events", nbins, xmin, xmax),
            var, "weight"
        )
        hsigABS = df_sig.Histo1D(
            (f"h_{var}_abssig", f"{var} (Signal) ABSw;{var};Normalized Events", nbins, xmin, xmax),
            var, "absweight"
        )
        hbkgABS = df_bkg.Histo1D(
            (f"h_{var}_absbkg", f"{var} (Background) ABSw;{var};Normalized Events", nbins, xmin, xmax),
            var, "absweight"
        )

        # Draw
        c = ROOT.TCanvas(f"c_{var}", "", 800, 600)
        hsig.SetLineColor(ROOT.kOrange)
        hsigABS.SetLineColor(ROOT.kRed)
        hbkg.SetLineColor(ROOT.kBlue)
        hbkgABS.SetLineColor(ROOT.kGreen+1)

        hsig.Scale(1.0 / hsig.Integral())
        hsigABS.Scale(1.0 / hsigABS.Integral())
        hbkg.Scale(1.0 / hbkg.Integral())
        hbkgABS.Scale(1.0 / hbkgABS.Integral())
        hsig.SetMaximum(1.5*max(hsig.GetMaximum(),hbkg.GetMaximum()))
        hsig.Draw("HIST")
        hsigABS.Draw("HIST SAME")
        hbkg.Draw("HIST SAME")
        hbkgABS.Draw("HIST SAME")

        # Explicit legend (middle right)
        leg = ROOT.TLegend(0.60, 0.40, 0.88, 0.65)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.04)

        leg.AddEntry(hsigABS.GetValue(), "sig", "l")
        leg.AddEntry(hbkg.GetValue(), "bkg", "l")
        leg.Draw()
        #c.BuildLegend()

        c.SaveAs(dir_map[category]+var+mesonCat+"_sig_bkg.png")

        print(f" → Saved {var}_sig_bkg.png")

def load_process_class(isSignal, variables, drawPlot=False):
#, target=0, weight=1.):
    
    # Heavy-lifting in C++ and remote access of data
    df = ROOT.RDataFrame(mytree)
    #df = df.Filter("(HiggsCandCorrMass>(125-20) and HiggsCandCorrMass<(125+20))","HiggsMass within reasonable range 125+-20")
    # FOR VH - VL - TTH - TTL
    #df = df.Filter("(HiggsCandCorrMass>(125-50) and HiggsCandCorrMass<(125+50))","HiggsMass within reasonable range 125+-50")
    df = df.Filter("(HCandMass>(125-25) and HCandMass<(125+25))","HiggsMass within reasonable range 125+-25")
    if isSignal and drawPlot: make_plots(df) # only make the plot once

    sig_ids = signal_map[mesonCat][category]
    sig_filter = " || ".join([f"mc == {mid}" for mid in sig_ids])
    bkg_ids = bkg_map[category]
    bkg_filter = " || ".join([f"mc == {mid}" for mid in bkg_ids])
    df = df.Filter(sig_filter if isSignal else bkg_filter)
    # Define the weight column (if not already present)
    df = df.Define("weight", "w_allSF")
#    df = df.Filter("((jetVBF1_Pt > 35 || jetVBF1_Pt > 35)) && jetVBF1_Pt > 25 && jetVBF1_Pt > 25","at least one > 35 GeV and both above 25")
#    df = df.Filter("((jetVBF1_Pt > 35 || jetVBF1_Pt > 35)) && jetVBF1_Pt > 25 && jetVBF1_Pt > 25","at least one > 35 GeV and both above 25")        
    nevts = df.Count().GetValue()
    print('case isSignal=',isSignal,' -- nevts',nevts)

    cols = df.AsNumpy(variables + ["weight"])
    # Push data to scipy ecosystem
    pdf = pd.DataFrame(cols)
    pdf['target'] = 1 if isSignal else 0
    return pdf

def diagnostic(bdt,test_data,test_labels,variables):

    # Probabilities (what you actually want)
    y_score = bdt.predict_proba(test_data)[:, 1]

    print("Score range:", y_score.min(), y_score.max())

    roc = roc_auc_score(test_labels, bdt.predict_proba(test_data)[:, 1])
    print(f"ROC AUC: {roc:.4f}")

    # Compute ROC curve and AUC
    fpr, tpr, thresholds = roc_curve(test_labels, y_score)
    roc_auc = auc(fpr, tpr)

    # Plot ROC curve
    plt.figure(figsize=(8,6))
    plt.plot(fpr, tpr, color='blue', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('XGBoost ROC Curve')
    plt.legend(loc="lower right")
    plt.grid(alpha=0.3)
    outfile = dir_map[category] + f"roc.png"
    plt.savefig(outfile, dpi=200)
    plt.close()
    print(f"📊 Saved roc to: {outfile}")

    # Set feature names for export
    bdt.get_booster().feature_names = variables
    print("Number of input vars:", len(variables))
    print("Booster feature names:", bdt.get_booster().feature_names)
    importance = bdt.get_booster().get_score(importance_type="gain")
    print("Non-zero features:", importance)

    importance_types = ["gain", "weight", "cover"]
    titles = {
        "gain": "XGBoost Feature Importance — Gain",
        "weight": "XGBoost Feature Importance — Weight (Split Count)",
        "cover": "XGBoost Feature Importance — Cover"
    }

    for imp in importance_types:
        plt.figure(figsize=(8, 6))
        xgb.plot_importance(bdt, importance_type=imp, max_num_features=20)
        plt.title(titles[imp])
#        plt.tight_layout()
        plt.subplots_adjust(left=0.35)

        outfile = dir_map[category] + f"feature_importance_{imp}.png"
        plt.savefig(outfile, dpi=200)
        plt.close()

        print(f"📊 Saved feature importance ({imp}) to: {outfile}")

def diagnostic636(bdt,dtest,test_labels,variables):

    """
    bdt        : xgboost.core.Booster
    dtest      : xgboost.DMatrix
    test_labels: numpy array (0/1)
    variables  : list of physics variable names
    """

    print('INSIDE DIAGNOSTIC')
    print("Number of input vars:", len(variables))
    print("Booster feature names:", bdt.feature_names)
    print("Booster attributes:", bdt.attributes())

    scores = bdt.predict(dtest)
    print("score range:", scores.min(), scores.max())
    print(type(bdt))

    y_pred = (scores > 0.5).astype(int)
    acc = accuracy_score(test_labels, y_pred)
    roc = roc_auc_score(test_labels, scores)
    print(f"Accuracy : {acc:.4f}")
    print(f"ROC AUC  : {roc:.4f}")

    # Compute ROC curve and AUC
    fpr, tpr, thresholds = roc_curve(test_labels, y_pred)
    roc_auc = auc(fpr, tpr)

    # Plot ROC curve
    plt.figure(figsize=(8,6))
    plt.plot(fpr, tpr, color='blue', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('XGBoost ROC Curve')
    plt.legend(loc="lower right")
    plt.grid(alpha=0.3)
    plt.show()

    importance_types = ["gain", "weight", "cover"]
    titles = {
        "gain": "XGBoost Feature Importance — Gain",
        "weight": "XGBoost Feature Importance — Weight (Split Count)",
        "cover": "XGBoost Feature Importance — Cover"
    }

    # Set feature names for export
    bdt.feature_names = variables

    for imp in importance_types:
        scores_imp = bdt.get_score(importance_type=imp)

        if not scores_imp:
            print(f"⚠️  No non-zero importances for {imp}")
            continue

        # Map f0 → real variable names

        mapped = {
            variables[int(k[0:])]: v
            for k, v in scores_imp.items()
        }

        print(f"\nNon-zero features ({imp}):")
        for k, v in sorted(mapped.items(), key=lambda x: -x[1]):
            print(f"  {k:25s} : {v:.4e}")

        # Plot
        plt.figure(figsize=(8, 6))
        xgb.plot_importance(
            bdt,
            importance_type=imp,
            max_num_features=20,
            xlabel=imp,
            show_values=False,
        )
        plt.title(titles[imp])
        plt.tight_layout()

        outfile = dir_map[category] + f"feature_importance_{imp}.png"
        plt.savefig(outfile, dpi=200)
        plt.close()

        print(f"📊 Saved feature importance ({imp}) to: {outfile}")

def _test_XGB_class(label):

    # -------------------------------------------------
    # Main training dataset (no Higgs mass)
    # -------------------------------------------------

    variables = variables_map[category] #+ variables_notUseful_map[category]

    sig_df = load_process_class(True,  variables, False)
    bkg_df = load_process_class(False, variables, False)

    data = pd.concat([sig_df, bkg_df], ignore_index=True).sample(frac=1, random_state=42)
    data["event"] = np.arange(len(data))
    #data["event"] = data.index

    # -------------------------------------------------
    # Build balanced, positive-definite weights
    # -------------------------------------------------

    # Build balanced weight
    data["weight_balanced"] = data["weight"]

    # Fix negative weights (MANDATORY)
    n_neg = (data["weight_balanced"] <= 0).sum()
    if n_neg:
        print(f"⚠️  Found {n_neg} non-positive weights → taking absolute value")
        
        data["weight_balanced"] = data["weight_balanced"].abs()
        data["weight_balanced"] = data["weight_balanced"].clip(lower=1e-6)

    # -------------------------------------------------
    # Split even odd
    # -------------------------------------------------

    train_mask = (data["event"] % 2 == 0)
    test_mask  = ~train_mask

    train_data   = data.loc[train_mask, variables]
    test_data    = data.loc[test_mask,  variables]

    train_labels = data.loc[train_mask, "target"]
    test_labels  = data.loc[test_mask,  "target"]

    train_weights = data.loc[train_mask, "weight_balanced"]
    test_weights  = data.loc[test_mask,  "weight_balanced"]

    # -------------------------------------------------
    # Correlation check with HiggsCandCorrMass
    # -------------------------------------------------

    variables_hm = variables_map[category] + ["HCandMass"]

    sig_df_hm = load_process_class(True,  variables_hm, True)
    bkg_df_hm = load_process_class(False, variables_hm, False)

    data_hm = pd.concat([sig_df_hm, bkg_df_hm], ignore_index=True)
    data_hm["event"] = np.arange(len(data_hm))

    train_mask_hm = (data_hm["event"] % 2 == 0)

    train_hm = data_hm.loc[train_mask_hm, variables_hm].to_numpy()

    corr = np.corrcoef(train_hm, rowvar=False)
    corr_df = pd.DataFrame(corr, index=variables_hm, columns=variables_hm).round(3)

    out_corr = f"output/correlation_matrix{category}.txt"
    corr_df.to_csv(out_corr, sep="\t")

    print(f"📐 Correlation matrix saved to {out_corr}")

    # -------------------------------------------------
    # Class imbalance (MANDATORY)
    # -------------------------------------------------

    n_sig = np.sum(train_labels == 1)
    n_bkg = np.sum(train_labels == 0)

    if n_sig == 0 or n_bkg == 0:
        raise RuntimeError("ERROR: signal or background count is zero")

    scale_pos_weight = n_bkg / n_sig
    paramsClass["scale_pos_weight"] = scale_pos_weight

    print(f"Signal events     : {n_sig}")
    print(f"Background events : {n_bkg}")
    print(f"scale_pos_weight  : {scale_pos_weight:.3f}")

    # -------------------------
    # some diagnostic i.e. of the variables
    # -------------------------    
    print("Train/test sizes:", train_data.shape, test_data.shape)
    print("Train signal fraction:", train_labels.mean(),
          "Test signal fraction:", test_labels.mean())

    print("Train labels unique:", np.unique(train_labels))
    print("Test labels unique:", np.unique(test_labels))
    print("Train positive fraction:", train_labels.mean())
    print("Test positive fraction:", test_labels.mean())

    # Basic statistics
    for v in variables:
        col = train_data[v]  # Series
        print(
            f"{v:20s}",
            "min =", col.min(),
            "max =", col.max(),
            "std =", col.std(),
    )

    '''
    KS      Interpretation
    < 0.05  No separation
    0.05–0.15  Very weak
    0.15–0.30  Some discrimination
    > 0.30  Strong
    '''

    # KS test between signal and background
    for v in variables:
        sig = train_data.loc[train_labels == 1, v]
        bkg = train_data.loc[train_labels == 0, v]

        ks_stat = scipy.stats.ks_2samp(sig, bkg).statistic

        if ks_stat < 0.05: comment='No separation'
        if ks_stat > 0.05 and ks_stat < 0.15: comment='very weak'
        if ks_stat > 0.15 and ks_stat < 0.30: comment='some discrimination'
        if ks_stat > 0.30: comment='STRONG'

        print(
            f"{v:20s}",
            "⟨sig⟩ =", sig.mean(),
            "⟨bkg⟩ =", bkg.mean(),
            "KS =", ks_stat,
            "==> ",comment
        )

    # Sanity checks
    print("Any NaNs in train_data:", train_data.isna().any().any())
    print("Any infs in train_data:", np.isinf(train_data.values).any())

    print("Mean train weight:", train_weights.mean())
    print("Min train weight :", train_weights.min())
    print("Max train weight :", train_weights.max())

    train_data = train_data.to_numpy()
    test_data = test_data.to_numpy()
    train_weights = train_weights.to_numpy()

    train_labels = train_labels.to_numpy()
    test_labels = test_labels.to_numpy()
    test_weights = test_weights.to_numpy()

    ###
    # -------------------------

    print("Start training in 632")

    eval_set = [(train_data, train_labels), (test_data, test_labels)]
    bdt = xgb.XGBClassifier(**paramsClass)
    bdt.fit(train_data, train_labels, sample_weight=train_weights, verbose=True, eval_set=eval_set)

    print("Training complete.")

    fOutName = f"output/classification_model_{category}.root"
    model_name = f"bdt_model_{category}"
    print("variables",variables)
    print("Export model ",model_name)

    ROOT.TMVA.Experimental.SaveXGBoost(bdt, model_name, fOutName, num_inputs=len(variables))
    print(f"output written to {fOutName} with name {model_name}")

    ###
    # -------------------------

    '''
    print("Start training in 636")

    feature_names = [f"f{i}" for i in range(len(variables))]

    dtrain = xgb.DMatrix(
        train_data,
        label=train_labels,
        weight=train_weights,
        feature_names=feature_names,
#        feature_names=variables,
    )

    dtest = xgb.DMatrix(
        test_data,
        label=test_labels,
        weight=test_weights,
        feature_names=feature_names,
#        feature_names=variables,
    )

    print('CALLING xgb.train')
    bdt = xgb.train(
        paramsClass,
        dtrain,
        num_boost_round=200,
        evals=[(dtrain, "train"), (dtest, "test")],
        verbose_eval=True
    )

#    eval_set = [(train_data, train_labels), (test_data, test_labels)]
#    eval_weights = [train_weights, test_weights]
    
#    bdt = xgb.XGBClassifier(**paramsClass)
#    bdt.fit(
#        train_data,
#        train_labels,
#        sample_weight=train_weights,
#        eval_set=[(test_data, test_labels)],
#        sample_weight_eval_set=[test_weights],
#        verbose=True
#    )

    print("Training complete.")

    fOutName = f"output/classification_model_{category}_TEST.root"
    model_name = f"bdt_model_{category}"
    print("variables",variables)
    print("Export model ",model_name)
    print(bdt.feature_names)

    clf = xgb.XGBClassifier()
    clf._Booster = bdt

   # Attach metadata ROOT expects
    #    tmva_features = [f"f{i}" for i in range(len(variables))]
    clf._Booster.feature_names = feature_names
    #    clf._Booster.feature_names = variables 
    clf._Booster.objective = "binary:logistic"
    
    ROOT.TMVA.Experimental.SaveXGBoost(
        clf,
        model_name,
        fOutName,
        num_inputs=len(variables)
    )
    
    print(f"✅ XGBoost model exported to {fOutName} as {model_name}")
    '''

#    ROOT.TMVA.Experimental.SaveXGBoost(bdt, model_name, fOutName, num_inputs=len(variables))
#    print(f"output written to {fOutName} with name {model_name}")

    # append the variables

    variables_ = ROOT.TList()
    for var in variables:
        print(var)
        variables_.Add(ROOT.TObjString(var))
    fOut = ROOT.TFile(fOutName, "UPDATE")
    fOut.WriteObject(variables_, "variables")
    print('FILE SAVED')
    
    ## call diagnostic
#    diagnostic636(bdt,dtest,test_labels,variables)
    diagnostic(bdt,test_data,test_labels,variables)    

if __name__ == "__main__":

    _test_XGB_class("default")
