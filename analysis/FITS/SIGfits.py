from ROOT import *

from prepareFits import getHisto

gSystem.Load("libHiggsAnalysisCombinedLimit.so")

blinded=False

xlowRange = 90.
xhighRange = 170.

x = RooRealVar('mh', 'm_{#gamma,meson}', xlowRange, xhighRange)

x.setRange("full", xlowRange  , xhighRange)
x.setRange("left", xlowRange, 115)
x.setRange("right", 135, xhighRange)

def  fitSig(tag , mesonCat, year):

    # Create a empty workspace
    w = RooWorkspace("w", "workspace")

    mcSig = ["WH", "ZH", "VBFH", "ZinvH", "ggH"]
    for sig in mcSig:

        # 1l: ZH, WH
        if (sig == "VBFH" and tag == "_Wcat"): continue
        if (sig == "ggH" and tag == "_Wcat"): continue
#        if (sig == "VBFHlow" and tag == "_Wcat"): continue
        if (sig == "ZinvH" and tag == "_Wcat"): continue
        # MET: WH, ZinvH
        if (sig == "ZH" and tag == "_Zinvcat"): continue
        if (sig == "VBFH" and tag == "_Zinvcat"): continue
#        if (sig == "VBFHlow" and tag == "_Zinvcat"): continue
        if (sig == "ggH" and tag == "_Zinvcat"): continue
        ##
        if (sig == "WH" and tag == "_Zcat"): continue
        if (sig == "VBFH" and tag == "_Zcat"): continue
        if (sig == "ggH" and tag == "_Zcat"): continue
#        if (sig == "VBFHlow" and tag == "_Zcat"): continue
        if (sig == "ZinvH" and tag == "_Zcat"): continue
        ##
        if (sig == "ZH" and tag == "_GFcat"): continue
        if (sig == "WH" and tag == "_GFcat"): continue
        if (sig == "ZinvH" and tag == "_GFcat"): continue
#        if (sig == "VBFH" and tag == "_GFcat"): continue
#        if (sig == "VBFHlow" and tag == "_GFcat"): continue
        ## ggh cat: GF and VBF
        if (sig == "WH" and tag == "_VBFcat"): continue
        if (sig == "ZH" and tag == "_VBFcat"): continue
        if (sig == "ZinvH" and tag == "_VBFcat"): continue
        if (sig == "ggH" and tag == "_VBFcat"): continue
#        if (sig == "VBFHlow" and tag == "_VBFcat"): continue
        ##
        if (sig == "WH" and tag == "_VBFcatlow"): continue
        if (sig == "ZH" and tag == "_VBFcatlow"): continue
        if (sig == "ZinvH" and tag == "_VBFcatlow"): continue
#        if (sig == "VBFH" and tag == "_VBFcatlow"): continue
        if (sig == "ggH" and tag == "_VBFcatlow"): continue

        print(tag, ' ', sig)

        data_full = getHisto(4, 200*10, 0. , 200., True, tag, mesonCat, True, sig)

        x = RooRealVar('mh', 'm_{#gamma,meson}', xlowRange, xhighRange)

        x.setRange("full", xlowRange, xhighRange)

        data = RooDataHist('datahist'+tag+'_'+sig, 'data', RooArgList(x), data_full)

        # -----------------------------------------------------------------------------

        cb_mu = RooRealVar('cb_mu'+tag+'_'+sig, 'cb_mu', 125., 125-10. , 125+10.)
        cb_sigma = RooRealVar('cb_sigma'+tag+'_'+sig, 'cb_sigma', 0., 3.)
        cb_alphaL = RooRealVar('cb_alphaL'+tag+'_'+sig, 'cb_alphaL', 0., 5.)
        cb_alphaR = RooRealVar('cb_alphaR'+tag+'_'+sig, 'cb_alphaR', 0., 5.)
        cb_nL = RooRealVar('cb_nL'+tag+'_'+sig, 'cb_nL', 2., 50.)
        cb_nR = RooRealVar('cb_nR'+tag+'_'+sig, 'cb_nR', 0., 5.)

        pdf_crystalball = RooDoubleCBFast('crystal_ball'+tag+'_'+sig, 'crystal_ball', x, cb_mu, cb_sigma, cb_alphaL, cb_nL, cb_alphaR, cb_nR)
        model = pdf_crystalball

        # -----------------------------------------------------------------------------

        model.fitTo(data,RooFit.Minimizer("Minuit2"),RooFit.Strategy(2),RooFit.Range("full"))

        # Here we will plot the results
        canvas = TCanvas("canvas", "canvas", 800, 800)

        plotFrameWithNormRange = x.frame(RooFit.Title("mH_"+tag+"_"+mesonCat+"_"+str(year)))

        # Plot only the blinded data, and then plot the PDF over the full range as well as both sidebands
        data.plotOn(plotFrameWithNormRange)
        model.plotOn(plotFrameWithNormRange, RooFit.LineColor(2), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineStyle(10))
        model.paramOn(plotFrameWithNormRange, RooFit.Layout(0.65,0.99,0.75))
        plotFrameWithNormRange.getAttText().SetTextSize(0.02);

        plotFrameWithNormRange.Draw()

        canvas.Draw()

        canvas.SaveAs("WS/signal_"+tag+"_"+mesonCat+"_"+str(year)+".png")

        # -----------------------------------------------------------------------------

        norm_SR = data_full.Integral(data_full.FindBin(xlowRange), data_full.FindBin(xhighRange))
        Sig_norm = RooRealVar(model.GetName()+ "_norm", model.GetName()+ "_norm", norm_SR) # no range means contants

        # -----------------------------------------------------------------------------
        # -----------------------------------------------------------------------------
        # Create workspace, import data and model

        cb_mu.setConstant()
        cb_sigma.setConstant()
        cb_alphaL.setConstant()
        cb_alphaR.setConstant()
        cb_nL.setConstant()
        cb_nR.setConstant()
        Sig_norm.setConstant()

        # Import model and all its components into the workspace
        getattr(w,'import')(model)

        getattr(w,'import')(Sig_norm)
        print('integral signal = ',Sig_norm.Print())

        # Import data into the workspace
        getattr(w,'import')(data)

        # Print workspace contents
        w.Print()

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Save workspace in file

    # Save the workspace into a ROOT file
    w.writeToFile("WS/Signal"+tag+"_"+mesonCat+"_"+str(year)+"_workspace.root")

def  fitBkg(tag , mesonCat, year):

    data_full = getHisto(4, 250, 0. , 250., True, tag, mesonCat, False, -1 )
    data = RooDataHist('datahist'+tag, 'data', RooArgList(x), data_full)

    blindedData = data.reduce(RooFit.CutRange("left,right"))

    # -----------------------------------------------------------------------------
    # BERN law
    bern_c0 = RooRealVar('bern_c0'+tag, 'bern_c0', 0.2, 0., 0.5)
    bern_c1 = RooRealVar('bern_c1'+tag, 'bern_c1', 0.1, 0., 1.)
    bern_c2 = RooRealVar('bern_c2'+tag, 'bern_c2', 0.1, 0., 1.)
    bern_c3 = RooRealVar('bern_c3'+tag, 'bern_c3', 0.5, 0., 5.)
    bern_c4 = RooRealVar('bern_c4'+tag, 'bern_c4', 0.5, 0., 5.)
    bern_c5 = RooRealVar('bern_c5'+tag, 'bern_c5', 1e-2, 0., 0.1)

    pdf_bern0 = RooBernstein('bern0'+tag, 'bern0', x,
                             RooArgList(bern_c0))
    pdf_bern1 = RooBernstein('bern1'+tag, 'bern1', x,
                             RooArgList(bern_c0, bern_c1))
    pdf_bern2 = RooBernstein('bern2'+tag, 'bern2', x,
                             RooArgList(bern_c0, bern_c1, bern_c2))
    pdf_bern3 = RooBernstein('bern3'+tag, 'bern3', x,
                             RooArgList(bern_c0, bern_c1, bern_c2, bern_c3))
    pdf_bern4 = RooBernstein('bern4'+tag, 'bern4', x,
                             RooArgList(bern_c0, bern_c1, bern_c2, bern_c3,bern_c4))
    pdf_bern5 = RooBernstein('bern5'+tag, 'bern5', x,
                             RooArgList(bern_c0, bern_c1, bern_c2, bern_c3,bern_c4, bern_c5))

#    if blinded: pdf_bern4.selectNormalizationRange(RooFit.CutRange("left,right"))

    # -----------------------------------------------------------------------------
    # GAUS law
    gauss_mu = RooRealVar('gauss_mu'+tag, 'gauss_mu', 30., 0, 50.)
    gauss_sigma = RooRealVar('gauss_sigma'+tag, 'gaus_sigma', 15, 10., 50)
    pdf_gauss = RooGaussian('gauss'+tag, 'gauss', x , gauss_mu, gauss_sigma)

    # -----------------------------------------------------------------------------
    # POW law
    # str expressions of formulae
    formula_pow1 = 'TMath::Power(@0, @1)'
    formula_pow2 = '(1.-@1)*TMath::Power(@0,@2) + @1*TMath::Power(@0,@3)'
    formula_pow3 = '(1.-@1-@2)*TMath::Power(@0,@3) + @1*TMath::Power(@0,@4) + @2*TMath::Power(@0,@5)'

    # Variables
    pow_frac1 = RooRealVar('frac1', 'frac1', 0.01, 0., 1.)
    pow_frac2 = RooRealVar('frac2', 'frac2', 0.01, 0., 1.)
    pow_p1 = RooRealVar('p1', 'p1', -2.555, -10., 0.)
    pow_p2 = RooRealVar('p2', 'p2', -8., -10., 0.)
    pow_p3 = RooRealVar('p3', 'p3', -10., -10., 0.)

    # Power Law PDFs
    pdf_pow1 = RooGenericPdf('pow1', 'pow1', formula_pow1,
                             RooArgList(x, pow_p1))
    pdf_pow2 = RooGenericPdf('pow2', 'pow2', formula_pow2,
                             RooArgList(x, pow_frac1, pow_p1, pow_p2))
    pdf_pow3 = RooGenericPdf('pow3', 'pow3', formula_pow3,
                            RooArgList(x, pow_frac1, pow_frac2, pow_p1, pow_p2, pow_p3))

    # -----------------------------------------------------------------------------
    # EXP law
    exp_p1 = RooRealVar('exp_p1'+tag, 'exp_p1', -0.1, -10, 0)
    exp_p2 = RooRealVar('exp_p2', 'exp_p2', -1e-2, -10, 0)
    exp_p3 = RooRealVar('exp_p3', 'exp_p3', -1e-3, -10, 0)
    exp_c1 = RooRealVar('exp_c1', 'exp_c1', 0., 1.)
    exp_c2 = RooRealVar('exp_c2', 'exp_c2', 0., 1.)
    exp_c3 = RooRealVar('exp_c3', 'exp_c3', 0., 1.)

    pdf_exp1 = RooExponential('exp1'+tag+'_bkg', 'exp1', x, exp_p1)
    pdf_single_exp2 = RooExponential('single_exp2', 'single_exp2', x, exp_p2)
    pdf_single_exp3 = RooExponential('single_exp3', 'single_exp3', x, exp_p3)

    pdf_exp3 = RooAddPdf('exp3', 'exp3',
                         RooArgList(pdf_exp1, pdf_single_exp2, pdf_single_exp3),
                         RooArgList(exp_c1, exp_c2, exp_c3))

    pdf_exp1_conv_gauss = RooFFTConvPdf('exp1_conv_gauss', 'exp1 (X) gauss', x, pdf_exp1, pdf_gauss)

#   Set #bins to be used for FFT sampling to 10000
    x.setBins(10000, "cache");
    # Construct landau (x) gauss
#    model = RooFFTConvPdf ("bxg", "bernstein (X) gauss", x, pdf_bern5, pdf_gauss);
#    model = RooFFTConvPdf ("bxg", "bernstein (X) gauss", x, pdf_bern1, pdf_gauss);

    if tag=='_VBFcat' or tag=='_VBFcatlow' or tag=='_Zinvcat':
        model = RooFFTConvPdf ('bxg'+tag+'_bkg', "bernstein (X) gauss", x, pdf_bern2, pdf_gauss);
    elif tag=='_GFcat':
        model = RooFFTConvPdf ('bxg'+tag+'_bkg', "bernstein (X) gauss", x, pdf_bern3, pdf_gauss);
    elif tag=='_Wcat' or tag=='_Zcat':
        model = pdf_exp1
#    model = RooFFTConvPdf ("bxg", "bernstein (X) gauss", x, pdf_exp3, pdf_gauss);
#    model = RooFFTConvPdf ("bxg", "bernstein (X) gauss", x, pdf_pow1, pdf_gauss);
#    model = pdf_gauss
#    model = pdf_bern2
#    model = pdf_bern
#    model = pdf_pow3
#    model = pdf_exp1
#    model = pdf_exp1_conv_gauss

    # -----------------------------------------------------------------------------

    if blinded: model.fitTo(blindedData,RooFit.Minimizer("Minuit2"),RooFit.Strategy(2),RooFit.Range("full"))
    else: model.fitTo(data,RooFit.Minimizer("Minuit2"),RooFit.Strategy(2),RooFit.Range("full"))

    # -----------------------------------------------------------------------------

    norm_SR = data_full.Integral(data_full.FindBin(xlowRange), data_full.FindBin(xhighRange))
    BKG_norm = RooRealVar(model.GetName()+ "_norm", model.GetName()+ "_norm", norm_SR, 0.5*norm_SR, 2*norm_SR)

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Make plot out of the frame

    # Here we will plot the results
    canvas = TCanvas("canvas", "canvas", 800, 800)
    #canvas.Divide(2, 1)

    plotFrameWithNormRange = x.frame(RooFit.Title("mH_"+tag+"_"+mesonCat+"_"+str(year)+" (isoCR) "))

    # Plot only the blinded data, and then plot the PDF over the full range as well as both sidebands

    if blinded:
        blindedData.plotOn(plotFrameWithNormRange)
        model.plotOn(plotFrameWithNormRange, RooFit.LineColor(4), RooFit.Range("left"), RooFit.NormRange("left,right"))
        model.plotOn(plotFrameWithNormRange, RooFit.LineColor(3), RooFit.Range("right"), RooFit.NormRange("left,right"))
#        model.plotOn(plotFrameWithNormRange, RooFit.LineColor(2), RooFit.Range("full"), RooFit.NormRange("left,right"), RooFit.LineStyle(10))
    else:
        data.plotOn(plotFrameWithNormRange)
        model.plotOn(plotFrameWithNormRange, RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(2), RooFit.LineStyle(10))
#        model.plotOn(plotFrameWithNormRange, RooFit.Components("exp3"), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kBlue)) ;
#        model.plotOn(plotFrameWithNormRange, RooFit.Components("bern4"), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kBlue)) ;
#        model.plotOn(plotFrameWithNormRange, RooFit.Components("pow3"), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kBlue)) ;
#        model.plotOn(plotFrameWithNormRange, RooFit.Components("gauss"), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kGreen)) ;
#        model.plotOn(plotFrameWithNormRange, RooFit.LineColor(4), RooFit.Range("full"), RooFit.NormRange("full"))
#        model.plotOn(plotFrameWithNormRange, RooFit.LineColor(3), RooFit.Range("full"), RooFit.NormRange("full"))

    model.paramOn(plotFrameWithNormRange, RooFit.Layout(0.75,0.99,0.95))
    plotFrameWithNormRange.getAttText().SetTextSize(0.02);

    plotFrameWithNormRange.Draw()
    hresid = plotFrameWithNormRange.residHist()
##    RooHist *hresid = plotFrameWithNormRange->residHist();

    canvas.Draw()

    if blinded: canvas.SaveAs("WS/bkg_"+tag+"_"+mesonCat+"_"+str(year)+"REDUCED.png")
    else: canvas.SaveAs("WS/bkg_"+tag+"_"+mesonCat+"_"+str(year)+".png")

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Create workspace, import data and model

    # Create a empty workspace
    w = RooWorkspace("w", "workspace")

    # Import model and all its components into the workspace
    getattr(w,'import')(model)

    # Import model_norm
    getattr(w,'import')(BKG_norm)
    print("integral BKG",BKG_norm.Print())

    # Import data into the workspace
    getattr(w,'import')(data)

    # Print workspace contents
    w.Print()

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Save workspace in file

    # Save the workspace into a ROOT file
    w.writeToFile("WS/Bkg"+tag+"_"+mesonCat+"_"+str(year)+"_workspace.root")

def makePlot():

    f = TFile("Signal_Wcat__RhoCat_2018_workspace.root")

    # Retrieve workspace from file
    w = f.Get("w")
    w.Print()

    x = w.var("mh")
    model = w.pdf("crystal_ball")
    data = w.data("datahist")

    model.fitTo(data)

    frame = x.frame()
    canvas = TCanvas("canvas", "canvas", 800, 800)
    model.plotOn(frame)
    frame.Draw();
    canvas.Draw()
    canvas.SaveAs("bkgTEST.png")

if __name__ == "__main__":

    fitSig('_Wcat','_RhoCat',2018)
    fitBkg('_Wcat','_RhoCat',2018)

    fitSig('_Zcat','_RhoCat',2018)
    fitBkg('_Zcat','_RhoCat',2018)

    fitSig('_Wcat','_PhiCat',2018)
    fitBkg('_Wcat','_PhiCat',2018)

    fitSig('_Zcat','_PhiCat',2018)
    fitBkg('_Zcat','_PhiCat',2018)

    fitSig('_VBFcatlow','_RhoCat',2018)
    fitBkg('_VBFcatlow','_RhoCat',2018)

    fitSig('_VBFcatlow','_PhiCat',2018)
    fitBkg('_VBFcatlow','_PhiCat',2018)

    fitSig('_VBFcat','_PhiCat',2018)
    fitBkg('_VBFcat','_PhiCat',2018)

    fitSig('_VBFcat','_RhoCat',2018)
    fitBkg('_VBFcat','_RhoCat',2018)

    fitSig('_Zinvcat','_PhiCat',2018)
    fitSig('_Zinvcat','_RhoCat',2018)

    fitBkg('_Zinvcat','_PhiCat',2018)
    fitBkg('_Zinvcat','_RhoCat',2018)

    fitSig('_GFcat','_PhiCat',2018)
    fitBkg('_GFcat','_PhiCat',2018)

    fitSig('_GFcat','_RhoCat',2018)
    fitBkg('_GFcat','_RhoCat',2018)

    exit()

#    fitBkg('_VBFcat','_PhiCat',2018)
#    makePlot()
