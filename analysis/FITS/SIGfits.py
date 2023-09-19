from ROOT import *

from prepareFits import *

gROOT.SetBatch()
gSystem.Load("libHiggsAnalysisCombinedLimit.so")

blinded=False

doMultiPdf=True
doCR=False
#histoEnum = 43
histoEnum = 44
#histoEnum = 4

workspaceName = 'WS_AUG31'
if histoEnum == 43: workspaceName = 'WSmva_AUG31'
if histoEnum == 44: workspaceName = 'WSmvaLowBin_AUG31'
if histoEnum == 43 and doCR: workspaceName = 'WSmvaCR'

if histoEnum == 44: binBDT="_bin1"

## for GF
xlowRange = 110.
xhighRange = 160.

## for Vcat
#xlowRange = 100.
#xhighRange = 150.

## VBF VBFcat
#xlowRange = 100.
#xhighRange = 170.

x = RooRealVar('mh', 'm_{#gamma,meson}', xlowRange, xhighRange)

x.setRange("full", xlowRange, xhighRange)
x.setRange("left", xlowRange, 115)
x.setRange("right", 135, xhighRange)

def  fitSig(tag , mesonCat, year):

    if mesonCat == '_RhoCat': MVAbin = MVAbinRho
    if mesonCat == '_PhiCat': MVAbin = MVAbinPhi
    if mesonCat == '_K0StarCat': MVAbin = MVAbinRho

    # Create a empty workspace
    w = RooWorkspace("w", "workspace")

    mcSig = ["WH", "WHl", "ZH", "ZHl", "TTH", "VBFH", "ZinvH", "ggH"]
    for sig in mcSig:

        # 1l-2l: ZH, WH, TTH
        if (sig == "VBFH" and tag == "_Vcat"): continue
        if (sig == "ggH" and tag == "_Vcat"): continue
        if (sig == "ZinvH" and tag == "_Vcat"): continue
        if (sig == "WHl" and tag == "_Vcat"): continue
        # 1l: ZHl, WH, TTH
        if (sig == "ZH" and tag == "_Wcat"): continue
        if (sig == "VBFH" and tag == "_Wcat"): continue
        if (sig == "ggH" and tag == "_Wcat"): continue
        if (sig == "ZinvH" and tag == "_Wcat"): continue
        # MET: WHl, ZinvH, (ZHl is very minor)
        if (sig == "ZHl" and tag == "_Zinvcat"): continue
        if (sig == "ZH" and tag == "_Zinvcat"): continue
        if (sig == "WH" and tag == "_Zinvcat"): continue
        if (sig == "VBFH" and tag == "_Zinvcat"): continue
        if (sig == "ggH" and tag == "_Zinvcat"): continue
        if (sig == "TTH" and tag == "_Zinvcat"): continue
        ##
        if (sig == "WH" and tag == "_Zcat"): continue
        if (sig == "VBFH" and tag == "_Zcat"): continue
        if (sig == "ggH" and tag == "_Zcat"): continue
        if (sig == "ZinvH" and tag == "_Zcat"): continue
        if (sig == "TTH" and tag == "_Zcat"): continue
        ## ggh cat: GF and VBF (in principle can all the other 3 categories as well)
        if (sig == "ZH" and tag == "_GFcat"): continue
        if (sig == "WH" and tag == "_GFcat"): continue
        if (sig == "ZHl" and tag == "_GFcat"): continue
        if (sig == "WHl" and tag == "_GFcat"): continue
        if (sig == "ZinvH" and tag == "_GFcat"): continue
        if (sig == "TTH" and tag == "_GFcat"): continue
        ##
        if (sig == "WH" and tag == "_VBFcat"): continue
        if (sig == "ZH" and tag == "_VBFcat"): continue
        if (sig == "WHl" and tag == "_VBFcat"): continue
        if (sig == "ZHl" and tag == "_VBFcat"): continue
        if (sig == "ZinvH" and tag == "_VBFcat"): continue
        if (sig == "ggH" and tag == "_VBFcat"): continue
        if (sig == "TTH" and tag == "_VBFcat"): continue
        ##
        if (sig == "WH" and tag == "_VBFcatlow"): continue
        if (sig == "ZH" and tag == "_VBFcatlow"): continue
        if (sig == "WHl" and tag == "_VBFcatlow"): continue
        if (sig == "ZHl" and tag == "_VBFcatlow"): continue
        if (sig == "ZinvH" and tag == "_VBFcatlow"): continue
        if (sig == "ggH" and tag == "_VBFcatlow"): continue
        if (sig == "TTH" and tag == "_VBFcatlow"): continue

        print(tag, ' ', sig)

        foo = -1 * histoEnum if doCR else histoEnum
        data_full = getHisto(foo, int(xhighRange - xlowRange)*10, xlowRange, xhighRange, True, tag, mesonCat, True, sig)

        data = RooDataHist('datahist'+tag+'_'+sig, 'data', RooArgList(x), data_full)

        # -----------------------------------------------------------------------------

        cb_mu = RooRealVar('cb_mu'+binBDT+mesonCat+tag+'_'+sig, 'cb_mu', 125., 125-10. , 125+10.)
        cb_sigma = RooRealVar('cb_sigma'+binBDT+mesonCat+tag+'_'+sig, 'cb_sigma', 0., 3.)
        cb_alphaL = RooRealVar('cb_alphaL'+binBDT+mesonCat+tag+'_'+sig, 'cb_alphaL', 0., 5.)
        cb_alphaR = RooRealVar('cb_alphaR'+binBDT+mesonCat+tag+'_'+sig, 'cb_alphaR', 0., 5.)
        cb_nL = RooRealVar('cb_nL'+binBDT+mesonCat+tag+'_'+sig, 'cb_nL', 0., 5.)
        cb_nR = RooRealVar('cb_nR'+binBDT+mesonCat+tag+'_'+sig, 'cb_nR', 0., 5.)

        pdf_crystalball = RooDoubleCBFast('crystal_ball'+binBDT+mesonCat+tag+'_'+sig, 'crystal_ball', x, cb_mu, cb_sigma, cb_alphaL, cb_nL, cb_alphaR, cb_nR)
        model = pdf_crystalball

        # -----------------------------------------------------------------------------

        model.fitTo(data,RooFit.Minimizer("Minuit2"),RooFit.Strategy(2),RooFit.Range("full"))

        # Here we will plot the results
        canvas = TCanvas("canvas", "canvas", 800, 800)

        titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- "
        if histoEnum == 43:
            if doCR:
                titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- MVA<"+str(MVAbin[tag])
                if tag == '_GFcat': titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- MVA<"+str(MVAbin[tag])
                if tag == '_Zinvcat': titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- MVA<"+str(MVAbin[tag])
            else:
                titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- MVA>"+str(MVAbin[tag])
                if tag == '_GFcat': titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- MVA>"+str(MVAbin[tag])
                if tag == '_Zinvcat': titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- MVA>"+str(MVAbin[tag])
        elif histoEnum == 44:
            titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- "+str(MVAbinLow[tag])+"<MVA<"+str(MVAbin[tag])
            if tag == '_GFcat': titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- "+str(MVAbinLow[tag])+"<MVA<"+str(MVAbin[tag])
            if tag == '_Zinvcat': titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- "+str(MVAbinLow[tag])+"<MVA<"+str(MVAbin[tag])

        plotFrameWithNormRange = x.frame(RooFit.Title(titleSTR))

        # Plot only the blinded data, and then plot the PDF over the full range as well as both sidebands
        data.plotOn(plotFrameWithNormRange)
        model.plotOn(plotFrameWithNormRange, RooFit.LineColor(2), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineStyle(10))
        model.paramOn(plotFrameWithNormRange, RooFit.Layout(0.65,0.99,0.75))
        plotFrameWithNormRange.getAttText().SetTextSize(0.02);

        plotFrameWithNormRange.Draw()

        canvas.Draw()

        if histoEnum == 43:
            if doCR: canvas.SaveAs(workspaceName+"/signal_"+sig+"_"+mesonCat+tag+"_"+str(year)+"_lowMVA.png")
            else: canvas.SaveAs(workspaceName+"/signal_"+sig+"_"+mesonCat+tag+"_"+str(year)+"_withMVA.png")
        elif histoEnum == 44: canvas.SaveAs(workspaceName+"/signal_"+sig+"_"+mesonCat+tag+"_"+str(year)+"_withMVAbin1.png")
        else : canvas.SaveAs(workspaceName+"/signal_"+sig+"_"+mesonCat+tag+"_"+str(year)+".png")

        # -----------------------------------------------------------------------------

        binLow = data_full.GetBin(1) #contains the first bin with low-edge
        binUp = data_full.GetBin(int(xhighRange-xlowRange)*10)  # second to last bin contains the upper-edge

        norm_SR = data_full.Integral(binLow, binUp)
        if doCR:
            Sig_norm = RooRealVar(model.GetName()+ "_normCR", model.GetName()+ "_normCR", norm_SR) # no range means contants
        else:
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
        print('INSIDE fitSig: integral signal = ',Sig_norm.Print())

        # Import data into the workspace
        getattr(w,'import')(data)

        # Print workspace contents
        w.Print()

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Save workspace in file

    if not blinded:
        if histoEnum == 43:
            if doCR: w.writeToFile(workspaceName+"/Signal"+tag+"_"+mesonCat+"_"+str(year)+"_workspace.root")
            else: w.writeToFile(workspaceName+"/Signal"+tag+"_"+mesonCat+"_"+str(year)+"_workspace.root")
        elif histoEnum == 44: w.writeToFile(workspaceName+"/Signal"+tag+"_"+mesonCat+"_"+str(year)+"_bin1_workspace.root")
        else: w.writeToFile(workspaceName+"/Signal"+tag+"_"+mesonCat+"_"+str(year)+"_workspace.root")

def  fitBkg(tag , mesonCat, year):

    if mesonCat == '_RhoCat': MVAbin = MVAbinRho
    if mesonCat == '_PhiCat': MVAbin = MVAbinPhi
    if mesonCat == '_K0StarCat': MVAbin = MVAbinRho

    x.setBins(int(115-xlowRange), "left")
    x.setBins(int(xhighRange-135), "right")

    foo = -1 * histoEnum if doCR else histoEnum
    data_full = getHisto(foo, int(xhighRange - xlowRange), xlowRange, xhighRange, True, tag, mesonCat, False, -1 )
    data = RooDataHist('datahist'+binBDT+mesonCat+tag, 'data', RooArgList(x), data_full)
    blindedData = data.reduce(RooFit.CutRange("left,right"))

    data_reduced_Manual = data_full.Clone()
    for k in range(115, 135):
        data_reduced_Manual.SetBinContent(data_reduced_Manual.FindBin(k),0)
        data_reduced_Manual.SetBinError(data_reduced_Manual.FindBin(k),0.)
    data_reduced = RooDataHist('datahistReduce'+binBDT+mesonCat+tag, 'dataReduced', RooArgList(x), data_reduced_Manual)

    # -----------------------------------------------------------------------------
    # BERN law
    bern_c0 = RooRealVar('bern_c0'+binBDT+mesonCat+tag, 'bern_c0', 0.5, 0., 1.)
    bern_c1 = RooRealVar('bern_c1'+binBDT+mesonCat+tag, 'bern_c1', 0.1, 0., 1.)
    bern_c2 = RooRealVar('bern_c2'+binBDT+mesonCat+tag, 'bern_c2', 0.1, 0., 1.)
    bern_c3 = RooRealVar('bern_c3'+binBDT+mesonCat+tag, 'bern_c3', 0.01, 0., 0.1) # limit this for the GF
    bern_c4 = RooRealVar('bern_c4'+binBDT+mesonCat+tag, 'bern_c4', 0.5, 0., 5.)
    bern_c5 = RooRealVar('bern_c5'+binBDT+mesonCat+tag, 'bern_c5', 1e-2, 0., 0.1)

    pdf_bern0 = RooBernstein('bern0'+binBDT+mesonCat+tag, 'bern0', x,
                             RooArgList(bern_c0))
    pdf_bern1 = RooBernstein('bern1'+binBDT+mesonCat+tag, 'bern1', x,
                             RooArgList(bern_c0, bern_c1))
    pdf_bern2 = RooBernstein('bern2'+binBDT+mesonCat+tag, 'bern2', x,
                             RooArgList(bern_c0, bern_c1, bern_c2))
    pdf_bern3 = RooBernstein('bern3'+binBDT+mesonCat+tag, 'bern3', x,
                             RooArgList(bern_c0, bern_c1, bern_c2, bern_c3))
    pdf_bern4 = RooBernstein('bern4'+binBDT+mesonCat+tag, 'bern4', x,
                             RooArgList(bern_c0, bern_c1, bern_c2, bern_c3, bern_c4))
    pdf_bern5 = RooBernstein('bern5'+binBDT+mesonCat+tag, 'bern5', x,
                             RooArgList(bern_c0, bern_c1, bern_c2, bern_c3, bern_c4, bern_c5))

#    if blinded: pdf_bern4.selectNormalizationRange(RooFit.CutRange("left,right"))

    # -----------------------------------------------------------------------------
    # chebychev law
    chebychev_c0 = RooRealVar('chebychev_c0'+binBDT+mesonCat+tag, 'chebychev_c0', 1.08, -1.1, 10.)
    chebychev_c1 = RooRealVar('chebychev_c1'+binBDT+mesonCat+tag, 'chebychev_c1', 0.4, -1., 1.)
    chebychev_c2 = RooRealVar('chebychev_c2'+binBDT+mesonCat+tag, 'chebychev_c2', 0.01, -0.1, 0.1)
    chebychev_c3 = RooRealVar('chebychev_c3'+binBDT+mesonCat+tag, 'chebychev_c3', 0., -1., 1.) # limit this for the GF

    pdf_chebychev1 = RooChebychev("chebychev1"+binBDT+mesonCat+tag, "chebychev1",x,
                                  RooArgList(chebychev_c0,chebychev_c1))

    pdf_chebychev2 = RooChebychev("chebychev2"+binBDT+mesonCat+tag, "chebychev2",x,
                                  RooArgList(chebychev_c0,chebychev_c1,chebychev_c2))

    pdf_chebychev3 = RooChebychev("chebychev3"+binBDT+mesonCat+tag,"chebychev3",x,
                                  RooArgList(chebychev_c0,chebychev_c1,chebychev_c2,chebychev_c3))

    # -----------------------------------------------------------------------------
    # GAUS law
    gauss_mu = RooRealVar('gauss_mu'+binBDT+mesonCat+tag, 'gauss_mu', 10., 0, 30.)
    gauss_sigma = RooRealVar('gauss_sigma'+binBDT+mesonCat+tag, 'gaus_sigma', 10, 2, 30) #starting value, min max
    pdf_gauss = RooGaussian('gauss'+binBDT+mesonCat+tag, 'gauss', x , gauss_mu, gauss_sigma)

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
    exp_p1 = RooRealVar('exp_p1'+binBDT+mesonCat+tag, 'exp_p1', -0.1, -10, 0)
    exp_p2 = RooRealVar('exp_p2'+binBDT+mesonCat+tag, 'exp_p2', -1e-2, -10, 0)
    exp_p3 = RooRealVar('exp_p3', 'exp_p3', -1e-3, -10, 0)
    exp_c1 = RooRealVar('exp_c1', 'exp_c1', 0., 1.)
    exp_c2 = RooRealVar('exp_c2', 'exp_c2', 0., 1.)
    exp_c3 = RooRealVar('exp_c3', 'exp_c3', 0., 1.)

    pdf_exp1 = RooExponential('exp1'+binBDT+mesonCat+tag, 'exp1', x, exp_p1)
    pdf_single_exp2 = RooExponential('single_exp2', 'single_exp2', x, exp_p2)
    pdf_single_exp3 = RooExponential('single_exp3', 'single_exp3', x, exp_p3)

    pdf_exp2 = RooAddPdf('exp2', 'exp2',
                         RooArgList(pdf_exp1, pdf_single_exp2),
                         RooArgList(exp_c1, exp_c2))

    pdf_exp3 = RooAddPdf('exp3', 'exp3',
                         RooArgList(pdf_exp1, pdf_single_exp2, pdf_single_exp3),
                         RooArgList(exp_c1, exp_c2, exp_c3))

    pdf_exp1_conv_gauss = RooFFTConvPdf('exp1_conv_gauss', 'exp1 (X) gauss', x, pdf_exp1, pdf_gauss)

#   Set #bins to be used for FFT sampling to 10000
#    x.setBins(10000, "cache");
    # Construct landau (x) gauss
#    model = RooFFTConvPdf ("bxg", "bernstein (X) gauss", x, pdf_bern5, pdf_gauss);
#    model = RooFFTConvPdf ("bxg", "bernstein (X) gauss", x, pdf_bern1, pdf_gauss);

    storedPdfs = RooArgList("store_"+binBDT+mesonCat+tag)

    if tag=='_VBFcatlow':
#        model = RooFFTConvPdf ('bxg'+mesonCat+tag+'_bkg', "bernstein (X) gauss", x, pdf_bern2, pdf_gauss);
        model = pdf_bern1
        model2 = pdf_chebychev1
    elif tag=='_VBFcat':
        if histoEnum == 4:
            model = pdf_bern2
            model2 = pdf_chebychev2
        else:
            model = pdf_bern1
            model2 = pdf_chebychev1
    elif tag=='_GFcat':
#        model = RooFFTConvPdf ('bxg'+mesonCat+tag+'_bkg', "bernstein (X) gauss", x, pdf_bern3, pdf_gauss);
        model = pdf_bern3
        model2 = pdf_chebychev3
    elif tag=='_Wcat' or tag=='_Zcat' or tag=='_Zinvcat' or tag=='_Vcat':
        model = pdf_exp1
        model2 = pdf_chebychev2
#    model = RooFFTConvPdf ("bxg", "bernstein (X) gauss", x, pdf_exp3, pdf_gauss);
#    model = RooFFTConvPdf ("bxg", "bernstein (X) gauss", x, pdf_pow1, pdf_gauss);
#    model = pdf_gauss
#    model = pdf_bern2
#    model = pdf_bern
#    model = pdf_pow3
#    model = pdf_exp1
#    model = pdf_exp1_conv_gauss

    if blinded: model.fitTo(blindedData,RooFit.Minimizer("Minuit2"),RooFit.Strategy(2),RooFit.Range("full"))
    else: fitresults = model.fitTo(data,RooFit.Minimizer("Minuit2"),RooFit.Strategy(2),RooFit.Range("full"),RooFit.Save(kTRUE))

    if doMultiPdf:
        storedPdfs.add(model)
        if blinded: model2.fitTo(blindedData,RooFit.Minimizer("Minuit2"),RooFit.Strategy(2),RooFit.Range("full"))
        else: fitresults2 = model2.fitTo(data,RooFit.Minimizer("Minuit2"),RooFit.Strategy(2),RooFit.Range("full"),RooFit.Save(kTRUE))
        storedPdfs.add(model2)  # extra PDF

    # -----------------------------------------------------------------------------


    binLow = data_full.GetBin(1) #contains the first bin with low-edge
    binUp = data_full.GetBin(int(xhighRange-xlowRange))  # second to last bin contains the upper-edge

    norm_range = data_full.Integral( binLow, binUp )
    print("--------------------------")
    print("NORM BKG",norm_range)
    print(' binX1 = ',data_full.GetXaxis().GetBinLowEdge(binLow)," - ",data_full.GetXaxis().GetBinUpEdge(binLow))
    print(' binX2 = ',data_full.GetXaxis().GetBinLowEdge(binUp)," - ",data_full.GetXaxis().GetBinUpEdge(binUp))
    print("--------------------------")

    if doCR:
        BKG_norm = RooRealVar(model.GetName()+ "_normCR", model.GetName()+ "_normCR", norm_range, 0.5*norm_range, 2*norm_range)
    else:
        if doMultiPdf:
            BKG_norm = RooRealVar("multipdf"+binBDT+mesonCat+tag+"_bkg"+"_norm",model.GetName()+"_norm", norm_range, 0.5*norm_range, 2*norm_range)
        else:
            BKG_norm = RooRealVar(model.GetName()+ "_norm", model.GetName()+ "_norm", norm_range, 0.5*norm_range, 2*norm_range)

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Make plot out of the frame

    # Here we will plot the results
    canvas = TCanvas("canvas", "canvas", 800, 800)
    #canvas.Divide(2, 1)

    titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- "
    if histoEnum == 43:
        if doCR:
            titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- MVA<"+str(MVAbin[tag])
            if tag == '_GFcat': titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- MVA<"+str(MVAbin[tag])
            if tag == '_Zinvcat': titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- MVA<"+str(MVAbin[tag])
        else:
            titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- MVA>"+str(MVAbin[tag])
            if tag == '_GFcat': titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- MVA>"+str(MVAbin[tag])
            if tag == '_Zinvcat': titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- MVA>"+str(MVAbin[tag])
    elif histoEnum == 44:
        titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- "+str(MVAbinLow[tag])+"<MVA<"+str(MVAbin[tag])
        if tag == '_GFcat': titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- "+str(MVAbinLow[tag])+"<MVA<"+str(MVAbin[tag])
        if tag == '_Zinvcat': titleSTR = "mH"+mesonCat+tag+"_"+str(year)+" -- "+str(MVAbinLow[tag])+"<MVA<"+str(MVAbin[tag])

#    plotFrameWithNormRange = x.frame(RooFit.Title("mH_"+tag+"_"+mesonCat+"_"+str(year)+" -- MVA>0.3"))
    plotFrameWithNormRange = x.frame(RooFit.Title(titleSTR))

    # Plot only the blinded data, and then plot the PDF over the full range as well as both sidebands

    if blinded:
        data_reduced.plotOn(plotFrameWithNormRange,RooFit.MarkerColor(kWhite),RooFit.LineColor(kWhite))

        model.plotOn(plotFrameWithNormRange, RooFit.Components(model.GetName()), RooFit.LineColor(kRed), RooFit.Range("left"), RooFit.NormRange("left,right"))
        model.plotOn(plotFrameWithNormRange, RooFit.Components(model.GetName()), RooFit.LineColor(kRed), RooFit.Range("right"), RooFit.NormRange("left,right"))
        model2.plotOn(plotFrameWithNormRange, RooFit.Components(model2.GetName()), RooFit.LineColor(kBlue), RooFit.Range("left"), RooFit.NormRange("left,right"))
        model2.plotOn(plotFrameWithNormRange, RooFit.Components(model2.GetName()), RooFit.LineColor(kBlue), RooFit.Range("right"), RooFit.NormRange("left,right"))

        data_reduced.plotOn(plotFrameWithNormRange, RooFit.Binning("left"))
        data_reduced.plotOn(plotFrameWithNormRange, RooFit.Binning("right"))

    else:
        data.plotOn(plotFrameWithNormRange)
#        model.plotOn(plotFrameWithNormRange, RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(2), RooFit.LineStyle(10))
        model.plotOn(plotFrameWithNormRange, RooFit.Components(model.GetName()), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kRed)) ;
        model2.plotOn(plotFrameWithNormRange, RooFit.Components(model2.GetName()), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kBlue)) ;
#        model.plotOn(plotFrameWithNormRange, RooFit.Components("exp3"), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kBlue)) ;
#        model.plotOn(plotFrameWithNormRange, RooFit.Components("bern4"), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kBlue)) ;
#        model.plotOn(plotFrameWithNormRange, RooFit.Components("pow3"), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kBlue)) ;
#        model.plotOn(plotFrameWithNormRange, RooFit.Components("gauss"), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kGreen)) ;
#        model.plotOn(plotFrameWithNormRange, RooFit.LineColor(4), RooFit.Range("full"), RooFit.NormRange("full"))
#        model.plotOn(plotFrameWithNormRange, RooFit.LineColor(3), RooFit.Range("full"), RooFit.NormRange("full"))
        name1 = model.GetName()+"_Norm[mh]_Comp["+model.GetName()+"]_Range[full]_NormRange[full]"
        name2 = model2.GetName()+"_Norm[mh]_Comp["+model2.GetName()+"]_Range[full]_NormRange[full]"
        chi2_1 = plotFrameWithNormRange.chiSquare(name1,"h_"+data.GetName(),fitresults.floatParsFinal().getSize())
        chi2_2 = plotFrameWithNormRange.chiSquare(name2,"h_"+data.GetName(),fitresults2.floatParsFinal().getSize())
#        plotFrameWithNormRange.Print("v")
        print('--------------------')
        print(model2.GetName(),"    chi2/ndof=",round(chi2_2,2)," ndof",fitresults2.floatParsFinal().getSize())
        print(model.GetName(),"    chi2/ndof=",round(chi2_1,2)," ndof",fitresults.floatParsFinal().getSize())
        print('--------------------')

#    model.paramOn(plotFrameWithNormRange, RooFit.Layout(0.6,0.99,0.95))
#    plotFrameWithNormRange.getAttText().SetTextSize(0.02);

    plotFrameWithNormRange.Draw()
#    hresid = plotFrameWithNormRange.residHist()


    offsetY = 0.60*data_full.GetMaximum()
    if tag == '_GFcat': offsetY = 0
    if (tag == '_VBFcat' or tag == '_VBFcatlow') and mesonCat == '_RhoCat' and histoEnum == 4: offsetY = 0
    if (tag == '_VBFcat' or tag == '_VBFcatlow') and mesonCat == '_PhiCat' and histoEnum == 4: offsetY = 0.80*data_full.GetMaximum()
    if (tag == '_VBFcat') and (histoEnum == 43 or histoEnum == 44): offsetY = 0.90*data_full.GetMaximum()
    latex = TLatex()
    latex.SetTextColor(kRed)
    latex.SetTextSize(0.04)
    latex.DrawLatex(111 ,offsetY + 0.10*data_full.GetMaximum(), model.GetName())
    latex.SetTextColor(kBlue)
    latex.DrawLatex(111 ,offsetY + 0.20*data_full.GetMaximum(), model2.GetName())

    canvas.Draw()

    if blinded:
        if histoEnum == 43:
            if doCR: canvas.SaveAs(workspaceName+"/bkg_"+mesonCat+tag+"_"+str(year)+"_lowMVA.png")
            else: canvas.SaveAs(workspaceName+"/bkg_"+mesonCat+tag+"_"+str(year)+"_REDUCED_withMVA.png")
        elif histoEnum == 44: canvas.SaveAs(workspaceName+"/bkg_"+mesonCat+tag+"_"+str(year)+"_REDUCED_withMVA_bin1.png")
        else: canvas.SaveAs(workspaceName+"/bkg_"+mesonCat+tag+"_"+str(year)+"REDUCED.png")
    else:
        if histoEnum == 43:
            if doCR: canvas.SaveAs(workspaceName+"/bkg_"+mesonCat+tag+"_"+str(year)+"_lowMVA.png")
            else: canvas.SaveAs(workspaceName+"/bkg_"+mesonCat+tag+"_"+str(year)+"_withMVA.png")
        elif histoEnum == 44: canvas.SaveAs(workspaceName+"/bkg_"+mesonCat+tag+"_"+str(year)+"_withMVA_bin1.png")
        else: canvas.SaveAs(workspaceName+"/bkg_"+mesonCat+tag+"_"+str(year)+".png")

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Create workspace, import data and model

    # Create a empty workspace
    w = RooWorkspace("w", "workspace")

    if doMultiPdf:
        pdf_cat = RooCategory("pdfindex"+binBDT+mesonCat+tag,"pdfindex"+binBDT+mesonCat+tag)
        pdf_bkg = RooMultiPdf("multipdf"+binBDT+mesonCat+tag+"_bkg","multipdf",pdf_cat,storedPdfs)
        getattr(w,'import')(pdf_bkg)
    else:
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

    if not blinded:
        if histoEnum == 43:
            if doCR: w.writeToFile(workspaceName+"/Bkg"+tag+"_"+mesonCat+"_"+str(year)+"_workspace.root")
            else: w.writeToFile(workspaceName+"/Bkg"+tag+"_"+mesonCat+"_"+str(year)+"_workspace.root")
        elif histoEnum == 44: w.writeToFile(workspaceName+"/Bkg"+tag+"_"+mesonCat+"_"+str(year)+"_bin1_workspace.root")
        else: w.writeToFile(workspaceName+"/Bkg"+tag+"_"+mesonCat+"_"+str(year)+"_workspace.root")

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

#    with open("example_VBFcatlow_RhoCat_2018_mHwindow.txt", "w") as f:
#        f.write("Events used in VBFcatlow RhoCat 2018 \n")
#        f.write("run:luminosityBlock:event \n")

    if not blinded: fitSig('_GFcat','_RhoCat',2018)
    fitBkg('_GFcat','_RhoCat',2018)

    if not blinded: fitSig('_GFcat','_PhiCat',2018)
    fitBkg('_GFcat','_PhiCat',2018)

    if not blinded: fitSig('_GFcat','_K0StarCat',2018)
    fitBkg('_GFcat','_K0StarCat',2018)

    '''

    if not blinded: fitSig('_VBFcatlow','_RhoCat',2018)
    fitBkg('_VBFcatlow','_RhoCat',2018)

    if not blinded: fitSig('_VBFcatlow','_PhiCat',2018)
    fitBkg('_VBFcatlow','_PhiCat',2018)

    if not blinded: fitSig('_VBFcatlow','_K0StarCat',2018)
    fitBkg('_VBFcatlow','_K0StarCat',2018)

    if not blinded: fitSig('_VBFcat','_RhoCat','Run2')
    fitBkg('_VBFcat','_RhoCat','Run2')

    if not blinded: fitSig('_VBFcat','_PhiCat','Run2')
    fitBkg('_VBFcat','_PhiCat','Run2')

    if not blinded: fitSig('_VBFcat','_K0StarCat','Run2')
    fitBkg('_VBFcat','_K0StarCat','Run2')

    if histoEnum == 4:

        if not blinded: fitSig('_Vcat','_RhoCat','Run2')
        fitBkg('_Vcat','_RhoCat','Run2')

        if not blinded: fitSig('_Vcat','_PhiCat','Run2')
        fitBkg('_Vcat','_PhiCat','Run2')

        if not blinded: fitSig('_Vcat','_K0StarCat','Run2')
        fitBkg('_Vcat','_K0StarCat','Run2')
    '''

    exit()

    fitSig('_Wcat','_RhoCat','Run2')
    fitBkg('_Wcat','_RhoCat','Run2')

    fitSig('_Wcat','_PhiCat','Run2')
    fitBkg('_Wcat','_PhiCat','Run2')

    fitSig('_Zcat','_RhoCat','Run2')
    fitBkg('_Zcat','_RhoCat','Run2')

    fitSig('_Zcat','_PhiCat','Run2')
    fitBkg('_Zcat','_PhiCat','Run2')

    fitSig('_Zinvcat','_RhoCat',2018)
    fitBkg('_Zinvcat','_RhoCat',2018)

    fitSig('_Zinvcat','_PhiCat',2018)
    fitBkg('_Zinvcat','_PhiCat',2018)

    exit()

#    makePlot()
