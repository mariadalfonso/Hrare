import ROOT

from prepareFits import *

ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
ROOT.gInterpreter.Declare('#include "HiggsAnalysis/CombinedLimit/interface/RooDoubleCBFast.h"')
#ROOT.gInterpreter.Declare('#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"')

#ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
#ROOT.gROOT.ProcessLine('.L RooDoubleCBFast.cc+')

blinded=False

doMultiPdf=False
doCR=False
#histoEnum = 43
#histoEnum = 44
histoEnum = 4

workspaceName = 'WS'
#if histoEnum == 43: workspaceName = 'WSmva_JUNE13'
#if histoEnum == 44: workspaceName = 'WSmvaLowBin_JUNE13'
#if histoEnum == 43: workspaceName = 'WSmva_NOV1_neuIso'
#if histoEnum == 44: workspaceName = 'WSmvaLowBin_NOV1_neuIso'
#if histoEnum == 43 and doCR: workspaceName = 'WSmvaCR'
#if histoEnum == 4: binBDT=""
#if histoEnum == 43: binBDT=""
#if histoEnum == 44: binBDT="_bin1"

## for GF
xlowRange = 50
xhighRange = 250

## for Vcat
#xlowRange = 100.
#xhighRange = 150.

## VBF VBFcat
#xlowRange = 50.
#xhighRange = 150.

x = ROOT.RooRealVar('mh'+'GFcat', 'm_{c,c,J#psi}', xlowRange, xhighRange)
#x = RooRealVar('mh'+'VBFcat', 'm_{#gamma,meson}', xlowRange, xhighRange)
#x = RooRealVar('mh'+'VBFcatlow', 'm_{#gamma,meson}', xlowRange, xhighRange)
#x = RooRealVar('mh'+'Vcat', 'm_{#gamma,meson}', xlowRange, xhighRange)

x.setRange("full", xlowRange, xhighRange)
x.setRange("left", xlowRange, 115)
x.setRange("right", 135, xhighRange)

def  fitSig(tag , year):

    # Create a empty workspace
    w = ROOT.RooWorkspace("w", "workspace")

    print('WORKSPACE DONE')
    
    data_full = getHisto(histoEnum, int(xhighRange - xlowRange), xlowRange, xhighRange, True,tag, True)
#    data_full = getHisto(histoEnum, 100, 0. , 250., True,tag, True)
    print('getHisto  DONE')
    
    data = ROOT.RooDataHist('datahist', 'data', ROOT.RooArgList(x), data_full)
#    data = ROOT.RooDataHist(data_full) 

    # -----------------------------------------------------------------------------
    
    cb_mu = ROOT.RooRealVar('cb_mu'+tag, 'cb_mu', 125., 125-10. , 125+10.)
    cb_sigma = ROOT.RooRealVar('cb_sigma'+tag, 'cb_sigma', 15, 14., 17.)
    cb_alphaL = ROOT.RooRealVar('cb_alphaL'+tag, 'cb_alphaL', 1.1, 1., 1.2)
    cb_alphaR = ROOT.RooRealVar('cb_alphaR'+tag, 'cb_alphaR', 0., 5.)
    cb_nL = ROOT.RooRealVar('cb_nL'+tag, 'cb_nL', 0., 5.)
    cb_nR = ROOT.RooRealVar('cb_nR'+tag, 'cb_nR', 0., 10.)

    pdf_crystalball = ROOT.RooDoubleCBFast('crystal_ball'+tag, 'crystal_ball', x, cb_mu, cb_sigma, cb_alphaL, cb_nL, cb_alphaR, cb_nR)
#    pdf_crystalball = ROOT.RooDoubleCB('crystal_ball'+tag, 'crystal_ball', x, cb_mu, cb_sigma, cb_alphaL, cb_nL, cb_alphaR, cb_nR)    
    model = pdf_crystalball

    # -----------------------------------------------------------------------------

    model.fitTo(data,ROOT.RooFit.Minimizer("Minuit2"),ROOT.RooFit.Strategy(2),ROOT.RooFit.Range("full"))
    
    # Here we will plot the results
    if True:
        # Create canvas with two pads
        canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
        pad1 = ROOT.TPad("pad1", "Top Pad", 0, 0.35, 1, 1)
        pad2 = ROOT.TPad("pad2", "Bottom Pad", 0, 0, 1, 0.35)

        pad1.SetBottomMargin(0.03)
        pad2.SetTopMargin(0.05)
        pad2.SetBottomMargin(0.3)

        pad1.Draw()
        pad2.Draw()

        # ========== TOP PAD ==========
        pad1.cd()


        titleSTR = "mH"+tag+"_"+str(year)+" -- "
        '''
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
        '''

#        titleSTR = ""
        plotFrameWithNormRange = x.frame(ROOT.RooFit.Title(titleSTR))

        # Plot only the blinded data, and then plot the PDF over the full range as well as both sidebands
        data.plotOn(plotFrameWithNormRange, ROOT.RooFit.Name("fit_data"))
        model.plotOn(plotFrameWithNormRange, ROOT.RooFit.LineColor(2), ROOT.RooFit.Range("full"), ROOT.RooFit.NormRange("full"), ROOT.RooFit.LineStyle(10), ROOT.RooFit.Name("fit_curve"))
        model.paramOn(plotFrameWithNormRange, ROOT.RooFit.Layout(0.50,0.99,0.85))
        plotFrameWithNormRange.getAttText().SetTextSize(0.02);

        plotFrameWithNormRange.Draw()

        # ========== BOTTOM PAD (ratio plot) ==========
        pad2.cd()

        # Compute residuals = (data - fit) / fit
        # Get RooHist of data and RooCurve of model
        plotFrameWithNormRange.Print("v")
        dataHist = plotFrameWithNormRange.getHist("fit_data")
        fitCurve = plotFrameWithNormRange.getCurve("fit_curve")

        # Create a new histogram for the ratio
        ratioHist = ROOT.TH1D("ratio", "", dataHist.GetN(), xlowRange, xhighRange)

        for i in range(dataHist.GetN()):
            x_val = dataHist.GetX()[i]
            y_data = dataHist.GetY()[i]
            y_fit = fitCurve.Eval(x_val)

            if y_fit != 0:
                ratio = (y_data - y_fit) / y_fit
            else:
                ratio = 0

            ratioHist.SetBinContent(i + 1, ratio)
            ratioHist.SetBinError(i + 1, dataHist.GetErrorY(i) / y_fit if y_fit != 0 else 0)

        ratioHist.GetYaxis().SetTitle("Ratio")
        ratioHist.GetYaxis().SetTitleSize(0.08)
        ratioHist.GetYaxis().SetLabelSize(0.08)
        ratioHist.GetYaxis().SetTitleOffset(0.4)
        ratioHist.GetXaxis().SetTitle("m_{c,c,J#psi}")
        ratioHist.GetXaxis().SetTitleSize(0.1)
        ratioHist.GetXaxis().SetLabelSize(0.08)
        ratioHist.SetLineColor(ROOT.kBlack)
        ratioHist.SetMarkerStyle(20)
        ratioHist.Draw("EP")

        # Draw horizontal line at 0
        line = ROOT.TLine(xlowRange, 0.0, xhighRange, 0.0)
        line.SetLineColor(ROOT.kRed)
        line.SetLineWidth(2)
        line.SetLineStyle(2)
        line.Draw()

        # ========== Save the canvas ==========
        canvas.Draw()
        canvas.SaveAs(workspaceName+"/signal_"+tag+"_"+str(year)+".png")
        chi2_ndf = plotFrameWithNormRange.chiSquare()
        print("Chi² / ndf =", chi2_ndf)

        # -----------------------------------------------------------------------------

    '''
        if histoEnum == 43:
            if doCR: canvas.SaveAs(workspaceName+"/signal_"+sig+"_"+mesonCat+tag+"_"+str(year)+"_lowMVA.png")
            else: canvas.SaveAs(workspaceName+"/signal_"+sig+"_"+mesonCat+tag+"_"+str(year)+"_withMVA.png")
        elif histoEnum == 44: canvas.SaveAs(workspaceName+"/signal_"+sig+"_"+mesonCat+tag+"_"+str(year)+"_withMVAbin1.png")
        else : canvas.SaveAs(workspaceName+"/signal_"+sig+"_"+mesonCat+tag+"_"+str(year)+".png")
    '''
        
    # -----------------------------------------------------------------------------

    binLow = data_full.GetBin(1) #contains the first bin with low-edge
    binUp = data_full.GetBin(int(xhighRange-xlowRange)*10)  # second to last bin contains the upper-edge

    norm_SR = data_full.Integral(binLow, binUp)

    '''
        # this is to account the h.c.
        if mesonCat == '_K0StarCat': norm_SR = norm_SR * 2
        # this is for the trigger SF
        if tag == '_VBFcat': norm_SR = 0.95 * norm_SR

        if doCR:
            Sig_norm = RooRealVar(model.GetName()+ "_normCR", model.GetName()+ "_normCR", norm_SR) # no range means contants
        else:
            Sig_norm = RooRealVar(model.GetName()+ "_norm", model.GetName()+ "_norm", norm_SR) # no range means contants
        '''

    Sig_norm = ROOT.RooRealVar(model.GetName()+ "_norm", model.GetName()+ "_norm", norm_SR) # no range means contants            

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

    w.writeToFile(workspaceName+"/Signal_workspace.root")
    '''
    if not blinded:
        if histoEnum == 43:
            if doCR: w.writeToFile(workspaceName+"/Signal"+tag+"_"+mesonCat+"_"+str(year)+"_workspace.root")
            else: w.writeToFile(workspaceName+"/Signal"+tag+"_"+mesonCat+"_"+str(year)+"_workspace.root")
        elif histoEnum == 44: w.writeToFile(workspaceName+"/Signal"+tag+"_"+mesonCat+"_"+str(year)+"_bin1_workspace.root")
        else: w.writeToFile(workspaceName+"/Signal"+tag+"_"+mesonCat+"_"+str(year)+"_workspace.root")
    '''

def  fitBkg(tag, year):

    x.setBins(int(115-xlowRange), "left")
    x.setBins(int(xhighRange-135), "right")

    foo = -1 * histoEnum if doCR else histoEnum
#    data_full = getHisto(foo, int(xhighRange - xlowRange), xlowRange, xhighRange, True, tag, False, -1 )
    data_full = getHisto(histoEnum, 100, 0. , 250., True,tag, False)    
    
    data = ROOT.RooDataHist('datahist'+tag, 'data', ROOT.RooArgList(x), data_full)
    blindedData = data.reduce(ROOT.RooFit.CutRange("left,right"))

    data_reduced_Manual = data_full.Clone()
    for k in range(115, 135):
        data_reduced_Manual.SetBinContent(data_reduced_Manual.FindBin(k),0)
        data_reduced_Manual.SetBinError(data_reduced_Manual.FindBin(k),0.)
    data_reduced = ROOT.RooDataHist('datahistReduce'+tag, 'dataReduced', ROOT.RooArgList(x), data_reduced_Manual)

    # -----------------------------------------------------------------------------
    # BERN law
    bern_c0 = ROOT.RooRealVar('bern_c0'+tag, 'bern_c0', 0.5, 0., 1.)
    bern_c1 = ROOT.RooRealVar('bern_c1'+tag, 'bern_c1', 0.1, 0., 1.)
    bern_c2 = ROOT.RooRealVar('bern_c2'+tag, 'bern_c2', 0.1, 0., 1.)
    bern_c3 = ROOT.RooRealVar('bern_c3'+tag, 'bern_c3', 0.01, 0., 0.1) # limit this for the GF
    bern_c4 = ROOT.RooRealVar('bern_c4'+tag, 'bern_c4', 0.5, 0., 5.)
    bern_c5 = ROOT.RooRealVar('bern_c5'+tag, 'bern_c5', 1e-2, 0., 0.1)

    pdf_bern0 = ROOT.RooBernstein('bern0'+tag, 'bern0', x,
                                  ROOT.RooArgList(bern_c0))
    pdf_bern1 = ROOT.RooBernstein('bern1'+tag, 'bern1', x,
                                  ROOT.RooArgList(bern_c0, bern_c1))
    pdf_bern2 = ROOT.RooBernstein('bern2'+tag, 'bern2', x,
                                  ROOT.RooArgList(bern_c0, bern_c1, bern_c2))
    pdf_bern3 = ROOT.RooBernstein('bern3'+tag, 'bern3', x,
                                  ROOT.RooArgList(bern_c0, bern_c1, bern_c2, bern_c3))
    pdf_bern4 = ROOT.RooBernstein('bern4'+tag, 'bern4', x,
                                  ROOT.RooArgList(bern_c0, bern_c1, bern_c2, bern_c3, bern_c4))
    pdf_bern5 = ROOT.RooBernstein('bern5'+tag, 'bern5', x,
                                  ROOT.RooArgList(bern_c0, bern_c1, bern_c2, bern_c3, bern_c4, bern_c5))

#    if blinded: pdf_bern4.selectNormalizationRange(ROOT.RooFit.CutRange("left,right"))

    # -----------------------------------------------------------------------------
    # chebychev law
    chebychev_c0 = ROOT.RooRealVar('chebychev_c0'+tag, 'chebychev_c0', 1.08, -1.1, 10.)
    chebychev_c1 = ROOT.RooRealVar('chebychev_c1'+tag, 'chebychev_c1', 0.4, -1., 1.)
    chebychev_c2 = ROOT.RooRealVar('chebychev_c2'+tag, 'chebychev_c2', 0.01, -0.1, 0.1)
    chebychev_c3 = ROOT.RooRealVar('chebychev_c3'+tag, 'chebychev_c3', 0., -1., 1.) # limit this for the GF

    pdf_chebychev1 = ROOT.RooChebychev("chebychev1"+tag, "chebychev1",x,
                                       ROOT.RooArgList(chebychev_c0,chebychev_c1))

    pdf_chebychev2 = ROOT.RooChebychev("chebychev2"+tag, "chebychev2",x,
                                       ROOT.RooArgList(chebychev_c0,chebychev_c1,chebychev_c2))

    pdf_chebychev3 = ROOT.RooChebychev("chebychev3"+tag,"chebychev3",x,
                                       ROOT.RooArgList(chebychev_c0,chebychev_c1,chebychev_c2,chebychev_c3))

    # -----------------------------------------------------------------------------
    # GAUS law
    gauss_mu = ROOT.RooRealVar('gauss_mu', 'gauss_mu', 120, 100, 140)
    gauss_sigma = ROOT.RooRealVar('gauss_sigma', 'gauss_sigma', 30, 20, 40) #starting value, min max
    pdf_gauss = ROOT.RooGaussian('gauss', 'gauss', x , gauss_mu, gauss_sigma)

    # -----------------------------------------------------------------------------
    # POW law
    # str expressions of formulae
    formula_pow1 = 'TMath::Power(@0, @1)'
    formula_pow2 = '(1.-@1)*TMath::Power(@0,@2) + @1*TMath::Power(@0,@3)'
    formula_pow3 = '(1.-@1-@2)*TMath::Power(@0,@3) + @1*TMath::Power(@0,@4) + @2*TMath::Power(@0,@5)'

    # Variables
    pow_frac1 = ROOT.RooRealVar('frac1', 'frac1', 0.01, 0., 1.)
    pow_frac2 = ROOT.RooRealVar('frac2', 'frac2', 0.01, 0., 1.)
    pow_p1 = ROOT.RooRealVar('p1', 'p1', -2.555, -10., 0.)
    pow_p2 = ROOT.RooRealVar('p2', 'p2', -8., -10., 0.)
    pow_p3 = ROOT.RooRealVar('p3', 'p3', -10., -10., 0.)

    # Power Law PDFs
    pdf_pow1 = ROOT.RooGenericPdf('pow1', 'pow1', formula_pow1,
                                  ROOT.RooArgList(x, pow_p1))
    pdf_pow2 = ROOT.RooGenericPdf('pow2', 'pow2', formula_pow2,
                                  ROOT.RooArgList(x, pow_frac1, pow_p1, pow_p2))
    pdf_pow3 = ROOT.RooGenericPdf('pow3', 'pow3', formula_pow3,
                                  ROOT.RooArgList(x, pow_frac1, pow_frac2, pow_p1, pow_p2, pow_p3))
    
    # -----------------------------------------------------------------------------
    # EXP law
    # these two from Davit
    exp_p1 = ROOT.RooRealVar('exp_p1', 'exp_p1', -0.0207, -0.022, -0.018)
    exp_p2 = ROOT.RooRealVar('exp_p2', 'exp_p2', -1e-2, -10, 10)
    # old param
#    exp_p1 = ROOT.RooRealVar('exp_p1'+tag, 'exp_p1', -0.1, -10, 0)
#    exp_p2 = ROOT.RooRealVar('exp_p2'+tag, 'exp_p2', -1e-2, -10, 0)
    exp_p3 = ROOT.RooRealVar('exp_p3', 'exp_p3', -1e-3, -10, 0)
    exp_c1 = ROOT.RooRealVar('exp_c1', 'exp_c1', 0., 1.)
    exp_c2 = ROOT.RooRealVar('exp_c2', 'exp_c2', 0., 1.)
    exp_c3 = ROOT.RooRealVar('exp_c3', 'exp_c3', 0., 1.)

    pdf_exp1 = ROOT.RooExponential('exp1'+tag, 'exp1', x, exp_p1)
    pdf_single_exp2 = ROOT.RooExponential('single_exp2', 'single_exp2', x, exp_p2)
    pdf_single_exp3 = ROOT.RooExponential('single_exp3', 'single_exp3', x, exp_p3)

    pdf_exp2 = ROOT.RooAddPdf('exp2', 'exp2',
                         ROOT.RooArgList(pdf_exp1, pdf_single_exp2),
                         ROOT.RooArgList(exp_c1, exp_c2))

    pdf_exp3 = ROOT.RooAddPdf('exp3', 'exp3',
                         ROOT.RooArgList(pdf_exp1, pdf_single_exp2, pdf_single_exp3),
                         ROOT.RooArgList(exp_c1, exp_c2, exp_c3))

    pdf_exp1_conv_gauss = ROOT.RooFFTConvPdf('exp1_conv_gauss', 'exp1 (X) gauss', x, pdf_exp1, pdf_gauss)

#   Set #bins to be used for FFT sampling to 10000
#    x.setBins(10000, "cache");
    # Construct landau (x) gauss
#    model = RooFFTConvPdf ("bxg", "bernstein (X) gauss", x, pdf_bern5, pdf_gauss);
#    model = RooFFTConvPdf ("bxg", "bernstein (X) gauss", x, pdf_bern1, pdf_gauss);

    storedPdfs = ROOT.RooArgList("store_"+tag)

    if tag=='_VBFcatlow':
#        model = RooFFTConvPdf ('bxg'+tag+'_bkg', "bernstein (X) gauss", x, pdf_bern2, pdf_gauss);
        model2 = pdf_chebychev1
        if histoEnum == 43 or histoEnum == 44: model = pdf_bern1
        if histoEnum == 4: model = pdf_bern2
    elif tag=='_VBFcat':
        if histoEnum == 4:
            model = pdf_bern2
            model2 = pdf_chebychev1
        else:
            model = pdf_bern1
            model2 = pdf_chebychev1
    elif tag=='_GFcat':
        if histoEnum == 4:
            model = pdf_exp1_conv_gauss
#            model = pdf_bern1
            model2 = pdf_chebychev1
        elif histoEnum == 43 or histoEnum == 44:
            model = pdf_bern2
            model2 = pdf_chebychev1
    elif tag=='_Wcat' or tag=='_Zcat' or tag=='_Zinvcat' or tag=='_Vcat':
        model = pdf_exp1
        model2 = pdf_chebychev1
#    model = RooFFTConvPdf ("bxg", "bernstein (X) gauss", x, pdf_exp3, pdf_gauss);
#    model = RooFFTConvPdf ("bxg", "bernstein (X) gauss", x, pdf_pow1, pdf_gauss);
#    model = pdf_gauss
#    model = pdf_bern2
#    model = pdf_bern
#    model = pdf_pow3
#    model = pdf_exp1
#    model = pdf_exp1_conv_gauss

    if blinded: model.fitTo(blindedData,ROOT.RooFit.Minimizer("Minuit2"),ROOT.RooFit.Strategy(2),ROOT.RooFit.Range("full"))
    else: fitresults = model.fitTo(data,ROOT.RooFit.Minimizer("Minuit2"),ROOT.RooFit.Strategy(2),ROOT.RooFit.Range("full"),ROOT.RooFit.Save(ROOT.kTRUE))

    if doMultiPdf:
        storedPdfs.add(model)
        if blinded: model2.fitTo(blindedData,ROOT.RooFit.Minimizer("Minuit2"),ROOT.RooFit.Strategy(2),ROOT.RooFit.Range("full"))
        else: fitresults2 = model2.fitTo(data,ROOT.RooFit.Minimizer("Minuit2"),ROOT.RooFit.Strategy(2),ROOT.RooFit.Range("full"),ROOT.RooFit.Save(ROOT.kTRUE))
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
        BKG_norm = ROOT.RooRealVar(model.GetName()+ "_normCR", model.GetName()+ "_normCR", norm_range, 0.5*norm_range, 2*norm_range)
    else:
        if doMultiPdf:
            BKG_norm = ROOT.RooRealVar("multipdf"+tag+"_bkg"+"_norm",model.GetName()+"_norm", norm_range, 0.5*norm_range, 2*norm_range)
        else:
            BKG_norm = ROOT.RooRealVar(model.GetName()+ "_norm", model.GetName()+ "_norm", norm_range, 0.5*norm_range, 2*norm_range)

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Make plot out of the frame

    # Here we will plot the results
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
    #canvas.Divide(2, 1)

    titleSTR = "mH"+tag+"_"+str(year)+" -- "
    if histoEnum == 43:
        if doCR:
            titleSTR = "mH"+tag+"_"+str(year)+" -- MVA<"+str(MVAbin[tag])
            if tag == '_GFcat': titleSTR = "mH"+tag+"_"+str(year)+" -- MVA<"+str(MVAbin[tag])
            if tag == '_Zinvcat': titleSTR = "mH"+tag+"_"+str(year)+" -- MVA<"+str(MVAbin[tag])
        else:
            titleSTR = "mH"+tag+"_"+str(year)+" -- MVA>"+str(MVAbin[tag])
            if tag == '_GFcat': titleSTR = "mH"+tag+"_"+str(year)+" -- MVA>"+str(MVAbin[tag])
            if tag == '_Zinvcat': titleSTR = "mH"+tag+"_"+str(year)+" -- MVA>"+str(MVAbin[tag])
    elif histoEnum == 44:
        titleSTR = "mH"+tag+"_"+str(year)+" -- "+str(MVAbinLow[tag])+"<MVA<"+str(MVAbin[tag])
        if tag == '_GFcat': titleSTR = "mH"+tag+"_"+str(year)+" -- "+str(MVAbinLow[tag])+"<MVA<"+str(MVAbin[tag])
        if tag == '_Zinvcat': titleSTR = "mH"+tag+"_"+str(year)+" -- "+str(MVAbinLow[tag])+"<MVA<"+str(MVAbin[tag])

#    plotFrameWithNormRange = x.frame(RooFit.Title("mH_"+tag+"_"+"_"+str(year)+" -- MVA>0.3"))
    plotFrameWithNormRange = x.frame(ROOT.RooFit.Title(titleSTR))

    # Plot only the blinded data, and then plot the PDF over the full range as well as both sidebands

    if blinded:
        data_reduced.plotOn(plotFrameWithNormRange, ROOT.RooFit.MarkerColor(ROOT.kWhite), ROOT.RooFit.LineColor(ROOT.kWhite))

        model.plotOn(plotFrameWithNormRange, ROOT.RooFit.Components(model.GetName()), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Range("left"), ROOT.RooFit.NormRange("left,right"))
        model.plotOn(plotFrameWithNormRange, ROOT.RooFit.Components(model.GetName()), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Range("right"), ROOT.RooFit.NormRange("left,right"))
        model2.plotOn(plotFrameWithNormRange, ROOT.RooFit.Components(model2.GetName()), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Range("left"), ROOT.RooFit.NormRange("left,right"))
        model2.plotOn(plotFrameWithNormRange, ROOT.RooFit.Components(model2.GetName()), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Range("right"), ROOT.RooFit.NormRange("left,right"))

        data_reduced.plotOn(plotFrameWithNormRange, ROOT.RooFit.Binning("left"))
        data_reduced.plotOn(plotFrameWithNormRange, ROOT.RooFit.Binning("right"))

    else:
        data.plotOn(plotFrameWithNormRange)
#        model.plotOn(plotFrameWithNormRange, RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(2), RooFit.LineStyle(10))
        model.plotOn(plotFrameWithNormRange, ROOT.RooFit.Components(model.GetName()), ROOT.RooFit.Range("full"), ROOT.RooFit.NormRange("full"), ROOT.RooFit.LineColor(ROOT.kRed)) ;
        model2.plotOn(plotFrameWithNormRange, ROOT.RooFit.Components(model2.GetName()), ROOT.RooFit.Range("full"), ROOT.RooFit.NormRange("full"), ROOT.RooFit.LineColor(ROOT.kBlue)) ;
#        model.plotOn(plotFrameWithNormRange, RooFit.Components("exp3"), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kBlue)) ;
#        model.plotOn(plotFrameWithNormRange, RooFit.Components("bern4"), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kBlue)) ;
#        model.plotOn(plotFrameWithNormRange, RooFit.Components("pow3"), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kBlue)) ;
#        model.plotOn(plotFrameWithNormRange, RooFit.Components("gauss"), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kGreen)) ;
#        model.plotOn(plotFrameWithNormRange, RooFit.LineColor(4), RooFit.Range("full"), RooFit.NormRange("full"))
#        model.plotOn(plotFrameWithNormRange, RooFit.LineColor(3), RooFit.Range("full"), RooFit.NormRange("full"))
        name1 = model.GetName()+"_Norm[mh]_Comp["+model.GetName()+"]_Range[full]_NormRange[full]"
        name2 = model2.GetName()+"_Norm[mh]_Comp["+model2.GetName()+"]_Range[full]_NormRange[full]"
        chi2_1 = plotFrameWithNormRange.chiSquare(name1,"h_"+data.GetName(),fitresults.floatParsFinal().getSize())
        if doMultiPdf: chi2_2 = plotFrameWithNormRange.chiSquare(name2,"h_"+data.GetName(),fitresults2.floatParsFinal().getSize())
#        plotFrameWithNormRange.Print("v")
        print('--------------------')
        if doMultiPdf: print(model2.GetName(),"    chi2/ndof=",round(chi2_2,2)," ndof",fitresults2.floatParsFinal().getSize())
        print(model.GetName(),"    chi2/ndof=",round(chi2_1,2)," ndof",fitresults.floatParsFinal().getSize())
        print('--------------------')

        if histoEnum == 4: fileToWrite="preselection_"+tag+"_"+"_"+str(year)+".txt"
        if histoEnum == 43: fileToWrite="bin0_"+tag+"_"+"_"+str(year)+".txt"
        if histoEnum == 44: fileToWrite="bin1_"+tag+"_"+"_"+str(year)+".txt"
        with open(fileToWrite, "a") as f:
            str1 = model.GetName()+"    chi2/ndof="+str(round(chi2_1,2))+" ndof"+str(fitresults.floatParsFinal().getSize())+"\n"
            f.write(str1)
            if doMultiPdf: str2 = model2.GetName()+"    chi2/ndof="+str(round(chi2_2,2))+" ndof"+str(fitresults2.floatParsFinal().getSize())+"\n"
            if doMultiPdf: f.write(str2)

#    model.paramOn(plotFrameWithNormRange, RooFit.Layout(0.6,0.99,0.95))
#    plotFrameWithNormRange.getAttText().SetTextSize(0.02);

    plotFrameWithNormRange.Draw()
#    hresid = plotFrameWithNormRange.residHist()


    offsetY = 0.60*data_full.GetMaximum()
    if tag == '_GFcat': offsetY = 0
    if (tag == '_VBFcat' or tag == '_VBFcatlow') and histoEnum == 4: offsetY = 0
    if (tag == '_VBFcat' or tag == '_VBFcatlow') and histoEnum == 4: offsetY = 0.80*data_full.GetMaximum()
    if (tag == '_VBFcat') and (histoEnum == 43 or histoEnum == 44): offsetY = 0.90*data_full.GetMaximum()
    latex = ROOT.TLatex()
    latex.SetTextColor(ROOT.kRed)
    latex.SetTextSize(0.04)
    latex.DrawLatex(111 ,offsetY + 0.10*data_full.GetMaximum(), model.GetName())
    latex.SetTextColor(ROOT.kBlue)
    latex.DrawLatex(111 ,offsetY + 0.20*data_full.GetMaximum(), model2.GetName())

    canvas.Draw()

    if blinded:
        if histoEnum == 43:
            if doCR: canvas.SaveAs(workspaceName+"/bkg_"+tag+"_"+str(year)+"_lowMVA.png")
            else: canvas.SaveAs(workspaceName+"/bkg_"+tag+"_"+str(year)+"_REDUCED_withMVA.png")
        elif histoEnum == 44: canvas.SaveAs(workspaceName+"/bkg_"+tag+"_"+str(year)+"_REDUCED_withMVA_bin1.png")
        else: canvas.SaveAs(workspaceName+"/bkg_"+tag+"_"+str(year)+"REDUCED.png")
    else:
        if histoEnum == 43:
            if doCR: canvas.SaveAs(workspaceName+"/bkg_"+tag+"_"+str(year)+"_lowMVA.png")
            else: canvas.SaveAs(workspaceName+"/bkg_"+tag+"_"+str(year)+"_withMVA.png")
        elif histoEnum == 44: canvas.SaveAs(workspaceName+"/bkg_"+tag+"_"+str(year)+"_withMVA_bin1.png")
        else: canvas.SaveAs(workspaceName+"/bkg_"+tag+"_"+str(year)+".png")

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Create workspace, import data and model

    # Create a empty workspace
    w = ROOT.RooWorkspace("w", "workspace")

    if doMultiPdf:
        pdf_cat = ROOT.RooCategory("pdfindex"+tag,"pdfindex"+tag)
        pdf_bkg = ROOT.RooMultiPdf("multipdf"+tag+"_bkg","multipdf",pdf_cat,storedPdfs)
        getattr(w,'import')(pdf_bkg)
    else:
        # Import model and all its components into the workspace
        getattr(w,'import')(model)

    # Import model_norm
    getattr(w,'import')(BKG_norm)
    print("integral BKG",BKG_norm.Print())

    # Import data into the workspace
    getattr(w,'import')(data)

#    print('del = ',w.TestBit(TObject.kNotDeleted))
#    print('onHeap = ', w.TestBit(TObject.kIsOnHeap))

    # Print workspace contents
    w.Print()

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Save workspace in file

    if not blinded:
        if histoEnum == 43:
            if doCR: w.writeToFile(workspaceName+"/Bkg"+tag+"_"+"_"+str(year)+"_workspace.root")
            else: w.writeToFile(workspaceName+"/Bkg"+tag+"_"+"_"+str(year)+"_workspace.root")
        elif histoEnum == 44: w.writeToFile(workspaceName+"/Bkg"+tag+"_"+str(year)+"_bin1_workspace.root")
        else: w.writeToFile(workspaceName+"/Bkg"+tag+"_"+"_"+str(year)+"_workspace.root")

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
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
    model.plotOn(frame)
    frame.Draw();
    canvas.Draw()
    canvas.SaveAs("bkgTEST.png")

if __name__ == "__main__":

    fitSig('_GFcat',2018)
    fitBkg('_GFcat',2018)
