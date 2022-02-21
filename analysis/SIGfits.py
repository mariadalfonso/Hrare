from ROOT import *

from prepareFits import getHisto

blinded=False

def  fitSig(tag , mesonCat, year ):

    data_full = getHisto(4, 200*10, 0. , 200., True, tag, mesonCat,True)
    
    xlowRange = 125-25.
    xhighRange = 125+25.
    x = RooRealVar('x', 'm_{#gamma,meson}', xlowRange, xhighRange)

    x.setRange("full", xlowRange, xhighRange)
    data = RooDataHist('datahist', 'data', RooArgList(x), data_full)
    
    xlow = 115
    xhigh = 135

    cb_mu = RooRealVar('cb_mu', 'cb_mu', xlow, 1e-4, xhigh)
    cb_sigmaL = RooRealVar('cb_sigmaL', 'cb_sigmaL',0., 2.)
    cb_sigmaR = RooRealVar('cb_sigmaR', 'cb_sigmaR',0., 2.)
    cb_alphaL = RooRealVar('cb_alphaL', 'cb_alphaL',0., 5.)
    cb_alphaR = RooRealVar('cb_alphaR', 'cb_alphaR',0., 5.)
    cb_nL = RooRealVar('cb_nL', 'cb_nL',0., 5.)
    cb_nR = RooRealVar('cb_nR', 'cb_nR',0., 5.)

    pdf_crystalball = RooCrystalBall('crystal_ball', 'crystal_ball', x, cb_mu, cb_sigmaL, cb_sigmaR, cb_alphaL, cb_nL, cb_alphaR, cb_nR)
    
    pdf_crystalball.fitTo(data)

    # Here we will plot the results
    canvas = TCanvas("canvas", "canvas", 800, 800)
    
    plotFrameWithNormRange = x.frame(RooFit.Title("mH_"+tag+"_"+mesonCat+"_"+str(year)))
    
    # Plot only the blinded data, and then plot the PDF over the full range as well as both sidebands
    data.plotOn(plotFrameWithNormRange)
    pdf_crystalball.plotOn(plotFrameWithNormRange, RooFit.LineColor(2), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineStyle(10))
    pdf_crystalball.paramOn(plotFrameWithNormRange, RooFit.Layout(0.65,0.99,0.75))
    plotFrameWithNormRange.getAttText().SetTextSize(0.02);

    plotFrameWithNormRange.Draw()
    
    canvas.Draw()
 
    canvas.SaveAs("signal_"+tag+"_"+mesonCat+"_"+str(year)+".png")


def  fitBkg(tag , mesonCat, year ):

    data_full = getHisto(4, 250, 0. , 250., True, tag, mesonCat, False)

    xlowRange = 60.
    xhighRange = 250.
    x = RooRealVar('x', 'm_{#gamma,meson}', xlowRange, xhighRange)

    x.setRange("sr", 90, 250)
    x.setRange("left", 90, 115)
    x.setRange("right", 135, 200)

    data = RooDataHist('datahist', 'data', RooArgList(x), data_full)

    blindedData = data.reduce(RooFit.CutRange("left,right"))

    #-------------------
    # BERN law
    bern_c0 = RooRealVar('bern_c0', 'bern_c0', 0.2, 0., 0.5)
    bern_c1 = RooRealVar('bern_c1', 'bern_c1', 0.1, 0., 1.)
    bern_c2 = RooRealVar('bern_c2', 'bern_c2', 0.1, 0., 1.)
    bern_c3 = RooRealVar('bern_c3', 'bern_c3', 0.5, 0., 5.)
    bern_c4 = RooRealVar('bern_c4', 'bern_c4', 0.5, 0., 5.)
    bern_c5 = RooRealVar('bern_c5', 'bern_c5', 1e-2, 0., 0.1)

    pdf_bern0 = RooBernstein('bern0', 'bern0', x,
                             RooArgList(bern_c0))
    pdf_bern1 = RooBernstein('bern1', 'bern1', x,
                             RooArgList(bern_c0, bern_c1))
    pdf_bern2 = RooBernstein('bern2', 'bern2', x,
                             RooArgList(bern_c0, bern_c1, bern_c2))
    pdf_bern3 = RooBernstein('bern3', 'bern3', x,
                             RooArgList(bern_c0, bern_c1, bern_c2, bern_c3))
    pdf_bern4 = RooBernstein('bern4', 'bern4', x,
                             RooArgList(bern_c0, bern_c1, bern_c2, bern_c3,
                                        bern_c4))
    pdf_bern5 = RooBernstein('bern5', 'bern5', x,
                             RooArgList(bern_c0, bern_c1, bern_c2, bern_c3,
                                        bern_c4, bern_c5))

    #-------------------
    # GAUSS law
    gauss_mu = RooRealVar('gauss_mu', 'gauss_mu', 30., 0, 50.)
    gauss_sigma = RooRealVar('gauss_sigma', 'gaus_sigma', 15, 10., 50)
    pdf_gauss = RooGaussian('gauss', 'gauss', x , gauss_mu, gauss_sigma)

    #-------------------
    # POWER law
    formula_pow1 = 'TMath::Power(@0, @1)'
    formula_pow2 = '(1.-@1)*TMath::Power(@0,@2) + @1*TMath::Power(@0,@3)'
    formula_pow3 = '(1.-@1-@2)*TMath::Power(@0,@3) + @1*TMath::Power(@0,@4) + @2*TMath::Power(@0,@5)'

    pow_frac1 = RooRealVar('frac1', 'frac1', 0.01, 0., 1.)
    pow_frac2 = RooRealVar('frac2', 'frac2', 0.01, 0., 1.)
    pow_p1 = RooRealVar('p1', 'p1', -2.555, -10., 0.)
    pow_p2 = RooRealVar('p2', 'p2', -8., -10., 0.)
    pow_p3 = RooRealVar('p3', 'p3', -10., -10., 0.)

    pdf_pow1 = RooGenericPdf('pow1', 'pow1', formula_pow1,
                             RooArgList(x, pow_p1))
    pdf_pow2 = RooGenericPdf('pow2', 'pow2', formula_pow2,
                             RooArgList(x, pow_frac1, pow_p1, pow_p2))
    pdf_pow3 = RooGenericPdf('pow3', 'pow3', formula_pow3,
                            RooArgList(x, pow_frac1, pow_frac2, pow_p1, pow_p2, pow_p3))

    #-------------------
    # EXP law
    exp_p1 = RooRealVar('exp_p1', 'exp_p1', -0.1, -10, 0)
    exp_p2 = RooRealVar('exp_p2', 'exp_p2', -1e-2, -10, 0)
    exp_p3 = RooRealVar('exp_p3', 'exp_p3', -1e-3, -10, 0)
    exp_c1 = RooRealVar('exp_c1', 'exp_c1', 0., 1.)
    exp_c2 = RooRealVar('exp_c2', 'exp_c2', 0., 1.)
    exp_c3 = RooRealVar('exp_c3', 'exp_c3', 0., 1.)

    pdf_exp1 = RooExponential('exp1', 'exp1', x, exp_p1)
    pdf_single_exp2 = RooExponential('single_exp2', 'single_exp2', x, exp_p2)
    pdf_single_exp3 = RooExponential('single_exp3', 'single_exp3', x, exp_p3)

    pdf_exp3 = RooAddPdf('exp3', 'exp3',
                         RooArgList(pdf_exp1, pdf_single_exp2, pdf_single_exp3),
                         RooArgList(exp_c1, exp_c2, exp_c3))

    pdf_exp1_conv_gauss = RooFFTConvPdf('exp1_conv_gauss', 'exp1 (X) gauss', x, pdf_exp1, pdf_gauss)

    x.setBins(10000, "cache");
    pdf_final = RooFFTConvPdf ("bxg", "bernstein (X) gauss", x, pdf_bern2, pdf_gauss);
#    pdf_final = RooFFTConvPdf ("bxg", "bernstein (X) gauss", x, pdf_exp3, pdf_gauss);
#    pdf_final = RooFFTConvPdf ("bxg", "bernstein (X) gauss", x, pdf_pow1, pdf_gauss);
#    pdf_final = pdf_gauss
#    pdf_final = pdf_bern2
#    pdf_final = pdf_bern
#    pdf_final = pdf_pow3
#    pdf_final = pdf_exp1
#    pdf_final = pdf_exp1_conv_gauss

    if blinded: pdf_final.fitTo(blindedData,RooFit.Minimizer("Minuit2"),RooFit.Strategy(2),RooFit.Range("sr"))
    else: pdf_final.fitTo(data,RooFit.Minimizer("Minuit2"),RooFit.Strategy(2),RooFit.Range("sr"))

    # --------------------------------
    # Here we will plot the results
    canvas = TCanvas("canvas", "canvas", 800, 800)

    plotFrameWithNormRange = x.frame(RooFit.Title("mH_"+tag+"_"+mesonCat+"_"+str(year)+" (isoCR) "))

    # Plot only the blinded data, and then plot the PDF over the full range as well as both sidebands

    if blinded:
        blindedData.plotOn(plotFrameWithNormRange)
        pdf_final.plotOn(plotFrameWithNormRange, RooFit.LineColor(4), RooFit.Range("left"), RooFit.NormRange("left,right"))
        pdf_final.plotOn(plotFrameWithNormRange, RooFit.LineColor(3), RooFit.Range("right"), RooFit.NormRange("left,right"))
#        pdf_final.plotOn(plotFrameWithNormRange, RooFit.LineColor(2), RooFit.Range("full"), RooFit.NormRange("left,right"), RooFit.LineStyle(10))
    else:
        data.plotOn(plotFrameWithNormRange)
        pdf_final.plotOn(plotFrameWithNormRange, RooFit.Range("sr"), RooFit.NormRange("sr"), RooFit.LineColor(2), RooFit.LineStyle(10))
#        pdf_final.plotOn(plotFrameWithNormRange, RooFit.Components("exp3"), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kBlue)) ;
#        pdf_final.plotOn(plotFrameWithNormRange, RooFit.Components("bern4"), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kBlue)) ;
#        pdf_final.plotOn(plotFrameWithNormRange, RooFit.Components("pow3"), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kBlue)) ;
#        pdf_final.plotOn(plotFrameWithNormRange, RooFit.Components("gauss"), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineColor(kGreen)) ;

    pdf_final.paramOn(plotFrameWithNormRange, RooFit.Layout(0.75,0.99,0.95))
    plotFrameWithNormRange.getAttText().SetTextSize(0.02);

    plotFrameWithNormRange.Draw()
    hresid = plotFrameWithNormRange.residHist()
##    RooHist *hresid = plotFrameWithNormRange->residHist();

    canvas.Draw()

    if blinded: canvas.SaveAs("bkg_"+tag+"_"+mesonCat+"_"+str(year)+"REDUCED.png")
    else: canvas.SaveAs("bkg_"+tag+"_"+mesonCat+"_"+str(year)+".png")

if __name__ == "__main__":

    fitSig('_Wcat','_PhiCat',2018)
    fitSig('_VBFcat','_PhiCat',2018)
    fitSig('_VBFcat','_RhoCat',2018)
    fitBkg('_VBFcat','_RhoCat',2018)
    fitBkg('_VBFcat','_PhiCat',2018)
    fitBkg('_Wcat','_PhiCat',2018)
    fitBkg('_Wcat','_RhoCat',2018)
