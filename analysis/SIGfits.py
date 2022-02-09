from ROOT import *

from prepareFits import getHisto


def  fitSig(tag , mesonCat, year ):

    
    data_full = getHisto(4, 200*10, 0. , 200., True, tag, mesonCat)
    
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

    # Left edge of box starts at 75% of Xaxis)
    plotFrameWithNormRange.Draw()
    
    canvas.Draw()
 
    canvas.SaveAs("signal_"+tag+"_"+mesonCat+"_"+str(year)+".png")


if __name__ == "__main__":

    fitSig('_Wcat','_PhiCat',2018)
    fitSig('_VBFcat','_PhiCat',2018)
    fitSig('_VBFcat','_RhoCat',2018)    
