from ROOT import *

from prepareFits import *

gROOT.SetBatch()
gSystem.Load("libHiggsAnalysisCombinedLimit.so")

xlowRange = 100.
xhighRange = 170.

#Define the observable --------------------------
x = ROOT.RooRealVar("mH", "mH", xlowRange, xhighRange, "GeV")
x.setRange("full", xlowRange, xhighRange)
#x.setRange("left", xlowRange, 115)
#x.setRange("right", 135, xhighRange)
x.setBins(int(xhighRange - xlowRange))

multicanvas = TCanvas("canvas", "canvas", 1600, 1600)
multicanvas.Divide(4,3)

canvasCorr = TCanvas("canvas2", "canvas2", 2000, 2000)                                                                                                                                   
canvasCorr.Divide(2,2)

canvasFits = TCanvas("canvas3", "canvas3", 2000, 2000)                                                                                                                                   
canvasFits.Divide(2,1)

def  checkBias(tag , mesonCat, year):

    # -----------------------------------------------------------------------------
    # SIGNAL law
    # USED in AN
    cb_mu = RooRealVar('cb_mu'+mesonCat+tag, 'cb_mu', 125., 125-10. , 125+10.)
    cb_sigma = RooRealVar('cb_sigma'+mesonCat+tag, 'cb_sigma', 0., 3.)
    cb_alphaL = RooRealVar('cb_alphaL'+mesonCat+tag, 'cb_alphaL', 0., 5.)
    cb_alphaR = RooRealVar('cb_alphaR'+mesonCat+tag, 'cb_alphaR', 0., 5.)
    cb_nL = RooRealVar('cb_nL'+mesonCat+tag, 'cb_nL', 0., 30.)
    cb_nR = RooRealVar('cb_nR'+mesonCat+tag, 'cb_nR', 0., 5.)

    pdf_crystalball = RooDoubleCBFast('crystal_ball'+mesonCat+tag+'_'+"ggH", 'crystal_ball', x, cb_mu, cb_sigma, cb_alphaL, cb_nL, cb_alphaR, cb_nR)

    # -----------------------------------------------------------------------------
    # BERN law

    bern_c0 = RooRealVar('bern_c0'+mesonCat+tag, 'bern_c0', 0.2, 0., 0.5)
    bern_c1 = RooRealVar('bern_c1'+mesonCat+tag, 'bern_c1', 0.1, 0., 1.)
    bern_c2 = RooRealVar('bern_c2'+mesonCat+tag, 'bern_c2', 0.1, 0., 1.)
    bern_c3 = RooRealVar('bern_c3'+mesonCat+tag, 'bern_c3', 0.01, 0., 0.1) # limit this for the GF                                                                                                                     
    bern_c4 = RooRealVar('bern_c4'+mesonCat+tag, 'bern_c4', 0.5, 0., 5.)
    bern_c5 = RooRealVar('bern_c5'+mesonCat+tag, 'bern_c5', 1e-2, 0., 0.1)

    pdf_bern0 = RooBernstein('bern0'+mesonCat+tag, 'bern0', x,
                             RooArgList(bern_c0))
    pdf_bern1 = RooBernstein('bern1'+mesonCat+tag, 'bern1', x,
                             RooArgList(bern_c0, bern_c1))
    pdf_bern2 = RooBernstein('bern2'+mesonCat+tag, 'bern2', x,
                             RooArgList(bern_c0, bern_c1, bern_c2))
    pdf_bern3 = RooBernstein('bern3'+mesonCat+tag, 'bern3', x,
                             RooArgList(bern_c0, bern_c1, bern_c2, bern_c3))
    pdf_bern4 = RooBernstein('bern4'+mesonCat+tag, 'bern4', x,
                             RooArgList(bern_c0, bern_c1, bern_c2, bern_c3, bern_c4))
    pdf_bern5 = RooBernstein('bern5'+mesonCat+tag, 'bern5', x,
                             RooArgList(bern_c0, bern_c1, bern_c2, bern_c3, bern_c4, bern_c5))

    #    if blinded: pdf_bern4.selectNormalizationRange(RooFit.CutRange("left,right"))
    # ---------------------------------------------------------------------------------
    # chebychev law

    chebychev_c0 = RooRealVar('chebychev_c0'+mesonCat+tag, 'chebychev_c0', 1.0, -2, 10.)
    chebychev_c1 = RooRealVar('chebychev_c1'+mesonCat+tag, 'chebychev_c1', 0.4, -1., 1.)
    chebychev_c2 = RooRealVar('chebychev_c2'+mesonCat+tag, 'chebychev_c2', 0.01, -0.1, 0.1)
    chebychev_c3 = RooRealVar('chebychev_c3'+mesonCat+tag, 'chebychev_c3', 0., -1., 1.) # limit this for the GF

    pdf_chebychev2 = RooChebychev("chebychev2"+mesonCat+tag, "chebychev2",x,
                                  RooArgList(chebychev_c0,chebychev_c1,chebychev_c2))

    pdf_chebychev3 = RooChebychev("chebychev3"+mesonCat+tag,"chebychev3",x,
                                  RooArgList(chebychev_c0,chebychev_c1,chebychev_c2,chebychev_c3))

    #
    # ---------------------------------------------------------------------------------
    #

    print('-------------------------------')
    print('HELLO check bias for ', mesonCat)
    print('-------------------------------')

    ## BKG
    dataBKG_full = getHisto(43, 250, 0., 250., True, tag, mesonCat, False, -1 )
    dataBKG = RooDataHist('datahist'+mesonCat+tag, 'data', RooArgList(x), dataBKG_full)

    N = dataBKG_full.Integral(dataBKG_full.FindBin(xlowRange), dataBKG_full.FindBin(xhighRange))
    Nbkg = ROOT.RooRealVar("Nbkg", "N of bkg events",N)    
    print('-------------------------------')
    print('INTEGRAL BKG : ',N)
    print('-------------------------------')    

    pdf1=pdf_chebychev3
    pdf2=pdf_bern3
    
    if tag=='_GFcat':
        fitresult1 = pdf1.fitTo(dataBKG,RooFit.Minimizer("Minuit2"),RooFit.Strategy(2),RooFit.Range("full"),RooFit.Save(kTRUE),RooFit.Minos(kTRUE))
        fitresult2 = pdf2.fitTo(dataBKG,RooFit.Minimizer("Minuit2"),RooFit.Strategy(2),RooFit.Range("full"),RooFit.Save(kTRUE),RooFit.Minos(kTRUE))       

    myFrame=x.frame()
    dataBKG.plotOn(myFrame)
    pdf1.plotOn(myFrame, RooFit.LineColor(ROOT.kGreen+1), RooFit.LineStyle(10))    
    pdf2.plotOn(myFrame, RooFit.LineColor(kMagenta), RooFit.LineStyle(10))
    multicanvas.cd(1)
    myFrame.Draw()
#    myFrame.Print('v')

    ROOT.gPad.SetLeftMargin(0.30)
    hcorr1 = fitresult1.correlationHist()
    hcorr2 = fitresult2.correlationHist()    
    canvasCorr.cd(1)
    hcorr1.Draw('colz')
    canvasCorr.cd(2)
    hcorr2.Draw('colz')    
    canvasCorr.SaveAs("~/public_html/Hrare/BIAS/corrSig.png")   

    name1 = pdf1.GetName()+"_Norm[mH]_Range[fit_nll_"+pdf1.GetName()+"_"+dataBKG.GetName()+"]_NormRange[fit_nll_"+pdf1.GetName()+"_"+dataBKG.GetName()+"]"
    name2 = pdf2.GetName()+"_Norm[mH]_Range[fit_nll_"+pdf2.GetName()+"_"+dataBKG.GetName()+"]_NormRange[fit_nll_"+pdf2.GetName()+"_"+dataBKG.GetName()+"]"    
    chi2_1 = myFrame.chiSquare(name1,"h_"+dataBKG.GetName(),fitresult1.floatParsFinal().getSize())
    chi2_2 = myFrame.chiSquare(name2,"h_"+dataBKG.GetName(),fitresult2.floatParsFinal().getSize())
        
    print('--------------------')
    print(pdf1.GetName(),"    chi2/ndof=",round(chi2_2,2)," ndof",fitresult1.floatParsFinal().getSize())
    print(pdf2.GetName(),"    chi2/ndof=",round(chi2_1,2)," ndof",fitresult2.floatParsFinal().getSize())
    print('--------------------')

    ##SIG
    dataSIG_full = getHisto(43 , 250, 0., 250., True, tag, mesonCat, True, "ggH")
    Ns = dataSIG_full.Integral(dataSIG_full.FindBin(xlowRange), dataSIG_full.FindBin(xhighRange)) * 1./1000
    print('-------------------------------')
    print('INTEGRAL SIGNALG : ',Ns)
    print('-------------------------------')
    dataSIG = RooDataHist('datahist'+tag+'_'+"ggH", 'data', RooArgList(x), dataSIG_full)
    fitresultsSig = pdf_crystalball.fitTo(dataSIG,RooFit.Minimizer("Minuit2"),RooFit.Strategy(2),RooFit.Range("full"),RooFit.Save(kTRUE),RooFit.Minos(kTRUE))


    # Plot
    canvasFits.cd(1)
    plotFrameSig = x.frame()
    dataSIG.plotOn(plotFrameSig)
    pdf_crystalball.plotOn(plotFrameSig, RooFit.LineColor(2), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineStyle(10))
    pdf_crystalball.paramOn(plotFrameSig, RooFit.Layout(0.65,0.99,0.75)) #double xmin, double xmax=0.99, double ymin=0.95
    plotFrameSig.getAttText().SetTextSize(0.02);
    plotFrameSig.Draw()
    canvasFits.cd(2)
#    hresid = plotFrameSig.residHist()
#    hresid.Draw()
    hpull = plotFrameSig.pullHist()
    hpull.Draw()
    canvasFits.SaveAs("~/public_html/Hrare/BIAS/fitSig.png")

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(1)
    hcorr = fitresultsSig.correlationHist()
    cor = fitresultsSig.correlationMatrix()
    cov = fitresultsSig.covarianceMatrix()

    canvasCorr.cd(3)
    ROOT.gPad.SetLeftMargin(0.30)
    hcorr.GetYaxis().SetTitleOffset(1.4)
    hcorr.Draw('colz')
    canvasCorr .SaveAs("~/public_html/Hrare/BIAS/corrSig.png")   
    
    Ns_ = ROOT.RooRealVar("Ns_","xsec*lumi*BRmeson",Ns)
    B_R_ = ROOT.RooRealVar("B_R_","branching_ratio",10.,9.,11.)
    Nsig_ggH = ROOT.RooFormulaVar("Nsig_ggH","@0*@1",ROOT.RooArgList(Ns_,B_R_))
    print('-------------------------------')
    print('INTEGRAL SIGNALG injected : ',Nsig_ggH.evaluate())
    print('INTEGRAL BKG injected : ',Nbkg)    
#    print('-------------------------------')
    argListN = ROOT.RooArgList(Nsig_ggH,Nbkg)
#    print('-------------------------------')
    print "arg list N size = ",argListN.getSize()
    print('-------------------------------')

    genPDF = ROOT.RooAddPdf("totPDF_chebychev",pdf1.GetName(),ROOT.RooArgList(pdf_crystalball, pdf1), ROOT.RooArgList(Nsig_ggH,Nbkg))
    fitPDF = ROOT.RooAddPdf("totPDF_bernstein",pdf2.GetName(),ROOT.RooArgList(pdf_crystalball, pdf2), ROOT.RooArgList(Nsig_ggH,Nbkg))    

    multicanvas.Draw()
    multicanvas.SaveAs("~/public_html/Hrare/BIAS/toy1.png")   
    
    #
    # -------------------------------------------------------------------------------------
    #

    mcstudy = RooMCStudy(
        genPDF,
        RooArgSet(x),
        RooFit.Binned(1),
        RooFit.Silence(1),
        RooFit.FitModel(fitPDF),
        RooFit.Extended(1),
        RooFit.FitOptions(RooFit.Save(kTRUE),
                          RooFit.Minimizer("Minuit2"),RooFit.Strategy(2),
                          RooFit.Range("full"),
                          RooFit.PrintEvalErrors(0))) # -1 switches off printing.

    print('RooMCStudy set up')
    #Generate and fit 100 experiments of 40k events each
    mcstudy.generateAndFit(100,40000,kTRUE)
    print('gen and fit done')

#    multicanvas.cd(2)
#    leg2 = ROOT.TLegend(0.12,0.75,0.44,0.97) #left positioning
#    leg2.SetHeader(" ")
#    leg2.AddEntry(0,"gen pdf: "+genPDF.GetTitle(),"")
#    leg2.AddEntry(0,"fit pdf: "+fitPDF.GetTitle(),"")

    multicanvas.cd(5)
    param1 = mcstudy.plotParam(cb_mu)
    param1.Draw()

    multicanvas.cd(6)
    frame2 = mcstudy.plotParam(cb_sigma)
    frame2.Draw()

    multicanvas.cd(7)
    frame3 = mcstudy.plotParam(cb_alphaL)
    frame3.Draw()

    multicanvas.cd(8)
    frame4 = mcstudy.plotParam(cb_alphaR)
    frame4.Draw()

    multicanvas.cd(9)
    frame5 = mcstudy.plotParam(cb_nL)
    frame5.Draw()

    multicanvas.cd(10)
    frame6 = mcstudy.plotParam(cb_nR)
    frame6.Draw()

    ##
    
    genData0 = mcstudy.genData(0) 
    genData0.Print()
    multicanvas.cd(2)
    frame9=x.frame()
    genData0.plotOn(frame9)
    mcstudy.fitResult(0).Print("v");
    frame9.Draw()

    multicanvas.cd(3)
    frame7=x.frame()
    genData25 = mcstudy.genData(25) 
    genData25.Print()
    genData25.plotOn(frame7)
    mcstudy.fitResult(25).Print("v");
    frame7.Draw()

    multicanvas.cd(4)
    frame8=x.frame()
    genData47 = mcstudy.genData(47) 
    genData47.Print()
    genData47.plotOn(frame8)
    mcstudy.fitResult(47).Print("v");
    frame8.Draw()    

    multicanvas.Draw()
    multicanvas.SaveAs("~/public_html/Hrare/BIAS/toy1.png")       

    multicanvas.cd(12)
    BRpull_frame = mcstudy.plotPull(B_R_, ROOT.RooFit.Range(-5,5), ROOT.RooFit.Bins(100), ROOT.RooFit.FitGauss(1))
    BRpull_frame.SetTitle("")
    BRpull_frame.SetTitleOffset(1.5,"y")
    BRpull_frame.SetXTitle("Gen: "+genPDF.GetTitle()+" Vs Fit: "+fitPDF.GetTitle())
    BRpull_frame.SetMaximum(1.1*BRpull_frame.GetMaximum())
    BRpull_frame.Draw()

    multicanvas.cd(11)    
    BRpar_frame = mcstudy.plotParam(B_R_, ROOT.RooFit.Range(8,12), ROOT.RooFit.Bins(40), ROOT.RooFit.FitGauss(1))
    BRpar_frame.SetTitle("")
    BRpar_frame.SetTitleOffset(1.5,"y")
    BRpar_frame.SetXTitle("Gen: "+genPDF.GetTitle()+" Vs Fit: "+fitPDF.GetTitle())
    BRpar_frame.SetMaximum(1.1*BRpar_frame.GetMaximum())
    BRpar_frame.Draw()
    
    multicanvas.Draw()

    multicanvas.SaveAs("~/public_html/Hrare/BIAS/toy1.png")       

if __name__ == "__main__":

    checkBias('_GFcat','_RhoCat',2018)

    
