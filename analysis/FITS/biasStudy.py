from ROOT import *

from prepareFits import *

gROOT.SetBatch()
gSystem.Load("libHiggsAnalysisCombinedLimit.so")

xlowRange = 100.
xhighRange = 170.

#Define the observable --------------------------
x = ROOT.RooRealVar("mH", "mH", xlowRange, xhighRange, "GeV")
x.setRange("full", xlowRange, xhighRange)
x.setBins(int(xhighRange - xlowRange))

multicanvas = TCanvas("canvas", "canvas", 1600, 1600)
multicanvas.Divide(4,4)

canvasCorr = TCanvas("canvas2", "canvas2", 2000, 2000)                                                                                                                                   
canvasCorr.Divide(2,2)

canvasFits = TCanvas("canvas3", "canvas3", 2000, 2000)                                                                                                                                   
canvasFits.Divide(2,1)

canvasPull = TCanvas("pull", "pull", 800, 800)

def  checkBias(tag , mesonCat, year, myBRvalue):

    # -----------------------------------------------------------------------------
    # SIGNAL law

    # USED in AN
    cb_mu = RooRealVar('cb_mu'+mesonCat+tag, 'cb_mu', 125., 125-10. , 125+10.)
    cb_sigma = RooRealVar('cb_sigma'+mesonCat+tag, 'cb_sigma', 0., 3.)
    cb_alphaL = RooRealVar('cb_alphaL'+mesonCat+tag, 'cb_alphaL', 0., 5.)
    cb_alphaR = RooRealVar('cb_alphaR'+mesonCat+tag, 'cb_alphaR', 0., 5.)
    cb_nL = RooRealVar('cb_nL'+mesonCat+tag, 'cb_nL', 0., 30.)
    cb_nR = RooRealVar('cb_nR'+mesonCat+tag, 'cb_nR', 0., 15.)

    pdf_crystalball = RooDoubleCBFast('crystal_ball'+mesonCat+tag+'_'+"ggH", 'crystal_ball', x, cb_mu, cb_sigma, cb_alphaL, cb_nL, cb_alphaR, cb_nR)

    cb_mu_2 = RooRealVar('cb_mu_2'+mesonCat+tag, 'cb_mu_2', 125., 125-10. , 125+10.)
    cb_sigma_2 = RooRealVar('cb_sigma_2'+mesonCat+tag, 'cb_sigma_2', 0., 3.)
    cb_alphaL_2 = RooRealVar('cb_alphaL_2'+mesonCat+tag, 'cb_alphaL_2', 0., 5.)
    cb_alphaR_2 = RooRealVar('cb_alphaR_2'+mesonCat+tag, 'cb_alphaR_2', 0., 5.)
    cb_nL_2 = RooRealVar('cb_nL_2'+mesonCat+tag, 'cb_nL_2', 0., 30.)
    cb_nR_2 = RooRealVar('cb_nR_2'+mesonCat+tag, 'cb_nR_2', 0., 15.)
    pdf_crystalball2 = RooDoubleCBFast('crystal_ball2'+mesonCat+tag+'_'+"ggH", 'crystal_ball2', x, cb_mu_2, cb_sigma_2, cb_alphaL_2, cb_nL_2, cb_alphaR_2, cb_nR_2)

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

    chebychev_c0 = RooRealVar('chebychev_c0'+mesonCat+tag, 'chebychev_c0', 0., -2, 10.)
    chebychev_c1 = RooRealVar('chebychev_c1'+mesonCat+tag, 'chebychev_c1', 0.4, -1., 1.)
    chebychev_c2 = RooRealVar('chebychev_c2'+mesonCat+tag, 'chebychev_c2', 0.01, -0.1, 0.1)
    chebychev_c3 = RooRealVar('chebychev_c3'+mesonCat+tag, 'chebychev_c3', 0., -1., 1.) # limit this for the GF
#    chebychev_c0 = RooRealVar('chebychev_c0'+mesonCat+tag, 'chebychev_c0', 1.0, -2, 10.)
#    chebychev_c1 = RooRealVar('chebychev_c1'+mesonCat+tag, 'chebychev_c1', 0.4, -1., 1.)
#    chebychev_c2 = RooRealVar('chebychev_c2'+mesonCat+tag, 'chebychev_c2', 0.01, -0.1, 0.1)
#    chebychev_c3 = RooRealVar('chebychev_c3'+mesonCat+tag, 'chebychev_c3', 0., -1., 1.) # limit this for the GF

    pdf_chebychev1 = RooChebychev("chebychev1"+mesonCat+tag, "chebychev1",x,
                                  RooArgList(chebychev_c0,chebychev_c1))

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
    dataBKG_full = getHisto(43, int(xhighRange - xlowRange), xlowRange , xhighRange , True, tag, mesonCat, False, -1 )
    dataBKG = RooDataHist('datahist'+mesonCat+tag, 'data', RooArgList(x), dataBKG_full)

    N = dataBKG_full.Integral(dataBKG_full.FindBin(xlowRange), dataBKG_full.FindBin(xhighRange))
    Nbkg = ROOT.RooRealVar("Nbkg", "N of bkg events",N)    
    print('-------------------------------')
    print('INTEGRAL BKG : ',N)
    print('-------------------------------')    

    if tag=='_GFcat':
        pdf1 = pdf_chebychev3
        pdf2 = pdf_bern3
    elif tag=='_VBFcat' or tag=='_VBFcatlow':
        pdf1 = pdf_bern1
        pdf2 = pdf_chebychev1

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
    canvasCorr.SaveAs("~/public_html/Hrare/BIAS/corrSig"+mesonCat+tag+".png")

    name1 = pdf1.GetName()+"_Norm[mH]_Range[fit_nll_"+pdf1.GetName()+"_"+dataBKG.GetName()+"]_NormRange[fit_nll_"+pdf1.GetName()+"_"+dataBKG.GetName()+"]"
    name2 = pdf2.GetName()+"_Norm[mH]_Range[fit_nll_"+pdf2.GetName()+"_"+dataBKG.GetName()+"]_NormRange[fit_nll_"+pdf2.GetName()+"_"+dataBKG.GetName()+"]"    
    chi2_1 = myFrame.chiSquare(name1,"h_"+dataBKG.GetName(),fitresult1.floatParsFinal().getSize())
    chi2_2 = myFrame.chiSquare(name2,"h_"+dataBKG.GetName(),fitresult2.floatParsFinal().getSize())
        
    print('--------------------')
    print(pdf1.GetName(),"    chi2/ndof=",round(chi2_2,2)," ndof",fitresult1.floatParsFinal().getSize())
    print(pdf2.GetName(),"    chi2/ndof=",round(chi2_1,2)," ndof",fitresult2.floatParsFinal().getSize())
    print('--------------------')

    ##SIG
    if tag=='_GFcat': name = "ggH"
    elif tag=='_VBFcat' or tag=='_VBFcatlow': name = "VBFH"

    dataSIG_full = getHisto(43 , int(xhighRange - xlowRange), xlowRange, xhighRange, True, tag, mesonCat, True, name)
    Ns = dataSIG_full.Integral(dataSIG_full.FindBin(xlowRange), dataSIG_full.FindBin(xhighRange)) * 1./1000
    print('-------------------------------')
    print('INTEGRAL SIGNALG : ',Ns)
    print('-------------------------------')
    dataSIG = RooDataHist('datahist'+tag+'_'+"ggH", 'data', RooArgList(x), dataSIG_full)
    fitresultsSig = pdf_crystalball.fitTo(dataSIG,RooFit.Minimizer("Minuit2"),RooFit.Strategy(2),RooFit.Range("full"),RooFit.Save(kTRUE),RooFit.Minos(kTRUE))
    fitresultsSig2 = pdf_crystalball2.fitTo(dataSIG,RooFit.Minimizer("Minuit2"),RooFit.Strategy(2),RooFit.Range("full"),RooFit.Save(kTRUE),RooFit.Minos(kTRUE))

    # Plot only the blinded data, and then plot the PDF over the full range as well as both sidebands

    canvasFits.cd(1)
    plotFrameSig = x.frame()
    dataSIG.plotOn(plotFrameSig)
    pdf_crystalball.plotOn(plotFrameSig, RooFit.LineColor(2), RooFit.Range("full"), RooFit.NormRange("full"), RooFit.LineStyle(10))
    pdf_crystalball.paramOn(plotFrameSig, RooFit.Layout(0.65,0.99,0.75)) #double xmin, double xmax=0.99, double ymin=0.95
    plotFrameSig.getAttText().SetTextSize(0.02)
    plotFrameSig.Draw()
    canvasFits.cd(2)
#    hresid = plotFrameSig.residHist()
#    hresid.Draw()
    hpull = plotFrameSig.pullHist()
    hpull.Draw()
    canvasFits.SaveAs("~/public_html/Hrare/BIAS/fitSig"+mesonCat+tag+".png")

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
    
    Ns_ = ROOT.RooRealVar("Ns_","xsec*lumi*BRmeson",Ns) # intial is 1
    B_R_ = ROOT.RooRealVar("B_R_","branching_ratio",myBRvalue,-3.,3.) #
    if tag=='_VBFcat' or tag=='_VBFcatlow':
        B_R_ = ROOT.RooRealVar("B_R_","branching_ratio",myBRvalue,-10.,10.) #
    Nsig_ggH = ROOT.RooFormulaVar("Nsig_ggH","@0*@1",ROOT.RooArgList(Ns_,B_R_))
    B_R_const = ROOT.RooRealVar("B_R_","branching_ratio",myBRvalue) #
    Nsig_ggH_const = ROOT.RooFormulaVar("Nsig_ggH","@0*@1",ROOT.RooArgList(Ns_,B_R_const))
    print('-------------------------------')
    print('INTEGRAL SIGNALG injected : ',Nsig_ggH.evaluate())
    print('INTEGRAL BKG injected : ',Nbkg)    
#    print('-------------------------------')
    argListN = ROOT.RooArgList(Nsig_ggH,Nbkg)
#    print('-------------------------------')
    print("arg list N size = ",argListN.getSize())
    print('-------------------------------')

    # set parameter used in the generation as constant
    cb_mu.setConstant(kTRUE)
    cb_sigma.setConstant(kTRUE)
    cb_alphaL.setConstant(kTRUE)
    cb_alphaR.setConstant(kTRUE)
    cb_nR.setConstant(kTRUE)
    cb_nL.setConstant(kTRUE)

    genPDF = ROOT.RooAddPdf("totPDF_chebychev",pdf1.GetName(),ROOT.RooArgList(pdf_crystalball, pdf1), ROOT.RooArgList(Nsig_ggH_const,Nbkg))
    fitPDF = ROOT.RooAddPdf("totPDF_bernstein",pdf2.GetName(),ROOT.RooArgList(pdf_crystalball, pdf2), ROOT.RooArgList(Nsig_ggH,Nbkg))
#    fitPDF = ROOT.RooAddPdf("totPDF_bernstein",pdf2.GetName(),ROOT.RooArgList(pdf_crystalball, pdf1), ROOT.RooArgList(Nsig_ggH,Nbkg))

    latex = TLatex()
    latex.SetTextColor(kGreen+1)
    latex.SetTextSize(0.04)
    latex.DrawLatex(110 ,0.10*dataBKG_full.GetMaximum(), pdf_bern3.GetName())
    latex.SetTextColor(kMagenta)
    latex.DrawLatex(110 ,0.20*dataBKG_full.GetMaximum(), pdf_chebychev3.GetName())

    latex.SetTextColor(kRed)
    latex.DrawLatex(110 ,0.30*dataBKG_full.GetMaximum(), genPDF.GetName())
    latex.SetTextColor(kBlue)
    latex.DrawLatex(110 ,0.40*dataBKG_full.GetMaximum(), fitPDF.GetName())

    multicanvas.Draw()
    multicanvas.SaveAs("~/public_html/Hrare/BIAS/toy1"+mesonCat+tag+".png")

    #
    # -------------------------------------------------------------------------------------
    #

    ntot = int(N+Nsig_ggH_const.evaluate())
    print('will generate ',ntot)

    mcstudy = RooMCStudy(
        genPDF,
        RooArgSet(x),
#        RooFit.Binned(1),
        RooFit.Silence(1),
        RooFit.FitModel(fitPDF),
        RooFit.Extended(1),
        RooFit.FitOptions(RooFit.Save(kTRUE),
                          RooFit.Minimizer("Minuit2"),RooFit.Strategy(2),
                          RooFit.Range("full"),
                          RooFit.PrintEvalErrors(0))) # -1 switches off printing.

    print('RooMCStudy set up')

    # Generate and fit 100 experiments of ntot events each
    mcstudy.generateAndFit(5000,ntot,kTRUE)
    print('gen and fit done')

    # -----------------------------------------------------------------------------
    # plot parameter of fit

    if False:
        multicanvas.cd(5)
        param1 = mcstudy.plotParam(cb_mu_2)
        param1.Draw()

        multicanvas.cd(6)
        frame2 = mcstudy.plotParam(cb_sigma_2)
        frame2.Draw()

        multicanvas.cd(7)
        frame3 = mcstudy.plotParam(cb_alphaL_2)
        frame3.Draw()

        multicanvas.cd(8)
        frame4 = mcstudy.plotParam(cb_alphaR_2)
        frame4.Draw()

        multicanvas.cd(9)
        frame5 = mcstudy.plotParam(cb_nL_2)
        frame5.Draw()

        multicanvas.cd(10)
        frame6 = mcstudy.plotParam(cb_nR_2)
        frame6.Draw()

    if tag=='_GFcat' and False:
        multicanvas.cd(13)
        param13 = mcstudy.plotParam(bern_c0)
        param13.Draw()

        multicanvas.cd(14)
        frame14 = mcstudy.plotParam(bern_c1)
        frame14.Draw()

        multicanvas.cd(15)
        frame15 = mcstudy.plotParam(bern_c2)
        frame15.Draw()

        multicanvas.cd(16)
        frame16 = mcstudy.plotParam(bern_c3)
        frame16.Draw()

    # -----------------------------------------------------------------------------
    # SIGNAL law

    multicanvas.cd(2)
    leg2 = ROOT.TLegend(0.12,0.75,0.44,0.97) #left positioning
    leg2.SetHeader(" ")
    leg2.SetBorderSize(0)
#    leg2.SetNColumns(1)
#    leg2.SetFillColorAlpha(0,0.)
#    leg2.SetLineColor(1)
#    leg2.SetLineStyle(1)
#    leg2.SetLineWidth(1)
#    leg2.SetFillStyle(1001)
    leg2.AddEntry(0,"gen pdf: "+genPDF.GetTitle(),"")
    leg2.AddEntry(0,"fit pdf: "+fitPDF.GetTitle(),"")

    ##
    
    genData0 = mcstudy.genData(0) 
    genData0.Print()
    multicanvas.cd(2)
    frame9=x.frame()
    genData0.plotOn(frame9)
    mcstudy.fitResult(0).Print("v")
    frame9.Draw()
    leg2.Draw()

    multicanvas.cd(3)
    frame7=x.frame()
    genData25 = mcstudy.genData(25) 
    genData25.Print()
    genData25.plotOn(frame7)
    mcstudy.fitResult(25).Print("v")
    frame7.Draw()

    multicanvas.cd(4)
    frame8=x.frame()
    genData47 = mcstudy.genData(47) 
    genData47.Print()
    genData47.plotOn(frame8)
    mcstudy.fitResult(47).Print("v")
    frame8.Draw()    

#    corrHist000 = mcstudy.fitResult(0).correlationHist("c027")
#    corrHist000.plotOn(frame9)

#    multicanvas.cd(10)
#    corrHist027 = mcstudy.fitResult(27).correlationHist("c027")
#    corrHist027.Draw("colz")

#    multicanvas.cd(11)
#    corrHist093 = mcstudy.fitResult(93).correlationHist("c093")
#    corrHist093.Draw("colz")

#    multicanvas.cd(12)
#    corrHist150 = mcstudy.fitResult(150).correlationHist("c150")
#    corrHist150.Draw("colz")

    multicanvas.Draw()
    multicanvas.SaveAs("~/public_html/Hrare/BIAS/toy1"+mesonCat+tag+".png")

    multicanvas.cd(12)
    BRpull_frame = mcstudy.plotPull(B_R_, ROOT.RooFit.Range(-5,5), ROOT.RooFit.Bins(100), ROOT.RooFit.FitGauss(1))
    BRpull_frame.SetTitle("")
    BRpull_frame.SetTitleOffset(1.5,"y")
    BRpull_frame.SetXTitle("Gen: "+genPDF.GetTitle()+" Vs Fit: "+fitPDF.GetTitle())
    BRpull_frame.SetMaximum(1.1*BRpull_frame.GetMaximum())
    BRpull_frame.Draw()

    multicanvas.cd(11)    
    BRpar_frame = mcstudy.plotParam(B_R_, ROOT.RooFit.Range(-5.,5.), ROOT.RooFit.Bins(200), ROOT.RooFit.FitGauss(1))
#    BRpar_frame = mcstudy.plotParam(B_R_, ROOT.RooFit.Range(0.9,1.1), ROOT.RooFit.Bins(20), ROOT.RooFit.FitGauss(1))
    BRpar_frame.SetTitle(mesonCat+tag)
    BRpar_frame.SetTitleOffset(1.5,"y")
    BRpar_frame.SetXTitle("Gen: "+genPDF.GetTitle()+" Vs Fit: "+fitPDF.GetTitle())
    BRpar_frame.SetMaximum(1.1*BRpar_frame.GetMaximum())
    BRpar_frame.Draw()

    multicanvas.Draw()

    multicanvas.SaveAs("~/public_html/Hrare/BIAS/toy1"+mesonCat+tag+".png")

    canvasPull.cd()
    BRpull_frame.Draw()
    latex2 = TLatex()
    latex2.SetTextSize(0.04)
    latex2.DrawLatex(-3 ,0.7*BRpull_frame.GetMaximum(), "BR="+str(myBRvalue))

    canvasPull.SaveAs("~/public_html/Hrare/BIAS/pull_toy"+mesonCat+tag+"BR"+str(myBRvalue)+".png")


if __name__ == "__main__":

    checkBias('_VBFcat','_PhiCat',2018, 0.001)
    checkBias('_VBFcatlow','_PhiCat',2018, 0.001)
    checkBias('_GFcat','_PhiCat',2018, 0.001)

    checkBias('_VBFcat','_RhoCat',2018, 0.001)
    checkBias('_VBFcatlow','_RhoCat',2018, 0.001)
    checkBias('_GFcat','_RhoCat',2018, 0.001)
