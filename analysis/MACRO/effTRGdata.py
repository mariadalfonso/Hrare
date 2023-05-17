import ROOT
import os
from array import array
import math

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gROOT.SetBatch()

ROOT.ROOT.EnableImplicitMT()
RDataFrame = ROOT.RDataFrame

#mesonCat = '_PhiCat'
mesonCat = '_RhoCat'
#anaCat = '_GFcat'
#anaCat = '_VBFcat'
#anaCat = '_VBFcatlow'
anaCat = '_Wcat'
#anaCat = '_Zcat'

year = '12016'

mytree = ROOT.TChain('events')
if year == '2018':
   mytree.Add('../TRGstudy/2018/outname_mc-1_Wcat'+mesonCat+'_2018_triggerData.root')
   mytree.Add('../TRGstudy/2018/outname_mc-2_Wcat'+mesonCat+'_2018_triggerData.root')
   mytree.Add('../TRGstudy/2018/outname_mc-3_Wcat'+mesonCat+'_2018_triggerData.root')
   mytree.Add('../TRGstudy/2018/outname_mc-4_Wcat'+mesonCat+'_2018_triggerData.root')
   mytree.Add('../TRGstudy/2018/outname_mc-11_Wcat'+mesonCat+'_2018_triggerData.root')
   mytree.Add('../TRGstudy/2018/outname_mc-12_Wcat'+mesonCat+'_2018_triggerData.root')
   mytree.Add('../TRGstudy/2018/outname_mc-13_Wcat'+mesonCat+'_2018_triggerData.root')
   mytree.Add('../TRGstudy/2018/outname_mc-14_Wcat'+mesonCat+'_2018_triggerData.root')   
if year == '2017':
   mytree.Add('../TRGstudy/2017/outname_mc-6_Wcat'+mesonCat+'_2017_triggerData.root')
   mytree.Add('../TRGstudy/2017/outname_mc-16_Wcat'+mesonCat+'_2017_triggerData.root')   
if year == '12016':
#   mytree.Add('../TRGstudy/12016/outname_mc-1_Wcat'+mesonCat+'_12016_triggerData.root')
   mytree.Add('../TRGstudy/12016/outname_mc-2_Wcat'+mesonCat+'_12016_triggerData.root')
   mytree.Add('../TRGstudy/12016/outname_mc-3_Wcat'+mesonCat+'_12016_triggerData.root')
   mytree.Add('../TRGstudy/12016/outname_mc-4_Wcat'+mesonCat+'_12016_triggerData.root')
   mytree.Add('../TRGstudy/12016/outname_mc-5_Wcat'+mesonCat+'_12016_triggerData.root')
   mytree.Add('../TRGstudy/12016/outname_mc-6_Wcat'+mesonCat+'_12016_triggerData.root')
#   mytree.Add('../TRGstudy/12016/outname_mc-11_Wcat'+mesonCat+'_12016_triggerData.root')
   mytree.Add('../TRGstudy/12016/outname_mc-12_Wcat'+mesonCat+'_12016_triggerData.root')
   mytree.Add('../TRGstudy/12016/outname_mc-13_Wcat'+mesonCat+'_12016_triggerData.root')
   mytree.Add('../TRGstudy/12016/outname_mc-14_Wcat'+mesonCat+'_12016_triggerData.root')
   mytree.Add('../TRGstudy/12016/outname_mc-15_Wcat'+mesonCat+'_12016_triggerData.root')
   mytree.Add('../TRGstudy/12016/outname_mc-16_Wcat'+mesonCat+'_12016_triggerData.root')         
   
mytreeMC = ROOT.TChain('events')
if year == '2018':
   mytreeMC.Add('../TRGstudy/2018/outname_mc2_Wcat'+mesonCat+'_2018_triggerData.root')
   mytreeMC.Add('../TRGstudy/2018/outname_mc1_Wcat'+mesonCat+'_2018_triggerData.root')
   mytreeMC.Add('../TRGstudy/2018/outname_mc31_Wcat'+mesonCat+'_2018_triggerData.root')
   mytreeMC.Add('../TRGstudy/2018/outname_mc32_Wcat'+mesonCat+'_2018_triggerData.root')
   mytreeMC.Add('../TRGstudy/2018/outname_mc33_Wcat'+mesonCat+'_2018_triggerData.root')
   mytreeMC.Add('../TRGstudy/2018/outname_mc34_Wcat'+mesonCat+'_2018_triggerData.root')
   mytreeMC.Add('../TRGstudy/2018/outname_mc35_Wcat'+mesonCat+'_2018_triggerData.root')
   mytreeMC.Add('../TRGstudy/2018/outname_mc36_Wcat'+mesonCat+'_2018_triggerData.root')
if year == '2017':
   mytreeMC.Add('../TRGstudy/2017/outname_mc2_Wcat'+mesonCat+'_2017_triggerData.root')
   mytreeMC.Add('../TRGstudy/2017/outname_mc1_Wcat'+mesonCat+'_2017_triggerData.root')
   mytreeMC.Add('../TRGstudy/2017/outname_mc31_Wcat'+mesonCat+'_2017_triggerData.root')
   mytreeMC.Add('../TRGstudy/2017/outname_mc32_Wcat'+mesonCat+'_2017_triggerData.root')
   mytreeMC.Add('../TRGstudy/2017/outname_mc33_Wcat'+mesonCat+'_2017_triggerData.root')
   mytreeMC.Add('../TRGstudy/2017/outname_mc34_Wcat'+mesonCat+'_2017_triggerData.root')
   mytreeMC.Add('../TRGstudy/2017/outname_mc35_Wcat'+mesonCat+'_2017_triggerData.root')
   mytreeMC.Add('../TRGstudy/2017/outname_mc36_Wcat'+mesonCat+'_2017_triggerData.root')
if year == '12016':
   mytreeMC.Add('../TRGstudy/12016/outname_mc2_Wcat'+mesonCat+'_12016_triggerData.root')
   mytreeMC.Add('../TRGstudy/12016/outname_mc1_Wcat'+mesonCat+'_12016_triggerData.root')
##mytreeMC.Add('../TRGstudy/12016/outname_mc31_Wcat'+mesonCat+'_12016_triggerData.root')
   mytreeMC.Add('../TRGstudy/12016/outname_mc32_Wcat'+mesonCat+'_12016_triggerData.root')
   mytreeMC.Add('../TRGstudy/12016/outname_mc33_Wcat'+mesonCat+'_12016_triggerData.root')
   mytreeMC.Add('../TRGstudy/12016/outname_mc34_Wcat'+mesonCat+'_12016_triggerData.root')
   mytreeMC.Add('../TRGstudy/12016/outname_mc35_Wcat'+mesonCat+'_12016_triggerData.root')
   mytreeMC.Add('../TRGstudy/12016/outname_mc36_Wcat'+mesonCat+'_12016_triggerData.root')
   #mytreeMC.Add('../TRGstudy/22016/outname_mc2_Wcat'+mesonCat+'_22016_triggerData.root')
   #mytreeMC.Add('../TRGstudy/22016/outname_mc1_Wcat'+mesonCat+'_22016_triggerData.root')
   #mytreeMC.Add('../TRGstudy/22016/outname_mc31_Wcat'+mesonCat+'_22016_triggerData.root')
   #mytreeMC.Add('../TRGstudy/22016/outname_mc32_Wcat'+mesonCat+'_22016_triggerData.root')
   #mytreeMC.Add('../TRGstudy/22016/outname_mc33_Wcat'+mesonCat+'_22016_triggerData.root')
   #mytreeMC.Add('../TRGstudy/22016/outname_mc34_Wcat'+mesonCat+'_22016_triggerData.root')
   #mytreeMC.Add('../TRGstudy/22016/outname_mc35_Wcat'+mesonCat+'_22016_triggerData.root')
   #mytreeMC.Add('../TRGstudy/22016/outname_mc36_Wcat'+mesonCat+'_22016_triggerData.root')

# Create the plot                                                                                                                                                                                           
def createCanvasPads(doLog):

   # Create canvas with pad                                                                                                                                                                                 
   c = ROOT.TCanvas("c", "", 600, 600)
   pad = ROOT.TPad("upper_pad", "", 0, 0, 1, 1)

   if doLog : pad.SetLogy(1)
   pad.SetTickx(False)
   pad.SetTicky(False)
   pad.Draw()
   pad.cd()

   ydiv = 0.2
   pad1 = ROOT.TPad("upper_pad", "", 0.,ydiv,1.,1.)
   if doLog : pad1.SetLogy(1)
   pad2 = ROOT.TPad("lower_pad", "", 0.,0.,1.,ydiv)

   pad1.Draw()
   pad2.Draw()

   return c, pad1, pad2


def createGraph(itemTRG, item, tree):

   nEntries = tree.GetEntries()
   print('nEntries=',nEntries)

   xmin = 20.
   xmax = 200.
   nBins = int((xmax-xmin)/5)

   h_num = ROOT.TH1F('h_num', ' ', nBins, xmin, xmax)
   h_den = ROOT.TH1F('h_den', ' ', nBins, xmin, xmax)   

   df = RDataFrame(tree)

   offline=False
   if item == 1:
      var = "goodPhotons_pt[0]"
      if year == '2018': offline = "true"
      if year == '2017' or year == '12016': offline = "abs(goodPhotons_eta[0])<1.4"

   h_num = df.Define("var","{}".format(var)).Define("offlineSel","{}".format(offline)).Define("weight","w*lumiIntegrated").Filter("triggerMu17_Photon30>0 and triggerMu>0 and offlineSel").Histo1D(("h_num","h_num",nBins,xmin,xmax),"var","weight")
   h_den = df.Define("var","{}".format(var)).Define("offlineSel","{}".format(offline)).Define("weight","w*lumiIntegrated").Filter("triggerMu>0 and offlineSel").Histo1D(("h_num","h_num",nBins,xmin,xmax),"var","weight")   

   print('h (passed) = ', h_num.Integral() )
   print('h (total) = ', h_den.Integral() )   
   print('nBins=',nBins, ' nnum=' ,h_num.GetNbinsX(), ' nden=' ,h_den.GetNbinsX(), )

   n = h_num.GetNbinsX()
   x, ex, eff, eff_error_low, eff_error_high = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )

   for myBin in range( n ):

      x.append(h_num.GetXaxis().GetBinLowEdge(myBin))
      #ex.append(h_num.GetBinWidth())
      ex.append(1)
      passed = h_num.Integral(h_num.FindBin(h_num.GetXaxis().GetBinLowEdge(myBin)),h_num.GetNbinsX()+1)
      total = h_den.Integral(h_den.FindBin(h_den.GetXaxis().GetBinLowEdge(myBin)),h_den.GetNbinsX()+1)
      #print('pt=',myBin,'  passed=',passed,'  total=',total,'  eff=',passed/total)
      eff.append(passed/total)
      err_low_ = (passed/total) - ROOT.TEfficiency.ClopperPearson(total , passed,  0.683,  False)
      err_high_ = ROOT.TEfficiency.ClopperPearson(total, passed,  0.683,  True) - (passed/total)
      eff_error_low.append(err_low_)
      eff_error_high.append(err_high_)

   print(len(x))
#   graph = ROOT.TGraphAsymmErrors(len(x),x,eff,ex,ex,eff_error_low,eff_error_high)
#   graph = ROOT.TGraph(len(x),x,eff)
   graph = ROOT.TGraphErrors(n,x,eff,ex,eff_error_low)

   graph.SetMarkerStyle(20)
   graph.GetXaxis().SetTitle("p_{T}^{#gamma}")
   graph.GetHistogram().GetYaxis().SetRangeUser(0.5,1.5)   

   return graph


def plot():

   mg = ROOT.TMultiGraph();
   yearString = year
   if year == '2018': yearString = '2018 B-D'
   if year == '2017': yearString = '2017 F'
   if year == '12016': yearString = '2016 B-F'
   mg.SetTitle(yearString)

   #sigmoid = ROOT.TF1("sigmoid","[0]/(1+exp([1]*([2]-x)))",0.5,1.);   
   
   #triggerAna data
   print("-- DATA")
   gr_num1 = createGraph(0, 1, mytree) # Hmass
   gr_num1.SetMarkerColor(ROOT.kBlue)
   gr_num1.SetLineColor(ROOT.kBlue)
   mg.Add(gr_num1);
   #sigmoid.SetParameter(0, 0.98); 
   #sigmoid.SetParameter(1, 0.1); 
   #sigmoid.SetParameter(2, 20); 
   #gr_num1.Fit("sigmoid")
   #fit_data_l=gr_num1.GetFunction("sigmoid")
   #fit_data_l.SetLineColor(ROOT.kBlue)
   #fit_data_l.SetLineStyle(ROOT.kDashed)
   
   #triggerAna MC
   print("-- MC")   
   gr_num1_mc = createGraph(0, 1, mytreeMC) # Hmass
   gr_num1_mc.SetMarkerColor(ROOT.kRed)   
   gr_num1_mc.SetLineColor(ROOT.kRed)
   mg.Add(gr_num1_mc);
   #sigmoid.SetParameter(0, 0.98); 
   #sigmoid.SetParameter(1, 0.1); 
   #sigmoid.SetParameter(2, 20); 
   #gr_num1_mc.Fit("sigmoid")
   #fit_mc_l = gr_num1_mc.GetFunction("sigmoid")
   #fit_mc_l.SetLineColor(ROOT.kRed)
   #fit_mc_l.SetLineStyle(ROOT.kDashed)
   
   doLog=False
   c, pad1, pad2 = createCanvasPads(doLog)

   pad1.cd() 

   mg.Draw("A P E");
   mg.GetHistogram().GetYaxis().SetRangeUser(0.5,1.5)
   
   line = ROOT.TLine(75, 0.5, 75, 1.)
   line.SetLineColor(ROOT.kGray)
   line.SetLineStyle(3)   
   line.Draw()

   line2 = ROOT.TLine(35, 0.5, 35, 1.)
   line2.SetLineColor(ROOT.kGray)
   line2.SetLineStyle(3)   
   line2.Draw()

   latex = ROOT.TLatex()
   latex.SetTextSize(0.04)
   latex.SetTextColor(ROOT.kBlue)
   latex.DrawLatex(130, 1.3,  "Single Muon")
   latex.SetTextColor(ROOT.kRed)
   latex.DrawLatex(130, 1.2,  "MC V+jets")
   latex.SetTextColor(ROOT.kBlack)
   if year == '2018': latex.DrawLatex(50, 1.25,  "WP=90")
   else: latex.DrawLatex(50, 1.25,  "WP=90 barrel")
   
   if year == '2018': c.SaveAs('~/public_html/Hrare/TRG/PhotonFromData'+mesonCat+year+'.png')
   else: c.SaveAs('~/public_html/Hrare/TRG/PhotonFromData'+mesonCat+year+'_WP90_barrel.png')
   
   pad2.cd() 
   
   ratio=ROOT.TF1("ratio","fit_mc_l/fit_data_l",20.,180)

   ratio.GetYaxis().SetRangeUser(0.8,1.2)
   ratio.GetXaxis().SetTitleOffset(4.)
   ratio.GetXaxis().SetTitleSize(0.15)
   ratio.GetXaxis().SetLabelSize(0.12)
   
   ratio.GetYaxis().SetTitleOffset(0.3)
   ratio.GetYaxis().SetTitleSize(0.15)
   ratio.GetYaxis().SetLabelSize(0.12)
   ratio.GetYaxis().SetTitle("SF(#epsilon^{data}/#epsilon^{MC})");
   
   ratio.Draw("hist")

   line3 = ROOT.TLine(20, 1., 200, 1.)
   line3.SetLineColor(ROOT.kGray)
   line3.SetLineStyle(1)
   line3.Draw()
   
if __name__ == "__main__":

   plot()
   
   exit()
