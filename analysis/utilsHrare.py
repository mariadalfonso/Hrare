import ROOT
import os
from subprocess import call,check_output

def findDataset(name):

    DASclient = "dasgoclient -query '%(query)s'"
    cmd= DASclient%{'query':'file dataset=%s'%name}
    print(cmd)
    check_output(cmd,shell=True)
    fileList=[ 'root://xrootd-cms.infn.it//'+x for x in check_output(cmd,shell=True).split() ]

    files_ROOT = ROOT.vector('string')()
    for f in fileList: files_ROOT.push_back(f)

    return files_ROOT

def findDIR(directory):

    rootFiles = ROOT.vector('string')()
    for root, directories, filenames in os.walk(directory):
        for f in filenames:

            filePath = os.path.join(os.path.abspath(root), f)
            if "failed/" in filePath: continue
            if "log/" in filePath: continue
            rootFiles.push_back(filePath)

    return rootFiles




def plot(h,filename,doLogX,color):

   ROOT.gStyle.SetOptStat(1);
   ROOT.gStyle.SetTextFont(42)
   c = ROOT.TCanvas("c", "", 800, 700)
   if doLogX: c.SetLogx();
#  c.SetLogy()

   h.SetTitle("")
   h.GetXaxis().SetTitleSize(0.04)
   h.GetYaxis().SetTitleSize(0.04)
   h.SetLineColor(color)

   h.Draw()

   label = ROOT.TLatex(); label.SetNDC(True)
   label.DrawLatex(0.175, 0.740, "#eta")
   label.DrawLatex(0.205, 0.775, "#rho,#omega")
   label.DrawLatex(0.270, 0.740, "#phi")
   label.SetTextSize(0.040); label.DrawLatex(0.100, 0.920, "#bf{CMS Simulation}")
   label.SetTextSize(0.030); label.DrawLatex(0.630, 0.920, "#sqrt{s} = 13 TeV, L_{int} = X fb^{-1}")

   print("saving file: {} ".format(filename))

   c.SaveAs(filename)
