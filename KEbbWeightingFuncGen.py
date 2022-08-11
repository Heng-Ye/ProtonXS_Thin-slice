import ROOT as RT
import sys
import math

from argparse import ArgumentParser as ap
import numpy as np
from array import array
from math import sqrt
from math import exp

#Run as:python KEbbWeightingFuncGen.py -d ./rw/kebb_reweight_data.root -c ./rw/kebb_reweight_mc.root -o ./make_plots/plots_beamxy_bmrw_calo_hyper

#--------------------------------------------------------------------------------------
parser = ap()

RT.gROOT.LoadMacro("~/protoDUNEStyle.C")
RT.gROOT.SetBatch(); #When running in batch mode, PyROOT does NOT display any graphics
RT.gStyle.SetOptStat(00000)
#RT.gStyle.SetErrorX(1.e-4)
RT.gStyle.SetTitleAlign(23)
RT.gStyle.SetTitleX(.5)
RT.gStyle.SetLineWidth(1)
tt = RT.TLatex();
tt.SetNDC();

#Read files -----------------------------------------------------------------------------------------
parser.add_argument("-d", type=str, help='Data File', default = "")
parser.add_argument("-c", type=str, help='MC File', default = "")
parser.add_argument("-o", type=str, help='Output folder', default = "")

args=parser.parse_args()

if not (args.d and args.o and args.c):
  print("--> Please provide all files (data, MC, output_file_folder)...")
  exit()
if (args.d): print('Read data: '+args.d)
if (args.c): print('MC: '+args.c)
if (args.o): print('Output folder: '+args.o)

#Read MC histograms -------------------
f_mc=RT.TFile(args.c, "OPEN")
mu_len_ref_mc=f_mc.Get("mu_len_ref")
sg_len_ref_mc=f_mc.Get("sg_len_ref")
mu_len_exp_mc=f_mc.Get("mu_len_exp")
sg_len_exp_mc=f_mc.Get("sg_len_exp")

#Read data histograms -----------------------
f_data=RT.TFile(args.d, "OPEN")
mu_len_ref_data=f_data.Get("mu_len_ref")
sg_len_ref_data=f_data.Get("sg_len_ref")
mu_len_exp_data=f_data.Get("mu_len_exp")
sg_len_exp_data=f_data.Get("sg_len_exp")

#f_data=RT.TFile(args.d, "OPEN")
#kebeam_data=f_data.Get("h1d_kebeam")

c0=RT.TCanvas("c0","",1200,1600)
c0.Divide(1,2)
c0.cd(1)

f2d_mu=RT.TH2D("f2d_mu","", 130, 0, 130, 500, 0, 500)
f2d_mu.SetTitle(";Length [cm]; #mu of KE_{bb}[MeV];")
f2d_mu.GetXaxis().CenterTitle()
f2d_mu.Draw()
mu_len_exp_mc.SetLineColor(2)
mu_len_exp_mc.SetMarkerColor(2)

mu_len_ref_data.SetLineColor(1)
mu_len_ref_data.SetMarkerColor(1)

mu_len_exp_data.SetLineColor(4)
mu_len_exp_data.SetMarkerColor(4)

mu_len_exp_mc.Draw("p same")
mu_len_ref_data.Draw("p same")
mu_len_exp_data.Draw("p same")

leg=RT.TLegend(0.44,0.65,.9,0.87)
leg.SetFillStyle(0)
leg.AddEntry(mu_len_ref_data, "Data(fit)", "p")
leg.AddEntry(mu_len_exp_data, "Data(measurement)", "p")
leg.AddEntry(mu_len_exp_mc, "MC(measurement)", "p")
leg.Draw()



c0.cd(2)
f2d_sg=RT.TH2D("f2d_sg","", 130, 0, 130, 50, 30, 80)
f2d_sg.SetTitle(";Length [cm]; #sigma of KE_{bb}[MeV];")
f2d_sg.GetXaxis().CenterTitle()
f2d_sg.Draw()

sg_len_exp_mc.SetLineColor(2)
sg_len_exp_mc.SetMarkerColor(2)

sg_len_ref_data.SetLineColor(1)
sg_len_ref_data.SetMarkerColor(1)

sg_len_exp_data.SetLineColor(4)
sg_len_exp_data.SetMarkerColor(4)

sg_len_exp_mc.Draw("p same")
sg_len_ref_data.Draw("p same")
sg_len_exp_data.Draw("p same")

leg2=RT.TLegend(0.1,0.65,.4,0.87)
leg2.SetFillStyle(0)
leg2.AddEntry(mu_len_ref_data, "Data(fit)", "p")
leg2.AddEntry(mu_len_exp_data, "Data(measurement)", "p")
leg2.AddEntry(mu_len_exp_mc, "MC(measurement)", "p")
leg2.Draw()


c0.Print(args.o+'/kebb_mu_sg.eps')

