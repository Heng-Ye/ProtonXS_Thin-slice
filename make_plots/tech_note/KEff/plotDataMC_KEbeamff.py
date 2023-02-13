import ROOT as RT
import sys
import math

from argparse import ArgumentParser as ap
import numpy as np
from array import array
from math import sqrt
from math import exp


def fitg(x, p):
  m=p[0];
  s=p[1];
  n=p[2];
  g=n*math.exp(-(x[0]-m)*(x[0]-m)/(2*s*s));
  return g


def VNFit(h, pre_mean, n_sigma):
  pre_max=h.GetMaximum()
  pre_rms=h.GetRMS()
  print('pre_max:',pre_max)
  print('pre_rms:',pre_rms)
  print('pre_max:',pre_max)
  
  #pre-fit
  gg=RT.TF1("gg", fitg, pre_mean-n_sigma*pre_rms, pre_mean+n_sigma*pre_rms, 3)
  gg.SetParameter(0, pre_mean)
  gg.SetParameter(1, pre_rms)
  gg.SetParameter(2, pre_max)
  h.Fit("gg","remn")

  #pos-fit
  g=RT.TF1("g", fitg, gg.GetParameter(0)-n_sigma*gg.GetParameter(1), gg.GetParameter(0)+n_sigma*gg.GetParameter(1), 3)
  g.SetParameter(0, gg.GetParameter(0))
  g.SetParameter(1, gg.GetParameter(1))
  g.SetParameter(2, gg.GetParameter(2))

  g.SetLineStyle(2)
  g.SetLineWidth(2)

  h.Fit("g","remn");
  return g

#Run the code using the following commend -------------------------------------------------------------------------------------------------------------
#python plotDataMC_KEbeamff.py -d /dune/data2/users/hyliao/protonana/v09_39_01/KEFF_EDEPT/proton_beamxy_beammom_bkg_runAll.root -c ../mc_proton_beamxy_beammom_nobmrw_edepteloss.root -o plots_beamxy_nobmrw_EdeptEloss

#python plotDataMC_KEbeamff.py -d /dune/data2/users/hyliao/protonana/v09_39_01/KEFF_EDEPT/proton_beamxy_beammom_bkg_runAll.root -c ../mc_proton_beamxy_beammom_nobmrw_edepteloss.root -crw ../mc_proton_beamxy_beammom_bmrw_edepteloss.root -o plots_beamxy_nobmrw_EdeptEloss

#--------------------------------------------------------------------------------------
parser = ap()

RT.gROOT.LoadMacro("~/protoDUNEStyle.C")
RT.gROOT.SetBatch(); #When running in batch mode, PyROOT does NOT display any graphics
RT.gStyle.SetOptStat(00000)
#RT.gStyle.SetErrorX(1.e-4)
RT.gStyle.SetTitleAlign(23)
RT.gStyle.SetTitleX(.5)
#RT.gStyle.SetLineWidth(1)
tt = RT.TLatex();
tt.SetNDC();

#Read files -----------------------------------------------------------------------------------------
parser.add_argument("-d", type=str, help='Data File', default = "")
parser.add_argument("-c", type=str, help='MC File', default = "")
parser.add_argument("-crw", type=str, help='MCrw File', default = "")
parser.add_argument("-o", type=str, help='Output folder', default = "")

args=parser.parse_args()
if not (args.d and args.o and args.c):
  print("--> Please provide all files (data, MC, MR_rw, output_file_folder). Thank you!")
  exit()

if (args.d): print('Read data: '+args.d)
if (args.c): print('MC: '+args.c)
if (args.crw): print('MCrw: '+args.c)
if (args.o): print('Output folder: '+args.o)

#Read data histograms ----------------------------------
f_data=RT.TFile(args.d, "OPEN")
#f_data.ls()

#beamff
keffbeam_stop_data=f_data.Get("h1d_keffbeam_stop")
n_ffbeam_stop_data=keffbeam_stop_data.Integral()
print('n_ffbeam_stop_data=',n_ffbeam_stop_data)

#Read MC histograms ----------------------------------
f_mc=RT.TFile(args.c, "OPEN")

#ff-truth
keff_stop_mc=f_mc.Get("h1d_keff_stop");
keffbeam_stop_mc=f_mc.Get("h1d_keffbeam_stop")
n_ff_stop_mc=keff_stop_mc.Integral()
n_ffbeam_stop_mc=keffbeam_stop_mc.Integral()
print('n_ff_stop_mc=',n_ff_stop_mc)
print('n_ffbeam_stop_mc=',n_ffbeam_stop_mc)

#Read MC(rw) histograms ----------------------------------
f_mc_rw=RT.TFile(args.crw, "OPEN")

#ff-truth
keff_stop_mc_rw=f_mc_rw.Get("h1d_keff_stop");
keffbeam_stop_mc_rw=f_mc_rw.Get("h1d_keffbeam_stop")
keff_stop_mc_rw.SetName("keff_stop_mc_rw")
keffbeam_stop_mc_rw.SetName("keffbeam_stop_mc_rw")

n_ff_stop_mc_rw=keff_stop_mc_rw.Integral()
n_ffbeam_stop_mc_rw=keffbeam_stop_mc_rw.Integral()
print('n_ff_stop_mc_rw=',n_ff_stop_mc_rw)
print('n_ffbeam_stop_mc_rw=',n_ffbeam_stop_mc_rw)

#normalization [mc] ------------------------------------------------
keffbeam_stop_mc.Scale(n_ffbeam_stop_data/n_ffbeam_stop_mc)
keff_stop_mc.Scale(n_ffbeam_stop_data/n_ff_stop_mc)

#normalization [mcrw] ------------------------------------------------
keffbeam_stop_mc_rw.Scale(n_ffbeam_stop_data/n_ffbeam_stop_mc_rw)
keff_stop_mc_rw.Scale(n_ffbeam_stop_data/n_ff_stop_mc_rw)

#fitting to extract mu and sigma ----------------------------
#data
fit_keffbeam_stop_data=VNFit(keffbeam_stop_data, 410, 3)

#mc
fit_keff_stop_mc=VNFit(keff_stop_mc, 400, 3)
fit_keffbeam_stop_mc=VNFit(keffbeam_stop_mc, 400, 3)

#mc(rw)
fit_keff_stop_mc_rw=VNFit(keff_stop_mc_rw, 400, 3)
fit_keffbeam_stop_mc_rw=VNFit(keffbeam_stop_mc_rw, 400, 3)

#get fitted mu and sigma [data] -----------------
#data
n_data_stop=1
mu_stop_data=[]
er_mu_stop_data=[]
sigma_stop_data=[]
er_sigma_stop_data=[]

for i in range(n_data_stop):
  fit=fit_keffbeam_stop_data
  if i==0: 
    fit=fit_keffbeam_stop_data
  m=fit.GetParameter(0)
  er_m=fit.GetParError(0)
  s=fit.GetParameter(1)
  er_s=fit.GetParError(1)
  mu_stop_data.append(m)
  er_mu_stop_data.append(er_m)
  sigma_stop_data.append(s)
  er_sigma_stop_data.append(er_s)
  print("Data i=",i," m=",m, "s=",s,"\n")

#mc
n_mc_stop=4
mu_stop_mc=[]
er_mu_stop_mc=[]
sigma_stop_mc=[]
er_sigma_stop_mc=[]

for i in range(n_mc_stop):
  fit=fit_keff_stop_mc
  if i==0: 
    fit=fit_keff_stop_mc
  if i==1: 
    fit=fit_keffbeam_stop_mc
  if i==2: 
    fit=fit_keff_stop_mc_rw
  if i==3: 
    fit=fit_keffbeam_stop_mc_rw

  m=fit.GetParameter(0)
  er_m=fit.GetParError(0)
  s=fit.GetParameter(1)
  er_s=fit.GetParError(1)
  mu_stop_mc.append(m)
  er_mu_stop_mc.append(er_m)
  sigma_stop_mc.append(s)
  er_sigma_stop_mc.append(er_s)
  print("MC i=",i," m=",m, "s=",s,"\n")

#Plots #########################################################################################################################################################
#[1]KEBeamff
c0=RT.TCanvas("c0","",1200,900)
c0.Divide(1,1)
c0.cd(1)

f2d_keff=RT.TH2D("f2d_keff","", 370, 220, 600, 600, 0, keffbeam_stop_data.GetBinContent(keffbeam_stop_data.GetMaximumBin())+400)
f2d_keff.SetTitle("Stopping Protons;Proton Kinetic Energy at TPC FF [MeV];")
f2d_keff.GetXaxis().CenterTitle()
f2d_keff.Draw()

keffbeam_stop_data.SetLineColor(1)
keffbeam_stop_data.SetMarkerColor(1)
fit_keffbeam_stop_data.SetLineColor(1)
keffbeam_stop_data.SetLineWidth(1)

keff_stop_mc.SetLineColor(2)
keff_stop_mc.SetMarkerColor(2)
fit_keff_stop_mc.SetLineColor(2)
keff_stop_mc.SetLineWidth(1)

keffbeam_stop_mc.SetLineColor(4)
keffbeam_stop_mc.SetMarkerColor(4)
fit_keffbeam_stop_mc.SetLineColor(4)
keffbeam_stop_mc.SetLineWidth(1)


fit_keffbeam_stop_data.Draw("same")
fit_keff_stop_mc.Draw("same")
fit_keffbeam_stop_mc.Draw("same")

keffbeam_stop_mc.Draw("ep same")
keff_stop_mc.Draw(" same")
keffbeam_stop_data.Draw(" same")

leg0=RT.TLegend(0.14,0.65,.9,0.9)
leg0.SetFillStyle(0)
txt0=[]
txt0.append("Data: #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_data[0],er_mu_stop_data[0],sigma_stop_data[0],er_sigma_stop_data[0]))
txt0.append("MC(truth): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[0],er_mu_stop_mc[0],sigma_stop_mc[0],er_sigma_stop_mc[0]))
txt0.append("MC(spec): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[1],er_mu_stop_mc[1],sigma_stop_mc[1],er_sigma_stop_mc[1]))

print(txt0[0])
leg0.AddEntry(keffbeam_stop_data, txt0[0], "ep")
leg0.AddEntry(keff_stop_mc, txt0[1], "l")
leg0.AddEntry(keffbeam_stop_mc, txt0[2], "l")
#leg0.SetNColumns(2);

leg0.Draw()

c0.Print(args.o+'/kebeamff_data_mc.eps')

"""
#Fitted result (nobmrw) 
#keffbeam_stop_data
   1  p0           4.09837e+02   2.93549e-01  -6.91597e-05   2.09147e-05
   2  p1           4.28448e+01   2.02807e-01   1.95520e-04  -1.15037e-04
   3  p2           4.11054e+02   3.35511e+00   3.35511e+00   4.19084e-05
#keff_stop_mc
   1  p0           4.15916e+02   2.41791e-01  -2.65188e-06   1.04777e-06
   2  p1           3.99959e+01   1.55280e-01   9.25748e-05  -5.88385e-05
   3  p2           4.38820e+02   3.09842e+00   3.09842e+00   3.72791e-05
#fit_keffbeam_stop_mc
   1  p0           4.15248e+02   2.24060e-01   1.38021e-05  -4.12830e-06
   2  p1           3.72018e+01   1.58229e-01   1.60618e-04  -1.49393e-04
   3  p2           4.75150e+02   3.44224e+00   3.44224e+00   3.04656e-05
"""

#[2]KEBeamff (MCRW using truth)
c1=RT.TCanvas("c0","",1200,900)
c1.Divide(1,1)
c1.cd(1)

f2d_keff1=RT.TH2D("f2d_keff1","", 370, 220, 600, 600, 0, keffbeam_stop_data.GetBinContent(keffbeam_stop_data.GetMaximumBin())+400)
f2d_keff1.SetTitle("Stopping Protons;Proton Kinetic Energy at TPC FF [MeV];")
f2d_keff1.GetXaxis().CenterTitle()
f2d_keff1.Draw()

keff_stop_mc_rw.SetLineColor(2)
keff_stop_mc_rw.SetMarkerColor(2)
fit_keff_stop_mc_rw.SetLineColor(2)
keff_stop_mc_rw.SetLineWidth(1)

keffbeam_stop_mc_rw.SetLineColor(4)
keffbeam_stop_mc_rw.SetMarkerColor(4)
fit_keffbeam_stop_mc_rw.SetLineColor(4)
keffbeam_stop_mc_rw.SetLineWidth(1)


fit_keffbeam_stop_data.Draw("same")
fit_keff_stop_mc_rw.Draw("same")
fit_keffbeam_stop_mc_rw.Draw("same")

keffbeam_stop_mc_rw.Draw("ep same")
keff_stop_mc_rw.Draw(" same")
keffbeam_stop_data.Draw(" same")

leg1=RT.TLegend(0.14,0.65,.9,0.9)
leg1.SetFillStyle(0)
txt1=[]
txt1.append("Data: #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_data[0],er_mu_stop_data[0],sigma_stop_data[0],er_sigma_stop_data[0]))
txt1.append("MC(truth): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[2],er_mu_stop_mc[2],sigma_stop_mc[2],er_sigma_stop_mc[2]))
txt1.append("MC(spec): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[3],er_mu_stop_mc[3],sigma_stop_mc[3],er_sigma_stop_mc[3]))

print(txt1[0])
leg1.AddEntry(keffbeam_stop_data, txt1[0], "ep")
leg1.AddEntry(keff_stop_mc_rw, txt1[1], "l")
leg1.AddEntry(keffbeam_stop_mc_rw, txt1[2], "l")
leg1.Draw()

c1.Print(args.o+'/kebeamff_data_mc_bmrw1_ff.eps')


