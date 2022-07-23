import ROOT as RT
import sys
import math

from argparse import ArgumentParser as ap
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

#--------------------------------------------------------------------------------------
parser = ap()

RT.gROOT.LoadMacro("~/protoDUNEStyle.C")
RT.gROOT.SetBatch(); #When running in batch mode, PyROOT does NOT display any graphics
RT.gStyle.SetOptStat(00000)
#RT.gStyle.SetErrorX(1.e-4)
RT.gStyle.SetTitleAlign(23)
RT.gStyle.SetTitleX(.5)
tt = RT.TLatex();
tt.SetNDC();

#Read files -----------------------------------------------------------------------------------------
parser.add_argument("-d", type=str, help='Data File', default = "")
parser.add_argument("-c", type=str, help='MC File', default = "")
parser.add_argument("-o", type=str, help='Output folder', default = "")

#fdata='/dune/data2/users/hyliao/protonana/v09_39_01/KEHY_old/proton_beamxy_beammom_runTMP_new.root'
#fmc=
#fmc_bmrw=
#outpath=
args=parser.parse_args()

if not (args.d and args.o and args.c):
  print("--> Please provide all files (data, MC, MR_rw, output_file_folder). Thank you!")
  exit()

if (args.d): print('Read data: '+args.d)
if (args.c): print('MC data: '+args.c)
if (args.o): print('Output folder: '+args.o)


#Read data histograms ----------------------------------
f_data=RT.TFile(args.d, "OPEN")
#f_data.ls()

#beam
kebeam_data=f_data.Get("h1d_kebeam")
kebeam_stop_data=f_data.Get("h1d_kebeam_stop")

#ff
kehy_data=f_data.Get("h1d_kehy")
kehy_stop_data=f_data.Get("h1d_kehy_stop")
kehy_inel_data=f_data.Get("h1d_kehy_inel")

#E-dept
kerange_stop_data=f_data.Get("h1d_kerange_stop")
kecalo_stop_data=f_data.Get("h1d_kecalo_stop")

#ke-end[calo]
kend_calo_stop_data=f_data.Get("h1d_kend_calo_stop")
kend_calo_el_data=f_data.Get("h1d_kend_calo_el")
kend_calo_inel_data=f_data.Get("h1d_kend_calo_inel")

#ke-end[bb]
kend_bb_stop_data=f_data.Get("h1d_kend_bb_stop")
kend_bb_el_data=f_data.Get("h1d_kend_bb_el")
kend_bb_inel_data=f_data.Get("h1d_kend_bb_inel")



#get entries 
n_beam_stop_data=kebeam_stop_data.Integral()

n_hy_stop_data=kehy_stop_data.Integral()
n_hy_inel_data=kehy_inel_data.Integral()

n_range_stop_data=kerange_stop_data.Integral()
n_calo_stop_data=kecalo_stop_data.Integral()

n_end_bb_stop_data=kend_bb_stop_data.Integral()
n_end_bb_el_data=kend_bb_el_data.Integral()
n_end_bb_inel_data=kend_bb_inel_data.Integral()

n_end_calo_stop_data=kend_calo_stop_data.Integral()
n_end_calo_el_data=kend_calo_el_data.Integral()
n_end_calo_inel_data=kend_calo_inel_data.Integral()

#Read MC histograms ----------------------------------
f_mc=RT.TFile(args.c, "OPEN")
f_mc.ls()

#truth
ke0_mc=f_mc.Get("h1d_ke0");
ke0_stop_mc=f_mc.Get("h1d_ke0_stop");

#beam
kebeam_mc=f_mc.Get("h1d_kebeam")
kebeam_stop_mc=f_mc.Get("h1d_kebeam_stop")

#ff-truth
keff_mc=f_mc.Get("h1d_keff");
keff_stop_mc=f_mc.Get("h1d_keff_stop");
keff_inel_mc=f_mc.Get("h1d_keff_inel");


#ff-hyper
kehy_mc=f_mc.Get("h1d_kehy")
kehy_stop_mc=f_mc.Get("h1d_kehy_stop")
kehy_inel_mc=f_mc.Get("h1d_kehy_inel")

#ff-hyper
kehy_mc=f_mc.Get("h1d_kehy")
kehy_stop_mc=f_mc.Get("h1d_kehy_stop")
kehy_inel_mc=f_mc.Get("h1d_kehy_inel")

#E-dept
kerange_stop_mc=f_mc.Get("h1d_kerange_stop")
kecalo_stop_mc=f_mc.Get("h1d_kecalo_stop")

#ke-end[calo]
kend_calo_stop_mc=f_mc.Get("h1d_kend_calo_stop")
kend_calo_el_mc=f_mc.Get("h1d_kend_calo_el")
kend_calo_inel_mc=f_mc.Get("h1d_kend_calo_inel")

#ke-end[bb]
kend_bb_stop_mc=f_mc.Get("h1d_kend_bb_stop")
kend_bb_el_mc=f_mc.Get("h1d_kend_bb_el")
kend_bb_inel_mc=f_mc.Get("h1d_kend_bb_inel")

#ke-end[truth]
kend_true_stop_mc=f_mc.Get("h1d_kend_true_stop");
kend_true_el_mc=f_mc.Get("h1d_kend_true_el");
kend_true_inel_mc=f_mc.Get("h1d_kend_true_inel");

n_beam_stop_mc=kebeam_stop_mc.Integral()
n_0_stop_mc=ke0_stop_mc.Integral()

n_ff_stop_mc=keff_stop_mc.Integral()
n_ff_inel_mc=keff_inel_mc.Integral()

n_hy_stop_mc=kehy_stop_mc.Integral()
n_hy_inel_mc=kehy_inel_mc.Integral()

n_range_stop_mc=kerange_stop_mc.Integral()
n_calo_stop_mc=kecalo_stop_mc.Integral()

n_end_true_stop_mc=kend_true_stop_mc.Integral()
n_end_true_el_mc=kend_true_el_mc.Integral()
n_end_true_inel_mc=kend_true_inel_mc.Integral()

n_end_calo_stop_mc=kend_calo_stop_mc.Integral()
n_end_calo_el_mc=kend_calo_el_mc.Integral()
n_end_calo_inel_mc=kend_calo_inel_mc.Integral()

n_end_bb_stop_mc=kend_bb_stop_mc.Integral()
n_end_bb_el_mc=kend_bb_el_mc.Integral()
n_end_bb_inel_mc=kend_bb_inel_mc.Integral()

#normalization [mc]
kebeam_stop_mc.Scale(n_beam_stop_data/n_beam_stop_mc)

ke0_stop_mc.Scale(n_beam_stop_data/n_0_stop_mc)
keff_stop_mc.Scale(n_hy_stop_data/n_ff_stop_mc)
keff_inel_mc.Scale(n_beam_stop_data/n_ff_inel_mc)

kehy_stop_mc.Scale(n_hy_stop_data/n_hy_stop_mc)
kehy_inel_mc.Scale(n_hy_inel_data/n_hy_inel_mc)

kerange_stop_mc.Scale(n_range_stop_data/n_range_stop_mc)
kecalo_stop_mc.Scale(n_calo_stop_data/n_calo_stop_mc)

kend_true_stop_mc.Scale(n_end_calo_stop_data/n_end_true_stop_mc)
kend_true_el_mc.Scale(n_end_calo_el_data/n_end_true_el_mc)
kend_true_inel_mc.Scale(n_end_calo_inel_data/n_end_true_inel_mc)

kend_calo_stop_mc.Scale(n_end_calo_stop_data/n_end_calo_stop_mc)
kend_calo_el_mc.Scale(n_end_calo_el_data/n_end_calo_el_mc)
kend_calo_inel_mc.Scale(n_end_calo_inel_data/n_end_calo_inel_mc)

kend_bb_stop_mc.Scale(n_end_bb_stop_data/n_end_bb_stop_mc)
kend_bb_el_mc.Scale(n_end_bb_el_data/n_end_bb_el_mc)
kend_bb_inel_mc.Scale(n_end_bb_inel_data/n_end_bb_inel_mc)


#fitting to extract mu and sigma ----------------------------
fit_kebeam_stop_data=VNFit(kebeam_stop_data, 430, 3)
fit_kebeam_stop_data.SetName("fit_kebeam_stop_data")
fit_kebeam_stop_data.SetLineStyle(2)
m_stop_data=fit_kebeam_stop_data.GetParameter(0)
err_m_stop_data=fit_kebeam_stop_data.GetParError(0)
s_stop_data=fit_kebeam_stop_data.GetParameter(1)
err_s_stop_data=fit_kebeam_stop_data.GetParError(1)

fit_kebeam_stop_mc=VNFit(kebeam_stop_mc, 430, 3)
fit_kebeam_stop_mc.SetName("fit_kebeam_stop_mc")
fit_kebeam_stop_mc.SetLineStyle(2)
m_stop_mc=fit_kebeam_stop_mc.GetParameter(0)
err_m_stop_mc=fit_kebeam_stop_mc.GetParError(0)
s_stop_mc=fit_kebeam_stop_mc.GetParameter(1)
err_s_stop_mc=fit_kebeam_stop_mc.GetParError(1)

fit_ke0_stop_mc=VNFit(ke0_stop_mc, 430, 3)
fit_ke0_stop_mc.SetName("fit_ke0_stop_mc")
fit_ke0_stop_mc.SetLineStyle(2)
m0_stop_mc=fit_ke0_stop_mc.GetParameter(0)
err_m0_stop_mc=fit_ke0_stop_mc.GetParError(0)
s0_stop_mc=fit_ke0_stop_mc.GetParameter(1)
err_s0_stop_mc=fit_ke0_stop_mc.GetParError(1)


#Plots
c0=RT.TCanvas("c0","",1200,900)
c0.Divide(1,1)
c0.cd(1)

f2d_ke=RT.TH2D("f2d_ke","", 370, 250, 620, 600, 0, kebeam_stop_data.GetBinContent(kebeam_stop_data.GetMaximumBin())+300)
f2d_ke.SetTitle("Stopping Protons;Initial Proton Kinetic Energy [MeV];")
f2d_ke.GetXaxis().CenterTitle()
f2d_ke.Draw()
kebeam_stop_data.SetLineColor(1)
fit_kebeam_stop_data.SetLineColor(1)
fit_ke0_stop_mc.SetLineColor(3)

kebeam_stop_mc.SetLineColor(2)
ke0_stop_mc.SetLineColor(3)

fit_kebeam_stop_mc.Draw("same")
kebeam_stop_mc.Draw("hist same")

fit_ke0_stop_mc.Draw("same")
ke0_stop_mc.Draw("hist same")

kebeam_stop_data.Draw("ep same")
fit_kebeam_stop_data.Draw("same")

leg0=RT.TLegend(0.14,0.65,.9,0.9)
leg0.SetFillStyle(0)
txt0=[]
txt0.append("Data: #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(m_stop_data,err_m_stop_data,s_stop_data,err_s_stop_data))
txt0.append("MC(truth): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(m0_stop_mc,err_m0_stop_mc,s0_stop_mc,err_s0_stop_mc))
txt0.append("MC(spec): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(m_stop_mc,err_m_stop_mc,s_stop_mc,err_s_stop_mc))

print(txt0[0])
leg0.AddEntry(kebeam_stop_data, txt0[0], "ep")
leg0.AddEntry(ke0_stop_mc, txt0[1], "l")
leg0.AddEntry(kebeam_stop_mc, txt0[2], "l")
#leg0.SetNColumns(2);

leg0.Draw()
c0.Print(args.o+'/ke_ini.eps')


