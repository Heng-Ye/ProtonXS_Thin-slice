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
#python plotDataMC_KE.py -d /dune/data2/users/hyliao/protonana/v09_39_01/KEHY/proton_beamxy_beammom_runAll.root -c ../mc_proton_beamxy_beammom_bmrw_hyper.root -drw /dune/data2/users/hyliao/protonana/v09_39_01/KEHY_BMRW/proton_beamxy_beammom_bmrw_runAll.root -crw ../mc_proton_beamxy_beammom_calo_afterkerw_old.root -o plots_beamxy_bmrw_calo_hyper

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
parser.add_argument("-drw", type=str, help='Data RW File', default = "")
parser.add_argument("-crw", type=str, help='MC RW File', default = "")
parser.add_argument("-o", type=str, help='Output folder', default = "")

args=parser.parse_args()

if not (args.d and args.o and args.c):
  print("--> Please provide all files (data, MC, MR_rw, output_file_folder). Thank you!")
  exit()

if (args.d): print('Read data: '+args.d)
if (args.c): print('MC: '+args.c)
if (args.drw): print('Read RW data: '+args.drw)
if (args.crw): print('MC RW: '+args.crw)
if (args.o): print('Output folder: '+args.o)

#Const. E-loss Values -------------------------------------------
const_eloss_mc=47.0058/1.00097
const_eloss_data=45.6084/0.99943

#Read data histograms ----------------------------------
f_data=RT.TFile(args.d, "OPEN")
#f_data.ls()

#beam
kebeam_data=f_data.Get("h1d_kebeam")
kebeam_stop_data=f_data.Get("h1d_kebeam_stop")

#ff-const E-loss from kebeam
keffbeam_data=f_data.Get("h1d_keffbeam")
keffbeam_stop_data=f_data.Get("h1d_keffbeam_stop")
keffbeam_inel_data=f_data.Get("h1d_keffbeam_inel")
keffbeam_el_data=f_data.Get("h1d_keffbeam_el")

#ff
kehy_data=f_data.Get("h1d_kehy")
kehy_stop_data=f_data.Get("h1d_kehy_stop")
kehy_inel_data=f_data.Get("h1d_kehy_inel")
kehy_el_data=f_data.Get("h1d_kehy_el")

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

n_ffbeam_stop_data=keffbeam_stop_data.Integral()
n_ffbeam_inel_data=keffbeam_inel_data.Integral()

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


#Read data rw histograms ----------------------------------
f_data_rw=RT.TFile(args.drw, "OPEN")
#f_data.ls()

#beam
kebeam_data_rw=f_data_rw.Get("h1d_kebeam")
kebeam_stop_data_rw=f_data_rw.Get("h1d_kebeam_stop")

#ff-const E-loss from kebeam
keffbeam_data_rw=f_data_rw.Get("h1d_keffbeam")
keffbeam_stop_data_rw=f_data_rw.Get("h1d_keffbeam_stop")
keffbeam_inel_data_rw=f_data_rw.Get("h1d_keffbeam_inel")
keffbeam_el_data_rw=f_data_rw.Get("h1d_keffbeam_el")

#ff
kehy_data_rw=f_data_rw.Get("h1d_kehy")
kehy_stop_data_rw=f_data_rw.Get("h1d_kehy_stop")
kehy_inel_data_rw=f_data_rw.Get("h1d_kehy_inel")
kehy_el_data_rw=f_data_rw.Get("h1d_kehy_el")

#E-dept
kerange_stop_data_rw=f_data_rw.Get("h1d_kerange_stop")
kecalo_stop_data_rw=f_data_rw.Get("h1d_kecalo_stop")

#ke-end[calo]
kend_calo_stop_data_rw=f_data_rw.Get("h1d_kend_calo_stop")
kend_calo_el_data_rw=f_data_rw.Get("h1d_kend_calo_el")
kend_calo_inel_data_rw=f_data_rw.Get("h1d_kend_calo_inel")

#ke-end[bb]
kend_bb_stop_data_rw=f_data_rw.Get("h1d_kend_bb_stop")
kend_bb_el_data_rw=f_data_rw.Get("h1d_kend_bb_el")
kend_bb_inel_data_rw=f_data_rw.Get("h1d_kend_bb_inel")

#get entries 
n_beam_stop_data_rw=kebeam_stop_data_rw.Integral()

n_ffbeam_stop_data_rw=keffbeam_stop_data_rw.Integral()
n_ffbeam_inel_data_rw=keffbeam_inel_data_rw.Integral()

n_hy_stop_data_rw=kehy_stop_data_rw.Integral()
n_hy_inel_data_rw=kehy_inel_data_rw.Integral()

n_range_stop_data_rw=kerange_stop_data_rw.Integral()
n_calo_stop_data_rw=kecalo_stop_data_rw.Integral()

n_end_bb_stop_data_rw=kend_bb_stop_data_rw.Integral()
n_end_bb_el_data_rw=kend_bb_el_data_rw.Integral()
n_end_bb_inel_data_rw=kend_bb_inel_data_rw.Integral()

n_end_calo_stop_data_rw=kend_calo_stop_data_rw.Integral()
n_end_calo_el_data_rw=kend_calo_el_data_rw.Integral()
n_end_calo_inel_data_rw=kend_calo_inel_data_rw.Integral()

n_ffbeam_el_data_rw=keffbeam_el_data_rw.Integral()


#Read MC histograms ----------------------------------
f_mc=RT.TFile(args.c, "OPEN")
#f_mc.ls()

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
#keff_el_mc=f_mc.Get("h1d_keff_el");

#ff-const E-loss from kebeam
keffbeam_mc=f_mc.Get("h1d_keffbeam")
keffbeam_stop_mc=f_mc.Get("h1d_keffbeam_stop")
keffbeam_inel_mc=f_mc.Get("h1d_keffbeam_inel")
keffbeam_el_mc=f_mc.Get("h1d_keffbeam_el")







#ff-hyper
kehy_mc=f_mc.Get("h1d_kehy")
kehy_stop_mc=f_mc.Get("h1d_kehy_stop")
kehy_inel_mc=f_mc.Get("h1d_kehy_inel")
#kehy_el_mc=f_mc.Get("h1d_kehy_el")

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

n_ffbeam_stop_mc=keffbeam_stop_mc.Integral()
n_ffbeam_inel_mc=keffbeam_inel_mc.Integral()

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


#Read MC histograms ----------------------------------
f_mc_rw=RT.TFile(args.crw, "OPEN")
#f_mc_rw.ls()

#truth
ke0_mc_rw=f_mc_rw.Get("h1d_ke0");
ke0_stop_mc_rw=f_mc_rw.Get("h1d_ke0_stop");

#beam
kebeam_mc_rw=f_mc_rw.Get("h1d_kebeam")
kebeam_stop_mc_rw=f_mc_rw.Get("h1d_kebeam_stop")

#ff-truth
keff_mc_rw=f_mc_rw.Get("h1d_keff");
keff_stop_mc_rw=f_mc_rw.Get("h1d_keff_stop");
#keff_inel_mc_rw=f_mc_rw.Get("h1d_keff_inel");

#ff-const E-loss from kebeam
#keffbeam_mc_rw=f_mc_rw.Get("h1d_keffbeam")
keffbeam_stop_mc_rw=f_mc_rw.Get("h1d_keffbeam_stop")
keffbeam_el_mc_rw=f_mc_rw.Get("h1d_keffbeam_el")
keffbeam_el_inel_mc_rw=f_mc_rw.Get("h1d_keffbeam_el_inel")
keffbeam_el_el_mc_rw=f_mc_rw.Get("h1d_keffbeam_el_el")
keffbeam_el_midcosmic_mc_rw=f_mc_rw.Get("h1d_keffbeam_el_midcosmic")
keffbeam_el_midpi_mc_rw=f_mc_rw.Get("h1d_keffbeam_el_midpi")
keffbeam_el_midp_mc_rw=f_mc_rw.Get("h1d_keffbeam_el_midp")
keffbeam_el_midmu_mc_rw=f_mc_rw.Get("h1d_keffbeam_el_midmu")
keffbeam_el_mideg_mc_rw=f_mc_rw.Get("h1d_keffbeam_el_mideg")
keffbeam_el_midother_mc_rw=f_mc_rw.Get("h1d_keffbeam_el_midother")

keffbeam_inel_mc_rw=f_mc_rw.Get("h1d_keffbeam_inel")
keffbeam_inel_inel_mc_rw=f_mc_rw.Get("h1d_keffbeam_inel_inel")
keffbeam_inel_el_mc_rw=f_mc_rw.Get("h1d_keffbeam_inel_el")
keffbeam_inel_midcosmic_mc_rw=f_mc_rw.Get("h1d_keffbeam_inel_midcosmic")
keffbeam_inel_midpi_mc_rw=f_mc_rw.Get("h1d_keffbeam_inel_midpi")
keffbeam_inel_midp_mc_rw=f_mc_rw.Get("h1d_keffbeam_inel_midp")
keffbeam_inel_midmu_mc_rw=f_mc_rw.Get("h1d_keffbeam_inel_midmu")
keffbeam_inel_mideg_mc_rw=f_mc_rw.Get("h1d_keffbeam_inel_mideg")
keffbeam_inel_midother_mc_rw=f_mc_rw.Get("h1d_keffbeam_inel_midother")


#ff-hyper
kehy_mc_rw=f_mc_rw.Get("h1d_kehy")
kehy_stop_mc_rw=f_mc_rw.Get("h1d_kehy_stop")
kehy_inel_mc_rw=f_mc_rw.Get("h1d_kehy_inel")
#kehy_el_mc_rw=f_mc_rw.Get("h1d_kehy_el")

#E-dept
kerange_stop_mc_rw=f_mc_rw.Get("h1d_kerange_stop")
kecalo_stop_mc_rw=f_mc_rw.Get("h1d_kecalo_stop")

#ke-end[calo]
kend_calo_stop_mc_rw=f_mc_rw.Get("h1d_kend_calo_stop")
kend_calo_el_mc_rw=f_mc_rw.Get("h1d_kend_calo_el")
kend_calo_el_el_mc_rw=f_mc_rw.Get("h1d_kend_calo_el_el")
kend_calo_el_inel_mc_rw=f_mc_rw.Get("h1d_kend_calo_el_inel")
kend_calo_el_midcosmic_mc_rw=f_mc_rw.Get("h1d_kend_calo_el_midcosmic")
kend_calo_el_midpi_mc_rw=f_mc_rw.Get("h1d_kend_calo_el_midpi")
kend_calo_el_midp_mc_rw=f_mc_rw.Get("h1d_kend_calo_el_midp")
kend_calo_el_midmu_mc_rw=f_mc_rw.Get("h1d_kend_calo_el_midmu")
kend_calo_el_mideg_mc_rw=f_mc_rw.Get("h1d_kend_calo_el_mideg")
kend_calo_el_midother_mc_rw=f_mc_rw.Get("h1d_kend_calo_el_midother")

kend_calo_inel_mc_rw=f_mc_rw.Get("h1d_kend_calo_inel")

#ke-end[bb]
kend_bb_stop_mc_rw=f_mc_rw.Get("h1d_kend_bb_stop")
kend_bb_el_mc_rw=f_mc_rw.Get("h1d_kend_bb_el")
kend_bb_inel_mc_rw=f_mc_rw.Get("h1d_kend_bb_inel")

#ke-end[truth]
kend_true_stop_mc_rw=f_mc_rw.Get("h1d_kend_true_stop");
kend_true_el_mc_rw=f_mc_rw.Get("h1d_kend_true_el");
kend_true_inel_mc_rw=f_mc_rw.Get("h1d_kend_true_inel");

n_beam_stop_mc_rw=kebeam_stop_mc_rw.Integral()
n_0_stop_mc_rw=ke0_stop_mc_rw.Integral()

n_ffbeam_stop_mc_rw=keffbeam_stop_mc_rw.Integral()
n_ffbeam_inel_mc_rw=keffbeam_inel_mc_rw.Integral()

n_ff_stop_mc_rw=keff_stop_mc_rw.Integral()
#n_ff_inel_mc_rw=keff_inel_mc_rw.Integral()

n_hy_stop_mc_rw=kehy_stop_mc_rw.Integral()
n_hy_inel_mc_rw=kehy_inel_mc_rw.Integral()

n_range_stop_mc_rw=kerange_stop_mc_rw.Integral()
n_calo_stop_mc_rw=kecalo_stop_mc_rw.Integral()

n_end_true_stop_mc_rw=kend_true_stop_mc_rw.Integral()
n_end_true_el_mc_rw=kend_true_el_mc_rw.Integral()
n_end_true_inel_mc_rw=kend_true_inel_mc_rw.Integral()

n_end_calo_stop_mc_rw=kend_calo_stop_mc_rw.Integral()
n_end_calo_el_mc_rw=kend_calo_el_mc_rw.Integral()
n_end_calo_inel_mc_rw=kend_calo_inel_mc_rw.Integral()

n_end_bb_stop_mc_rw=kend_bb_stop_mc_rw.Integral()
n_end_bb_el_mc_rw=kend_bb_el_mc_rw.Integral()
n_end_bb_inel_mc_rw=kend_bb_inel_mc_rw.Integral()

n_ffbeam_el_mc_rw=keffbeam_el_mc_rw.Integral()



#normalization [mc] ------------------------------------------------
kebeam_stop_mc.Scale(n_beam_stop_data/n_beam_stop_mc)

ke0_stop_mc.Scale(n_beam_stop_data/n_0_stop_mc)

keff_stop_mc.Scale(n_hy_stop_data/n_ff_stop_mc)
keff_inel_mc.Scale(n_beam_stop_data/n_ff_inel_mc)

keffbeam_stop_mc.Scale(n_ffbeam_stop_data/n_ffbeam_stop_mc)
keffbeam_inel_mc.Scale(n_ffbeam_inel_data/n_ffbeam_inel_mc)

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

keffbeam_stop_mc.Scale(n_hy_stop_data/n_ffbeam_stop_mc)
keffbeam_inel_mc.Scale(n_hy_inel_data/n_ffbeam_inel_mc)

#rw case
kebeam_stop_mc_rw.Scale(n_beam_stop_data/n_beam_stop_mc_rw)

ke0_stop_mc_rw.Scale(n_beam_stop_data/n_0_stop_mc_rw)

keff_stop_mc_rw.Scale(n_hy_stop_data/n_ff_stop_mc_rw)
#keff_inel_mc_rw.Scale(n_beam_stop_data/n_ff_inel_mc_rw)

keffbeam_stop_mc_rw.Scale(n_ffbeam_stop_data_rw/n_ffbeam_stop_mc_rw)
keffbeam_inel_mc_rw.Scale(n_ffbeam_inel_data_rw/n_ffbeam_inel_mc_rw)

kehy_stop_mc_rw.Scale(n_hy_stop_data/n_hy_stop_mc_rw)
kehy_inel_mc_rw.Scale(n_hy_inel_data/n_hy_inel_mc_rw)

kerange_stop_mc_rw.Scale(n_range_stop_data/n_range_stop_mc_rw)
kecalo_stop_mc_rw.Scale(n_calo_stop_data/n_calo_stop_mc_rw)

kend_true_stop_mc_rw.Scale(n_end_calo_stop_data/n_end_true_stop_mc_rw)
kend_true_el_mc_rw.Scale(n_end_calo_el_data/n_end_true_el_mc_rw)
kend_true_inel_mc_rw.Scale(n_end_calo_inel_data/n_end_true_inel_mc_rw)

kend_calo_stop_mc_rw.Scale(n_end_calo_stop_data/n_end_calo_stop_mc_rw)
kend_calo_el_mc_rw.Scale(n_end_calo_el_data/n_end_calo_el_mc_rw)
kend_calo_el_el_mc_rw.Scale(n_end_calo_el_data/n_end_calo_el_mc_rw)
kend_calo_el_inel_mc_rw.Scale(n_end_calo_el_data/n_end_calo_el_mc_rw)
kend_calo_el_midcosmic_mc_rw.Scale(n_end_calo_el_data/n_end_calo_el_mc_rw)
kend_calo_el_midpi_mc_rw.Scale(n_end_calo_el_data/n_end_calo_el_mc_rw)
kend_calo_el_midp_mc_rw.Scale(n_end_calo_el_data/n_end_calo_el_mc_rw)
kend_calo_el_midmu_mc_rw.Scale(n_end_calo_el_data/n_end_calo_el_mc_rw)
kend_calo_el_mideg_mc_rw.Scale(n_end_calo_el_data/n_end_calo_el_mc_rw)
kend_calo_el_midother_mc_rw.Scale(n_end_calo_el_data/n_end_calo_el_mc_rw)

kend_calo_inel_mc_rw.Scale(n_end_calo_inel_data/n_end_calo_inel_mc_rw)

kend_bb_stop_mc_rw.Scale(n_end_bb_stop_data/n_end_bb_stop_mc_rw)
kend_bb_el_mc_rw.Scale(n_end_bb_el_data/n_end_bb_el_mc_rw)
kend_bb_inel_mc_rw.Scale(n_end_bb_inel_data/n_end_bb_inel_mc_rw)

keffbeam_el_mc_rw.Scale(n_ffbeam_el_data_rw/n_ffbeam_el_mc_rw)
keffbeam_el_inel_mc_rw.Scale(n_ffbeam_el_data_rw/n_ffbeam_el_mc_rw)
keffbeam_el_el_mc_rw.Scale(n_ffbeam_el_data_rw/n_ffbeam_el_mc_rw)
keffbeam_el_midcosmic_mc_rw.Scale(n_ffbeam_el_data_rw/n_ffbeam_el_mc_rw)
keffbeam_el_midpi_mc_rw.Scale(n_ffbeam_el_data_rw/n_ffbeam_el_mc_rw)
keffbeam_el_midp_mc_rw.Scale(n_ffbeam_el_data_rw/n_ffbeam_el_mc_rw)
keffbeam_el_midmu_mc_rw.Scale(n_ffbeam_el_data_rw/n_ffbeam_el_mc_rw)
keffbeam_el_mideg_mc_rw.Scale(n_ffbeam_el_data_rw/n_ffbeam_el_mc_rw)
keffbeam_el_midother_mc_rw.Scale(n_ffbeam_el_data_rw/n_ffbeam_el_mc_rw)


keffbeam_inel_mc_rw.Scale(n_ffbeam_inel_data_rw/n_ffbeam_inel_mc_rw)
keffbeam_inel_inel_mc_rw.Scale(n_ffbeam_inel_data_rw/n_ffbeam_inel_mc_rw)
keffbeam_inel_el_mc_rw.Scale(n_ffbeam_inel_data_rw/n_ffbeam_inel_mc_rw)
keffbeam_inel_midcosmic_mc_rw.Scale(n_ffbeam_inel_data_rw/n_ffbeam_inel_mc_rw)
keffbeam_inel_midpi_mc_rw.Scale(n_ffbeam_inel_data_rw/n_ffbeam_inel_mc_rw)
keffbeam_inel_midp_mc_rw.Scale(n_ffbeam_inel_data_rw/n_ffbeam_inel_mc_rw)
keffbeam_inel_midmu_mc_rw.Scale(n_ffbeam_inel_data_rw/n_ffbeam_inel_mc_rw)
keffbeam_inel_mideg_mc_rw.Scale(n_ffbeam_inel_data_rw/n_ffbeam_inel_mc_rw)
keffbeam_inel_midother_mc_rw.Scale(n_ffbeam_inel_data_rw/n_ffbeam_inel_mc_rw)


#fitting to extract mu and sigma ----------------------------
#Before rw
#Data
#initial KE
fit_kebeam_stop_data=VNFit(kebeam_stop_data, 430, 3)

#FF_constEloss
fit_keffbeam_stop_data=VNFit(keffbeam_stop_data, 410, 3)
fit_keffbeam_inel_data=VNFit(keffbeam_inel_data, 410, 3)
fit_keffbeam_el_data=VNFit(keffbeam_el_data, 410, 3)

#FF
fit_kehy_stop_data=VNFit(kehy_stop_data, 410, 3)
fit_kehy_inel_data=VNFit(kehy_inel_data, 410, 3)
fit_kehy_el_data=VNFit(kehy_el_data, 410, 3)

#EDept
fit_kerange_stop_data=VNFit(kerange_stop_data, 400, 3)
fit_kecalo_stop_data=VNFit(kecalo_stop_data, 400, 3)

#End_Calo
fit_kend_calo_stop_data=VNFit(kend_calo_stop_data, 20, 3)
fit_kend_calo_el_data=VNFit(kend_calo_el_data, 20, 3)
fit_kend_calo_inel_data=VNFit(kend_calo_inel_data, 20, 3)

#End_bb
fit_kend_bb_stop_data=VNFit(kend_bb_stop_data, 20, 3)
fit_kend_bb_el_data=VNFit(kend_bb_el_data, 20, 3)
fit_kend_bb_inel_data=VNFit(kend_bb_inel_data, 20, 3)

#MC
#initial KE
fit_ke0_stop_mc=VNFit(ke0_stop_mc, 430, 3)
fit_kebeam_stop_mc=VNFit(kebeam_stop_mc, 430, 3)

#FF-truth
fit_keff_stop_mc=VNFit(keff_stop_mc, 430, 3)
fit_keff_inel_mc=VNFit(keff_inel_mc,430, 3)
#fit_keff_el_mc=VNFit(keff_el_mc,430, 3)

#FF-const E-loss
fit_keffbeam_mc=VNFit(keffbeam_mc, 400, 3)
fit_keffbeam_stop_mc=VNFit(keffbeam_stop_mc, 400, 3)
fit_keffbeam_inel_mc=VNFit(keffbeam_inel_mc, 400, 3)
fit_keffbeam_el_mc=VNFit(keffbeam_el_mc, 400, 3)

#FF-HY
fit_kehy_stop_mc=VNFit(kehy_stop_mc, 400, 3)
fit_kehy_inel_mc=VNFit(kehy_inel_mc, 400, 3)
#fit_kehy_el_mc=VNFit(kehy_el_mc, 400, 3)

#EDept
fit_kerange_stop_mc=VNFit(kerange_stop_mc, 420, 3)
fit_kecalo_stop_mc=VNFit(kecalo_stop_mc, 420, 3)

#End_true
fit_kend_true_stop_mc=VNFit(kend_true_stop_mc, 0, 1)
fit_kend_true_el_mc=VNFit(kend_true_el_mc, 20, 3)
fit_kend_true_inel_mc=VNFit(kend_true_inel_mc, 20, 3)

#End_calo
fit_kend_calo_stop_mc=VNFit(kend_calo_stop_mc , 20, 3)
fit_kend_calo_el_mc=VNFit(kend_calo_el_mc , 20, 3)
fit_kend_calo_inel_mc=VNFit(kend_calo_inel_mc , 20, 3)
fit_kend_bb_stop_mc=VNFit(kend_bb_stop_mc , 20, 3)
fit_kend_bb_el_mc=VNFit(kend_bb_el_mc , 20, 3)
fit_kend_bb_inel_mc=VNFit(kend_bb_inel_mc , 20, 3)

#after rw
#Data
#initial KE
fit_kebeam_stop_data_rw=VNFit(kebeam_stop_data_rw, 430, 3)

#FF_constEloss
fit_keffbeam_stop_data_rw=VNFit(keffbeam_stop_data_rw, 410, 3)
fit_keffbeam_inel_data_rw=VNFit(keffbeam_inel_data_rw, 410, 3)
fit_keffbeam_el_data_rw=VNFit(keffbeam_el_data_rw, 410, 3)

#FF
fit_kehy_stop_data_rw=VNFit(kehy_stop_data_rw, 410, 3)
fit_kehy_inel_data_rw=VNFit(kehy_inel_data_rw, 410, 3)
fit_kehy_el_data_rw=VNFit(kehy_el_data_rw, 410, 3)

#EDept
fit_kerange_stop_data_rw=VNFit(kerange_stop_data_rw, 400, 3)
fit_kecalo_stop_data_rw=VNFit(kecalo_stop_data_rw, 400, 3)

#End_Calo
fit_kend_calo_stop_data_rw=VNFit(kend_calo_stop_data_rw, 20, 3)
fit_kend_calo_el_data_rw=VNFit(kend_calo_el_data_rw, 20, 3)
fit_kend_calo_inel_data_rw=VNFit(kend_calo_inel_data_rw, 20, 3)

#End_bb
fit_kend_bb_stop_data_rw=VNFit(kend_bb_stop_data_rw, 20, 3)
fit_kend_bb_el_data_rw=VNFit(kend_bb_el_data_rw, 20, 3)
fit_kend_bb_inel_data_rw=VNFit(kend_bb_inel_data_rw, 20, 3)

#MC
#initial KE
fit_ke0_stop_mc_rw=VNFit(ke0_stop_mc_rw, 430, 3)
fit_kebeam_stop_mc_rw=VNFit(kebeam_stop_mc_rw, 430, 3)

#FF-truth
fit_keff_stop_mc_rw=VNFit(keff_stop_mc_rw, 430, 3)

#FF-const E-loss
#fit_keffbeam_mc_rw=VNFit(keffbeam_mc_rw, 400, 3)
fit_keffbeam_stop_mc_rw=VNFit(keffbeam_stop_mc_rw, 400, 3)
fit_keffbeam_inel_mc_rw=VNFit(keffbeam_inel_mc_rw, 400, 3)
fit_keffbeam_el_mc_rw=VNFit(keffbeam_el_mc_rw, 400, 3)

#FF-HY
fit_kehy_stop_mc_rw=VNFit(kehy_stop_mc_rw, 400, 3)
fit_kehy_inel_mc_rw=VNFit(kehy_inel_mc_rw, 400, 3)
#fit_kehy_el_mc_rw=VNFit(kehy_el_mc_rw, 400, 3)

#EDept
fit_kerange_stop_mc_rw=VNFit(kerange_stop_mc_rw, 420, 3)
fit_kecalo_stop_mc_rw=VNFit(kecalo_stop_mc_rw, 420, 3)

#End_true
fit_kend_true_stop_mc_rw=VNFit(kend_true_stop_mc_rw, 0, 1)
fit_kend_true_el_mc_rw=VNFit(kend_true_el_mc_rw, 20, 3)
fit_kend_true_inel_mc_rw=VNFit(kend_true_inel_mc_rw, 20, 3)

#End_calo
fit_kend_calo_stop_mc_rw=VNFit(kend_calo_stop_mc_rw , 20, 3)
fit_kend_calo_el_mc_rw=VNFit(kend_calo_el_mc_rw , 20, 3)
fit_kend_calo_inel_mc_rw=VNFit(kend_calo_inel_mc_rw , 20, 3)
fit_kend_bb_stop_mc_rw=VNFit(kend_bb_stop_mc_rw , 20, 3)
fit_kend_bb_el_mc_rw=VNFit(kend_bb_el_mc_rw , 20, 3)
fit_kend_bb_inel_mc_rw=VNFit(kend_bb_inel_mc_rw , 20, 3)


#get fitted mu and sigma [mc] ------------------------
#stopping protons
n_mc_stop=11

mu_stop_mc=[]
er_mu_stop_mc=[]
sigma_stop_mc=[]
er_sigma_stop_mc=[]

for i in range(n_mc_stop):
  fit=fit_ke0_stop_mc
  if i==0: 
    fit=fit_ke0_stop_mc
  if i==1: 
    fit=fit_kebeam_stop_mc
  if i==2:
    fit=fit_keff_stop_mc
  if i==3:
    fit=fit_keffbeam_stop_mc
  if i==4:
    fit=fit_kehy_stop_mc
  if i==5:
    fit=fit_kerange_stop_mc
  if i==6:
    fit=fit_kecalo_stop_mc
  if i==7:
    fit=fit_kend_calo_stop_mc
  if i==8:
    fit=fit_kend_bb_stop_mc
  if i==9:
    fit=fit_kend_true_stop_mc
  if i==10:
    fit=fit_keffbeam_stop_mc_rw

  m=fit.GetParameter(0)
  er_m=fit.GetParError(0)
  s=fit.GetParameter(1)
  er_s=fit.GetParError(1)
  #mu_sg_mc.append([m, er_m, s, er_s]) 
  mu_stop_mc.append(m)
  er_mu_stop_mc.append(er_m)
  sigma_stop_mc.append(s)
  er_sigma_stop_mc.append(er_s)
  print("i=",i," m=",m, "s=",s,"\n")
    
#print("n_mc_stop=",n_mc_stop,"\n")
#print("mu_sg_mc::",mu_sg_mc)
#print("len(mu_sg_mc):",mu_sg_mc)
#print("len(mu_sg_mc[0]):",len(mu_sg_mc[0]))
#print("const_eloss_mc:",const_eloss_mc)


#el protons(mc)
n_mc_el_rw=2

mu_el_mc_rw=[]
er_mu_el_mc_rw=[]
sigma_el_mc_rw=[]
er_sigma_el_mc_rw=[]

for i in range(n_mc_el_rw):
  fit=fit_keffbeam_el_mc_rw
  if i==0: 
    fit=fit_keffbeam_el_mc_rw
  if i==1: 
    fit=fit_kend_calo_el_mc_rw
  m=fit.GetParameter(0)
  er_m=fit.GetParError(0)
  s=fit.GetParameter(1)
  er_s=fit.GetParError(1)
  mu_el_mc_rw.append(m)
  er_mu_el_mc_rw.append(er_m)
  sigma_el_mc_rw.append(s)
  er_sigma_el_mc_rw.append(er_s)
    
#inel protons(mc)
n_mc_inel_rw=1

mu_inel_mc_rw=[]
er_mu_inel_mc_rw=[]
sigma_inel_mc_rw=[]
er_sigma_inel_mc_rw=[]

for i in range(n_mc_inel_rw):
  fit=fit_keffbeam_inel_mc_rw
  if i==0: 
    fit=fit_keffbeam_inel_mc_rw
  m=fit.GetParameter(0)
  er_m=fit.GetParError(0)
  s=fit.GetParameter(1)
  er_s=fit.GetParError(1)
  mu_inel_mc_rw.append(m)
  er_mu_inel_mc_rw.append(er_m)
  sigma_inel_mc_rw.append(s)
  er_sigma_inel_mc_rw.append(er_s)

	
#get fitted mu and sigma [data] -----------------
n_data_stop=8

mu_stop_data=[]
er_mu_stop_data=[]
sigma_stop_data=[]
er_sigma_stop_data=[]

for i in range(n_mc_stop):
  fit=fit_kebeam_stop_data
  if i==0: 
    fit=fit_kebeam_stop_data
  if i==1: 
    fit=fit_keffbeam_stop_data
  if i==2:
    fit=fit_kehy_stop_data
  if i==3:
    fit=fit_kerange_stop_data
  if i==4:
    fit=fit_kecalo_stop_data
  if i==5:
    fit=fit_kend_calo_stop_data
  if i==6:
    fit=fit_kend_bb_stop_data
  if i==7: 
    fit=fit_keffbeam_stop_data_rw

  m=fit.GetParameter(0)
  er_m=fit.GetParError(0)
  s=fit.GetParameter(1)
  er_s=fit.GetParError(1)
  mu_stop_data.append(m)
  er_mu_stop_data.append(er_m)
  sigma_stop_data.append(s)
  er_sigma_stop_data.append(er_s)



#el protons(data)
n_data_el_rw=2

mu_el_data_rw=[]
er_mu_el_data_rw=[]
sigma_el_data_rw=[]
er_sigma_el_data_rw=[]

for i in range(n_data_el_rw):
  fit=fit_keffbeam_el_data_rw
  if i==0: 
    fit=fit_keffbeam_el_data_rw
  if i==1: 
    fit=fit_kend_calo_el_data_rw
  m=fit.GetParameter(0)
  er_m=fit.GetParError(0)
  s=fit.GetParameter(1)
  er_s=fit.GetParError(1)
  mu_el_data_rw.append(m)
  er_mu_el_data_rw.append(er_m)
  sigma_el_data_rw.append(s)
  er_sigma_el_data_rw.append(er_s)


#el protons(data)
n_data_inel_rw=1

mu_inel_data_rw=[]
er_mu_inel_data_rw=[]
sigma_inel_data_rw=[]
er_sigma_inel_data_rw=[]

for i in range(n_data_inel_rw):
  fit=fit_keffbeam_inel_data_rw
  if i==0: 
    fit=fit_keffbeam_inel_data_rw
  m=fit.GetParameter(0)
  er_m=fit.GetParError(0)
  s=fit.GetParameter(1)
  er_s=fit.GetParError(1)
  mu_inel_data_rw.append(m)
  er_mu_inel_data_rw.append(er_m)
  sigma_inel_data_rw.append(s)
  er_sigma_inel_data_rw.append(er_s)


#Plots #########################################################################################################################################################
#KEbeam ------------------------------------------------------------------------------------------------------------------
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

kebeam_stop_data.SetLineWidth(1)
kebeam_stop_mc.SetLineWidth(1)
ke0_stop_mc.SetLineWidth(1)

fit_kebeam_stop_mc.Draw("same")
kebeam_stop_mc.Draw("hist same")

fit_ke0_stop_mc.Draw("same")
ke0_stop_mc.Draw("hist same")

kebeam_stop_data.Draw("ep same")
fit_kebeam_stop_data.Draw("same")

leg0=RT.TLegend(0.14,0.65,.9,0.9)
leg0.SetFillStyle(0)
txt0=[]
txt0.append("Data: #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_data[0],er_mu_stop_data[0],sigma_stop_data[0],er_sigma_stop_data[0]))
txt0.append("MC(truth): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[0],er_mu_stop_mc[0],sigma_stop_mc[0],er_sigma_stop_mc[0]))
txt0.append("MC(spec): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[1],er_mu_stop_mc[1],sigma_stop_mc[1],er_sigma_stop_mc[1]))

print(txt0[0])
leg0.AddEntry(kebeam_stop_data, txt0[0], "ep")
leg0.AddEntry(ke0_stop_mc, txt0[1], "l")
leg0.AddEntry(kebeam_stop_mc, txt0[2], "l")
#leg0.SetNColumns(2);

leg0.Draw()
c0.Print(args.o+'/ke_ini.eps')

#KEFF --------------------------------------------------------------------------------------------------------------------------
c1=RT.TCanvas("c1","",1200,900)
c1.Divide(1,1)
c1.cd(1)

f2d_keff=RT.TH2D("f2d_keff","", 370, 220, 600, 600, 0, keffbeam_stop_data.GetBinContent(keffbeam_stop_data.GetMaximumBin())+400)
f2d_keff.SetTitle("Stopping Protons;Proton Kinetic Energy at TPC FF [MeV];")
f2d_keff.GetXaxis().CenterTitle()
f2d_keff.Draw()

keffbeam_stop_data.SetLineColor(1)
fit_keffbeam_stop_data.SetLineColor(1)

keff_stop_mc.SetLineColor(3)
fit_keff_stop_mc.SetLineColor(3)

keffbeam_stop_mc.SetLineColor(2)
fit_keffbeam_stop_mc.SetLineColor(2)

kehy_stop_data.SetLineColor(4)
kehy_stop_data.SetMarkerColor(4)
fit_kehy_stop_data.SetLineColor(4)

kehy_stop_mc.SetLineColor(6)
kehy_stop_mc.SetMarkerColor(6)
fit_kehy_stop_mc.SetLineColor(6)

#kebeff_stop_mc.SetLineColor(3)

keff_stop_mc.SetLineWidth(1)
keffbeam_stop_mc.SetLineWidth(1)
kehy_stop_data.SetLineWidth(1)
kehy_stop_mc.SetLineWidth(1)
keffbeam_stop_data.SetLineWidth(1)


fit_keff_stop_mc.Draw("same")
keff_stop_mc.Draw("hist same")

fit_keffbeam_stop_mc.Draw("same")
keffbeam_stop_mc.Draw("hist same")

fit_kehy_stop_data.Draw("same")
kehy_stop_data.Draw("ep same")

fit_kehy_stop_mc.Draw("same")
kehy_stop_mc.Draw("hist same")

keffbeam_stop_data.Draw("ep same")
fit_keffbeam_stop_data.Draw("same")

leg1=RT.TLegend(0.14,0.65,.9,0.9)
leg1.SetFillStyle(0)
txt1=[]
txt1.append("Data(spec): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_data[1],er_mu_stop_data[1],sigma_stop_data[1],er_sigma_stop_data[1]))
txt1.append("Data(fit): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_data[2],er_mu_stop_data[2],sigma_stop_data[2],er_sigma_stop_data[2]))
txt1.append("MC(truth): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[2],er_mu_stop_mc[2],sigma_stop_mc[2],er_sigma_stop_mc[2]))
txt1.append("MC(spec): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[3],er_mu_stop_mc[3],sigma_stop_mc[3],er_sigma_stop_mc[3]))
txt1.append("MC(fit): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[4],er_mu_stop_mc[4],sigma_stop_mc[4],er_sigma_stop_mc[4]))

print(txt1[0])
leg1.AddEntry(keffbeam_stop_data, txt1[0], "ep")
leg1.AddEntry(kehy_stop_data, txt1[1], "ep")
leg1.AddEntry(keff_stop_mc, txt1[2], "l")
leg1.AddEntry(keffbeam_stop_mc, txt1[3], "l")
leg1.AddEntry(kehy_stop_mc, txt1[4], "l")

#leg1.SetNColumns(2);

leg1.Draw()
c1.Print(args.o+'/ke_ff.eps')


#Ecalo --------------------------------------------------------------------------------------------------------------------------
c2=RT.TCanvas("c2","",1200,900)
c2.Divide(1,1)
c2.cd(1)

f2d_kecalo=RT.TH2D("f2d_kecalo","", 500, 100, 600, 600, 0, kecalo_stop_data.GetBinContent(kecalo_stop_data.GetMaximumBin())+400)
f2d_kecalo.SetTitle("Stopping Protons;Proton Kinetic Energy [MeV];")
f2d_kecalo.GetXaxis().CenterTitle()
f2d_kecalo.Draw()

kecalo_stop_data.SetLineColor(1)
fit_kecalo_stop_data.SetLineColor(1)

kecalo_stop_mc.SetLineColor(2)
fit_kecalo_stop_mc.SetLineColor(2)

kecalo_stop_mc.SetLineWidth(1)
kecalo_stop_data.SetLineWidth(1)

fit_kecalo_stop_mc.Draw("same")
kecalo_stop_mc.Draw("hist same")
fit_kecalo_stop_data.Draw("same")
kecalo_stop_data.Draw("ep same")

leg2=RT.TLegend(0.14,0.65,.9,0.9)
leg2.SetFillStyle(0)
txt2=[]
txt2.append("Data(calo): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_data[4],er_mu_stop_data[4],sigma_stop_data[4],er_sigma_stop_data[4]))
txt2.append("MC(calo): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[6],er_mu_stop_mc[6],sigma_stop_mc[6],er_sigma_stop_mc[6]))

print(txt2[0])
leg2.AddEntry(kecalo_stop_data, txt2[0], "ep")
leg2.AddEntry(kecalo_stop_mc, txt2[1], "l")

#leg1.SetNColumns(2);

leg2.Draw()
c2.Print(args.o+'/ke_calo.eps')


#Erange --------------------------------------------------------------------------------------------------------------------------
c3=RT.TCanvas("c3","",1200,900)
c3.Divide(1,1)
c3.cd(1)

f2d_kerange=RT.TH2D("f2d_kerange","", 500, 100, 600, 600, 0, kerange_stop_data.GetBinContent(kerange_stop_data.GetMaximumBin())+400)
f2d_kerange.SetTitle("Stopping Protons;Proton Kinetic Energy [MeV];")
f2d_kerange.GetXaxis().CenterTitle()
f2d_kerange.Draw()

kerange_stop_data.SetLineColor(1)
fit_kerange_stop_data.SetLineColor(1)

kerange_stop_mc.SetLineColor(2)
fit_kerange_stop_mc.SetLineColor(2)

kerange_stop_mc.SetLineWidth(1)
kerange_stop_data.SetLineWidth(1)

fit_kerange_stop_data.Draw("same")
fit_kerange_stop_mc.Draw("same")
kerange_stop_mc.Draw("hist same")
kerange_stop_data.Draw("ep same")

leg3=RT.TLegend(0.14,0.65,.9,0.9)
leg3.SetFillStyle(0)
txt3=[]
txt3.append("Data(range): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_data[3],er_mu_stop_data[3],sigma_stop_data[3],er_sigma_stop_data[3]))
txt3.append("MC(range): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[6],er_mu_stop_mc[5],sigma_stop_mc[5],er_sigma_stop_mc[5]))

print(txt3[0])
leg3.AddEntry(kerange_stop_data, txt2[0], "ep")
leg3.AddEntry(kerange_stop_mc, txt2[1], "l")

leg3.Draw()
c3.Print(args.o+'/ke_range.eps')

#KEend_calo --------------------------------------------------------------------------------------------------------------------------
#stop
c4=RT.TCanvas("c4","",1200,900)
c4.Divide(1,1)
c4.cd(1)

f2d_kend_calo=RT.TH2D("f2d_kend_calo","", 450, -150, 300, 600, 0, kend_calo_stop_data.GetBinContent(kend_calo_stop_data.GetMaximumBin())+800)
f2d_kend_calo.SetTitle("Stopping Protons;Proton Kinetic Energy at Track End [MeV];")
f2d_kend_calo.GetXaxis().CenterTitle()
f2d_kend_calo.Draw()

kend_calo_stop_data.SetLineColor(1)
fit_kend_calo_stop_data.SetLineColor(1)

kend_calo_stop_mc.SetLineColor(2)
fit_kend_calo_stop_mc.SetLineColor(2)

kend_true_stop_mc.SetLineColor(3)
fit_kend_true_stop_mc.SetLineColor(3)

kend_true_stop_mc.SetLineWidth(1)
kend_calo_stop_mc.SetLineWidth(1)
kend_calo_stop_data.SetLineWidth(1)

fit_kend_calo_stop_data.Draw("same")
fit_kend_calo_stop_mc.Draw("same")
fit_kend_true_stop_mc.Draw("same")
kend_true_stop_mc.Draw("hist same")
kend_calo_stop_mc.Draw("hist same")
kend_calo_stop_data.Draw("ep same")

leg4=RT.TLegend(0.14,0.65,.9,0.9)
leg4.SetFillStyle(0)
txt4=[]
txt4.append("Data(calo): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_data[5],er_mu_stop_data[5],sigma_stop_data[5],er_sigma_stop_data[5]))
txt4.append("MC(calo): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[7],er_mu_stop_mc[7],sigma_stop_mc[7],er_sigma_stop_mc[7]))
txt4.append("MC(truth): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[9],er_mu_stop_mc[9],sigma_stop_mc[9],er_sigma_stop_mc[9]))

print(txt4[0])
leg4.AddEntry(kend_calo_stop_data, txt4[0], "ep")
leg4.AddEntry(kend_calo_stop_mc, txt4[1], "l")
leg4.AddEntry(kend_true_stop_mc, txt4[2], "l")

leg4.Draw()
c4.Print(args.o+'/kend_calo_stop.eps')






#el_rw
c4_el=RT.TCanvas("c4_el","",1200,900)
c4_el.Divide(1,1)
c4_el.cd(1)

f2d_kend_calo_el=RT.TH2D("f2d_kend_calo_el","", 420, -120, 300, 600, 0, kend_calo_el_data.GetBinContent(kend_calo_el_data.GetMaximumBin())+250)
f2d_kend_calo_el.SetTitle("Elastic-scattering Protons;Proton Kinetic Energy at Track End [MeV];")
f2d_kend_calo_el.GetXaxis().CenterTitle()
f2d_kend_calo_el.Draw()

kend_calo_el_data.SetLineColor(1)
fit_kend_calo_el_data.SetLineColor(1)
kend_calo_el_data_rw.SetLineColor(4)
fit_kend_calo_el_data.SetLineColor(4)

kend_calo_el_mc.SetLineColor(2)
fit_kend_calo_el_mc.SetLineColor(2)
kend_calo_el_mc_rw.SetLineColor(7)
fit_kend_calo_el_mc_rw.SetLineColor(7)

kend_true_el_mc.SetLineColor(3)
fit_kend_true_el_mc.SetLineColor(3)

kend_true_el_mc.SetLineWidth(1)
kend_calo_el_mc_rw.SetLineWidth(1)
kend_calo_el_data_rw.SetLineWidth(1)


#fit_kend_calo_el_data.Draw("same")
#fit_kend_calo_el_mc.Draw("same")
#fit_kend_true_el_mc.Draw("same")
kend_true_el_mc.Draw("hist same")
#kend_calo_el_mc.Draw("hist same")
#kend_calo_el_data.Draw("ep same")
kend_calo_el_mc_rw.Draw("hist same")
kend_calo_el_data_rw.Draw("ep same")

leg4_el=RT.TLegend(0.54,0.65,.9,0.9)
leg4_el.SetFillStyle(0)
txt4_el=[]
#txt4_el.append("Data(calo)".format())
txt4_el.append("Data(calo)[rw]".format())
#txt4_el.append("MC(calo)".format())
txt4_el.append("MC(calo)[rw]".format())
txt4_el.append("MC(truth)".format())

print(txt4_el[0])
#leg4_el.AddEntry(kend_calo_el_data, txt4_el[0], "ep")
leg4_el.AddEntry(kend_calo_el_data_rw, txt4_el[0], "ep")
#leg4_el.AddEntry(kend_calo_el_mc, txt4_el[2], "l")
leg4_el.AddEntry(kend_calo_el_mc_rw, txt4_el[1], "l")
leg4_el.AddEntry(kend_true_el_mc, txt4_el[2], "l")

leg4_el.Draw()
c4_el.Print(args.o+'/kend_calo_el_rw.eps')



#inel
c4_inel=RT.TCanvas("c4_inel","",1200,900)
c4_inel.Divide(1,1)
c4_inel.cd(1)

f2d_kend_calo_inel=RT.TH2D("f2d_kend_calo_inel","", 750, -150, 600, 600, 0, kend_calo_inel_data.GetBinContent(kend_calo_inel_data.GetMaximumBin())+250)
f2d_kend_calo_inel.SetTitle("Inelastic-scattering Protons;Proton Kinetic Energy at Track End [MeV];")
f2d_kend_calo_inel.GetXaxis().CenterTitle()
f2d_kend_calo_inel.Draw()

kend_calo_inel_data.SetLineColor(1)
fit_kend_calo_inel_data.SetLineColor(1)

kend_calo_inel_mc.SetLineColor(2)
fit_kend_calo_inel_mc.SetLineColor(2)

kend_true_inel_mc.SetLineColor(3)
fit_kend_true_inel_mc.SetLineColor(3)

#fit_kend_calo_inel_data.Draw("same")
#fit_kend_calo_inel_mc.Draw("same")
#fit_kend_true_inel_mc.Draw("same")
kend_true_inel_mc.Draw("hist same")
kend_calo_inel_mc.Draw("hist same")
kend_calo_inel_data.Draw("ep same")

leg4_inel=RT.TLegend(0.54,0.65,.9,0.9)
leg4_inel.SetFillStyle(0)
txt4_inel=[]
txt4_inel.append("Data(calo)".format())
txt4_inel.append("MC(calo)".format())
txt4_inel.append("MC(truth)".format())

print(txt4_inel[0])
leg4_inel.AddEntry(kend_calo_inel_data, txt4_inel[0], "ep")
leg4_inel.AddEntry(kend_calo_inel_mc, txt4_inel[1], "l")
leg4_inel.AddEntry(kend_true_inel_mc, txt4_inel[2], "l")

leg4_inel.Draw()
c4_inel.Print(args.o+'/kend_calo_inel.eps')


#inel_rw
c4_inel=RT.TCanvas("c4_inel","",1200,900)
c4_inel.Divide(1,1)
c4_inel.cd(1)

f2d_kend_calo_inel=RT.TH2D("f2d_kend_calo_inel","", 420, -150, 600, 600, 0, kend_calo_inel_data.GetBinContent(kend_calo_inel_data.GetMaximumBin())+250)
f2d_kend_calo_inel.SetTitle("Elastic-scattering Protons;Proton Kinetic Energy at Track End [MeV];")
f2d_kend_calo_inel.GetXaxis().CenterTitle()
f2d_kend_calo_inel.Draw()

kend_calo_inel_data.SetLineColor(1)
fit_kend_calo_inel_data.SetLineColor(1)
kend_calo_inel_data_rw.SetLineColor(4)
fit_kend_calo_inel_data.SetLineColor(4)

kend_calo_inel_mc.SetLineColor(2)
fit_kend_calo_inel_mc.SetLineColor(2)
kend_calo_inel_mc_rw.SetLineColor(7)
fit_kend_calo_inel_mc_rw.SetLineColor(7)

kend_true_inel_mc.SetLineColor(3)
fit_kend_true_inel_mc.SetLineColor(3)

#fit_kend_calo_inel_data.Draw("same")
#fit_kend_calo_inel_mc.Draw("same")
#fit_kend_true_inel_mc.Draw("same")
kend_true_inel_mc.Draw("hist same")
#kend_calo_inel_mc.Draw("hist same")
#kend_calo_inel_data.Draw("ep same")
kend_calo_inel_mc_rw.Draw("hist same")
kend_calo_inel_data_rw.Draw("ep same")

leg4_inel=RT.TLegend(0.54,0.65,.9,0.9)
leg4_inel.SetFillStyle(0)
txt4_inel=[]
#txt4_inel.append("Data(calo)".format())
txt4_inel.append("Data(calo)[rw]".format())
#txt4_inel.append("MC(calo)".format())
txt4_inel.append("MC(calo)[rw]".format())
txt4_inel.append("MC(truth)".format())

print(txt4_inel[0])
#leg4_inel.AddEntry(kend_calo_inel_data, txt4_inel[0], "ep")
leg4_inel.AddEntry(kend_calo_inel_data_rw, txt4_inel[0], "ep")
#leg4_inel.AddEntry(kend_calo_inel_mc, txt4_inel[2], "l")
leg4_inel.AddEntry(kend_calo_inel_mc_rw, txt4_inel[1], "l")
leg4_inel.AddEntry(kend_true_inel_mc, txt4_inel[2], "l")

leg4_inel.Draw()
c4_inel.Print(args.o+'/kend_calo_inel_rw.eps')



#KEend_bb --------------------------------------------------------------------------------------------------------------------------
c5=RT.TCanvas("c5","",1200,900)
c5.Divide(1,1)
c5.cd(1)

f2d_kend_bb=RT.TH2D("f2d_kend_bb","", 450, -150, 300, 600, 0, kend_bb_stop_data.GetBinContent(kend_bb_stop_data.GetMaximumBin())+800)
f2d_kend_bb.SetTitle("Stopping Protons;Proton Kinetic Energy at Track End [MeV];")
f2d_kend_bb.GetXaxis().CenterTitle()
f2d_kend_bb.Draw()

kend_bb_stop_data.SetLineColor(1)
fit_kend_bb_stop_data.SetLineColor(1)

kend_bb_stop_mc.SetLineColor(2)
fit_kend_bb_stop_mc.SetLineColor(2)

kend_true_stop_mc.SetLineColor(3)
fit_kend_true_stop_mc.SetLineColor(3)

fit_kend_bb_stop_data.Draw("same")
fit_kend_bb_stop_mc.Draw("same")
fit_kend_true_stop_mc.Draw("same")
kend_true_stop_mc.Draw("hist same")
kend_bb_stop_mc.Draw("hist same")
kend_bb_stop_data.Draw("ep same")

leg5=RT.TLegend(0.14,0.65,.9,0.9)
leg5.SetFillStyle(0)
txt5=[]
txt5.append("Data(bb): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_data[6],er_mu_stop_data[6],sigma_stop_data[6],er_sigma_stop_data[6]))
txt5.append("MC(bb): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[8],er_mu_stop_mc[8],sigma_stop_mc[8],er_sigma_stop_mc[8]))
txt5.append("MC(truth): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[9],er_mu_stop_mc[9],sigma_stop_mc[9],er_sigma_stop_mc[9]))

print(txt5[0])
leg5.AddEntry(kend_bb_stop_data, txt5[0], "ep")
leg5.AddEntry(kend_bb_stop_mc, txt5[1], "l")
leg5.AddEntry(kend_true_stop_mc, txt5[2], "l")

leg5.Draw()
c5.Print(args.o+'/kend_bb.eps')

#inel
c5_inel=RT.TCanvas("c5_inel","",1200,900)
c5_inel.Divide(1,1)
c5_inel.cd(1)

f2d_kend_bb_inel=RT.TH2D("f2d_kend_bb_inel","", 750, -150, 600, 600, 0, kend_bb_inel_data.GetBinContent(kend_bb_inel_data.GetMaximumBin())+250)
f2d_kend_bb_inel.SetTitle("Inelastic-scattering Protons;Proton Kinetic Energy at Track End [MeV];")
f2d_kend_bb_inel.GetXaxis().CenterTitle()
f2d_kend_bb_inel.Draw()

kend_bb_inel_data.SetLineColor(1)
fit_kend_bb_inel_data.SetLineColor(1)

kend_bb_inel_mc.SetLineColor(2)
fit_kend_bb_inel_mc.SetLineColor(2)

kend_true_inel_mc.SetLineColor(3)
fit_kend_true_inel_mc.SetLineColor(3)

#fit_kend_bb_inel_data.Draw("same")
#fit_kend_bb_inel_mc.Draw("same")
#fit_kend_true_inel_mc.Draw("same")
kend_true_inel_mc.Draw("hist same")
kend_bb_inel_mc.Draw("hist same")
kend_bb_inel_data.Draw("ep same")

leg5_inel=RT.TLegend(0.54,0.65,.9,0.9)
leg5_inel.SetFillStyle(0)
txt5_inel=[]
txt5_inel.append("Data(bb)".format())
txt5_inel.append("MC(bb)".format())
txt5_inel.append("MC(truth)".format())

print(txt5_inel[0])
leg5_inel.AddEntry(kend_bb_inel_data, txt5_inel[0], "ep")
leg5_inel.AddEntry(kend_bb_inel_mc, txt5_inel[1], "l")
leg5_inel.AddEntry(kend_true_inel_mc, txt5_inel[2], "l")

leg5_inel.Draw()
c5_inel.Print(args.o+'/kend_bb_inel.eps')

#el_rw
c5_el=RT.TCanvas("c5_el","",1200,900)
c5_el.Divide(1,1)
c5_el.cd(1)

f2d_kend_bb_el=RT.TH2D("f2d_kend_bb_el","", 420, -120, 300, 600, 0, kend_bb_el_data.GetBinContent(kend_bb_el_data.GetMaximumBin())+250)
f2d_kend_bb_el.SetTitle("Elastic-scattering Protons;Proton Kinetic Energy at Track End [MeV];")
f2d_kend_bb_el.GetXaxis().CenterTitle()
f2d_kend_bb_el.Draw()

kend_bb_el_data.SetLineColor(1)
fit_kend_bb_el_data.SetLineColor(1)
kend_bb_el_data_rw.SetLineColor(4)
kend_bb_el_data_rw.SetMarkerColor(4)
fit_kend_bb_el_data.SetLineColor(4)

kend_bb_el_mc.SetLineColor(2)
fit_kend_bb_el_mc.SetLineColor(2)
kend_bb_el_mc_rw.SetLineColor(7)
fit_kend_bb_el_mc_rw.SetLineColor(7)

kend_true_el_mc.SetLineColor(3)
fit_kend_true_el_mc.SetLineColor(3)

#fit_kend_bb_el_data.Draw("same")
#fit_kend_bb_el_mc.Draw("same")
#fit_kend_true_el_mc.Draw("same")
kend_true_el_mc.Draw("hist same")
#kend_bb_el_mc.Draw("hist same")
#kend_bb_el_data.Draw("ep same")
kend_bb_el_mc_rw.Draw("hist same")
kend_bb_el_data_rw.Draw("ep same")

leg5_el=RT.TLegend(0.54,0.65,.9,0.9)
leg5_el.SetFillStyle(0)
txt5_el=[]
#txt5_el.append("Data(bb)".format())
txt5_el.append("Data(bb)[rw]".format())
#txt5_el.append("MC(bb)".format())
txt5_el.append("MC(bb)[rw]".format())
txt5_el.append("MC(truth)".format())

print(txt5_el[0])
#leg5_el.AddEntry(kend_bb_el_data, txt5_el[0], "ep")
leg5_el.AddEntry(kend_bb_el_data_rw, txt5_el[0], "ep")
#leg5_el.AddEntry(kend_bb_el_mc, txt5_el[2], "l")
leg5_el.AddEntry(kend_bb_el_mc_rw, txt5_el[1], "l")
leg5_el.AddEntry(kend_true_el_mc, txt5_el[2], "l")

leg5_el.Draw()
c5_el.Print(args.o+'/kend_bb_el_rw.eps')


#inel_rw
c5_inel=RT.TCanvas("c5_inel","",1200,900)
c5_inel.Divide(1,1)
c5_inel.cd(1)

f2d_kend_bb_inel=RT.TH2D("f2d_kend_bb_inel","", 420, -110, 600, 600, 0, kend_bb_inel_data.GetBinContent(kend_bb_inel_data.GetMaximumBin())+250)
f2d_kend_bb_inel.SetTitle("Elastic-scattering Protons;Proton Kinetic Energy at Track End [MeV];")
f2d_kend_bb_inel.GetXaxis().CenterTitle()
f2d_kend_bb_inel.Draw()

kend_bb_inel_data.SetLineColor(1)
fit_kend_bb_inel_data.SetLineColor(1)
kend_bb_inel_data_rw.SetLineColor(4)
fit_kend_bb_inel_data.SetLineColor(4)

kend_bb_inel_mc.SetLineColor(2)
fit_kend_bb_inel_mc.SetLineColor(2)
kend_bb_inel_mc_rw.SetLineColor(7)
fit_kend_bb_inel_mc_rw.SetLineColor(7)

kend_true_inel_mc.SetLineColor(3)
fit_kend_true_inel_mc.SetLineColor(3)

#fit_kend_bb_inel_data.Draw("same")
#fit_kend_bb_inel_mc.Draw("same")
#fit_kend_true_inel_mc.Draw("same")
kend_true_inel_mc.Draw("hist same")
#kend_bb_inel_mc.Draw("hist same")
#kend_bb_inel_data.Draw("ep same")
kend_bb_inel_mc_rw.Draw("hist same")
kend_bb_inel_data_rw.Draw("ep same")

leg5_inel=RT.TLegend(0.54,0.65,.9,0.9)
leg5_inel.SetFillStyle(0)
txt5_inel=[]
#txt5_inel.append("Data(bb)".format())
txt5_inel.append("Data(bb)[rw]".format())
#txt5_inel.append("MC(bb)".format())
txt5_inel.append("MC(bb)[rw]".format())
txt5_inel.append("MC(truth)".format())

print(txt5_inel[0])
#leg5_inel.AddEntry(kend_bb_inel_data, txt5_inel[0], "ep")
leg5_inel.AddEntry(kend_bb_inel_data_rw, txt5_inel[0], "ep")
#leg5_inel.AddEntry(kend_bb_inel_mc, txt5_inel[2], "l")
leg5_inel.AddEntry(kend_bb_inel_mc_rw, txt5_inel[1], "l")
leg5_inel.AddEntry(kend_true_inel_mc, txt5_inel[2], "l")

leg5_inel.Draw()
c5_inel.Print(args.o+'/kend_bb_inel_rw.eps')



#inel_rw
c5x_inel=RT.TCanvas("c5x_inel","",1200,900)
c5x_inel.Divide(1,1)
c5x_inel.cd(1)

f2d_keffbeam_inel=RT.TH2D("f2d_keffbeam_inel","", 420, 200, 600, 600, 0, keffbeam_inel_data.GetBinContent(keffbeam_inel_data.GetMaximumBin())+250)
f2d_keffbeam_inel.SetTitle("Elastic-scattering Protons;Proton Kinetic Energy at Track End [MeV];")
f2d_keffbeam_inel.GetXaxis().CenterTitle()
f2d_keffbeam_inel.Draw()

keffbeam_inel_data.SetLineColor(1)
fit_keffbeam_inel_data.SetLineColor(1)
keffbeam_inel_data_rw.SetLineColor(4)
fit_keffbeam_inel_data.SetLineColor(4)

keffbeam_inel_mc.SetLineColor(2)
fit_keffbeam_inel_mc.SetLineColor(2)
keffbeam_inel_mc_rw.SetLineColor(7)
fit_keffbeam_inel_mc_rw.SetLineColor(7)

keff_inel_mc.SetLineColor(3)

#fit_keffbeam_inel_data.Draw("same")
#fit_keffbeam_inel_mc.Draw("same")
#fit_keff_true_inel_mc.Draw("same")
keff_inel_mc.Draw("hist same")
#keffbeam_inel_mc.Draw("hist same")
#keffbeam_inel_data.Draw("ep same")
keffbeam_inel_mc_rw.Draw("hist same")
keffbeam_inel_data_rw.Draw("ep same")

leg5_inel=RT.TLegend(0.54,0.65,.9,0.9)
leg5_inel.SetFillStyle(0)
txt5_inel=[]
#txt5_inel.append("Data(bb)".format())
txt5_inel.append("Data(bb)[rw]".format())
#txt5_inel.append("MC(bb)".format())
txt5_inel.append("MC(bb)[rw]".format())
txt5_inel.append("MC(truth)".format())

print(txt5_inel[0])
#leg5_inel.AddEntry(keffbeam_inel_data, txt5_inel[0], "ep")
leg5_inel.AddEntry(keffbeam_inel_data_rw, txt5_inel[0], "ep")
#leg5_inel.AddEntry(keffbeam_inel_mc, txt5_inel[2], "l")
leg5_inel.AddEntry(keffbeam_inel_mc_rw, txt5_inel[1], "l")
leg5_inel.AddEntry(keff_inel_mc, txt5_inel[2], "l")

leg5_inel.Draw()
c5x_inel.Print(args.o+'/keffbeam_inel_rw.eps')



#KEFFbeam --------------------------------------------------------------------------------------------------------------------------
#stop
c1_rw=RT.TCanvas("c1_rw","",1200,900)
c1_rw.Divide(1,1)
c1_rw.cd(1)

f2d_keff_rw=RT.TH2D("f2d_keff_rw","", 370, 220, 600, 600, 0, keffbeam_stop_data.GetBinContent(keffbeam_stop_data.GetMaximumBin())+400)
f2d_keff_rw.SetTitle("Stopping Protons;Proton Kinetic Energy at TPC FF [MeV];")
f2d_keff_rw.GetXaxis().CenterTitle()
f2d_keff_rw.Draw()

keffbeam_stop_data_rw.SetLineColor(8)
keffbeam_stop_data_rw.SetMarkerColor(8)
keffbeam_stop_data_rw.SetLineWidth(1)
fit_keffbeam_stop_data_rw.SetLineColor(8)

keffbeam_stop_mc_rw.SetLineColor(7)
keffbeam_stop_mc_rw.SetMarkerColor(7)
keffbeam_stop_mc_rw.SetLineWidth(1)
fit_keffbeam_stop_mc_rw.SetLineColor(7)

keff_stop_mc.SetLineWidth(1)
kehy_stop_mc.SetLineWidth(1)
kehy_stop_data.SetLineWidth(1)
keffbeam_stop_data.SetLineWidth(1)
keffbeam_stop_data_rw.SetLineWidth(1)
keffbeam_stop_mc.SetLineWidth(1)
keffbeam_stop_mc_rw.SetLineWidth(1)

fit_keff_stop_mc.Draw("same")
keff_stop_mc.Draw("hist same")

fit_keffbeam_stop_mc.Draw("same")
keffbeam_stop_mc.Draw("hist same")

fit_kehy_stop_data.Draw("same")
kehy_stop_data.Draw("ep same")

fit_kehy_stop_mc.Draw("same")
kehy_stop_mc.Draw("hist same")

fit_keffbeam_stop_data.Draw("same")
keffbeam_stop_data.Draw("ep same")

fit_keffbeam_stop_mc_rw.Draw("same")
keffbeam_stop_mc_rw.Draw("hist same")

fit_keffbeam_stop_data_rw.Draw("same")
keffbeam_stop_data_rw.Draw("ep same")


leg1_rw=RT.TLegend(0.14,0.65,.9,0.9)
leg1_rw.SetFillStyle(0)
txt1_rw=[]
txt1_rw.append("Data(spec): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_data[1],er_mu_stop_data[1],sigma_stop_data[1],er_sigma_stop_data[1]))
txt1_rw.append("Data(fit): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_data[2],er_mu_stop_data[2],sigma_stop_data[2],er_sigma_stop_data[2]))
txt1_rw.append("Data(rw): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_data[7],er_mu_stop_data[7],sigma_stop_data[7],er_sigma_stop_data[7]))
txt1_rw.append("MC(truth): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[2],er_mu_stop_mc[2],sigma_stop_mc[2],er_sigma_stop_mc[2]))
txt1_rw.append("MC(spec): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[3],er_mu_stop_mc[3],sigma_stop_mc[3],er_sigma_stop_mc[3]))
txt1_rw.append("MC(fit): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[4],er_mu_stop_mc[4],sigma_stop_mc[4],er_sigma_stop_mc[4]))
txt1_rw.append("MC(rw): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[10],er_mu_stop_mc[10],sigma_stop_mc[10],er_sigma_stop_mc[10]))

print(txt1_rw[0])
leg1_rw.AddEntry(keffbeam_stop_data, txt1_rw[0], "ep")
leg1_rw.AddEntry(kehy_stop_data, txt1_rw[1], "ep")
leg1_rw.AddEntry(keffbeam_stop_data_rw, txt1_rw[2], "l")
leg1_rw.AddEntry(keff_stop_mc, txt1_rw[3], "l")
leg1_rw.AddEntry(keffbeam_stop_mc, txt1_rw[4], "l")
leg1_rw.AddEntry(kehy_stop_mc, txt1_rw[5], "l")
leg1_rw.AddEntry(keffbeam_stop_mc_rw, txt1_rw[6], "l")

#leg1.SetNColumns(2);

leg1_rw.Draw()
c1_rw.Print(args.o+'/ke_ff_rw.eps')





#el -----------------------------------------------------------------------------------------------------------------------------------------
c1_el_rw=RT.TCanvas("c1_el_rw","",1200,900)
c1_el_rw.Divide(1,1)
c1_el_rw.cd(1)

f2d_keff_el_rw=RT.TH2D("f2d_keff_el_rw","", 370, 220, 600, 600, 0, keffbeam_el_data.GetBinContent(keffbeam_el_data.GetMaximumBin())+400)
f2d_keff_el_rw.SetTitle("Elastic-scattering Protons;Proton Kinetic Energy at TPC FF [MeV];")
f2d_keff_el_rw.GetXaxis().CenterTitle()
f2d_keff_el_rw.Draw()

keffbeam_el_data_rw.SetLineColor(1)
keffbeam_el_data_rw.SetMarkerColor(1)

keffbeam_el_data_rw.SetLineWidth(1)
fit_keffbeam_el_data_rw.SetLineColor(1)

keffbeam_el_mc_rw.SetLineColor(1)
keffbeam_el_mc_rw.SetMarkerColor(1)

keffbeam_el_mc_rw.SetLineWidth(1)
fit_keffbeam_el_mc_rw.SetLineColor(1)

keffbeam_el_inel_mc_rw.SetFillColor(2)
keffbeam_el_el_mc_rw.SetFillColor(4)
keffbeam_el_midcosmic_mc_rw.SetFillColor(5)
keffbeam_el_midpi_mc_rw.SetFillColor(6)
keffbeam_el_midp_mc_rw.SetFillColor(3)
keffbeam_el_midmu_mc_rw.SetFillColor(28)
keffbeam_el_mideg_mc_rw.SetFillColor(30)
keffbeam_el_midother_mc_rw.SetFillColor(15)

keffbeam_el_inel_mc_rw.SetLineColor(2)
keffbeam_el_el_mc_rw.SetLineColor(4)
keffbeam_el_midcosmic_mc_rw.SetLineColor(5)
keffbeam_el_midpi_mc_rw.SetLineColor(6)
keffbeam_el_midp_mc_rw.SetLineColor(3)
keffbeam_el_midmu_mc_rw.SetLineColor(28)
keffbeam_el_mideg_mc_rw.SetLineColor(30)
keffbeam_el_midother_mc_rw.SetLineColor(15)

hs_keffbeam_el_rw=RT.THStack("hs_keffbeam_el_rw","");
hs_keffbeam_el_rw.Add(keffbeam_el_inel_mc_rw)
hs_keffbeam_el_rw.Add(keffbeam_el_el_mc_rw)
hs_keffbeam_el_rw.Add(keffbeam_el_midp_mc_rw)
hs_keffbeam_el_rw.Add(keffbeam_el_midcosmic_mc_rw)
hs_keffbeam_el_rw.Add(keffbeam_el_midpi_mc_rw)
hs_keffbeam_el_rw.Add(keffbeam_el_midmu_mc_rw)
hs_keffbeam_el_rw.Add(keffbeam_el_mideg_mc_rw)
hs_keffbeam_el_rw.Add(keffbeam_el_midother_mc_rw)

keffbeam_el_mc_rw.SetLineWidth(1)
keffbeam_el_data_rw.SetLineWidth(1)

hs_keffbeam_el_rw.Draw("hist same")
#fit_keffbeam_stop_mc_rw.Draw("same")
keffbeam_el_mc_rw.Draw("hist same")
#fit_keffbeam_stop_data_rw.Draw("same")
keffbeam_el_data_rw.Draw("ep same")


tleg11_rw=RT.TLegend(0.14,0.75,.9,0.87)
tleg11_rw.SetFillStyle(0)
txt11_rw=[]
txt11_rw.append("Data(rw): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_el_data_rw[0],er_mu_el_data_rw[0],sigma_el_data_rw[0],er_sigma_el_data_rw[0]))
txt11_rw.append("MC(rw): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_el_mc_rw[0],er_mu_el_mc_rw[0],sigma_el_mc_rw[0],er_sigma_el_mc_rw[0]))
tleg11_rw.AddEntry(keffbeam_el_data_rw, txt11_rw[0], "ep")
tleg11_rw.AddEntry(0, "", "")
tleg11_rw.AddEntry(keffbeam_el_mc_rw, txt11_rw[1], "l")
tleg11_rw.AddEntry(0, "", "")
tleg11_rw.SetNColumns(2);


leg11_rw=RT.TLegend(0.14,0.62,.88,0.73)
leg11_rw.SetFillStyle(0)
leg11_rw.AddEntry(keffbeam_el_inel_mc_rw, "Inel", "f")
leg11_rw.AddEntry(keffbeam_el_el_mc_rw, "El", "f")
leg11_rw.AddEntry(keffbeam_el_midcosmic_mc_rw, "misID:cosmic", "f")
leg11_rw.AddEntry(keffbeam_el_midp_mc_rw, "misID:p", "f")
leg11_rw.AddEntry(keffbeam_el_midpi_mc_rw, "misID:#pi", "f")
leg11_rw.AddEntry(keffbeam_el_midmu_mc_rw, "misID:#mu", "f")
leg11_rw.AddEntry(keffbeam_el_mideg_mc_rw, "misID:e/#gamma", "f")
leg11_rw.AddEntry(keffbeam_el_midother_mc_rw, "misID:other", "f")
leg11_rw.SetNColumns(4);

tleg11_rw.Draw()
leg11_rw.Draw()

c1_el_rw.Print(args.o+'/ke_ff_el_rw.eps')


#KEffbeam_el: data/mc ratio ---------------------------------------------------
R_keffbeam_el=keffbeam_el_data_rw.Clone()
R_keffbeam_el.Divide(R_keffbeam_el, keffbeam_el_mc_rw)

r1_el_rw=RT.TCanvas("r1_el_rw","",1200,900)
r1_el_rw.Divide(1,1)
r1_el_rw.cd(1)

r2d_keff_el_rw=RT.TH2D("r2d_keff_el_rw","", 370, 220, 600, 10, 0, 3)
r2d_keff_el_rw.SetTitle("Elastic-scattering Protons;Proton Kinetic Energy at TPC FF [MeV]; Data/MC")
r2d_keff_el_rw.GetXaxis().CenterTitle()
r2d_keff_el_rw.Draw("")
R_keffbeam_el.Draw("ep same")

line_keff_el_rw=RT.TLine(220, 1, 600, 1)
line_keff_el_rw.SetLineColor(1)
line_keff_el_rw.SetLineStyle(2)
line_keff_el_rw.Draw("same")


r1_el_rw.Print(args.o+'/ratio_ff_el_rw.eps')








#inel:keff -----------------------------------------------------------------------------------------------------------------------------------------
c1_inel_rw=RT.TCanvas("c1_inel_rw","",1200,900)
c1_inel_rw.Divide(1,1)
c1_inel_rw.cd(1)

f2d_keff_inel_rw=RT.TH2D("f2d_keff_inel_rw","", 370, 220, 600, 600, 0, keffbeam_inel_data.GetBinContent(keffbeam_inel_data.GetMaximumBin())+400)
f2d_keff_inel_rw.SetTitle("Inelastic-scattering Protons;Proton Kinetic Energy at TPC FF [MeV];")
f2d_keff_inel_rw.GetXaxis().CenterTitle()
f2d_keff_inel_rw.Draw()

keffbeam_inel_data_rw.SetLineColor(1)
keffbeam_inel_data_rw.SetMarkerColor(1)

keffbeam_inel_data_rw.SetLineWidth(1)
fit_keffbeam_inel_data_rw.SetLineColor(1)

keffbeam_inel_mc_rw.SetLineColor(1)
keffbeam_inel_mc_rw.SetMarkerColor(1)

keffbeam_inel_mc_rw.SetLineWidth(1)
fit_keffbeam_inel_mc_rw.SetLineColor(1)

keffbeam_inel_inel_mc_rw.SetFillColor(2)
keffbeam_inel_el_mc_rw.SetFillColor(4)
keffbeam_inel_midcosmic_mc_rw.SetFillColor(5)
keffbeam_inel_midpi_mc_rw.SetFillColor(6)
keffbeam_inel_midp_mc_rw.SetFillColor(3)
keffbeam_inel_midmu_mc_rw.SetFillColor(28)
keffbeam_inel_mideg_mc_rw.SetFillColor(30)
keffbeam_inel_midother_mc_rw.SetFillColor(15)

keffbeam_inel_inel_mc_rw.SetLineColor(2)
keffbeam_inel_el_mc_rw.SetLineColor(4)
keffbeam_inel_midcosmic_mc_rw.SetLineColor(5)
keffbeam_inel_midpi_mc_rw.SetLineColor(6)
keffbeam_inel_midp_mc_rw.SetLineColor(3)
keffbeam_inel_midmu_mc_rw.SetLineColor(28)
keffbeam_inel_mideg_mc_rw.SetLineColor(30)
keffbeam_inel_midother_mc_rw.SetLineColor(15)

hs_keffbeam_inel_rw=RT.THStack("hs_keffbeam_inel_rw","");
hs_keffbeam_inel_rw.Add(keffbeam_inel_inel_mc_rw)
hs_keffbeam_inel_rw.Add(keffbeam_inel_el_mc_rw)
hs_keffbeam_inel_rw.Add(keffbeam_inel_midp_mc_rw)
hs_keffbeam_inel_rw.Add(keffbeam_inel_midcosmic_mc_rw)
hs_keffbeam_inel_rw.Add(keffbeam_inel_midpi_mc_rw)
hs_keffbeam_inel_rw.Add(keffbeam_inel_midmu_mc_rw)
hs_keffbeam_inel_rw.Add(keffbeam_inel_mideg_mc_rw)
hs_keffbeam_inel_rw.Add(keffbeam_inel_midother_mc_rw)

keffbeam_inel_mc_rw.SetLineWidth(1)
keffbeam_inel_data_rw.SetLineWidth(1)

hs_keffbeam_inel_rw.Draw("hist same")
#fit_keffbeam_stop_mc_rw.Draw("same")
keffbeam_inel_mc_rw.Draw("hist same")
#fit_keffbeam_stop_data_rw.Draw("same")
keffbeam_inel_data_rw.Draw("ep same")


itleg11_rw=RT.TLegend(0.14,0.75,.9,0.87)
itleg11_rw.SetFillStyle(0)
itxt11_rw=[]
itxt11_rw.append("Data(rw): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_inel_data_rw[0],er_mu_inel_data_rw[0],sigma_inel_data_rw[0],er_sigma_inel_data_rw[0]))
itxt11_rw.append("MC(rw): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_inel_mc_rw[0],er_mu_inel_mc_rw[0],sigma_inel_mc_rw[0],er_sigma_inel_mc_rw[0]))
itleg11_rw.AddEntry(keffbeam_inel_data_rw, itxt11_rw[0], "ep")
itleg11_rw.AddEntry(0, "", "")
itleg11_rw.AddEntry(keffbeam_inel_mc_rw, itxt11_rw[1], "l")
itleg11_rw.AddEntry(0, "", "")
itleg11_rw.SetNColumns(2);


ileg11_rw=RT.TLegend(0.14,0.62,.88,0.73)
ileg11_rw.SetFillStyle(0)
ileg11_rw.AddEntry(keffbeam_inel_inel_mc_rw, "Inel", "f")
ileg11_rw.AddEntry(keffbeam_inel_el_mc_rw, "El", "f")
ileg11_rw.AddEntry(keffbeam_inel_midcosmic_mc_rw, "misID:cosmic", "f")
ileg11_rw.AddEntry(keffbeam_inel_midp_mc_rw, "misID:p", "f")
ileg11_rw.AddEntry(keffbeam_inel_midpi_mc_rw, "misID:#pi", "f")
ileg11_rw.AddEntry(keffbeam_inel_midmu_mc_rw, "misID:#mu", "f")
ileg11_rw.AddEntry(keffbeam_inel_mideg_mc_rw, "misID:e/#gamma", "f")
ileg11_rw.AddEntry(keffbeam_inel_midother_mc_rw, "misID:other", "f")
ileg11_rw.SetNColumns(4);

itleg11_rw.Draw()
ileg11_rw.Draw()

c1_inel_rw.Print(args.o+'/ke_ff_inel_rw.eps')


#KEffbeam_el: data/mc ratio ---------------------------------------------------
R_keffbeam_inel=keffbeam_inel_data_rw.Clone()
R_keffbeam_inel.Divide(R_keffbeam_inel, keffbeam_inel_mc_rw)

r1_inel_rw=RT.TCanvas("r1_inel_rw","",1200,900)
r1_inel_rw.Divide(1,1)
r1_inel_rw.cd(1)

r2d_keff_inel_rw=RT.TH2D("r2d_keff_inel_rw","", 370, 220, 600, 10, 0, 3)
r2d_keff_inel_rw.SetTitle("Inelastic-scattering Protons;Proton Kinetic Energy at TPC FF [MeV]; Data/MC")
r2d_keff_inel_rw.GetXaxis().CenterTitle()
r2d_keff_inel_rw.Draw("")
R_keffbeam_inel.Draw("ep same")

line_keff_inel_rw=RT.TLine(220, 1, 600, 1)
line_keff_inel_rw.SetLineColor(1)
line_keff_inel_rw.SetLineStyle(2)
line_keff_inel_rw.Draw("same")


r1_inel_rw.Print(args.o+'/ratio_ff_inel_rw.eps')










#KEnd(stop)
c4_rw=RT.TCanvas("c4","",1200,900)
c4_rw.Divide(1,1)
c4_rw.cd(1)

f2d_kend_calo_rw=RT.TH2D("f2d_kend_calo_rw","", 450, -150, 300, 600, 0, kend_calo_stop_data_rw.GetBinContent(kend_calo_stop_data_rw.GetMaximumBin())+800)
f2d_kend_calo_rw.SetTitle("Stopping Protons;Proton Kinetic Energy at Track End [MeV];")
f2d_kend_calo_rw.GetXaxis().CenterTitle()
f2d_kend_calo_rw.Draw()

kend_calo_stop_data_rw.SetLineColor(8)
kend_calo_stop_data_rw.SetMarkerColor(8)
fit_kend_calo_stop_data.SetLineColor(8)

kend_calo_stop_mc_rw.SetLineColor(7)
kend_calo_stop_mc_rw.SetMarkerColor(7)
fit_kend_calo_stop_mc_rw.SetLineColor(7)

fit_kend_calo_stop_data.Draw("same")
fit_kend_calo_stop_mc.Draw("same")
fit_kend_true_stop_mc.Draw("same")
kend_true_stop_mc.Draw("hist same")
kend_calo_stop_mc.Draw("hist same")
kend_calo_stop_data.Draw("ep same")

fit_kend_calo_stop_data_rw.Draw("same")
fit_kend_calo_stop_mc_rw.Draw("same")
kend_calo_stop_mc_rw.Draw("hist same")
kend_calo_stop_data_rw.Draw("ep same")



leg4_rw=RT.TLegend(0.14,0.65,.9,0.9)
leg4_rw.SetFillStyle(0)
txt4_rw=[]
txt4_rw.append("Data(calo): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_data[5],er_mu_stop_data[5],sigma_stop_data[5],er_sigma_stop_data[5]))
txt4_rw.append("MC(calo): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[7],er_mu_stop_mc[7],sigma_stop_mc[7],er_sigma_stop_mc[7]))
txt4_rw.append("MC(truth): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[9],er_mu_stop_mc[9],sigma_stop_mc[9],er_sigma_stop_mc[9]))

print(txt4_rw[0])
leg4_rw.AddEntry(kend_calo_stop_data, txt4_rw[0], "ep")
leg4_rw.AddEntry(kend_calo_stop_mc, txt4_rw[1], "l")
leg4_rw.AddEntry(kend_true_stop_mc, txt4_rw[2], "l")

leg4_rw.Draw()
c4_rw.Print(args.o+'/kend_calo_stop_rw.eps')





#el -----------------------------------------------------------------------------------------------------------------------------------------
c4_el_rw=RT.TCanvas("c4_el_rw","",1200,900)
c4_el_rw.Divide(1,1)
c4_el_rw.cd(1)

f2d_kend_calo_el_rw=RT.TH2D("f2d_kend_calo_rw","", 450, -150, 300, 600, 0, kend_calo_el_data_rw.GetBinContent(kend_calo_el_data_rw.GetMaximumBin())+800)
f2d_kend_calo_el_rw.SetTitle("Elastic-scattering Protons;Proton Kinetic Energy at Track End [MeV];")
f2d_kend_calo_el_rw.GetXaxis().CenterTitle()
f2d_kend_calo_el_rw.Draw()


kend_calo_el_data_rw.SetLineColor(1)
kend_calo_el_data_rw.SetMarkerColor(1)

kend_calo_el_data_rw.SetLineWidth(1)
fit_kend_calo_el_data_rw.SetLineColor(1)

kend_calo_el_mc_rw.SetLineColor(1)
kend_calo_el_mc_rw.SetMarkerColor(1)

kend_calo_el_mc_rw.SetLineWidth(1)
fit_kend_calo_el_mc_rw.SetLineColor(1)

kend_calo_el_inel_mc_rw.SetFillColor(2)
kend_calo_el_el_mc_rw.SetFillColor(4)
kend_calo_el_midcosmic_mc_rw.SetFillColor(5)
kend_calo_el_midpi_mc_rw.SetFillColor(6)
kend_calo_el_midp_mc_rw.SetFillColor(3)
kend_calo_el_midmu_mc_rw.SetFillColor(28)
kend_calo_el_mideg_mc_rw.SetFillColor(30)
kend_calo_el_midother_mc_rw.SetFillColor(15)

kend_calo_el_inel_mc_rw.SetLineColor(2)
kend_calo_el_el_mc_rw.SetLineColor(4)
kend_calo_el_midcosmic_mc_rw.SetLineColor(5)
kend_calo_el_midpi_mc_rw.SetLineColor(6)
kend_calo_el_midp_mc_rw.SetLineColor(3)
kend_calo_el_midmu_mc_rw.SetLineColor(28)
kend_calo_el_mideg_mc_rw.SetLineColor(30)
kend_calo_el_midother_mc_rw.SetLineColor(15)

hs_kend_calo_el_rw=RT.THStack("hs_kend_calo_el_rw","");
hs_kend_calo_el_rw.Add(kend_calo_el_inel_mc_rw)
hs_kend_calo_el_rw.Add(kend_calo_el_el_mc_rw)
hs_kend_calo_el_rw.Add(kend_calo_el_midp_mc_rw)
hs_kend_calo_el_rw.Add(kend_calo_el_midcosmic_mc_rw)
hs_kend_calo_el_rw.Add(kend_calo_el_midpi_mc_rw)
hs_kend_calo_el_rw.Add(kend_calo_el_midmu_mc_rw)
hs_kend_calo_el_rw.Add(kend_calo_el_mideg_mc_rw)
hs_kend_calo_el_rw.Add(kend_calo_el_midother_mc_rw)

kend_calo_el_mc_rw.SetLineWidth(1)
kend_calo_el_data_rw.SetLineWidth(1)

hs_kend_calo_el_rw.Draw("hist same")
#fit_kend_calo_stop_mc_rw.Draw("same")
kend_calo_el_mc_rw.Draw("hist same")
#fit_kend_calo_stop_data_rw.Draw("same")
kend_calo_el_data_rw.Draw("ep same")


tleg111_rw=RT.TLegend(0.14,0.75,.9,0.87)
tleg111_rw.SetFillStyle(0)
txt111_rw=[]
txt111_rw.append("Data(rw): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_el_data_rw[1],er_mu_el_data_rw[1],sigma_el_data_rw[1],er_sigma_el_data_rw[1]))
txt111_rw.append("MC(rw): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_el_mc_rw[1],er_mu_el_mc_rw[1],sigma_el_mc_rw[1],er_sigma_el_mc_rw[1]))
tleg111_rw.AddEntry(kend_calo_el_data_rw, txt111_rw[0], "ep")
tleg111_rw.AddEntry(0, "", "")
tleg111_rw.AddEntry(kend_calo_el_mc_rw, txt111_rw[1], "l")
tleg111_rw.AddEntry(0, "", "")
tleg111_rw.SetNColumns(2);


leg111_rw=RT.TLegend(0.14,0.62,.88,0.73)
leg111_rw.SetFillStyle(0)
leg111_rw.AddEntry(kend_calo_el_inel_mc_rw, "Inel", "f")
leg111_rw.AddEntry(kend_calo_el_el_mc_rw, "El", "f")
leg111_rw.AddEntry(kend_calo_el_midcosmic_mc_rw, "misID:cosmic", "f")
leg111_rw.AddEntry(kend_calo_el_midp_mc_rw, "misID:p", "f")
leg111_rw.AddEntry(kend_calo_el_midpi_mc_rw, "misID:#pi", "f")
leg111_rw.AddEntry(kend_calo_el_midmu_mc_rw, "misID:#mu", "f")
leg111_rw.AddEntry(kend_calo_el_mideg_mc_rw, "misID:e/#gamma", "f")
leg111_rw.AddEntry(kend_calo_el_midother_mc_rw, "misID:other", "f")
leg111_rw.SetNColumns(4);

tleg111_rw.Draw()
leg111_rw.Draw()

c4_el_rw.Print(args.o+'/kend_calo_el_rw.eps')


#KEffbeam_el: data/mc ratio -----------------------------------------------------------------------------
R_kend_calo_el_rw=kend_calo_el_data_rw.Clone()
R_kend_calo_el_rw.Divide(R_kend_calo_el_rw, kend_calo_el_mc_rw)

re1_el_rw=RT.TCanvas("re1_el_rw","",1200,900)
re1_el_rw.Divide(1,1)
re1_el_rw.cd(1)

r2d_kend_el_rw=RT.TH2D("r2d_kend_el_rw","", 450, -150, 300, 10, 0, 3)
r2d_kend_el_rw.SetTitle("Elastic-scattering Protons;Proton Kinetic Energy at Track End [MeV]; Data/MC")
r2d_kend_el_rw.GetXaxis().CenterTitle()
r2d_kend_el_rw.Draw("")
R_kend_calo_el_rw.Draw("ep same")

line_kend_el_rw=RT.TLine(-150, 1, 300, 1)
line_kend_el_rw.SetLineColor(1)
line_kend_el_rw.SetLineStyle(2)
line_kend_el_rw.Draw("same")


re1_el_rw.Print(args.o+'/ratio_end_el_rw.eps')




'''
#Construct the weighting func ----------------------------------------
#denominator: KEfit_stop_data
#Nominator: [1] data: (KEbeam-constE)
#           [2] MC: (KEbeam-constE)

#Weighting func (data)
mu_data_denom=mu_stop_data[2]
sigma_data_denom=sigma_stop_data[2]

mu_data_nom=mu_stop_data[1]
sigma_data_nom=sigma_stop_data[1]

#Weighting func (MC)
mu_mc_nom=mu_stop_mc[3]
sigma_mc_nom=sigma_stop_mc[3]

#Gaussian
print("(data) mu_data_denom:",mu_data_denom)
print("(data) sigma_data_denom:",sigma_data_denom)

print("(data) mu_data_nom:",mu_data_nom)
print("(data) sigma_data_nom:",sigma_data_nom)

print("(mc) mu_nom:",mu_mc_nom)
print("(mc) sigma_nom:",sigma_mc_nom)


#FF
#test=[]
#ex=[0,0,0]

#test1 = [1, 2, 3]
#test.append([4, 5, 6])

#test = array.array('d', [1, 2, 3])
#test.append([1, 2, 3]) 
#test.append([4, 5, 6]) 
#print(test)
#print(len(test))
#print(len(test[0]))

gr_stop_mc = RT.TGraphErrors(len(mu_stop_mc), array('d', mu_stop_mc[:]), array('d', sigma_stop_mc[:]), array('d', er_mu_stop_mc[:]), array('d', er_sigma_stop_mc[:]))

cx_stop_mc=RT.TCanvas("cx_stop_mc","",1200,900)
cx_stop_mc.Divide(1,1)
cx_stop_mc.cd(1)
gr_stop_mc.Draw("ap")
cx_stop_mc.Print('test.eps')
'''

