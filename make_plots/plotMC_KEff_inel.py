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

#read data --------------------------------------------------------------------------------------------------------------------------------------------------
file_mc="../mc_proton_beamxy_beammom_nobmrw_studyKEconst.root"
#file_mc="../mc_proton_beamxy_beammom_nobmrw_studyKEfit.root"
#file_mc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/PDSPProd4a_MC_6GeV_gen_datadriven_reco1_sce_off_v1/mc_sceoff_KEFit.root"

f_mc=RT.TFile(file_mc, "OPEN")
trklen_dkeff_inel=f_mc.Get("h2d_trklen_dkeff_inel")
trklen_dkeff_el=f_mc.Get("h2d_trklen_dkeff_el")
out_path="./plots_KEffFit_inel/"

#Plots #########################################################################################################################################################
#KEbeam ------------------------------------------------------------------------------------------------------------------
c0=RT.TCanvas("c0","",1200,900)
c0.Divide(1,1)
c0.cd(1)

f2d_0=RT.TH2D("f2d_0","", 1400, 0, 140, 1000, -700, 300)
f2d_0.SetTitle("Inelastic-scattering Candidates; Reconstructed Track Length [cm]; KE_{FF}(Fit)-KE_{FF}(truth) [MeV]")
f2d_0.GetXaxis().CenterTitle()
f2d_0.GetYaxis().CenterTitle()
f2d_0.Draw()
trklen_dkeff_inel.Draw("colz same")

ll=RT.TLine(0,0,140,0)
ll.SetLineColor(2)
ll.SetLineStyle(2)
ll.SetLineWidth(2)
ll.Draw("same")
c0.Print(out_path+'/trklen_kefit_minus_keffconst_inel.eps')
#c0.Print(out_path+'/trklen_kefit_minus_keff_inel.eps')
#c0.Print(out_path+'/trklen_kefit_minus_keff_inel_sceoff.eps')

#el
c1=RT.TCanvas("c0","",1200,900)
c1.Divide(1,1)
c1.cd(1)
f2d_1=RT.TH2D("f2d_1","", 1400, 0, 140, 1000, -700, 300)
f2d_1.SetTitle("Elastic-scattering Candidates; Reconstructed Track Length [cm]; KE_{FF}(Fit)-KE_{FF}(truth) [MeV]")
f2d_1.GetXaxis().CenterTitle()
f2d_1.GetYaxis().CenterTitle()
f2d_1.Draw()
trklen_dkeff_el.Draw("colz same")
ll.Draw("same")
#c1.Print(out_path+'/trklen_kefit_minus_keff_el.eps')
c1.Print(out_path+'/trklen_kefit_minus_keffconst_el.eps')

f2d_2=RT.TH2D("f2d_1","", 1400, 0, 140, 2000, -100, 100)
f2d_2.SetTitle("Elastic-scattering Candidates; Reconstructed Track Length [cm]; KE_{FF}(Fit)-KE_{FF}(truth) [MeV]")
f2d_2.GetXaxis().CenterTitle()
f2d_2.GetYaxis().CenterTitle()
f2d_2.Draw()
trklen_dkeff_el.Draw("colz same")
ll.Draw("same")
#c1.Print(out_path+'/trklen_kefit_minus_keff_el_zoom.eps')
c1.Print(out_path+'/trklen_kefit_minus_keffconst_el_zoom.eps')


