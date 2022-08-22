import ROOT as RT
import sys
import math

from argparse import ArgumentParser as ap
import numpy as np
from array import array
from math import sqrt
from math import exp


def fitg(x, p):
  m=p[0]
  s=p[1]
  n=p[2]
  g=n*math.exp(-(x[0]-m)*(x[0]-m)/(2*s*s))
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

  h.Fit("g","remn")
  return g


#MC File ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
file_mc='../mc_proton_beamxy_beammom_nobmrw.root'
out_path='./keend_study'

#Plt Style -------------------------------------------------------------------------------
parser = ap()

RT.gROOT.LoadMacro("~/protoDUNEStyle.C")
RT.gROOT.SetBatch(); #PyROOT does NOT display any graphics in batch mode
RT.gStyle.SetOptStat(00000)
#RT.gStyle.SetErrorX(1.e-4)
RT.gStyle.SetTitleAlign(23)
RT.gStyle.SetTitleX(.5)
RT.gStyle.SetLineWidth(1)
tt = RT.TLatex();
tt.SetNDC();

#read files --------------------
f_mc=RT.TFile(file_mc, "OPEN")

#end ----------------------------------------------------------------------
#calo
kend_calo_el=f_mc.Get("h1d_kend_calo_el")
kend_calo_inel=f_mc.Get("h1d_kend_calo_inel")

#bb
kend_bb_el=f_mc.Get("h1d_kend_bb_el")
kend_bb_inel=f_mc.Get("h1d_kend_bb_inel")

#truth
kend_true_el=f_mc.Get("h1d_kend_true_el")
kend_true_inel=f_mc.Get("h1d_kend_true_inel")


#[1]------------------------------------------------------------------------------------------------------
c0_end_mc=RT.TCanvas("c0_end_mc","",1200,900)
c0_end_mc.Divide(1,1)
c0_end_mc.cd(1).SetLogy(0)

#f2d_mc=RT.TH2D("f2d_mc","", 700, -100, 600, 900, 0.01, 2000) #logy
f2d_mc=RT.TH2D("f2d_mc","", 700, -100, 600, 900, 0., 1000) #liny
f2d_mc.SetTitle("Reconstructed Elastic-scattering Protons; Kinetic Energy at Track End [MeV]; Events")
f2d_mc.Draw("")

kend_calo_el.SetLineColor(2)
kend_calo_el.SetMarkerColor(2)
kend_bb_el.SetLineColor(4)
kend_bb_el.SetMarkerColor(4)
kend_true_el.SetLineColor(3)
kend_true_el.SetMarkerColor(3)

kend_true_el.Draw("hist same")
kend_calo_el.Draw("hist same")
kend_bb_el.Draw("hist same")

leg0=RT.TLegend(0.7,0.7,.87,0.9)
leg0.SetFillStyle(0)
leg0.AddEntry(kend_calo_el, "Calo.", "l")
leg0.AddEntry(kend_bb_el, "BB.", "l")
leg0.AddEntry(kend_true_el, "Truth", "l")
#leg0.SetNColumns(2);
leg0.Draw()


#KEffbeam_KEhy_stop.Draw("colz same")
#ll=RT.TLine(200,200,600,600)
#ll.SetLineColor(2)
#ll.SetLineStyle(2)
#ll.Draw()
c0_end_mc.Print(out_path+'/keend_el.eps')


#[2]------------------------------------------------------------------------------------------------------
c1_end_mc=RT.TCanvas("c1_end_mc","",1200,900)
c1_end_mc.Divide(1,1)
c1_end_mc.cd(1).SetLogy(0)

f2d_mc=RT.TH2D("f2d_mc","", 700, -100, 600, 900, 0., 400) #liny
f2d_mc.SetTitle("Reconstructed Inelastic-scattering Protons; Kinetic Energy at Track End [MeV]; Events")
f2d_mc.Draw("")

kend_calo_inel.SetLineColor(2)
kend_calo_inel.SetMarkerColor(2)
kend_bb_inel.SetLineColor(4)
kend_bb_inel.SetMarkerColor(4)
kend_true_inel.SetLineColor(3)
kend_true_inel.SetMarkerColor(3)

kend_true_inel.Draw("hist same")
kend_calo_inel.Draw("hist same")
kend_bb_inel.Draw("hist same")

leg1=RT.TLegend(0.7,0.7,.87,0.9)
leg1.SetFillStyle(0)
leg1.AddEntry(kend_calo_inel, "Calo.", "l")
leg1.AddEntry(kend_bb_inel, "BB.", "l")
leg1.AddEntry(kend_true_inel, "Truth", "l")
#leg0.SetNColumns(2);
leg1.Draw()


#KEffbeam_KEhy_stop.Draw("colz same")
#ll=RT.TLine(200,200,600,600)
#ll.SetLineColor(2)
#ll.SetLineStyle(2)
#ll.Draw()
c1_end_mc.Print(out_path+'/keend_inel.eps')

