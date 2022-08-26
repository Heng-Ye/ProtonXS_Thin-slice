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
  pre_mean=h.GetMean()
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


def VRFit(h, fmin, fmax):
  pre_mean=h.GetMean()
  pre_max=h.GetMaximum()
  pre_rms=h.GetRMS()
  n_sigma=1
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
  #g=RT.TF1("g", fitg, gg.GetParameter(0)-n_sigma*gg.GetParameter(1), gg.GetParameter(0)+n_sigma*gg.GetParameter(1), 3)
  g=RT.TF1("g", fitg, fmin, fmax, 3)
  g.SetParameter(0, gg.GetParameter(0))
  g.SetParameter(1, gg.GetParameter(1))
  g.SetParameter(2, gg.GetParameter(2))

  g.SetLineStyle(2)
  g.SetLineWidth(2)

  h.Fit("g","remn")
  return g


#MC File ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#file_mc='../mc_proton_beamxy_beammom_nobmrw.root'
file_mc='../mc_proton_beamxy_beammom_nobmrw_KEffeqFit.root'
out_path='./keend_summary'

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

#set colors
kend_calo_inel.SetLineColor(2)
kend_calo_inel.SetMarkerColor(2)

kend_bb_inel.SetLineColor(2)
kend_bb_inel.SetMarkerColor(2)

kend_calo_el.SetLineColor(4)
kend_calo_el.SetMarkerColor(4)

kend_bb_el.SetLineColor(4)
kend_bb_el.SetMarkerColor(4)

kend_true_inel.SetLineColor(3)
kend_true_inel.SetMarkerColor(3)

kend_true_el.SetLineColor(7)
kend_true_el.SetMarkerColor(7)

kend_calo_inel.SetLineWidth(1)
kend_calo_el.SetLineWidth(1)
kend_bb_inel.SetLineWidth(1)
kend_bb_el.SetLineWidth(1)
kend_true_inel.SetLineWidth(1)
kend_true_el.SetLineWidth(1)


#Summary plot
c0_end_mc=RT.TCanvas("c0_end_mc","",1200,900)
c0_end_mc.Divide(1,1)
c0_end_mc.cd(1).SetLogy(0)

f2d_mc=RT.TH2D("f2d_mc","", 700, -100, 600, 900, 0, 900)
#f2d_mc.SetTitle("KE (Bethe-Bloch) with KE_{FF}=(KE_{beam}-#DeltaE)*R; Kinetic Energy at Track End [MeV]; Events")
f2d_mc.SetTitle("KE (Bethe-Bloch) with KE_{FF}=KE(Fit); Kinetic Energy at Track End [MeV]; Events")
f2d_mc.Draw("")

kend_true_inel.Draw("hist same")
kend_true_el.Draw("hist same")
kend_bb_inel.Draw("hist same")
kend_bb_el.Draw("hist same")

leg=RT.TLegend(0.5,0.6,.91,0.87)
leg.SetFillStyle(0)
leg.AddEntry(kend_bb_inel, "Inel.(reco.)", "l")
leg.AddEntry(kend_true_inel, "Inel.(truth)", "l")
leg.AddEntry(kend_bb_el, "El.(reco.)", "l")
leg.AddEntry(kend_true_el, "El.(truth)", "l")
leg.Draw()
c0_end_mc.Print(out_path+'/kend_bb_summary_KEffIsKEFit.eps')


#f2d_mc.SetTitle("KE (Calo) with KE_{FF}=(KE_{beam}-#DeltaE)*R; Kinetic Energy at Track End [MeV]; Events")
f2d_mc.SetTitle("KE (Calo) with KE_{FF}=KE(Fit); Kinetic Energy at Track End [MeV]; Events")
f2d_mc.Draw("")
kend_true_inel.Draw("hist same")
kend_true_el.Draw("hist same")
kend_calo_el.Draw("hist same")
kend_calo_inel.Draw("hist same")
leg.Draw()
#c0_end_mc.Print(out_path+'/kend_calo_summary.eps')
c0_end_mc.Print(out_path+'/kend_calo_summary_KEffIsKEFit.eps')





