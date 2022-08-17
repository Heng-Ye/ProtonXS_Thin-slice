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


file_mc='/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/mc_kecalo_nobmrw_beamxy_ElossTuneHyper.root'
file_data='/dune/data2/users/hyliao/protonana/v09_39_01/KEHY_ELoss/data_kehy_beamxy_runAll.root'
out_path='./keff_study'

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


#ff --------------------------------------------------
keff_truth_mc=f_mc.Get("h1d_keff_truth")
keff_truth_stop_mc=f_mc.Get("h1d_keff_truth_stop")
keff_truth_el_mc=f_mc.Get("h1d_keff_truth_el")
keff_truth_inel_mc=f_mc.Get("h1d_keff_truth_inel")

n_keff_truth_mc=keff_truth_mc.Integral()
n_keff_truth_stop_mc=keff_truth_stop_mc.Integral()

#hy --------------------------------------------------
kehy_mc=f_mc.Get("h1d_kehy")
kehy_stop_mc=f_mc.Get("h1d_kehy_stop")
n_kehy_mc=kehy_mc.Integral()
n_kehy_stop_mc=kehy_stop_mc.Integral()


#normalization ------------------------------------------------
keff_truth_stop_mc.Scale(n_keff_truth_mc/n_keff_truth_stop_mc)
#kehy_mc.Scale(n_keff_truth_mc/n_kehy_mc)
kehy_stop_mc.Scale(n_keff_truth_mc/n_kehy_stop_mc)


#fit------------------------------------------------------------
fit_keff_truth_mc=VNFit(keff_truth_mc, 430, 3)
fit_keff_truth_stop_mc=VNFit(keff_truth_stop_mc, 430, 3)
fit_kehy_stop_mc=VNFit(kehy_stop_mc, 430, 3)
fit_kehy_mc=VNFit(kehy_mc, 430, 3)


#save fits --------------------------------
n_mc_stop=4

mu_stop_mc=[]
er_mu_stop_mc=[]
sigma_stop_mc=[]
er_sigma_stop_mc=[]

for i in range(n_mc_stop):
  fit=fit_keff_truth_mc
  if i==0: 
    fit=fit_keff_truth_mc
  if i==1: 
    fit=fit_keff_truth_stop_mc
  if i==2: 
    fit=fit_kehy_mc
  if i==3: 
    fit=fit_kehy_stop_mc

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




#KEff(MC): Comparsion of cut effect ------------------------------------------------------------------------------------------------------------------
c0_ff_mc=RT.TCanvas("c0_ff_mc","",1200,900)
c0_ff_mc.Divide(1,1)
c0_ff_mc.cd(1)

f2d_ff_mc=RT.TH2D("f2d_ff_mc","", 400, 220, 620, 600, 0, keff_truth_mc.GetBinContent(keff_truth_mc.GetMaximumBin())+5000)
f2d_ff_mc.SetTitle(";Proton Kinetic Energy at TPC FF[MeV];")
f2d_ff_mc.GetXaxis().CenterTitle()
f2d_ff_mc.Draw()

keff_truth_mc.SetLineColor(1)
keff_truth_stop_mc.SetLineColor(2)
keff_truth_stop_mc.SetMarkerColor(2)
kehy_mc.SetLineColor(3)
kehy_mc.SetMarkerColor(3)
kehy_stop_mc.SetLineColor(4)
kehy_stop_mc.SetMarkerColor(4)

fit_keff_truth_mc.SetLineColor(1)
fit_keff_truth_mc.SetLineStyle(2)
fit_keff_truth_stop_mc.SetLineColor(2)
fit_keff_truth_stop_mc.SetLineStyle(2)
fit_kehy_mc.SetLineColor(3)
fit_kehy_mc.SetLineStyle(2)
fit_kehy_stop_mc.SetLineColor(4)
fit_kehy_stop_mc.SetLineStyle(2)






keff_truth_mc.Draw("hist same")
fit_keff_truth_mc.Draw("same")
keff_truth_stop_mc.Draw("hist same")
fit_keff_truth_stop_mc.Draw("same")
kehy_mc.Draw("hist same")
fit_kehy_mc.Draw("same")
kehy_stop_mc.Draw("hist same")
fit_kehy_stop_mc.Draw("same")

leg_mc=RT.TLegend(0.14,0.55,.9,0.9)
leg_mc.SetFillStyle(0)
txt_mc=[]
txt_mc.append("MC Truth (before stopping proton cut): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[0],er_mu_stop_mc[0],sigma_stop_mc[0],er_sigma_stop_mc[0]))
txt_mc.append("MC Truth (stopping proton cut): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[1],er_mu_stop_mc[1],sigma_stop_mc[1],er_sigma_stop_mc[1]))
txt_mc.append("MC Fit (before proton cut): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[2],er_mu_stop_mc[2],sigma_stop_mc[2],er_sigma_stop_mc[2]))
txt_mc.append("MC Fit (stopping proton cut): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_stop_mc[3],er_mu_stop_mc[3],sigma_stop_mc[3],er_sigma_stop_mc[3]))

leg_mc.AddEntry(keff_truth_mc, txt_mc[0], "l")
leg_mc.AddEntry(keff_truth_stop_mc, txt_mc[1], "l")
leg_mc.AddEntry(kehy_mc, txt_mc[2], "l")
leg_mc.AddEntry(kehy_stop_mc, txt_mc[3], "l")

leg_mc.Draw()


c0_ff_mc.Print(out_path+'/keff_before_after_stop_cut.eps')


