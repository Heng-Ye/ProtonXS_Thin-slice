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
file_mc='../mc_keff.root'
file_mc_2d='../mc_proton_beamxy_beammom_nobmrw.root'
out_path='./keff_study'
#file_mc='../mc_keff_newlikelihood.root'
#out_path='./keff_study_likelihood'

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

#read files -------------------------
f_mc=RT.TFile(file_mc, "OPEN")
f_mc_2d=RT.TFile(file_mc_2d, "OPEN")

#ff ------------------------------------------------------------------------
KEffbeam_KEhy_stop=f_mc.Get("h2d_KEffbeam_KEhy_stop")
KEffbeam_KEhy_inel=f_mc.Get("h2d_KEffbeam_KEhy_inel")

KEffbeam_dKEhy_stop=f_mc.Get("h2d_KEffbeam_dKEhy_stop")
KEffbeam_dKEhy_inel=f_mc.Get("h2d_KEffbeam_dKEhy_inel")

ratio_KEffbeam_KEhy_stop=f_mc.Get("h1d_ratio_KEffbeam_KEhy_stop")

keffbeam_keff_stop=f_mc_2d.Get("h2d_keffbeam_keff_stop")
keffbeam_keff_el=f_mc_2d.Get("h2d_keffbeam_keff_el")
keffbeam_keff_inel=f_mc_2d.Get("h2d_keffbeam_keff_inel")

kehy_keff_stop=f_mc_2d.Get("h2d_kehy_keff_stop")
kehy_keff_el=f_mc_2d.Get("h2d_kehy_keff_el")
kehy_keff_inel=f_mc_2d.Get("h2d_kehy_keff_inel")




#[1]
c0_ff_mc=RT.TCanvas("c0_ff_mc","",1200,900)
c0_ff_mc.Divide(1,1)
c0_ff_mc.cd(1)

f2d_mc=RT.TH2D("f2d_mc","", 400,200,600, 400,200,600)
f2d_mc.SetTitle("Stopping Protons; KE_{beam}-#DeltaE [MeV]; KE(fit) [MeV]")
f2d_mc.Draw("")

KEffbeam_KEhy_stop.Draw("colz same")
ll=RT.TLine(200,200,600,600)
ll.SetLineColor(2)
ll.SetLineStyle(2)
ll.Draw()
c0_ff_mc.Print(out_path+'/keffbeam_kehy_stop.eps')

#[2]
f2d_mc2=RT.TH2D("f2d_mc","", 400,200,600, 500,100,600)
f2d_mc2.SetTitle("Inelastic-scattering Protons; KE_{beam}-#DeltaE [MeV]; KE(fit) [MeV]")
f2d_mc2.Draw("")

KEffbeam_KEhy_inel.Draw("colz same")
ll2=RT.TLine(200,200,600,600)
ll2.SetLineColor(2)
ll2.SetLineStyle(2)
ll2.Draw()
c0_ff_mc.Print(out_path+'/keffbeam_kehy_inel.eps')

#[3]
f2dx_mc=RT.TH2D("f2d_mc","", 400,200,600, 1200,-600,600)
f2dx_mc.SetTitle("Stopping Protons; KE_{beam}-#DeltaE [MeV]; KE(fit)-KE(truth) [MeV]")
f2dx_mc.Draw()
KEffbeam_dKEhy_stop.Draw("colz same")

ll3=RT.TLine(200,0,600,0)
ll3.SetLineColor(2)
ll3.SetLineStyle(2)
ll3.Draw()
c0_ff_mc.Print(out_path+'/keffbeam_dkehy_stop.eps')

#[4]
#f2dx_mc=RT.TH2D("f2d_mc","", 400,200,600, 1200,-600,600)
f2dx_mc.SetTitle("Inelastic-scattering Protons; KE_{beam}-#DeltaE [MeV]; KE(fit)-KE(truth) [MeV]")
f2dx_mc.Draw()
KEffbeam_dKEhy_inel.Draw("colz same")

ll4=RT.TLine(200,0,600,0)
ll4.SetLineColor(2)
ll4.SetLineStyle(2)
ll4.Draw()
c0_ff_mc.Print(out_path+'/keffbeam_dkehy_inel.eps')

#[5]Ratio of (KEHY/(KEbeam-dE))
c0_r_mc=RT.TCanvas("c0_r_mc","",1200,900)
c0_r_mc.Divide(1,1)
c0_r_mc.cd(1)
f2d_ratio_KEffbeam_KEhy_stop=RT.TH2D("f2d_ratio_KEffbeam_KEhy_stop","", 10,0.5,1.5, 100,0,3000)
f2d_ratio_KEffbeam_KEhy_stop.SetTitle("Stopping Protons; KE(fit)/(KE_{beam}-#DeltaE); Events")
f2d_ratio_KEffbeam_KEhy_stop.Draw()
ratio_KEffbeam_KEhy_stop.Draw("same")

fit_ratio_KEffbeam_KEhy_stop=VNFit(ratio_KEffbeam_KEhy_stop, 1, 3)
fit_ratio_KEffbeam_KEhy_stop.SetLineColor(2)
fit_ratio_KEffbeam_KEhy_stop.SetMarkerStyle(2)
fit_ratio_KEffbeam_KEhy_stop.Draw("same")

mu=fit_ratio_KEffbeam_KEhy_stop.GetParameter(0)
er_mu=fit_ratio_KEffbeam_KEhy_stop.GetParError(0)
sigma=fit_ratio_KEffbeam_KEhy_stop.GetParameter(1)
er_sigma=fit_ratio_KEffbeam_KEhy_stop.GetParError(1)

print('mu:',mu)
print('sigma:',sigma)

leg0=RT.TLegend(0.1,0.7,.92,0.95)
leg0.SetFillStyle(0)
txt0=[]
txt0.append("Data: #mu={:.4f}#pm{:.4f} MeV, #sigma={:.4f}#pm{:.4f} MeV".format(mu,er_mu,sigma,er_sigma))
leg0.AddEntry(ratio_KEffbeam_KEhy_stop, txt0[0], "l")
#leg0.SetNColumns(2);
leg0.Draw()

c0_r_mc.Print(out_path+'/ratio_kehy_keffbeam_stop.eps')


#[6]KEff(truth) vs (KEbeam-dE)*R --------------------------------------------------
c1_ff_mc=RT.TCanvas("c1_ff_mc","",1200,900)
c1_ff_mc.Divide(1,1)
c1_ff_mc.cd(1)

f2d1_mc=RT.TH2D("f2d1_mc","", 400, 200, 600, 400, 200, 600)
f2d1_mc.SetTitle("Stopping Protons; (KE_{beam}-#DeltaE)*R [MeV]; KE(truth) [MeV]")
f2d1_mc.Draw("")

keffbeam_keff_stop.Draw("colz same")
ll1=RT.TLine(200,200,600,600)
ll1.SetLineColor(2)
ll1.SetLineStyle(2)
ll1.Draw()
c1_ff_mc.Print(out_path+'/keffbeam_keff_stop.eps')

f2d1_mc.SetTitle("Elastic-scattering Protons; (KE_{beam}-#DeltaE)*R [MeV]; KE(truth) [MeV]")
f2d1_mc.Draw("")
keffbeam_keff_el.Draw("colz same")
ll1.Draw()
c1_ff_mc.Print(out_path+'/keffbeam_keff_el.eps')

profx_keffbeam_keff_el=keffbeam_keff_el.ProfileX()
f2d1_mc.Draw("")
keffbeam_keff_el.Draw("colz same")
fit_profx_keffbeam_keff_el=profx_keffbeam_keff_el.Fit("pol1","","same",320,510)

ll2=RT.TLine(320,320,510,510)
ll2.SetLineColor(2)
ll2.SetLineStyle(2)
ll2.Draw()

c1_ff_mc.Print(out_path+'/keffbeam_keff_el_with_fit.eps')








f2d1_mc.SetTitle("Inelastic-scattering Protons; (KE_{beam}-#DeltaE)*R [MeV]; KE(truth) [MeV]")
f2d1_mc.Draw("")
keffbeam_keff_inel.Draw("colz same")
ll1.Draw()
c1_ff_mc.Print(out_path+'/keffbeam_keff_inel.eps')


profx_keffbeam_keff_inel=keffbeam_keff_inel.ProfileX()
f2d1_mc.Draw("")
keffbeam_keff_inel.Draw("colz same")
fit_profx_keffbeam_keff_inel=profx_keffbeam_keff_inel.Fit("pol1","","same",320,510)

ll2.Draw()


c1_ff_mc.Print(out_path+'/keffbeam_keff_inel_with_fit.eps')






