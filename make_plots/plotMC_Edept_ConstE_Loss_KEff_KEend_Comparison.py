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


#MC File ---------------------------------------------------------------------
file_mc='../mc_proton_beamxy_beammom_nobmrw_ElasFitted_v09_39_01.root'
out_path='./keff_keend_comparison/ElasFitted'
#file_mc='../mc_proton_beamxy_beammom_nobmrw_ElasTrue_v09_39_01.root'
#out_path='./keff_keend_comparison/ElasTrue'
#file_mc='../mc_proton_beamxy_beammom_nobmrw_AllFitted_v09_39_01.root'
#out_path='./keff_keend_comparison/AllFitted'
#file_mc='../mc_proton_beamxy_beammom_nobmrw_AllTrue_v09_39_01.root'
#out_path='./keff_keend_comparison/AllTrue'

#read files -------------------------------------
f_mc=RT.TFile(file_mc, "OPEN")

#ff_truth
keff_el=f_mc.Get("h1d_keff_el")
keff_inel=f_mc.Get("h1d_keff_inel")

#ff_reco
keffbeam_el=f_mc.Get("h1d_keffbeam_el")
keffbeam_inel=f_mc.Get("h1d_keffbeam_inel")

#end_bb
kend_bb_el=f_mc.Get("h1d_kend_bb_el")
kend_bb_inel=f_mc.Get("h1d_kend_bb_inel")

#end_true
kend_true_inel=f_mc.Get("h1d_kend_true_inel")
kend_true_el=f_mc.Get("h1d_kend_true_el")

#calo
kend_calo_el=f_mc.Get("h1d_kend_calo_el")
kend_calo_inel=f_mc.Get("h1d_kend_calo_inel")

#fit
pre_mean_ff=400
n_s=3
fit_keff_el=VNFit(keff_el, pre_mean_ff, n_s)
fit_keff_inel=VNFit(keff_inel, pre_mean_ff, n_s)
fit_keffbeam_el=VNFit(keffbeam_el, pre_mean_ff, n_s)
fit_keffbeam_inel=VNFit(keffbeam_inel, pre_mean_ff, n_s)


#MC File (const-E) ------------------------------------------------------------
file_mc2='../mc_proton_beamxy_beammom_nobmrw_constE_v09_39_01.root'

#read files -------------------------------------
f_mc2=RT.TFile(file_mc2, "OPEN")

#ff_truth
keff_el2=f_mc2.Get("h1d_keff_el")
keff_inel2=f_mc2.Get("h1d_keff_inel")

#ff_reco
keffbeam_el2=f_mc2.Get("h1d_keffbeam_el")
keffbeam_inel2=f_mc2.Get("h1d_keffbeam_inel")

#end_bb
kend_bb_el2=f_mc2.Get("h1d_kend_bb_el")
kend_bb_inel2=f_mc2.Get("h1d_kend_bb_inel")

#end_true
kend_true_inel2=f_mc2.Get("h1d_kend_true_inel")
kend_true_el2=f_mc2.Get("h1d_kend_true_el")

#calo
kend_calo_el2=f_mc2.Get("h1d_kend_calo_el")
kend_calo_inel2=f_mc2.Get("h1d_kend_calo_inel")

#fit2
fit_keff_el2=VNFit(keff_el2, pre_mean_ff, n_s)
fit_keff_inel2=VNFit(keff_inel2, pre_mean_ff, n_s)
fit_keffbeam_el2=VNFit(keffbeam_el2, pre_mean_ff, n_s)
fit_keffbeam_inel2=VNFit(keffbeam_inel2, pre_mean_ff, n_s)



#set colors ------------------------
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

kend_true_el.SetLineColor(3)
kend_true_el.SetMarkerColor(3)

keff_el.SetLineColor(3)
keff_inel.SetLineColor(3)

fit_keff_el.SetLineColor(3)
fit_keff_el.SetMarkerColor(3)

fit_keff_inel.SetLineColor(3)
fit_keff_inel.SetMarkerColor(3)

fit_keff_el2.SetLineColor(3)
fit_keff_el2.SetMarkerColor(3)

fit_keff_inel2.SetLineColor(3)
fit_keff_inel2.SetMarkerColor(3)


keffbeam_el.SetLineColor(4)
keffbeam_inel.SetLineColor(2)

fit_keffbeam_el.SetLineColor(4)
fit_keffbeam_inel.SetLineColor(4)
fit_keffbeam_el.SetMarkerColor(4)
fit_keffbeam_inel.SetMarkerColor(4)

fit_keffbeam_inel.SetLineColor(2)
fit_keffbeam_inel.SetMarkerColor(2)




kend_calo_inel.SetLineWidth(1)
kend_calo_el.SetLineWidth(1)
kend_bb_inel.SetLineWidth(1)
kend_bb_el.SetLineWidth(1)
kend_true_inel.SetLineWidth(1)
kend_true_el.SetLineWidth(1)
keffbeam_el.SetLineWidth(1)
keffbeam_inel.SetLineWidth(1)
keff_el.SetLineWidth(1)
keff_inel.SetLineWidth(1)

fit_keff_el.SetLineStyle(2)
fit_keff_inel.SetLineStyle(2)
fit_keff_el2.SetLineStyle(2)
fit_keff_inel2.SetLineStyle(2)



kend_calo_inel2.SetLineColor(6)
kend_calo_inel2.SetMarkerColor(6)

kend_bb_inel2.SetLineColor(4)
kend_bb_inel2.SetMarkerColor(4)

kend_calo_el2.SetLineColor(6)
kend_calo_el2.SetMarkerColor(6)

kend_bb_el2.SetLineColor(2)
kend_bb_el2.SetMarkerColor(2)

keffbeam_el2.SetLineColor(2)
keffbeam_inel2.SetLineColor(4)

kend_calo_inel2.SetLineWidth(1)
kend_calo_el2.SetLineWidth(1)
kend_bb_inel2.SetLineWidth(1)
kend_bb_el2.SetLineWidth(1)
kend_true_inel2.SetLineWidth(1)
kend_true_el2.SetLineWidth(1)
keffbeam_el2.SetLineWidth(1)
keffbeam_inel2.SetLineWidth(1)
keff_el2.SetLineWidth(1)
keff_inel2.SetLineWidth(1)

mu_el=[]
er_mu_el=[]
sigma_el=[]
er_sigma_el=[]

for i in range(3):
  fit=fit_keff_el
  if i==0: 
    fit=fit_keff_el
  if i==1: 
    fit=fit_keffbeam_el
  if i==2:
    fit=fit_keffbeam_el2
  m=fit.GetParameter(0)
  er_m=fit.GetParError(0)
  s=fit.GetParameter(1)
  er_s=fit.GetParError(1)

  mu_el.append(m)
  er_mu_el.append(er_m)
  sigma_el.append(s)
  er_sigma_el.append(er_s)
  print("i=",i," m=",m, "s=",s,"\n")

mu_inel=[]
er_mu_inel=[]
sigma_inel=[]
er_sigma_inel=[]

for i in range(3):
  fit=fit_keff_inel
  if i==0: 
    fit=fit_keff_inel
  if i==1: 
    fit=fit_keffbeam_inel
  if i==2:
    fit=fit_keffbeam_inel2
  m=fit.GetParameter(0)
  er_m=fit.GetParError(0)
  s=fit.GetParameter(1)
  er_s=fit.GetParError(1)

  mu_inel.append(m)
  er_mu_inel.append(er_m)
  sigma_inel.append(s)
  er_sigma_inel.append(er_s)
  print("j=",i," m=",m, "s=",s,"\n")


#Summary plot --------------------------------------------------------------------------------------------------------

#el,ff -------------------------------------------------------------------------------------------------------
c0_ff_mc=RT.TCanvas("c0_ff_mc","",1200,900)
c0_ff_mc.Divide(1,1)
c0_ff_mc.cd(1).SetLogy(0)

f2d_ff=RT.TH2D("f2d_ff","", 400, 200, 600, 900, 0, 900)
f2d_ff.SetTitle("Reco Elastic Scattering Proton Candidates; Kinetic Energy at TPC Front Face [MeV]; Events")
f2d_ff.Draw("")

keff_el.Draw("hist same")
fit_keff_el.Draw("same")
keffbeam_el2.Draw("hist same")
fit_keffbeam_el2.Draw("same")
keffbeam_el.Draw("hist same")
fit_keffbeam_el.Draw("same")

txt_el=[]
txt_el.append("El.(reco.) [Edept.]: #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_el[1], er_mu_el[1], sigma_el[1], er_sigma_el[1]))
txt_el.append("El.(reco.) [Const.]: #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_el[2], er_mu_el[2], sigma_el[2], er_sigma_el[2]))
txt_el.append("El.(truth): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_el[0], er_mu_el[0], sigma_el[0], er_sigma_el[0]))

leg_ff=RT.TLegend(0.16, 0.65, .9, 0.87)
leg_ff.SetFillStyle(0)
leg_ff.AddEntry(keffbeam_el, txt_el[0], "l")
leg_ff.AddEntry(keffbeam_el2, txt_el[1], "l")
leg_ff.AddEntry(keff_el, txt_el[2], "l")
#leg_ff.AddEntry(keffbeam_el, "El.(reco.) [KE_{FF}=E_{Edept}]", "l")
#leg_ff.AddEntry(keffbeam_el2, "El.(reco.) [KE_{FF}=E_{Const}]", "l")
#leg_ff.AddEntry(keff_el, "El.(truth)", "l")
leg_ff.Draw()
c0_ff_mc.Print(out_path+'/keff_el_comparison.eps')

#inel,ff ------------------------------------------------------------------------------------------------------------
c0_ff_inel_mc=RT.TCanvas("c0_ff_inel_mc","",1200,900)
c0_ff_inel_mc.Divide(1,1)
c0_ff_inel_mc.cd(1).SetLogy(0)

f2d_ff_inel=RT.TH2D("f2d_ff_inel","", 400, 200, 600, 1000, 0, 1000)
f2d_ff_inel.SetTitle("Reco Inelastic Scattering Proton Candidates; Kinetic Energy at TPC Front Face [MeV]; Events")
f2d_ff_inel.Draw("")

keff_inel.Draw("hist same")
fit_keff_inel.Draw("same")
keffbeam_inel2.Draw("hist same")
fit_keffbeam_inel2.Draw("same")
keffbeam_inel.Draw("hist same")
fit_keffbeam_inel.Draw("same")

txt_inel=[]
txt_inel.append("Inel.(reco.) [Edept.]: #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_inel[1], er_mu_inel[1], sigma_inel[1], er_sigma_inel[1]))
txt_inel.append("Inel.(reco.) [Const.]: #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_inel[2], er_mu_inel[2], sigma_inel[2], er_sigma_inel[2]))
txt_inel.append("Inel.(truth): #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_inel[0], er_mu_inel[0], sigma_inel[0], er_sigma_inel[0]))
leg_inel_ff2=RT.TLegend(0.16, 0.68, .9, 0.9)
leg_inel_ff2.SetFillStyle(0)
leg_inel_ff2.AddEntry(keffbeam_inel, txt_inel[0], "l")
leg_inel_ff2.AddEntry(keffbeam_inel2, txt_inel[1], "l")
leg_inel_ff2.AddEntry(keff_inel, txt_inel[2], "l")

#leg_inel_ff2.AddEntry(keffbeam_inel, "Inel.(reco.) [KE_{FF}=E_{Edept}]", "l")
#leg_inel_ff2.AddEntry(keffbeam_inel2, "Inel.(reco.) [KE_{FF}=E_{Const}]", "l")
#leg_inel_ff2.AddEntry(keff_inel, "Inel.(truth)", "l")
leg_inel_ff2.Draw()
c0_ff_inel_mc.Print(out_path+'/keff_inel_comparison.eps')


#el,end ----------------------------------------------------------------------------------------------------
c0_end_mc=RT.TCanvas("c0_end_mc","",1200,900)
c0_end_mc.Divide(1,1)
c0_end_mc.cd(1).SetLogy(0)

f2d_mc=RT.TH2D("f2d_mc","", 700, -100, 600, 900, 0, 900)
f2d_mc.SetTitle("Reco Elastic Scattering Proton Candidates; Kinetic Energy at Track End [MeV]; Events")
f2d_mc.Draw("")

kend_true_el.Draw("hist same")
kend_bb_el2.Draw("hist same")
kend_bb_el.Draw("hist same")

leg=RT.TLegend(0.3,0.6,.8,0.87)
leg.SetFillStyle(0)
leg.AddEntry(kend_bb_el, "El.(reco.) [Edept.]", "l")
leg.AddEntry(kend_bb_el2, "El.(reco.) [Const.]", "l")
leg.AddEntry(kend_true_el, "El.(truth)", "l")
leg.Draw()
c0_end_mc.Print(out_path+'/kend_el_comparison.eps')

#inel,end ----------------------------------------------------------------------------------------------------
c0_end_inel_mc=RT.TCanvas("c0_end_inel_mc","",1200,900)
c0_end_inel_mc.Divide(1,1)
c0_end_inel_mc.cd(1).SetLogy(0)

f2d_inel_mc=RT.TH2D("f2d_inel_mc","", 700, -100, 600, 900, 0, 900)
f2d_inel_mc.SetTitle("Reco Inelastic Scattering Proton Candidates; Kinetic Energy at Track End [MeV]; Events")
f2d_inel_mc.Draw("")

kend_true_inel.Draw("hist same")
kend_bb_inel2.Draw("hist same")
kend_bb_inel.Draw("hist same")

leg_inel=RT.TLegend(0.3,0.6,.8,0.87)
leg_inel.SetFillStyle(0)
leg_inel.AddEntry(kend_bb_inel, "Inel.(reco.) [Edept.]", "l")
leg_inel.AddEntry(kend_bb_inel2, "Inel.(reco.) [Const.]", "l")
leg_inel.AddEntry(kend_true_inel, "Inel.(truth)", "l")
leg_inel.Draw()
c0_end_inel_mc.Print(out_path+'/kend_inel_comparison.eps')



