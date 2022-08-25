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
file_mc='../mc_proton_beamxy_beammom_nobmrw.root'
#file_mc='../mc_proton_beamxy_beammom_nobmrw_old.root'
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

#2D scatter plots
kebb_keendtruth_stop=f_mc.Get("h2d_kebb_keendtruth_stop")
kebb_keendtruth_el=f_mc.Get("h2d_kebb_keendtruth_el")
kebb_keendtruth_inel=f_mc.Get("h2d_kebb_keendtruth_inel")

kebb_kebbtruth_stop=f_mc.Get("h2d_kebb_kebbtruth_stop")
kebb_kebbtruth_el=f_mc.Get("h2d_kebb_kebbtruth_el")
kebb_kebbtruth_inel=f_mc.Get("h2d_kebb_kebbtruth_inel")

kecalo_keendtruth_stop=f_mc.Get("h2d_kecalo_keendtruth_stop")
kecalo_keendtruth_el=f_mc.Get("h2d_kecalo_keendtruth_el")
kecalo_keendtruth_inel=f_mc.Get("h2d_kecalo_keendtruth_inel")

kebbtruelength_keendtruth_stop=f_mc.Get("h2d_kebbtruelength_keendtruth_stop")
kebbtruelength_keendtruth_el=f_mc.Get("h2d_kebbtruelength_keendtruth_el")
kebbtruelength_keendtruth_inel=f_mc.Get("h2d_kebbtruelength_keendtruth_inel")

kend_bbtruelength_stop=f_mc.Get("h1d_kend_bbtruelength_stop")
kend_bbtruelength_el=f_mc.Get("h1d_kend_bbtruelength_el")
kend_bbtruelength_inel=f_mc.Get("h1d_kend_bbtruelength_inel")

ratio_rangetrue_rangereco_stop=f_mc.Get("h1d_ratio_rangetrue_rangereco_stop")
ratio_rangetrue_rangereco_el=f_mc.Get("h1d_ratio_rangetrue_rangereco_el")
ratio_rangetrue_rangereco_inel=f_mc.Get("h1d_ratio_rangetrue_rangereco_inel")

recorange_truerange_stop=f_mc.Get("h2d_recorange_truerange_stop")
recorange_truerange_el=f_mc.Get("h2d_recorange_truerange_el")
recorange_truerange_inel=f_mc.Get("h2d_recorange_truerange_inel")

ratio_rangehy_rangereco_stop=f_mc.Get("h1d_ratio_rangehy_rangereco_stop")
ratio_rangehy_rangereco_el=f_mc.Get("h1d_ratio_rangehy_rangereco_el")

recorange_rangehy_el=f_mc.Get("h2d_recorange_rangehy_el")
recorange_rangehy_stop=f_mc.Get("h2d_recorange_rangehy_stop")


#[1]------------------------------------------------------------------------------------------------------
c0_end_mc=RT.TCanvas("c0_end_mc","",1200,900)
c0_end_mc.Divide(1,1)
c0_end_mc.cd(1).SetLogy(0)

#f2d_mc=RT.TH2D("f2d_mc","", 700, -100, 600, 900, 0.01, 2000) #logy
#f2d_mc=RT.TH2D("f2d_mc","", 800, -200, 600, 900, 0., 1000) #liny
f2d_mc=RT.TH2D("f2d_mc","", 800, -200, 600, 900, 0., 1500) #liny+bbtruth_length

f2d_mc.SetTitle("Reconstructed Elastic-scattering Protons; Kinetic Energy at Track End [MeV]; Events")
f2d_mc.Draw("")

kend_calo_el.SetLineColor(2)
kend_calo_el.SetMarkerColor(2)
kend_bb_el.SetLineColor(4)
kend_bb_el.SetMarkerColor(4)
kend_true_el.SetLineColor(3)
kend_true_el.SetMarkerColor(3)
kend_bbtruelength_el.SetLineColor(51)
kend_bbtruelength_el.SetMarkerColor(51)

kend_true_el.Draw("hist same")
kend_calo_el.Draw("hist same")
kend_bb_el.Draw("hist same")
kend_bbtruelength_el.Draw("hist same")

leg0=RT.TLegend(0.45,0.6,.87,0.87)
leg0.SetFillStyle(0)
leg0.AddEntry(kend_calo_el, "Calo.", "l")
leg0.AddEntry(kend_bb_el, "BB", "l")
leg0.AddEntry(kend_true_el, "Truth", "l")
leg0.AddEntry(kend_bbtruelength_el, "BB with Range(truth)", "l")

#leg0.SetNColumns(2);
leg0.Draw()


#KEffbeam_KEhy_stop.Draw("colz same")
#ll=RT.TLine(200,200,600,600)
#ll.SetLineColor(2)
#ll.SetLineStyle(2)
#ll.Draw()
#c0_end_mc.Print(out_path+'/keend_el.eps')
c0_end_mc.Print(out_path+'/keend__el.eps')


#[2]------------------------------------------------------------------------------------------------------
c1_end_mc=RT.TCanvas("c1_end_mc","",1200,900)
c1_end_mc.Divide(1,1)
c1_end_mc.cd(1).SetLogy(0)

f2d_mc=RT.TH2D("f2d_mc","", 800, -200, 600, 900, 0., 400) #liny
f2d_mc.SetTitle("Reconstructed Inelastic-scattering Protons; Kinetic Energy at Track End [MeV]; Events")
f2d_mc.Draw("")

kend_calo_inel.SetLineColor(2)
kend_calo_inel.SetMarkerColor(2)
kend_bb_inel.SetLineColor(4)
kend_bb_inel.SetMarkerColor(4)
kend_true_inel.SetLineColor(3)
kend_true_inel.SetMarkerColor(3)
kend_bbtruelength_inel.SetLineColor(51)
kend_bbtruelength_inel.SetMarkerColor(51)


kend_true_inel.Draw("hist same")
kend_calo_inel.Draw("hist same")
kend_bb_inel.Draw("hist same")
kend_bbtruelength_inel.Draw("hist same")


#leg1=RT.TLegend(0.7,0.7,.87,0.9)
leg1=RT.TLegend(0.45,0.6,.87,0.87)
leg1.SetFillStyle(0)
leg1.AddEntry(kend_calo_inel, "Calo.", "l")
leg1.AddEntry(kend_bb_inel, "BB.", "l")
#leg1.AddEntry(kend_true_inel, "Truth", "l")
leg1.AddEntry(kend_bbtruelength_el, "BB with Range(truth)", "l")
#leg0.SetNColumns(2);
leg1.Draw()


#KEffbeam_KEhy_stop.Draw("colz same")
#ll=RT.TLine(200,200,600,600)
#ll.SetLineColor(2)
#ll.SetLineStyle(2)
#ll.Draw()
#c1_end_mc.Print(out_path+'/keend_inel.eps')
c1_end_mc.Print(out_path+'/keend__inel.eps')


#[3]--------------------------------------------------------------------------------------------------------------------------------------
#stop
c2_end_mc=RT.TCanvas("c2_end_mc","",1200,900)
c2_end_mc.Divide(1,1)
c2_end_mc.cd(1).SetLogy(0)

f2d2_mc=RT.TH2D("f2d_mc","", 800, -200, 600, 800, -200, 600) #liny
f2d2_mc.SetTitle("Reconstructed Stopping Protons; Kinetic Energy at Track End (BB) [MeV]; Kinetic Energy at Track End (Truth) [MeV]")
f2d2_mc.Draw("")

kebb_keendtruth_stop.SetMarkerSize(.3)
kebb_keendtruth_stop.Draw("same")

ll=RT.TLine(-200,-200,600,600)
ll.SetLineColor(2)
ll.SetLineStyle(2)
ll.Draw()
c2_end_mc.Print(out_path+'/kebb_kend_truth_stop.eps')

#el
f2d2_mc.SetTitle("Reconstructed Elastic Protons; Kinetic Energy at Track End (BB) [MeV]; Kinetic Energy at Track End (Truth) [MeV]")
f2d2_mc.Draw("")

kebb_keendtruth_el.SetMarkerSize(.3)
kebb_keendtruth_el.Draw("same")
ll.Draw()
c2_end_mc.Print(out_path+'/kebb_kend_truth_el.eps')

#inel
f2d2_mc.SetTitle("Reconstructed Inelastic Protons; Kinetic Energy at Track End (BB) [MeV]; Kinetic Energy at Track End (Truth) [MeV]")
f2d2_mc.Draw("")

kebb_keendtruth_inel.SetMarkerSize(.3)
kebb_keendtruth_inel.Draw("same")
ll.Draw()
c2_end_mc.Print(out_path+'/kebb_kend_truth_inel.eps')

#[4]------------------------------------------------------------------------------------------------------
f2d2_mc.SetTitle("Reconstructed Stopping Protons; (KE_{beam}-#DeltaE)*R-E_{dept} [MeV]; KE at Track End (Truth) [MeV]")
f2d2_mc.Draw("")
kecalo_keendtruth_stop.SetMarkerSize(.3)
kecalo_keendtruth_stop.Draw("same")
ll.Draw()
c2_end_mc.Print(out_path+'/kecalo_keend_truth_stop.eps')

f2d2_mc.SetTitle("Reconstructed Elastic-scattering Protons; (KE_{beam}-#DeltaE)*R-E_{dept} [MeV]; KE at Track End (Truth) [MeV]")
f2d2_mc.Draw("")
kecalo_keendtruth_el.SetMarkerSize(.3)
kecalo_keendtruth_el.Draw("same")
ll.Draw()
c2_end_mc.Print(out_path+'/kecalo_keend_truth_el.eps')

f2d2_mc.SetTitle("Reconstructed Inelastic-scattering Protons; (KE_{beam}-#DeltaE)*R-E_{dept} [MeV]; KE at Track End (Truth) [MeV]")
f2d2_mc.Draw("")
kecalo_keendtruth_inel.SetMarkerSize(.3)
kecalo_keendtruth_inel.Draw("same")
ll.Draw()
c2_end_mc.Print(out_path+'/kecalo_keend_truth_inel.eps')

f2d2_mc.SetTitle("Reconstructed Elastic-scattering Protons; KE_{BB} at Track End (True Range) [MeV]; Kinetic Energy at Track End (Truth) [MeV]")
f2d2_mc.Draw()
kebbtruelength_keendtruth_el.SetMarkerSize(.3)
kebbtruelength_keendtruth_el.Draw("same")
ll.Draw()
c2_end_mc.Print(out_path+'/kebbtruelength_keend_truth_el.eps')





#[5]--------------------------------------------------------------------------------------------------------------------------------------
#el
c3_mc=RT.TCanvas("c3_mc","",1200,900)
c3_mc.Divide(1,1)
c3_mc.cd(1).SetLogy(0)

f2d3_mc=RT.TH2D("f2d_mc","", 140,0,140,140,0,140) #liny
f2d3_mc.SetTitle("Reconstructed Stopping Protons; Track Length (Reco) [cm]; Track Length (Truth) [cm]")
f2d3_mc.Draw("")

recorange_truerange_stop.SetMarkerSize(.3)
recorange_truerange_stop.Draw("same")

lll=RT.TLine(0,0,140,140)
lll.SetLineColor(2)
lll.SetLineStyle(2)
lll.Draw()
c3_mc.Print(out_path+'/recorange_truerange_stop.eps')


f2d3_mc.SetTitle("Reconstructed Stopping Protons; Track Length (Reco) [cm]; Track Length (Fit) [cm]")
f2d3_mc.Draw("")
recorange_rangehy_stop.SetMarkerSize(.3)
recorange_rangehy_stop.Draw("same")
lll.Draw()
c3_mc.Print(out_path+'/recorange_hyrange_stop.eps')

#[6]ratio ---------------------------------------------------------------------------------------------------
c4_mc=RT.TCanvas("c4_mc","",1200,900)
c4_mc.Divide(1,1)
c4_mc.cd(1).SetLogy(1)

f2dl_mc=RT.TH2D("f2dl_mc","", 60,-.1,1.7,100,1,25000) #liny
f2dl_mc.SetTitle("Stopping Protons; Range Ratio; Events")
f2dl_mc.Draw("")

#ratio_rangetrue_rangereco_el.SetLineColor(2)
#ratio_rangehy_rangereco_el.SetLineColor(4)

ratio_rangetrue_rangereco_stop.SetLineColor(2)
ratio_rangehy_rangereco_stop.SetLineColor(4)

ratio_rangetrue_rangereco_stop.Draw("hist same")
ratio_rangehy_rangereco_stop.Draw("hist same")

#fit ratio
fit_ratio_rangetrue_rangereco_stop=VRFit(ratio_rangetrue_rangereco_stop, 0.9,1.1)
fit_ratio_rangehy_rangereco_stop=VRFit(ratio_rangehy_rangereco_stop, 0.9,1.1)
fit_ratio_rangetrue_rangereco_stop.SetLineColor(2)
fit_ratio_rangetrue_rangereco_stop.SetLineStyle(2)

fit_ratio_rangehy_rangereco_stop.SetLineColor(4)
fit_ratio_rangehy_rangereco_stop.SetLineStyle(2)

n_mc=2

mu_mc=[]
er_mu_mc=[]
sigma_mc=[]
er_sigma_mc=[]

for i in range(n_mc):
  fit=fit_ratio_rangetrue_rangereco_stop
  if i==0: 
    fit=fit_ratio_rangetrue_rangereco_stop
  if i==1: 
    fit=fit_ratio_rangehy_rangereco_stop
  m=fit.GetParameter(0)
  er_m=fit.GetParError(0)
  s=fit.GetParameter(1)
  er_s=fit.GetParError(1)
  mu_mc.append(m)
  er_mu_mc.append(er_m)
  sigma_mc.append(s)
  er_sigma_mc.append(er_s)
  print("i=",i," m=",m, "s=",s,"\n")


fit_ratio_rangetrue_rangereco_stop.Draw("same")
fit_ratio_rangehy_rangereco_stop.Draw("same")

leg1l=RT.TLegend(0.2,0.7,.7,0.87)
leg1l.SetFillStyle(0)
txt_mc=[]
txt_mc.append("Range(truth)/Range(reco): #mu={:.2f}#pm{:.2f}, #sigma={:.2f}#pm{:.2f}".format(mu_mc[0],er_mu_mc[0],sigma_mc[0],er_sigma_mc[0]))
txt_mc.append("Range(Fit)/Range(reco): #mu={:.2f}#pm{:.2f}, #sigma={:.2f}#pm{:.2f}".format(mu_mc[1],er_mu_mc[1],sigma_mc[1],er_sigma_mc[1]))

leg1l.AddEntry(ratio_rangetrue_rangereco_stop, txt_mc[0], "l")
leg1l.AddEntry(ratio_rangehy_rangereco_stop, txt_mc[1], "l")
leg1l.Draw()

c4_mc.Print(out_path+'/ratios_range_stop.eps')

