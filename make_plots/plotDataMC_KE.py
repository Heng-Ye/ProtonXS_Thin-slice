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

#ff-const E-loss from kebeam
keffbeam_mc=f_mc.Get("h1d_keffbeam")
keffbeam_stop_mc=f_mc.Get("h1d_keffbeam_stop")
keffbeam_inel_mc=f_mc.Get("h1d_keffbeam_inel")

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


#fitting to extract mu and sigma ----------------------------
#Data
#initial KE
fit_kebeam_stop_data=VNFit(kebeam_stop_data, 430, 3)

#FF_constEloss
fit_keffbeam_stop_data=VNFit(keffbeam_stop_data, 410, 3)
fit_keffbeam_inel_data=VNFit(keffbeam_inel_data, 410, 3)

#FF
fit_kehy_stop_data=VNFit(kehy_stop_data, 410, 3)
fit_kehy_inel_data=VNFit(kehy_inel_data, 410, 3)

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

#FF-const E-loss
fit_keffbeam_mc=VNFit(keffbeam_mc, 400, 3)
fit_keffbeam_stop_mc=VNFit(keffbeam_stop_mc, 400, 3)
fit_keffbeam_inel_mc=VNFit(keffbeam_inel_mc, 400, 3)

#FF-HY
fit_kehy_stop_mc=VNFit(kehy_stop_mc, 400, 3)
fit_kehy_inel_mc=VNFit(kehy_inel_mc, 400, 3)

#EDept
fit_kerange_stop_mc=VNFit(kerange_stop_mc, 420, 3)
fit_kecalo_stop_mc=VNFit(kecalo_stop_mc, 420, 3)

#End_true
fit_kend_true_stop_mc=VNFit(kend_true_stop_mc, 20, 3)
fit_kend_true_el_mc=VNFit(kend_true_el_mc, 20, 3)
fit_kend_true_inel_mc=VNFit(kend_true_inel_mc, 20, 3)

#End_calo
fit_kend_calo_stop_mc=VNFit(kend_calo_stop_mc , 20, 3)
fit_kend_calo_el_mc=VNFit(kend_calo_el_mc , 20, 3)
fit_kend_calo_inel_mc=VNFit(kend_calo_inel_mc , 20, 3)
fit_kend_bb_stop_mc=VNFit(kend_bb_stop_mc , 20, 3)
fit_kend_bb_el_mc=VNFit(kend_bb_el_mc , 20, 3)
fit_kend_bb_inel_mc=VNFit(kend_bb_inel_mc , 20, 3)

#get fitted mu and sigma [mc] ------------------------
n_mc_stop=10

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
    fit=fit_kend_bb_el_mc
  if i==9:
    fit=fit_kend_true_stop_mc
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
    

print("n_mc_stop=",n_mc_stop,"\n")
#print("mu_sg_mc::",mu_sg_mc)
#print("len(mu_sg_mc):",mu_sg_mc)
#print("len(mu_sg_mc[0]):",len(mu_sg_mc[0]))
print("const_eloss_mc:",const_eloss_mc)
	
#get fitted mu and sigma [data] -----------------
n_data_stop=7

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
    fit=fit_kend_bb_el_data
  m=fit.GetParameter(0)
  er_m=fit.GetParError(0)
  s=fit.GetParameter(1)
  er_s=fit.GetParError(1)
  mu_stop_data.append(m)
  er_mu_stop_data.append(er_m)
  sigma_stop_data.append(s)
  er_sigma_stop_data.append(er_s)

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









#FF
test=[]
ex=[0,0,0]

#test1 = [1, 2, 3]
#test.append([4, 5, 6])

#test = array.array('d', [1, 2, 3])
test.append([1, 2, 3]) 
test.append([4, 5, 6]) 
print(test)
print(len(test))
print(len(test[0]))



gr_stop_mc = RT.TGraphErrors(len(mu_stop_mc), array('d', mu_stop_mc[:]), array('d', sigma_stop_mc[:]), array('d', er_mu_stop_mc[:]), array('d', er_sigma_stop_mc[:]))

cx_stop_mc=RT.TCanvas("cx_stop_mc","",1200,900)
cx_stop_mc.Divide(1,1)
cx_stop_mc.cd(1)
gr_stop_mc.Draw("ap")
cx_stop_mc.Print('test.eps')





