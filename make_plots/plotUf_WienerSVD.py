import ROOT as RT
import sys
import math

from argparse import ArgumentParser as ap
import numpy as np
from array import array
from math import sqrt
from math import exp

#--------------------------------------------------------------------------------------
parser = ap()
RT.gROOT.LoadMacro("~/protoDUNEStyle.C")
RT.gROOT.SetBatch(); #When running in batch mode, PyROOT will NOT display any graphics
RT.gStyle.SetOptStat(00000)
#RT.gStyle.SetErrorX(1.e-4)
RT.gStyle.SetTitleAlign(23)
RT.gStyle.SetTitleX(.5)
RT.gStyle.SetLineWidth(1)
tt = RT.TLatex();
tt.SetNDC();

#Read files -----------------------------------------------------------------------------------------
path_b='../Wiener_SVD_files/MC_nobmrw/'
path_uf=path_b

#before uf files
arg_b_int=path_b+'input_wiener_svd_int.root'
arg_b_inc=path_b+'input_wiener_svd_inc.root'
arg_b_inc_st=path_b+'input_wiener_svd_inc_st.root'

#uf files
arg_uf_int=path_uf+'uf_wiener_svd_int.root'
arg_uf_inc=path_uf+'uf_wiener_svd_inc.root'
arg_uf_inc_st=path_uf+'uf_wiener_svd_inc_st.root'

#uf using Bayesian unfold
#arg_uf_roo_int='../RooUnfold_files/MC_nobmrw/input_uf_roounfold_int.root'
#arg_uf_roo_inc='../RooUnfold_files/MC_nobmrw/input_uf_roounfold_inc.root'
#arg_uf_roo_inc_st='../RooUnfold_files/MC_nobmrw/input_uf_roounfold_inc_st.root'

arg_uf_roo_int='../RooUnfold_files/MC_nobmrw/SVD/input_uf_roounfold_int.root'
arg_uf_roo_inc='../RooUnfold_files/MC_nobmrw/SVD/input_uf_roounfold_inc.root'
arg_uf_roo_inc_st='../RooUnfold_files/MC_nobmrw/SVD/input_uf_roounfold_inc_st.root'

#plot output
plot_out='./plots_Wiener_SVD/'

#Load histograms -------------------------------------------
#before unfolding
#int
f_b_int=RT.TFile(arg_b_int, "OPEN")
true_int=f_b_int.Get("htrue_signal")
reco_int=f_b_int.Get("hmeas")
reco_nobkgsub_int=f_b_int.Get("hmeas_before_bkgsub")

n_true_int=true_int.Integral()
n_reco_int=reco_int.Integral()

#inc
f_b_inc=RT.TFile(arg_b_inc, "OPEN")
true_inc=f_b_inc.Get("htrue_signal")
reco_inc=f_b_inc.Get("hmeas")
reco_nobkgsub_inc=f_b_inc.Get("hmeas_before_bkgsub")
n_true_inc=true_inc.Integral()
n_reco_inc=reco_inc.Integral()

#inc_st
f_b_inc_st=RT.TFile(arg_b_inc_st, "OPEN")
true_inc_st=f_b_inc_st.Get("htrue_signal")
reco_inc_st=f_b_inc_st.Get("hmeas")
reco_nobkgsub_inc_st=f_b_inc_st.Get("hmeas_before_bkgsub")
n_true_inc_st=true_inc_st.Integral()
n_reco_inc_st=reco_inc_st.Integral()

#after Weiner SVD unfolding
#int
f_uf_int=RT.TFile(arg_uf_int, "OPEN")
unf_int=f_uf_int.Get("unf")
n_uf_int=unf_int.Integral()
print('n_true_int:',n_true_int)
print('n_reco_int:',n_reco_int)
print('n_uf_int:',n_uf_int)
print('n_reco_int/n_uf_int:',n_reco_int/n_uf_int)
#unf_int.Scale(n_reco_int/n_uf_int)

#inc
f_uf_inc=RT.TFile(arg_uf_inc, "OPEN")
unf_inc=f_uf_inc.Get("unf")
n_uf_inc=unf_inc.Integral()
print('n_true_inc:',n_true_inc)
print('n_reco_inc:',n_reco_inc)
print('n_uf_inc:',n_uf_inc)
print('n_reco_inc/n_uf_inc:',n_reco_inc/n_uf_inc)
#unf_inc.Scale(n_reco_inc/n_uf_inc)

#inc_st
f_uf_inc_st=RT.TFile(arg_uf_inc_st, "OPEN")
unf_inc_st=f_uf_inc_st.Get("unf")
n_uf_inc_st=unf_inc_st.Integral()
print('n_true_inc_st:',n_true_inc_st)
print('n_reco_inc_st:',n_reco_inc_st)
print('n_uf_inc_st:',n_uf_inc_st)
print('n_reco_inc/n_uf_inc:',n_reco_inc_st/n_uf_inc_st)
#unf_inc_st.Scale(n_reco_inc_st/n_uf_inc_st)

#after Roounfold -----------------------------------------
#int
f_uf_roo_int=RT.TFile(arg_uf_roo_int, "OPEN")
uf_bayes_int=f_uf_roo_int.Get("data_int_uf")
n_uf_bayes_int=uf_bayes_int.Integral()
#uf_bayes_int.Scale(n_reco_int/n_uf_bayes_int)

#inc
f_uf_roo_inc=RT.TFile(arg_uf_roo_inc, "OPEN")
uf_bayes_inc=f_uf_roo_inc.Get("data_inc_uf")
n_uf_bayes_inc=uf_bayes_inc.Integral()
#uf_bayes_inc.Scale(n_reco_inc/n_uf_bayes_inc)

#inc_st
f_uf_roo_inc_st=RT.TFile(arg_uf_roo_inc_st, "OPEN")
uf_bayes_inc_st=f_uf_roo_inc_st.Get("data_st_inc_uf")
n_uf_bayes_inc_st=uf_bayes_inc_st.Integral()
#uf_bayes_inc_st.Scale(n_reco_inc_st/n_uf_bayes_inc_st)


#set colors ------------------------
#int
true_int.SetLineColor(2)
true_int.SetMarkerColor(2)

unf_int.SetLineColor(3)
unf_int.SetMarkerColor(3)

uf_bayes_int.SetLineColor(4)
uf_bayes_int.SetMarkerColor(4)
#uf_bayes_int.SetMarkerStyle(24)

reco_int.SetLineColor(13)
reco_int.SetMarkerColor(13)
reco_int.SetMarkerStyle(24)


#inc
true_inc.SetLineColor(2)
true_inc.SetMarkerColor(2)

unf_inc.SetLineColor(3)
unf_inc.SetMarkerColor(3)

reco_inc.SetLineColor(13)
reco_inc.SetMarkerColor(13)
reco_inc.SetMarkerStyle(24)

uf_bayes_inc.SetLineColor(4)
uf_bayes_inc.SetMarkerColor(4)
#uf_bayes_inc.SetMarkerStyle(24)


#inc_st
true_inc_st.SetLineColor(2)
true_inc_st.SetMarkerColor(2)

unf_inc_st.SetLineColor(3)
unf_inc_st.SetMarkerColor(3)

reco_inc_st.SetLineColor(13)
reco_inc_st.SetMarkerColor(13)
reco_inc_st.SetMarkerStyle(24)

uf_bayes_inc_st.SetLineColor(4)
uf_bayes_inc_st.SetMarkerColor(4)
#uf_bayes_inc_st.SetMarkerStyle(24)


#Plots ------------------------------------------------------------------------------------------------------
#int
c0_int=RT.TCanvas("c0_int","",1200,900)
c0_int.Divide(1,1)
c0_int.cd(1).SetLogy(0)

f2d_int=RT.TH2D("f2d_ff","", 42,-1,41,100,-300,5000)
f2d_int.SetTitle("Interacting Protons; Slice ID; ")
f2d_int.Draw("")

true_int.Draw("hist same")
reco_int.Scale(n_uf_int/n_reco_int)
reco_int.Draw("ep same")
uf_bayes_int.Draw("ep same")
unf_int.Draw("ep same")

leg_int=RT.TLegend(0.16, 0.65, .9, 0.87)
leg_int.SetFillStyle(0)
leg_int.AddEntry(true_int, 'Truth (Validation set)', "l")
#leg_int.AddEntry(reco_int, 'Reco before unfolding (Test set)', "ep")
leg_int.AddEntry(unf_int, 'Reco after Wiener SVD unfolding (Test set)', "ep")
#leg_int.AddEntry(uf_bayes_int, 'Reco after Bayesian unfolding (Test set)', "ep")
leg_int.AddEntry(uf_bayes_int, 'Reco after SVD unfolding (Test set)', "ep")
leg_int.Draw()
c0_int.Print(plot_out+'spec_int.eps')

#inc
c0_inc=RT.TCanvas("c0_inc","",1200,900)
c0_inc.Divide(1,1)
c0_inc.cd(1).SetLogy(0)

f2d_inc=RT.TH2D("f2d_inc","", 42,-1,41,100,-500,5000)
f2d_inc.SetTitle("Incident Protons; Slice ID; ")
f2d_inc.Draw("")

true_inc.Draw("hist same")
reco_inc.Scale(n_uf_inc/n_reco_inc)
reco_inc.Draw("ep same")
uf_bayes_inc.Draw("ep same")
unf_inc.Draw("ep same")

leg_inc=RT.TLegend(0.16, 0.65, .9, 0.87)
leg_inc.SetFillStyle(0)
leg_inc.AddEntry(true_inc, 'Truth (Validation Set)', "l")
#leg_inc.AddEntry(reco_inc, 'Reco (before unfolding)', "l")
leg_inc.AddEntry(unf_inc, 'Reco after Wiener SVD Unfolding (Test set)', "ep")
#leg_inc.AddEntry(uf_bayes_inc, 'Reco after Bayesian Unfolding (Test set)', "ep")
leg_inc.AddEntry(uf_bayes_inc, 'Reco after SVD Unfolding (Test set)', "ep")
leg_inc.Draw()
c0_inc.Print(plot_out+'spec_inc.eps')

#inc_st
c0_inc_st=RT.TCanvas("c0_inc_st","",1200,900)
c0_inc_st.Divide(1,1)
c0_inc_st.cd(1).SetLogy(0)

f2d_inc_st=RT.TH2D("f2d_inc_st","", 42,-1,41,100,-800,15000)
f2d_inc_st.SetTitle("Incident Protons; Start Slice ID; ")
f2d_inc_st.Draw("")

true_inc_st.Draw("hist same")
reco_inc_st.Scale(n_uf_inc_st/n_reco_inc_st)
reco_inc_st.Draw("ep same")
uf_bayes_inc_st.Draw("ep same")
unf_inc_st.Draw("ep same")

leg_inc_st=RT.TLegend(0.16, 0.65, .9, 0.87)
leg_inc_st.SetFillStyle(0)
leg_inc_st.AddEntry(true_inc, 'Truth (Validation Set)', "l")
#leg_inc_st.AddEntry(reco_inc, 'Reco (before unfolding)', "l")
leg_inc_st.AddEntry(unf_inc, 'Reco after Wiener SVD unfolding (Test set)', "ep")
#leg_inc_st.AddEntry(uf_bayes_inc_st, 'Reco after Bayesian unfolding (Test set)', "ep")
leg_inc_st.AddEntry(uf_bayes_inc_st, 'Reco after SVD unfolding (Test set)', "ep")
leg_inc_st.Draw()
c0_inc_st.Print(plot_out+'spec_inc_st.eps')

