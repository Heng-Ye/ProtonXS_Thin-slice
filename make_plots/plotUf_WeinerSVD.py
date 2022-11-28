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
RT.gROOT.SetBatch(); #When running in batch mode, PyROOT does NOT display any graphics
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

#plot output
plot_out='./plots_Weiner_SVD/'

#parser.add_argument("-d_int", type=str, help='Data INT File', default = "")
#args=parser.parse_args()

#if not (args.d and args.o and args.c):
#  print("--> Please provide all files (data, MC, MR_rw, output_file_folder). Thank you!")
#  exit()
#if (args.d): print('Read data: '+args.d)


#Load histograms --------------------------------------
#before unfolding
f_b_int=RT.TFile(arg_b_int, "OPEN")
true_int=f_b_int.Get("htrue_signal")
reco_int=f_b_int.Get("hmeas")
reco_nobkgsub_int=f_b_int.Get("hmeas_before_bkgsub")

f_b_inc=RT.TFile(arg_b_inc, "OPEN")
true_inc=f_b_inc.Get("htrue_signal")
reco_inc=f_b_inc.Get("hmeas")
reco_nobkgsub_inc=f_b_inc.Get("hmeas_before_bkgsub")

f_b_inc_st=RT.TFile(arg_b_inc_st, "OPEN")
true_inc_st=f_b_inc_st.Get("htrue_signal")
reco_inc_st=f_b_inc_st.Get("hmeas")
reco_nobkgsub_inc_st=f_b_inc_st.Get("hmeas_before_bkgsub")


#after Weiner SVD unfolding
f_uf_int=RT.TFile(arg_uf_int, "OPEN")
unf_int=f_uf_int.Get("unf")

f_uf_inc=RT.TFile(arg_uf_inc, "OPEN")
unf_inc=f_uf_inc.Get("unf")

f_uf_inc_st=RT.TFile(arg_uf_inc_st, "OPEN")
unf_inc_st=f_uf_inc_st.Get("unf")



#set colors ------------------------
#
true_int.SetLineColor(2)
true_int.SetMarkerColor(2)

unf_int.SetLineColor(4)
unf_int.SetMarkerColor(4)

reco_int.SetLineColor(3)
reco_int.SetMarkerColor(3)

#
true_inc.SetLineColor(2)
true_inc.SetMarkerColor(2)

unf_inc.SetLineColor(4)
unf_inc.SetMarkerColor(4)

reco_inc.SetLineColor(3)
reco_inc.SetMarkerColor(3)

#
true_inc_st.SetLineColor(2)
true_inc_st.SetMarkerColor(2)

unf_inc_st.SetLineColor(4)
unf_inc_st.SetMarkerColor(4)

reco_inc_st.SetLineColor(3)
reco_inc_st.SetMarkerColor(3)



#el,ff -------------------------------------------------------------------------------------------------------
c0_int=RT.TCanvas("c0_int","",1200,900)
c0_int.Divide(1,1)
c0_int.cd(1).SetLogy(0)

f2d_int=RT.TH2D("f2d_ff","", 42,-1,41,100,-300,4000)
f2d_int.SetTitle("Interacting Protons; Slice ID; ")
f2d_int.Draw("")

true_int.Draw("hist same")
reco_int.Draw("ep same")
unf_int.Draw("ep same")

#txt_int=[]
#txt_el.append("El.(reco.) [Const.]: #mu={:.1f}#pm{:.1f} MeV, #sigma={:.1f}#pm{:.1f} MeV".format(mu_el[2], er_mu_el[2], sigma_el[2], er_sigma_el[2]))

leg_int=RT.TLegend(0.16, 0.65, .9, 0.87)
leg_int.SetFillStyle(0)
leg_int.AddEntry(true_int, 'Truth', "l")
leg_int.AddEntry(reco_int, 'Reco (before Weiner SVD unfolding)', "l")
leg_int.AddEntry(unf_int, 'Reco (after Weiner SVD unfolding)', "l")
leg_int.Draw()
c0_int.Print(plot_out+'spec_int.eps')




c0_inc=RT.TCanvas("c0_inc","",1200,900)
c0_inc.Divide(1,1)
c0_inc.cd(1).SetLogy(0)

f2d_inc=RT.TH2D("f2d_inc","", 42,-1,41,100,-300,4000)
f2d_inc.SetTitle("Incident Protons; Slice ID; ")
f2d_inc.Draw("")

true_inc.Draw("hist same")
reco_inc.Draw("ep same")
unf_inc.Draw("ep same")

leg_inc=RT.TLegend(0.16, 0.65, .9, 0.87)
leg_inc.SetFillStyle(0)
leg_inc.AddEntry(true_inc, 'Truth', "l")
leg_inc.AddEntry(reco_inc, 'Reco (before Weiner SVD unfolding)', "l")
leg_inc.AddEntry(unf_inc, 'Reco (after Weiner SVD unfolding)', "l")
leg_inc.Draw()
c0_inc.Print(plot_out+'spec_inc.eps')



c0_inc_st=RT.TCanvas("c0_inc_st","",1200,900)
c0_inc_st.Divide(1,1)
c0_inc_st.cd(1).SetLogy(0)

f2d_inc_st=RT.TH2D("f2d_inc_st","", 42,-1,41,100,-800,15000)
f2d_inc_st.SetTitle("Incident Protons; Start Slice ID; ")
f2d_inc_st.Draw("")

true_inc_st.Draw("hist same")
reco_inc_st.Draw("ep same")
unf_inc_st.Draw("ep same")

leg_inc_st=RT.TLegend(0.16, 0.65, .9, 0.87)
leg_inc_st.SetFillStyle(0)
leg_inc_st.AddEntry(true_inc, 'Truth', "l")
leg_inc_st.AddEntry(reco_inc, 'Reco (before Weiner SVD unfolding)', "l")
leg_inc_st.AddEntry(unf_inc, 'Reco (after Weiner SVD unfolding)', "l")
leg_inc_st.Draw()
c0_inc_st.Print(plot_out+'spec_inc_st.eps')

