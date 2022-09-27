import ROOT as RT
import sys
import math

from argparse import ArgumentParser as ap
import numpy as np
from array import array
from math import sqrt
from math import exp

#Functions ----------------------------------------
def fitg(x, p):
  m=p[0]
  s=p[1]
  n=p[2]
  g=n*math.exp(-(x[0]-m)*(x[0]-m)/(2*s*s))
  return g

def govg(x, par):
  #g1
  m1=par[0]
  s1=par[1]
  g1=-(x[0]-m1)*(x[0]-m1)/(2*s1*s1)

  #g2
  m2=par[2]
  s2=par[3]
  g2=-(x[0]-m2)*(x[0]-m2)/(2*s2*s2)

  #g2/g1
  g_ov_g=0
  g_ov_g=math.exp(g2-g1)
  #if (m1==m2 and s1==s2) g_ov_g=1
  return g_ov_g


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

#Gaussians for weighting function -----------------------------------------------------------------------------------------
#mu_denom_data=411.05145837595467; //new with event-by-event R corr (const E-loss using stopping protons)
#sg_denom_data=47.48714821962207; //new with event-by-event R corr (const E-loss using stopping protons)

mu_denom_data=411.06645442311424 #new with event-by-event (KEHY(fit) using stopping protons)
sg_denom_data=47.076122305960645 #new with event-by-event (KEHY(fit) using stopping protons)

#double mu_nom_mc=416.224743039812; //for mc(KEbeam-const) with R=1 (R=Ratio of KEff(Fit)/(KEbeam-dE))
#double sg_nom_mc=42.786018962508784; //

mu_nom_mc=416.1620092367158 #for mc [KE(Fit)]
sg_nom_mc=40.48356757740762

#weighting func. (ke) --------------------------
xmin=0
xmax=800
kerw=RT.TF1("kerw", govg, xmin, xmax, 4)
kerw.SetParameter(0, mu_nom_mc)
kerw.SetParameter(1, sg_nom_mc)
kerw.SetParameter(2, mu_denom_data)
kerw.SetParameter(3, sg_denom_data)

#MC ---------------------------------------
g_mc=RT.TF1("g_mc", fitg, xmin, xmax, 3)
g_mc.SetName("g_mc")
g_mc.SetParameter(0, mu_nom_mc)
g_mc.SetParameter(1, sg_nom_mc)
g_mc.SetParameter(2, 1)

#Data ----------------------------------------
g_data=RT.TF1("g_data", fitg, xmin, xmax, 3)
g_data.SetName("g_data")
g_data.SetParameter(0, mu_denom_data)
g_data.SetParameter(1, sg_denom_data)
g_data.SetParameter(2, 1)


#KEff(MC): Comparsion of cut effect --------------------------------
#wf ---------------------------------------------------
c0_wf=RT.TCanvas("c0_wf","",1200,900)
c0_wf.Divide(1,1)
c0_wf.cd(1)

f2d_wf=RT.TH2D("f2d_wf","", 800, 0, 800, 12, 0, 120)
f2d_wf.SetTitle("Weighting Function; Proton Kinetic Energy [MeV];")
f2d_wf.GetXaxis().CenterTitle()
f2d_wf.Draw()

kerw.SetLineColor(2)
kerw.SetLineStyle(2)
kerw.Draw("same")

c0_wf.Print('wf/wf.eps')

#gaussian: data and mc ---------------------------------
c1_g=RT.TCanvas("c0_wf","",1200,900)
c1_g.Divide(1,1)
c1_g.cd(1)

f2d_g=RT.TH2D("f2d_g","", 800, 0, 800, 12, 0, 1.2)
f2d_g.SetTitle("; Proton Kinetic Energy [MeV];")
f2d_g.GetXaxis().CenterTitle()
f2d_g.Draw()

g_data.SetLineColor(4)
g_mc.SetLineColor(3)

g_data.SetLineStyle(2)
g_mc.SetLineStyle(2)

g_data.Draw("same")
g_mc.Draw("same")

leg_mc=RT.TLegend(0.14,0.7,.7,0.9)
leg_mc.SetFillStyle(0)
txt_mc=[]
txt_mc.append("MC")
txt_mc.append("Data")
leg_mc.AddEntry(g_mc, txt_mc[0], "l")
leg_mc.AddEntry(g_data, txt_mc[1], "l")
leg_mc.Draw()

c1_g.Print('wf/gaussians.eps')

