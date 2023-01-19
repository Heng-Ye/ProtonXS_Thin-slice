#include "TVector3.h"
#include <TH2.h>
#include <TH1.h>
#include "TH2D.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TLine.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH3F.h>
#include <TString.h>
#include <TProfile2D.h>
#include <THStack.h>
#include "TGraph.h" 
#include "TGraphSmooth.h" 
#include "TParameter.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TSystem.h"
#include "string"
#include "vector"
#include "TSpline.h"
#include "TH3F.h"
#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>

#include <stdio.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>

void plot_xs_sys_uncertainties() {
	//plot output path ---------------------------//
	//TString str_fig_out="./plot_xs_comparison/";
	//TString str_fig_out="./plot_xs_comparison_BB/";
	//TString str_fig_out="./plot_xs_comparison_BMRW/";
	TString str_fig_out="./plot_xs_comparison_El/";
	//TString str_fig_out="./plot_xs_comparison_MisIDP/";

	//Load XS files ----------------------------------------------------------------------------------------------------//
	//reco xs
        //TString str="./xs_files/xs_Eslice_dE20MeV_40slcs_bmrw_data_cen_sysbmrw.root"; //KEff shift due to errors of fit 
        //TString str="./xs_files/xs_Eslice_dE20MeV_40slcs_bmrw_data_cen_sysbeloss.root"; //do not use this one (over-counting!)
        //TString str="./xs_files/xs_Eslice_dE20MeV_40slcs_bmrw_data_cen_sysBB.root";

        //TString str="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_SYS_BB_newslcid.root";
        TString str="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_SYS_ElScale_newslcid.root";
        //TString str="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_SYS_MisIDPScale_newslcid.root";
	TFile *f_in = TFile::Open(str.Data()); //xs
	TGraphErrors *reco_xs=(TGraphErrors *)f_in->Get("reco_xs"); //xs, stat
	TGraphAsymmErrors *reco_xs_all=(TGraphAsymmErrors *)f_in->Get("reco_xs_all"); //xs, err_bmrw, err_fiberpos

        //Geant4 -------------------------------------------------------------------------------------//
        TFile f_xs("/dune/data2/users/hyliao/GeantReweight/xs_cascade/proton_cross_section.root");
        TGraph *total_inel_KE = (TGraph*)f_xs.Get("inel_KE");

	//Neutrino Gen -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        TFile f_other("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/make_plots/model_prediction_screen_grab/model_pred_screengrab.root");
        TGraph *xs_GENIE_ha2018 = (TGraph*)f_other.Get("xs_GENIE_ha2018");
        TGraph *xs_NEUT_2019 = (TGraph*)f_other.Get("xs_NEUT_2019");
        TGraph *xs_NuWRO_2019 = (TGraph*)f_other.Get("xs_NuWRO_2019");
        TGraph *xs_GENIE_hN2018 = (TGraph*)f_other.Get("xs_GENIE_hN2018");
        TGraph *xs_GENIE_INCL_pp = (TGraph*)f_other.Get("xs_GENIE_INCL_pp");
	
        //plot style ---------------------------------------------------------------------------//
        gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
        gROOT->SetStyle("protoDUNEStyle");
        gROOT->ForceStyle();
        gStyle->SetTitleX(0.5);
        gStyle->SetTitleAlign(23);
        gStyle->SetOptStat(0);

	xs_GENIE_ha2018->SetMarkerColor(1); xs_GENIE_ha2018->SetLineColor(1); 
	xs_NEUT_2019->SetMarkerColor(7); xs_NEUT_2019->SetLineColor(7);
	xs_GENIE_INCL_pp->SetMarkerColor(6); xs_GENIE_INCL_pp->SetLineColor(6);
	xs_GENIE_hN2018->SetMarkerColor(kOrange-3); xs_GENIE_hN2018->SetLineColor(kOrange-3);
	xs_NuWRO_2019->SetMarkerColor(4); xs_NuWRO_2019->SetLineColor(4);	

	xs_GENIE_ha2018->SetLineStyle(2);
	xs_NEUT_2019->SetLineStyle(2);
	xs_GENIE_INCL_pp->SetLineStyle(2);
	xs_GENIE_hN2018->SetLineStyle(2);
	xs_NuWRO_2019->SetLineStyle(2);

        //gStyle->SetPadLeftMargin(0.13);
        gStyle->SetPadRightMargin(0.0);
 	gStyle->SetPadLeftMargin(0.16);
    	gStyle->SetPadRightMargin(0.04);//0.02

        TCanvas *c_xs = new TCanvas("c_xs", "c_xs", 1400, 900);
	c_xs->Divide(1,1);
	c_xs->cd(1);

        float ymax=1500;
        float xmax=550;
        float xmin=0;

        TH2D *f2d_xs=new TH2D("f2d_xs","",450,xmin,xmax,ymax,0,ymax);
        f2d_xs->GetXaxis()->SetTitle("Proton Kinetic Energy [MeV]");
        f2d_xs->GetYaxis()->SetTitleOffset(1.3);
        f2d_xs->GetYaxis()->SetTitle("P-Ar inelastic cross section [mb]");
        f2d_xs->Draw("");

        //gr_truexs->SetLineWidth(2);
        //gr_truexs->SetMarkerColor(3);
        //gr_truexs->SetLineColor(3);
        //gr_truexs->SetMarkerStyle(21);
        //gr_truexs->SetMarkerSize(1.5);

        total_inel_KE->SetLineColor(4);
        total_inel_KE->Draw("c same");

	reco_xs_all->SetMarkerColor(2);
	reco_xs_all->SetLineColor(2);
	reco_xs_all->SetMarkerStyle(20);
	reco_xs_all->Draw("p same");
	reco_xs->Draw("p same");



/*
	xs_GENIE_ha2018->Draw("c same");
	xs_NEUT_2019->Draw("c same");
	xs_GENIE_INCL_pp->Draw("c same");
	xs_GENIE_hN2018->Draw("c same");
	xs_NuWRO_2019->Draw("c same");

        //gr_truexs->Draw("p same");

        gr_recoxs->SetLineWidth(2);
        gr_recoxs->SetMarkerColor(1);
        gr_recoxs->SetLineColor(1);
        gr_recoxs->SetMarkerStyle(20);
        gr_recoxs->SetMarkerSize(1.5);
        gr_recoxs->Draw("p same");

	double ke_ff_mean_stop=4.12972e+02;
	double ke_ff_sigma_stop=4.08988e+01;
	TLine l_mean(ke_ff_mean_stop,0,ke_ff_mean_stop,ymax);
	TLine l_left(ke_ff_mean_stop-3.*ke_ff_sigma_stop,0,ke_ff_mean_stop-3.*ke_ff_sigma_stop,ymax);
	TLine l_right(ke_ff_mean_stop+3.*ke_ff_sigma_stop,0,ke_ff_mean_stop+3.*ke_ff_sigma_stop,ymax);
	l_mean.SetLineColor(1);
	l_left.SetLineColor(1);
        l_right.SetLineColor(1);

	l_left.SetLineStyle(2);	
	l_right.SetLineStyle(2);	

	//l_mean.Draw("same");
	//l_left.Draw("same");
	//l_right.Draw("same");
*/

        //TLegend *leg_xs = new TLegend(0.6,0.6,0.9,0.85);
        //TLegend *leg_xs = new TLegend(0.4,0.6,0.9,0.85);
        TLegend *leg_xs = new TLegend(0.2,0.7,0.9,0.85);
        leg_xs->SetFillStyle(0);
	leg_xs->SetNColumns(2);
        leg_xs->AddEntry(reco_xs, "Data(stat.)", "pe");
        leg_xs->AddEntry(reco_xs_all, "Data(stat.+sys.)", "pe");
        //leg_xs->AddEntry(gr_truexs, "MC Truth", "pe");
        leg_xs->AddEntry(total_inel_KE, "Geant4", "l");

/*
        leg_xs->AddEntry(xs_GENIE_ha2018, "GENIE hA2018", "l");
        leg_xs->AddEntry(xs_GENIE_hN2018, "GENIE hN2018", "l");
        leg_xs->AddEntry(xs_NEUT_2019, "NEUT 2019", "l");
        leg_xs->AddEntry(xs_NuWRO_2019, "NuWRO 2019", "l");
        leg_xs->AddEntry(xs_GENIE_INCL_pp, "GENIE INCL++", "l");
*/
	
        //leg_xs->AddEntry(gr_recoxs, "MC Reco", "pe");
        leg_xs->Draw();

   	//pDUNE Logo
        float logo_y=ymax+15;
        TLatex **txt_pdune1=new TLatex*[1];
        txt_pdune1[0]=new TLatex(xmin, logo_y, Form("#bf{DUNE:ProtoDUNE-SP}"));
        txt_pdune1[0]->SetTextColor(1);
        txt_pdune1[0]->Draw();
        //
        //Beam Logo
        TLatex **txt_p1=new TLatex*[1];
        //txt_p1[0]=new TLatex(xmax-6.3, logo_y, Form("Protons (1 GeV/c)")); //x:0-20
        //txt_p1[0]=new TLatex(xmax-120, logo_y, Form("Protons (1 GeV/c)")); //x:0-40
        //txt_p1[0]=new TLatex(142, logo_y, Form("Protons (1 GeV/c)")); //x:0-40
        txt_p1[0]=new TLatex(390, logo_y, Form("Protons (1 GeV/c)")); //x:0-40

        txt_p1[0]->SetTextColor(1);
        txt_p1[0]->SetTextSize(0.05);
        txt_p1[0]->Draw();

	c_xs->Print(Form("%sxs_sys_data.eps",str_fig_out.Data()));
	//c_xs->Print(Form("%sxs_data_noreco.eps",outpath.Data()));

/*
        TH2D *f2d_xs_ext=new TH2D("f2d_xs_ext","",450,xmin,600,ymax,0,ymax);
        f2d_xs_ext->GetXaxis()->SetTitle("Proton Kinetic Energy [MeV]");
        f2d_xs_ext->GetYaxis()->SetTitleOffset(1.3);
        f2d_xs_ext->GetYaxis()->SetTitle("P-Ar inelastic cross section [mb]");
	f2d_xs_ext->Draw();
        total_inel_KE->Draw("c same");
        gr_truexs->Draw("p same");
        gr_recoxs->Draw("p same");
        leg_xs->Draw();
        txt_pdune1[0]->Draw();
        TLatex **txt2_p1=new TLatex*[1];
        txt2_p1[0]=new TLatex(600-180, logo_y, Form("Protons (1 GeV/c)")); //x:0-40
        txt2_p1[0]->Draw();
	c_xs->Print(Form("%sxs_data_largerrange.eps",outpath.Data()));
*/


}
