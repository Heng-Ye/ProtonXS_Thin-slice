#include "TVector3.h"
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
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

#include "../headers/BasicParameters.h"
#include "../headers/BasicFunctions.h"


void plot_Advanced_KEff() {

	TString outpath="./plot_Advanced_KEff/";
        TString fmc="../mc_ke_advancedKEff.root";

	TString str_end=Form(""); //all protons

        //plot style --------------------------------------------------------------------//
        gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
        gROOT->SetStyle("protoDUNEStyle");
        gROOT->ForceStyle();
        gStyle->SetTitleX(0.5);
        gStyle->SetTitleAlign(23);
        gStyle->SetOptStat(0);

	//gStyle->SetPalette(53);

	//read mc [after bmrw] ------------------------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fmc.Data());

	TH1D *h1d_keff_truth_el=(TH1D *)f_mc->Get("h1d_keff_truth_el");
	TH1D *h1d_keff_truth_inel=(TH1D *)f_mc->Get("h1d_keff_truth_inel");

	TH1D *h1d_keff_reco_constErange_el=(TH1D *)f_mc->Get("h1d_keff_reco_constErange_el");
	TH1D *h1d_keff_reco_constErange_inel=(TH1D *)f_mc->Get("h1d_keff_reco_constErange_inel");

	TH1D *h1d_keff_reco_constEcalo_el=(TH1D *)f_mc->Get("h1d_keff_reco_constEcalo_el");
	TH1D *h1d_keff_reco_constEcalo_inel=(TH1D *)f_mc->Get("h1d_keff_reco_constEcalo_inel");

	TH1D *h1d_keff_reco_Edept_range_el=(TH1D *)f_mc->Get("h1d_keff_reco_Edept_range_el");
	TH1D *h1d_keff_reco_Edept_calo_el=(TH1D *)f_mc->Get("h1d_keff_reco_Edept_calo_el");

	TH1D *h1d_keff_reco_Edept_range_inel=(TH1D *)f_mc->Get("h1d_keff_reco_Edept_range_inel");
	//TH1D *h1d_keff_reco_Edept_calo_inel=(TH1D *)f_mc->Get("h1d_keff_reco_Edept_calo_inel");





	TCanvas *c_el=new TCanvas(Form("c_el"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	//gStyle->SetPadRightMargin(0.16); 
	//gStyle->SetPadLeftMargin(0.12); 

	c_el->Divide(1,1);
	c_el->cd(1)->SetLogy();


	TH2D *f2d_el=new TH2D("f2d_el","",1000,0,800,100,0.5,50000);
	f2d_el->SetTitle("Elastic-scattering Protons; Proton KE_{ff} [MeV]; Counts");
	f2d_el->Draw();
	
	h1d_keff_truth_el->SetLineColor(4);
	h1d_keff_reco_constErange_el->SetLineColor(1);
	h1d_keff_reco_Edept_range_el->SetLineColor(2);

	h1d_keff_truth_el->Draw("same");
	h1d_keff_reco_constErange_el->Draw("same");
	//h1d_keff_reco_Edept_calo_el
	h1d_keff_reco_Edept_range_el->Draw("same");

	TLegend *leg = new TLegend(0.15,0.7,0.5,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(h1d_keff_truth_el, "Truth","l");
	leg->AddEntry(h1d_keff_reco_constErange_el, "Reco: Const E-loss (using Erange)","l");
	leg->AddEntry(h1d_keff_reco_Edept_range_el, "Reco: CSDA Estimation","l");
	leg->Draw();



	TCanvas *c_inel=new TCanvas(Form("c_inel"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	//gStyle->SetPadRightMargin(0.16); 
	//gStyle->SetPadLeftMargin(0.12); 

	c_inel->Divide(1,1);
	c_inel->cd(1)->SetLogy();


	TH2D *f2d_inel=new TH2D("f2d_el","",1000,0,800,100,0.5,50000);
	f2d_inel->SetTitle("Inelastic-scattering Protons; Proton KE_{ff} [MeV]; Counts");
	f2d_inel->Draw();
	
	h1d_keff_truth_inel->SetLineColor(4);
	h1d_keff_reco_constErange_inel->SetLineColor(1);
	h1d_keff_reco_Edept_range_inel->SetLineColor(2);
	//h1d_keff_reco_Edept_calo_inel->SetLineColor(3);

	h1d_keff_truth_inel->Draw("same");
	h1d_keff_reco_constErange_inel->Draw("same");
	//h1d_keff_reco_Edept_calo_inel->Draw("same");
	h1d_keff_reco_Edept_range_inel->Draw("same");

	TLegend *leg_inel = new TLegend(0.15,0.7,0.5,0.9);
	leg_inel->SetFillStyle(0);
	leg_inel->AddEntry(h1d_keff_truth_inel, "Truth","l");
	leg_inel->AddEntry(h1d_keff_reco_constErange_inel, "Reco: Const E-loss (using Erange)","l");
	leg_inel->AddEntry(h1d_keff_reco_Edept_range_inel, "Reco: Adapted E-loss Estimation","l");
	//leg_inel->AddEntry(h1d_keff_reco_Edept_calo_inel, "Reco: Adapted E-loss Estimation2","l");

	leg_inel->Draw();




}
