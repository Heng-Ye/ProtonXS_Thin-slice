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


void plot_Advanced_KEend() {

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

	TH1D *h1d_keend_truth_el=(TH1D *)f_mc->Get("h1d_keend_truth_el");
	TH1D *h1d_keend_truth_inel=(TH1D *)f_mc->Get("h1d_keend_truth_inel");

	TH1D *h1d_kebb_truth_el=(TH1D *)f_mc->Get("h1d_kebb_truth_el");
	TH1D *h1d_kebb_truth_inel=(TH1D *)f_mc->Get("h1d_kebb_truth_inel");

	TH1D *h1d_kebb_reco_constErange_el=(TH1D *)f_mc->Get("h1d_kebb_reco_constErange_el");
	TH1D *h1d_kebb_reco_constErange_inel=(TH1D *)f_mc->Get("h1d_kebb_reco_constErange_inel");

	TH1D *h1d_kebb_reco_Edept_range_el=(TH1D *)f_mc->Get("h1d_kebb_reco_Edept_range_el");
	TH1D *h1d_kebb_reco_Edept_range_inel=(TH1D *)f_mc->Get("h1d_kebb_reco_Edept_range_inel");

	TH1D *h1d_kecalo_reco_constEcalo_el=(TH1D *)f_mc->Get("h1d_kecalo_reco_constEcalo_el");
	TH1D *h1d_kecalo_reco_constEcalo_inel=(TH1D *)f_mc->Get("h1d_kecalo_reco_constEcalo_inel");

	TH1D *h1d_kecalo_reco_Edept_calo_el=(TH1D *)f_mc->Get("h1d_kecalo_reco_Edept_calo_el");
	TH1D *h1d_kecalo_reco_Edept_calo_inel=(TH1D *)f_mc->Get("h1d_kecalo_reco_Edept_calo_inel");

	h1d_keend_truth_el->SetLineColor(3);
	h1d_keend_truth_inel->SetLineColor(3);

	h1d_kebb_truth_el->SetLineColor(3);
	h1d_kebb_truth_inel->SetLineColor(3);

	h1d_kebb_reco_Edept_range_el->SetLineColor(4);
	h1d_kebb_reco_Edept_range_inel->SetLineColor(4);

	h1d_kebb_reco_constErange_el->SetLineColor(4);
	h1d_kebb_reco_constErange_inel->SetLineColor(4);

	h1d_kecalo_reco_constEcalo_el->SetLineColor(2);
	h1d_kecalo_reco_constEcalo_inel->SetLineColor(2);

	h1d_kecalo_reco_Edept_calo_el->SetLineColor(2);
	h1d_kecalo_reco_Edept_calo_inel->SetLineColor(2);


	

	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.08); 
	gStyle->SetPadLeftMargin(0.15); 

	TCanvas *c_el1=new TCanvas(Form("c_el1"),"",900, 600);
	c_el1->Divide(1,1);
	//c_el1->cd(1)->SetLogy();
	
	TH2D *f2d_el1=new TH2D("f2d_el1","",900,-400,500,100,0.,3000);
	f2d_el1->SetTitle("Elastic-scattering Protons; Proton Kinetic Energy [MeV]; Counts");
	f2d_el1->GetYaxis()->SetTitleOffset(1.2);
	f2d_el1->Draw();
	h1d_keend_truth_el->Draw("same");
	//h1d_kebb_truth_el->Draw("same");
	h1d_kebb_reco_constErange_el->Draw("same");
	h1d_kecalo_reco_constEcalo_el->Draw("same");
	//h1d_kebb_reco_Edept_range_el->Draw("same");
	//h1d_kecalo_reco_Edept_calo_el->Draw("same");

	TLegend *leg_el1 = new TLegend(0.493318,0.642361,0.883073,0.875);
	leg_el1->SetFillStyle(0);
	leg_el1->AddEntry(h1d_keend_truth_el, "Truth ","l");
	leg_el1->AddEntry(h1d_kebb_reco_constErange_el, "Reco (Bethe-Bloch: Const. E-loss)","l");
	leg_el1->AddEntry(h1d_kecalo_reco_constEcalo_el, "Reco (Calo: Const. E-loss)","l");
	leg_el1->Draw();
	c_el1->Print(Form("%sKEend_el_1.eps",outpath.Data()));

	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.08); 
	gStyle->SetPadLeftMargin(0.15); 

	TCanvas *c_el2=new TCanvas(Form("c_el2"),"",900, 600);
	//gStyle->SetTitleX(0.5); 
	//gStyle->SetTitleAlign(23); 
	//gStyle->SetPadRightMargin(0.16); 
	//gStyle->SetPadLeftMargin(0.15); 

	c_el2->Divide(1,1);
	//c_el2->cd(1)->SetLogy();
	TH2D *f2d_el2=new TH2D("f2d_el2","",900,-400,500,100,0.,3000);
	f2d_el2->SetTitle("Elastic-scattering Protons; Proton Kinetic Energy [MeV]; Counts");
	f2d_el2->GetYaxis()->SetTitleOffset(1.2);
	f2d_el2->Draw();
	h1d_keend_truth_el->Draw("same");
	h1d_kebb_reco_Edept_range_el->Draw("same");
	h1d_kecalo_reco_Edept_calo_el->Draw("same");

	TLegend *leg_el2 = new TLegend(0.553318,0.642361,0.883073,0.875);
	leg_el2->SetFillStyle(0);
	leg_el2->AddEntry(h1d_keend_truth_el, "Truth ","l");
	//leg_el->AddEntry(h1d_kebb_truth_el, "Truth (Bethe-Bloch: KE_{ff}(truth), Range(truth))","l");
	leg_el2->AddEntry(h1d_kebb_reco_Edept_range_el, "Reco (Bethe-Bloch: Edept.)","l");
	leg_el2->AddEntry(h1d_kecalo_reco_Edept_calo_el, "Reco (Calo: Edept.)","l");
	leg_el2->Draw();
	c_el2->Print(Form("%sKEend_el_2.eps",outpath.Data()));



	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.08); 
	gStyle->SetPadLeftMargin(0.15); 

	TCanvas *c_el3=new TCanvas(Form("c_el3"),"",900, 600);
	c_el3->Divide(1,1);
	c_el3->cd(1)->SetLogy();
	TH2D *f2d_el3=new TH2D("f2d_el3","",900,-100,100,100,0.,40000);
	f2d_el3->SetTitle("Elastic-scattering Protons; Proton Kinetic Energy [MeV]; Counts");
	f2d_el3->GetYaxis()->SetTitleOffset(1.2);
	f2d_el3->Draw();
	h1d_keend_truth_el->Draw("same");
	h1d_kebb_reco_Edept_range_el->Draw("same");
	h1d_kecalo_reco_Edept_calo_el->Draw("same");

	TLegend *leg_el3 = new TLegend(0.169333,0.64881,0.5,0.880882);
	leg_el3->SetFillStyle(0);
	leg_el3->AddEntry(h1d_keend_truth_el, "Truth ","l");
	//leg_el->AddEntry(h1d_kebb_truth_el, "Truth (Bethe-Bloch: KE_{ff}(truth), Range(truth))","l");
	leg_el3->AddEntry(h1d_kebb_reco_Edept_range_el, "Reco (Bethe-Bloch: Edept.)","l");
	leg_el3->AddEntry(h1d_kecalo_reco_Edept_calo_el, "Reco (Calo: Edept.)","l");
	leg_el3->Draw();
	c_el3->Print(Form("%sKEend_el_3_logy.eps",outpath.Data()));


	//inel
	TCanvas *c_inel1=new TCanvas(Form("c_inel1"),"",900, 600);
	c_inel1->Divide(1,1);
	//c_inel1->cd(1)->SetLogy();
	
	TH2D *f2d_inel1=new TH2D("f2d_inel1","",1000,-400,600,100,0.,1000);
	f2d_inel1->SetTitle("Inelastic-scattering Protons; Proton Kinetic Energy [MeV]; Counts");
	f2d_inel1->GetYaxis()->SetTitleOffset(1.2);
	f2d_inel1->Draw();
	h1d_keend_truth_inel->Draw("same");
	h1d_kecalo_reco_constEcalo_inel->Draw("same");
	h1d_kebb_reco_constErange_inel->Draw("same");

	TLegend *leg_inel1 = new TLegend(0.2,0.6,0.6,0.88);
	leg_inel1->SetFillStyle(0);
	leg_inel1->AddEntry(h1d_keend_truth_inel, "Truth ","l");
	leg_inel1->AddEntry(h1d_kebb_reco_constErange_inel, "Reco (Bethe-Bloch: Const. E-loss)","l");
	leg_inel1->AddEntry(h1d_kecalo_reco_constEcalo_inel, "Reco (Calo: Const. E-loss)","l");
	leg_inel1->Draw();
	c_inel1->Print(Form("%sKEend_inel_1.eps",outpath.Data()));



	TCanvas *c_inel2=new TCanvas(Form("c_inel2"),"",900, 600);
	//gStyle->SetTitleX(0.5); 
	//gStyle->SetTitleAlign(23); 
	//gStyle->SetPadRightMargin(0.16); 
	//gStyle->SetPadLeftMargin(0.15); 

	c_inel2->Divide(1,1);
	//c_inel2->cd(1)->SetLogy();
	TH2D *f2d_inel2=new TH2D("f2d_inel2","",1000,-400,600,100,0.,1000);

	f2d_inel2->SetTitle("Inelastic-scattering Protons; Proton Kinetic Energy [MeV]; Counts");
	f2d_inel2->GetYaxis()->SetTitleOffset(1.2);
	f2d_inel2->Draw();
	h1d_keend_truth_inel->Draw("same");
	h1d_kebb_reco_Edept_range_inel->Draw("same");
	h1d_kecalo_reco_Edept_calo_inel->Draw("same");

	TLegend *leg_inel2 = new TLegend(0.2,0.6,0.6,0.88);
	leg_inel2->SetFillStyle(0);
	leg_inel2->AddEntry(h1d_keend_truth_inel, "Truth ","l");
	//leg_el->AddEntry(h1d_kebb_truth_el, "Truth (Bethe-Bloch: KE_{ff}(truth), Range(truth))","l");
	leg_inel2->AddEntry(h1d_kebb_reco_Edept_range_inel, "Reco (Bethe-Bloch: Edept.)","l");
	leg_inel2->AddEntry(h1d_kecalo_reco_Edept_calo_inel, "Reco (Calo: Edept.)","l");
	leg_inel2->Draw();
	c_inel2->Print(Form("%sKEend_inel_2.eps",outpath.Data()));









/*
	//inel
	TCanvas *c_inel=new TCanvas(Form("c_inel"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	//gStyle->SetPadRightMargin(0.16); 
	//gStyle->SetPadLeftMargin(0.12); 

	c_inel->Divide(1,1);
	//c_inel->cd(1)->SetLogy();
	
	TH2D *f2d_inel=new TH2D("f2d_inel","",1000,-100,800,100,0.5,1300);
	f2d_inel->SetTitle("Inelastic-scattering Protons; Proton KE_{ff} [MeV]; Counts");
	f2d_inel->Draw();
	h1d_keend_truth_inel->Draw("same");
	h1d_kebb_truth_inel->Draw("same");
	h1d_kebb_reco_constErange_inel->Draw("same");
	h1d_kebb_reco_Edept_range_inel->Draw("same");
	h1d_kecalo_reco_constEcalo_inel->Draw("same");
	h1d_kecalo_reco_Edept_calo_inel->Draw("same");

	TLegend *leg_inel = new TLegend(0.15,0.7,0.5,0.9);
	leg_inel->SetFillStyle(0);
	leg_inel->AddEntry(h1d_keend_truth_inel, "Truth (MC Trajs.)","l");
	leg_inel->AddEntry(h1d_kebb_truth_inel, "Truth (Bethe-Bloch: KE_{ff}(truth), Range(truth))","l");
	leg_inel->AddEntry(h1d_kebb_reco_constErange_inel, "Reco (Bethe-Bloch: Const. E-loss)","l");
	leg_inel->AddEntry(h1d_kebb_reco_Edept_range_inel, "Reco (Bethe-Bloch: Edept.)","l");
	leg_inel->AddEntry(h1d_kecalo_reco_constEcalo_inel, "Reco (Calo: Const. E-loss)","l");
	leg_inel->AddEntry(h1d_kecalo_reco_Edept_calo_inel, "Reco (Calo: Const. Edept.)","l");
	leg_inel->Draw();
*/



}
