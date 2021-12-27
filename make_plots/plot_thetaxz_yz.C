//#include "../headers/BasicAnaFunc.h"
#include <vector>
#include <iostream>
#include "TGraphAsymmErrors.h"
#include "../headers/SliceParams.h"
#include "../headers/TemplateFitter.h"

void Pal2()
{
   static Int_t  colors[50];
   static Bool_t initialized = kFALSE;
 
   Double_t Red[3]    = { 1.00, 0.50, 0.00};
   Double_t Green[3]  = { 0.50, 0.00, 1.00};
   Double_t Blue[3]   = { 1.00, 0.00, 0.50};
   Double_t Length[3] = { 0.00, 0.50, 1.00 };
 
   if(!initialized){
      Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,50);
      for (int i=0; i<50; i++) colors[i] = FI+i;
      initialized = kTRUE;
      return;
   }
   gStyle->SetPalette(50,colors);
}

void plot_thetaxz_yz() {

	TString file=Form("../mc_truelen.root");
	TFile *fin = TFile::Open(file.Data());

	TProfile2D *h2d_thetaXZ_thetaYZ_BQ_inel=(TProfile2D *)fin->Get("h2d_thetaXZ_thetaYZ_BQ_inel");
	TProfile2D *h2d_thetaXZ_thetaYZ_BQ_el=(TProfile2D *)fin->Get("h2d_thetaXZ_thetaYZ_BQ_el");
	TProfile2D *h2d_thetaXZ_thetaYZ_BQ_misidp=(TProfile2D *)fin->Get("h2d_thetaXZ_thetaYZ_BQ_misidp");

	TProfile2D *h2d_thetaXZ_thetaYZ_Pos_inel=(TProfile2D *)fin->Get("h2d_thetaXZ_thetaYZ_Pos_inel");
	TProfile2D *h2d_thetaXZ_thetaYZ_Pos_el=(TProfile2D *)fin->Get("h2d_thetaXZ_thetaYZ_Pos_el");
	TProfile2D *h2d_thetaXZ_thetaYZ_Pos_misidp=(TProfile2D *)fin->Get("h2d_thetaXZ_thetaYZ_Pos_misidp");


	h2d_thetaXZ_thetaYZ_BQ_inel->SetTitle("Inel. [BQ Cut]; #Theta_{XZ} [deg]; #Theta_{YZ} [deg]");
	h2d_thetaXZ_thetaYZ_BQ_el->SetTitle("El. [BQ Cut]; #Theta_{XZ} [deg]; #Theta_{YZ} [deg]");
	h2d_thetaXZ_thetaYZ_BQ_misidp->SetTitle("MisID:P [BQ Cut]; #Theta_{XZ} [deg]; #Theta_{YZ} [deg]");

	h2d_thetaXZ_thetaYZ_Pos_inel->SetTitle("Inel. [Pos Cut]; #Theta_{XZ} [deg]; #Theta_{YZ} [deg]");
	h2d_thetaXZ_thetaYZ_Pos_el->SetTitle("El. [Pos Cut]; #Theta_{XZ} [deg]; #Theta_{YZ} [deg]");
	h2d_thetaXZ_thetaYZ_Pos_misidp->SetTitle("MisID:P [Pos Cut]; #Theta_{XZ} [deg]; #Theta_{YZ} [deg]");


	//config -------------------------------------------------------//
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleYOffset(1.2);
	//gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	gStyle->SetOptStat(0);
	gStyle->SetPadLeftMargin(0.15);
        gStyle->SetFrameLineWidth(1);

	//gStyle->SetPalette(60);
	//gStyle->SetPalette(kLightTemperature);
	//gStyle->SetPalette(kTemperatureMap);
	//gStyle->SetPalette(kColorPrintableOnGrey);
	//gStyle->SetPalette(kBlueGreenYellow);
	//gStyle->SetPalette(kFruitPunch);
	//gStyle->SetPalette(kGreenPink);
	//gStyle->SetPalette(kPastel);
	//gStyle->SetPalette(kRedBlue);
	//gStyle->SetPalette(kCool);
	//gStyle->SetPalette(kWaterMelon);
	//gStyle->SetPalette(kBlackBody);
	//gStyle->SetPadTopMargin(0.2);
	//gStyle->SetPadRightMargin(0.17);
	//config -------------------------------------------------------//


	//inel
	//gStyle->SetPalette(kBlueRedYellow);
	//gStyle->SetPalette(kRainBow);
	gStyle->SetPalette(kBlackBody);
	//gStyle->SetPalette(kCMYK);
	//gStyle->SetPalette(kGreenPink);
	//gStyle->SetPalette(kPastel);
	//gStyle->SetPalette(kTemperatureMap);
	//gStyle->SetPalette(kColorPrintableOnGrey);
   	//TExec *ex2 = new TExec("ex2","Pal2();");
   	//ex2->Draw();

	TCanvas *c_ = new TCanvas("c_", "c_", 1200, 900);
	c_->Divide(1,1);
	//TH2D *f2d=new TH2D("f2d","",15,-1,150,15,-1,150);
	//f2d->SetTitle(";#Theta_{XZ}[deg];#Theta_{YZ}[deg]");
	//f2d->Draw();
	
	//TProfile2D *h2d=(TProfile2D *)h2d_thetaXZ_thetaYZ_Pos_misidp->Clone();
	//TProfile2D *h2d=(TProfile2D *)h2d_thetaXZ_thetaYZ_Pos_inel->Clone();
	//TProfile2D *h2d=(TProfile2D *)h2d_thetaXZ_thetaYZ_Pos_el->Clone();
	
	//TProfile2D *h2d=(TProfile2D *)h2d_thetaXZ_thetaYZ_BQ_misidp->Clone();
	//TProfile2D *h2d=(TProfile2D *)h2d_thetaXZ_thetaYZ_BQ_inel->Clone();
	TProfile2D *h2d=(TProfile2D *)h2d_thetaXZ_thetaYZ_BQ_el->Clone();


	h2d->GetXaxis()->SetRangeUser(-1,150);
	h2d->GetYaxis()->SetRangeUser(-1,150);
	h2d->Draw("colz");

	//c_->Print("./plots_bkgstudy/thetaxz_thetayz_Pos_misidp.eps");
	//c_->Print("./plots_bkgstudy/thetaxz_thetayz_Pos_inel.eps");
	//c_->Print("./plots_bkgstudy/thetaxz_thetayz_Pos_el.eps");

	//c_->Print("./plots_bkgstudy/thetaxz_thetayz_BQ_misidp.eps");
	//c_->Print("./plots_bkgstudy/thetaxz_thetayz_BQ_inel.eps");
	c_->Print("./plots_bkgstudy/thetaxz_thetayz_BQ_el.eps");




/*
	TLine *line=new TLine(0,1,120,1);
	line->SetLineStyle(2);

	//inel
	TCanvas *c_ = new TCanvas("c_", "c_", 1200, 900);
	c_->Divide(1,1);
	c_->SetLogy();
	h1d_truetrklen_NoCut->SetLineColor(15);
	h1d_truetrklen_Pan->SetLineColor(6);
	h1d_truetrklen_CaloSz->SetLineColor(1);
	h1d_truetrklen_Pos->SetLineColor(2);
	h1d_truetrklen_BQ->SetLineColor(3);
	h1d_truetrklen_RecoInel->SetLineColor(4);
	h1d_truetrklen_BeamMatch->SetLineColor(7);

	TH2D *f2d=new TH2D("f2d","",150, 0, 150, 200, 0.1, 20000);
	f2d->GetXaxis()->SetTitle("True Track Length [cm]");
	f2d->GetYaxis()->SetTitle("Counts");
	f2d->Draw();
	h1d_truetrklen_NoCut->Draw("hist same");

	h1d_truetrklen_Pan->Draw("hist same");
	h1d_truetrklen_CaloSz->Draw("hist same");
	h1d_truetrklen_Pos->Draw("hist same");
	h1d_truetrklen_BQ->Draw("hist same");
	h1d_truetrklen_RecoInel->Draw("hist same");
	h1d_truetrklen_BeamMatch->Draw("hist same");

	TLegend *leg = new TLegend(0.6,0.65,0.85,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(h1d_truetrklen_NoCut, "True Inel.","l");
	leg->AddEntry(h1d_truetrklen_Pan, "Pandora Slice","l");
	leg->AddEntry(h1d_truetrklen_CaloSz, "Calo Size","l");
	leg->AddEntry(h1d_truetrklen_Pos, "Position","l");
	leg->AddEntry(h1d_truetrklen_BQ, "BQ (Pos+Cos#Theta)","l");
	leg->AddEntry(h1d_truetrklen_RecoInel, "Reco Inel.","l");
	leg->AddEntry(h1d_truetrklen_BeamMatch, "Beam Match","l");
	leg->Draw();


	TString fig_out=Form("plots_bkgstudy/truelen_cuts.eps");
	c_->Print(fig_out.Data());


	TCanvas *c_1 = new TCanvas("c_1", "c_1", 1200, 900);
	c_1->Divide(1,1);

	TH2D *f2d1=new TH2D("f2d1","",30,0,120,150,0,1.5);
	//TH2D *f2d1=new TH2D("f2d1","",50, 0, 50,150,0,1.5);
	f2d1->GetXaxis()->SetTitle("True Track Length [cm]");
	f2d1->GetYaxis()->SetTitle("Efficiency");
	f2d1->Draw();

	Eff_NoCut->SetLineColor(15);
	Eff_NoCut->SetMarkerColor(15);

	Eff_Pan->SetLineColor(6);
	Eff_Pan->SetMarkerColor(6);

	Eff_CaloSz->SetLineColor(1);
	Eff_CaloSz->SetMarkerColor(1);

	Eff_Pos->SetLineColor(2);
	Eff_Pos->SetMarkerColor(2);

	Eff_BQ->SetLineColor(3);
	Eff_BQ->SetMarkerColor(3);

	Eff_RecoInel->SetLineColor(4);
	Eff_RecoInel->SetMarkerColor(4);

	Eff_BeamMatch->SetLineColor(7);
	Eff_BeamMatch->SetMarkerColor(7);

	Eff_NoCut->SetMarkerStyle(20);
	Eff_Pan->SetMarkerStyle(20);
	Eff_CaloSz->SetMarkerStyle(20);
	Eff_Pos->SetMarkerStyle(20);
	Eff_BQ->SetMarkerStyle(20);
	Eff_RecoInel->SetMarkerStyle(20);
	Eff_BeamMatch->SetMarkerStyle(20);


	Eff_NoCut->Draw("p same");
	Eff_Pan->Draw("p same");
	Eff_CaloSz->Draw("p same");
	Eff_Pos->Draw("p same");
	Eff_BQ->Draw("p same");
	Eff_RecoInel->Draw("p same");
	Eff_BeamMatch->Draw("p same");

	TLegend *leg1 = new TLegend(0.2,0.65,0.85,0.9);
	leg1->SetFillStyle(0);
	leg1->AddEntry(h1d_truetrklen_NoCut, "True Inel.","l");
	leg1->AddEntry(h1d_truetrklen_Pan, "Pandora Slice","l");
	leg1->AddEntry(h1d_truetrklen_CaloSz, "Calo Size","l");
	leg1->AddEntry(h1d_truetrklen_Pos, "Position","l");
	leg1->AddEntry(h1d_truetrklen_BQ, "BQ (Pos+Cos#Theta)","l");
	leg1->AddEntry(h1d_truetrklen_RecoInel, "Reco Inel.","l");
	leg1->AddEntry(h1d_truetrklen_BeamMatch, "Beam Match","l");
	leg1->SetNColumns(3);
	leg1->Draw();

	//TString fig_out1=Form("plots_bkgstudy/eff_truelen_cuts.eps");
	//TString fig_out1=Form("plots_bkgstudy/eff_truelen_cuts_zoom1.eps");
	TString fig_out1=Form("plots_bkgstudy/eff_truelen_cuts_zoom3.eps");
	c_1->Print(fig_out1.Data());

	




	gPad->Update();
	for (int j=0; j<nn_cos; ++j) {
        	float tmp_min=(float)j*dcos;
                float tmp_max=tmp_min+dcos;

		//TPaveLabel* title = new TPaveLabel(0.,0.96,0.15,0.99,Form("cos#Theta:%.1f-%.1f",tmp_min,tmp_max));
		//title->SetLineColor(0);
		//title->SetFillColor(0);
		//title->SetBorderSize(0);
		//title->SetTextColor(4);
  		//title->Draw();

		TString tit_inel=Form("cos#Theta:%.2f-%.2f; Reco Track Length [cm]; Reco Track Length/True Track Length",tmp_min,tmp_max);
		tp2d_recorange_eff_inel[j]->SetTitle(tit_inel.Data());

		c_inel->cd(nn_cos-j);
		tp2d_recorange_eff_inel[j]->GetZaxis()->SetTitle("#chi^{2}PID");
		tp2d_recorange_eff_inel[j]->SetMarkerSize(0.08);
		tp2d_recorange_eff_inel[j]->SetMarkerColor(2);
		tp2d_recorange_eff_inel[j]->Draw("");
		line->Draw();
	}

	TString fig_out_inel=Form("plots_bkgstudy/tp2d_range_reco_eff_inel.eps");
	c_inel->Print(fig_out_inel.Data());

	for (int j=0; j<nn_cos; ++j) {
        	float tmp_min=(float)j*dcos;
                float tmp_max=tmp_min+dcos;

		TString tit_inel=Form("cos#Theta:%.2f-%.2f; Reco Track Length [cm]; Reco Track Length/True Track Length",tmp_min,tmp_max);

		//zoom
		if (j==nn_cos-1) {
			TCanvas *c_inel_zoom = new TCanvas("c_inel_zoom", "c_inel_zoom", 1200, 900);
			c_inel_zoom->Divide(1,1);
			c_inel_zoom->cd(1);
			TH2D *f2d_zoom_inel=new TH2D("f2d_zoom_inel","",20,0,20,40,-.1,2.4);
			f2d_zoom_inel->SetTitle(tit_inel.Data());
			f2d_zoom_inel->Draw();
			//tp2d_recorange_eff_inel[j]->SetMarkerSize(0.5);
			tp2d_recorange_eff_inel[j]->Draw("colz same");

			TLine *linex=new TLine(0,1,20,1);
			linex->SetLineStyle(2);
			linex->Draw();

			TString fig_out_inel_zoom=Form("plots_bkgstudy/tp2d_range_reco_eff_inel_zoom.eps");
			c_inel_zoom->Print(fig_out_inel_zoom.Data());
		}

	}


	//el
	TCanvas *c_el = new TCanvas("c_el", "c_el", 1200, 900);
	c_el->Divide(5,5);

	gPad->Update();
	for (int j=0; j<nn_cos; ++j) {
        	float tmp_min=(float)j*dcos;
                float tmp_max=tmp_min+dcos;

		TString tit_el=Form("cos#Theta:%.2f-%.2f; Reco Track Length [cm]; Reco Track Length/True Track Length",tmp_min,tmp_max);
		tp2d_recorange_eff_el[j]->SetTitle(tit_el.Data());

		c_el->cd(nn_cos-j);
		tp2d_recorange_eff_el[j]->GetZaxis()->SetTitle("#chi^{2}PID");
		tp2d_recorange_eff_el[j]->SetMarkerSize(0.08);
		tp2d_recorange_eff_el[j]->SetMarkerColor(2);
		tp2d_recorange_eff_el[j]->Draw("");
		line->Draw();
	}

	TString fig_out_el=Form("plots_bkgstudy/tp2d_range_reco_eff_el.eps");
	c_el->Print(fig_out_el.Data());

	for (int j=0; j<nn_cos; ++j) {
        	float tmp_min=(float)j*dcos;
                float tmp_max=tmp_min+dcos;

		TString tit_el=Form("cos#Theta:%.2f-%.2f; Reco Track Length [cm]; Reco Track Length/True Track Length",tmp_min,tmp_max);

		//zoom
		if (j==nn_cos-1) {
			TCanvas *c_el_zoom = new TCanvas("c_el_zoom", "c_el_zoom", 1200, 900);
			c_el_zoom->Divide(1,1);
			c_el_zoom->cd(1);
			TH2D *f2d_zoom_el=new TH2D("f2d_zoom_el","",20,0,20,40,-.1,1.2);
			f2d_zoom_el->SetTitle(tit_el.Data());
			f2d_zoom_el->Draw();
			//tp2d_recorange_eff_inel[j]->SetMarkerSize(0.5);
			tp2d_recorange_eff_el[j]->Draw("colz same");

			TLine *linex=new TLine(0,1,20,1);
			linex->SetLineStyle(2);
			linex->Draw();

			TString fig_out_el_zoom=Form("plots_bkgstudy/tp2d_range_reco_eff_el_zoom.eps");
			c_el_zoom->Print(fig_out_el_zoom.Data());
		}

	}


	//misidp
	TCanvas *c_misidp = new TCanvas("c_misidp", "c_misidp", 1200, 900);
	c_misidp->Divide(5,5);

	gPad->Update();
	for (int j=0; j<nn_cos; ++j) {
        	float tmp_min=(float)j*dcos;
                float tmp_max=tmp_min+dcos;

		TString tit_misidp=Form("cos#Theta:%.2f-%.2f; Reco Track Length [cm]; Reco Track Length/True Track Length",tmp_min,tmp_max);
		tp2d_recorange_eff_misidp[j]->SetTitle(tit_misidp.Data());

		c_misidp->cd(nn_cos-j);
		tp2d_recorange_eff_misidp[j]->GetZaxis()->SetTitle("#chi^{2}PID");
		tp2d_recorange_eff_misidp[j]->SetMarkerSize(0.08);
		tp2d_recorange_eff_misidp[j]->SetMarkerColor(2);
		tp2d_recorange_eff_misidp[j]->Draw("");
		line->Draw();
	}

	TString fig_out_misidp=Form("plots_bkgstudy/tp2d_range_reco_eff_misidp.eps");
	c_misidp->Print(fig_out_misidp.Data());

	for (int j=0; j<nn_cos; ++j) {
        	float tmp_min=(float)j*dcos;
                float tmp_max=tmp_min+dcos;

		TString tit_misidp=Form("cos#Theta:%.2f-%.2f; Reco Track Length [cm]; Reco Track Length/True Track Length",tmp_min,tmp_max);

		//zoom
		if (j==nn_cos-1) {
			TCanvas *c_misidp_zoom = new TCanvas("c_misidp_zoom", "c_misidp_zoom", 1200, 900);
			c_misidp_zoom->Divide(1,1);
			c_misidp_zoom->cd(1);
			TH2D *f2d_zoom_misidp=new TH2D("f2d_zoom_misidp","",20,0,20,40,-1,3);
			f2d_zoom_misidp->SetTitle(tit_misidp.Data());
			f2d_zoom_misidp->Draw();
			//tp2d_recorange_eff_inel[j]->SetMarkerSize(0.5);
			tp2d_recorange_eff_misidp[j]->Draw("colz same");

			TLine *linex=new TLine(0,1,20,1);
			linex->SetLineStyle(2);
			linex->Draw();

			TString fig_out_misidp_zoom=Form("plots_bkgstudy/tp2d_range_reco_eff_misidp_zoom.eps");
			c_misidp_zoom->Print(fig_out_misidp_zoom.Data());
		}

	}

	//ntrklen vs chi2pid-------------------------------------------------------------------------------------------------------------------------------//
	//inel
	TCanvas *c1_inel = new TCanvas("c1_inel", "c1_inel", 1200, 900);
	c1_inel->Divide(5,5);

	gPad->Update();
	for (int j=0; j<nn_cos; ++j) {
        	float tmp_min=(float)j*dcos;
                float tmp_max=tmp_min+dcos;

		TString tit_inel=Form("cos#Theta:%.2f-%.2f; Reco Track Length/CSDA; #chi^{2}PID",tmp_min,tmp_max);
		tp2d_ntrklen_chi2_inel[j]->SetTitle(tit_inel.Data());

		c1_inel->cd(nn_cos-j);
		tp2d_ntrklen_chi2_inel[j]->GetZaxis()->SetTitle("cos#Theta");
		tp2d_ntrklen_chi2_inel[j]->SetMarkerSize(0.01);
		tp2d_ntrklen_chi2_inel[j]->SetMarkerColor(2);
		tp2d_ntrklen_chi2_inel[j]->Draw("");
	}

	TString fig1_out_inel=Form("plots_bkgstudy/tp2d_ntrklen_chi2pid_inel.eps");
	c1_inel->Print(fig1_out_inel.Data());


	//el
	TCanvas *c1_el = new TCanvas("c1_el", "c1_el", 1200, 900);
	c1_el->Divide(5,5);

	gPad->Update();
	for (int j=0; j<nn_cos; ++j) {
        	float tmp_min=(float)j*dcos;
                float tmp_max=tmp_min+dcos;

		TString tit_el=Form("cos#Theta:%.2f-%.2f; Reco Track Length/CSDA; #chi^{2}PID",tmp_min,tmp_max);
		tp2d_ntrklen_chi2_el[j]->SetTitle(tit_el.Data());

		c1_el->cd(nn_cos-j);
		tp2d_ntrklen_chi2_el[j]->GetZaxis()->SetTitle("cos#Theta");
		tp2d_ntrklen_chi2_el[j]->SetMarkerSize(0.01);
		tp2d_ntrklen_chi2_el[j]->SetMarkerColor(2);
		tp2d_ntrklen_chi2_el[j]->Draw("");
	}

	TString fig1_out_el=Form("plots_bkgstudy/tp2d_ntrklen_chi2pid_el.eps");
	c1_el->Print(fig1_out_el.Data());

	//misidp
	TCanvas *c1_misidp = new TCanvas("c1_misidp", "c1_misidp", 1200, 900);
	c1_misidp->Divide(5,5);

	gPad->Update();
	for (int j=0; j<nn_cos; ++j) {
        	float tmp_min=(float)j*dcos;
                float tmp_max=tmp_min+dcos;

		TString tit_misidp=Form("cos#Theta:%.2f-%.2f; Reco Track Length/CSDA; #chi^{2}PID",tmp_min,tmp_max);
		tp2d_ntrklen_chi2_misidp[j]->SetTitle(tit_misidp.Data());

		c1_misidp->cd(nn_cos-j);
		tp2d_ntrklen_chi2_misidp[j]->GetZaxis()->SetTitle("cos#Theta");
		tp2d_ntrklen_chi2_misidp[j]->SetMarkerSize(0.01);
		tp2d_ntrklen_chi2_misidp[j]->SetMarkerColor(2);
		tp2d_ntrklen_chi2_misidp[j]->Draw("");
	}

	TString fig1_out_misidp=Form("plots_bkgstudy/tp2d_ntrklen_chi2pid_misidp.eps");
	c1_misidp->Print(fig1_out_misidp.Data());
*/


/*
	for (int j=0; j<nn_cos; ++j) {
        	float tmp_min=(float)j*dcos;
                float tmp_max=tmp_min+dcos;

		TCanvas *c_all_ = new TCanvas("c_all_", "c_all_", 4600,2400);
		c_all_->Divide(3,2);

		gPad->Update();
		//gPad->SetTopMargin(5);

		//TPaveLabel* title = new TPaveLabel(0.35,0.96,0.65,0.99,Form("cos#Theta:%.1f-%.1f",tmp_min,tmp_max));
		TPaveLabel* title = new TPaveLabel(0.,0.96,0.15,0.99,Form("cos#Theta:%.1f-%.1f",tmp_min,tmp_max));
		title->SetLineColor(0);
		title->SetFillColor(0);
		title->SetBorderSize(0);
		title->SetTextColor(4);
  		title->Draw();

		tp2d_recorange_eff_el[j]->SetTitle("El.; Reco Track Length [cm]; Reco Track Length/True Track Length");
		tp2d_recorange_eff_inel[j]->SetTitle("Inel.; Reco Track Length [cm]; Reco Track Length/True Track Length");
		tp2d_recorange_eff_misidp[j]->SetTitle("MisID:P; Reco Track Length [cm]; Reco Track Length/True Track Length");

		c_all_->cd(1);
		tp2d_recorange_eff_inel[j]->GetZaxis()->SetTitle("#chi^{2}PID");
		tp2d_recorange_eff_inel[j]->SetMarkerSize(1);
		tp2d_recorange_eff_inel[j]->Draw("colz");
		line->Draw();

		c_all_->cd(2);
		tp2d_recorange_eff_el[j]->GetZaxis()->SetTitle("#chi^{2}PID");
		tp2d_recorange_eff_el[j]->Draw("colz");
		line->Draw();

		c_all_->cd(3);
		tp2d_recorange_eff_misidp[j]->GetZaxis()->SetTitle("#chi^{2}PID");
		tp2d_recorange_eff_misidp[j]->Draw("colz");
		line->Draw();


		//ntrklen vs chi2PID --------------------------------------------------------------------------------------------//
		tp2d_ntrklen_chi2_el[j]->SetTitle("El.; Reco Track Length/CSDA; #chi^{2}PID");
		tp2d_ntrklen_chi2_inel[j]->SetTitle("Inel. ; Reco Track Length/CSDA; #chi^{2}PID");
		tp2d_ntrklen_chi2_misidp[j]->SetTitle("MisID:P; Reco Track Length/CSDA; #chi^{2}PID");

		c_all_->cd(4);
		tp2d_ntrklen_chi2_inel[j]->GetZaxis()->SetTitle("cos#Theta");
		tp2d_ntrklen_chi2_inel[j]->Draw("colz");

		c_all_->cd(5);
		tp2d_ntrklen_chi2_el[j]->GetZaxis()->SetTitle("cos#Theta");
		tp2d_ntrklen_chi2_el[j]->Draw("colz");

		c_all_->cd(6);
		tp2d_ntrklen_chi2_misidp[j]->GetZaxis()->SetTitle("cos#Theta");
		tp2d_ntrklen_chi2_misidp[j]->Draw("colz");

		TString fig_out_all_eps=Form("%s_cosTheta_%d.eps",fout_path.Data(),j);
		c_all_->Print(fig_out_all_eps.Data());

	
		delete c_all_;
	}
*/






}
