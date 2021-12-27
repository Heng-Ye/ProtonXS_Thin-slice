//#include "../headers/BasicAnaFunc.h"
#include <vector>
#include <iostream>
#include "TGraphAsymmErrors.h"
#include "../headers/SliceParams.h"
#include "../headers/TemplateFitter.h"

void plot_vtx_misidp_eff() {

	TString file_mc=Form("../mc_vtx_rangereco.root");
	TFile *fmc = TFile::Open(file_mc.Data());

	TH2D *h2d=(TH2D *)fmc->Get(Form("h2d_trklen_dr_RecoInelMIDP"));
	TH2D *h2d_el=(TH2D *)fmc->Get(Form("h2d_trklen_dr_RecoInelMIDP_el"));
	TH2D *h2d_inel=(TH2D *)fmc->Get(Form("h2d_trklen_dr_RecoInelMIDP_inel"));
	TH2D *h2d_misidp=(TH2D *)fmc->Get(Form("h2d_trklen_dr_RecoInelMIDP_misidp"));

	//figout
	TString fig_out="plots_vtx/recolen_dr_recoelmidp.eps";
	TString fig_out_inel="plots_vtx/recolen_dr_recoelmidp_inel.eps";
	TString fig_out_el="plots_vtx/recolen_dr_recoelmidp_el.eps";
	TString fig_out_misidp="plots_vtx/recolen_dr_recoelmidp_misidp.eps";

       	//int n_b=300;
        //double b_min=0;
        //double b_max=150;
       	int n_x=h2d->GetNbinsX();
        double x_min=h2d->GetXaxis()->GetXmin();
        double x_max=h2d->GetXaxis()->GetXmax();
	double bin_x=(x_max-x_min)/(double)n_x;
	cout<<"n_x:"<<n_x<<endl;
	cout<<"x_min:"<<x_min<<" - x_max:"<<x_max<<endl;

        //int n_dr=600;
        //double dr_min=-150;
        //double dr_max=150;
        int n_y=h2d->GetNbinsY();
        double y_min=h2d->GetYaxis()->GetXmin();
        double y_max=h2d->GetYaxis()->GetXmax();
	double bin_y=(y_max-y_min)/(double)n_y;
	cout<<"\nn_y:"<<n_y<<endl;
	cout<<"y_min:"<<y_min<<" - y_max:"<<y_max<<endl;

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

	//inel ---------------------------------------------------------------------------------------//
	TCanvas *c_inel = new TCanvas("c_inel", "c_inel", 1200, 900);
	c_inel->Divide(1,1);
	c_inel->cd(1);
	h2d->SetTitle("Reco InEl+MisID:P-Rich Cut; Reco Track Length [cm]; #Delta L");
	//h2d->GetYaxis()->SetRangeUser(-10,10);
	h2d->Draw("colz");
	c_inel->Print(fig_out.Data());


	h2d_inel->SetTitle("Reco InEl+MisID:P-Rich Cut [True Inel.]; Reco Track Length [cm]; #Delta L");
	//h2d_inel->GetYaxis()->SetRangeUser(-10,10);
	h2d_inel->Draw("colz");
	c_inel->Print(fig_out_inel.Data());

	h2d_el->SetTitle("Reco InEl+MisID:P-Rich Cut [True El.]; Reco Track Length [cm]; #Delta L");
	//h2d_el->GetYaxis()->SetRangeUser(-10,10);
	h2d_el->Draw("colz");
	c_inel->Print(fig_out_el.Data());

	h2d_misidp->SetTitle("Reco InEl+MisID:P-Rich Cut [MisID:P]; Reco Track Length [cm]; #Delta L");
	//h2d_misidp->GetYaxis()->SetRangeUser(-10,10);
	h2d_misidp->Draw("colz");
	c_inel->Print(fig_out_misidp.Data());

	//h2d_inel->GetYaxis()->SetRangeUser(-10,10);
	//h2d_inel->GetXaxis()->SetRangeUser(0,110);
	//h2d_inel->Draw("colz");
	//c_inel->Print(fig_out_inel_zoom.Data());












}
