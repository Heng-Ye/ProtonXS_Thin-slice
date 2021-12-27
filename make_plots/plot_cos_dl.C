//#include "../headers/BasicAnaFunc.h"
#include <vector>
#include <iostream>
#include "TGraphAsymmErrors.h"
#include "../headers/SliceParams.h"
#include "../headers/TemplateFitter.h"

void plot_cos_dl() {

	//TString file_mc=Form("../mc_vtx_rangereco_old.root");
	TString file_mc=Form("../mc_vtx_rangereco.root");
	TFile *fmc = TFile::Open(file_mc.Data());

	TH2D *h2d=(TH2D *)fmc->Get(Form("h2d_cos_dr_PosRecoInel"));
	TH2D *h2d_el=(TH2D *)fmc->Get(Form("h2d_cos_dr_PosRecoInel_el"));
	TH2D *h2d_inel=(TH2D *)fmc->Get(Form("h2d_cos_dr_PosRecoInel_inel"));
	TH2D *h2d_misidp=(TH2D *)fmc->Get(Form("h2d_cos_dr_PosRecoInel_misidp"));

	//figout
	TString fig_out="plots_misidp/cos_dr_posrecoinel.eps";
	TString fig_out_inel="plots_misidp/cos_dr_posrecoinel_inel.eps";
	TString fig_out_el="plots_misidp/cos_dr_posrecoinel_el.eps";
	TString fig_out_misidp="plots_misidp/cos_dr_posrecoinel_misidp.eps";


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

	gStyle->SetPalette(kBird);
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

	TCanvas *c_ = new TCanvas("c_", "c_", 1200, 900);
	c_->Divide(1,1);
	TH2D *f2d=new TH2D("","",10,0,1, 300, -150,150);


	TLine *line1=new TLine(0.96,-150,0.96,150);
	line1->SetLineColor(2);
	line1->SetLineStyle(2);

	TLine *line2=new TLine(0.9,-150,0.9,150);
	line2->SetLineColor(3);
	line2->SetLineStyle(2);

	f2d->SetTitle(Form("Pos+RecoInel Cut;cos#Theta; #Delta L [cm]"));
	f2d->Draw();
        h2d->Draw("colz same");
	line1->Draw("same");
	line2->Draw("same");
	c_->Print(fig_out.Data());

	f2d->SetTitle(Form("Pos+RecoInel Cut [Inel.];cos#Theta; #Delta L [cm]"));
	f2d->Draw();
        h2d_inel->Draw("colz same");
	line1->Draw("same");
	line2->Draw("same");
	c_->Print(fig_out_inel.Data());

	f2d->SetTitle(Form("Pos+RecoInel [El.];cos#Theta; #Delta L [cm]"));
	f2d->Draw();
        h2d_el->Draw("colz same");
	line1->Draw("same");
	line2->Draw("same");
	c_->Print(fig_out_el.Data());

	f2d->SetTitle(Form("Pos+RecoInel [MisID:P];cos#Theta; #Delta L [cm]"));
	f2d->Draw();
        h2d_misidp->Draw("colz same");
	line1->Draw("same");
	line2->Draw("same");
	c_->Print(fig_out_misidp.Data());


}
