//#include "../headers/BasicAnaFunc.h"
#include <vector>
#include <iostream>
#include "../headers/SliceParams.h"
#include "../headers/TemplateFitter.h"

void plot_recotrue_length() {

	TString file=Form("../prod4a_misidpstudy_dx4cm_25slcs.root");

	//TString str=Form("recotrklen_truetrklen");
	//TString fout_path=Form("plots_bkgstudy/range_reco_true_");

	TString str=Form("ntrklen_chi2");
	TString fout_path=Form("plots_bkgstudy/ntrklen_chi2");

	//config
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	//gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	gStyle->SetOptStat(0);
	//gStyle->SetPadLeftMargin(0.13);

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
	gStyle->SetPalette(kBlackBody);


	//get data -------------------------------------------------------------------------------------//
	TFile *fin = TFile::Open(file.Data());

	TProfile2D *h2d_inel=(TProfile2D *)fin->Get(Form("h2d_%s_inel",str.Data()));
	TProfile2D *h2d_el=(TProfile2D *)fin->Get(Form("h2d_%s_el",str.Data()));
	TProfile2D *h2d_misidp=(TProfile2D *)fin->Get(Form("h2d_%s_misidp",str.Data()));

	//draw fig
	//inel
	TCanvas *c_inel = new TCanvas("c_inel", "c_inel", 900,700);
	h2d_inel->SetTitle("Inel.");
	//c_inel->SetLogz();
	h2d_inel->Draw("colz");
	h2d_inel->GetZaxis()->SetTitle("cos#Theta");
	TLine *line=new TLine(0,0,120,120);
	line->SetLineStyle(2);
	line->SetLineColor(1);
	//line->Draw();
	c_inel->Print(Form("%s_inel.eps",fout_path.Data()));

	//el
	TCanvas *c_el = new TCanvas("c_el", "c_el", 900,700);
	h2d_el->Draw("colz");
	//line->Draw();
	c_el->Print(Form("%s_el.eps",fout_path.Data()));

	//misid:p
	TCanvas *c_misidp = new TCanvas("c_misidp", "c_misidp", 900,700);
	//c_misidp->SetLogz();

	h2d_misidp->Draw("colz");
	//line->Draw();
	c_misidp->Print(Form("%s_misidp.eps",fout_path.Data()));







}
