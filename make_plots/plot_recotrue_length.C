//#include "../headers/BasicAnaFunc.h"
#include <vector>
#include <iostream>
#include "../headers/SliceParams.h"
#include "../headers/TemplateFitter.h"

void plot_recotrue_length() {

	TString file=Form("../prod4a_misidpstudy_dx4cm_25slcs.root");

	TString str=Form("recotrklen_truetrklen");
	TString fout_path=Form("plots_bkgstudy/range_reco_true");
	//TString fout_path2=Form("plots_bkgstudy/ntrklen_chi2pid");

	//TString str=Form("ntrklen_chi2");
	//TString fout_path=Form("plots_bkgstudy/ntrklen_chi2");

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
	gStyle->SetPadRightMargin(0.15);

	//get data -------------------------------------------------------------------------------------//
	TFile *fin = TFile::Open(file.Data());

	TProfile2D *h2d_inel=(TProfile2D *)fin->Get(Form("h2d_%s_inel",str.Data()));
	TProfile2D *h2d_el=(TProfile2D *)fin->Get(Form("h2d_%s_el",str.Data()));
	TProfile2D *h2d_misidp=(TProfile2D *)fin->Get(Form("h2d_%s_misidp",str.Data()));


        int nn_cos=11;
        float dcos=0.1;

	TProfile2D *tp2d_range_reco_true_inel[nn_cos];
	TProfile2D *tp2d_range_reco_true_el[nn_cos];
	TProfile2D *tp2d_range_reco_true_misidp[nn_cos];

	TProfile2D *tp2d_ntrklen_chi2_inel[nn_cos];
	TProfile2D *tp2d_ntrklen_chi2_el[nn_cos];
	TProfile2D *tp2d_ntrklen_chi2_misidp[nn_cos];

	for (int j=0; j<nn_cos; ++j) {
		tp2d_range_reco_true_inel[j]=(TProfile2D *)fin->Get(Form("tp2d_range_reco_true_inel_%d",j));
		tp2d_range_reco_true_el[j]=(TProfile2D *)fin->Get(Form("tp2d_range_reco_true_el_%d",j));
		tp2d_range_reco_true_misidp[j]=(TProfile2D *)fin->Get(Form("tp2d_range_reco_true_misidp_%d",j));

		tp2d_ntrklen_chi2_inel[j]=(TProfile2D *)fin->Get(Form("tp2d_ntrklen_chi2_inel_%d",j));
		tp2d_ntrklen_chi2_el[j]=(TProfile2D *)fin->Get(Form("tp2d_ntrklen_chi2_el_%d",j));
		tp2d_ntrklen_chi2_misidp[j]=(TProfile2D *)fin->Get(Form("tp2d_ntrklen_chi2_misidp_%d",j));
	}






	//draw fig
	//inel
	TCanvas *c_inel = new TCanvas("c_inel", "c_inel", 900,700);
	h2d_inel->SetTitle("Inel.");
	//c_inel->SetLogz();
	h2d_inel->GetZaxis()->SetTitle("cos#Theta");
	h2d_inel->Draw("colz");
	TLine *line=new TLine(0,0,120,120);
	line->SetLineStyle(2);
	line->SetLineColor(1);
	line->Draw();
	c_inel->Print(Form("%s_inel.eps",fout_path.Data()));

	//el
	TCanvas *c_el = new TCanvas("c_el", "c_el", 900,700);
	h2d_el->GetZaxis()->SetTitle("cos#Theta");
	h2d_el->Draw("colz");
	line->Draw();
	c_el->Print(Form("%s_el.eps",fout_path.Data()));

	//misid:p
	TCanvas *c_misidp = new TCanvas("c_misidp", "c_misidp", 900,700);
	//c_misidp->SetLogz();
	h2d_misidp->GetZaxis()->SetTitle("cos#Theta");
	h2d_misidp->Draw("colz");
	line->Draw();
	c_misidp->Print(Form("%s_misidp.eps",fout_path.Data()));


/*
	//all-in-one
	TCanvas *c_all = new TCanvas("c_all", "c_all", 900,700);
	h2d_inel->SetMarkerColor(2);
	h2d_inel->SetMarkerStyle(20);
	h2d_inel->SetMarkerSize(.3);

	h2d_el->SetMarkerColor(4);
	//h2d_el->SetMarkerStyle(24);
	h2d_el->SetMarkerSize(.3);
	
	h2d_misidp->SetMarkerColor(3);
	h2d_misidp->SetMarkerStyle(4);
	h2d_misidp->SetMarkerSize(.4);
	
	c_all->cd(1);
	h2d_inel->Draw("p");
	h2d_el->Draw("psame");
	h2d_misidp->Draw("psame");

	c_all->Print(Form("%s_all.png",fout_path.Data()));
	c_all->Print(Form("%s_all.eps",fout_path.Data()));
*/

	//cosTheta - slice-by-slice ----------------------------//
	for (int j=0; j<nn_cos; ++j) {
        	float tmp_min=(float)j*dcos;
                float tmp_max=tmp_min+dcos;

		TCanvas *c_all_ = new TCanvas("c_all_", "c_all_", 2000,1200);
		c_all_->Divide(4,2);
		//c_all_->SetTitle(Form("cos#Theta:%.1f-%.1f",tmp_min,tmp_max));
		TPaveLabel* title = new TPaveLabel(0.35,0.96,0.65,0.99,Form("cos#Theta:%.1f-%.1f",tmp_min,tmp_max));
		title->SetLineColor(0);
		title->SetTextColor(4);
		title->SetFillColor(0);
		title->SetBorderSize(0);
  		title->Draw();

		tp2d_range_reco_true_el[j]->SetMarkerColor(4);
		tp2d_range_reco_true_inel[j]->SetMarkerColor(2);
		tp2d_range_reco_true_misidp[j]->SetMarkerColor(3);

		tp2d_range_reco_true_el[j]->SetMarkerSize(.2);
		tp2d_range_reco_true_inel[j]->SetMarkerSize(.2);
		tp2d_range_reco_true_misidp[j]->SetMarkerSize(.2);

		tp2d_range_reco_true_el[j]->SetTitle("El.");
		tp2d_range_reco_true_inel[j]->SetTitle("Inel.");
		tp2d_range_reco_true_misidp[j]->SetTitle("MisID:P");
		
		c_all_->cd(1);
		tp2d_range_reco_true_el[j]->Draw();
		line->Draw();

		c_all_->cd(2);
		tp2d_range_reco_true_inel[j]->Draw();
		line->Draw();

		c_all_->cd(3);
		tp2d_range_reco_true_misidp[j]->Draw();
		line->Draw();

		c_all_->cd(3);
		tp2d_range_reco_true_misidp[j]->Draw();
		line->Draw();

		c_all_->cd(4);
		tp2d_range_reco_true_inel[j]->Draw();
		tp2d_range_reco_true_el[j]->Draw("same");
		tp2d_range_reco_true_misidp[j]->Draw("same");
		line->Draw();

		//ntrklen vs chi2pid
		tp2d_ntrklen_chi2_el[j]->SetMarkerColor(4);
		tp2d_ntrklen_chi2_inel[j]->SetMarkerColor(2);
		tp2d_ntrklen_chi2_misidp[j]->SetMarkerColor(3);

		tp2d_ntrklen_chi2_el[j]->SetMarkerSize(.2);
		tp2d_ntrklen_chi2_inel[j]->SetMarkerSize(.2);
		tp2d_ntrklen_chi2_misidp[j]->SetMarkerSize(.2);

		tp2d_ntrklen_chi2_el[j]->SetTitle("El.");
		tp2d_ntrklen_chi2_inel[j]->SetTitle("Inel.");
		tp2d_ntrklen_chi2_misidp[j]->SetTitle("MisID:P");
		
		c_all_->cd(5);
		tp2d_range_reco_true_el[j]->Draw();

		c_all_->cd(6);
		tp2d_ntrklen_chi2_inel[j]->Draw();

		c_all_->cd(7);
		tp2d_ntrklen_chi2_misidp[j]->Draw();

		c_all_->cd(8);
		tp2d_ntrklen_chi2_inel[j]->Draw();
		tp2d_ntrklen_chi2_el[j]->Draw("same");
		tp2d_ntrklen_chi2_misidp[j]->Draw("same");








		TString fig_out_all;
		if (j==0) fig_out_all=Form("%s_cosTheta.pdf(",fout_path.Data());
		if (j>0&&j<nn_cos-1) fig_out_all=Form("%s_cosTheta.pdf",fout_path.Data());
		if (j==nn_cos-1) fig_out_all=Form("%s_cosTheta.pdf)",fout_path.Data());

		c_all_->Print(fig_out_all.Data());

		delete c_all_;
	}







}
