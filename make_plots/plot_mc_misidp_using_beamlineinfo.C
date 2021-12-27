//#include "../headers/BasicAnaFunc.h"
#include <vector>
#include <iostream>
#include "../headers/SliceParams.h"
#include "../headers/TemplateFitter.h"

void plot_mc_misidp_using_beamlineinfo() {

	TString file_mc=Form("../prod4a_misidpstudy_dxycut_dx4cm_25slcs.root");
	TFile *fin_mc = TFile::Open(file_mc.Data());

	TString file_data=Form("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_run5387_prod4reco2_bkgstudy_dx4cm_25slcs_largerbin.root");
	TFile *fin_data = TFile::Open(file_data.Data());


	TH1D *h1d_dxy_misidp_lenle0=(TH1D*)fin_mc->Get(Form("h1d_dxy_BQ_misidp_lenle0"));
	TH1D *h1d_dxy_misidp_lengt0=(TH1D*)fin_mc->Get(Form("h1d_dxy_BQ_misidp_lengt0"));
	TH1D *h1d_dxy_inel=(TH1D*)fin_mc->Get(Form("h1d_dxy_BQ_inel"));
	TH1D *h1d_dxy_el=(TH1D*)fin_mc->Get(Form("h1d_dxy_BQ_el"));

	TH1D *h1d_dxy_data=(TH1D*)fin_data->Get(Form("h1d_dxy_BQ"));


	h1d_dxy_misidp_lengt0->SetLineColor(3);
	h1d_dxy_misidp_lenle0->SetLineColor(6);
	h1d_dxy_inel->SetLineColor(2);
	h1d_dxy_el->SetLineColor(4);

	//config -------------------------------------------------------//
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	//gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleYOffset(1.2);
	gStyle->SetTitleXOffset(2.4);
	gStyle->SetTitleAlign(23); 
	gStyle->SetTitleX(0.5);
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
	//get data -------------------------------------------------------------------------------------//

	TCanvas *c_ = new TCanvas("c_", "c_", 1000, 1200);
	c_->Divide(1,3);
	c_->cd(1);
	c_->cd(1)->SetLogy();
	//gStyle->SetTitleX(0.99);
	//gStyle->SetOptStat(0);
	h1d_dxy_misidp_lengt0->GetXaxis()->SetTitleOffset(1.45);
	h1d_dxy_misidp_lengt0->SetTitle("BQ; #sqrt{(X_{beam}-X_{TPC})^{2}+(Y_{beam}-Y_{TPC})^{2}} [cm];Counts");
	h1d_dxy_misidp_lengt0->Draw("hist");
	h1d_dxy_misidp_lenle0->Draw("hist same");


	TLegend *leg = new TLegend(0.7,0.65,0.85,0.84);
	leg->SetFillStyle(0);
	leg->AddEntry(h1d_dxy_misidp_lengt0, "MisID:P_lenGT0","l");
	leg->AddEntry(h1d_dxy_misidp_lenle0, "MisID:P_lenLE0","l");
	leg->Draw();

	c_->cd(2);
	c_->cd(2)->SetLogy();
	h1d_dxy_inel->SetTitle("; #sqrt{(X_{beam}-X_{TPC})^{2}+(Y_{beam}-Y_{TPC})^{2}} [cm]; Counts");
	h1d_dxy_inel->GetXaxis()->SetTitleOffset(1.45);
	h1d_dxy_inel->Draw("hist same");
	h1d_dxy_el->Draw("hist same");

	TLegend *leg2 = new TLegend(0.7,0.65,0.85,0.84);
	leg2->SetFillStyle(0);
	leg2->AddEntry(h1d_dxy_inel, "Inel.","l");
	leg2->AddEntry(h1d_dxy_el, "El.","l");
	leg2->Draw();

	//c_->Print("./plots_bkgstudy/h1d_misidp_cut.eps");


	//TCanvas *c_d = new TCanvas("c_d", "c_d", 1200, 400);
	//c_d->Divide(1,1);
	//c_d->cd(1);
	//c_d->cd(1)->SetLogy();
	c_->cd(3)->SetLogy();

	h1d_dxy_data->GetXaxis()->SetTitleOffset(1.45);
	h1d_dxy_data->SetTitle("BQ; #sqrt{(X_{beam}-X_{TPC})^{2}+(Y_{beam}-Y_{TPC})^{2}} [cm];Counts");
	h1d_dxy_data->SetLineColor(2);
	h1d_dxy_data->Draw("hist");

	TLegend *leg3 = new TLegend(0.7,0.65,0.85,0.84);
	leg3->SetFillStyle(0);
	leg3->AddEntry(h1d_dxy_data, "Data","l");
	leg3->Draw();
	//c_d->Print("./plots_bkgstudy/h1d_data.eps");
	c_->Print("./plots_bkgstudy/h1d_dxy_all.eps");


/*
	TFile *fin = TFile::Open(file.Data());

        //int nn_cos=10;
	//const int nnn_cos=10;
        //float dcos=0.1;

        int nn_cos=25;
	const int nnn_cos=25;
        float dcos=0.04;

	TProfile2D *tp2d_recorange_eff_el[nnn_cos];
	TProfile2D *tp2d_recorange_eff_inel[nnn_cos];
	TProfile2D *tp2d_recorange_eff_misidp[nnn_cos];

	TProfile2D *tp2d_ntrklen_chi2_el[nnn_cos];
	TProfile2D *tp2d_ntrklen_chi2_inel[nnn_cos];
	TProfile2D *tp2d_ntrklen_chi2_misidp[nnn_cos];
	for (int j=0; j<nn_cos; ++j) {
		tp2d_recorange_eff_el[j]=(TProfile2D *)fin->Get(Form("tp2d_recorange_eff_el_%d",j));
		tp2d_recorange_eff_inel[j]=(TProfile2D *)fin->Get(Form("tp2d_recorange_eff_inel_%d",j));
		tp2d_recorange_eff_misidp[j]=(TProfile2D *)fin->Get(Form("tp2d_recorange_eff_misidp_%d",j));

		tp2d_ntrklen_chi2_el[j]=(TProfile2D *)fin->Get(Form("tp2d_ntrklen_chi2_el_%d",j));
		tp2d_ntrklen_chi2_inel[j]=(TProfile2D *)fin->Get(Form("tp2d_ntrklen_chi2_inel_%d",j));
		tp2d_ntrklen_chi2_misidp[j]=(TProfile2D *)fin->Get(Form("tp2d_ntrklen_chi2_misidp_%d",j));
	}

	//cosTheta - slice-by-slice ----------------------------//
	gStyle->SetPadBottomMargin(0.1);
	//gStyle->SetPadLeftMargin(0.15);
	//gStyle->SetPadRightMargin(0);
	//gStyle->SetPadTopMargin(-0.16);

	TLine *line=new TLine(0,1,120,1);
	line->SetLineStyle(2);

	//inel
	TCanvas *c_inel = new TCanvas("c_inel", "c_inel", 1200, 900);
	c_inel->Divide(5,5);

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
