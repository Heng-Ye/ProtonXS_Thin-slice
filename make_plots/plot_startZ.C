//#include "../headers/BasicAnaFunc.h"
#include <vector>
#include <iostream>
#include "TGraphAsymmErrors.h"
#include "../headers/SliceParams.h"
#include "../headers/TemplateFitter.h"

void plot_startZ() {

	TString file_mc1=Form("../prod4a_thinslice_dx4cm_25slcs_1stHitKEff.root");
	TFile *fmc1 = TFile::Open(file_mc1.Data());
	TH1D *zst_mc=(TH1D *)fmc1->Get(Form("reco_startZ_sce"));
	zst_mc->SetLineColor(4);

	TString file_mc2=Form("../mc_truelen_fixtruelencalc.root");
	TFile *fmc2 = TFile::Open(file_mc2.Data());
	TH2D *h2d_mc=(TH2D *)fmc2->Get(Form("h2d_trueEndZ_ketrue_NoCut_inel"));
	TH1D *zEnd_true=h2d_mc->ProjectionX();
	zEnd_true->SetLineColor(6); 	

	TString file_data=Form("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_run5387_prod4a_thinslice_dx4cm_25slcs.root");
	TFile *fdata = TFile::Open(file_data.Data());
	TH1D *zst_data=(TH1D *)fdata->Get(Form("reco_startZ_sce"));
	zst_data->SetLineColor(2);
	
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

	TLine *line=new TLine(0,1,120,1);
	line->SetLineStyle(2);

	//inel
	TCanvas *c_ = new TCanvas("c_", "c_", 1200, 900);
	c_->Divide(1,2);
	c_->cd(1);
	zst_mc->Scale((double)zst_data->Integral()/(double)zst_mc->Integral());
	zst_mc->SetTitle("CaloSz Cut; Reco. StartZ [cm]; Counts");
	zst_mc->GetXaxis()->SetRangeUser(-4,9);
	zst_mc->Draw("hist");
	zst_data->Draw("hist same");

        TLegend *leg0 = new TLegend(0.5,0.65,.85,0.88);
        leg0->SetFillStyle(0);
        leg0->AddEntry(zst_data, "Data (after SCE corr.)", "l");
        leg0->AddEntry(zst_mc, "MC (after SCE corr.)", "l");
        leg0->Draw();

	c_->cd(2);
	zEnd_true->SetTitle("Inelastic Scattering Protons; True StartZ [cm]; Counts");
	zEnd_true->GetXaxis()->SetRangeUser(-30,50);
	zEnd_true->Draw("hist");


	TString fig_out=Form("startZ_afterSCE.eps");
	c_->Print(fig_out.Data());







}
