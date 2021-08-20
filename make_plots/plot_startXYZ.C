#include "THStack.h"

void plot_startXYZ(TString fin_mc, TString fin_data, TString fout_path) {

	//pDUNE style
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	//gStyle->SetOptStat(0);

	//mc -------------------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fin_mc.Data());
	TH1D *reco_startX_sce_mc=(TH1D *)f_mc->Get(Form("reco_startX_sce"));
	TH1D *reco_startY_sce_mc=(TH1D *)f_mc->Get(Form("reco_startY_sce"));
	TH1D *reco_startZ_sce_mc=(TH1D *)f_mc->Get(Form("reco_startZ_sce"));

	reco_startX_sce_mc->SetLineColor(2); reco_startX_sce_mc->SetMarkerColor(2); 
	reco_startY_sce_mc->SetLineColor(2); reco_startY_sce_mc->SetMarkerColor(2); 
	reco_startZ_sce_mc->SetLineColor(2); reco_startZ_sce_mc->SetMarkerColor(2); 

	//data -------------------------------------------------------------------------//
	TFile *f_data = TFile::Open(fin_data.Data());
	TH1D *reco_startX_sce_data=(TH1D *)f_data->Get(Form("reco_startX_sce"));
	TH1D *reco_startY_sce_data=(TH1D *)f_data->Get(Form("reco_startY_sce"));
	TH1D *reco_startZ_sce_data=(TH1D *)f_data->Get(Form("reco_startZ_sce"));

	reco_startX_sce_data->SetLineColor(1); reco_startX_sce_data->SetMarkerColor(1); 
	reco_startY_sce_data->SetLineColor(1); reco_startY_sce_data->SetMarkerColor(1); 
	reco_startZ_sce_data->SetLineColor(1); reco_startZ_sce_data->SetMarkerColor(1); 

	int n_data=reco_startZ_sce_data->Integral();
	int n_mc=reco_startZ_sce_mc->Integral();

	reco_startX_sce_mc->Scale((double)n_data/(double)n_mc);
	reco_startY_sce_mc->Scale((double)n_data/(double)n_mc);
	reco_startZ_sce_mc->Scale((double)n_data/(double)n_mc);


	TCanvas *c_x=new TCanvas(Form("cx"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_x->Divide(1,1);
	c_x->cd(1);
	reco_startX_sce_mc->SetTitle("Calo Size; Start X (after SCE correction) [cm];");
	reco_startX_sce_mc->Draw("hist");
	reco_startX_sce_data->Draw("ep same");	

	TLegend *legx = new TLegend(0.1,0.6,0.6,0.85);
	legx->SetFillStyle(0);
	legx->AddEntry(reco_startX_sce_data, "Data","ep");
	legx->AddEntry(reco_startX_sce_mc, "MC","l");
	legx->Draw();
	c_x->Print(Form("%s/startX_sce.eps",fout_path.Data()));




	TCanvas *c_y=new TCanvas(Form("cy"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_y->Divide(1,1);
	c_y->cd(1);
	reco_startY_sce_mc->SetTitle("Calo Size; Start Y (after SCE correction) [cm];");
	reco_startY_sce_mc->Draw("hist");
	reco_startY_sce_data->Draw("ep same");	

	TLegend *legy = new TLegend(0.1,0.6,0.6,0.85);
	legy->SetFillStyle(0);
	legy->AddEntry(reco_startX_sce_data, "Data","ep");
	legy->AddEntry(reco_startX_sce_mc, "MC","l");
	legy->Draw();
	c_y->Print(Form("%s/startY_sce.eps",fout_path.Data()));



	TCanvas *c_z=new TCanvas(Form("cz"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_z->Divide(1,1);
	c_z->cd(1);
	reco_startZ_sce_mc->SetTitle("Calo Size; Start Z (after SCE correction) [cm];");
	reco_startZ_sce_mc->Draw("hist");
	reco_startZ_sce_data->Draw("ep same");	

	TLegend *legz = new TLegend(0.75,0.75,0.85,0.85);
	legz->SetFillStyle(0);
	legz->AddEntry(reco_startZ_sce_data, "Data","ep");
	legz->AddEntry(reco_startZ_sce_mc, "MC","l");
	legz->Draw();
	c_z->Print(Form("%s/startZ_sce.eps",fout_path.Data()));






/*
	TH2D *f2d=new TH2D("f2d",Form("%s",str_obs.Data()),100, -10, 10,11200,0,11200);
	f2d->GetXaxis()->SetTitle(Form("#Delta%s/#sigma_{%s}",str_obs.Data(), str_obs.Data()));

   	TPad *pad1 = new TPad(Form("pad1"), Form("pad1"), 0, 0., 1, 1.);
    	pad1->SetGridx();
    	pad1->SetGridy();
    	pad1->Draw();
    	pad1->cd();

	f2d->GetXaxis()->SetNdivisions(32);
   	//f2d->GetYaxis()->SetNdivisions(64);
	f2d->Draw("same");
	hs->Draw("hist same");

	TLegend *leg = new TLegend(0.2,0.6,0.8,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(h1d_el, "El","f");
	leg->AddEntry(h1d_inel, "Inel","f");

	leg->AddEntry(h1d_midcosmic,"misID:cosmic","f");
	leg->AddEntry(h1d_midp, "misID:p","f");
	leg->AddEntry(h1d_midpi, "misID:#pi","f");

	leg->AddEntry(h1d_midmu,"misID:#mu","f");
	leg->AddEntry(h1d_mideg, "misID:e/#gamma","f");
	leg->AddEntry(h1d_midother, "misID:other","f");

	leg->SetNColumns(3);
	leg->Draw();

	c_->Update();
	c_->Print(Form("%s/hdelta%s.eps",fout_path.Data(),str_obs.Data()));
*/

}
