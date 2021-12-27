#include "THStack.h"

void plot_trklen_different_config() {

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	//gStyle->SetOptStat(0);
	
	TString str_cut=Form("Pos");

	
	//MC SCE with SCE Corr
	TFile *fmc_sce_withCorr = TFile::Open("../mc_proton_after_bmrw_HD.root");
	TH1D *trklen_sce_withCorr=(TH1D *)fmc_sce_withCorr->Get(Form("h1d_trklen_%s",str_cut.Data()));
	TH1D *trklen_true=(TH1D *)fmc_sce_withCorr->Get(Form("h1d_truetrklen_%s",str_cut.Data()));

	//MC SCE No SCE Corr
	TFile *fmc_sce_NoCorr = TFile::Open("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1_noscecorr/mc_sceon_nocorr.root");
	TH1D *trklen_sce_noCorr=(TH1D *)fmc_sce_NoCorr->Get(Form("h1d_trklen_%s",str_cut.Data()));

	//MC SCE OFF
	TFile *fmc_sce_off = TFile::Open("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/PDSPProd4a_MC_6GeV_gen_datadriven_reco1_sce_off_v1/mc_sceoff.root");
	TH1D *trklen_sce_off=(TH1D *)fmc_sce_off->Get(Form("h1d_trklen_%s",str_cut.Data()));
	
	trklen_sce_withCorr->SetLineColor(1);
	trklen_sce_noCorr->SetLineColor(2);
	trklen_sce_off->SetLineColor(4);
	trklen_true->SetLineColor(3);

	int n_sce_withCorr=trklen_sce_withCorr->Integral();
	int n_sce_noCorr=trklen_sce_noCorr->Integral();
	int n_sce_off=trklen_sce_off->Integral();
	int n_true=trklen_true->Integral();

	trklen_sce_noCorr->Scale((double)n_sce_withCorr/(double)n_sce_noCorr);
	trklen_sce_off->Scale((double)n_sce_withCorr/(double)n_sce_off);
	trklen_true->Scale((double)n_sce_withCorr/(double)n_true);



	TCanvas *c_=new TCanvas(Form("c"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_->Divide(1,1);
	c_->cd(1);

	TH2D *f2d=new TH2D("f2d",Form("%s","MC Pos"),121,-1,120,2600,0,2600);
	f2d->GetXaxis()->SetTitle("Track Length [cm]");
	f2d->Draw();

	trklen_sce_withCorr->Draw("hist same");
	trklen_sce_noCorr->Draw("hist same");
	//trklen_sce_off->Draw("hist same");
	//trklen_true->Draw("hist same");

	TLegend *leg = new TLegend(0.2,0.6,0.8,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(trklen_sce_withCorr, "SCE ON with Correction","l");
	leg->AddEntry(trklen_sce_noCorr, "SCE ON with No Correction","l");
	//leg->AddEntry(trklen_sce_off, "SCE OFF","l");
	//leg->AddEntry(trklen_true, "True","l");

	//leg->SetNColumns(3);
	leg->Draw();

	//c_->Print(Form("./plots/mc_trklen_different_configs_%s_zoom.eps",str_cut.Data()));
	c_->Print(Form("./plots/mc_trklen_different_configs_%s_sceoncorr_ononcorr.eps",str_cut.Data()));
	//c_->Print(Form("./plots/mc_trklen_different_configs_%s_sceoncorr_sceoff.eps",str_cut.Data()));
	//c_->Print(Form("./plots/mc_trklen_different_configs_%s_sceoncorr_true.eps",str_cut.Data()));




/*
	//TH2D *f2d=new TH2D("f2d",Form("%s",str_cut.Data()),136,-4,132,250,0,250);
	TH2D *f2d=new TH2D("f2d",Form("%s",str_cut.Data()),136,-4,132,300,0,300);
	f2d->GetXaxis()->SetTitle("Track Length [cm]");
	f2d->Draw();
	//c_->cd(1)->SetLogy();
	hs->Draw("hist same");
	trklen_data->Draw("ep same");

	TLegend *leg = new TLegend(0.2,0.6,0.8,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(trklen_data, "Data", "ep");
	leg->AddEntry(trklen_inel, "Inel","f");
	leg->AddEntry(trklen_el, "El","f");

	leg->AddEntry(trklen_midcosmic,"misID:cosmic","f");
	leg->AddEntry(trklen_midp, "misID:p","f");
	leg->AddEntry(trklen_midpi, "misID:#pi","f");

	leg->AddEntry(trklen_midmu,"misID:#mu","f");
	leg->AddEntry(trklen_mideg, "misID:e/#gamma","f");
	leg->AddEntry(trklen_midother, "misID:other","f");

	leg->SetNColumns(3);
	leg->Draw();

	c_->Print(Form("%s/trklentrue_%s.eps",fout_path.Data(),str_cut.Data()));
*/


}
