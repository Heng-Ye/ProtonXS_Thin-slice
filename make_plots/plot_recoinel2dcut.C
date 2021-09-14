#include "THStack.h"
#include "../headers/BasicParameters.h"
#include "../headers/BasicFunctions.h"

void plot_recoinel2dcut() {

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);

	TFile *f_mc = TFile::Open("../mc_proton_studyRecoInelCut.root");
	TH2D* ntrklen_chi2pid_BQ=(TH2D*)f_mc->Get("ntrklen_chi2pid_BQ");
	TH2D* ntrklen_chi2pid_RecoInel=(TH2D*)f_mc->Get("ntrklen_chi2pid_RecoInel");
	TH2D* ntrklen_chi2pid_RecoInel2=(TH2D*)f_mc->Get("ntrklen_chi2pid_RecoInel2");
	TH2D* ntrklen_chi2pid_RecoInel3=(TH2D*)f_mc->Get("ntrklen_chi2pid_RecoInel3");
	ntrklen_chi2pid_BQ->SetTitle(";Proton Track Length/CSDA; #chi^{2} PID [a.u.]");
	ntrklen_chi2pid_RecoInel->SetTitle(";Proton Track Length/CSDA; #chi^{2} PID [a.u.]");
	ntrklen_chi2pid_RecoInel2->SetTitle(";Proton Track Length/CSDA; #chi^{2} PID [a.u.]");
	ntrklen_chi2pid_RecoInel3->SetTitle(";Proton Track Length/CSDA; #chi^{2} PID [a.u.]");


	TFile *f_data = TFile::Open("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_proton_bmrw.root");
	TH2D* ntrklen_chi2pid_BQ_data=(TH2D*)f_data->Get("ntrklen_chi2pid_BQ");
	TCanvas *c0=new TCanvas(Form("c0"),"",900, 600);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	c0->Divide(1,1);
	c0->cd(1);
	ntrklen_chi2pid_BQ_data->SetMarkerSize(.7);
	ntrklen_chi2pid_BQ_data->SetTitle(";Proton Track Length/CSDA; #chi^{2} PID [a.u.]");
	ntrklen_chi2pid_BQ_data->Draw("colz");
	c0->Print("./plots_recoinel/ntrklen_chi2_BQ_data_colz.eps");	


	TCanvas *c1=new TCanvas(Form("c1"),"",900, 600);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	c1->Divide(1,1);
	c1->cd(1);
	ntrklen_chi2pid_BQ->SetMarkerSize(.1);
	ntrklen_chi2pid_BQ->Draw("");
	c1->Print("./plots_recoinel/ntrklen_chi2_BQ.eps");	


	TCanvas *c2=new TCanvas(Form("c2"),"",900, 600);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	c2->Divide(1,1);
	c2->cd(1);
	ntrklen_chi2pid_RecoInel->Draw("colz");
	c2->Print("./plots_recoinel/ntrklen_chi2_defaultCut.eps");	

	TCanvas *c3=new TCanvas(Form("c3"),"",900, 600);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	c3->Divide(1,1);
	c3->cd(1);
	ntrklen_chi2pid_RecoInel2->Draw("colz");
	c3->Print("./plots_recoinel/ntrklen_chi2_Cut1.eps");	

	TCanvas *c4=new TCanvas(Form("c4"),"",900, 600);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	c4->Divide(1,1);
	c4->cd(1);
	ntrklen_chi2pid_RecoInel3->Draw("colz");
	c4->Print("./plots_recoinel/ntrklen_chi2_Cut2.eps");	


}
