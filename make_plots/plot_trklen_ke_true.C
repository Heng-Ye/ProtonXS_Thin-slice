#include "/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4_MC_1GeV_reco1_sce_datadriven/SliceParams.h"


void plot_trklen_ke_true(TString fin, TString outpath){

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	//gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 

	gStyle->SetOptStat(0);

	TFile *f0 = TFile::Open(fin.Data());
	TH2D *trklen_ke_true_inel=(TH2D *)f0->Get("trklen_ke_true_inel");
	TH2D *trklen_ke_true_el=(TH2D *)f0->Get("trklen_ke_true_el");

	TH2D* frame2d=new TH2D("frame2d","", 130, 0, 130, 620, -20, 600); //trklenend_2d
	frame2d->GetXaxis()->SetTitle("True Track Length [cm]");
	frame2d->GetYaxis()->SetTitle("True KE [MeV]");

	TCanvas *c1=new TCanvas("c1",""); 
	c1->Divide(1,1);
	c1->cd(1);
	frame2d->Draw();
	trklen_ke_true_el->Draw("colz same");
	c1->Print(Form("%s/trklen_ke_trueEL.eps",outpath.Data()));


	TCanvas *c2=new TCanvas("c2",""); 
	c2->Divide(1,1);
	c2->cd(1);
	frame2d->Draw();
	trklen_ke_true_inel->Draw("colz same");
	c2->Print(Form("%s/trklen_ke_trueInEL.eps",outpath.Data()));


}

