#include "THStack.h"

void plot_reco_cosineTheta(TString fin, TString fout_path) {

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	//gStyle->SetOptStat(0);

	TFile *f0 = TFile::Open(fin.Data());

	TH1D *cosine_inel=(TH1D *)f0->Get(Form("reco_cosineTheta_inel"));
	TH1D *cosine_el=(TH1D *)f0->Get(Form("reco_cosineTheta_el"));

	TH1D *cosine_midcosmic=(TH1D *)f0->Get(Form("reco_cosineTheta_midcosmic"));
	TH1D *cosine_midpi=(TH1D *)f0->Get(Form("reco_cosineTheta_midpi"));
	TH1D *cosine_midp=(TH1D *)f0->Get(Form("reco_cosineTheta_midp"));

	TH1D *cosine_midmu=(TH1D *)f0->Get(Form("reco_cosineTheta_midmu"));
	TH1D *cosine_mideg=(TH1D *)f0->Get(Form("reco_cosineTheta_mideg"));
	TH1D *cosine_midother=(TH1D *)f0->Get(Form("reco_cosineTheta_midother"));

	cosine_inel->SetFillColor(2); cosine_inel->SetLineColor(2);
	cosine_el->SetFillColor(4); cosine_el->SetLineColor(4);

	cosine_midp->SetFillColor(3); cosine_midp->SetLineColor(3);
	cosine_midcosmic->SetFillColor(5); cosine_midcosmic->SetLineColor(5);
	cosine_midpi->SetFillColor(6); cosine_midpi->SetLineColor(6);

	cosine_midmu->SetFillColor(28); cosine_midmu->SetLineColor(28);
	cosine_mideg->SetFillColor(30); cosine_mideg->SetLineColor(30);
	cosine_midother->SetFillColor(15); cosine_midother->SetLineColor(15);

	cosine_el->GetXaxis()->SetTitle("cos#Theta");

	THStack* hs=new THStack("hs","");
	hs->Add(cosine_inel);
	hs->Add(cosine_el);

	hs->Add(cosine_midcosmic);
	hs->Add(cosine_midp);
	hs->Add(cosine_midpi);

	hs->Add(cosine_midmu);
	hs->Add(cosine_mideg);
	hs->Add(cosine_midother);


	TCanvas *c_=new TCanvas(Form("c"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_->Divide(1,1);
	c_->cd(1);
	c_->cd(1)->SetLogy();
	TH2D *f2d=new TH2D("f2d",Form("%s","CaloSize"),100,0.9,1,92000,0,92000);
	f2d->GetXaxis()->SetTitle("cos#Theta");
	f2d->Draw();
	hs->Draw("hist same");

	TLegend *leg = new TLegend(0.2,0.6,0.8,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(cosine_el, "El","f");
	leg->AddEntry(cosine_inel, "Inel","f");

	leg->AddEntry(cosine_midcosmic,"misID:cosmic","f");
	leg->AddEntry(cosine_midp, "misID:p","f");
	leg->AddEntry(cosine_midpi, "misID:#pi","f");

	leg->AddEntry(cosine_midmu,"misID:#mu","f");
	leg->AddEntry(cosine_mideg, "misID:e/#gamma","f");
	leg->AddEntry(cosine_midother, "misID:other","f");

	leg->SetNColumns(3);
	leg->Draw();

	c_->Print(Form("%s/cosinereco_%s.eps",fout_path.Data(),"calosz"));

}
