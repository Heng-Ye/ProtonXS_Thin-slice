#include "THStack.h"

void plot_dz(TString fin, TString fout_path, TString str_cut) {

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	//gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	gStyle->SetOptStat(0);

	//TString str_cut=Form("NoCut");
	TFile *f0 = TFile::Open(fin.Data());

	TH1D *z_inel=(TH1D *)f0->Get(Form("dz_inel_%s",str_cut.Data()));
	TH1D *z_el=(TH1D *)f0->Get(Form("dz_el_%s",str_cut.Data()));
	TH1D *z_mcs=(TH1D *)f0->Get(Form("dz_mcs_%s",str_cut.Data()));

	TH1D *z_upinel=(TH1D *)f0->Get(Form("dz_upinel_%s",str_cut.Data()));
	TH1D *z_upel=(TH1D *)f0->Get(Form("dz_upel_%s",str_cut.Data()));
	TH1D *z_upmcs=(TH1D *)f0->Get(Form("dz_upmcs_%s",str_cut.Data()));

	TH1D *z_misinel=(TH1D *)f0->Get(Form("dz_misinel_%s",str_cut.Data()));
	TH1D *z_misel=(TH1D *)f0->Get(Form("dz_misel_%s",str_cut.Data()));
	TH1D *z_mismcs=(TH1D *)f0->Get(Form("dz_mismcs_%s",str_cut.Data()));

	z_inel->SetFillColor(2); z_inel->SetLineColor(2);
	z_el->SetFillColor(4); z_el->SetLineColor(4);
	z_mcs->SetFillColor(3); z_mcs->SetLineColor(3);

	z_misinel->SetFillColor(7); z_misinel->SetLineColor(7);
	z_misel->SetFillColor(6); z_misel->SetLineColor(6);
	z_mismcs->SetFillColor(5); z_mismcs->SetLineColor(5);

	z_upinel->SetFillColor(41); z_upinel->SetLineColor(41);
	z_upel->SetFillColor(30); z_upel->SetLineColor(30);
	z_upmcs->SetFillColor(48); z_upmcs->SetLineColor(48);

	z_mcs->GetXaxis()->SetTitle("Reco EndZ [cm]");

	THStack* hs=new THStack("hs","");
	hs->Add(z_mcs);
	hs->Add(z_el);
	hs->Add(z_inel);

	hs->Add(z_upmcs);
	hs->Add(z_upel);
	hs->Add(z_upinel);

	hs->Add(z_mismcs);
	hs->Add(z_misel);
	hs->Add(z_misinel);


	TCanvas *c_=new TCanvas(Form("c"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_->Divide(1,1);
	c_->cd(1);
	TH2D *f2d=new TH2D("f2d","",100,-100,100,1200,0,1200);
	f2d->GetXaxis()->SetTitle("EndZ (Reco-Truth) [cm]");
	//f2d->GetYaxis()->SetTitle("Counts");
	f2d->Draw();
	//c_->cd(1)->SetLogy();

	hs->Draw("same");

	TLegend *leg = new TLegend(0.5,0.6,0.8,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(z_mcs,"NoHS","f");
	leg->AddEntry(z_el, "El","f");
	leg->AddEntry(z_inel, "Inel","f");

	leg->AddEntry(z_mismcs,"mis:NoHS","f");
	leg->AddEntry(z_misel, "mis:El","f");
	leg->AddEntry(z_misinel, "mis:Inel","f");

	leg->AddEntry(z_upmcs,"up:NoHS","f");
	leg->AddEntry(z_upel, "up:El","f");
	leg->AddEntry(z_upinel, "up:Inel","f");

	leg->SetNColumns(3);
	leg->Draw();

	c_->Print(Form("%s/dz_%s.eps",fout_path.Data(),str_cut.Data()));

}
