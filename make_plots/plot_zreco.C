#include "THStack.h"

void plot_zreco(TString fin, TString fout_path, TString str_cut) {

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	gStyle->SetOptStat(0);

	TFile *f0 = TFile::Open(fin.Data());

	TH1D *z_inel=(TH1D *)f0->Get(Form("zend_reco_inel_%s",str_cut.Data()));
	TH1D *z_el=(TH1D *)f0->Get(Form("zend_reco_el_%s",str_cut.Data()));
	TH1D *z_mcs=(TH1D *)f0->Get(Form("zend_reco_mcs_%s",str_cut.Data()));

	TH1D *z_midcosmic=(TH1D *)f0->Get(Form("zend_reco_midcosmic_%s",str_cut.Data()));
	TH1D *z_midpi=(TH1D *)f0->Get(Form("zend_reco_midpi_%s",str_cut.Data()));
	TH1D *z_midp=(TH1D *)f0->Get(Form("zend_reco_midp_%s",str_cut.Data()));

	TH1D *z_midmu=(TH1D *)f0->Get(Form("zend_reco_midmu_%s",str_cut.Data()));
	TH1D *z_mideg=(TH1D *)f0->Get(Form("zend_reco_mideg_%s",str_cut.Data()));
	TH1D *z_midother=(TH1D *)f0->Get(Form("zend_reco_midother_%s",str_cut.Data()));

	z_inel->SetFillColor(2); z_inel->SetLineColor(2);
	z_el->SetFillColor(4); z_el->SetLineColor(4);
	z_mcs->SetFillColor(3); z_mcs->SetLineColor(3);

	z_midp->SetFillColor(7); z_midp->SetLineColor(7);
	z_midcosmic->SetFillColor(5); z_midcosmic->SetLineColor(5);
	z_midpi->SetFillColor(6); z_midpi->SetLineColor(6);

	z_midmu->SetFillColor(1); z_midmu->SetLineColor(1);
	z_mideg->SetFillColor(30); z_mideg->SetLineColor(30);
	z_midother->SetFillColor(15); z_midother->SetLineColor(15);

	z_mcs->GetXaxis()->SetTitle("Reco EndZ [cm]");

	THStack* hs=new THStack("hs","");
	hs->Add(z_mcs);
	hs->Add(z_el);
	hs->Add(z_inel);

	hs->Add(z_midcosmic);
	hs->Add(z_midp);
	hs->Add(z_midpi);

	hs->Add(z_midmu);
	hs->Add(z_mideg);
	hs->Add(z_midother);


	TCanvas *c_=new TCanvas(Form("c"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_->Divide(1,1);
	c_->cd(1);
	//TH2D *f2d=new TH2D("f2d","",200,-50,150,1200,0,1200);
	TH2D *f2d=new TH2D("f2d","",200,-50,150,3200,0,3200);
	f2d->GetXaxis()->SetTitle("Reco EndZ [cm]");
	//f2d->GetYaxis()->SetTitle("Counts");
	f2d->Draw();
	//c_->cd(1)->SetLogy();

	hs->Draw("same");

	TLegend *leg = new TLegend(0.2,0.6,0.8,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(z_mcs,"Nohs","f");
	leg->AddEntry(z_el, "El","f");
	leg->AddEntry(z_inel, "Inel","f");

	leg->AddEntry(z_midcosmic,"misID:cosmic","f");
	leg->AddEntry(z_midp, "misID:p","f");
	leg->AddEntry(z_midpi, "misID:#pi","f");

	leg->AddEntry(z_midmu,"misID:#mu","f");
	leg->AddEntry(z_mideg, "misID:e/#gamma","f");
	leg->AddEntry(z_midother, "misID:other","f");

	leg->SetNColumns(3);
	leg->Draw();

	c_->Print(Form("%s/zreco_%s.eps",fout_path.Data(),str_cut.Data()));

}
