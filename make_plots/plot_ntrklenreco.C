#include "THStack.h"

void plot_ntrklenreco(TString fin, TString fout_path, TString str_cut) {

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	//gStyle->SetOptStat(0);

	TFile *f0 = TFile::Open(fin.Data());

	TH1D *ntrklen_inel=(TH1D *)f0->Get(Form("ntrklen_inel_%s",str_cut.Data()));
	TH1D *ntrklen_el=(TH1D *)f0->Get(Form("ntrklen_el_%s",str_cut.Data()));

	TH1D *ntrklen_midcosmic=(TH1D *)f0->Get(Form("ntrklen_midcosmic_%s",str_cut.Data()));
	TH1D *ntrklen_midpi=(TH1D *)f0->Get(Form("ntrklen_midpi_%s",str_cut.Data()));
	TH1D *ntrklen_midp=(TH1D *)f0->Get(Form("ntrklen_midp_%s",str_cut.Data()));

	TH1D *ntrklen_midmu=(TH1D *)f0->Get(Form("ntrklen_midmu_%s",str_cut.Data()));
	TH1D *ntrklen_mideg=(TH1D *)f0->Get(Form("ntrklen_mideg_%s",str_cut.Data()));
	TH1D *ntrklen_midother=(TH1D *)f0->Get(Form("ntrklen_midother_%s",str_cut.Data()));

	ntrklen_inel->SetFillColor(2); ntrklen_inel->SetLineColor(2);
	ntrklen_el->SetFillColor(4); ntrklen_el->SetLineColor(4);

	ntrklen_midp->SetFillColor(3); ntrklen_midp->SetLineColor(3);
	ntrklen_midcosmic->SetFillColor(5); ntrklen_midcosmic->SetLineColor(5);
	ntrklen_midpi->SetFillColor(6); ntrklen_midpi->SetLineColor(6);

	ntrklen_midmu->SetFillColor(28); ntrklen_midmu->SetLineColor(28);
	ntrklen_mideg->SetFillColor(30); ntrklen_mideg->SetLineColor(30);
	ntrklen_midother->SetFillColor(15); ntrklen_midother->SetLineColor(15);

	ntrklen_el->GetXaxis()->SetTitle("Reco Track Length/CSDA [a.u.]");

	THStack* hs=new THStack("hs","");
	hs->Add(ntrklen_el);
	hs->Add(ntrklen_inel);

	hs->Add(ntrklen_midcosmic);
	hs->Add(ntrklen_midp);
	hs->Add(ntrklen_midpi);

	hs->Add(ntrklen_midmu);
	hs->Add(ntrklen_mideg);
	hs->Add(ntrklen_midother);


	TCanvas *c_=new TCanvas(Form("c"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_->Divide(1,1);
	c_->cd(1);
	TH2D *f2d=new TH2D("f2d",Form("%s",str_cut.Data()),100, -0.02, 1.2, 6200,0,6200);
	f2d->GetXaxis()->SetTitle("Reco Track Length/CSDA [a.u.]");
	f2d->Draw();
	hs->Draw("hist same");

	TLegend *leg = new TLegend(0.2,0.6,0.8,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(ntrklen_el, "El","f");
	leg->AddEntry(ntrklen_inel, "Inel","f");

	leg->AddEntry(ntrklen_midcosmic,"misID:cosmic","f");
	leg->AddEntry(ntrklen_midp, "misID:p","f");
	leg->AddEntry(ntrklen_midpi, "misID:#pi","f");

	leg->AddEntry(ntrklen_midmu,"misID:#mu","f");
	leg->AddEntry(ntrklen_mideg, "misID:e/#gamma","f");
	leg->AddEntry(ntrklen_midother, "misID:other","f");

	leg->SetNColumns(3);
	leg->Draw();

	c_->Print(Form("%s/ntrklenreco_%s.eps",fout_path.Data(),str_cut.Data()));

}
