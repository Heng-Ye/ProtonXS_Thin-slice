#include "THStack.h"

void plot_dtrklen(TString fin, TString fout_path, TString str_cut) {

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	gStyle->SetOptStat(0);

	TFile *f0 = TFile::Open(fin.Data());

	TH1D *dtrklen_inel=(TH1D *)f0->Get(Form("dtrklen_inel_%s",str_cut.Data()));
	TH1D *dtrklen_el=(TH1D *)f0->Get(Form("dtrklen_el_%s",str_cut.Data()));

	TH1D *dtrklen_midcosmic=(TH1D *)f0->Get(Form("dtrklen_midcosmic_%s",str_cut.Data()));
	TH1D *dtrklen_midpi=(TH1D *)f0->Get(Form("dtrklen_midpi_%s",str_cut.Data()));
	TH1D *dtrklen_midp=(TH1D *)f0->Get(Form("dtrklen_midp_%s",str_cut.Data()));

	TH1D *dtrklen_midmu=(TH1D *)f0->Get(Form("dtrklen_midmu_%s",str_cut.Data()));
	TH1D *dtrklen_mideg=(TH1D *)f0->Get(Form("dtrklen_mideg_%s",str_cut.Data()));
	TH1D *dtrklen_midother=(TH1D *)f0->Get(Form("dtrklen_midother_%s",str_cut.Data()));

	dtrklen_inel->SetFillColor(2); dtrklen_inel->SetLineColor(2);
	dtrklen_el->SetFillColor(4); dtrklen_el->SetLineColor(4);

	dtrklen_midp->SetFillColor(3); dtrklen_midp->SetLineColor(3);
	dtrklen_midcosmic->SetFillColor(5); dtrklen_midcosmic->SetLineColor(5);
	dtrklen_midpi->SetFillColor(6); dtrklen_midpi->SetLineColor(6);

	dtrklen_midmu->SetFillColor(28); dtrklen_midmu->SetLineColor(28);
	dtrklen_mideg->SetFillColor(30); dtrklen_mideg->SetLineColor(30);
	dtrklen_midother->SetFillColor(15); dtrklen_midother->SetLineColor(15);

	dtrklen_el->GetXaxis()->SetTitle("Reco Track Length - Truth Track Length [cm]");

	THStack* hs=new THStack("hs","");
	hs->Add(dtrklen_el);
	hs->Add(dtrklen_inel);

	hs->Add(dtrklen_midcosmic);
	hs->Add(dtrklen_midp);
	hs->Add(dtrklen_midpi);

	hs->Add(dtrklen_midmu);
	hs->Add(dtrklen_mideg);
	hs->Add(dtrklen_midother);


	TCanvas *c_=new TCanvas(Form("c"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_->Divide(1,1);
	c_->cd(1);
	TH2D *f2d=new TH2D("f2d",Form("%s",str_cut.Data()),200, -100, 100, 200000, 0, 200000);
	f2d->GetXaxis()->SetTitle("(Reco-Truth) Track Length [cm]");
	f2d->Draw();
	c_->cd(1)->SetLogy();

	hs->Draw("hist same");

	TLegend *leg = new TLegend(0.2,0.6,0.8,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(dtrklen_el, "El","f");
	leg->AddEntry(dtrklen_inel, "Inel","f");

	leg->AddEntry(dtrklen_midcosmic,"misID:cosmic","f");
	leg->AddEntry(dtrklen_midp, "misID:p","f");
	leg->AddEntry(dtrklen_midpi, "misID:#pi","f");

	leg->AddEntry(dtrklen_midmu,"misID:#mu","f");
	leg->AddEntry(dtrklen_mideg, "misID:e/#gamma","f");
	leg->AddEntry(dtrklen_midother, "misID:other","f");

	leg->SetNColumns(3);
	leg->Draw();

	c_->Print(Form("%s/dtrklen_%s.eps",fout_path.Data(),str_cut.Data()));

}
