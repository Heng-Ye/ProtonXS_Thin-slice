#include "THStack.h"

void plot_trklentrue(TString fin, TString fout_path, TString str_cut) {

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	//gStyle->SetOptStat(0);
	
	//mc ----------------------------------------------------------------------------------------//
	TFile *f0 = TFile::Open(fin.Data());

	TH1D *trklen_inel=(TH1D *)f0->Get(Form("h1d_truetrklen_inel_%s",str_cut.Data()));
	TH1D *trklen_el=(TH1D *)f0->Get(Form("h1d_truetrklen_el_%s",str_cut.Data()));

	TH1D *trklen_midcosmic=(TH1D *)f0->Get(Form("h1d_truetrklen_midcosmic_%s",str_cut.Data()));
	TH1D *trklen_midpi=(TH1D *)f0->Get(Form("h1d_truetrklen_midpi_%s",str_cut.Data()));
	TH1D *trklen_midp=(TH1D *)f0->Get(Form("h1d_truetrklen_midp_%s",str_cut.Data()));

	TH1D *trklen_midmu=(TH1D *)f0->Get(Form("h1d_truetrklen_midmu_%s",str_cut.Data()));
	TH1D *trklen_mideg=(TH1D *)f0->Get(Form("h1d_truetrklen_mideg_%s",str_cut.Data()));
	TH1D *trklen_midother=(TH1D *)f0->Get(Form("h1d_truetrklen_midother_%s",str_cut.Data()));

	trklen_inel->SetFillColor(2); trklen_inel->SetLineColor(2);
	trklen_el->SetFillColor(4); trklen_el->SetLineColor(4);

	trklen_midp->SetFillColor(3); trklen_midp->SetLineColor(3);
	trklen_midcosmic->SetFillColor(5); trklen_midcosmic->SetLineColor(5);
	trklen_midpi->SetFillColor(6); trklen_midpi->SetLineColor(6);

	trklen_midmu->SetFillColor(28); trklen_midmu->SetLineColor(28);
	trklen_mideg->SetFillColor(30); trklen_mideg->SetLineColor(30);
	trklen_midother->SetFillColor(15); trklen_midother->SetLineColor(15);

	int n_mc_inel=trklen_inel->Integral();
	int n_mc_el=trklen_el->Integral();
	int n_mc_midcosmic=trklen_midcosmic->Integral();
	int n_mc_midpi=trklen_midpi->Integral();
	int n_mc_midp=trklen_midpi->Integral();
	int n_mc_midmu=trklen_midmu->Integral();
	int n_mc_mideg=trklen_mideg->Integral();
	int n_mc_midother=trklen_midother->Integral();
	int n_mc=n_mc_inel+n_mc_el+n_mc_midcosmic+n_mc_midpi+n_mc_midp+n_mc_midmu+n_mc_mideg+n_mc_midother;

	//double norm_mc=(double)n_data/(double)n_mc;
	double norm_mc=1;
	trklen_inel->Scale(norm_mc);
	trklen_el->Scale(norm_mc);
	trklen_midcosmic->Scale(norm_mc);
	trklen_midpi->Scale(norm_mc);
	trklen_midp->Scale(norm_mc);
	trklen_midmu->Scale(norm_mc);
	trklen_mideg->Scale(norm_mc);
	trklen_midother->Scale(norm_mc);

	//--------------------------------------------------------------------------------------------//

	trklen_inel->GetXaxis()->SetTitle("Truth Track Length [cm]");

	THStack* hs=new THStack("hs","");
	hs->Add(trklen_inel);
	hs->Add(trklen_el);

	hs->Add(trklen_midcosmic);
	hs->Add(trklen_midp);
	hs->Add(trklen_midpi);

	hs->Add(trklen_midmu);
	hs->Add(trklen_mideg);
	hs->Add(trklen_midother);


	TCanvas *c_=new TCanvas(Form("c"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_->Divide(1,1);
	c_->cd(1);
	TH2D *f2d=new TH2D("f2d",Form("%s",str_cut.Data()),136,-4,132,1100,0,1100);
	f2d->GetXaxis()->SetTitle("Truth Track Length [cm]");
	f2d->Draw();
	//c_->cd(1)->SetLogy();
	hs->Draw("hist same");
	//trklen_data->Draw("ep same");

	TLegend *leg = new TLegend(0.2,0.6,0.8,0.9);
	leg->SetFillStyle(0);
	//leg->AddEntry(trklen_data, "Data", "ep");
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

	c_->Print(Form("%s/trklentruth_%s.eps",fout_path.Data(),str_cut.Data()));

}
