#include "THStack.h"

void plot_hdeltaXYZ(TString fin, TString fin_data, TString fout_path, TString str_obs) {

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	//gStyle->SetOptStat(0);

	//data --------------------------------------------------------------------------------------//
	TFile *f_data = TFile::Open(fin_data.Data());
	TH1D *h1d=(TH1D *)f_data->Get(Form("hdelta%s",str_obs.Data()));
	int n_data=h1d->Integral();

	//mc ----------------------------------------------------------------------------------------//
	TFile *f0 = TFile::Open(fin.Data());

	TH1D *h1d_inel=(TH1D *)f0->Get(Form("hdelta%s_inel",str_obs.Data()));
	TH1D *h1d_el=(TH1D *)f0->Get(Form("hdelta%s_el",str_obs.Data()));

	TH1D *h1d_midcosmic=(TH1D *)f0->Get(Form("hdelta%s_midcosmic",str_obs.Data()));
	TH1D *h1d_midpi=(TH1D *)f0->Get(Form("hdelta%s_midpi",str_obs.Data()));
	TH1D *h1d_midp=(TH1D *)f0->Get(Form("hdelta%s_midp",str_obs.Data()));

	TH1D *h1d_midmu=(TH1D *)f0->Get(Form("hdelta%s_midmu",str_obs.Data()));
	TH1D *h1d_mideg=(TH1D *)f0->Get(Form("hdelta%s_mideg",str_obs.Data()));
	TH1D *h1d_midother=(TH1D *)f0->Get(Form("hdelta%s_midother",str_obs.Data()));

	h1d_inel->SetFillColor(2); h1d_inel->SetLineColor(2);
	h1d_el->SetFillColor(4); h1d_el->SetLineColor(4);

	h1d_midp->SetFillColor(3); h1d_midp->SetLineColor(3);
	h1d_midcosmic->SetFillColor(5); h1d_midcosmic->SetLineColor(5);
	h1d_midpi->SetFillColor(6); h1d_midpi->SetLineColor(6);

	h1d_midmu->SetFillColor(28); h1d_midmu->SetLineColor(28);
	h1d_mideg->SetFillColor(30); h1d_mideg->SetLineColor(30);
	h1d_midother->SetFillColor(15); h1d_midother->SetLineColor(15);

	int n_inel=h1d_inel->Integral();
	int n_el=h1d_el->Integral();
	int n_midcosmic=h1d_midcosmic->Integral();
	int n_midpi=h1d_midpi->Integral();
	int n_midp=h1d_midp->Integral();
	int n_midmu=h1d_midmu->Integral();
	int n_mideg=h1d_mideg->Integral();
	int n_midother=h1d_midother->Integral();
	int n_mc=n_inel+n_el+n_midcosmic+n_midpi+n_midp+n_midmu+n_mideg+n_midother;
	double mc_scale=(double)n_data/(double)n_mc;
	h1d_inel->Scale(mc_scale);	
	h1d_el->Scale(mc_scale);	
	h1d_midp->Scale(mc_scale);	
	h1d_midcosmic->Scale(mc_scale);	
	h1d_midpi->Scale(mc_scale);	
	h1d_midmu->Scale(mc_scale);	
	h1d_mideg->Scale(mc_scale);	
	h1d_midother->Scale(mc_scale);	


	//h1d_el->GetXaxis()->SetTitle("");
	h1d_inel->GetXaxis()->SetTitle(Form("#Delta%s/#sigma_{%s}",str_obs.Data(), str_obs.Data()));

	THStack* hs=new THStack("hs","");
	hs->Add(h1d_inel);
	hs->Add(h1d_el);

	hs->Add(h1d_midcosmic);
	hs->Add(h1d_midp);
	hs->Add(h1d_midpi);

	hs->Add(h1d_midmu);
	hs->Add(h1d_mideg);
	hs->Add(h1d_midother);


	TCanvas *c_=new TCanvas(Form("c"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_->Divide(1,1);
	c_->cd(1);
	float ymax=1500;
	if (str_obs.Data()==Form("%s","XY")) ymax=2100;

	TH2D *f2d=new TH2D("f2d",Form("%s",str_obs.Data()),100, -10, 10,2000,0,ymax);
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
	h1d->Draw("ep same");

	TLegend *leg = new TLegend(0.2,0.6,0.8,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(h1d, "Data","ep");
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

}
