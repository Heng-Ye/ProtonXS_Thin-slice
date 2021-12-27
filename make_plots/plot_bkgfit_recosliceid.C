//#include "../headers/BasicFunctions.h"
#include <vector>
#include <iostream>
#include "THStack.h"
#include "../headers/TemplateFitter.h"

void plot_bkgfit_recosliceid(TString fin, TString fin_data, TString fout_path, TString fout_filename, TString rep, TString title1) {

	//TString rep="h_recosliceid_allevts_cuts";
	TString title=Form("%s; Reco SliceID; Events", title1.Data());
	//TString title=Form("%s; Reco SliceID; Events", "All Protons");

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	//gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	gStyle->SetOptStat(0);
	gStyle->SetPadLeftMargin(0.13);

	//get data ----------------------------------------------------//
	TFile *f_data = TFile::Open(fin_data.Data());
	TH1D *h1d_data=(TH1D *)f_data->Get(Form("%s",rep.Data()));
	int n_data=h1d_data->Integral();
	cout<<"n_data:"<<n_data<<endl;

	//get mc -----------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fin.Data());
	TH1D *h1d_inel=(TH1D *)f_mc->Get(Form("%s_inel",rep.Data()));
	TH1D *h1d_el=(TH1D *)f_mc->Get(Form("%s_el", rep.Data()));
	TH1D *h1d_midcosmic=(TH1D *)f_mc->Get(Form("%s_midcosmic", rep.Data()));
	TH1D *h1d_midpi=(TH1D *)f_mc->Get(Form("%s_midpi", rep.Data()));
	TH1D *h1d_midp=(TH1D *)f_mc->Get(Form("%s_midp", rep.Data()));
	TH1D *h1d_midmu=(TH1D *)f_mc->Get(Form("%s_midmu", rep.Data()));
	TH1D *h1d_mideg=(TH1D *)f_mc->Get(Form("%s_mideg",rep.Data()));
	TH1D *h1d_midother=(TH1D *)f_mc->Get(Form("%s_midother", rep.Data()));

	h1d_inel->Sumw2();
	h1d_el->Sumw2();
h1d_midcosmic->Sumw2();
h1d_midpi->Sumw2();
h1d_midp->Sumw2();
h1d_midmu->Sumw2();
h1d_mideg->Sumw2();
h1d_midother->Sumw2();


	TH1D *h1d_mc=(TH1D*)h1d_inel->Clone(); h1d_mc->Sumw2();
	h1d_mc->SetName("h1d_mc");
	h1d_mc->Add(h1d_el);
	h1d_mc->Add(h1d_midcosmic);
	h1d_mc->Add(h1d_midpi);
	h1d_mc->Add(h1d_midp);
	h1d_mc->Add(h1d_midmu);
	h1d_mc->Add(h1d_mideg);
	h1d_mc->Add(h1d_midother);
	

	TH1D *MC_2=(TH1D *)f_mc->Get(Form("%s_midp", rep.Data())); 
	TH1D *MC_1=(TH1D *)h1d_inel->Clone(); MC_1->Sumw2();
	MC_1->Add(h1d_el);
	MC_1->Add(h1d_midcosmic);
	MC_1->Add(h1d_midpi);
	MC_1->Add(h1d_midmu);
	MC_1->Add(h1d_mideg);
	MC_1->Add(h1d_midother);


	h1d_inel->SetFillColor(2); h1d_inel->SetLineColor(2);
	h1d_el->SetFillColor(4); h1d_el->SetLineColor(4);
	h1d_midp->SetFillColor(3); h1d_midp->SetLineColor(3);
	h1d_midcosmic->SetFillColor(5); h1d_midcosmic->SetLineColor(5);
	h1d_midpi->SetFillColor(6); h1d_midpi->SetLineColor(6);
	h1d_midmu->SetFillColor(28); h1d_midmu->SetLineColor(28);
	h1d_mideg->SetFillColor(30); h1d_mideg->SetLineColor(30);
	h1d_midother->SetFillColor(15); h1d_midother->SetLineColor(15);

	int n_mc_inel=h1d_inel->Integral();
	int n_mc_el=h1d_el->Integral();
	int n_mc_midcosmic=h1d_midcosmic->Integral();
	int n_mc_midpi=h1d_midpi->Integral();
	int n_mc_midp=h1d_midp->Integral();
	int n_mc_midmu=h1d_midmu->Integral();
	int n_mc_mideg=h1d_mideg->Integral();
	int n_mc_midother=h1d_midother->Integral();
	int n_mc=n_mc_inel+n_mc_el+n_mc_midcosmic+n_mc_midpi+n_mc_midp+n_mc_midmu+n_mc_mideg+n_mc_midother;

	double norm_mc=(double)n_data/(double)n_mc;
	//double norm_mc=1;


	cout<<"n_data:"<<n_data<<endl;
	cout<<"n_mc:"<<n_mc<<endl;
	cout<<"h1d_inel:"<<h1d_inel->Integral()<<endl;

	
	h1d_inel->Scale(norm_mc);
	h1d_el->Scale(norm_mc);
	h1d_midcosmic->Scale(norm_mc);
	h1d_midpi->Scale(norm_mc);
	h1d_midp->Scale(norm_mc);
	h1d_midmu->Scale(norm_mc);
	h1d_mideg->Scale(norm_mc);
	h1d_midother->Scale(norm_mc);
	h1d_mc->Scale(norm_mc);

	MC_1->Scale(norm_mc);
	MC_2->Scale(norm_mc);

	cout<<"h1d_inel:"<<h1d_inel->Integral()<<endl;

	//prepare stack histo.
	THStack* hs=new THStack("hs","");
	hs->Add(h1d_inel);
	hs->Add(h1d_el);
	hs->Add(h1d_midp);
	hs->Add(h1d_midcosmic);
	hs->Add(h1d_midpi);
	hs->Add(h1d_midmu);
	hs->Add(h1d_mideg);
	hs->Add(h1d_midother);

	//h1d_mc->GetXaxis()->SetTitle("#chi^{2}PID");
	
	cout<<"n_mc:"<<n_mc<<endl;
	//cout<<"hs->Integral():"<<hs->Integral()<<endl;


        //data/MC -----------------------------------------//
        TH1D *R=(TH1D*)h1d_data->Clone();
        R->Divide(h1d_data, h1d_mc);
        //-------------------------------------------------//

	TCanvas *c_ = new TCanvas("c_", "c_", 900,700);
        TPad *pad1 = new TPad("pad1","pad1",0,0.29,1.,1.);
        pad1->SetBottomMargin(0.03);
        pad1->SetGridx();
        pad1->SetGridy();
        pad1->Draw();
        pad1->cd();
        //pad1->SetLogy();



	float xmin=-2;
	float xmax=26;
        float ymax=90;

	TH2D *f2d=new TH2D("f2d",Form("%s",""),100,xmin,xmax,3000,0,ymax); //logy
	f2d->SetTitle(title.Data());
	f2d->GetYaxis()->SetTitleOffset(1.);
	f2d->Draw();
	hs->Draw("hist same");
	h1d_mc->Draw("same hist");
	h1d_data->Draw("ep same");

	TLegend *leg = new TLegend(0.2,0.65,0.8,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(h1d_data, "Data","ep");
	leg->AddEntry(h1d_inel, "Inel","f");
	leg->AddEntry(h1d_el, "El","f");

	leg->AddEntry(h1d_midcosmic,"misID:cosmic","f");
	leg->AddEntry(h1d_midp, "misID:p","f");
	leg->AddEntry(h1d_midpi, "misID:#pi","f");

	leg->AddEntry(h1d_midmu,"misID:#mu","f");
	leg->AddEntry(h1d_mideg, "misID:e/#gamma","f");
	leg->AddEntry(h1d_midother, "misID:other","f");

	leg->SetNColumns(3);
	leg->Draw();


        c_->cd();
        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
        pad2->SetTopMargin(0.02);
        pad2->SetBottomMargin(0.25);
        pad2->SetGridx();
        pad2->SetGridy();
        pad2->Draw();
        pad2->cd();

        //TH2D *f2d2=new TH2D("f2d2",Form("%s","Pos"),100,xmin,xmax,64,0,6.4); //liny
        TH2D *f2d2=new TH2D("f2d2",Form("%s","Pos"),100,xmin,xmax,64,0,2.1); //liny
        f2d2->SetTitle(";Reco SliceID;Data/MC");
        f2d2->GetXaxis()->SetLabelSize(0.1);
        f2d2->GetYaxis()->SetLabelSize(0.1);

        f2d2->GetXaxis()->SetTitleSize(0.1);
        f2d2->GetYaxis()->SetTitleSize(0.1);
        f2d2->GetYaxis()->SetTitleOffset(.5);
        f2d2->Draw();
        TLine *line1=new TLine(xmin,1,xmax,1);
        line1->SetLineColor(7);
        line1->SetLineWidth(2);
        line1->Draw("same");
        R->Draw("ep same");



	c_->Print(Form("%s/%s.eps",fout_path.Data(), fout_filename.Data()));


/*
	//best-fit	
  	TemplateFitter fitter;
		
	//Min.(h0-h1-s*h2)
	
	//h0:data
	//h1:mc(mc except misID:p)
	//h2:mc(misID:p)
	vector<double> scal_midp;
	vector<double> err_scal_midp;
			
	fitter.SetHistograms(h1d_data, MC_1, MC_2);
    	//fitter.SetFitRange(4,-1);
    	fitter.Fit();
	scal_midp.push_back(fitter.GetPar());
    	err_scal_midp.push_back(fitter.GetParError());

	cout<<"scal_midp:"<<scal_midp.at(0)<<" +- "<<err_scal_midp.at(0)<<endl;
*/


}
