//#include "../headers/BasicAnaFunc.h"
#include <vector>
#include <iostream>
#include "../headers/SliceParams.h"
#include "../headers/TemplateFitter.h"

void plot_bkgfit_ke() {

	TString fin=Form("../mc_proton_beamxy_beammom_bmrw_by_kebeamff_bkg.root");
	TString fin_data=Form("/dune/data2/users/hyliao/protonana/v09_39_01/KEHY_KEBEAMFF_BKG/proton_beamxy_beammom_bkg_runAll.root");

	TString fout_path=Form("./plots_bkgstudy_ke/KEend_misidp.eps");
	TString observable="h1d_kend_bb_misidp";
	TString title="Misidentified secondary proton candidates";

	//TString fout_path=Form("./plots_bkgstudy_ke/KEend_el.eps");
	//TString observable="h1d_kend_bb_el";
	//TString title="Elastic-scattering proton candidates";

	TString ratio_title="; Proton Kinetic Energy at Interaction [MeV]; Data/MC";
	bool El_or_Misidp=0; 
	El_or_Misidp=1; //set 1 for misid:p to be subtracted
	//El_or_Misidp=0; //set 0 for el to be subtracted

	//plot style -------------------------------------------------//
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	//gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	gStyle->SetOptStat(0);
	gStyle->SetPadLeftMargin(0.13);

	//get data -------------------------------------------------------------------------------------//
	TFile *f_data = TFile::Open(fin_data.Data());
	TH1D *h0_data=(TH1D *)f_data->Get(Form("%s",observable.Data()));

	//get mc -----------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fin.Data());

	TH1D *h0_mc=(TH1D *)f_mc->Get(Form("%s",observable.Data()));
	TH1D *h0_mc_inel=(TH1D *)f_mc->Get(Form("%s_inel",observable.Data()));
	TH1D *h0_mc_el=(TH1D *)f_mc->Get(Form("%s_el",observable.Data()));
	TH1D *h0_mc_midcosmic=(TH1D *)f_mc->Get(Form("%s_midcosmic",observable.Data()));
	TH1D *h0_mc_midpi=(TH1D *)f_mc->Get(Form("%s_midpi",observable.Data()));
	TH1D *h0_mc_midp=(TH1D *)f_mc->Get(Form("%s_midp",observable.Data()));
	TH1D *h0_mc_midmu=(TH1D *)f_mc->Get(Form("%s_midmu",observable.Data()));
	TH1D *h0_mc_mideg=(TH1D *)f_mc->Get(Form("%s_mideg",observable.Data()));
	TH1D *h0_mc_midother=(TH1D *)f_mc->Get(Form("%s_midother",observable.Data()));
	
	TH1D *h_data=new TH1D("h_data","",h0_data->GetNbinsX(),h0_data->GetXaxis()->GetXmin(),h0_data->GetXaxis()->GetXmax());
	TH1D *h_mc=new TH1D("h_mc","",h0_mc->GetNbinsX(),h0_mc->GetXaxis()->GetXmin(),h0_mc->GetXaxis()->GetXmax());
	TH1D *h_mc_inel=new TH1D("h_mc_inel","",h0_mc->GetNbinsX(),h0_mc->GetXaxis()->GetXmin(),h0_mc->GetXaxis()->GetXmax());
	TH1D *h_mc_el=new TH1D("h_mc_el","",h0_mc->GetNbinsX(),h0_mc->GetXaxis()->GetXmin(),h0_mc->GetXaxis()->GetXmax());
	TH1D *h_mc_midcosmic=new TH1D("h_mc_midcosmic","",h0_mc->GetNbinsX(),h0_mc->GetXaxis()->GetXmin(),h0_mc->GetXaxis()->GetXmax());
	TH1D *h_mc_midpi=new TH1D("h_mc_midpi","",h0_mc->GetNbinsX(),h0_mc->GetXaxis()->GetXmin(),h0_mc->GetXaxis()->GetXmax());
	TH1D *h_mc_midp=new TH1D("h_mc_midp","",h0_mc->GetNbinsX(),h0_mc->GetXaxis()->GetXmin(),h0_mc->GetXaxis()->GetXmax());
	TH1D *h_mc_midmu=new TH1D("h_mc_midmu","",h0_mc->GetNbinsX(),h0_mc->GetXaxis()->GetXmin(),h0_mc->GetXaxis()->GetXmax());
	TH1D *h_mc_mideg=new TH1D("h_mc_mideg","",h0_mc->GetNbinsX(),h0_mc->GetXaxis()->GetXmin(),h0_mc->GetXaxis()->GetXmax());
	TH1D *h_mc_midother=new TH1D("h_mc_midother","",h0_mc->GetNbinsX(),h0_mc->GetXaxis()->GetXmin(),h0_mc->GetXaxis()->GetXmax());

	for (int i=0; i<h0_data->GetNbinsX(); ++i) {
		float ebin=h0_data->GetXaxis()->GetBinCenter(i);
		if (ebin>-50) {
			h_data->SetBinContent(i,h0_data->GetBinContent(i));

			h_mc->SetBinContent(i,h0_mc->GetBinContent(i));
			h_mc_inel->SetBinContent(i,h0_mc_inel->GetBinContent(i));
			h_mc_el->SetBinContent(i,h0_mc_el->GetBinContent(i));
			h_mc_midp->SetBinContent(i,h0_mc_midp->GetBinContent(i));
			h_mc_midcosmic->SetBinContent(i,h0_mc_midcosmic->GetBinContent(i));
			h_mc_midpi->SetBinContent(i,h0_mc_midpi->GetBinContent(i));
			h_mc_midmu->SetBinContent(i,h0_mc_midmu->GetBinContent(i));
			h_mc_mideg->SetBinContent(i,h0_mc_mideg->GetBinContent(i));
			h_mc_midother->SetBinContent(i,h0_mc_midother->GetBinContent(i));
		}	

	}

	std::cout<<"nbin of h_data:"<<h_data->GetNbinsX()<<std::endl;
	std::cout<<"xmin [data]:"<<h_data->GetXaxis()->GetXmin()<<std::endl;
	std::cout<<"xmax [data]:"<<h_data->GetXaxis()->GetXmax()<<std::endl;

	std::cout<<"nbin of h_mc:"<<h_mc->GetNbinsX()<<std::endl;
	std::cout<<"xmin [mc]:"<<h_mc->GetXaxis()->GetXmin()<<std::endl;
	std::cout<<"xmax [mc]:"<<h_mc->GetXaxis()->GetXmax()<<std::endl;
/*
	float range_xmin=0;
	float range_xmax=1000;
	h_data->GetXaxis()->SetRange(range_xmin,range_xmax);
	h_mc->GetXaxis()->SetRange(range_xmin,range_xmax);
	h_mc_inel->GetXaxis()->SetRange(range_xmin,range_xmax);
	h_mc_el->GetXaxis()->SetRange(range_xmin,range_xmax);
	h_mc_midcosmic->GetXaxis()->SetRange(range_xmin,range_xmax);
	h_mc_midpi->GetXaxis()->SetRange(range_xmin,range_xmax);
	h_mc_midp->GetXaxis()->SetRange(range_xmin,range_xmax);
	h_mc_midmu->GetXaxis()->SetRange(range_xmin,range_xmax);
	h_mc_mideg->GetXaxis()->SetRange(range_xmin,range_xmax);
	h_mc_midother->GetXaxis()->SetRange(range_xmin,range_xmax);
*/

	h_mc_inel->SetFillColor(2);
	h_mc_el->SetFillColor(4);
	h_mc_midp->SetFillColor(3);
	h_mc_midcosmic->SetFillColor(5);
	h_mc_midpi->SetFillColor(6);
	h_mc_midmu->SetFillColor(28);
	h_mc_mideg->SetFillColor(30);
	h_mc_midother->SetFillColor(15);

	h_mc_inel->SetLineColor(2);
	h_mc_el->SetLineColor(4);
	h_mc_midp->SetLineColor(3);
	h_mc_midcosmic->SetLineColor(5);
	h_mc_midpi->SetLineColor(6);
	h_mc_midmu->SetLineColor(28);
	h_mc_mideg->SetLineColor(30);
	h_mc_midother->SetLineColor(15);

	int n_data=h_data->Integral();
	cout<<"n_data:"<<n_data<<endl;

	int n_mc=h_mc->Integral();
	std::cout<<"n_mc:"<<n_mc<<std::endl;

	int n_mc_inel=h_mc_inel->Integral();
	int n_mc_el=h_mc_el->Integral();
	int n_mc_midp=h_mc_midp->Integral();
	int n_mc_midcosmic=h_mc_midcosmic->Integral();
	int n_mc_midpi=h_mc_midpi->Integral();
	int n_mc_midmu=h_mc_midmu->Integral();
	int n_mc_mideg=h_mc_mideg->Integral();
	int n_mc_midother=h_mc_midother->Integral();
	//int n_mc=n_mc_inel+n_mc_el+n_mc_midp+n_mc_midcosmic+n_mc_midpi+n_mc_midmu+n_mc_mideg+n_mc_midother;

	//std::cout<<"n_mc_inel+n_mc_el+n_mc_midp+n_mc_midcosmic+n_mc_midpi+n_mc_midmu+n_mc_mideg+n_mc_midother="<<n_mc_inel+n_mc_el+n_mc_midp+n_mc_midcosmic+n_mc_midpi+n_mc_midmu+n_mc_mideg+n_mc_midother<<std::endl;

	double norm_mc=(double)n_data/(double)n_mc;
	cout<<"n_mc:"<<n_mc<<endl;

	//Rebin
	int n_rb=5;
	h_data->Rebin(n_rb);
	h_mc->Rebin(n_rb);
	h_mc_inel->Rebin(n_rb);
	h_mc_el->Rebin(n_rb);
	h_mc_midp->Rebin(n_rb);
	h_mc_midcosmic->Rebin(n_rb);
	h_mc_midpi->Rebin(n_rb);
	h_mc_midmu->Rebin(n_rb);
	h_mc_mideg->Rebin(n_rb);
	h_mc_midother->Rebin(n_rb);

	TH1D *MC1=(TH1D *)h_mc->Clone();
	TH1D *MC2;
	if (El_or_Misidp==true) MC2=(TH1D *)h_mc_midp->Clone();
	if (El_or_Misidp==false) MC2=(TH1D *)h_mc_el->Clone();
	MC1->Add(MC2,-1);
	//MC1->Add(h_mc_midp,-1);
	MC2->Scale(norm_mc);
	MC1->Scale(norm_mc);

	//stack histogram --------------------
	THStack* hs=new THStack("hs","");
	hs->Add(h_mc_inel);
	hs->Add(h_mc_el);
	hs->Add(h_mc_midp);
	hs->Add(h_mc_midcosmic);
	hs->Add(h_mc_midpi);
	hs->Add(h_mc_midmu);
	hs->Add(h_mc_mideg);
	hs->Add(h_mc_midother);

	//data/mc ratio --------------------------
       	TH1D *R_raw=(TH1D*)h_data->Clone();
       	R_raw->Divide(h_data, h_mc);
	R_raw->SetTitle(Form("%s",ratio_title.Data()));
	
	//draw fig
	TCanvas *c_ = new TCanvas("c_", "c_", 900,700);
        TPad *pad1 = new TPad("pad1","pad1",0,0.29,1.,1.);
        pad1->SetBottomMargin(0.03);
        pad1->SetGridx();
        pad1->SetGridy();
        pad1->SetLogy(0);
        pad1->Draw();
        pad1->cd();

	float xmin=-100;
	float xmax=550;
       	float ymax=1500; //for log y of misid:p
       	//float ymax=18000; //for liny of el

	//TString title=Form("SliceID %d",k-1);
	TH2D *f2d=new TH2D("f2d",Form("%s",""),100,xmin,xmax,3000,0,ymax); //logy
	f2d->SetTitle(title.Data());
	f2d->GetYaxis()->SetTitleOffset(1.);
	f2d->Draw();
	hs->Draw("same hist");
	h_data->Draw("ep same");

	TLegend *leg = new TLegend(0.2,0.65,0.8,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(h_data, "Data","ep");
	leg->AddEntry(h_mc_inel, "Inel","f");
	leg->AddEntry(h_mc_el, "El","f");

	leg->AddEntry(h_mc_midcosmic,"misID:cosmic","f");
	leg->AddEntry(h_mc_midp, "misID:p","f");
	leg->AddEntry(h_mc_midpi, "misID:#pi","f");

	leg->AddEntry(h_mc_midmu,"misID:#mu","f");
	leg->AddEntry(h_mc_mideg, "misID:e/#gamma","f");
	leg->AddEntry(h_mc_midother, "misID:other","f");

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

       	TH2D *f2d2=new TH2D("f2d2",Form("%s",""),100,xmin,xmax,64,-1,5.); //liny for fake data
       	f2d2->SetTitle(Form("%s",ratio_title.Data()));
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
       	R_raw->Draw("ep same");

	//fitting to derive the scaling factor --------------------
  	TemplateFitter fitter;
		
	//Min.(h0-h1-s*h2)
	//h0:data
	//h1:mc(mc except misID:p)
	//h2:mc(misID:p)			
	fitter.SetHistograms(h_data, MC1, MC2);
    	//fitter.SetFitRange(4,-1);
    	fitter.Fit();

	vector<double> scal_fit;
	vector<double> err_scal_fit;
	scal_fit.push_back(fitter.GetPar());
    	err_scal_fit.push_back(fitter.GetParError());

        TLine *line_fit=new TLine(xmin, scal_fit.at(0), xmax, scal_fit.at(0));
        TLine *line_fit_up=new TLine(xmin, scal_fit.at(0)+err_scal_fit.at(0), xmax, scal_fit.at(0)+err_scal_fit.at(0));
        TLine *line_fit_dn=new TLine(xmin, scal_fit.at(0)-err_scal_fit.at(0), xmax, scal_fit.at(0)-err_scal_fit.at(0));

        line_fit->SetLineColor(2);
        line_fit->SetLineWidth(2);

        line_fit_up->SetLineColor(2);
        line_fit_up->SetLineWidth(2);
        line_fit_up->SetLineStyle(2);

        line_fit_dn->SetLineColor(2);
        line_fit_dn->SetLineWidth(2);
        line_fit_dn->SetLineStyle(2);

        line_fit->Draw("same");
	line_fit_up->Draw("same");
	line_fit_dn->Draw("same");


	//Best-fit -------------------------------------------------------------
	//[1]misid:p protons:
	//KEend_bb
	//Best fit = 1.03137 error = 0.0772546

 	//EXT PARAMETER                                   STEP         FIRST   
  	//NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   	//1  corr_fact    1.03137e+00   7.72546e-02   8.45675e-05  -1.29152e-01


	c_->Print(Form("%s",fout_path.Data()));






}
