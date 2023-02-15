//#include "../headers/BasicAnaFunc.h"
#include <vector>
#include <iostream>
#include "../headers/SliceParams.h"
#include "../headers/TemplateFitter.h"
#include "RooFit.h"
#include "RooBifurGauss.h"

Double_t Quintic(Double_t *x, Double_t *par) {
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0]+par[3]*pow(x[0],3)+par[4]*pow(x[0],4)+par[5]*pow(x[0],5);
}

Double_t MSE(Double_t *x, Double_t *par) {
	return par[0]*(par[1]+TMath::Exp(par[2]*(x[0]-par[3])))*1./(1.+TMath::Exp((x[0]-par[4])/par[5]));
}

Double_t Rayleigh(Double_t *x, Double_t *par) {
	double m=par[0];
	double s=par[1];
	double n=par[2];

	double g=n*TMath::Exp(-(x[0]-m)*(x[0]-m)/(2*s*s))*(x[0]/pow(s,2));
	return g;
}

Double_t Weibull(Double_t *x, Double_t *par) {
	double a=par[0]; //scale parameter, or characteristic life
	double b=par[1]; //shape parameter (or slope)
	double c=par[2]; //location parameter (or failure free life)

	double g=b/a*pow((x[0]-c)/a,b-1)*TMath::Exp(-pow((x[0]-c)/a,b));
	return g;
}

Double_t ExtremeVal(Double_t *x, Double_t *par) {
	double a=par[0];
	double b=par[1];

	double g=1./b*TMath::Exp(((a-x[0])/b)-TMath::Exp((a-x[0])/b));
	return g;
}

Double_t fdistribution_pdf(Double_t *x, Double_t *par) {
	double n=par[0];
	double m=par[1];
	double x0=par[2];

	// Inlined to enable clad-auto-derivation for this function.
	// function is defined only for both n and m > 0
	if (n < 0 || m < 0)
		return std::numeric_limits<double>::quiet_NaN();
	if ((x[0]-x0) < 0)
		return 0.0;

	return std::exp((n/2) * std::log(n) + (m/2) * std::log(m) + ROOT::Math::lgamma((n+m)/2) - ROOT::Math::lgamma(n/2) - ROOT::Math::lgamma(m/2)
			+ (n/2 -1) * std::log(x[0]-x0) - ((n+m)/2) * std::log(m +  n*(x[0]-x0)) );

}

Double_t ReversedLognormal(Double_t *x, Double_t *par) {
	double delta=par[0];
	double eta=par[1];
	double gamma=par[2];
	return delta/((eta-x[0])*sqrt(2.*TMath::Pi()))*std::exp(-0.5*pow(gamma+delta*TMath::Log(eta-x[0]),2));
}


TF1* MIDP_SHAPE_LOGNORM(TH1D* h, float xmin, float xmax) {
	TF1 * f0 = new TF1("f0","[0]*ROOT::Math::lognormal_pdf(x,[1],[2],[3])",xmin,xmax);
	
	//set initial parameters
	double p[4];
	p[0]=h->GetEntries()*h->GetXaxis()->GetBinWidth(1); //area of hist
		
   	// find median of histogram 
   	double prob[] = {0.8}; 
   	double q[1]; 
   	h->GetQuantiles(1,q,prob);
   	double median = q[0];
   	// find mode of histogram 
   	double  mode = h->GetBinCenter( h->GetMaximumBin());
	std::cout << "histogram mode is " << mode  << " median is " << median << std::endl;
   	if (median < 0) { 
		median=-median;
      		//Error("lognormal","not valid histogram median");
      		//return;
   	}

   	// m is log(median)
   	p[1] = std::log(median); 

   	// s2 is  log(median) - log(mode) 
   	p[2] = std::sqrt( std::log(median/mode) ); 

	p[3] = 330;

   	f0->SetParameters(p); 	
	h->Fit("f0","remn");

	const int n_iter=20;
	TF1 **fx=new TF1*[n_iter];	
	for (int i=0; i<n_iter; ++i) {
		//fx[i] = new TF1(Form("fx_%d",i),ReversedLognormal,xmin,xmax,3);
		fx[i] = new TF1(Form("fx_%d",i),"[0]*ROOT::Math::lognormal_pdf(x,[1],[2],[3])",xmin,xmax); 
		if (i>0) for (int k=0; k<4; ++k) { fx[i]->SetParameter(k, fx[i-1]->GetParameter(k)); }
		else for (int k=0; k<4; ++k) { fx[i]->SetParameter(k, f0->GetParameter(k)); }
		h->Fit(fx[i],"remn");
	}

	return fx[n_iter-1];
}


TF1* MIDP_SHAPE_Weibull(TH1D* h, float xmin, float xmax) {
	//TF1 * f0 = new TF1("f0",ExtremeVal,xmin,xmax, 2);
	//f0->SetParameter(0,-350);

	TF1 *f0= new TF1("f0","ROOT::Math::fdistribution_pdf(x,[0],[1])",xmin,xmax);
	h->Fit("f0","remn");

	const int n_iter=20;
	TF1 **fx=new TF1*[n_iter];	
	for (int i=0; i<n_iter; ++i) {
		fx[i] = new TF1(Form("fx_%d",i),ExtremeVal,xmin,xmax,2);
		if (i>0) for (int k=0; k<2; ++k) { fx[i]->SetParameter(k, fx[i-1]->GetParameter(k)); }
		else for (int k=0; k<4; ++k) { fx[i]->SetParameter(k, f0->GetParameter(k)); }
		h->Fit(fx[i],"remn");
	}

	return fx[n_iter-1];

	/*
	   f0->SetParameter(2,350);
	   f0->SetParameter(0,80);
	//f0->SetParameter(1,0.15);
	f0->SetParLimits(0,0,999999999999);
	f0->SetParLimits(1,0,999999999999);
	h->Fit("f0","remn");

	TF1 * f1 = new TF1("f1",Weibull,xmin,xmax,3);
	for (int i=0; i<3; ++i) {
	f1->SetParameter(i, f0->GetParameter(i));
	}
	h->Fit("f1","remn");

	TF1 * f2 = new TF1("f2",Weibull,xmin,xmax,3);
	for (int i=0; i<3; ++i) {
	f2->SetParameter(i, f1->GetParameter(i));
	}
	h->Fit("f2","remn");

	return f2;
	*/

}


//Quintic regression ------------------------------------------------------------------------------------------------------//
TF1* MIDP_SHAPE(TH1D* h, float xmin, float xmax) {
	TF1 * f0 = new TF1("f0","2.*gaus(x,[0],[1],[2])*ROOT::Math::normal_cdf([3]*x,1,0)",xmin,xmax);
	f0->SetParameters(3.89529e+01,3.35908e+02,7.42532e+01,1.67572e+02);
	h->Fit("f0","remn");

	const int n_iter=5;
	TF1 **fx=new TF1*[n_iter];	
	for (int i=0; i<n_iter; ++i) {
		fx[i] = new TF1(Form("fx_%d",i),"2.*gaus(x,[0],[1],[2])*ROOT::Math::normal_cdf([3]*x,1,0)",xmin,xmax);
		if (i>0) for (int k=0; k<4; ++k) { fx[i]->SetParameter(k, fx[i-1]->GetParameter(k)); }
		else for (int k=0; k<4; ++k) { fx[i]->SetParameter(k, f0->GetParameter(k)); }
		h->Fit(fx[i],"remn");
	}

	return fx[n_iter-1];
}	

void plot_bkgfit_ke() {

	//Load Data/MC files --------------------------------------------------------------------------------------------------------------//
	//TString fin=Form("../mc_proton_beamxy_beammom_bmrw_by_kebeamff_bkg.root");
	//TString fin_data=Form("/dune/data2/users/hyliao/protonana/v09_39_01/KEHY_KEBEAMFF_BKG/proton_beamxy_beammom_bkg_runAll.root");

	TString fin=Form("../mc_proton_beamxy_beammom_bmrw_edepteloss.root");
	TString fin_data=Form("/dune/data2/users/hyliao/protonana/v09_39_01/KEFF_EDEPT/proton_beamxy_beammom_bkg_runAll.root");

	//Histogram Title -------------------------------------------------------------//
	TString ratio_title="; Proton Kinetic Energy at Interaction [MeV]; Data/MC";

	//Background Subtraction Setting ------------------------------//
	int proton_bkg_subtract=0; 
	//proton_bkg_subtract=1; //set 1 for misid:p to be subtracted
	//proton_bkg_subtract=0; //set 0 for el to be subtracted
	//proton_bkg_subtract=2; //set 2 for inel to be subtracted


	//KEend ---------------------------------------------------------------//
	//un-block to plot the specific interaction channels
	TString fout_path=Form("./plots_bkgstudy_ke/KEend_misidp.eps");
	TString observable="h1d_kend_bb_misidp";
	proton_bkg_subtract=1; //set 1 for misid:p to be subtracted
	TString title="Misidentified secondary proton candidates";

	//TString fout_path=Form("./plots_bkgstudy_ke/KEend_el.eps");
	//TString observable="h1d_kend_bb_el";
	//proton_bkg_subtract=0; //set 0 for el to be subtracted
	//TString title="Elastic-scattering proton candidates";

	//TString fout_path=Form("./plots_bkgstudy_ke/KEend_inel.eps");
	//TString observable="h1d_kend_bb_inel";
	//proton_bkg_subtract=2; //set 2 for inel to be subtracted
	//TString title="Inelastic-scattering proton candidates";

	//KEffbeam ---------------------------------------------------------//
	//TString fout_path=Form("./plots_bkgstudy_ke/KEffbeam_inel.eps");
	//TString observable="h1d_keffbeam_inel";
	//proton_bkg_subtract=2; //set 2 for inel to be subtracted
	//TString title="Inelastic-scattering proton candidates";

	//TString fout_path=Form("./plots_bkgstudy_ke/KEffbeam_el.eps");
	//TString observable="h1d_keffbeam_el";
	//proton_bkg_subtract=0; //set 0 for el to be subtracted
	//TString title="Elastic-scattering proton candidates";

	//TString fout_path=Form("./plots_bkgstudy_ke/KEffbeam_misidp.eps");
	//TString observable="h1d_keffbeam_misidp";
	//proton_bkg_subtract=1; //set 1 for misid:p to be subtracted
	//TString title="Misidentified secondary proton candidates";



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

	//get mc -----------------------------------------------------------------------------------------------------------------------//
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

	//Rebin ----------------------//
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

	//MC scaling
	h_mc->Scale(norm_mc);
	h_mc_inel->Scale(norm_mc);
	h_mc_el->Scale(norm_mc);
	h_mc_midp->Scale(norm_mc);
	h_mc_midcosmic->Scale(norm_mc);
	h_mc_midpi->Scale(norm_mc);
	h_mc_midmu->Scale(norm_mc);
	h_mc_mideg->Scale(norm_mc);
	h_mc_midother->Scale(norm_mc);

	//for (int i=0; i<h_data->GetNbinsX(); ++i) {
	//cout<<h_data->GetBinCenter(i)<<" "<<h_data->GetBinContent(i)<<endl;
	//}

	//Prepare MCs to be fitted ----------------------------- 
	//Min.(h0-h1-s*h2)
	//h0:data
	//h1:mc(mc except misID:p)
	//h2:mc(misID:p)
	TH1D *MC1=(TH1D *)h_mc->Clone();
	TH1D *MC2;
	if (proton_bkg_subtract==1) MC2=(TH1D *)h_mc_midp->Clone();
	if (proton_bkg_subtract==0) MC2=(TH1D *)h_mc_el->Clone();
	if (proton_bkg_subtract==2) MC2=(TH1D *)h_mc_inel->Clone();
	MC1->Add(MC2,-1);

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



	//data/mc ratio ----------------------------------
	TH1D *R_raw=(TH1D*)h_data->Clone();
	R_raw->Divide(h_data, h_mc);
	R_raw->SetTitle(Form("%s",ratio_title.Data()));

	for (int i=0; i<h_data->GetNbinsX(); ++i) {
		cout<<"\nE["<<i<<"]="<<h_data->GetBinCenter(i)<<endl;	
		cout<<"h_data["<<i<<"]="<<h_data->GetBinContent(i)<<" Err="<<h_data->GetBinError(i)<<endl;	
		cout<<"h_mc["<<i<<"]="<<h_mc->GetBinContent(i)<<" Err="<<h_mc->GetBinError(i)<<endl;	

		//cross-check err
		double tmp_err_r1=pow(h_data->GetBinError(i)/h_data->GetBinContent(i),2);
		double tmp_err_r2=pow(h_mc->GetBinError(i)/h_mc->GetBinContent(i),2);
		double tmp_r=(h_data->GetBinContent(i)/h_mc->GetBinContent(i))*sqrt(tmp_err_r1+tmp_err_r2);
		cout<<"R_raw["<<i<<"]="<<R_raw->GetBinContent(i)<<" Err="<<R_raw->GetBinError(i)<<" tmp_r="<<tmp_r<<endl;	

	}

	//draw fig
	TCanvas *c_ = new TCanvas("c_", "c_", 900,700);
	TPad *pad1 = new TPad("pad1","pad1",0,0.29,1.,1.);
	pad1->SetBottomMargin(0.03);
	pad1->SetGridx();
	pad1->SetGridy();
	pad1->SetLogy(0);
	pad1->Draw();
	pad1->cd();

	//KEend ----------------------------------//
	//float xmin=-100;
	float xmin=-10;
	float xmax=550;
	float ymax=150; //for lin y of misid:p
	//float ymax=1550; //for lin y of inel
	//float ymax=2200; //for liny of el (after rebin=5)
	//float ymax=18000; //for liny of el

	//KEff ----------------------------------//
	//float xmin=250;
	//float xmax=600;
	//float ymax=3200; //for lin y of inel,el
	//float ymax=300; //for lin y of misidp 


	//TString title=Form("SliceID %d",k-1);
	TH2D *f2d=new TH2D("f2d",Form("%s",""),100,xmin,xmax,3000,0,ymax); //logy
	f2d->SetTitle(title.Data());
	f2d->GetYaxis()->SetTitleOffset(1.);
	f2d->Draw();
	hs->Draw("same hist");
	h_data->Draw("ep same");

	TF1 *misidp_shape_data=MIDP_SHAPE_LOGNORM(h_data, 0, 540);
	misidp_shape_data->SetLineColor(1);
	misidp_shape_data->SetLineStyle(2);
	misidp_shape_data->Draw("same");

	//h_mc->SetLineColor(1);
	//h_mc->Draw("same hist");
	TF1 *misidp_shape_mc=MIDP_SHAPE_LOGNORM(h_mc, 0, 540);
	misidp_shape_mc->SetLineColor(4);
	misidp_shape_mc->SetLineStyle(2);
	misidp_shape_mc->Draw("same");


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


	//Best-fit[const-E-loss] ==============================================================
	//KEend_bb ------------
	//[1]misid:p protons:
	//NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
	//1  corr_fact    9.98892e-01   6.96991e-02   7.88407e-05  -1.64638e-04
	//EXTERNAL ERROR MATRIX.    NDIM=  25    NPAR=  1    ERR DEF=1
	//4.859e-03 
	//Best fit = 0.998892 error = 0.0696991
	//-----------------------------------------------------------------------
	//[2]el protons:
	//EXT PARAMETER                                   STEP         FIRST   
	//NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
	//1  corr_fact    9.99491e-01   1.66378e-02   3.34628e-05   4.41530e-04
	//EXTERNAL ERROR MATRIX.    NDIM=  25    NPAR=  1    ERR DEF=1
	//2.768e-04 
	//Best fit = 0.999491 error = 0.0166378
	//------------------------------------------------------------------------
	//[3]inel protons:
	//EXT PARAMETER                                   STEP         FIRST   
	//NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
	//1  corr_fact    9.99526e-01   1.26930e-02   1.65193e-05   3.29811e-04
	//EXTERNAL ERROR MATRIX.    NDIM=  25    NPAR=  1    ERR DEF=1
	//1.611e-04 
	//Best fit = 0.999526 error = 0.012693
	//=========================================================================
	//KEff=KEbeam-constE--------------
	//[1]inel
	//EXT PARAMETER                                   STEP         FIRST   
	//NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
	//1  corr_fact    9.96384e-01   1.26871e-02   2.55850e-05   1.92184e-02
	//EXTERNAL ERROR MATRIX.    NDIM=  25    NPAR=  1    ERR DEF=1
	//1.610e-04 
	//Best fit = 0.996384 error = 0.0126871
	//-------------------------------------------------------------------------
	//[2]el
	//EXT PARAMETER                                   STEP         FIRST   
	//NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
	//1  corr_fact    9.23894e-01   9.46050e-02   2.36701e-04   8.31353e-02
	//EXTERNAL ERROR MATRIX.    NDIM=  25    NPAR=  1    ERR DEF=1
	//8.953e-03 
	//Best fit = 0.923894 error = 0.094605
	//-------------------------------------------------------------------------
	//[3]misid:p
	//EXT PARAMETER                                   STEP         FIRST   
	//NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
	//1  corr_fact    1.24328e+00   5.00595e-01   4.98813e-04   8.84935e-05
	//EXTERNAL ERROR MATRIX.    NDIM=  25    NPAR=  1    ERR DEF=1
	//2.525e-01 
	//Best fit = 1.24328 error = 0.500595
	//=========================================================================

	//Best-fit[Edept-E-loss] ==============================================================
	//KEend_bb ------------
	//[1]misid:p protons:
	//EXT PARAMETER                                   STEP         FIRST   
	//NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
	//1  corr_fact    9.82934e-01   7.05165e-02   8.38410e-05  -4.35257e-02
	//EXTERNAL ERROR MATRIX.    NDIM=  25    NPAR=  1    ERR DEF=1
	//4.974e-03 
	//Best fit = 0.982934 error = 0.0705165
	//[2]el protons:
	//NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
	//1  corr_fact    9.96484e-01   1.66868e-02   2.97235e-05   1.88657e-02
	//EXTERNAL ERROR MATRIX.    NDIM=  25    NPAR=  1    ERR DEF=1
	//2.785e-04 
	//Best fit = 0.996484 error = 0.0166868
	//[3]inel protons:
	//EXT PARAMETER                                   STEP         FIRST   
	//NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
	//1  corr_fact    9.99085e-01   1.28888e-02   1.21266e-05   8.91422e-04
	//EXTERNAL ERROR MATRIX.    NDIM=  25    NPAR=  1    ERR DEF=1
	//1.661e-04 
	//Best fit = 0.999085 error = 0.0128888
	//=========================================================================
	//KEff=KEbeam-Edept--------------
	//[1]inel
	//EXT PARAMETER                                   STEP         FIRST   
	//NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
	//1  corr_fact    9.95818e-01   1.28686e-02   2.24885e-05   1.85718e-02
	//EXTERNAL ERROR MATRIX.    NDIM=  25    NPAR=  1    ERR DEF=1
	//1.656e-04 
	//Best fit = 0.995818 error = 0.0128686
	//[2]el protons:
	//EXT PARAMETER                                   STEP         FIRST   
	//NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
	//1  corr_fact    9.93652e-01   1.44007e-02   2.62226e-05   7.34494e-02
	//EXTERNAL ERROR MATRIX.    NDIM=  25    NPAR=  1    ERR DEF=1
	//2.074e-04 
	//Best fit = 0.993652 error = 0.0144007
	//[3]misid:p
	//NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
	//1  corr_fact    9.97398e-01   7.04523e-02   5.73121e-05  -8.24530e-04
	//EXTERNAL ERROR MATRIX.    NDIM=  25    NPAR=  1    ERR DEF=1
	//4.964e-03 
	//Best fit = 0.997398 error = 0.0704523
	//=========================================================================


	c_->Print(Form("%s",fout_path.Data()));


}
