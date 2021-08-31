#include "THStack.h"
#include "../headers/BasicParameters.h"
#include "../headers/BasicFunctions.h"

void plot_2d_ntrklen_chi2pid() {
	//TString rep="Pos";	
	TString rep="BQ";	

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);

	TFile *f_mc = TFile::Open("../mc_proton_studyRecoInelCut.root");
	TH2D* h2d=(TH2D*)f_mc->Get(Form("ntrklen_chi2pid_%s", rep.Data()));
	TH2D* h2d_inel=(TH2D*)f_mc->Get(Form("ntrklen_chi2pid_%s_inel", rep.Data()));
	TH2D* h2d_el=(TH2D*)f_mc->Get(Form("ntrklen_chi2pid_%s_el", rep.Data()));
	TH2D* h2d_midp=(TH2D*)f_mc->Get(Form("ntrklen_chi2pid_%s_midp", rep.Data()));

	h2d_inel->SetMarkerColor(2);
	h2d_el->SetMarkerColor(4);
	h2d_midp->SetMarkerColor(3);




	TCanvas *c0=new TCanvas(Form("c0"),"",900, 600);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	c0->Divide(1,1);
	c0->cd(1);
	h2d_inel->SetTitle(";Proton Track Length/CSDA; #chi^{2} PID [a.u.]");
	h2d_inel->SetMarkerSize(0.2);
	h2d_inel->Draw("");
	c0->Print(Form("./plots_recoinel/ntrklen_chi2_%s_inel.eps", rep.Data()));


	TCanvas *c1=new TCanvas(Form("c1"),"",900, 600);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	c1->Divide(1,1);
	c1->cd(1);
	h2d_el->SetTitle(";Proton Track Length/CSDA; #chi^{2} PID [a.u.]");
	h2d_el->SetMarkerSize(0.2);
	h2d_el->Draw("");
	c1->Print(Form("./plots_recoinel/ntrklen_chi2_%s_el.eps", rep.Data()));	



	TCanvas *c2=new TCanvas(Form("c2"),"",900, 600);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	c2->Divide(1,1);
	c2->cd(1);
	h2d_midp->SetTitle(";Proton Track Length/CSDA; #chi^{2} PID [a.u.]");
	h2d_midp->SetMarkerSize(0.2);
	h2d_midp->Draw("");
	c2->Print(Form("./plots_recoinel/ntrklen_chi2_%s_midp.eps", rep.Data()));	



	TCanvas *c3=new TCanvas(Form("c3"),"",900, 600);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	c3->Divide(1,1);
	c3->cd(1);
	h2d_inel->SetTitle(";Proton Track Length/CSDA; #chi^{2} PID [a.u.]");
	h2d_inel->Draw("");
	h2d_midp->SetMarkerSize(0.4);
	h2d_midp->Draw("same");
	c3->Print(Form("./plots_recoinel/ntrklen_chi2_%s_inel_midp.eps", rep.Data()));	



	TCanvas *c0_=new TCanvas(Form("c0_"),"",900, 600);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	c0_->Divide(1,1);
	c0_->cd(1);
   	c0_->cd(1)->SetLogz();
	//gStyle->SetPalette(kDarkBodyRadiator);
	//gStyle->SetPalette(kRainBow);
	//gStyle->SetPalette(kColorPrintableOnGrey);
	//gStyle->SetPalette(kBlueRedYellow);
	gStyle->SetPalette(kBlackBody); //good
	//gStyle->SetPalette(kLightTemperature);
	//gStyle->SetPalette(kCool);

	
	h2d->SetTitle(";Proton Track Length/CSDA [a.u.]; #chi^{2} PID [a.u.]");
	//h2d->SetMarkerSize(0.2);
	//h2d->RebinX();
	h2d->RebinY(4);
	h2d->Draw("colz");

	//draw cuts
	double mean_norm_trklen_csda=9.01289e-01; //prod4a (spec)
	double sigma_norm_trklen_csda=7.11431e-02; //prod4a (spec)
	double min_norm_trklen_csda=mean_norm_trklen_csda-2.*sigma_norm_trklen_csda;
	double max_norm_trklen_csda=mean_norm_trklen_csda+3.*sigma_norm_trklen_csda;

	//new inel cut using chi^2 info
	double pid_1=7.5;
	double pid_2=10.;

	TLine *l1=new TLine(0,pid_1,min_norm_trklen_csda,pid_1);
	TLine *l2=new TLine(min_norm_trklen_csda,pid_2,max_norm_trklen_csda,pid_2);
	l1->SetLineStyle(2);
	l2->SetLineStyle(2);
	l1->SetLineWidth(2);
	l2->SetLineWidth(2);
	l1->SetLineColor(1);
	l2->SetLineColor(1);
	l1->Draw("same");
	l2->Draw("same");


        //pDUNE Logo
        float logo_y=151;
	float xmin=-0.02;
	float xmax=1.2;
        TLatex **txt_pdune1=new TLatex*[1];
        txt_pdune1[0]=new TLatex(xmin, logo_y, Form("#bf{DUNE:ProtoDUNE-SP}"));
        txt_pdune1[0]->SetTextColor(1);
        txt_pdune1[0]->Draw();
        //
        //Beam Logo
        TLatex **txt_p1=new TLatex*[1];
        //txt_p1[0]=new TLatex(xmax-6.3, logo_y, Form("Protons (1 GeV/c)")); //x:0-20
        txt_p1[0]=new TLatex(xmax-0.35, logo_y, Form("Protons (1 GeV/c)")); //x:0-40
        txt_p1[0]->SetTextColor(1);
        txt_p1[0]->SetTextSize(0.05);
        txt_p1[0]->Draw();

	c0_->Print(Form("./plots_recoinel/ntrklen_chi2_%s_mc_colz.eps", rep.Data()));




}
