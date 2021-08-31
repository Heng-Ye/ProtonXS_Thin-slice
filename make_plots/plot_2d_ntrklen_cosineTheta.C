#include "THStack.h"
#include "../headers/BasicParameters.h"
#include "../headers/BasicFunctions.h"

void plot_2d_ntrklen_cosineTheta() {
	TString rep="Pos";
	TString title=";Proton Track Length/CSDA range[a.u.]; cos#theta[a.u.]";	
	TString fig_out="ntrklen_cosTheta_";

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);

	TFile *f_mc = TFile::Open("../mc_proton_studyRecoInelCut.root");
	//TH2D* h2d=(TH2D*)f_mc->Get(Form("ntrklen_chi2pid_%s", rep.Data()));
	TH2D* h2d_inel=(TH2D*)f_mc->Get(Form("ntrklen_cosineTheta_%s_inel", rep.Data()));
	TH2D* h2d_el=(TH2D*)f_mc->Get(Form("ntrklen_cosineTheta_%s_el", rep.Data()));
	TH2D* h2d_midp=(TH2D*)f_mc->Get(Form("ntrklen_cosineTheta_%s_midp", rep.Data()));

	h2d_inel->SetMarkerColor(2);
	h2d_el->SetMarkerColor(4);
	h2d_midp->SetMarkerColor(3);


	TCanvas *c0=new TCanvas(Form("c0"),"",900, 600);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	c0->Divide(1,1);
	c0->cd(1);
	h2d_inel->SetTitle(title.Data());
	h2d_inel->SetMarkerSize(0.2);
	h2d_inel->Draw("");
	c0->Print(Form("./plots_recoinel/%s_%s_inel.eps", fig_out.Data(), rep.Data()));


	TCanvas *c1=new TCanvas(Form("c1"),"",900, 600);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	c1->Divide(1,1);
	c1->cd(1);
	h2d_el->SetTitle(title.Data());
	h2d_el->SetMarkerSize(0.2);
	h2d_el->Draw("");
	c1->Print(Form("./plots_recoinel/%s_%s_el.eps", fig_out.Data(), rep.Data()));	



	TCanvas *c2=new TCanvas(Form("c2"),"",900, 600);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	c2->Divide(1,1);
	c2->cd(1);
	h2d_midp->SetTitle(title.Data());
	h2d_midp->SetMarkerSize(0.2);
	h2d_midp->Draw("");
	c2->Print(Form("./plots_recoinel/%s_%s_midp.eps", fig_out.Data(), rep.Data()));	








	TCanvas *c3=new TCanvas(Form("c3"),"",900, 600);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	c3->Divide(1,1);
	c3->cd(1);
	h2d_inel->SetTitle(";Proton Track Length/CSDA; #chi^{2} PID [a.u.]");
	h2d_inel->Draw("");
	h2d_midp->SetMarkerSize(0.4);
	h2d_midp->Draw("same");
	c3->Print(Form("./plots_recoinel/%s_%s_inel_midp.eps", fig_out.Data(), rep.Data()));	





}
