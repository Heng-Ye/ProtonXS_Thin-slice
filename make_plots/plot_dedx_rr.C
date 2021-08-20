#include "THStack.h"
#include "../headers/BasicParameters.h"
#include "../headers/BasicFunctions.h"

void plot_dedx_rr() {

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);

	TFile *f_mc = TFile::Open("../prod4a_thinslice_dx4cm_25slcs.root");
	TH2D* rr_dedx_recostop=(TH2D*)f_mc->Get("rr_dedx_recostop");
	//TGraph *gr_predict_dedx_resrange=(TGraph *)f_mc->Get("gr_predict_dedx_resrange");
	gr_predict_dedx_resrange->SetLineColor(2);
	//gr_predict_dedx_resrange->SetLineWidth(1);
	gr_predict_dedx_resrange->SetMarkerColor(2);
	gr_predict_dedx_resrange->SetMarkerSize(.2);

	TCanvas *c1=new TCanvas(Form("c1"),"",900, 600);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	c1->Divide(1,1);
	c1->cd(1);
	rr_dedx_recostop->SetTitle("Stopping protons; Residual Range[cm]; dE/dx [MeV/cm]");	
	rr_dedx_recostop->Draw("colz");
	gr_predict_dedx_resrange->Draw("cp same");

	TLegend *leg=new TLegend(0.522519,0.758562,0.808461,0.86135);
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	TLegendEntry* ll[2];
	ll[0]=leg->AddEntry(gr_predict_dedx_resrange,"Expectation","l"); ll[0]->SetTextColor(1); ll[0]->SetLineColor(2);			
	leg->Draw();
	c1->Print("./plots_recoinel/dedx_rr_mc_recostop.eps");	


}
