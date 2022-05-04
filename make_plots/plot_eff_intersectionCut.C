#include "TVector3.h"
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include "TH2D.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TLine.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH3F.h>
#include <TString.h>
#include <TProfile2D.h>
#include <THStack.h>
#include "TGraph.h" 
#include "TGraphSmooth.h" 
#include "TParameter.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TSystem.h"
#include "string"
#include "vector"
#include "TSpline.h"
#include "TH3F.h"
#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>

#include <stdio.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "../headers/BasicParameters.h"
#include "../headers/BasicFunctions.h"


void plot_KEend_All() {

	TString outpath="./plot_Advanced_KEff/";
        TString fmc="../mc_ke_advancedKEff_nobmrw.root";
	//TString str=Form("h1d_kebb_reco"); //all protons
	//TString str=Form("h1d_kebb_recoEl"); //El
	
	//TString str=Form("h1d_kebb_recoInel"); //Inel
	//TString str_truth=Form("h1d_keend_truth_inel");
	//TString str_title=Form("Reco. Inelastic-scattering Protons; Proton Kinetic Energy[MeV]; Counts");
	//TString str_figout=Form("%sreco_inel.eps",outpath.Data());

	TString str=Form("h1d_kebb_recoEl"); //El
	TString str_truth=Form("h1d_keend_truth_el");
	TString str_title=Form("Reco. Elastic-scattering Protons; Proton Kinetic Energy[MeV]; Counts");
	TString str_figout=Form("%sreco_el.eps",outpath.Data());


	//TString str=Form("h1d_kebb_reco"); //all
	//TString str_truth=Form("h1d_keend_truth");
	//TString str_title=Form("Reco. Incident Protons; Proton Kinetic Energy[MeV]; Counts");
	//TString str_figout=Form("%sreco_all.eps",outpath.Data());



        //plot style --------------------------------------------------------------------//
        gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
        gROOT->SetStyle("protoDUNEStyle");
        gROOT->ForceStyle();
        gStyle->SetTitleX(0.5);
        gStyle->SetTitleAlign(23);
        gStyle->SetOptStat(0);

	//gStyle->SetPalette(53);

	//read mc [after bmrw] ------------------------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fmc.Data());
	TH1D *h1d=(TH1D *)f_mc->Get(Form("%s",str.Data()));
	TH1D *h1d_inel=(TH1D *)f_mc->Get(Form("%s_inel",str.Data()));
	TH1D *h1d_el=(TH1D *)f_mc->Get(Form("%s_el",str.Data()));
	TH1D *h1d_midcosmic=(TH1D *)f_mc->Get(Form("%s_midcosmic",str.Data()));
	TH1D *h1d_midpi=(TH1D *)f_mc->Get(Form("%s_midpi",str.Data()));
	TH1D *h1d_midp=(TH1D *)f_mc->Get(Form("%s_midp",str.Data()));
	TH1D *h1d_midmu=(TH1D *)f_mc->Get(Form("%s_midmu",str.Data()));
	TH1D *h1d_mideg=(TH1D *)f_mc->Get(Form("%s_mideg",str.Data()));
	TH1D *h1d_midother=(TH1D *)f_mc->Get(Form("%s_midother",str.Data()));

	TH1D *h1d_truth=(TH1D *)f_mc->Get(Form("%s",str_truth.Data()));
	h1d_truth->SetLineColor(7);

	h1d_inel->SetFillColor(2); h1d_inel->SetLineColor(2);
	h1d_el->SetFillColor(4); h1d_el->SetLineColor(4);
	h1d_midp->SetFillColor(3); h1d_midp->SetLineColor(3);
	h1d_midcosmic->SetFillColor(5); h1d_midcosmic->SetLineColor(5);
	h1d_midpi->SetFillColor(6); h1d_midpi->SetLineColor(6);
	h1d_midmu->SetFillColor(28); h1d_midmu->SetLineColor(28);
	h1d_mideg->SetFillColor(30); h1d_mideg->SetLineColor(30);
	h1d_midother->SetFillColor(15); h1d_midother->SetLineColor(15);

	THStack* hs=new THStack("hs","");
	hs->Add(h1d_inel);
	hs->Add(h1d_el);
	hs->Add(h1d_midp);
	hs->Add(h1d_midcosmic);
	hs->Add(h1d_midpi);
	hs->Add(h1d_midmu);
	hs->Add(h1d_mideg);
	hs->Add(h1d_midother);

	TCanvas *c_ = new TCanvas("c_", "c_", 1200,800);
	c_->cd(1)->SetLogy();

	//TH2D *f2d=new TH2D("f2d",Form("%s",""),550,-50,500,150,0,1500); 
	//TH2D *f2d=new TH2D("f2d",Form("%s",""),100,-50,50,100,1,500000); 
	TH2D *f2d=new TH2D("f2d",Form("%s",""),100,-50,50,100,0.5,500000); 
	//TH2D *f2d=new TH2D("f2d",Form("%s",""),100,-50,50,100,0,200000); 

	f2d->SetTitle(Form("%s",str_title.Data()));
	f2d->GetYaxis()->SetTitleOffset(1.1);
	f2d->Draw();
	//f2d->GetXaxis()->SetLabelSize(0);
	
	hs->Draw("hist same");
	//h1d_mc->Draw("hist same");
	h1d->Draw("same");
	h1d_truth->Draw("same");

	TLegend *leg = new TLegend(0.2,0.6,0.8,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(h1d_truth, "Truth","l");
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

	c_->Print(Form("%s",str_figout.Data()));



}
