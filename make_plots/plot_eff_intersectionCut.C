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


void plot_eff_intersectionCut() {

        //TString fmc="../mc_kebbbkg_nobmrw_HD.root";
        TString fmc="../mc_kebkg_nobmrw.root";
        TString fdata="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_kebkg_withintersection.root";


	//TString str=Form("trklen_RecoEl"); //El
	//TString str=Form("trklen_RecoInEl"); //InEl
	TString str=Form("trklen_All"); //All
	TString outpath=Form("./plot_eff_intersection/%s.eps",str.Data());
	TString outpath2=Form("./plot_eff_intersection/%s_ba.eps",str.Data());


        //plot style -------------------------------------------------//
        gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
        gROOT->SetStyle("protoDUNEStyle");
        gROOT->ForceStyle();
        gStyle->SetTitleX(0.5);
        gStyle->SetTitleAlign(23);
        gStyle->SetOptStat(0);
	//-------------------------------------------------------------//

	//gStyle->SetPalette(53);

	//read mc [after bmrw] ------------------------------------------------------------------------------//
	//mc
	TFile *f_mc = TFile::Open(fmc.Data());
	TH1D *h1d_b=(TH1D *)f_mc->Get(Form("%s",str.Data()));
	TH1D *h1d_a=(TH1D *)f_mc->Get(Form("%s_intersection",str.Data()));
	TH1D *R=(TH1D *)h1d_a->Clone();
	R->Divide(h1d_b);

	//data
	TFile *f_data = TFile::Open(fdata.Data());
	TH1D *h1d_b_data=(TH1D *)f_data->Get(Form("%s",str.Data()));
	TH1D *h1d_a_data=(TH1D *)f_data->Get(Form("%s_intersection",str.Data()));
	TH1D *R_data=(TH1D *)h1d_a_data->Clone();
	R_data->Divide(h1d_b_data);


	TCanvas *c_ = new TCanvas("c_", "c_", 1200,800);
	//c_->cd(1)->SetLogy();
	R->SetTitle("; Reco track length [cm]; Survival Fraction");
	R->SetLineColor(2);
	R->SetMarkerColor(2);
	R->Draw("ep");
	R_data->Draw("ep same");

TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
leg->SetFillStyle(0);
leg->AddEntry(R_data, "Data","ep");
leg->AddEntry(R, "MC","ep");
leg->Draw();

	c_->Print(Form("%s",outpath.Data()));



	TCanvas *c_2 = new TCanvas("c_2", "c_2", 1200,800);
	//c_->cd(1)->SetLogy();
	h1d_b->Draw();
	h1d_a->SetLineColor(2);
	h1d_a->Draw("same");
	c_2->Print(Form("%s",outpath2.Data()));






/*
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

*/


}
