#include "TVector3.h"
#include <TH2.h>
#include <TH1.h>
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

TGraphErrors* egr(TGraphErrors*in, int n_st) {

	Double_t *x0 = in->GetX();
   	Double_t *y0 = in->GetY();
	Double_t *ex0 = in->GetEX();
   	Double_t *ey0 = in->GetEY();

	vector<double> x;
	vector<double> y;
	vector<double> ex;
	vector<double> ey;

  	for(int i=0; i<in->GetN()-n_st; i++) {
		x.push_back(x0[i]);
		y.push_back(y0[i]);
		ex.push_back(ex0[i]);
		ey.push_back(ey0[i]);
	}
	cout<<"input size="<<in->GetN()<<endl;
	cout<<"output size="<<x.size()<<endl;
	TGraphErrors* egr_out=new TGraphErrors(x.size(), &x.at(0), &y.at(0), &ex.at(0), &ey.at(0));

	return egr_out;
}



void makeXS_Compare_Eslice_Thinslice() {

	TString outpath="./plots_XS_comp/";

	//[0]root -l make_ThinSliceEdataXS.C to produce the root file that contains gr_recoxs & gr_truexs

	//Load XSs -----------------------------------------------------------------------------//
        //TString str_fmc_eslice="./xs_rootfiles/xs_Eslice_dE20MeV_40slcs_nobmrw_stslcplus0.5.root";
        TString str_fmc_eslice="./xs_rootfiles/xs_Eslice_dE20MeV_40slcs_nobmrw_stslcplus0.5_new.root";
        TString str_fmc_sslice="./xs_rootfiles/xs_thinslice_dx5cm_20slcs_nobmrw.root";
	TFile *fmc_eslice = TFile::Open(str_fmc_eslice.Data());
	TFile *fmc_sslice = TFile::Open(str_fmc_sslice.Data());

	//TGraphErrors* tmp_recoxs_eslice=(TGraphErrors*)fmc_eslice->Get("gr_recoxs");
	//TGraphErrors* tmp_truexs_eslice=(TGraphErrors*)fmc_eslice->Get("gr_truexs");
	TGraphErrors* recoxs_eslice=egr((TGraphErrors*)fmc_eslice->Get("gr_recoxs"),1); //skip 1st point
	TGraphErrors* truexs_eslice=egr((TGraphErrors*)fmc_eslice->Get("gr_truexs"),1); //skip 1st point


	TGraphErrors* recoxs_sslice=(TGraphErrors*)fmc_sslice->Get("gr_recoxs");
	TGraphErrors* truexs_sslice=(TGraphErrors*)fmc_sslice->Get("gr_truexs");

        TFile f_xs("/dune/data2/users/hyliao/GeantReweight/xs_cascade/proton_cross_section.root");
        TGraph *total_inel_KE = (TGraph*)f_xs.Get("inel_KE");


        //plot style --------------------------------------------------------------------//
        gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
        gROOT->SetStyle("protoDUNEStyle");
        gROOT->ForceStyle();
        gStyle->SetTitleX(0.5);
        gStyle->SetTitleAlign(23);
        gStyle->SetOptStat(0);


        //gStyle->SetPadLeftMargin(0.13);
        gStyle->SetPadRightMargin(0.0);
 	gStyle->SetPadLeftMargin(0.16);
    	gStyle->SetPadRightMargin(0.04);//0.02



        TCanvas *c_xs = new TCanvas("c_xs", "c_xs", 1400, 900);
	c_xs->Divide(1,1);
	c_xs->cd(1);

        float ymin=200;
        float ymax=1200;
        //float xmax=620;
        float xmax=420;
        float xmin=0;

        TH2D *f2d_xs=new TH2D("f2d_xs","",450,xmin,xmax,ymax,ymin,ymax);
        f2d_xs->GetXaxis()->SetTitle("Proton Kinetic Energy [MeV]");
        f2d_xs->GetYaxis()->SetTitleOffset(1.3);
        f2d_xs->GetYaxis()->SetTitle("P-Ar inelastic cross section [mb]");
        f2d_xs->Draw("");

        total_inel_KE->SetLineColor(2);
        total_inel_KE->SetLineWidth(3);

        total_inel_KE->Draw("c same");


        //truexs_eslice->SetLineWidth(2);
        truexs_eslice->SetMarkerColor(3);
        truexs_eslice->SetLineColor(3);
        truexs_eslice->SetMarkerStyle(20);
        truexs_eslice->SetMarkerSize(1.5);


        //recoxs_eslice->SetLineWidth(2);
        recoxs_eslice->SetMarkerColor(1);
        recoxs_eslice->SetLineColor(1);
        recoxs_eslice->SetMarkerStyle(20);
        recoxs_eslice->SetMarkerSize(1.5);


        truexs_sslice->SetLineWidth(1);
        truexs_sslice->SetMarkerColor(94);
        truexs_sslice->SetLineColor(94);
        truexs_sslice->SetMarkerStyle(24);
        truexs_sslice->SetMarkerSize(1.5);

        recoxs_sslice->SetLineWidth(1);
        recoxs_sslice->SetMarkerColor(4);
        recoxs_sslice->SetLineColor(4);
        recoxs_sslice->SetMarkerStyle(24);
        recoxs_sslice->SetMarkerSize(1.5);


        truexs_sslice->Draw("p same");
        recoxs_sslice->Draw("p same");
        truexs_eslice->Draw("pz same");
        recoxs_eslice->Draw("pz same");


        TLegend *leg_xs = new TLegend(0.3,0.6,0.95,0.85);
        leg_xs->SetFillStyle(0);
	leg_xs->SetNColumns(3);
        leg_xs->AddEntry(total_inel_KE, "Geant4", "l");
        leg_xs->AddEntry((TObject*)0, "", "");
        leg_xs->AddEntry((TObject*)0, "", "");
        leg_xs->AddEntry(truexs_eslice, "MC Truth (E-slice)", "pe");
        leg_xs->AddEntry(recoxs_eslice, "MC Reco (E-slice)", "pe");
        leg_xs->AddEntry((TObject*)0, "", "");
        leg_xs->AddEntry(truexs_sslice, "MC Truth (Thin-slice)", "pe");
        leg_xs->AddEntry(recoxs_sslice, "MC Reco (Thin-slice)", "pe");
        leg_xs->AddEntry((TObject*)0, "", "");
        leg_xs->Draw();


/*
        TLegend *leg_xs = new TLegend(0.3,0.6,0.95,0.85);
        leg_xs->SetFillStyle(0);
        leg_xs->AddEntry(total_inel_KE, "Geant4", "l");
        leg_xs->AddEntry(truexs_eslice, "MC Truth (E-slice)", "pe");
        leg_xs->AddEntry(recoxs_eslice, "MC Reco (E-slice)", "pe");
        leg_xs->Draw();
*/

        //pDUNE Logo
        float logo_y=ymax+9;
        TLatex **txt_pdune1=new TLatex*[1];
        txt_pdune1[0]=new TLatex(xmin, logo_y, Form("#bf{DUNE:ProtoDUNE-SP}"));
        txt_pdune1[0]->SetTextColor(1);
        txt_pdune1[0]->Draw();
        //
        //Beam Logo
        TLatex **txt_p1=new TLatex*[1];
        txt_p1[0]=new TLatex(xmax-142, logo_y, Form("Proton MC (1 GeV/c)"));
        txt_p1[0]->SetTextColor(1);
        txt_p1[0]->SetTextSize(0.05);
        txt_p1[0]->Draw();


	//c_xs->Print(Form("%sxs_data.eps",outpath.Data()));
	c_xs->Print(Form("%spArxs_eslie_thineslice.eps",outpath.Data()));
	c_xs->Print(Form("%spArxs_eslie_thineslice.pdf",outpath.Data()));
	//c_xs->Print(Form("%spArxs_eslie.pdf",outpath.Data()));
	//c_xs->Print(Form("%spArxs_eslie.eps",outpath.Data()));




/*

	TGraphErrors* recoxs_eslice=(TGraphErrors*)fmc_eslice->Get("gr_recoxs");
	TGraphErrors* truexs_eslice=



        truexs_eslice->SetLineWidth(2);
        truexs_eslice->SetMarkerColor(3);
        truexs_eslice->SetLineColor(3);
        truexs_eslice->SetMarkerStyle(21);
        truexs_eslice->SetMarkerSize(1.5);
        truexs_eslice->Draw("p same");

        gr_recoxs->SetLineWidth(2);
        gr_recoxs->SetMarkerColor(1);
        gr_recoxs->SetLineColor(1);
        gr_recoxs->SetMarkerStyle(20);
        gr_recoxs->SetMarkerSize(1.5);
        //gr_recoxs->Draw("p same");
*/



/*

	double ke_ff_mean_stop=4.12972e+02;
	double ke_ff_sigma_stop=4.08988e+01;
	TLine l_mean(ke_ff_mean_stop,0,ke_ff_mean_stop,ymax);
	TLine l_left(ke_ff_mean_stop-3.*ke_ff_sigma_stop,0,ke_ff_mean_stop-3.*ke_ff_sigma_stop,ymax);
	TLine l_right(ke_ff_mean_stop+3.*ke_ff_sigma_stop,0,ke_ff_mean_stop+3.*ke_ff_sigma_stop,ymax);
	l_mean.SetLineColor(1);
	l_left.SetLineColor(1);
        l_right.SetLineColor(1);

	l_left.SetLineStyle(2);	
	l_right.SetLineStyle(2);	

	//l_mean.Draw("same");
	//l_left.Draw("same");
	//l_right.Draw("same");

        TLegend *leg_xs = new TLegend(0.6,0.6,0.9,0.85);
        leg_xs->SetFillStyle(0);
        leg_xs->AddEntry(total_inel_KE, "Geant4", "l");
        leg_xs->AddEntry(gr_truexs, "MC Truth", "pe");
        //leg_xs->AddEntry(gr_recoxs, "MC Reco", "pe");
        leg_xs->Draw();


	//c_xs->Print(Form("%sxs_data.eps",outpath.Data()));
	c_xs->Print(Form("%sxs_data_noreco.eps",outpath.Data()));



        TFile f_xs("/dune/data2/users/hyliao/GeantReweight/xs_cascade/proton_cross_section.root");
        TGraph *total_inel_KE = (TGraph*)f_xs.Get("inel_KE");

        //gStyle->SetPadLeftMargin(0.13);
        gStyle->SetPadRightMargin(0.0);
 	gStyle->SetPadLeftMargin(0.16);
    	gStyle->SetPadRightMargin(0.04);//0.02
        TCanvas *c_xs = new TCanvas("c_xs", "c_xs", 1400, 900);
	c_xs->Divide(1,1);
	c_xs->cd(1);

        float ymax=1200;
        float xmax=620;
        float xmin=0;

        TH2D *f2d_xs=new TH2D("f2d_xs","",450,xmin,xmax,ymax,0,ymax);
        f2d_xs->GetXaxis()->SetTitle("Proton Kinetic Energy [MeV]");
        f2d_xs->GetYaxis()->SetTitleOffset(1.3);
        f2d_xs->GetYaxis()->SetTitle("P-Ar inelastic cross section [mb]");
        f2d_xs->Draw("");

        total_inel_KE->SetLineColor(2);
        total_inel_KE->Draw("c same");
        gr_truexs->SetLineWidth(2);
        gr_truexs->SetMarkerColor(3);
        gr_truexs->SetLineColor(3);
        gr_truexs->SetMarkerStyle(21);
        gr_truexs->SetMarkerSize(1.5);
        gr_truexs->Draw("p same");
*/



}
