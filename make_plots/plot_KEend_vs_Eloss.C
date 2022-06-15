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

//Gaussian fit
TF1* FitKEEnd(TH1D* h, Int_t col) {
        //pre-fit parameters
        float pre_mean=h->GetBinCenter(h->GetMaximumBin());
        float pre_max=h->GetBinContent(h->GetMaximumBin());
        float pre_rms=50.;
        cout<<"mean: "<<pre_mean<<endl;
        cout<<"rms: "<<pre_rms<<endl;
        cout<<"max: "<<pre_max<<endl;
        cout<<""<<endl;

        //1st fitting ---------------------------------------------------------//
        TF1 *gg=new TF1("gg", fitg, pre_mean-3*pre_rms, pre_mean+3*pre_rms, 3);
        gg->SetParameter(0,pre_mean);
        gg->SetParameter(1,pre_rms);
        gg->SetParameter(2,pre_max);
        //if (pre_rms>1.0e+06) { gg->SetParLimits(1,0,100); }

        //gg->SetLineColor(col);
        //gg->SetLineStyle(2);
        h->Fit("gg","remn");

        //2nd fitting -----------------------------------------------------------------------------------------------------------//
        //TF1 *g=new TF1("g", fitg, gg->GetParameter(0)-1.*gg->GetParameter(1), gg->GetParameter(0)+1.*gg->GetParameter(1), 3);
        TF1 *g=new TF1("g", fitg, gg->GetParameter(0)-40., gg->GetParameter(0)+20., 3);
        //TF1 *g=new TF1("g",fitg,0.3,0.5,3);

        //TF1 *g=new TF1("g",fitg,gg->GetParameter(0)-1,gg->GetParameter(0)+.5,3);
        g->SetParameter(0, gg->GetParameter(0));
        g->SetParameter(1, gg->GetParameter(1));
        g->SetParameter(2, gg->GetParameter(2));

        //g->SetParLimits(0,gg->GetParameter(0)-1*gg->GetParameter(1), gg->GetParameter(0)+1*gg->GetParameter(1));
	//double sss=gg->GetParameter(1); if (sss<0) sss=-sss;
        //g->SetParLimits(1,0,5.*sss);
        //g->SetParLimits(2,0,100.*sqrt(pre_max));

        g->SetLineColor(col);
        g->SetLineStyle(2);
        g->SetLineWidth(2);

       h->Fit("g","remn");
       return g;
}



void plot_KEend_vs_Eloss() {

	TString outpath="./plot_KEend_Eloss_bmrw/";
	//TString fin="../mc_kecalo_bmrw.root";
	//TString fin="../mc_kecalo_nobmrw_ElossTune.root";
	TString fin="../mc_kecalo_bmrw_ElossTune.root";
	
	//setup e-loss map 
	vector<double> Eloss;
	const int n_eloss=40;
        double Eloss_start=37.;
        double dEloss=0.5;
        for (int k=0; k<n_eloss; ++k) {
                double eloss_step=Eloss_start+(double)k*dEloss;
                Eloss.push_back(eloss_step);

                //h1d_KEend_tune_bmrw[k]=new TH1D(Form("h1d_KEend_tune_bmrw_%d",k),"", nke, kemin, kemax);
        }


	//const int n_eloss=2;
        //for (int k=0; k<n_eloss; ++k) {
		//if (k==0) Eloss.push_back(55);
		//if (k==1) Eloss.push_back(40);
	//}
	TH1D *h1d[n_eloss];


        //plot style --------------------------------------------------------------------//
        gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
        gROOT->SetStyle("protoDUNEStyle");
        gROOT->ForceStyle();
        gStyle->SetTitleX(0.5);
        gStyle->SetTitleAlign(23);
        gStyle->SetOptStat(0);

	//gStyle->SetPalette(53);

	//read mc [after bmrw] ------------------------------------------------------------------------------//
	TFile *f_ = TFile::Open(fin.Data());

	TCanvas *c_ = new TCanvas("c_", "c_", 1200,1800);
	c_->Divide(5,8);

	TF1 *fitke[n_eloss];
	vector<double> mu;
	vector<double> err_mu;
	vector<double> ex;
	TH2D* f2d[n_eloss];
        for (int k=0; k<n_eloss; ++k) {
                //double eloss_step=Eloss_start+(double)k*dEloss;
                //Eloss.push_back(eloss_step);

                h1d[k]=(TH1D*)f_->Get(Form("h1d_KEend_tune_bmrw_%d",k));

		//fit -------------------------------------------------------------------------------//
		fitke[k]=FitKEEnd(h1d[k],2);

		c_->cd(k+1);
		f2d[k]=new TH2D(Form("f2d_%d",k),"",200,-100,100, 100, 0,10000);
		f2d[k]->SetTitle(Form("E-loss=%.1f MeV; Proton KE [MeV]; Events",Eloss.at(k)));
		f2d[k]->Draw();
		h1d[k]->Draw("same");
		fitke[k]->Draw("same");

	
		//get mean & err_mean from the fit ------------//
		double mu_fit=fitke[k]->GetParameter(0);
		double err_mu_fit=fitke[k]->GetParError(0);
		mu.push_back(mu_fit);
		err_mu.push_back(err_mu_fit);
		ex.push_back(0);
        }
	c_->Print(Form("%sKEend_fit_MC.eps",outpath.Data()));
	
        gStyle->SetOptFit(0);
	TCanvas *c_2 = new TCanvas("c_2", "c_2", 1200,800);
	c_2->Divide(1,1);
	TGraphErrors *gr_de_mu = new TGraphErrors(mu.size(), &Eloss.at(0), &mu.at(0), &ex.at(0), &err_mu.at(0));
	//gr_de_mu->SetTitle("MC");
	//gr_de_mu->GetXaxis()->SetTitle("<#DeltaE> [MeV]");
	//gr_de_mu->GetYaxis()->SetTitle("#mu of KE_{end} [MeV]");
	TH2D* f2d_fit=new TH2D("f2d_fit","", 25, 35, 60, 25,-10,15);
	//f2d_fit->SetTitle("MC; <#DeltaE> [MeV]; #mu of KE_{end} [MeV]");	
	f2d_fit->SetTitle("MC (bmrw); <#DeltaE> [MeV]; #mu of KE_{end} [MeV]");	
	f2d_fit->Draw();

	gr_de_mu->Draw("p same");

	double fit_min=37.;
	double fit_max=58.;
        TF1 *line_best=new TF1("line_best", "[0]+[1]*x", fit_min, fit_max);
	double mm=(mu.at(0)-mu.at(-1+mu.size()))/(Eloss.at(0)-Eloss.at(-1+Eloss.size()));
	double bb=mu.at(0)-mm*Eloss.at(0);
        line_best->SetParameters(bb,mm);
        //line_best->SetParameter(1, mm);

        line_best->SetLineColor(2);
        line_best->SetLineStyle(1);
        line_best->SetLineWidth(2);

       	gr_de_mu->Fit("line_best","rem+");

	double p[2];
	double err_p[2];
	for (int i=0; i<2; ++i) {
		p[i]=line_best->GetParameter(i);
		err_p[i]=line_best->GetParError(i);

		cout<<"p["<<i<<"]:"<<p[i]<<endl;
		cout<<"err_p["<<i<<"]:"<<err_p[i]<<endl;
	}

	TLine *l0=new TLine (fit_min,0,fit_max,0);
	l0->SetLineColor(3);
        l0->SetLineStyle(3);
        l0->SetLineWidth(3);
	l0->Draw("same");

	TLegend *leg = new TLegend(0.4,0.6,0.85,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry((TObject*)0, Form("p0:%.2f#pm%.2f",p[0],err_p[0]),"");
	leg->AddEntry((TObject*)0, Form("p1:%.2f#pm%.2f",p[1],err_p[1]),"");
	leg->AddEntry((TObject*)0, Form("#chi^{2}/NDF=%.2f/%d",line_best->GetChisquare(), line_best->GetNDF()),"");
	leg->AddEntry((TObject*)0, Form("Extracted <#DeltaE>=%.2f MeV",-p[0]/p[1]),"");

	leg->Draw();


	c_2->Print(Form("%sKEend_mu_dE_MC.eps",outpath.Data()));

	


/*
	TCanvas *c_ = new TCanvas("c_", "c_", 1200,800);
	c_->cd(1)->SetLogy();
	TH2D *f2d=new TH2D("f2d",Form("%s",""),100,-50,50,100,0.5,500000); 
	f2d->SetTitle(Form("%s",str_title.Data()));
	f2d->GetYaxis()->SetTitleOffset(1.1);
	f2d->Draw();

	c_->Print(Form("%s",str_figout.Data()));
*/


}
