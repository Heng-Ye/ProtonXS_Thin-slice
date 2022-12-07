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
#include <utility>      // std::pair, std::make_pair
#include <algorithm>

#include "../headers/BasicParameters.h"
#include "../headers/BasicFunctions.h"

bool compare(std::pair<int, double> p1, std::pair<int, double> p2) {
    return p1.second<p2.second;
}

void plot_Sys_BeamMomentumReweighting() {
        //plot style --------------------------------------------------------------------//
        gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
        gROOT->SetStyle("protoDUNEStyle");
        gROOT->ForceStyle();
        gStyle->SetTitleX(0.5);
        gStyle->SetTitleAlign(23);
        gStyle->SetOptStat(0);

	//Data -----------------------------------------------------------------------------------------------------------//
	double mu_denom_data=411.05145837595467; //new with event-by-event R corr (const E-loss using stopping protons)
	double sg_denom_data=47.48714821962207; //new with event-by-event R corr (const E-loss using stopping protons)

	double err_mu_denom_data=0.4;
	double err_sg_denom_data=0.3;

	//MC ------------------------------------------------------------------------------------------------------//
	double mu_nom_mc=416.224743039812; //for mc(KEbeam-const) with R=1 (R=Ratio of KEff(Fit)/(KEbeam-dE))
	double sg_nom_mc=42.786018962508784; //

	//double mu_nom_mc=416.1620092367158; //for mc [KE(Fit)]
	//double sg_nom_mc=40.48356757740762; //

	//i= 0  m= 416.2247430398121 s= 42.786018962508784 [kebeam-dE]*1 (R=1)
	//i= 1  m= 416.1620092367158 s= 40.48356757740762 [KE(Fit)] 

	double err_mu_nom_mc=0.3;
	double err_sg_nom_mc=0.2;

	double mu_kemin=mu_nom_mc-5.*sg_nom_mc;
	double mu_kemax=mu_nom_mc+5.*sg_nom_mc;

	//list of sys
	vector< pair <double,double> > denom; //mu, sigma 
	vector< pair <double,double> > nom; //mu, sigma

	int cnt1=0;
	for (int i1=0; i1<=2; ++i1) {
		double sign=0;
		if (i1==0) sign=0;
		if (i1==1) sign=1;
		if (i1==2) sign=-1;
		
		double mu1=mu_denom_data+sign*err_mu_denom_data;
		double mu2=mu_nom_mc+sign*err_mu_nom_mc;

		for (int i2=0; i2<=2; ++i2) {
			double sign=0;
			if (i2==0) sign=0;
			if (i2==1) sign=1;
			if (i2==2) sign=-1;
	
			double sg1=sg_denom_data+err_sg_denom_data;
			double sg2=sg_nom_mc+err_sg_nom_mc;

			denom.push_back(make_pair(mu1,sg1));
			nom.push_back(make_pair(mu2,sg2));
			if (i1==0&&i2==0) std::cout<<"index of (0,)="<<cnt1<<std::endl;
			cnt1++;
		}	
	}

	cout<<"denom.size()="<<denom.size()<<endl;
	const int n_tot=(const int)pow(denom.size(),2);
	TF1 **kerw=new TF1*[n_tot];

	int cnt_=0;
	TCanvas *c_ = new TCanvas("c_", "c_", 1200, 900);
	c_->cd(1);
	double plot_ymin=0.9;
	double plot_ymax=2;
	//TH2D *f2d=new TH2D("f2d",Form("f2d"), 550, 150, 700, 10, ymin, ymax);
	TH2D *f2d=new TH2D("f2d",Form("f2d"), 550, 150, 700, 10, plot_ymin, plot_ymax);
	f2d->SetTitle("Beam Momentum Reweighting Systematics; Proton Kinetic Energy[MeV]; Weighting Function");
	f2d->Draw();
	vector< pair <int,double> > sys_limit;
	double val_ref=440;  //referenceto determine upper and lower curves	
	for (int j=0; j<denom.size(); ++j) {
		//cout <<"denom:"<< denom[j].first << " " << denom[j].second<<"| nom:"<< nom[j].first << " " << nom[j].second<< endl;

		for (int k=0; k<nom.size(); ++k) {
			kerw[cnt_]=new TF1(Form("kerw_%d",cnt_),govg, 0, 800, 4);	
			kerw[cnt_]->SetParameter(0, nom[k].first);
			kerw[cnt_]->SetParameter(1, nom[k].second);
			kerw[cnt_]->SetParameter(2, denom[j].first);
			kerw[cnt_]->SetParameter(3, denom[j].second);
			
			kerw[cnt_]->SetLineColor(3);
			kerw[cnt_]->Draw("same");

			sys_limit.push_back(make_pair(cnt_, kerw[cnt_]->Eval(val_ref)));
			
			cnt_++;
		}
					
	}
	kerw[0]->SetLineColor(2);
	kerw[0]->Draw("same");

	TLine **ver_lines=new TLine*[2]; 
	for (int k=0; k<2; ++k) {
        	if (k==0) ver_lines[k]=new TLine(mu_kemin, plot_ymin,mu_kemin,plot_ymax);
        	if (k==1) ver_lines[k]=new TLine(mu_kemax, plot_ymin,mu_kemax,plot_ymax);
		ver_lines[k]->SetLineStyle(2);
		ver_lines[k]->Draw("same");
	}
	
	TLegend *leg = new TLegend(0.4,0.7,0.6,0.8);
	leg->SetFillStyle(1);
	leg->AddEntry(kerw[0], "Central value","l");
	leg->AddEntry(kerw[1], "#pm1-#sigma","l");
	leg->Draw();

	c_->Print("./plot_bmrw_sys/bmrw_sys.eps");


	//one-by-one
	TCanvas *c_2 = new TCanvas("c_2", "c_2", 1200, 900);
	c_2->cd(1);
	double plot_ymin2=0.95;
	double plot_ymax2=1.10;
	double plot_xmin=390;
	double plot_xmax=470;
	TH2D *f2d2=new TH2D("f2d2",Form("f2d"), 90, plot_xmin, plot_xmax, 10, plot_ymin2, plot_ymax2);
	f2d2->SetTitle("Beam Momentum Reweighting Systematics; Proton Kinetic Energy[MeV]; Weighting Function");
	f2d2->Draw();

	//find maximum/minimum curve with given reference value -----------------------------------
	const auto max_p = max_element(sys_limit.begin(), sys_limit.end(), compare);
	const auto min_p = min_element(sys_limit.begin(), sys_limit.end(), compare);

    	int key_max = max_p->first;
    	double max_val = max_p->second;
    	int key_min = min_p->first;
    	double min_val = min_p->second;

	std::cout<<"max_val(ref_val="<<val_ref<<" MeV)="<<max_val<<" | key="<<key_max<<std::endl;
	std::cout<<"min_val(ref_val="<<val_ref<<" MeV)="<<min_val<<" | key="<<key_min<<std::endl;
	//-------------------------------------------------------------------------------------

	for (int k=0; k<cnt_; ++k) {
		kerw[0]->Draw("same");

		if (k==key_max) kerw[k]->SetLineColor(6);
		if (k==key_min) kerw[k]->SetLineColor(7);
		kerw[k]->Draw("same");

		//.push_back(make_pair(mu2,sg2));
		//	double tmp_up=kerw[cnt_]->Eval(440);
		//	double tmp_dn=kerw[cnt_]->Eval(440);

		//TLegend *leg_ = new TLegend(0.2,0.9-k*0.05,0.6,0.95-k*0.05);
		//leg_->SetFillStyle(0);
		//leg_->AddEntry((TObject*)0, Form("Index:%d",k),"");
		//leg_->AddEntry(kerw[0], "Central value","l");
		//leg_->AddEntry(kerw[1], Form("#mu1:%.2f#pm%.2f"),"l");
		//leg_->Draw();

		if (k==0) c_2->Print("./plot_bmrw_sys/bmrw_sys_indiv.pdf(");
		//else if (k>0&&k<cnt_-1) c_2->Print("./plot_bmrw_sys/bmrw_sys_indiv.pdf");
		else c_2->Print("./plot_bmrw_sys/bmrw_sys_indiv.pdf");
					
	}
	kerw[key_max]->Draw("same");
	kerw[key_min]->Draw("same");
        c_2->Print("./plot_bmrw_sys/bmrw_sys_indiv.pdf)");

	//Summary plot of beam mom. systematics --------------//
	TCanvas *c_3 = new TCanvas("c_3", "c_3", 1200, 900);
	//c_3->SetGrid(0);
   	//c_3->DrawFrame(0,0,2.2,12);

	c_3->cd(1);
	f2d2->Draw();

   	const Int_t n = 100;
	double dE=(plot_xmax-plot_xmin)/(double)n;
   	Double_t x[n], y[n],ymin[n], ymax[n];
   	Int_t i;
   	for (i=0;i<n;i++) {
		x[i]=plot_xmin+dE*(double)i;
		ymax[i]=kerw[key_max]->Eval(x[i]);
		ymin[i]=kerw[key_min]->Eval(x[i]);	
	}
	TGraph *grmin = new TGraph(n,x,ymin);
   	TGraph *grmax = new TGraph(n,x,ymax);

   	TGraph *grshade = new TGraph(2*n);
   	for (i=0;i<n;i++) {
      		grshade->SetPoint(i,x[i],ymax[i]);
      		grshade->SetPoint(n+i,x[n-i-1],ymin[n-i-1]);
   	}
   	grshade->SetFillColorAlpha(3, 0.2);
   	grshade->Draw("f");
   	//grmin->Draw("l");
   	//grmax->Draw("l");

	kerw[0]->Draw("same");
	kerw[key_max]->SetLineColor(1);
	kerw[key_min]->SetLineColor(1);
	//kerw[key_max]->SetLineStyle(3);
	//kerw[key_min]->SetLineStyle(3);
	kerw[key_max]->Draw("same");
	kerw[key_min]->Draw("same");

	TLegend *leg_3 = new TLegend(0.35,0.7,0.65,0.8);
	leg_3->SetFillStyle(0);
	leg_3->AddEntry(kerw[0], "Central value","l");
	leg_3->AddEntry(kerw[key_max], "#pm1-#sigma","l");
	leg_3->Draw();

	c_3->Print("./plot_bmrw_sys/bmrw_sys_summary.eps");




}
