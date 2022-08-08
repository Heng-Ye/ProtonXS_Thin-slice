#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include <utility>      // std::pair, std::make_pair
#include <string>       // std::string
#include <iostream>     // std::cout
#include "TString.h"

#include "TH1.h"
#include "TRandom2.h"
#include "TFile.h"

//#include "TF1.h"
//#include "./headers/BasicParameters.h"
//#include "./headers/BasicFunctions.h"
#include "./headers/BetheBloch.h"

using namespace std;

// number of distribution generated points
#define NGEN 100000

Double_t fitg(Double_t* x,Double_t *par) {
	double m=par[0];
	double s=par[1];
	double n=par[2];

	double g=n*TMath::Exp(-(x[0]-m)*(x[0]-m)/(2*s*s));
	//Double_t g=n/(s*sqrt(2*3.14159))*TMath::Exp(-(x[0]-m)*(x[0]-m)/(2*s*s));
	return g;
}

TF1* VNFit(TH1F* h, float pre_mean, float n_sigma) {
	//pre-fit parameters
	//float pre_mean=h->GetBinCenter(h->GetMaximumBin());
	float pre_max=h->GetMaximum();
	float pre_rms=h->GetRMS();
	cout<<"mean: "<<pre_mean<<endl;
	cout<<"rms: "<<pre_rms<<endl;
	cout<<"max: "<<pre_max<<endl;
	cout<<""<<endl;

	//1st fitting
	TF1 *gg=new TF1("gg",fitg,pre_mean-n_sigma*pre_rms,pre_mean+n_sigma*pre_rms,3);
	gg->SetParameter(0,pre_mean);
	gg->SetParameter(1,pre_rms);
	gg->SetParameter(2,pre_max);
	//if (pre_rms>1.0e+06) { gg->SetParLimits(1,0,100); }

	//gg->SetLineColor(col);
	gg->SetLineStyle(2);
	h->Fit("gg","remn");

	//2nd fitting
	TF1 *g=new TF1("g",fitg,gg->GetParameter(0)-n_sigma*gg->GetParameter(1),gg->GetParameter(0)+n_sigma*gg->GetParameter(1),3);
	//TF1 *g=new TF1("g",fitg,gg->GetParameter(0)-1,gg->GetParameter(0)+.5,3);
	g->SetParameter(0,gg->GetParameter(0));
	g->SetParameter(1,gg->GetParameter(1));
	g->SetParameter(2,gg->GetParameter(2));

	//g->SetParLimits(0,gg->GetParameter(0)-3*gg->GetParameter(1), gg->GetParameter(0)+3*gg->GetParameter(1));
	//double sss=gg->GetParameter(1); if (sss<0) sss=-sss;
	//g->SetParLimits(1,0,5.*sss);
	//g->SetParLimits(2,0,10.*sqrt(pre_max));

	//g->SetLineColor(col);
	g->SetLineStyle(2);
	g->SetLineWidth(2);

	h->Fit("g","remn");
	return g;
}


int main(int argc, char* argv[]) {

        //const. E-loss assumption --------------------------------------------------//     
        double const_eloss_mc=47.0058/1.00097; //const E-loss from fit (calo)
        //#p[0]:47.0058 err_p[0]:0.372157 p[1]:-1.00097 err_p[1]:0.00787403

        //New weighting func (using KEff_fit_stop at TPC FF as a reference) ------//
        double mu_denom_data=411.06602388610895;
        double sg_denom_data=47.075678784947826;

        double mu_nom_data=390.81237292943916; //for data
        double sg_nom_data=47.52091718691363; //for data

        double mu_nom_mc=388.560260293186; //for mc(KEbeam-const)
        double sg_nom_mc=43.13168235197187; //formc

        //weighting func. (ke)
        //TF1 *kerw=new TF1(Form("kerw"),govg,0,800,4);
        //kerw->SetParameter(0, mu_nom_mc);
        //kerw->SetParameter(1, sg_nom_mc);
        //kerw->SetParameter(2, mu_denom_data);
        //kerw->SetParameter(3, sg_denom_data);

        //Basic configure ------//
        BetheBloch BB;
        BB.SetPdgCode(2212);
        //----------------------//      

        //double kebb=-50; kebb=BB.KEAtLength(ke_ffbeam_MeV, range_reco);

	//double ke_fit=BB.KEAtLength(
   	int n=NGEN; //number of points in generator
	double nke=1000;
	double ke_min=0;
	double ke_max=1000;


	const int n_range=1;
	TH1F *gaus_kefit_data[n_range];
	TF1* Fit_gaus_kefit_data[n_range];

	//TH1F * gaus_kefit = new TH1F("gaus_kefit","", nke, ke_min, ke_max); 
	TFile *fout = new TFile("test.root","RECREATE");


	for (int j=0; j<n_range; ++j) {
		gaus_kefit_data[j]=new TH1F(Form("gaus_kefit_data_%d",j),Form("Range:%d [cm]",j), nke, ke_min, ke_max);
		#pragma omp parallel 
		{
  			printf("thread = %d\n", omp_get_thread_num());
		
			double range=0;
			for (int i = 0; i < n; ++i) {
      				double ke_val_fit=gRandom->Gaus(mu_nom_data,sg_nom_data);
       				gaus_kefit_data[j]->Fill(ke_val_fit);
   			}
		
        		//double kebb=-50; kebb=BB.KEAtLength(ke_ffbeam_MeV, range_reco);

               		//printf("Hello");
		}
		Fit_gaus_kefit_data[j]=VNFit(gaus_kefit_data[j], gaus_kefit_data[j]->GetMean(), 3);
		Fit_gaus_kefit_data[j]->SetLineColor(2);
		Fit_gaus_kefit_data[j]->SetLineStyle(2);
		Fit_gaus_kefit_data[j]->SetName(Form("Fit_gaus_kefit_data_%d",j));
		gaus_kefit_data[j]->Write();
		Fit_gaus_kefit_data[j]->Write();
	}
	fout->Close();

	return 0;
}
