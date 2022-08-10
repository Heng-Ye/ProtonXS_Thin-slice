#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include <utility>      // std::pair, std::make_pair
#include <string>       // std::string
#include <iostream>     // std::cout
#include "TString.h"

#include "TH1.h"
#include "TH2.h"
#include "TRandom2.h"
#include "TFile.h"
#include "TGraphErrors.h"


//#include "TF1.h"
//#include "./headers/BasicParameters.h"
//#include "./headers/BasicFunctions.h"
#include "./headers/BetheBloch.h"

using namespace std;

// number of distribution generated points
#define NGEN 50000000

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
	TF1 *g=new TF1("g",fitg,gg->GetParameter(0)-n_sigma*gg->GetParameter(1),gg->GetParameter(0)+n_sigma*gg->GetParameter(1),5);
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

	double Mu_ref=mu_denom_data;
	double Sg_ref=sg_denom_data;

	//double Mu_exp=mu_nom_mc;
	//double Sg_exp=sg_nom_mc;

	double Mu_exp=mu_nom_data;
	double Sg_exp=sg_nom_data;

	//Output file name -------------------
	//TString str_out="kebb_reweight_mc.root";
	TString str_out="kebb_reweight_data.root";

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

   	int n=NGEN; //number of points in generator

	//KE=KE(range) -------------------------------------------------------------------------------------------------------------------//
	float nke=100;
	float ke_min=0;
	float ke_max=1000;

	const int n_range=201; //# of range data points
	float range=0; //range start
	float d_range=.5; //step size [cm]
	float d_slice=0.01; //slice size [cm]

	TH1F *gaus_ref[n_range]; //served as truth
	TH1F *gaus_exp[n_range]; //measurement
	TF1* Fit_ref[n_range];
	TF1* Fit_exp[n_range];

	vector<float> mu_ref;
	vector<float> er_mu_ref;
	vector<float> sg_ref;
	vector<float> er_sg_ref;

	vector<float> mu_exp;
	vector<float> er_mu_exp;
	vector<float> sg_exp;
	vector<float> er_sg_exp;

	vector<float> range_vec; //vector to store range values
	vector<float> zero;
	for (int j=0; j<n_range; ++j) {
		float tmp_range=range+(float)j*d_range;
		range_vec.push_back(tmp_range);
		zero.push_back(0);
		cout<<"Create Gaussian at range="<<tmp_range<<"cm"<<endl;
		
		gaus_ref[j]=new TH1F(Form("gaus_ref_%d",j),Form("Length:%.1f [cm]",tmp_range), nke, ke_min, ke_max);
		gaus_exp[j]=new TH1F(Form("gaus_exp_%d",j),Form("Length:%.1f [cm]",tmp_range), nke, ke_min, ke_max);
	}

	TH2F *trklen_ke_ref=new TH2F("trklen_ke_ref","",200,0,200,1000,0,1000);
	TH2F *trklen_ke_exp=new TH2F("trklen_ke_exp","",200,0,200,1000,0,1000);
	TH1F *KEend_bb_ref=new TH1F("KEend_bb_ref","",1000,0,1000);
	TH1F *KEend_bb_exp=new TH1F("KEend_bb_exp","",1000,0,1000);

	//TH2F *KEff_KEend_bb_fit=new TH2F("KEff_KEend_bb_fit","",2000,0,200,10000,0,1000);

	TFile *fout = new TFile(str_out.Data(),"RECREATE");
	#pragma omp parallel 
	{ //parallel
  		printf("thread = %d\n", omp_get_thread_num());
		#pragma omp for
		for (int i = 0; i < n; ++i) { //event generator loop
			//keff
      			double keff_ref=gRandom->Gaus(Mu_ref, Sg_ref);
			double keff_exp=gRandom->Gaus(Mu_exp, Sg_exp);

			//travel distance (assuming full energy deposition)
			double range_ref=BB.RangeFromKESpline(keff_ref); //max. length
			double range_exp=BB.RangeFromKESpline(keff_exp); //max. length

			double kend_bb_ref=BB.KEAtLength(keff_ref, range_ref);
			double kend_bb_exp=BB.KEAtLength(keff_exp, range_exp);

			trklen_ke_ref->Fill(range_ref, keff_ref);
			trklen_ke_exp->Fill(range_exp, keff_exp);
			KEend_bb_ref->Fill(kend_bb_ref);
			KEend_bb_exp->Fill(kend_bb_exp);
			//KEff_KEend_bb_fit->Fill(kend_bb_fit,keff_fit);

			for (int j=0; j<n_range; ++j) { //event evolution as traveling inside tpc
				float range_cen=range_vec.at(j);
				float range_plus=range_cen+d_slice;
				float range_munus=range_cen-d_slice;

				if (range_ref>=range_cen) {
					double len_seg_ref=BB.KEAtLength(keff_ref, range_cen);
					gaus_ref[j]->Fill(len_seg_ref);
				}

				if (range_exp>=range_cen) {
					double len_seg_exp=BB.KEAtLength(keff_exp, range_cen);
					gaus_exp[j]->Fill(len_seg_exp);
				}
			} //event evolution as traveling inside tpc
   		} //event generator loop
	} //parallel

	for (int j=0; j<n_range; ++j) { //fit gaussian in each seg
		Fit_ref[j]=VNFit(gaus_ref[j], gaus_ref[j]->GetMean(), 3);
		Fit_ref[j]->SetLineColor(2);
		Fit_ref[j]->SetLineStyle(2);
		Fit_ref[j]->SetName(Form("Fit_ref_%d",j));

		mu_ref.push_back(Fit_ref[j]->GetParameter(0));
		er_mu_ref.push_back(Fit_ref[j]->GetParError(0));
		sg_ref.push_back(Fit_ref[j]->GetParameter(1));
		er_sg_ref.push_back(Fit_ref[j]->GetParError(1));

		Fit_exp[j]=VNFit(gaus_exp[j], gaus_exp[j]->GetMean(), 3);
		Fit_exp[j]->SetLineColor(2);
		Fit_exp[j]->SetLineStyle(2);
		Fit_exp[j]->SetName(Form("Fit_exp_%d",j));

		mu_exp.push_back(Fit_exp[j]->GetParameter(0));
		er_mu_exp.push_back(Fit_exp[j]->GetParError(0));
		sg_exp.push_back(Fit_exp[j]->GetParameter(1));
		er_sg_exp.push_back(Fit_exp[j]->GetParError(1));

		gaus_ref[j]->Write();
		Fit_ref[j]->Write();
		gaus_exp[j]->Write();
		Fit_exp[j]->Write();
		
	} //fit gaussian in each seg

	TGraphErrors *mu_len_ref=new TGraphErrors(range_vec.size(), &range_vec.at(0), &mu_ref.at(0), &zero.at(0), &er_mu_ref.at(0));
	TGraphErrors *sg_len_ref=new TGraphErrors(range_vec.size(), &range_vec.at(0), &sg_ref.at(0), &zero.at(0), &er_sg_ref.at(0));
	mu_len_ref->SetName("mu_len_ref");
	sg_len_ref->SetName("sg_len_ref");

	TGraphErrors *mu_len_exp=new TGraphErrors(range_vec.size(), &range_vec.at(0), &mu_exp.at(0), &zero.at(0), &er_mu_exp.at(0));
	TGraphErrors *sg_len_exp=new TGraphErrors(range_vec.size(), &range_vec.at(0), &sg_exp.at(0), &zero.at(0), &er_sg_exp.at(0));
	mu_len_exp->SetName("mu_len_exp");
	sg_len_exp->SetName("sg_len_exp");



	mu_len_ref->Write();
	sg_len_ref->Write();

	mu_len_exp->Write();
	sg_len_exp->Write();

	trklen_ke_ref->Write();
	trklen_ke_exp->Write();

	KEend_bb_ref->Write();
	KEend_bb_exp->Write();

	//KEff_KEend_bb_fit->Write();


	fout->Close();



	return 0;
}
