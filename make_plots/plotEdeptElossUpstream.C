#include "../headers/BasicParameters.h"
#include "../headers/BasicFunctions.h"
#include "TGraphErrors.h"

// Quadratic function
Double_t ElossFit(Double_t *x, Double_t *par) {
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

void plotEdeptElossUpstream() {
	//Input file names -------------------------------------------------//
	//TString fmc="../mc_proton_beamxy_beammom_edept_elossupstream_old.root";
	TString fmc="../mc_proton_beamxy_beammom_edept_elossupstream_old.root";
	TString fdata="/dune/data2/users/hyliao/protonana/v09_39_01/Diff_KEffbeam_KEFit_old/proton_beamxy_beammom_edept_eloss_runAll.root";
	TFile *f_mc = TFile::Open(fmc.Data());
	TFile *f_data = TFile::Open(fdata.Data());

        //plot style --------------------------------------------------------------------//
        gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
        gROOT->SetStyle("protoDUNEStyle");
        gROOT->ForceStyle();
        gStyle->SetTitleX(0.5);
        gStyle->SetTitleAlign(23);
        gStyle->SetOptStat(0);

	//Energy-dependent E-loss upstream ------------------------------------------------------------------//
	const int n_kebeam_slice=14; //14 energy slicing in total
	int d_kebeam=50; //50 MeV Slice
	int kebeam_min=0;
	vector<int> KEbeam_slice;
	vector<float> KEbeam_cen;
	vector<float> ex;

	TH1D **diff_kebeam_kefit_stop=new TH1D*[n_kebeam_slice];
	TH1D **diff_kebeam_kefit_el=new TH1D*[n_kebeam_slice];
	TH1D **diff_kebeam_keff_stop=new TH1D*[n_kebeam_slice];
	TH1D **diff_kebeam_keff_el=new TH1D*[n_kebeam_slice];

	TH1D **diff_kebeam_kefit_stop_data=new TH1D*[n_kebeam_slice];

	TCanvas *c_mc_diff_kebeam_kefit_stop=new TCanvas(Form("c_mc_diff_kebeam_kefit_stop"),"",1200, 1200);
	c_mc_diff_kebeam_kefit_stop->Divide(4,4);

	TCanvas *c_mc_diff_kebeam_kefit_stop_data=new TCanvas(Form("c_mc_diff_kebeam_kefit_stop_data"),"",1200, 1200);
	c_mc_diff_kebeam_kefit_stop_data->Divide(4,4);

	TF1 **fit_diff_kebeam_kefit_stop=new TF1*[n_kebeam_slice];
	TF1 **fit_diff_kebeam_kefit_el=new TF1*[n_kebeam_slice];
	TF1 **fit_diff_kebeam_keff_stop=new TF1*[n_kebeam_slice];
	TF1 **fit_diff_kebeam_keff_el=new TF1*[n_kebeam_slice];

	TF1 **fit_diff_kebeam_kefit_stop_data=new TF1*[n_kebeam_slice];

	//TF1* fit_kend_bb_inel_data; fit_kend_bb_inel_data=VFit(kend_bb_inel_data, 1);

	vector<float> mu_diff_kebeam_kefit_stop;
	vector<float> err_mu_diff_kebeam_kefit_stop;

	vector<float> mu_diff_kebeam_keff_stop;
	vector<float> err_mu_diff_kebeam_keff_stop;

	vector<float> mu_diff_kebeam_kefit_stop_data;
	vector<float> err_mu_diff_kebeam_kefit_stop_data;

	for (int ii=0; ii<n_kebeam_slice; ++ii) {
		kebeam_min+=d_kebeam;
		KEbeam_slice.push_back(kebeam_min);

		KEbeam_cen.push_back((float)kebeam_min-0.5*(float)d_kebeam);
		ex.push_back(0.5*(float)d_kebeam);

		diff_kebeam_kefit_stop[ii]=(TH1D *)f_mc->Get(Form("diff_kebeam_kefit_stop_%d",ii));
		diff_kebeam_kefit_el[ii]=(TH1D *)f_mc->Get(Form("diff_kebeam_kefit_el_%d",ii));
		diff_kebeam_keff_stop[ii]=(TH1D *)f_mc->Get(Form("diff_kebeam_keff_stop_%d",ii));
		diff_kebeam_keff_el[ii]=(TH1D *)f_mc->Get(Form("diff_kebeam_keff_el_%d",ii));

		diff_kebeam_kefit_stop_data[ii]=(TH1D *)f_data->Get(Form("diff_kebeam_kefit_stop_%d",ii));


		c_mc_diff_kebeam_kefit_stop->cd(ii+1);

		diff_kebeam_kefit_stop[ii]->Draw();
		fit_diff_kebeam_kefit_stop[ii]=VFit(diff_kebeam_kefit_stop[ii], 1);
		mu_diff_kebeam_kefit_stop.push_back(fit_diff_kebeam_kefit_stop[ii]->GetParameter(0));
		err_mu_diff_kebeam_kefit_stop.push_back(fit_diff_kebeam_kefit_stop[ii]->GetParError(0));

		diff_kebeam_kefit_el[ii]->SetLineColor(2);	
		diff_kebeam_kefit_el[ii]->Draw("same");
		fit_diff_kebeam_kefit_el[ii]=VFit(diff_kebeam_kefit_el[ii], 2);
		fit_diff_kebeam_kefit_el[ii]->Draw("same");
		
		diff_kebeam_keff_stop[ii]->SetLineColor(4);
		diff_kebeam_keff_stop[ii]->Draw("same");
		fit_diff_kebeam_keff_stop[ii]=VFit(diff_kebeam_keff_stop[ii], 4);
		fit_diff_kebeam_keff_stop[ii]->Draw("same");
		mu_diff_kebeam_keff_stop.push_back(fit_diff_kebeam_keff_stop[ii]->GetParameter(0));
		err_mu_diff_kebeam_keff_stop.push_back(fit_diff_kebeam_keff_stop[ii]->GetParError(0));

		
		diff_kebeam_keff_el[ii]->SetLineColor(3);
		diff_kebeam_keff_el[ii]->Draw("same");
		fit_diff_kebeam_keff_el[ii]=VFit(diff_kebeam_keff_el[ii], 3);
		fit_diff_kebeam_keff_el[ii]->Draw("same");

		c_mc_diff_kebeam_kefit_stop_data->cd(ii+1);
		diff_kebeam_kefit_stop_data[ii]->Draw();
		fit_diff_kebeam_kefit_stop_data[ii]=VFit(diff_kebeam_kefit_stop_data[ii], 2);
		fit_diff_kebeam_kefit_stop_data[ii]->Draw("same");
		mu_diff_kebeam_kefit_stop_data.push_back(fit_diff_kebeam_kefit_stop_data[ii]->GetParameter(0));
		err_mu_diff_kebeam_kefit_stop_data.push_back(fit_diff_kebeam_kefit_stop_data[ii]->GetParError(0));
			

	}
	//gg->SetParameter(0,pre_mean);
	//gg->SetParameter(1,pre_rms);
	//gg->SetParameter(2,pre_max);
	//c_hit_int->Print(Form("%sHit_matrix_int.eps",outpath.Data()));


	TCanvas *c_mc_kebeam_kebeam_kefit_stop=new TCanvas(Form("c_mc_kebeam_kebeam_kefit_stop"),"",900, 600);
	c_mc_kebeam_kebeam_kefit_stop->Divide(1,1);
	c_mc_kebeam_kebeam_kefit_stop->cd(1);
	TH2D *f2d=new TH2D("f2d","",350,300,650,65,0,65);
	f2d->GetXaxis()->SetTitle("KE_{Beam} [MeV]");
	f2d->GetYaxis()->SetTitle("KE_{Beam}-KE_{Fit} [MeV]");
        int st_plt=6;
	int n_plt=6;
	TGraphErrors *gr_mc_kebeam_kebeam_kefit_stop= new TGraphErrors(n_plt, &KEbeam_cen.at(st_plt), &mu_diff_kebeam_kefit_stop.at(st_plt), &ex.at(st_plt), &err_mu_diff_kebeam_kefit_stop.at(st_plt));
	TGraphErrors *gr_mc_kebeam_kebeam_keff_stop= new TGraphErrors(n_plt, &KEbeam_cen.at(st_plt), &mu_diff_kebeam_keff_stop.at(st_plt), &ex.at(st_plt), &err_mu_diff_kebeam_keff_stop.at(st_plt));
	TGraphErrors *gr_mc_kebeam_kebeam_kefit_stop_data= new TGraphErrors(n_plt, &KEbeam_cen.at(st_plt), &mu_diff_kebeam_kefit_stop_data.at(st_plt), &ex.at(st_plt), &err_mu_diff_kebeam_kefit_stop_data.at(st_plt));
	
	f2d->Draw();
	gr_mc_kebeam_kebeam_kefit_stop->Draw("p same");

	gr_mc_kebeam_kebeam_keff_stop->SetMarkerColor(2);
	gr_mc_kebeam_kebeam_keff_stop->SetLineColor(4);
	gr_mc_kebeam_kebeam_keff_stop->Draw("p same");

	gr_mc_kebeam_kebeam_kefit_stop_data->SetMarkerColor(4);
	gr_mc_kebeam_kebeam_kefit_stop_data->SetLineColor(4);
	gr_mc_kebeam_kebeam_kefit_stop_data->Draw("p same");

	
	// create a TF1 with the range from 0 to 3 and 6 parameters
	float fit_emin=300;
	float fit_emax=650;
   	TF1 *ElossFit_mc_kebeam_kebeam_kefit_stop = new TF1("ElossFit_mc_kebeam_kebeam_kefit_stop", ElossFit, fit_emin, fit_emax, 3);
   	TF1 *ElossFit_mc_kebeam_kebeam_keff_stop = new TF1("ElossFit_mc_kebeam_kebeam_keff_stop", ElossFit, fit_emin, fit_emax, 3);
   	TF1 *ElossFit_mc_kebeam_kebeam_kefit_stop_data = new TF1("ElossFit_mc_kebeam_kebeam_kefit_stop_data", ElossFit, fit_emin, fit_emax, 3);


	ElossFit_mc_kebeam_kebeam_kefit_stop->SetLineColor(2);
	ElossFit_mc_kebeam_kebeam_kefit_stop->SetLineStyle(2);

	ElossFit_mc_kebeam_kebeam_keff_stop->SetLineColor(4);
	ElossFit_mc_kebeam_kebeam_keff_stop->SetLineStyle(2);

	gr_mc_kebeam_kebeam_kefit_stop->Fit("ElossFit_mc_kebeam_kebeam_kefit_stop","rem+");
	gr_mc_kebeam_kebeam_kefit_stop->Draw("same");


	ElossFit_mc_kebeam_kebeam_kefit_stop_data->SetLineColor(4);
	ElossFit_mc_kebeam_kebeam_kefit_stop_data->SetLineStyle(2);
	gr_mc_kebeam_kebeam_kefit_stop_data->Fit("ElossFit_mc_kebeam_kebeam_kefit_stop_data","rem+");
	gr_mc_kebeam_kebeam_kefit_stop_data->Draw("same");




}
