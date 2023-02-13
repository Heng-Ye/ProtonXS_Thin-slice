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
	TString fmc="../mc_proton_beamxy_beammom_edept_elossupstream.root";
	TString fdata="/dune/data2/users/hyliao/protonana/v09_39_01/Diff_KEffbeam_KEFit/proton_beamxy_beammom_edept_eloss_runAll.root";
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

	//Plot results --------------------------------------------------------------------------------------------------------------------------------------------
	TCanvas *c_mc_kebeam_kebeam_kefit_stop=new TCanvas(Form("c_mc_kebeam_kebeam_kefit_stop"),"",900, 600);
	TPad *pad1 = new TPad("pad1","",0,0,1,1);
	pad1->Draw();
	pad1->cd();

	c_mc_kebeam_kebeam_kefit_stop->Divide(1,1);
	c_mc_kebeam_kebeam_kefit_stop->cd(1);

	TH2D *f2d=new TH2D("f2d","",300,300,600,64, 0,64);
	f2d->SetTitle("Stopping Protons; KE_{Beam} [MeV]; #mu (KE_{Beam}-KE_{Fit}) [MeV]");
	//f2d->GetXaxis()->SetTitle("KE_{Beam} [MeV]");
	//f2d->GetYaxis()->SetTitle("KE_{Beam}-KE_{Fit} [MeV]");

	int st_plt=6;
	int n_plt=6;
	TGraphErrors *gr_mc_kebeam_kebeam_kefit_stop= new TGraphErrors(n_plt, &KEbeam_cen.at(st_plt), &mu_diff_kebeam_kefit_stop.at(st_plt), &ex.at(st_plt), &err_mu_diff_kebeam_kefit_stop.at(st_plt));
	TGraphErrors *gr_mc_kebeam_kebeam_keff_stop= new TGraphErrors(n_plt, &KEbeam_cen.at(st_plt), &mu_diff_kebeam_keff_stop.at(st_plt), &ex.at(st_plt), &err_mu_diff_kebeam_keff_stop.at(st_plt));
	TGraphErrors *gr_data_kebeam_kebeam_kefit_stop= new TGraphErrors(n_plt, &KEbeam_cen.at(st_plt), &mu_diff_kebeam_kefit_stop_data.at(st_plt), &ex.at(st_plt), &err_mu_diff_kebeam_kefit_stop_data.at(st_plt));

	f2d->Draw();

	gr_mc_kebeam_kebeam_kefit_stop->SetMarkerColor(2);
	gr_mc_kebeam_kebeam_kefit_stop->SetLineColor(2);
	gr_mc_kebeam_kebeam_kefit_stop->Draw("p same");

	gr_mc_kebeam_kebeam_keff_stop->SetMarkerColor(2);
	gr_mc_kebeam_kebeam_keff_stop->SetLineColor(2);
	//gr_mc_kebeam_kebeam_keff_stop->Draw("p same");

	gr_data_kebeam_kebeam_kefit_stop->SetMarkerColor(1);
	gr_data_kebeam_kebeam_kefit_stop->SetLineColor(1);
	gr_data_kebeam_kebeam_kefit_stop->Draw("p same");


	//Fit the energy loss upstream --------------------------------------------------------------------------------------------------------------------------------------//
	float fit_emin=300;
	float fit_emax=600;
	TF1 *ElossFit_mc_kebeam_kebeam_kefit_stop = new TF1("ElossFit_mc_kebeam_kebeam_kefit_stop", ElossFit, fit_emin, fit_emax, 3);
	//TF1 *ElossFit_mc_kebeam_kebeam_keff_stop = new TF1("ElossFit_mc_kebeam_kebeam_keff_stop", ElossFit, fit_emin, fit_emax, 3);
	TF1 *ElossFit_data_kebeam_kebeam_kefit_stop = new TF1("ElossFit_data_kebeam_kebeam_kefit_stop", ElossFit, fit_emin, fit_emax, 3);

	ElossFit_mc_kebeam_kebeam_kefit_stop->SetLineColor(2);
	ElossFit_mc_kebeam_kebeam_kefit_stop->SetLineStyle(2);

	//ElossFit_mc_kebeam_kebeam_keff_stop->SetLineColor(4);
	//ElossFit_mc_kebeam_kebeam_keff_stop->SetLineStyle(2);

	ElossFit_data_kebeam_kebeam_kefit_stop->SetLineColor(1);
	ElossFit_data_kebeam_kebeam_kefit_stop->SetLineStyle(2);

	gr_mc_kebeam_kebeam_kefit_stop->Fit("ElossFit_mc_kebeam_kebeam_kefit_stop");
	ElossFit_mc_kebeam_kebeam_kefit_stop->Draw("same");
	gr_data_kebeam_kebeam_kefit_stop->Fit("ElossFit_data_kebeam_kebeam_kefit_stop");
	ElossFit_data_kebeam_kebeam_kefit_stop->Draw("same");

	pad1->Update(); //this will force the generation of the "stats" box
	TPaveStats *ps1 = (TPaveStats*)gr_data_kebeam_kebeam_kefit_stop->GetListOfFunctions()->FindObject("stats");
	ps1->SetX1NDC(0.12); ps1->SetX2NDC(0.37);
	ps1->SetY1NDC(0.55); ps1->SetY2NDC(0.75);

	//pad1->Update(); //this will force the generation of the "stats" box
	TPaveStats *ps2 = (TPaveStats*)gr_mc_kebeam_kebeam_kefit_stop->GetListOfFunctions()->FindObject("stats");
	ps2->SetX1NDC(0.37); ps2->SetX2NDC(0.62);
	ps2->SetY1NDC(0.55); ps2->SetY2NDC(0.75);
	ps2->SetTextColor(kRed);
	pad1->Modified();

	TLegend *leg = new TLegend(0.15,0.8,0.4,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(gr_data_kebeam_kebeam_kefit_stop, "Data","ep");
	leg->AddEntry(gr_mc_kebeam_kebeam_kefit_stop, "MC","ep");
	leg->SetNColumns(2);
	leg->Draw();

	vector< pair <double,double> > eloss_mc;
	vector< pair <double,double> > eloss_data;
	for (int i=0; i<3; ++i) {
		eloss_mc.push_back(std::make_pair(ElossFit_mc_kebeam_kebeam_kefit_stop->GetParameter(i), ElossFit_mc_kebeam_kebeam_kefit_stop->GetParError(i)));
		eloss_data.push_back(std::make_pair(ElossFit_data_kebeam_kebeam_kefit_stop->GetParameter(i), ElossFit_data_kebeam_kebeam_kefit_stop->GetParError(i)));
	}

	//systematics: fit 1-sigma up and down

	c_mc_kebeam_kebeam_kefit_stop->Print(Form("plot_edept_eloss_upstream/eloss_kebeam.eps"));





	/*
	// ========== For fit up/down shapes ========== //
	TH1D *hint = new TH1D("hint", "Fitted Gaussian with .95 conf.band", 1000., KE_low, KE_high);
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
	hint -> SetStats(false);
	hint -> SetFillColorAlpha(kRed, 0.3);
	hint -> SetMarkerSize(0);
	hint -> SetLineColor(kRed);
	if(!histname.Contains("ff_reco")) hint -> Draw("e3 same");

	TH1D *h_up = new TH1D("h_up", "", 1000., KE_low, KE_high);
	TH1D *h_down = new TH1D("h_down","", 1000., KE_low, KE_high);
	for(int i = 1; i < 1001; i++){
	double this_content = hint -> GetBinContent(i);
	double this_err = hint -> GetBinError(i);
	h_up -> SetBinContent(i, this_content + this_err);
	h_down -> SetBinContent(i, this_content - this_err);
	h_up -> SetBinError(i, (this_content + this_err) / 20.);
	h_down ->  SetBinError(i, (this_content - this_err) / 20.);
	}

	TF1 * h_up_pol_2 = new TF1("h_up_pol_2", "pol2", KE_low + 0., KE_high + 0.);
	h_up -> Fit(h_up_pol_2, "R0", "", KE_low + 0., KE_high + 0.);
	h_up_pol_2 -> SetLineColor(kBlue);
	if(!histname.Contains("ff_reco")) h_up_pol_2 -> Draw("lsame");  

	TF1 * h_down_pol_2 = new TF1("h_down_pol_2", "pol2", KE_low + 0., KE_high + 0.);
	h_down -> Fit(h_down_pol_2, "R0", "", KE_low + 0., KE_high + 0.);
	h_down_pol_2 -> SetLineColor(kCyan);
	if(!histname.Contains("ff_reco")) h_down_pol_2 ->Draw("lsame");

	TLegend * l_err = new TLegend(0.2, 0.30, 0.92, 0.35);
	l_err -> SetNColumns(2);
	l_err -> AddEntry(h_up_pol_2, "Fitted Up", "l");
	l_err -> AddEntry(h_down_pol_2, "Fitted Down", "l");
	l_err -> Draw("same");

	TString p0_up, p1_up, p2_up;
	p0_up = Form("p_{0} = %.3e", h_up_pol_2 -> GetParameter(0));
	p1_up = Form("p_{1} = %.3e", h_up_pol_2 -> GetParameter(1));
	p2_up = Form("p_{2} = %.3e", h_up_pol_2 -> GetParameter(2));
	TLatex latex_p0_up, latex_p1_up, latex_p2_up;
	latex_p0_up.SetNDC();
	latex_p1_up.SetNDC();
	latex_p2_up.SetNDC();
	latex_p0_up.SetTextSize(0.035);
	latex_p1_up.SetTextSize(0.035);
	latex_p2_up.SetTextSize(0.035);
	if(!histname.Contains("ff_reco")){
	latex_p0_up.DrawLatex(0.20, 0.26, p0_up);
	latex_p1_up.DrawLatex(0.20, 0.22, p1_up);
	latex_p2_up.DrawLatex(0.20, 0.18, p2_up);
	}

	TString p0_down, p1_down, p2_down;
	p0_down = Form("p_{0} = %.3e", h_down_pol_2 -> GetParameter(0));
	p1_down = Form("p_{1} = %.3e", h_down_pol_2 -> GetParameter(1));
	p2_down = Form("p_{2} = %.3e", h_down_pol_2 -> GetParameter(2));
	TLatex latex_p0_down, latex_p1_down, latex_p2_down;
	latex_p0_down.SetNDC();
	latex_p1_down.SetNDC();
	latex_p2_down.SetNDC();
	latex_p0_down.SetTextSize(0.035);
	latex_p1_down.SetTextSize(0.035);
	latex_p2_down.SetTextSize(0.035);
	if(!histname.Contains("ff_reco")){
	latex_p0_down.DrawLatex(0.60, 0.26, p0_down);
	latex_p1_down.DrawLatex(0.60, 0.22, p1_down);
	latex_p2_down.DrawLatex(0.60, 0.18, p2_down);
	}
	*/


}
