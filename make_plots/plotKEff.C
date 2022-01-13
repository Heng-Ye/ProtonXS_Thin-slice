#include "../headers/BasicParameters.h"
#include "../headers/BasicAnaFunc.h"

//void plotBeamKE_DataMC(TString fdata, TString fmc, TString outpath) {
void plotBeamKE_DataMC() {
	TString fdata="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_ke.root";
	TString fmc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/mc_ke.root";
	TString outpath="./plots_bmrw/";


	//TString rep="trklen";
	//TString x_axis_label="Proton Track Length [cm]";

	//plot style --------------------------------------------------------------------//
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	gStyle->SetOptStat(0);

	//read data ---------------------------------------------------------------------//
	TFile *f_data = TFile::Open(fdata.Data());
	TH1D *pbeam_data=(TH1D*)f_data->Get("h1d_kebeam");
	TF1 *pbeam_fit_data=(TF1*)f_data->Get("kebeam_fit");
	//TH1D* h1d=(TH1D*)f_data->Get(Form("h1d_%s_stop",rep.Data()));

	TH1D *pbeam_stop_data=(TH1D*)f_data->Get("h1d_kerange_stop");
	TF1 *pbeam_stop_fit_data=(TF1*)f_data->Get("kerange_stop_fit");
	int n_data=pbeam_data->Integral(); 
	int n_stop_data=pbeam_stop_data->Integral(); 
	cout<<"n_data:"<<n_data<<endl;
	cout<<"n_stop_data:"<<n_stop_data<<endl;

	pbeam_data->SetLineColor(1); 	      pbeam_data->SetMarkerColor(1);
	pbeam_fit_data->SetLineColor(1);      pbeam_fit_data->SetMarkerColor(1);
	pbeam_stop_data->SetLineColor(4);     pbeam_stop_data->SetMarkerColor(4);
	pbeam_stop_fit_data->SetLineColor(4); pbeam_stop_fit_data->SetMarkerColor(4);

	//read MC -----------------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fmc.Data());
	TH1D *pbeam_mc=(TH1D*)f_mc->Get("h1d_kebeam");
	TF1 *pbeam_fit_mc=(TF1*)f_mc->Get("kebeam_fit");

	TH1D *pbeam_stop_mc=(TH1D*)f_mc->Get("h1d_kerange_stop");
	TF1 *pbeam_stop_fit_mc=(TF1*)f_mc->Get("kerange_stop_fit");

	TH1D *pff_stop_mc=(TH1D*)f_mc->Get("h1d_keff_stop");
	TF1 *pff_stop_fit_mc=(TF1*)f_mc->Get("keff_stop_fit");

	TH1D *pbeam_bmrw_mc=(TH1D*)f_mc->Get("h1d_kebeam_bmrw");
	TF1 *pbeam_bmrw_fit_mc=(TF1*)f_mc->Get("kebeam_bmrw_fit");

	TH1D *pbeam_stop_bmrw_mc=(TH1D*)f_mc->Get("h1d_kerange_stop_bmrw");
	TF1 *pbeam_stop_bmrw_fit_mc=(TF1*)f_mc->Get("kerange_stop_bmrw_fit");

	TH1D *keff=(TH1D*)f_mc->Get("h1d_keff");
	TH1D *keff2=(TH1D*)f_mc->Get("h1d_keff2");

	TH1D *keff_recoinel=(TH1D*)f_mc->Get("h1d_keff_recoinel");
	TH1D *keff2_recoinel=(TH1D*)f_mc->Get("h1d_keff2_recoinel");

	TH1D *keff_stop=(TH1D*)f_mc->Get("h1d_keff_stop");
	TH1D *keff2_stop=(TH1D*)f_mc->Get("h1d_keff2_stop");

	int n_mc=pbeam_mc->Integral(); 
	int n_stop_mc=pbeam_stop_mc->Integral(); 
	int n_ff_mc=pff_stop_mc->Integral(); 
	cout<<"n_mc:"<<n_mc<<endl;
	cout<<"n_stop_mc:"<<n_stop_mc<<endl;

	pbeam_mc->Scale((double)n_data/(double)n_mc);
	pbeam_stop_mc->Scale((double)n_data/(double)n_mc);
	//pbeam_stop_mc->Scale((double)(n_stop_data)/(double)n_stop_mc);
	pff_stop_mc->Scale((double)n_data/(double)n_mc);
	//pbeam_fit_mc->Scale((double)n_data/(double)n_mc);
	//cout<<"pbeam_mc:"<<pbeam_mc->Integral()<<endl;
	
	double n_mc_bmrw=pbeam_bmrw_mc->Integral(); 
	
	pbeam_bmrw_mc->Scale((double)n_data/(double)n_mc_bmrw);
	pbeam_stop_bmrw_mc->Scale((double)n_stop_data/(double)pbeam_stop_bmrw_mc->Integral());
	//pbeam_stop_bmrw_mc->Scale((double)n_data/(double)n_mc_bmrw);

	pbeam_mc->SetLineColor(2); 	     pbeam_mc->SetMarkerColor(2);
	pbeam_fit_mc->SetLineColor(2);       pbeam_fit_mc->SetMarkerColor(2);

	pbeam_bmrw_mc->SetLineColor(2);	     pbeam_bmrw_mc->SetMarkerColor(2);
	
	pbeam_stop_mc->SetLineColor(3);      pbeam_stop_mc->SetMarkerColor(3);

	pbeam_stop_bmrw_mc->SetLineColor(6); pbeam_stop_bmrw_mc->SetMarkerColor(6);     	


	pbeam_stop_fit_mc->SetLineColor(2);  pbeam_stop_fit_mc->SetMarkerColor(2);
	pff_stop_mc->SetLineColor(6);        pff_stop_mc->SetMarkerColor(6);
	pff_stop_fit_mc->SetLineColor(6);    pff_stop_fit_mc->SetMarkerColor(6);

	//Proton Momentum --------------------------------------------------------------//
	TCanvas *c0=new TCanvas("c0","");
	c0->Divide(1,1);
	c0->cd(1);
	//TH2D* frame2d=new TH2D("frame2d","", 600, 600, 1400, 700, 0, 700); //zend_2d
	TH2D* frame2d=new TH2D("frame2d","", 800, 0, 800, 700, 0, 700); //zend_2d
	frame2d->SetTitle(";Proton Energy [MeV];");
	frame2d->GetXaxis()->CenterTitle();
	frame2d->Draw();
	pbeam_data->Draw("ep same");
	pbeam_stop_data->Draw("ep same");
	//pbeam_stop_fit_data->Draw("same");
	//pbeam_fit_mc->Draw("same");

	pbeam_mc->Draw("hist same");
        pbeam_stop_mc->Draw("hist same");
	pbeam_stop_bmrw_mc->Draw("hist same");

        //TF1* fit_pbeam_mc=VFit(pbeam_mc, 2);
	//fit_pbeam_mc->Draw("same");
        //TF1* fit_pbeam_stop_mc=VFit(pbeam_stop_mc, 3);
	//fit_pbeam_stop_mc->Draw("same");

	//pbeam_stop_fit_mc->Draw("same");
	//read rw histograms -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

	TLegend *leg0 = new TLegend(0.14,0.65,.6,0.85);
	leg0->SetFillStyle(0);
	leg0->AddEntry(pbeam_data, "Data (Protons before entering TPC)", "ep");
	leg0->AddEntry(pbeam_stop_data, "Data (Stopping Protons)", "ep");
	leg0->AddEntry(pbeam_mc, "MC (Protons before entering TPC)", "l");
	leg0->AddEntry(pbeam_stop_mc, "MC (Stopping Protons)", "l");
	leg0->AddEntry(pbeam_stop_bmrw_mc, "MC (Stopping Protons) [after BMRW]", "l");

	leg0->Draw();

        //pDUNE Logo
        TLatex **txt_pdune=new TLatex*[1];
        txt_pdune[0]=new TLatex(0.002, 706, Form("#bf{DUNE:ProtoDUNE-SP}"));
        txt_pdune[0]->SetTextColor(1);
        txt_pdune[0]->SetTextSize(0.06);
        txt_pdune[0]->Draw();
        //
        //Beam Logo
        TLatex **txt_p=new TLatex*[1];
        txt_p[0]=new TLatex(530,706, Form("Protons (1 GeV/c)"));
        txt_p[0]->SetTextColor(1);
        //txt_pdune[0]->SetTextSize(0.07);
        txt_p[0]->Draw();
	c0->Print(Form("%s/kebeam_data_mc.eps",outpath.Data()));




	//print out (mean,sigma) of KEs ----------------------------------------------------//
	//data
	//beam, stopping protons
	double mean_data_kebeam_stop=pbeam_fit_data->GetParameter(0);
	double err_mean_data_kebeam_stop=pbeam_fit_data->GetParError(0);
	double sigma_data_kebeam_stop=pbeam_fit_data->GetParameter(1);
	double err_sigma_data_kebeam_stop=pbeam_fit_data->GetParError(1);

	cout<<"double mean_data_kebeam_stop="<<mean_data_kebeam_stop<<";"<<endl;
	cout<<"double err_mean_data_kebeam_stop="<<err_mean_data_kebeam_stop<<";"<<endl;
	cout<<"double sigma_data_kebeam_stop="<<sigma_data_kebeam_stop<<";"<<endl;
	cout<<"double err_sigma_data_kebeam_stop="<<err_sigma_data_kebeam_stop<<";\n"<<endl;
	cout<<"mean_data_kebeam_stop="<<mean_data_kebeam_stop<<" #pm "<<err_mean_data_kebeam_stop<<endl;
	cout<<"sigma_data_kebeam_stop="<<sigma_data_kebeam_stop<<" #pm "<<err_sigma_data_kebeam_stop<<"\n\n"<<endl;

	//tpc, stopping protons
	double mean_data_kerange_stop=pbeam_stop_fit_data->GetParameter(0);
	double err_mean_data_kerange_stop=pbeam_stop_fit_data->GetParError(0);
	double sigma_data_kerange_stop=pbeam_stop_fit_data->GetParameter(1);
	double err_sigma_data_kerange_stop=pbeam_stop_fit_data->GetParError(1);

	cout<<"double mean_data_kerange_stop="<<mean_data_kerange_stop<<";"<<endl;
	cout<<"double err_mean_data_kerange_stop="<<err_mean_data_kerange_stop<<";"<<endl;
	cout<<"double sigma_data_kerange_stop="<<sigma_data_kerange_stop<<";"<<endl;
	cout<<"double err_sigma_data_kerange_stop="<<err_sigma_data_kerange_stop<<";\n"<<endl;
	cout<<"mean_data_kerange_stop="<<mean_data_kerange_stop<<" #pm "<<err_mean_data_kerange_stop<<endl;
	cout<<"sigma_data_kerange_stop="<<sigma_data_kerange_stop<<" #pm "<<err_sigma_data_kerange_stop<<"\n\n"<<endl;


	double Eloss_data=mean_data_kebeam_stop-mean_data_kerange_stop;
	double err_Eloss_data=sqrt(pow(err_mean_data_kebeam_stop,2)+pow(err_mean_data_kerange_stop,2));
	cout<<"Eloss_data="<<Eloss_data<<";"<<endl;
	cout<<"err_Elsigma_data_range_stoposs_data="<<err_Eloss_data<<";"<<endl;

	//mc
	//beam, stopping protons
	double mean_mc_kebeam_stop=pbeam_fit_mc->GetParameter(0);
	double err_mean_mc_kebeam_stop=pbeam_fit_mc->GetParError(0);
	double sigma_mc_kebeam_stop=pbeam_fit_mc->GetParameter(1);
	double err_sigma_mc_kebeam_stop=pbeam_fit_mc->GetParError(1);

	cout<<"double mean_mc_kebeam_stop="<<mean_mc_kebeam_stop<<";"<<endl;
	cout<<"double err_mean_mc_kebeam_stop="<<err_mean_mc_kebeam_stop<<";"<<endl;
	cout<<"double sigma_mc_kebeam_stop="<<sigma_mc_kebeam_stop<<";"<<endl;
	cout<<"double err_sigma_mc_kebeam_stop="<<err_sigma_mc_kebeam_stop<<";\n"<<endl;
	cout<<"mean_mc_kebeam_stop="<<mean_mc_kebeam_stop<<" #pm "<<err_mean_mc_kebeam_stop<<endl;
	cout<<"sigma_mc_kebeam_stop="<<sigma_mc_kebeam_stop<<" #pm "<<err_sigma_mc_kebeam_stop<<"\n\n"<<endl;


	//tpc, stopping protons
	double mean_mc_kerange_stop=pbeam_stop_fit_mc->GetParameter(0);
	double err_mean_mc_kerange_stop=pbeam_stop_fit_mc->GetParError(0);
	double sigma_mc_kerange_stop=pbeam_stop_fit_mc->GetParameter(1);
	double err_sigma_mc_kerange_stop=pbeam_stop_fit_mc->GetParError(1);

	cout<<"double mean_mc_kerange_stop="<<mean_mc_kerange_stop<<";"<<endl;
	cout<<"double err_mean_mc_kerange_stop="<<err_mean_mc_kerange_stop<<";"<<endl;
	cout<<"double sigma_mc_kerange_stop="<<sigma_mc_kerange_stop<<";"<<endl;
	cout<<"double err_sigma_mc_kerange_stop="<<err_sigma_mc_kerange_stop<<";\n"<<endl;
	cout<<"mean_mc_kerange_stop="<<mean_mc_kerange_stop<<" #pm "<<err_mean_mc_kerange_stop<<endl;
	cout<<"sigma_mc_kerange_stop="<<sigma_mc_kerange_stop<<" #pm "<<err_sigma_mc_kerange_stop<<"\n\n"<<endl;



	double Eloss_mc=mean_mc_kebeam_stop-mean_mc_kerange_stop;
	double err_Eloss_mc=sqrt(pow(err_mean_mc_kebeam_stop,2)+pow(err_mean_mc_kerange_stop,2));
	cout<<"Eloss_mc="<<Eloss_mc<<";"<<endl;
	cout<<"err_Eloss_mc="<<err_Eloss_mc<<";"<<endl;

	//mc, bmrw
	//beam, stopping protons
	double mean_mc_kebeam_stop_bmrw=pbeam_bmrw_fit_mc->GetParameter(0);
	double err_mean_mc_kebeam_stop_bmrw=pbeam_bmrw_fit_mc->GetParError(0);
	double sigma_mc_kebeam_stop_bmrw=pbeam_bmrw_fit_mc->GetParameter(1);
	double err_sigma_mc_kebeam_stop_bmrw=pbeam_bmrw_fit_mc->GetParError(1);

	cout<<"double mean_mc_kebeam_stop_bmrw="<<mean_mc_kebeam_stop_bmrw<<";"<<endl;
	cout<<"double err_mean_mc_kebeam_stop_bmrw="<<err_mean_mc_kebeam_stop_bmrw<<";"<<endl;
	cout<<"double sigma_mc_kebeam_stop_bmrw="<<sigma_mc_kebeam_stop_bmrw<<";"<<endl;
	cout<<"double err_sigma_mc_kebeam_stop_bmrw="<<err_sigma_mc_kebeam_stop_bmrw<<";\n"<<endl;
	cout<<"mean_mc_kebeam_stop_bmrw="<<mean_mc_kebeam_stop_bmrw<<" #pm "<<err_mean_mc_kebeam_stop_bmrw<<endl;
	cout<<"sigma_mc_kebeam_stop_bmrw="<<sigma_mc_kebeam_stop_bmrw<<" #pm "<<err_sigma_mc_kebeam_stop_bmrw<<"\n\n"<<endl;


	//tpc, stopping protons
	double mean_mc_kerange_stop_bmrw=pbeam_stop_bmrw_fit_mc->GetParameter(0);
	double err_mean_mc_kerange_stop_bmrw=pbeam_stop_bmrw_fit_mc->GetParError(0);
	double sigma_mc_kerange_stop_bmrw=pbeam_stop_bmrw_fit_mc->GetParameter(1);
	double err_sigma_mc_kerange_stop_bmrw=pbeam_stop_bmrw_fit_mc->GetParError(1);

	cout<<"double mean_mc_kerange_stop_bmrw="<<mean_mc_kerange_stop_bmrw<<";"<<endl;
	cout<<"double err_mean_mc_kerange_stop_bmrw="<<err_mean_mc_kerange_stop_bmrw<<";"<<endl;
	cout<<"double sigma_mc_kerange_stop_bmrw="<<sigma_mc_kerange_stop_bmrw<<";"<<endl;
	cout<<"double err_sigma_mc_kerange_stop_bmrw="<<err_sigma_mc_kerange_stop_bmrw<<";\n"<<endl;
	cout<<"mean_mc_kerange_stop_bmrw="<<mean_mc_kerange_stop_bmrw<<" #pm "<<err_mean_mc_kerange_stop_bmrw<<endl;
	cout<<"sigma_mc_kerange_stop_bmrw="<<sigma_mc_kerange_stop_bmrw<<" #pm "<<err_sigma_mc_kerange_stop_bmrw<<"\n\n"<<endl;


	double Eloss_bmrw_mc=mean_mc_kebeam_stop_bmrw-mean_mc_kerange_stop_bmrw;
	double err_Eloss_bmrw_mc=sqrt(pow(err_mean_mc_kebeam_stop_bmrw,2)+pow(err_mean_mc_kerange_stop_bmrw,2));
	cout<<"Eloss_bmrw_mc="<<Eloss_bmrw_mc<<";"<<endl;
	cout<<"err_Eloss_bmrw_mc="<<err_Eloss_bmrw_mc<<";"<<endl;

/*
	//Energy loss study
	//KEff=KEbeam-const_Eloss
	//[1]stopping protons
	TCanvas *c0_=new TCanvas("c0_","");
	c0_->Divide(1,1);
	c0_->cd(1);
	TH2D* frame2d_=new TH2D("frame2d_","", 800, 0, 800, 700, 0, 17000); //zend_2d
	frame2d_->SetTitle(";Proton Energy [MeV];");
	frame2d_->GetXaxis()->CenterTitle();
	//frame2d_->Draw();
	keff2->SetLineColor(2);
	keff->Draw("hist");
	keff2->Draw("hist same");
*/

/*
	TCanvas *c0_1=new TCanvas("c0_1","");
	c0_1->Divide(1,1);
	c0_1->cd(1);
	TH2D* frame2d_1=new TH2D("frame2d_1","", 800, 0, 800, 700, 0, 17000); //zend_2d
	frame2d_1->SetTitle(";Proton Energy [MeV];");
	frame2d_1->GetXaxis()->CenterTitle();
	frame2d_1->Draw();
	keff_stop->SetLineColor(1);
	keff2_stop->SetLineColor(2);
	keff2_stop->Scale((double)keff_stop->Integral()/(double)keff2_stop->Integral());
	keff_stop->Draw("hist same");
	keff2_stop->Draw("hist same");
*/


	TCanvas *c0_2=new TCanvas("c0_2","");
	c0_2->Divide(1,1);
	c0_2->cd(1);
	TH2D* frame2d_2=new TH2D("frame2d_2","", 800, 0, 800, 700, 0, 17000); //zend_2d
	frame2d_2->SetTitle(";Proton Energy [MeV];");
	frame2d_2->GetXaxis()->CenterTitle();
	frame2d_2->Draw();
	keff2_recoinel->SetLineColor(2);
	keff_recoinel->SetLineColor(1);
	keff_recoinel->Draw("hist same");
	keff2_recoinel->Draw("hist same");

	
	//[2]reco_inel protons

/*
	TH1D *keff=(TH1D*)f_mc->Get("h1d_keff");
	TH1D *keff2=(TH1D*)f_mc->Get("h1d_keff2");

	TH1D *keff_recoinel=(TH1D*)f_mc->Get("h1d_keff_recoinel");
	TH1D *keff2_recoinel=(TH1D*)f_mc->Get("h1d_keff2_recoinel");

	TH1D *keff_stop=(TH1D*)f_mc->Get("h1d_keff_stop");
	TH1D *keff2_stop=(TH1D*)f_mc->Get("h1d_keff2_stop");
*/




/*

	//Proton Momentum [only MC] --------------------------------------------------------------//
	TCanvas *c01=new TCanvas("c01","");
	c01->Divide(1,1);
	c01->cd(1);
	TH2D* frame2d01=new TH2D("frame2d01","", 600, 600, 1400, 700, 0, 700); //zend_2d
	frame2d01->SetTitle(";Proton Momentum [MeV/c];");
	frame2d01->GetXaxis()->CenterTitle();
	frame2d01->Draw();
	pbeam_mc->Draw("hist same");
	fit_pbeam_mc->Draw("hist same");

        pbeam_stop_mc->Draw("hist same");
	fit_pbeam_stop_mc->Draw("same");

	pff_stop_mc->Draw("hist same");
        TF1* fit_pff_stop_mc=VFit(pff_stop_mc, 6);
	fit_pff_stop_mc->Draw("same");


	TLegend *leg01 = new TLegend(0.14,0.65,.6,0.85);
	leg01->SetFillStyle(0);
	leg01->AddEntry(pbeam_mc, "MC (Protons before entering TPC)", "l");
	leg01->AddEntry(pff_stop_mc, "MC (Protons at TPC Front Face)", "l");
	leg01->AddEntry(pbeam_stop_mc, "MC (Stopping Protons)", "l");
	leg01->Draw();

        //pDUNE Logo
        //TLatex **txt_pdune=new TLatex*[1];
        //txt_pdune[0]=new TLatex(600.002, 706, Form("#bf{DUNE:ProtoDUNE-SP}"));
        //txt_pdune[0]->SetTextColor(1);
        txt_pdune[0]->Draw();
        //
        //Beam Logo
        //TLatex **txt_p=new TLatex*[1];
        //txt_p[0]=new TLatex(1145,706, Form("Protons (1 GeV/c)"));
        //txt_p[0]->SetTextColor(1);
        //txt_pdune[0]->SetTextSize(0.07);
        txt_p[0]->Draw();


	c01->Print(Form("%s/pbeam_only_mc.eps",outpath.Data()));




	//Beam (mu, sigma) -------------------------------------------------------------------------------------------------------------------//
	vector<float> mu_data;
	vector<float> err_mu_data;
	vector<float> sigma_data;
	vector<float> err_sigma_data;
	mu_data.push_back(pbeam_fit_data->GetParameter(0));		err_mu_data.push_back(pbeam_fit_data->GetParError(0));
	sigma_data.push_back(pbeam_fit_data->GetParameter(1));		err_sigma_data.push_back(pbeam_fit_data->GetParError(1));
	mu_data.push_back(pbeam_stop_fit_data->GetParameter(0));	err_mu_data.push_back(pbeam_stop_fit_data->GetParError(0));
	sigma_data.push_back(pbeam_stop_fit_data->GetParameter(1));	err_sigma_data.push_back(pbeam_stop_fit_data->GetParError(1));

	vector<float> mu_mc;
	vector<float> err_mu_mc;
	vector<float> sigma_mc;
	vector<float> err_sigma_mc;
	mu_mc.push_back(pbeam_fit_mc->GetParameter(0));		err_mu_mc.push_back(pbeam_fit_mc->GetParError(0));
	sigma_mc.push_back(pbeam_fit_mc->GetParameter(1));	err_sigma_mc.push_back(pbeam_fit_mc->GetParError(1));

	//mu_mc.push_back(pff_stop_fit_mc->GetParameter(0));	err_mu_mc.push_back(pff_stop_fit_mc->GetParError(0));
	//sigma_mc.push_back(pff_stop_fit_mc->GetParameter(1));	err_sigma_mc.push_back(pff_stop_fit_mc->GetParError(1));

	mu_mc.push_back(pbeam_stop_fit_mc->GetParameter(0));	err_mu_mc.push_back(pbeam_stop_fit_mc->GetParError(0));
	sigma_mc.push_back(pbeam_stop_fit_mc->GetParameter(1));	err_sigma_mc.push_back(pbeam_stop_fit_mc->GetParError(1));

	TGraphErrors *ms_data=new TGraphErrors(mu_data.size(), &mu_data.at(0), &sigma_data.at(0), &err_mu_data.at(0), &err_sigma_data.at(0));
	TGraphErrors *ms_mc=new TGraphErrors(mu_mc.size(), &mu_mc.at(0), &sigma_mc.at(0), &err_mu_mc.at(0), &err_sigma_mc.at(0));
	
	ms_mc->SetMarkerColor(2);
	ms_mc->SetLineColor(2);

	TCanvas *c1=new TCanvas("c1","");
	c1->Divide(1,1);
	c1->cd(1);
	TH2D* frame2d2=new TH2D("frame2d2","", 200, 930, 1040, 42, 48, 90);
	frame2d2->SetTitle(";#mu [MeV/c];#sigma [MeV/c]");
	frame2d2->GetXaxis()->CenterTitle();
	frame2d2->Draw();
	ms_data->Draw("p same");
	ms_mc->Draw("p same");

	float emin=930.;
	float emax=1040.;
	TF1 *fit_data = new TF1("fit_data", "[0]+[1]*x", emin, emax);
	fit_data->SetLineColor(1);
	fit_data->SetLineStyle(2);
	ms_data->Fit(fit_data,"remn");
	fit_data->Draw("same");
	
        TF1 *fit_mc = new TF1("fit_mc", "[0]+[1]*x", emin, emax);
        fit_mc->SetLineColor(2);
        fit_mc->SetLineStyle(2);
        ms_mc->Fit(fit_mc,"remn");
        fit_mc->Draw("same");

	TLegend *leg1 = new TLegend(0.14,0.65,.6,0.85);
	leg1->SetFillStyle(0);
	leg1->AddEntry(ms_data, Form("Data: #sigma=%.2f+%.2f*#mu",fit_data->GetParameter(0),fit_data->GetParameter(1)), "ep");
	leg1->AddEntry(ms_mc, Form("MC: #sigma=%.2f+%.2f*#mu",fit_mc->GetParameter(0),fit_mc->GetParameter(1)), "ep");
	leg1->Draw();

        //pDUNE Logo
        TLatex **txt_pdune1=new TLatex*[1];
        txt_pdune1[0]=new TLatex(930.002, 90.6, Form("#bf{DUNE:ProtoDUNE-SP}"));
        txt_pdune1[0]->SetTextColor(1);
        //txt_pdune1[0]->Draw();
        //
        //Beam Logo
        TLatex **txt_p1=new TLatex*[1];
        txt_p1[0]=new TLatex(1002,90.6, Form("Protons (1 GeV/c)"));
        txt_p1[0]->SetTextColor(1);
        txt_p1[0]->SetTextSize(0.05);
        //txt_p1[0]->Draw();


	c1->Print(Form("%s/pbeam_mus_igma_data_mc.eps",outpath.Data()));
*/




}
