#include "../headers/BasicParameters.h"
#include "../headers/BasicFunctions.h"

void plotDataMCKE_BMRW(TString fdata, TString fmc, TString fmc_bmrw, TString outpath) {

	TString rep="trklen";
	TString x_axis_label="Proton Track Length [cm]";

	//plot style --------------------------------------------------------------------//
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	gStyle->SetOptStat(0);

	//read data ----------------------------------------------------//
	TFile *f_data = TFile::Open(fdata.Data());
	TH1D *kebeam_data=(TH1D*)f_data->Get("h1d_kebeam");
	TH1D *pbeam_data=(TH1D*)f_data->Get("h1d_pbeam");

	TH1D *kebeam_stop_data=(TH1D*)f_data->Get("h1d_kebeam_stop");
	TH1D *pbeam_stop_data=(TH1D*)f_data->Get("h1d_pbeam_stop");

	TH1D *kerange_stop_data=(TH1D*)f_data->Get("h1d_kerange_stop");
	TH1D *prange_stop_data=(TH1D*)f_data->Get("h1d_prange_stop");

	TH1D *kecalo_stop_data=(TH1D*)f_data->Get("h1d_kecalo_stop");
	TH1D *pcalo_stop_data=(TH1D*)f_data->Get("h1d_pcalo_stop");

	float n_pbeam_stop_data=pbeam_stop_data->Integral();
	float n_prange_stop_data=prange_stop_data->Integral();
	float n_pcalo_stop_data=pcalo_stop_data->Integral();

	//read MC ------------------------------------------------------//
	TFile *f_mc = TFile::Open(fmc.Data());
	TH1D *kebeam_mc=(TH1D*)f_mc->Get("h1d_kebeam");
	TH1D *pbeam_mc=(TH1D*)f_mc->Get("h1d_pbeam");

	TH1D *kebeam_stop_mc=(TH1D*)f_mc->Get("h1d_kebeam_stop");
	TH1D *pbeam_stop_mc=(TH1D*)f_mc->Get("h1d_pbeam_stop");

	TH1D *ke0_mc=(TH1D*)f_mc->Get("h1d_ke0");
	TH1D *p0_mc=(TH1D*)f_mc->Get("h1d_p0");

	TH1D *ke0_stop_mc=(TH1D*)f_mc->Get("h1d_ke0_stop");
	TH1D *p0_stop_mc=(TH1D*)f_mc->Get("h1d_p0_stop");

	TH1D *keff_mc=(TH1D*)f_mc->Get("h1d_keff");
	TH1D *pff_mc=(TH1D*)f_mc->Get("h1d_pff");

	TH1D *keff_stop_mc=(TH1D*)f_mc->Get("h1d_keff_stop");
	TH1D *pff_stop_mc=(TH1D*)f_mc->Get("h1d_pff_stop");

	TH1D *kerange_stop_mc=(TH1D*)f_mc->Get("h1d_kerange_stop");
	TH1D *prange_stop_mc=(TH1D*)f_mc->Get("h1d_prange_stop");


	TH1D *kecalo_stop_mc=(TH1D*)f_mc->Get("h1d_kecalo_stop");
	TH1D *pcalo_stop_mc=(TH1D*)f_mc->Get("h1d_pcalo_stop");


	float n_pbeam_stop_mc=pbeam_stop_mc->Integral();
	float n_p0_stop_mc=p0_stop_mc->Integral();
	float n_prange_stop_mc=prange_stop_mc->Integral();
	float n_pcalo_stop_mc=pcalo_stop_mc->Integral();
	float n_pff_stop_mc=pff_stop_mc->Integral();



	//read MC[bmrw] ------------------------------------------------------//
	TFile *f_mc_bmrw = TFile::Open(fmc_bmrw.Data());
	TH1D *kebeam_mc_bmrw=(TH1D*)f_mc_bmrw->Get("h1d_kebeam");
	TH1D *pbeam_mc_bmrw=(TH1D*)f_mc_bmrw->Get("h1d_pbeam");

	TH1D *kebeam_stop_mc_bmrw=(TH1D*)f_mc_bmrw->Get("h1d_kebeam_stop");
	TH1D *pbeam_stop_mc_bmrw=(TH1D*)f_mc_bmrw->Get("h1d_pbeam_stop");

	TH1D *ke0_mc_bmrw=(TH1D*)f_mc_bmrw->Get("h1d_ke0");
	TH1D *p0_mc_bmrw=(TH1D*)f_mc_bmrw->Get("h1d_p0");

	TH1D *ke0_stop_mc_bmrw=(TH1D*)f_mc_bmrw->Get("h1d_ke0_stop");
	TH1D *p0_stop_mc_bmrw=(TH1D*)f_mc_bmrw->Get("h1d_p0_stop");

	TH1D *keff_mc_bmrw=(TH1D*)f_mc_bmrw->Get("h1d_keff");
	TH1D *pff_mc_bmrw=(TH1D*)f_mc_bmrw->Get("h1d_pff");

	TH1D *keff_stop_mc_bmrw=(TH1D*)f_mc_bmrw->Get("h1d_keff_stop");
	TH1D *pff_stop_mc_bmrw=(TH1D*)f_mc_bmrw->Get("h1d_pff_stop");

	TH1D *kerange_stop_mc_bmrw=(TH1D*)f_mc_bmrw->Get("h1d_kerange_stop");
	TH1D *prange_stop_mc_bmrw=(TH1D*)f_mc_bmrw->Get("h1d_prange_stop");


	TH1D *kecalo_stop_mc_bmrw=(TH1D*)f_mc_bmrw->Get("h1d_kecalo_stop");
	TH1D *pcalo_stop_mc_bmrw=(TH1D*)f_mc_bmrw->Get("h1d_pcalo_stop");


	float n_pbeam_stop_mc_bmrw=pbeam_stop_mc_bmrw->Integral();
	float n_p0_stop_mc_bmrw=p0_stop_mc_bmrw->Integral();
	float n_prange_stop_mc_bmrw=prange_stop_mc_bmrw->Integral();
	float n_pcalo_stop_mc_bmrw=pcalo_stop_mc_bmrw->Integral();
	float n_pff_stop_mc_bmrw=pff_stop_mc_bmrw->Integral();


	//Proton Momentum --------------------------------------------------------------//
	pbeam_stop_mc->Scale(n_pbeam_stop_data/n_pbeam_stop_mc);
	pbeam_stop_mc_bmrw->Scale(n_pbeam_stop_data/n_pbeam_stop_mc_bmrw);
	p0_stop_mc->Scale(n_pbeam_stop_data/n_p0_stop_mc);
	p0_stop_mc_bmrw->Scale(n_pbeam_stop_data/n_p0_stop_mc_bmrw);
	//Fit Gaussians on momenta ...
	//[1]
	TF1* fit_pbeam_stop_data; fit_pbeam_stop_data=VFit(pbeam_stop_data, 1);
	fit_pbeam_stop_data->SetName("fit_pbeam_stop_data");
	fit_pbeam_stop_data->SetLineStyle(2);
	double m_stop_data=fit_pbeam_stop_data->GetParameter(0); //Data prod4 reco2
	double err_m_stop_data=fit_pbeam_stop_data->GetParError(0);
	double s_stop_data=fit_pbeam_stop_data->GetParameter(1); //Data prod4 reco2
	double err_s_stop_data=fit_pbeam_stop_data->GetParError(1);

	//[2]spec
	TF1* fit_pbeam_stop_mc; fit_pbeam_stop_mc=VFit(pbeam_stop_mc, 2);
	fit_pbeam_stop_mc->SetName("fit_pbeam_stop_mc");
	fit_pbeam_stop_mc->SetLineStyle(2);
	double m_stop_mc=fit_pbeam_stop_mc->GetParameter(0); //
	double err_m_stop_mc=fit_pbeam_stop_mc->GetParError(0);
	double s_stop_mc=fit_pbeam_stop_mc->GetParameter(1); //
	double err_s_stop_mc=fit_pbeam_stop_mc->GetParError(1);

	//[3]
	TF1* fit_p0_stop_mc; fit_p0_stop_mc=VFit(p0_stop_mc, 3);
	fit_p0_stop_mc->SetName("fit_p0_stop_mc");
	fit_p0_stop_mc->SetLineStyle(2);
	double m0_stop_mc=fit_p0_stop_mc->GetParameter(0); //
	double err_m0_stop_mc=fit_p0_stop_mc->GetParError(0);
	double s0_stop_mc=fit_p0_stop_mc->GetParameter(1); //
	double err_s0_stop_mc=fit_p0_stop_mc->GetParError(1);

	//[2_bmrw]
	TF1* fit_pbeam_stop_mc_bmrw; fit_pbeam_stop_mc_bmrw=VFit(pbeam_stop_mc_bmrw, 4);
	fit_pbeam_stop_mc_bmrw->SetName("fit_pbeam_stop_mc_bmrw");
	fit_pbeam_stop_mc_bmrw->SetLineStyle(2);
	double m_stop_mc_bmrw=fit_pbeam_stop_mc_bmrw->GetParameter(0); //
	double err_m_stop_mc_bmrw=fit_pbeam_stop_mc_bmrw->GetParError(0);
	double s_stop_mc_bmrw=fit_pbeam_stop_mc_bmrw->GetParameter(1); //
	double err_s_stop_mc_bmrw=fit_pbeam_stop_mc_bmrw->GetParError(1);

	//[3_bmrw]
	TF1* fit_p0_stop_mc_bmrw; fit_p0_stop_mc_bmrw=VFit(p0_stop_mc_bmrw, 7);
	fit_p0_stop_mc_bmrw->SetName("fit_p0_stop_mc_bmrw");
	fit_p0_stop_mc_bmrw->SetLineStyle(2);
	double m0_stop_mc_bmrw=fit_p0_stop_mc_bmrw->GetParameter(0); //
	double err_m0_stop_mc_bmrw=fit_p0_stop_mc_bmrw->GetParError(0);
	double s0_stop_mc_bmrw=fit_p0_stop_mc_bmrw->GetParameter(1); //
	double err_s0_stop_mc_bmrw=fit_p0_stop_mc_bmrw->GetParError(1);


	TCanvas *c0=new TCanvas("c0","",1200,900);
	c0->Divide(1,1);
	c0->cd(1);
	TH2D* f2d_p=new TH2D("f2d_p","", 700, 600, 1300, 600, 0, pbeam_stop_data->GetBinContent(pbeam_stop_data->GetMaximumBin())+150); //
	f2d_p->SetTitle(";Initial Proton Momentum [MeV/c];");
	f2d_p->GetXaxis()->CenterTitle();
	f2d_p->Draw();
	pbeam_stop_data->SetLineColor(1);
	pbeam_stop_mc->SetLineColor(2);
	pbeam_stop_mc_bmrw->SetLineColor(4);
	p0_stop_mc->SetLineColor(3);
	p0_stop_mc_bmrw->SetLineColor(7);

	fit_pbeam_stop_data->Draw("same");
	fit_pbeam_stop_mc->Draw("same");
	fit_p0_stop_mc->Draw("same");
	fit_pbeam_stop_mc_bmrw->Draw("same");
	fit_p0_stop_mc_bmrw->Draw("same");

	pbeam_stop_mc->Draw("hist same");
	pbeam_stop_mc_bmrw->Draw("hist same");
	p0_stop_mc->Draw("hist same");
	p0_stop_mc_bmrw->Draw("hist same");
	pbeam_stop_data->Draw("ep same");

	TLegend *leg0 = new TLegend(0.14,0.65,.9,0.9);
	leg0->SetFillStyle(0);
        //leg0->AddEntry((TObject*)0, "", "");     
	leg0->AddEntry(pbeam_stop_data, Form("Data: #mu=%.1f#pm%.1f MeV/c, #sigma=%.1f#pm%.1f MeV/c",m_stop_data, err_m_stop_data, s_stop_data, err_s_stop_data), "ep");
	leg0->AddEntry(p0_stop_mc, Form("MC(truth): #mu=%.1f#pm%.1f MeV/c, #sigma=%.1f#pm%.1f MeV/c",m0_stop_mc,err_m0_stop_mc,s0_stop_mc,err_s0_stop_mc), "l");
	leg0->AddEntry(pbeam_stop_mc, Form("MC(spec): #mu=%.1f#pm%.1f MeV/c, #sigma=%.1f#pm%.1f MeV/c",m_stop_mc,err_m_stop_mc,s_stop_mc,err_s_stop_mc), "l");
	leg0->AddEntry(p0_stop_mc_bmrw, Form("MC(truth)+bmrw: #mu=%.1f#pm%.1f MeV/c, #sigma=%.1f#pm%.1f MeV/c",m0_stop_mc_bmrw,err_m0_stop_mc_bmrw,s0_stop_mc_bmrw,err_s0_stop_mc_bmrw), "l");
	leg0->AddEntry(pbeam_stop_mc_bmrw, Form("MC(spec)+bmrw: #mu=%.1f#pm%.1f MeV/c, #sigma=%.1f#pm%.1f MeV/c",m_stop_mc_bmrw,err_m_stop_mc_bmrw,s_stop_mc_bmrw,err_s_stop_mc_bmrw), "l");
        //leg0->SetNColumns(2);
	leg0->Draw();
	c0->Print(Form("%s/mom_ini_bmrw.eps",outpath.Data()));


	//Edept in TPC -------------------------------------------------------------------------------------------------------------------------//
	prange_stop_mc->Scale(n_prange_stop_data/n_prange_stop_mc);
	pcalo_stop_mc->Scale(n_pcalo_stop_data/n_pcalo_stop_mc);
	pff_stop_mc->Scale(n_prange_stop_data/n_pff_stop_mc);

	prange_stop_mc_bmrw->Scale(n_prange_stop_data/n_prange_stop_mc_bmrw);
	pcalo_stop_mc_bmrw->Scale(n_pcalo_stop_data/n_pcalo_stop_mc_bmrw);
	pff_stop_mc_bmrw->Scale(n_prange_stop_data/n_pff_stop_mc_bmrw);

	//Fit Gaussians on momenta ...
	//[1]
	TF1* fit_prange_stop_data; fit_prange_stop_data=VFit(prange_stop_data, 1);
	fit_prange_stop_data->SetName("fit_prange_stop_data");
	fit_prange_stop_data->SetLineStyle(2);
	double m_range_stop_data=fit_prange_stop_data->GetParameter(0); //Data prod4 reco2
	double err_m_range_stop_data=fit_prange_stop_data->GetParError(0);
	double s_range_stop_data=fit_prange_stop_data->GetParameter(1); //Data prod4 reco2
	double err_s_range_stop_data=fit_prange_stop_data->GetParError(1);

	//[2]
	TF1* fit_prange_stop_mc; fit_prange_stop_mc=VFit(prange_stop_mc, 2);
	fit_prange_stop_mc->SetName("fit_prange_stop_mc");
	fit_prange_stop_mc->SetLineStyle(2);
	double m_range_stop_mc=fit_prange_stop_mc->GetParameter(0); //
	double err_m_range_stop_mc=fit_prange_stop_mc->GetParError(0);
	double s_range_stop_mc=fit_prange_stop_mc->GetParameter(1); //
	double err_s_range_stop_mc=fit_prange_stop_mc->GetParError(1);

	//[2_bmrw]
	TF1* fit_prange_stop_mc_bmrw; fit_prange_stop_mc_bmrw=VFit(prange_stop_mc_bmrw, 4);
	fit_prange_stop_mc_bmrw->SetName("fit_prange_stop_mc_bmrw");
	fit_prange_stop_mc_bmrw->SetLineStyle(2);
	double m_range_stop_mc_bmrw=fit_prange_stop_mc_bmrw->GetParameter(0); //
	double err_m_range_stop_mc_bmrw=fit_prange_stop_mc_bmrw->GetParError(0);
	double s_range_stop_mc_bmrw=fit_prange_stop_mc_bmrw->GetParameter(1); //
	double err_s_range_stop_mc_bmrw=fit_prange_stop_mc_bmrw->GetParError(1);

	//[3]
	TF1* fit_pcalo_stop_data; fit_pcalo_stop_data=VFit(pcalo_stop_data, 1);
	fit_pcalo_stop_data->SetName("fit_pcalo_stop_data");
	fit_pcalo_stop_data->SetLineStyle(2);
	double m_calo_stop_data=fit_pcalo_stop_data->GetParameter(0); //
	double err_m_calo_stop_data=fit_pcalo_stop_data->GetParError(0);
	double s_calo_stop_data=fit_pcalo_stop_data->GetParameter(1); //
	double err_s_calo_stop_data=fit_pcalo_stop_data->GetParError(1);

	//[4]
	TF1* fit_pcalo_stop_mc; fit_pcalo_stop_mc=VFit(pcalo_stop_mc, 2);
	fit_pcalo_stop_mc->SetName("fit_pcalo_stop_mc");
	fit_pcalo_stop_mc->SetLineStyle(2);
	double m_calo_stop_mc=fit_pcalo_stop_mc->GetParameter(0); //
	double err_m_calo_stop_mc=fit_pcalo_stop_mc->GetParError(0);
	double s_calo_stop_mc=fit_pcalo_stop_mc->GetParameter(1); //
	double err_s_calo_stop_mc=fit_pcalo_stop_mc->GetParError(1);

	//[4_bmrw]
	TF1* fit_pcalo_stop_mc_bmrw; fit_pcalo_stop_mc_bmrw=VFit(pcalo_stop_mc_bmrw, 4);
	fit_pcalo_stop_mc_bmrw->SetName("fit_pcalo_stop_mc_bmrw");
	fit_pcalo_stop_mc_bmrw->SetLineStyle(2);
	double m_calo_stop_mc_bmrw=fit_pcalo_stop_mc_bmrw->GetParameter(0); //
	double err_m_calo_stop_mc_bmrw=fit_pcalo_stop_mc_bmrw->GetParError(0);
	double s_calo_stop_mc_bmrw=fit_pcalo_stop_mc_bmrw->GetParameter(1); //
	double err_s_calo_stop_mc_bmrw=fit_pcalo_stop_mc_bmrw->GetParError(1);


	//[5]
	TF1* fit_pff_stop_mc; fit_pff_stop_mc=VFit(pff_stop_mc, 7);
	fit_pff_stop_mc->SetName("fit_pff_stop_mc");
	fit_pff_stop_mc->SetLineStyle(2);
	double m_ff_stop_mc=fit_pff_stop_mc->GetParameter(0); //
	double err_m_ff_stop_mc=fit_pff_stop_mc->GetParError(0);
	double s_ff_stop_mc=fit_pff_stop_mc->GetParameter(1); //
	double err_s_ff_stop_mc=fit_pff_stop_mc->GetParError(1);

	//[5_bmrw]
	TF1* fit_pff_stop_mc_bmrw; fit_pff_stop_mc_bmrw=VFit(pff_stop_mc_bmrw, 3);
	fit_pff_stop_mc_bmrw->SetName("fit_pff_stop_mc_bmrw");
	fit_pff_stop_mc_bmrw->SetLineStyle(2);
	double m_ff_stop_mc_bmrw=fit_pff_stop_mc_bmrw->GetParameter(0); //
	double err_m_ff_stop_mc_bmrw=fit_pff_stop_mc_bmrw->GetParError(0);
	double s_ff_stop_mc_bmrw=fit_pff_stop_mc_bmrw->GetParameter(1); //
	double err_s_ff_stop_mc_bmrw=fit_pff_stop_mc_bmrw->GetParError(1);


	TCanvas *c1x=new TCanvas("c1x","",1200,900);
	c1x->Divide(1,1);
	c1x->cd(1);
	TH2D* f2d_px=new TH2D("f2d_px","", 100, 600, 1300, 600, 0, prange_stop_mc_bmrw->GetBinContent(prange_stop_mc_bmrw->GetMaximumBin())+150); //
	f2d_px->SetTitle("Stopping Protons; Proton Momentum [MeV/c];");
	f2d_px->GetXaxis()->CenterTitle();
	prange_stop_data->SetLineColor(1);
	pcalo_stop_data->SetLineColor(1);
	pcalo_stop_data->SetMarkerColor(1);
	pcalo_stop_mc->SetLineColor(2);
	pcalo_stop_mc_bmrw->SetLineColor(4);
	pff_stop_mc->SetLineColor(7);
	pff_stop_mc_bmrw->SetLineColor(3);
	prange_stop_mc->SetLineColor(2);
	prange_stop_mc_bmrw->SetLineColor(4);
	
	f2d_px->Draw();
	prange_stop_mc->Draw("hist same");
	//pcalo_stop_data->Draw("ep same");
	//pcalo_stop_mc->Draw("hist same");
	pff_stop_mc->Draw("hist same");
	prange_stop_mc_bmrw->Draw("hist same");
	pff_stop_mc_bmrw->Draw("hist same");

	fit_prange_stop_data->Draw("same");
	fit_prange_stop_mc->Draw("same");
	fit_prange_stop_mc_bmrw->Draw("same");
	fit_pff_stop_mc->Draw("same");
	fit_pff_stop_mc_bmrw->Draw("same");
	prange_stop_data->Draw("ep same");

	//p0_stop_mc->SetLineColor(3);
	//p0_stop_mc->Scale(n_pbeam_stop_data/n_p0_stop_mc);
	//p0_stop_mc->Draw("hist same");

	TLegend *leg0x = new TLegend(0.14,0.65,.9,0.9);
	leg0x->SetFillStyle(0);
	leg0x->AddEntry(prange_stop_data, Form("Data(range): #mu=%.1f#pm%.1f MeV/c, #sigma=%.1f#pm%.1f MeV/c",m_range_stop_data,err_m_range_stop_data,s_range_stop_data,err_s_range_stop_data), "ep");
	leg0x->AddEntry(pff_stop_mc, Form("MC(FF-truth): #mu=%.1f#pm%.1f MeV/c, #sigma=%.1f#pm%.1f MeV/c",m_ff_stop_mc,err_m_ff_stop_mc,s_ff_stop_mc,err_s_ff_stop_mc), "l");
	leg0x->AddEntry(pbeam_stop_mc, Form("MC(range): #mu=%.1f#pm%.1f MeV/c, #sigma=%.1f#pm%.1f MeV/c",m_range_stop_mc,err_m_range_stop_mc,s_range_stop_mc,err_s_range_stop_mc), "l");
	leg0x->AddEntry(pff_stop_mc_bmrw, Form("MC(FF-truth)+bmrw: #mu=%.1f#pm%.1f MeV/c, #sigma=%.1f#pm%.1f MeV/c",m_ff_stop_mc_bmrw,err_m_ff_stop_mc_bmrw,s_ff_stop_mc_bmrw,err_s_ff_stop_mc_bmrw), "l");
	leg0x->AddEntry(pbeam_stop_mc_bmrw, Form("MC(range)+bmrw: #mu=%.1f#pm%.1f MeV/c, #sigma=%.1f#pm%.1f MeV/c",m_range_stop_mc_bmrw,err_m_range_stop_mc_bmrw,s_range_stop_mc_bmrw,err_s_range_stop_mc_bmrw), "l");

	//leg0x->AddEntry(pcalo_stop_data, "MC(calo)", "l");
	leg0x->Draw();

	c1x->Print(Form("%s/mom_range_bmrw.eps",outpath.Data()));

	c1x->Draw();
	f2d_px->Draw();
	pcalo_stop_mc->Draw("hist same");
	pcalo_stop_mc_bmrw->Draw("hist same");
	//pff_stop_mc->Draw("hist same");
	//pff_stop_mc_bmrw->Draw("hist same");
	pcalo_stop_data->Draw("ep same");

	fit_pcalo_stop_data->Draw("same");
	fit_pcalo_stop_mc->Draw("same");
	fit_pcalo_stop_mc_bmrw->Draw("same");
	//fit_pff_stop_mc->Draw("same");
	//fit_pff_stop_mc_bmrw->Draw("same");

	TLegend *leg0xy = new TLegend(0.14,0.65,.8,0.85);
	leg0xy->SetFillStyle(0);
	leg0xy->AddEntry(pcalo_stop_data, Form("Data(calo): #mu=%.1f#pm%.1f MeV/c, #sigma=%.1f#pm%.1f MeV/c",m_calo_stop_data,err_m_calo_stop_data,s_calo_stop_data,err_s_calo_stop_data), "ep");
	//leg0xy->AddEntry(pff_stop_mc, "MC (FF-truth)", "l");
	//leg0x->AddEntry(pff_stop_mc_bmrw, "MC (FF-truth)+bmrw", "l");
	leg0xy->AddEntry(pcalo_stop_mc, Form("MC(calo): #mu=%.1f#pm%.1f MeV/c, #sigma=%.1f#pm%.1f MeV/c",m_calo_stop_mc,err_m_calo_stop_mc,s_calo_stop_mc,err_s_calo_stop_mc), "l");
	leg0xy->AddEntry(pcalo_stop_mc_bmrw, Form("MC(calo)+bmrw: #mu=%.1f#pm%.1f MeV/c, #sigma=%.1f#pm%.1f MeV/c",m_calo_stop_mc_bmrw,err_m_calo_stop_mc_bmrw,s_calo_stop_mc_bmrw,err_s_calo_stop_mc_bmrw), "l");
	leg0xy->Draw();
	c1x->Print(Form("%s/mom_calo_bmrw.eps",outpath.Data()));



	//print out the fitting result ------------------------------------------------------------------------------------------------------//
	cout<<"================================================"<<endl;
	cout<<"      mu            |  sigma  "<<endl;
	cout<<"------------------------------------------------"<<endl;
	cout<<"initial:"<<endl;
	cout<<"data(spec):"<<m_stop_data<<"+-"<<err_m_stop_data<<" | "<<s_stop_data<<" +-"<<err_s_stop_data<<endl;
	cout<<"MC(spec):"<<m_stop_mc<<"+-"<<err_m_stop_mc<<" | "<<s_stop_mc<<" +-"<<err_s_stop_mc<<endl;
	cout<<"MC(truth):"<<m0_stop_mc<<"+-"<<err_m0_stop_mc<<" | "<<s0_stop_mc<<" +-"<<err_s0_stop_mc<<endl;
	cout<<"------------------------------------------------"<<endl;
	cout<<"end (and ff):"<<endl;
	cout<<"data(range):"<<m_range_stop_data<<"+-"<<err_m_range_stop_data<<" | "<<s_range_stop_data<<" +-"<<err_s_range_stop_data<<endl;
	cout<<"mc(range):"<<m_range_stop_mc<<"+-"<<err_m_range_stop_mc<<" | "<<s_range_stop_mc<<" +-"<<err_s_range_stop_mc<<endl;
	cout<<"data(calo):"<<m_calo_stop_data<<"+-"<<err_m_calo_stop_data<<" | "<<s_calo_stop_data<<" +-"<<err_s_calo_stop_data<<endl;
	cout<<"mc(calo):"<<m_calo_stop_mc<<"+-"<<err_m_calo_stop_mc<<" | "<<s_calo_stop_mc<<" +-"<<err_s_calo_stop_mc<<endl;
	cout<<"mc(ff):"<<m_ff_stop_mc<<"+-"<<err_m_ff_stop_mc<<" | "<<s_ff_stop_mc<<" +-"<<err_s_ff_stop_mc<<endl;
	cout<<"------------------------------------------------"<<endl;
	cout<<"mc_bmrw(range):"<<m_range_stop_mc_bmrw<<"+-"<<err_m_range_stop_mc_bmrw<<" | "<<s_range_stop_mc_bmrw<<" +-"<<err_s_range_stop_mc_bmrw<<endl;
	cout<<"mc_bmrw(calo):"<<m_calo_stop_mc_bmrw<<"+-"<<err_m_calo_stop_mc_bmrw<<" | "<<s_calo_stop_mc_bmrw<<" +-"<<err_s_calo_stop_mc_bmrw<<endl;
	cout<<"mc_bmrw(ff):"<<m_ff_stop_mc_bmrw<<"+-"<<err_m_ff_stop_mc_bmrw<<" | "<<s_ff_stop_mc_bmrw<<" +-"<<err_s_ff_stop_mc_bmrw<<endl;


	//summarize the rsult -------------------//
	//data
	vector<double> m_data;
	vector<double> err_m_data;
	vector<double> s_data;
	vector<double> err_s_data;
	double err_m_sys_stop_data=sqrt(0.01*0.01+0.008*0.008); //1% for B-field & 0.8% for fiber-position
	double err_m_all_stop_data=sqrt(pow(err_m_sys_stop_data,2)+pow(err_m_stop_data,2));
	double err_m_all_range_stop_data=sqrt(pow(err_m_sys_stop_data,2)+pow(err_m_range_stop_data,2));
	double err_m_all_calo_stop_data=sqrt(pow(err_m_sys_stop_data,2)+pow(err_m_calo_stop_data,2));

	m_data.push_back(m_stop_data); err_m_data.push_back(err_m_all_stop_data); s_data.push_back(s_stop_data); err_s_data.push_back(err_s_stop_data);
	m_data.push_back(m_range_stop_data); err_m_data.push_back(err_m_all_range_stop_data); s_data.push_back(s_range_stop_data); err_s_data.push_back(err_s_range_stop_data);
	m_data.push_back(m_calo_stop_data); err_m_data.push_back(err_m_all_calo_stop_data); s_data.push_back(s_calo_stop_data); err_s_data.push_back(err_s_calo_stop_data);


	//mc
	vector<double> m_mc;
	vector<double> err_m_mc;
	vector<double> s_mc;
	vector<double> err_s_mc;
	m_mc.push_back(m0_stop_mc); err_m_mc.push_back(err_m0_stop_mc); s_mc.push_back(s0_stop_mc); err_s_mc.push_back(err_s0_stop_mc);
	m_mc.push_back(m_stop_mc); err_m_mc.push_back(err_m_stop_mc); s_mc.push_back(s_stop_mc); err_s_mc.push_back(err_s_stop_mc);
	
	m_mc.push_back(m_ff_stop_mc); err_m_mc.push_back(err_m_ff_stop_mc); s_mc.push_back(s_ff_stop_mc); err_s_mc.push_back(err_s_ff_stop_mc);
	m_mc.push_back(m_range_stop_mc); err_m_mc.push_back(err_m_range_stop_mc); s_mc.push_back(s_range_stop_mc); err_s_mc.push_back(err_s_range_stop_mc);
	m_mc.push_back(m_calo_stop_mc); err_m_mc.push_back(err_m_calo_stop_mc); s_mc.push_back(s_calo_stop_mc); err_s_mc.push_back(err_s_calo_stop_mc);

	TGraphErrors *ms_data=new TGraphErrors(m_data.size(), &m_data.at(0), &s_data.at(0), &err_m_data.at(0), &err_s_data.at(0));
	TGraphErrors *ms_mc=new TGraphErrors(m_mc.size(), &m_mc.at(0), &s_mc.at(0), &err_m_mc.at(0), &err_s_mc.at(0));
	ms_mc->SetMarkerColor(2);
	ms_mc->SetLineColor(2);


	//mc_bmrw
	vector<double> m_mc_bmrw;
	vector<double> err_m_mc_bmrw;
	vector<double> s_mc_bmrw;
	vector<double> err_s_mc_bmrw;
	m_mc_bmrw.push_back(m0_stop_mc_bmrw); err_m_mc_bmrw.push_back(err_m0_stop_mc_bmrw); s_mc_bmrw.push_back(s0_stop_mc_bmrw); err_s_mc_bmrw.push_back(err_s0_stop_mc_bmrw);
	m_mc_bmrw.push_back(m_stop_mc_bmrw); err_m_mc_bmrw.push_back(err_m_stop_mc_bmrw); s_mc_bmrw.push_back(s_stop_mc_bmrw); err_s_mc_bmrw.push_back(err_s_stop_mc_bmrw);
	
	m_mc_bmrw.push_back(m_ff_stop_mc_bmrw); err_m_mc_bmrw.push_back(err_m_ff_stop_mc_bmrw); s_mc_bmrw.push_back(s_ff_stop_mc_bmrw); err_s_mc_bmrw.push_back(err_s_ff_stop_mc_bmrw);
	m_mc_bmrw.push_back(m_range_stop_mc_bmrw); err_m_mc_bmrw.push_back(err_m_range_stop_mc_bmrw); s_mc_bmrw.push_back(s_range_stop_mc_bmrw); err_s_mc_bmrw.push_back(err_s_range_stop_mc_bmrw);
	m_mc_bmrw.push_back(m_calo_stop_mc_bmrw); err_m_mc_bmrw.push_back(err_m_calo_stop_mc_bmrw); s_mc_bmrw.push_back(s_calo_stop_mc_bmrw); err_s_mc_bmrw.push_back(err_s_calo_stop_mc_bmrw);

	TGraphErrors *ms_mc_bmrw=new TGraphErrors(m_mc_bmrw.size(), &m_mc_bmrw.at(0), &s_mc_bmrw.at(0), &err_m_mc_bmrw.at(0), &err_s_mc_bmrw.at(0));
	ms_mc->SetMarkerColor(2);
	ms_mc->SetLineColor(2);
	ms_mc_bmrw->SetMarkerColor(3);
	ms_mc_bmrw->SetLineColor(3);

	TCanvas *cms=new TCanvas("cms","",1200,900);
	cms->Divide(1,1);
	cms->cd(1);
	float emin=910.;
	float emax=1030.;
	float smin=40;
	float smax=120;
	TH2D* frame2d2=new TH2D("frame2d2","", 200, emin, emax, 42, smin, smax);
	frame2d2->SetTitle(";#mu [MeV/c];#sigma [MeV/c]");
	frame2d2->GetXaxis()->CenterTitle();
	frame2d2->Draw();
	ms_data->Draw("p same");
	ms_mc->Draw("p same");
	ms_mc_bmrw->Draw("p same");


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

        TF1 *fit_mc_bmrw = new TF1("fit_mc_bmrw", "[0]+[1]*x", emin, emax);
        fit_mc_bmrw->SetLineColor(3);
        fit_mc_bmrw->SetLineStyle(3);
        ms_mc_bmrw->Fit(fit_mc_bmrw,"remn");
        fit_mc_bmrw->Draw("same");


	TLegend *leg1 = new TLegend(0.14,0.65,.6,0.85);
	//leg1->SetFillStyle(0);
	leg1->AddEntry(ms_data, Form("Data: #sigma=%.2f%.2f*#mu",fit_data->GetParameter(0),fit_data->GetParameter(1)), "ep");
	leg1->AddEntry(ms_mc, Form("MC: #sigma=%.2f%.2f*#mu",fit_mc->GetParameter(0),fit_mc->GetParameter(1)), "ep");
	leg1->AddEntry(ms_mc_bmrw, Form("MC(bmrw): #sigma=%.2f%.2f*#mu",fit_mc_bmrw->GetParameter(0),fit_mc_bmrw->GetParameter(1)), "ep");
	leg1->Draw();

        //pDUNE Logo
        TLatex **txt_pdune1=new TLatex*[1];
        txt_pdune1[0]=new TLatex(emin+.002, emax+0.6, Form("#bf{DUNE:ProtoDUNE-SP}"));
        txt_pdune1[0]->SetTextColor(1);
        txt_pdune1[0]->Draw();
        //
        //Beam Logo
        TLatex **txt_p1=new TLatex*[1];
        txt_p1[0]=new TLatex(emax-200,emax+0.6, Form("Protons (1 GeV/c)"));
        txt_p1[0]->SetTextColor(1);
        txt_p1[0]->SetTextSize(0.05);
        txt_p1[0]->Draw();

	cms->Print(Form("%s/mu_sigma_data_mc.eps",outpath.Data()));








}
