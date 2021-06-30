#include "/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4_MC_1GeV_reco1_sce_datadriven/SliceParams.h"

Double_t rr2dedx(Double_t* x,Double_t *par) {
        double a=17.;
        double b=-0.42;

        return a*pow(x[0],b);
}

void plot_KEs(TString fin, TString outpath){

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	//gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 

	gStyle->SetOptStat(0);

	TFile *f0 = TFile::Open(fin.Data());
	TH2D *KE_ff_recostop=(TH2D *)f0->Get("KE_ff_recostop");
	TH2D *KE_range_recostop=(TH2D *)f0->Get("KE_range_recostop");
	TH2D *KE_calo_recostop=(TH2D *)f0->Get("KE_calo_recostop");
	TH2D *KE_rrange_recostop=(TH2D *)f0->Get("KE_rrange_recostop");
	TH2D *KE_rrange2_recostop=(TH2D *)f0->Get("KE_rrange2_recostop"); //dedx vs rr using LV
	TH2D *KE_simide_recostop=(TH2D *)f0->Get("KE_simide_recostop");
	TH2D *rr_dedx_recostop=(TH2D *)f0->Get("rr_dedx_recostop");

	TH1D *dKE_range_ff_recostop=(TH1D *)f0->Get("dKE_range_ff_recostop");
	TH1D *dKE_calo_ff_recostop=(TH1D *)f0->Get("dKE_calo_ff_recostop");
	TH1D *dKE_rrange_ff_recostop=(TH1D *)f0->Get("dKE_rrange_ff_recostop"); //ub
	TH1D *dKE_rrange2_ff_recostop=(TH1D *)f0->Get("dKE_rrange2_ff_recostop"); //LV

	dKE_range_ff_recostop->GetXaxis()->SetTitle("KE_{range}-KE_{ff} [MeV]");
	dKE_calo_ff_recostop->GetXaxis()->SetTitle("KE_{calo}-KE_{ff} [MeV]");
	dKE_rrange_ff_recostop->GetXaxis()->SetTitle("KE_{uB}-KE_{ff} [MeV]");
	dKE_rrange2_ff_recostop->GetXaxis()->SetTitle("KE_{Landau-Vavilov}-KE_{ff} [MeV]");


	TFile *fdedx_rr=new TFile("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/dedx_rr/dedx_rr.root");
	TGraph *gr_predict_dedx_resrange=(TGraph *)fdedx_rr->Get("gr_predict_dedx_resrange");

	KE_ff_recostop->SetLineColor(1);  
	KE_range_recostop->SetLineColor(2);  
	KE_calo_recostop->SetLineColor(4);  
	KE_rrange_recostop->SetLineColor(3); 
	KE_rrange2_recostop->SetLineColor(14);  
 
	KE_simide_recostop->SetLineColor(7);  

	int n_rb=2;
	//KE_ff_recostop->Rebin(n_rb);
	//KE_range_recostop->Rebin(n_rb);
	//KE_calo_recostop->Rebin(n_rb);
	//KE_rrange_recostop->Rebin(n_rb);
	//KE_simide_recostop->Rebin(n_rb);
	


	TH2D* frame2d=new TH2D("frame2d","", 400, 200, 600, 400, 0, 400); //zend_2d
	//frame2d->GetXaxis()->SetTitle("Proton KE [MeV]");
	//frame2d->GetYaxis()->SetTitle("Counts");
	frame2d->SetTitle("Reco Stopping Protons;Proton KE [MeV];Counts");


	TCanvas *c0=new TCanvas("c0",""); 
	c0->Divide(1,1);
	c0->cd(1);
	frame2d->Draw();
	KE_ff_recostop->Draw("same");  
	KE_range_recostop->Draw("same");  
	KE_calo_recostop->Draw("same");  
	//KE_rrange_recostop->Draw("same");  
	//KE_rrange2_recostop->Draw("same");  

	//KE_simide_recostop->Draw("same");  

	TLegend *leg5 = new TLegend(0.15,0.65,0.75,0.85);
	leg5->SetFillStyle(0);
	leg5->AddEntry(KE_ff_recostop, "KE_{ff}", "l");
	leg5->AddEntry(KE_range_recostop, "KE_{range}", "l");
	leg5->AddEntry(KE_calo_recostop, "KE_{calo}", "l");
	//leg5->AddEntry(KE_rrange_recostop, "KE_{Microboone approx: dE/dx=17 #times rr^{-0.42}}", "l");
	//leg5->AddEntry(KE_rrange2_recostop, "KE_{Landau-Vavilov}", "l");
	leg5->Draw();

	
	c0->Print(Form("%s/KEs.eps",outpath.Data()));


	TCanvas *c1=new TCanvas("c1",""); 
	c1->Divide(1,1);
	c1->cd(1);
	frame2d->Draw();
	KE_ff_recostop->Draw("same");  
	KE_simide_recostop->Draw("same");  

	TLegend *leg6 = new TLegend(0.15,0.65,0.75,0.85);
	leg6->SetFillStyle(0);
	leg6->AddEntry(KE_ff_recostop, "KE_{ff}", "l");
	leg6->AddEntry(KE_simide_recostop, "KE_{simIDE}", "l");
	leg6->Draw();

	c1->Print(Form("%s/KE_truths.eps",outpath.Data()));

	TCanvas *c2=new TCanvas("c2",""); 
	c2->Divide(1,1);
	c2->cd(1);
	rr_dedx_recostop->SetTitle("Reco Stop Cut; Residual Range [cm]; dE/dx [MeV/cm]");
	rr_dedx_recostop->Draw("colz");
	gr_predict_dedx_resrange->SetMarkerSize(0.3);
	gr_predict_dedx_resrange->SetLineColor(2);
	gr_predict_dedx_resrange->SetLineWidth(2);
	gr_predict_dedx_resrange->SetLineStyle(7);
	
	gr_predict_dedx_resrange->SetMarkerColor(2);
	gr_predict_dedx_resrange->Draw("p same");

	TF1 *fun=new TF1(Form("rr2dedx"),rr2dedx,0,120,0);
	fun->SetLineColor(1);
	fun->SetLineStyle(2);
	//fun->Draw("same");

	TLegend *legz = new TLegend(0.32,0.65,0.85,0.85);
	//legz->SetFillStyle(0);
	legz->AddEntry(gr_predict_dedx_resrange, "Landau-Vavilov", "l");
	//legz->AddEntry(fun, "Microboone approx: dE/dx=17 #times rr^{-0.42}", "l");
	legz->Draw();

	c2->Print(Form("%s/dedx_rr_recostop.eps",outpath.Data()));
	
	TCanvas *c3=new TCanvas("c3",""); 
	c3->Divide(1,1);
	c3->cd(1);
	TH2D* frame2d_zoom1=new TH2D("frame2d_zoom1","", 80, 40, 120, 8, 0, 8);
	//frame2d_zoom1->SetTitle("Reco Stop Cut; Residual Range [cm]; dE/dx [MeV/cm]")
	frame2d_zoom1->GetXaxis()->SetTitle("Residual Range[cm]");
	frame2d_zoom1->GetYaxis()->SetTitle("dE/dx[MeV/cm]");
	frame2d_zoom1->Draw();
	rr_dedx_recostop->Draw("colz same");
	gr_predict_dedx_resrange->Draw("cp same");
	//fun->Draw("same");
	legz->Draw();

	c3->Print(Form("%s/dedx_rr_recostop_zoom1.eps",outpath.Data()));
	

	TCanvas *c4=new TCanvas("c4",""); 
	c4->Divide(1,1);
	c4->cd(1);
	TH2D* frame2d_zoom2=new TH2D("frame2d_zoom2","", 40, 0, 40, 30, 0, 30);
	frame2d_zoom2->GetXaxis()->SetTitle("Residual Range[cm]");
	frame2d_zoom2->GetYaxis()->SetTitle("dE/dx[MeV/cm]");
	frame2d_zoom2->Draw();
	rr_dedx_recostop->Draw("colz same");
	gr_predict_dedx_resrange->Draw("cp same");
	//fun->Draw("same");
	legz->Draw();

	c4->Print(Form("%s/dedx_rr_recostop_zoom2.eps",outpath.Data()));
	

	TCanvas *c5=new TCanvas("c5",""); 
	c5->Divide(1,1);
	c5->cd(1);
	TH2D* frame2d_zoomx=new TH2D("frame2d_zoomx","", 120, 0, 120, 30, 0, 30);
	frame2d_zoomx->GetXaxis()->SetTitle("Residual Range[cm]");
	frame2d_zoomx->GetYaxis()->SetTitle("dE/dx[MeV/cm]");
	frame2d_zoomx->Draw();
	gr_predict_dedx_resrange->Draw("p same");
	//fun->Draw("same");

	legz->Draw();

	c5->Print(Form("%s/dedx_rr_theories.eps",outpath.Data()));



	gStyle->SetOptStat(1);
	TCanvas *c6=new TCanvas("c6",""); 
	c6->Divide(1,1);
	c6->cd(1);
	dKE_range_ff_recostop->Draw("hist");
	c6->Print(Form("%s/dke_range_ff.eps",outpath.Data()));
	

	TCanvas *c7=new TCanvas("c7",""); 
	c7->Divide(1,1);
	c7->cd(1);
	dKE_calo_ff_recostop->Draw("hist");
	c7->Print(Form("%s/dke_calo_ff.eps",outpath.Data()));

	TCanvas *c8=new TCanvas("c8",""); 
	c8->Divide(1,1);
	c8->cd(1);
	dKE_rrange_ff_recostop->Draw("hist");
	c8->Print(Form("%s/dke_ub_ff.eps",outpath.Data()));

	TCanvas *c9=new TCanvas("c9",""); 
	c9->Divide(1,1);
	c9->cd(1);
	dKE_rrange2_ff_recostop->Draw("hist");
	c9->Print(Form("%s/dke_ub_lv.eps",outpath.Data()));

	TCanvas *c10=new TCanvas("c10",""); 
	c10->Divide(1,1);
	c10->cd(1);
	gStyle->SetOptStat(0);
	dKE_range_ff_recostop->SetLineColor(2);
	dKE_calo_ff_recostop->SetLineColor(4); 
	dKE_rrange_ff_recostop->SetLineColor(3); 
	dKE_rrange2_ff_recostop->SetLineColor(14);

	TH2D* frame2d_zoomxx=new TH2D("frame2d_zoomxx","", 200, -100, 100, 800, 0, 800);
	frame2d_zoomxx->GetXaxis()->SetTitle("#DeltaKE [MeV]");
	frame2d_zoomxx->Draw();
	dKE_range_ff_recostop->Draw("hist same"); 
	dKE_calo_ff_recostop->Draw("hist same"); 
	dKE_rrange_ff_recostop->Draw("hist same"); 
	dKE_rrange2_ff_recostop->Draw("hist same"); 

	TLegend *leg55 = new TLegend(0.5,0.65,0.87,0.88);
	leg55->SetFillStyle(0);
	leg55->AddEntry(dKE_range_ff_recostop, Form("KE_{range}-KE_{ff}: %.2f MeV",dKE_range_ff_recostop->GetMean()),"l");
	leg55->AddEntry(dKE_calo_ff_recostop,Form("KE_{calo}-KE_{ff}: %.2f MeV",dKE_calo_ff_recostop->GetMean()),"l");
	leg55->AddEntry(dKE_rrange_ff_recostop,Form("KE_{uB}-KE_{ff}: %.2f MeV",dKE_rrange_ff_recostop->GetMean()),"l");
	leg55->AddEntry(dKE_rrange2_ff_recostop,Form("KE_{Landau-Vavilov}-KE_{ff}: %.2f MeV",dKE_rrange2_ff_recostop->GetMean()),"l");
	leg55->Draw();

	c10->Print(Form("%s/dke_summary.eps",outpath.Data()));
	

}

