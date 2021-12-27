//#include "../headers/BasicAnaFunc.h"
#include <vector>
#include <iostream>
#include "TGraphAsymmErrors.h"
#include "../headers/SliceParams.h"
#include "../headers/TemplateFitter.h"

void plot_vtx_eff() {

	TString file_mc=Form("../mc_vtx_rangereco.root");
	TFile *fmc = TFile::Open(file_mc.Data());

	TH2D *h2d=(TH2D *)fmc->Get(Form("h2d_trklen_dr_RecoInel"));
	TH2D *h2d_el=(TH2D *)fmc->Get(Form("h2d_trklen_dr_RecoInel_el"));
	TH2D *h2d_inel=(TH2D *)fmc->Get(Form("h2d_trklen_dr_RecoInel_inel"));
	TH2D *h2d_misidp=(TH2D *)fmc->Get(Form("h2d_trklen_dr_RecoInel_misidp"));

	TH2D *h2d_2060=(TH2D *)fmc->Get(Form("h2d_trklen_dr_RecoInel_2060_inel"));
	TH2D *h2d_0010=(TH2D *)fmc->Get(Form("h2d_trklen_dr_RecoInel_0010_inel"));

	TH2D *h2d_good=(TH2D *)fmc->Get(Form("h2d_trklen_dr_RecoInel_GOOD"));
	TH2D *h2d_bad=(TH2D *)fmc->Get(Form("h2d_trklen_dr_RecoInel_BAD"));
	TH2D *h2d_ugly=(TH2D *)fmc->Get(Form("h2d_trklen_dr_RecoInel_UGLY"));

	TH2D *h2d_true_inel=(TH2D *)fmc->Get(Form("h2d_truetrklen_dr_RecoInel_inel"));


	//figout
	TString fig_out="plots_vtx/recolen_dr_recoinel.eps";
	TString fig_out_inel="plots_vtx/recolen_dr_recoinel_inel.eps";
	TString fig_out_el="plots_vtx/recolen_dr_recoinel_el.eps";
	TString fig_out_misidp="plots_vtx/recolen_dr_recoinel_misidp.eps";

	TString fig_out_inel_zoom="plots_vtx/recolen_dr_recoinel_inel_zoom.eps";
	TString fig_out_inel_zoom1d="plots_vtx/dr_recoinel_inel_zoom.eps.eps";
	TString fig_out_inel_zoom1d_2="plots_vtx/dr_recoinel_inel_zoom2.eps.eps";

	TString fig_out_all="plots_vtx/recolen_good_bad_ugly.eps";


	TString fig_out_misidp_zoom="recolen_dr_BQ_misidp_zoom.eps";
	TString fig_out_misidp_zoomy="recolen_dr_BQ_misidp_zoomy.eps";

	TString fig_out_el_zoom="recolen_dr_BQ_el_zoom.eps";
	TString fig_out_el_zoomy="recolen_dr_BQ_el_zoomy.eps";

       	//int n_b=300;
        //double b_min=0;
        //double b_max=150;
       	int n_x=h2d->GetNbinsX();
        double x_min=h2d->GetXaxis()->GetXmin();
        double x_max=h2d->GetXaxis()->GetXmax();
	double bin_x=(x_max-x_min)/(double)n_x;
	cout<<"n_x:"<<n_x<<endl;
	cout<<"x_min:"<<x_min<<" - x_max:"<<x_max<<endl;

        //int n_dr=600;
        //double dr_min=-150;
        //double dr_max=150;
        int n_y=h2d->GetNbinsY();
        double y_min=h2d->GetYaxis()->GetXmin();
        double y_max=h2d->GetYaxis()->GetXmax();
	double bin_y=(y_max-y_min)/(double)n_y;
	cout<<"\nn_y:"<<n_y<<endl;
	cout<<"y_min:"<<y_min<<" - y_max:"<<y_max<<endl;

	//config -------------------------------------------------------//
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleYOffset(1.2);
	//gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	gStyle->SetOptStat(0);
	gStyle->SetPadLeftMargin(0.15);
        gStyle->SetFrameLineWidth(1);

	//gStyle->SetPalette(60);
	//gStyle->SetPalette(kLightTemperature);
	//gStyle->SetPalette(kTemperatureMap);
	//gStyle->SetPalette(kColorPrintableOnGrey);
	//gStyle->SetPalette(kBlueGreenYellow);
	//gStyle->SetPalette(kFruitPunch);
	//gStyle->SetPalette(kGreenPink);
	//gStyle->SetPalette(kPastel);
	//gStyle->SetPalette(kRedBlue);
	//gStyle->SetPalette(kCool);
	//gStyle->SetPalette(kWaterMelon);
	//gStyle->SetPalette(kBlackBody);
	//gStyle->SetPadTopMargin(0.2);
	//gStyle->SetPadRightMargin(0.17);
	//config -------------------------------------------------------//

	//inel ---------------------------------------------------------------------------------------//
	TCanvas *c_inel = new TCanvas("c_inel", "c_inel", 1200, 900);
	c_inel->Divide(1,1);
	c_inel->cd(1);
	h2d->SetTitle("Reco Inel. Cut; Reco Track Length [cm]; #Delta L");
	h2d->GetYaxis()->SetRangeUser(-120,100);
	h2d->Draw("colz");
	c_inel->Print(fig_out.Data());


	h2d_inel->SetTitle("Reco Inel. Cut [True Inel.]; Reco Track Length [cm]; #Delta L");
	h2d_inel->GetYaxis()->SetRangeUser(-120,100);
	h2d_inel->Draw("colz");
	c_inel->Print(fig_out_inel.Data());

	h2d_el->SetTitle("Reco Inel. Cut [True El.]; Reco Track Length [cm]; #Delta L");
	h2d_el->GetYaxis()->SetRangeUser(-120,100);
	h2d_el->Draw("colz");
	c_inel->Print(fig_out_el.Data());

	h2d_misidp->SetTitle("Reco Inel. Cut [MisID:P]; Reco Track Length [cm]; #Delta L");
	h2d_misidp->GetYaxis()->SetRangeUser(-120,100);
	h2d_misidp->Draw("colz");
	c_inel->Print(fig_out_misidp.Data());

	h2d_inel->GetYaxis()->SetRangeUser(-10,10);
	h2d_inel->GetXaxis()->SetRangeUser(0,110);
	h2d_inel->Draw("colz");
	c_inel->Print(fig_out_inel_zoom.Data());



	TCanvas *c1D = new TCanvas("c1D", "c1D", 1200, 900);
	c1D->Divide(1,1);
	c1D->cd(1);
	//c1D->cd(1)->SetLogy();

	TH2D *f2d_y=new TH2D("","",10,-100,20, 100,0.001,.2);
	f2d_y->SetTitle(Form("Reco Inel. Cut [True Inel.]; #Delta L [cm]"));
	f2d_y->Draw();
	
	TH1D *h1d_2060=h2d_2060->ProjectionY("h1d_2060");
	h1d_2060->SetLineColor(4);
	h1d_2060->Scale(1./(double)h1d_2060->Integral());
	h1d_2060->Draw("hist same");

	TH1D *h1d_0010=(TH1D *)h2d_0010->ProjectionY("h1d_0010");
	h1d_0010->SetLineColor(2);
	h1d_0010->Scale(1./(double)h1d_0010->Integral());
	h1d_0010->Draw("hist same");

        TLegend *leg0 = new TLegend(0.2,0.65,.5,0.88);
        leg0->SetFillStyle(0);
        leg0->AddEntry(h1d_2060, Form("20-60 cm"), "l");
        leg0->AddEntry(h1d_0010, Form("0-10 cm"), "l");
        leg0->Draw();

	c1D->Print(fig_out_inel_zoom1d.Data());

	TH2D *f2d_y2=new TH2D("","",10,-20,15, 100,0.001,.2);
	f2d_y2->SetTitle(Form("Reco Inel. Cut [True Inel.]; #Delta L [cm]"));
	f2d_y2->Draw();
	h1d_2060->Draw("hist same");
	h1d_0010->Draw("hist same");
        leg0->Draw();
	c1D->Print(fig_out_inel_zoom1d_2.Data());
	

	TCanvas *c1Dx = new TCanvas("c1Dx", "c1Dx", 1200, 900);
	c1Dx->Divide(1,1);
	c1Dx->cd(1);

	TH2D *f2dx=new TH2D("","",140,0,120,100,0,1600);
	f2dx->SetTitle(Form("Reco Inel. Cut [True Inel.]; Reco Track Length [cm]; Counts"));
	f2dx->Draw();

	TH1D *h1d=h2d->ProjectionX("h1d");
	TH1D *h1d_good=h2d_good->ProjectionX("h1d_good");
	TH1D *h1d_bad=h2d_bad->ProjectionX("h1d_bad");
	TH1D *h1d_ugly=h2d_ugly->ProjectionX("h1d_ugly");
	TH1D *h1d_true_inel=h2d_true_inel->ProjectionX("h1d_true_inel");

	int n_rb=4;

	h1d->Rebin(n_rb);
	h1d_good->Rebin(n_rb);
	h1d_bad->Rebin(n_rb);
	h1d_ugly->Rebin(n_rb);
	h1d_true_inel->Rebin(n_rb);

	
	h1d->SetLineColor(1);
	h1d_good->SetLineColor(2);
	h1d_bad->SetLineColor(4);
	h1d_ugly->SetLineColor(3);
	h1d_true_inel->SetLineColor(15);

	TH1D *h1d_all=(TH1D*)h1d_good->Clone("h1d_all");
	h1d_all->Add(h1d_bad);
	h1d_all->Add(h1d_ugly);
	h1d_all->SetLineColor(1);

	//h1d->Draw("same");
	h1d_all->Draw("same");
	h1d_good->Draw("same");
	h1d_bad->Draw("same");
	h1d_ugly->Draw("same");
	h1d_true_inel->Scale((double)h1d_good->Integral()/(double)h1d_true_inel->Integral());
	h1d_true_inel->Draw("hist same");

        TLegend *leg = new TLegend(0.45,0.65,.85,0.88);
        leg->SetFillStyle(0);
        leg->AddEntry(h1d, Form("True Inel."), "l");
        leg->AddEntry(h1d_good, Form("True Inel. with |#DeltaL|#leq3 cm"), "l");
        leg->AddEntry(h1d_bad, Form("True Inel. with #DeltaL<-3 cm"), "l");
        leg->AddEntry(h1d_ugly, Form("True Inel. with #DeltaL>3 cm"), "l");
        leg->AddEntry(h1d_true_inel, Form("True Inel. with true track length"), "l");
        leg->Draw();

	c1Dx->Print(fig_out_all.Data());





/*
	double xmin=0;
	double xmax=12;
	double ymin=-20;
	double ymax=10;

	int xmin_bin=(int)(xmin-x_min)/bin_x;
	int xmax_bin=(int)(xmax-x_min)/bin_x;	
	cout<<"xmin_bin-xmax_bin:"<<xmin_bin<<" - "<<xmax_bin<<endl;

	double xmin2=20;
	double xmax2=60;
	double ymin2=-20;
	double ymax2=10;

	int xmin_bin2=(int)(xmin2-x_min)/bin_x;
	int xmax_bin2=(int)(xmax2-x_min)/bin_x;	
	cout<<"xmin_bin2-xmax_bin2:"<<xmin_bin2<<" - "<<xmax_bin2<<endl;

	//inel ---------------------------------------------------------------------------------------//
	TCanvas *c_inel = new TCanvas("c_inel", "c_inel", 1200, 900);
	c_inel->Divide(1,1);
	c_inel->cd(1);
	h2d_inel->SetTitle("; Reco Track Length[cm]; #Delta L");
	h2d_inel->GetYaxis()->SetRangeUser(-100,50);
	h2d_inel->Draw("colz");

	//c_inel->Print(fig_out_inel.Data());
	//h2d_inel->GetXaxis()->SetRangeUser(xmin,xmax);
	//h2d_inel->GetYaxis()->SetRangeUser(ymin,ymax);
	//h2d_inel->Draw("colz");
	//c_inel->Print(fig_out_inel_zoom.Data());

	TCanvas *c_y = new TCanvas("c_y", "c_y", 1200, 900);
	c_y->Divide(1,1);
	c_y->cd(1);
	c_y->cd(1)->SetLogy();
	//TH2D *f2d_y=new TH2D("","",10,-20,10, 100,0,2000);
	TH2D *f2d_y=new TH2D("","",10,-100,10, 100,0,2000);
	f2d_y->SetTitle(Form("; Track Length(reco)-Track Length(truth) [cm]"));
	f2d_y->Draw();
	
	TH1D *h1d_inel=(TH1D *)h2d_inel->ProjectionY("h1d_inel", xmin_bin, xmax_bin, "");
	h1d_inel->SetName("h1d_inel");
	//h1d_inel->SetTitle(Form("; Track Length(reco)-Track Length(truth) [cm]"));
	h1d_inel->SetLineColor(2);
	h1d_inel->Draw("same");

	TH1D *h1d_inel2=(TH1D *)h2d_inel->ProjectionY("h1d_inel2", xmin_bin2, xmax_bin2, "");
	h1d_inel2->SetName("h1d_inel2");
	h1d_inel2->SetLineColor(4);
	h1d_inel2->Scale((double)h1d_inel->Integral()/(double)h1d_inel2->Integral());
	h1d_inel2->Draw("hist same");

        TLegend *leg0 = new TLegend(0.5,0.65,.85,0.88);
        leg0->SetFillStyle(0);
        leg0->AddEntry(h1d_inel, Form("%.0f-%.0f cm",xmin,xmax), "l");
        leg0->AddEntry(h1d_inel2, Form("%.0f-%.0f cm",xmin2,xmax2), "l");
        leg0->Draw();

	c_y->Print(fig_out_inel_zoomy.Data());

	//MisID:P ---------------------------------------------------------------------------------------//
	TCanvas *c_misidp = new TCanvas("c_misidp", "c_misidp", 1200, 900);
	c_misidp->Divide(1,1);
	c_misidp->cd(1);
	h2d_misidp->SetTitle("; Reco Track Length[cm]; #Delta L");
	h2d_misidp->Draw("colz");
	c_misidp->Print(fig_out_misidp.Data());
	h2d_misidp->GetXaxis()->SetRangeUser(xmin,xmax);
	h2d_misidp->GetYaxis()->SetRangeUser(ymin,ymax);
	h2d_misidp->Draw("colz");
	c_misidp->Print(fig_out_misidp_zoom.Data());

	TCanvas *c_y_misidp = new TCanvas("c_y_misidp", "c_y_misidp", 1200, 900);
	c_y_misidp->Divide(1,1);
	c_y_misidp->cd(1);
	c_y_misidp->cd(1)->SetLogy();
	TH2D *f2d_y_misidp=new TH2D("","",40,-20,10, 100,0,200);
	f2d_y_misidp->SetTitle(Form("; Track Length(reco)-Track Length(truth) [cm]"));
	f2d_y_misidp->Draw();
	
	TH1D *h1d_misidp=(TH1D *)h2d_misidp->ProjectionY("h1d_misidp", xmin_bin, xmax_bin, "");
	h1d_misidp->SetName("h1d_misidp");
	h1d_misidp->SetLineColor(2);
	h1d_misidp->Draw("same");

	TH1D *h1d_misidp2=(TH1D *)h2d_misidp->ProjectionY("h1d_misidp2", xmin_bin2, xmax_bin2, "");
	h1d_misidp2->SetName("h1d_misidp2");
	h1d_misidp2->SetLineColor(4);
	h1d_misidp2->Scale((double)h1d_misidp->Integral()/(double)h1d_misidp2->Integral());
	h1d_misidp2->Draw("hist same");

        TLegend *leg0_misidp = new TLegend(0.5,0.65,.85,0.88);
        leg0_misidp->SetFillStyle(0);
        leg0_misidp->AddEntry(h1d_misidp, Form("%.0f-%.0f cm",xmin,xmax), "l");
        leg0_misidp->AddEntry(h1d_misidp2, Form("%.0f-%.0f cm",xmin2,xmax2), "l");
        leg0_misidp->Draw();

	c_y_misidp->Print(fig_out_misidp_zoomy.Data());


	//El. ---------------------------------------------------------------------------------------//
	TCanvas *c_el = new TCanvas("c_el", "c_el", 1200, 900);
	c_el->Divide(1,1);
	c_el->cd(1);
	h2d_el->SetTitle("; Reco Track Length[cm]; #Delta L");
	h2d_el->Draw("colz");
	c_el->Print(fig_out_el.Data());
	h2d_el->GetXaxis()->SetRangeUser(xmin,xmax);
	h2d_el->GetYaxis()->SetRangeUser(ymin,ymax);
	h2d_el->Draw("colz");
	c_el->Print(fig_out_el_zoom.Data());

	TCanvas *c_y_el = new TCanvas("c_y_el", "c_y_el", 1200, 900);
	c_y_el->Divide(1,1);
	c_y_el->cd(1);
	c_y_el->cd(1)->SetLogy();
	TH2D *f2d_y_el=new TH2D("","",40,-20,10, 100,0,200);
	f2d_y_el->SetTitle(Form("; Track Length(reco)-Track Length(truth) [cm]"));
	f2d_y_el->Draw();
	
	TH1D *h1d_el=(TH1D *)h2d_el->ProjectionY("h1d_el", xmin_bin, xmax_bin, "");
	h1d_el->SetName("h1d_el");
	h1d_el->SetLineColor(2);
	h1d_el->Draw("same");

	TH1D *h1d_el2=(TH1D *)h2d_el->ProjectionY("h1d_el2", xmin_bin2, xmax_bin2, "");
	h1d_el2->SetName("h1d_el2");
	h1d_el2->SetLineColor(4);
	h1d_el2->Scale((double)h1d_el->Integral()/(double)h1d_el2->Integral());
	h1d_el2->Draw("hist same");

        TLegend *leg0_el = new TLegend(0.5,0.65,.85,0.88);
        leg0_el->SetFillStyle(0);
        leg0_el->AddEntry(h1d_el, Form("%.0f-%.0f cm",xmin,xmax), "l");
        leg0_el->AddEntry(h1d_el2, Form("%.0f-%.0f cm",xmin2,xmax2), "l");
        leg0_el->Draw();

	c_y_el->Print(fig_out_el_zoomy.Data());





      	for (int binx=1; binx<=n_x; binx++) {
        	for (int biny=1;biny<=n_y;biny++) {
            		//if (GetBinContent(binx,biny) > threshold) return binx;
         	}
      	}
	
	h2d_inel
	for (int i=1; i <= hist->GetNbinsX(); i++) {
   		if (hist->GetBinContent(i)) {
      			minXbin = i;
      			break;
   		}
	}

       	int n_b=h2d->GetNbinsX();
        double b_min=h2d->GetXaxis()->GetXmin();
        double b_max=h2d->GetXaxis()->GetXmax();
	cout<<"n_b:"<<n_b<<endl;
	cout<<"b_min:"<<b_min<<" - b_max:"<<b_max<<endl;

        //int n_dr=600;
        //double dr_min=-150;
        //double dr_max=150;
        int n_dr=h2d->GetNbinsY();
        double dr_min=h2d->GetYaxis()->GetXmin();
        double dr_max=h2d->GetYaxis()->GetXmax();










	TH1D *zst_mc=(TH1D *)fmc1->Get(Form("reco_startZ_sce"));
	zst_mc->SetLineColor(4);

	TString file_mc2=Form("../mc_truelen_fixtruelencalc.root");
	TFile *fmc2 = TFile::Open(file_mc2.Data());
	TH2D *h2d_mc=(TH2D *)fmc2->Get(Form("h2d_trueEndZ_ketrue_NoCut_inel"));
	TH1D *zEnd_true=h2d_mc->ProjectionX();
	zEnd_true->SetLineColor(6); 	

	TString file_data=Form("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_run5387_prod4a_thinslice_dx4cm_25slcs.root");
	TFile *fdata = TFile::Open(file_data.Data());
	TH1D *zst_data=(TH1D *)fdata->Get(Form("reco_startZ_sce"));
	zst_data->SetLineColor(2);
	
	//config -------------------------------------------------------//
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleYOffset(1.2);
	//gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	gStyle->SetOptStat(0);
	gStyle->SetPadLeftMargin(0.15);
        gStyle->SetFrameLineWidth(1);

	//gStyle->SetPalette(60);
	//gStyle->SetPalette(kLightTemperature);
	//gStyle->SetPalette(kTemperatureMap);
	//gStyle->SetPalette(kColorPrintableOnGrey);
	//gStyle->SetPalette(kBlueGreenYellow);
	//gStyle->SetPalette(kFruitPunch);
	//gStyle->SetPalette(kGreenPink);
	//gStyle->SetPalette(kPastel);
	//gStyle->SetPalette(kRedBlue);
	//gStyle->SetPalette(kCool);
	//gStyle->SetPalette(kWaterMelon);
	//gStyle->SetPalette(kBlackBody);
	//gStyle->SetPadTopMargin(0.2);
	//gStyle->SetPadRightMargin(0.17);
	//config -------------------------------------------------------//

	TLine *line=new TLine(0,1,120,1);
	line->SetLineStyle(2);

	//inel
	TCanvas *c_ = new TCanvas("c_", "c_", 1200, 900);
	c_->Divide(1,2);
	c_->cd(1);
	zst_mc->Scale((double)zst_data->Integral()/(double)zst_mc->Integral());
	zst_mc->SetTitle("CaloSz Cut; Reco. StartZ [cm]; Counts");
	zst_mc->GetXaxis()->SetRangeUser(-4,9);
	zst_mc->Draw("hist");
	zst_data->Draw("hist same");

        TLegend *leg0 = new TLegend(0.5,0.65,.85,0.88);
        leg0->SetFillStyle(0);
        leg0->AddEntry(zst_data, "Data (after SCE corr.)", "l");
        leg0->AddEntry(zst_mc, "MC (after SCE corr.)", "l");
        leg0->Draw();

	c_->cd(2);
	zEnd_true->SetTitle("Inelastic Scattering Protons; True StartZ [cm]; Counts");
	zEnd_true->GetXaxis()->SetRangeUser(-30,50);
	zEnd_true->Draw("hist");


	TString fig_out=Form("startZ_afterSCE.eps");
	c_->Print(fig_out.Data());




*/


}
