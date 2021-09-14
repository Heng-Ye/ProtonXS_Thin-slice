#include "THStack.h"

void study_bkg_using_chi2pid(TString fin_mc, TString fin_data, TString rep, TString fout_path) {

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	//gStyle->SetOptStat(0);
	
	//data --------------------------------------------------------------------------------------//
	TFile *f_data = TFile::Open(fin_data.Data());
	TH1D *h1d_data=(TH1D *)f_data->Get(Form("h1d_%s",rep.Data()));
	int n_data=h1d_data->Integral();

	//mc ----------------------------------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fin_mc.Data());
	TH1D *h1d_inel=(TH1D *)f_mc->Get(Form("h1d_%s_inel",rep.Data()));
	TH1D *h1d_el=(TH1D *)f_mc->Get(Form("h1d_%s_el", rep.Data()));
	TH1D *h1d_midcosmic=(TH1D *)f_mc->Get(Form("h1d_%s_midcosmic", rep.Data()));
	TH1D *h1d_midpi=(TH1D *)f_mc->Get(Form("h1d_%s_midpi", rep.Data()));
	TH1D *h1d_midp=(TH1D *)f_mc->Get(Form("h1d_%s_midp", rep.Data()));
	TH1D *h1d_midmu=(TH1D *)f_mc->Get(Form("h1d_%s_midmu", rep.Data()));
	TH1D *h1d_mideg=(TH1D *)f_mc->Get(Form("h1d_%s_mideg",rep.Data()));
	TH1D *h1d_midother=(TH1D *)f_mc->Get(Form("h1d_%s_midother", rep.Data()));

	TH1D *h1d_mc=(TH1D*)h1d_inel->Clone(); h1d_mc->Sumw2();
	h1d_mc->Add(h1d_el);
	h1d_mc->Add(h1d_midcosmic);
	h1d_mc->Add(h1d_midpi);
	h1d_mc->Add(h1d_midp);
	h1d_mc->Add(h1d_midmu);
	h1d_mc->Add(h1d_mideg);
	h1d_mc->Add(h1d_midother);

	h1d_inel->SetFillColor(2); h1d_inel->SetLineColor(2);
	h1d_el->SetFillColor(4); h1d_el->SetLineColor(4);
	h1d_midp->SetFillColor(3); h1d_midp->SetLineColor(3);
	h1d_midcosmic->SetFillColor(5); h1d_midcosmic->SetLineColor(5);
	h1d_midpi->SetFillColor(6); h1d_midpi->SetLineColor(6);
	h1d_midmu->SetFillColor(28); h1d_midmu->SetLineColor(28);
	h1d_mideg->SetFillColor(30); h1d_mideg->SetLineColor(30);
	h1d_midother->SetFillColor(15); h1d_midother->SetLineColor(15);

	int n_mc_inel=h1d_inel->Integral();
	int n_mc_el=h1d_el->Integral();
	int n_mc_midcosmic=h1d_midcosmic->Integral();
	int n_mc_midpi=h1d_midpi->Integral();
	int n_mc_midp=h1d_midp->Integral();
	int n_mc_midmu=h1d_midmu->Integral();
	int n_mc_mideg=h1d_mideg->Integral();
	int n_mc_midother=h1d_midother->Integral();
	int n_mc=n_mc_inel+n_mc_el+n_mc_midcosmic+n_mc_midpi+n_mc_midp+n_mc_midmu+n_mc_mideg+n_mc_midother;
	double norm_mc=(double)n_data/(double)n_mc;

	h1d_inel->Scale(norm_mc);
	h1d_el->Scale(norm_mc);
	h1d_midcosmic->Scale(norm_mc);
	h1d_midpi->Scale(norm_mc);
	h1d_midp->Scale(norm_mc);
	h1d_midmu->Scale(norm_mc);
	h1d_mideg->Scale(norm_mc);
	h1d_midother->Scale(norm_mc);
	h1d_mc->Scale(norm_mc);
	h1d_mc->GetXaxis()->SetTitle("#chi^{2}PID");
	
	cout<<"n_mc:"<<n_mc<<endl;
	cout<<"h1d_mc->Integral():"<<h1d_mc->Integral()<<endl;

	//rebin histograms
	int n_rb=2;
	h1d_inel->Rebin(n_rb);
	h1d_el->Rebin(n_rb);
	h1d_midcosmic->Rebin(n_rb);
	h1d_midpi->Rebin(n_rb);
	h1d_midp->Rebin(n_rb);
	h1d_midmu->Rebin(n_rb);
	h1d_mideg->Rebin(n_rb);
	h1d_midother->Rebin(n_rb);
	h1d_mc->Rebin(n_rb);
	h1d_data->Rebin(n_rb);

	h1d_data->Scale(1./(double)n_rb);
	h1d_mc->Scale(1./(double)n_rb);
	h1d_midother->Scale(1./(double)n_rb);
	h1d_mideg->Scale(1./(double)n_rb);
	h1d_midmu->Scale(1./(double)n_rb);
	h1d_midp->Scale(1./(double)n_rb);
	h1d_midpi->Scale(1./(double)n_rb);
	h1d_midcosmic->Scale(1./(double)n_rb);
	h1d_el->Scale(1./(double)n_rb);
	h1d_inel->Scale(1./(double)n_rb);



	//data/MC -----------------------------------------//
	TH1D *R=(TH1D*)h1d_data->Clone();
	R->Divide(h1d_data, h1d_mc);
	//-------------------------------------------------//

	THStack* hs=new THStack("hs","");
	hs->Add(h1d_inel);
	hs->Add(h1d_el);
	hs->Add(h1d_midp);
	hs->Add(h1d_midcosmic);
	hs->Add(h1d_midpi);
	hs->Add(h1d_midmu);
	hs->Add(h1d_mideg);
	hs->Add(h1d_midother);

	TCanvas *c_ = new TCanvas("c_", "c_", 800,800);
	TPad *pad1 = new TPad("pad1","pad1",0,0.3,1.,1.);
	pad1->SetBottomMargin(0.03);
	pad1->SetGridx();
	pad1->SetGridy();
	pad1->Draw();
	pad1->cd();
	pad1->SetLogy();

	float xmin=0.;
	//float xmax=20;
	float xmax=40;
        float ymax=3000; //logy
        //float ymax=550; //logy

	TH2D *f2d=new TH2D("f2d",Form("%s",""),100,xmin,xmax,3000,0,ymax); //logy
	f2d->Draw();
	f2d->GetXaxis()->SetLabelSize(0);
	
	
	hs->Draw("hist same");
	//h1d_mc->Draw("hist same");
	h1d_data->Draw("ep same");

	TLegend *leg = new TLegend(0.2,0.6,0.8,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(h1d_data, "Data","ep");
	leg->AddEntry(h1d_inel, "Inel","f");
	leg->AddEntry(h1d_el, "El","f");

	leg->AddEntry(h1d_midcosmic,"misID:cosmic","f");
	leg->AddEntry(h1d_midp, "misID:p","f");
	leg->AddEntry(h1d_midpi, "misID:#pi","f");

	leg->AddEntry(h1d_midmu,"misID:#mu","f");
	leg->AddEntry(h1d_mideg, "misID:e/#gamma","f");
	leg->AddEntry(h1d_midother, "misID:other","f");

	leg->SetNColumns(3);
	leg->Draw();

        //pDUNE Logo
        float logo_y=ymax+150;
        TLatex **txt_pdune1=new TLatex*[1];
        txt_pdune1[0]=new TLatex(xmin, logo_y, Form("#bf{DUNE:ProtoDUNE-SP}"));
        txt_pdune1[0]->SetTextColor(1);
        txt_pdune1[0]->Draw();
        //
        //Beam Logo
        TLatex **txt_p1=new TLatex*[1];
        //txt_p1[0]=new TLatex(xmax-6.3, logo_y, Form("Protons (1 GeV/c)")); //x:0-20
        txt_p1[0]=new TLatex(xmax-12.3, logo_y, Form("Protons (1 GeV/c)")); //x:0-40
        txt_p1[0]->SetTextColor(1);
        txt_p1[0]->SetTextSize(0.05);
        txt_p1[0]->Draw();


	c_->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
	pad2->SetTopMargin(0.02);
	pad2->SetBottomMargin(0.25);
	pad2->SetGridx();
	pad2->SetGridy();
	pad2->Draw();
	pad2->cd();

	TH2D *f2d2=new TH2D("f2d2",Form("%s","BQ"),100,xmin,xmax,64,0,6.4); //liny
	f2d2->SetTitle(";#chi^{2}PID;Data/MC");
	f2d2->GetXaxis()->SetLabelSize(0.1);
	f2d2->GetYaxis()->SetLabelSize(0.1);

	f2d2->GetXaxis()->SetTitleSize(0.1);
	f2d2->GetYaxis()->SetTitleSize(0.1);
	f2d2->GetYaxis()->SetTitleOffset(0.3);


	f2d2->Draw();
	TLine *line1=new TLine(xmin,1,xmax,1);
	line1->SetLineColor(7);
	line1->SetLineWidth(2);
	line1->Draw("same");
	R->Draw("ep same");

	//c_->Print(Form("%s/bkg_study_%s_allrange_logy.eps",fout_path.Data(),rep.Data()));
	c_->Print(Form("%s/bkg_study_%s_0-40_rb_logy.eps",fout_path.Data(),rep.Data()));
	//c_->Print(Form("%s/bkg_study_%s_allrange_liny.eps",fout_path.Data(),rep.Data()));





/*
	//h1d_inel->SetTitle("; Proton Track Length/CSDA;");
	//h1d_inel->SetTitle("; #chi^{2} PID;");
	h1d_inel->SetTitle("; Proton Track Length [cm];");
	THStack* h1s=new THStack("h1s","");
	h1s->Add(h1d_inel);
	h1s->Add(h1d_el);
	h1s->Add(h1d_midp);
	h1s->Add(h1d_midcosmic);
	h1s->Add(h1d_midpi);
	h1s->Add(h1d_midmu);
	h1s->Add(h1d_mideg);
	h1s->Add(h1d_midother);

	//2D dists -------------------------------------------------------------------------------//
	//[1d distributions]
	//ntrklen
	//TCanvas *c1d_0=new TCanvas(Form("c1d0"),"",900, 600);
	TCanvas *c1d_0=new TCanvas(Form("c1d0"),"",1600, 900);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	c1d_0->Divide(1,1);
	c1d_0->cd(1);
	//c1d_0->cd(1)->SetLogy();
	c1d_0->cd(1)->SetGrid();

	//TH2D *f2d_ntrklen=new TH2D("f2d_ntrklen","",100,-0.02,1.2,100,0,5600);
	//TH2D *f2d_ntrklen=new TH2D("f2d_ntrklen","",150, 0,150,1000,0,6600);
	//TH2D *f2d_ntrklen=new TH2D("f2d_ntrklen","",150, 0,150,600,0,600);
	//TH2D *f2d_ntrklen=new TH2D("f2d_ntrklen","",20, 0,20,2500,0,2500);
	//TH2D *f2d_ntrklen=new TH2D("f2d_ntrklen","",55, 20,75,2500,0,2500);
	//TH2D *f2d_ntrklen=new TH2D("f2d_ntrklen","",75, 75,150,2500,0,2500);
	//TH2D *f2d_ntrklen=new TH2D("f2d_ntrklen","",40, 2,6,1300,0,1300); //median dedx
	TH2D *f2d_ntrklen=new TH2D("f2d_ntrklen","",40, 2,10,1300,0,1300);
	//f2d_ntrklen->GetXaxis()->SetTitle("Proton Track Length [cm]");
	//f2d_ntrklen->GetXaxis()->SetTitle("#chi^{2} PID");
	//f2d_ntrklen->GetXaxis()->SetTitle("Median dE/dx [MeV/cm]");
	f2d_ntrklen->GetXaxis()->SetTitle("Energy Deposition/Range [MeV/cm]");
   	f2d_ntrklen->GetXaxis()->SetNdivisions(-502);
	f2d_ntrklen->Draw();
	h1s->Draw("hist same");
	h1d_data->SetLineWidth(2);
	//h1d_data->Rebin(5);
	//h1d_data->Scale(1./5.);
	h1d_data->Draw("ep same");
	h1d_recoinel->Draw("hist same");
		
	TLegend *leg = new TLegend(0.3,0.6,0.85,0.9);
	leg->SetFillStyle(0);
	//leg->AddEntry(trklen_data, "Data", "ep");
	leg->AddEntry(h1d_data, "Data","ep");
	leg->AddEntry(h1d_inel, "Inel","f");
	leg->AddEntry(h1d_el, "El","f");
	leg->AddEntry(h1d_midcosmic,"misID:cosmic","f");
	leg->AddEntry(h1d_midp, "misID:p","f");
	leg->AddEntry(h1d_midpi, "misID:#pi","f");
	leg->AddEntry(h1d_midmu,"misID:#mu","f");
	leg->AddEntry(h1d_mideg, "misID:e/#gamma","f");
	leg->AddEntry(h1d_midother, "misID:other","f");
	leg->AddEntry(h1d_recoinel, "Reco Inel.(Default)","l");
	leg->SetNColumns(3);
	leg->Draw();
	//c1d_0->Print(Form("%s/%s_%s.eps", fout_path.Data(), rep.Data()));
	//c1d_0->Print(Form("%s/%s_%s_zoom_logy.eps", fout_path.Data(), rep.Data()));
	//c1d_0->Print(Form("%s/%s_%s_zoom2_logy.eps", fout_path.Data(), rep.Data()));
	c1d_0->Print(Form("%s/%s_%s_all_logy.eps", fout_path.Data(), rep.Data()));
	

*/



}
