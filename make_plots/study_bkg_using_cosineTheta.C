#include "THStack.h"


double err(double n1, double er_n1, double n2, double er_n2) {
	double r=n1/n2;
	double er_r=sqrt(pow((er_n1/n1),2)+pow((er_n2/n2),2));
	er_r*=r;

	return er_r;
}


void study_bkg_using_cosineTheta(TString fin, TString fin_data, TString rep, TString fout_path) {

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	//gStyle->SetOptStat(0);

	//data ---------------------------------------------------------//
	TFile *f_data = TFile::Open(fin_data.Data());
	TH1D *cosine_data=(TH1D *)f_data->Get(Form("%s",rep.Data()));
	int n_data=cosine_data->Integral();

	//mc -----------------------------------------------------------------------------//
	TFile *f0 = TFile::Open(fin.Data());

	TH1D *cosine_mc=(TH1D *)f0->Get(Form("%s",rep.Data()));
	TH1D *cosine_inel=(TH1D *)f0->Get(Form("%s_inel",rep.Data()));
	TH1D *cosine_el=(TH1D *)f0->Get(Form("%s_el",rep.Data()));

	TH1D *cosine_midcosmic=(TH1D *)f0->Get(Form("%s_midcosmic",rep.Data()));
	TH1D *cosine_midpi=(TH1D *)f0->Get(Form("%s_midpi",rep.Data()));
	TH1D *cosine_midp=(TH1D *)f0->Get(Form("%s_midp",rep.Data()));

	TH1D *cosine_midmu=(TH1D *)f0->Get(Form("%s_midmu",rep.Data()));
	TH1D *cosine_mideg=(TH1D *)f0->Get(Form("%s_mideg",rep.Data()));
	TH1D *cosine_midother=(TH1D *)f0->Get(Form("%s_midother",rep.Data()));

	cosine_inel->SetFillColor(2); cosine_inel->SetLineColor(2);
	cosine_el->SetFillColor(4); cosine_el->SetLineColor(4);

	cosine_midp->SetFillColor(3); cosine_midp->SetLineColor(3);
	cosine_midcosmic->SetFillColor(5); cosine_midcosmic->SetLineColor(5);
	cosine_midpi->SetFillColor(6); cosine_midpi->SetLineColor(6);

	cosine_midmu->SetFillColor(28); cosine_midmu->SetLineColor(28);
	cosine_mideg->SetFillColor(30); cosine_mideg->SetLineColor(30);
	cosine_midother->SetFillColor(15); cosine_midother->SetLineColor(15);

	int n_mc_inel=cosine_inel->Integral();
	int n_mc_el=cosine_el->Integral();
	int n_midcosmic=cosine_midcosmic->Integral();
	int n_midpi=cosine_midpi->Integral();
	int n_midp=cosine_midp->Integral();
	int n_midmu=cosine_midmu->Integral();
	int n_mideg=cosine_mideg->Integral();
	int n_midother=cosine_midother->Integral();
	int n_mc=n_mc_inel+n_mc_el+n_midcosmic+n_midpi+n_midp+n_midmu+n_mideg+n_midother;
	double mc_scale=(double)n_data/(double)n_mc;
	cosine_inel->Scale(mc_scale);
	cosine_el->Scale(mc_scale);
	cosine_midcosmic->Scale(mc_scale);
	cosine_midpi->Scale(mc_scale);
	cosine_midp->Scale(mc_scale);
	cosine_midmu->Scale(mc_scale);
	cosine_mideg->Scale(mc_scale);
	cosine_midother->Scale(mc_scale);
	cosine_mc->Scale(mc_scale);
	cosine_el->GetXaxis()->SetTitle("cos#Theta");

	//data/MC -----------------------------------------//
	TH1D *R=(TH1D*)cosine_data->Clone();
	R->Divide(cosine_data,cosine_mc);
	//-------------------------------------------------//


	THStack* hs=new THStack("hs","");
	hs->Add(cosine_inel);
	hs->Add(cosine_el);

	hs->Add(cosine_midcosmic);
	hs->Add(cosine_midp);
	hs->Add(cosine_midpi);

	hs->Add(cosine_midmu);
	hs->Add(cosine_mideg);
	hs->Add(cosine_midother);




	TCanvas *c_ = new TCanvas("c_", "c_", 800,800);
	TPad *pad1 = new TPad("pad1","pad1",0,0.3,1.,1.);
	pad1->SetBottomMargin(0.03);
	pad1->SetGridx();
	pad1->SetGridy();
	pad1->Draw();
	pad1->cd();
	pad1->SetLogy();

	float xmin=0.4;
	float xmax=1.;
	float ymax=92000;
	TH2D *f2d=new TH2D("f2d",Form("%s",""),100,xmin,xmax,92000,0,ymax); //logy

	f2d->Draw();
	f2d->GetXaxis()->SetLabelSize(0);
	hs->Draw("hist same");
	//cosine_mc->Draw("hist same");
	cosine_data->Draw("ep same");

	TLegend *leg = new TLegend(0.2,0.6,0.8,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(cosine_data, "Data","ep");
	leg->AddEntry(cosine_inel, "Inel","f");
	leg->AddEntry(cosine_el, "El","f");

	leg->AddEntry(cosine_midcosmic,"misID:cosmic","f");
	leg->AddEntry(cosine_midp, "misID:p","f");
	leg->AddEntry(cosine_midpi, "misID:#pi","f");

	leg->AddEntry(cosine_midmu,"misID:#mu","f");
	leg->AddEntry(cosine_mideg, "misID:e/#gamma","f");
	leg->AddEntry(cosine_midother, "misID:other","f");

	leg->SetNColumns(3);
	leg->Draw();

        //pDUNE Logo
        float logo_y=ymax+9000;
        TLatex **txt_pdune1=new TLatex*[1];
        txt_pdune1[0]=new TLatex(xmin, logo_y, Form("#bf{DUNE:ProtoDUNE-SP}"));
        txt_pdune1[0]->SetTextColor(1);
        txt_pdune1[0]->Draw();
        //
        //Beam Logo
        TLatex **txt_p1=new TLatex*[1];
        txt_p1[0]=new TLatex(xmax-.185, logo_y, Form("Protons (1 GeV/c)"));
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

	TH2D *f2d2=new TH2D("f2d2",Form("%s","Pos"),100,xmin,xmax,64,0,6.4); //liny
	f2d2->SetTitle(";cos#Theta;Data/MC");
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



	//c_->Print(Form("%s/cosinereco_%s.eps",fout_path.Data(),"calosz"));
	//c_->Print(Form("%s/cosinereco_%s_allrange.eps",fout_path.Data(),"calosz"));
	//c_->Print(Form("%s/cosinereco_%s_allrange_liny.eps",fout_path.Data(),"calosz"));
	c_->Print(Form("%s/bkg2_study_%s_allrange_logy.eps",fout_path.Data(),rep.Data()));







	/*

	   TCanvas *c_=new TCanvas(Form("c"),"",900, 600);
	//TCanvas *c_=new TCanvas(Form("c"),"",1600, 10000);
	//	c_->SetWindowSize(900, 100);
	//c_-> Divide(9,57);
	//


	//allC = new TCanvas(“allC”,“All Plots”);
	//allC -> SetCanvasSize(1600, 10000);
	//allC -> SetWindowSize(1600, 800);
	//allC -> Divide(9,57);


	c_->Divide(1,2,0,1);
	c_->cd(1);
	c_->cd(1)->SetLogy();
	//c_->cd(1)->SetLogy(0);
	c_->cd(1)->SetGridx();
	c_->cd(1)->SetGridy();
	c_->cd(1)->SetBottomMargin(1);
	//gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	//TH2D *f2d=new TH2D("f2d",Form("%s","CaloSize"),100,0.9,1,92000,0,92000);
	//TH2D *f2d=new TH2D("f2d",Form("%s","Pos"),100,0.4,1.,92000,0,92000); //logy
	//TH2D *f2d=new TH2D("f2d",Form("%s","Pos"),100,0.4,1.,50,0,50); //liny
	float xmin=0.4;
	float xmax=1.;
	TH2D *f2d=new TH2D("f2d",Form("%s","Pos"),100,xmin,xmax,92000,0,92000); //logy
	//TH2D *f2d=new TH2D("f2d",Form("%s","CaloSize"),100,0.4,1.,50,0,50); //liny
	//TH2D *f2d=new TH2D("f2d",Form("%s","CaloSize"),100,0.,1.,2000,0,2000);
	//f2d->GetXaxis()->SetTitle("cos#Theta");

	//TPad *pad1 = new TPad(Form("pad1"), Form("pad1"), 0, 0., 1, 1.);
	//pad1->SetGridx();
	//pad1->SetGridy();
	//pad1->cd();
	//pad1->cd(1)->SetLogy();
	//pad1->Draw();

	//f2d->GetXaxis()->SetNdivisions(32);
	//f2d->Draw("same");


	f2d->Draw();
	hs->Draw("hist same");
	//cosine_mc->Draw("hist same");
	cosine_data->Draw("ep same");

	TLegend *leg = new TLegend(0.2,0.6,0.8,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(cosine_data, "Data","ep");
	leg->AddEntry(cosine_inel, "Inel","f");
	leg->AddEntry(cosine_el, "El","f");

	leg->AddEntry(cosine_midcosmic,"misID:cosmic","f");
	leg->AddEntry(cosine_midp, "misID:p","f");
	leg->AddEntry(cosine_midpi, "misID:#pi","f");

	leg->AddEntry(cosine_midmu,"misID:#mu","f");
	leg->AddEntry(cosine_mideg, "misID:e/#gamma","f");
	leg->AddEntry(cosine_midother, "misID:other","f");

	leg->SetNColumns(3);
	leg->Draw();

	c_->cd(2);
	TH2D *f2d2=new TH2D("f2d2",Form("%s","Pos"),100,xmin,xmax,50,0,5); //liny
	f2d2->SetTitle(";cos#Theta;Data/MC");
	f2d2->Draw();
	R->Draw("ep same");

	//c_->Print(Form("%s/cosinereco_%s.eps",fout_path.Data(),"calosz"));
	//c_->Print(Form("%s/cosinereco_%s_allrange.eps",fout_path.Data(),"calosz"));
	//c_->Print(Form("%s/cosinereco_%s_allrange_liny.eps",fout_path.Data(),"calosz"));
	c_->Print(Form("%s/bkg2_study_%s_allrange_logy.eps",fout_path.Data(),rep.Data()));
	//c_->Print(Form("%s/%s_allrange_liny.eps",fout_path.Data(),rep.Data()));
	*/


}
