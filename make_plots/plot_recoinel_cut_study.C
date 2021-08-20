#include "THStack.h"

void plot_recoinel_cut_study(TString fin_mc, TString fin_data, TString fout_path, TString rep, TString str_cut) {

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	//gStyle->SetOptStat(0);
	
	//data --------------------------------------------------------------------------------------//
	TFile *f_data = TFile::Open(fin_data.Data());
	//TH1D *trklen_data=(TH1D *)f_data->Get(Form("trklen_reco_%s",str_cut.Data()));
	//trklen_data->SetMarkerColor(1);
	//trklen_data->SetLineColor(1);
	//int n_data=trklen_data->Integral();
	
	//mc ----------------------------------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fin_mc.Data());
	//2d
/*
	TH2D *h2d = (TH2D*)f_mc->Get(Form("ntrklen_chi2pid_%s",str_cut.Data()));
	TH2D *h2d_inel = (TH2D*)f_mc->Get(Form("ntrklen_chi2pid_%s_inel",str_cut.Data()));
	TH2D *h2d_el = (TH2D*)f_mc->Get(Form("ntrklen_chi2pid_%s_el",str_cut.Data()));
	TH2D *h2d_midp = (TH2D*)f_mc->Get(Form("ntrklen_chi2pid_midp"));
	
	h2d->SetTitle(";Proton Track Length/CSDA [a.u.]; #chi^{2} PID [a.u.]");
	h2d_inel->SetTitle("True Inel.;Proton Track Length/CSDA [a.u.]; #chi^{2} PID [a.u.]");
	h2d_el->SetTitle("True El.;Proton Track Length/CSDA [a.u.]; #chi^{2} PID [a.u.]");
	h2d_midp->SetTitle("MisID:P;Proton Track Length/CSDA [a.u.]; #chi^{2} PID [a.u.]");

	h2d->SetMarkerSize(0.1);

	h2d_inel->SetMarkerSize(0.2);
	h2d_inel->SetMarkerColor(2);

	h2d_el->SetMarkerSize(0.2);
	h2d_el->SetMarkerColor(4);

	h2d_midp->SetMarkerSize(0.2);
	h2d_midp->SetMarkerColor(3);
*/
	
	//1d: ntrklen
	TH1D *ntrklen_inel=(TH1D *)f_mc->Get(Form("h1d_%s_%s_inel",rep.Data(), str_cut.Data()));
	TH1D *ntrklen_el=(TH1D *)f_mc->Get(Form("h1d_%s_%s_el", rep.Data(), str_cut.Data()));

	TH1D *ntrklen_midcosmic=(TH1D *)f_mc->Get(Form("h1d_%s_%s_midcosmic", rep.Data(), str_cut.Data()));
	TH1D *ntrklen_midpi=(TH1D *)f_mc->Get(Form("h1d_%s_%s_midpi", rep.Data(), str_cut.Data()));
	TH1D *ntrklen_midp=(TH1D *)f_mc->Get(Form("h1d_%s_%s_midp", rep.Data(), str_cut.Data()));

	TH1D *ntrklen_midmu=(TH1D *)f_mc->Get(Form("h1d_%s_%s_midmu", rep.Data(), str_cut.Data()));
	TH1D *ntrklen_mideg=(TH1D *)f_mc->Get(Form("h1d_%s_%s_mideg",rep.Data(), str_cut.Data()));
	TH1D *ntrklen_midother=(TH1D *)f_mc->Get(Form("h1d_%s_%s_midother", rep.Data(), str_cut.Data()));

	ntrklen_inel->SetFillColor(2); ntrklen_inel->SetLineColor(2);
	ntrklen_el->SetFillColor(4); ntrklen_el->SetLineColor(4);

	ntrklen_midp->SetFillColor(3); ntrklen_midp->SetLineColor(3);
	ntrklen_midcosmic->SetFillColor(5); ntrklen_midcosmic->SetLineColor(5);
	ntrklen_midpi->SetFillColor(6); ntrklen_midpi->SetLineColor(6);

	ntrklen_midmu->SetFillColor(28); ntrklen_midmu->SetLineColor(28);
	ntrklen_mideg->SetFillColor(30); ntrklen_mideg->SetLineColor(30);
	ntrklen_midother->SetFillColor(15); ntrklen_midother->SetLineColor(15);

	//ntrklen_inel->SetTitle("; Proton Track Length/CSDA;");
	//ntrklen_inel->SetTitle("; #chi^{2} PID;");
	ntrklen_inel->SetTitle("; Proton Track Length [cm];");
	THStack* hs_ntrklen=new THStack("hs_ntrklen","");
	hs_ntrklen->Add(ntrklen_inel);
	hs_ntrklen->Add(ntrklen_el);
	hs_ntrklen->Add(ntrklen_midp);
	hs_ntrklen->Add(ntrklen_midcosmic);
	hs_ntrklen->Add(ntrklen_midpi);
	hs_ntrklen->Add(ntrklen_midmu);
	hs_ntrklen->Add(ntrklen_mideg);
	hs_ntrklen->Add(ntrklen_midother);

	//2D dists -------------------------------------------------------------------------------//
/*
	//[1]
	TCanvas *c2d_0=new TCanvas(Form("c2d0"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c2d_0->Divide(1,1);
	c2d_0->cd(1);
	h2d->Draw("");
	c2d_0->Print(Form("%s/ntrklen_chi2pid_%s.eps",fout_path.Data(),str_cut.Data()));
	//h2d->Draw("colz");
	//c2d_0->Print(Form("%s/ntrklen_chi2pid_%s_colz.eps",fout_path.Data(),str_cut.Data()));

	//[2]
	TCanvas *c2d_1=new TCanvas(Form("c2d1"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c2d_1->Divide(1,1);
	c2d_1->cd(1);
	h2d_inel->Draw("");
	c2d_1->Print(Form("%s/ntrklen_chi2pid_%s_inel.eps",fout_path.Data(),str_cut.Data()));
	//h2d_inel->Draw("colz");
	//c2d_1->Print(Form("%s/ntrklen_chi2pid_%s_inel_colz.eps",fout_path.Data(),str_cut.Data()));

	//[2]
	TCanvas *c2d_2=new TCanvas(Form("c2d2"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c2d_2->Divide(1,1);
	c2d_2->cd(1);
	h2d_el->Draw("");
	c2d_2->Print(Form("%s/ntrklen_chi2pid_%s_el.eps",fout_path.Data(),str_cut.Data()));
	//h2d_el->Draw("colz");
	//c2d_2->Print(Form("%s/ntrklen_chi2pid_%s_el_colz.eps",fout_path.Data(),str_cut.Data()));

	//[3]
	TCanvas *c2d_3=new TCanvas(Form("c2d3"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c2d_3->Divide(1,1);
	c2d_3->cd(1);
	h2d_midp->Draw("");
	c2d_3->Print(Form("%s/ntrklen_chi2pid_%s_midp.eps",fout_path.Data(),str_cut.Data()));
	//h2d_midp->Draw("colz");
	//c2d_3->Print(Form("%s/ntrklen_chi2pid_%s_midp_colz.eps",fout_path.Data(),str_cut.Data()));

	//[2]&[3]
	TCanvas *c2d_4=new TCanvas(Form("c2d4"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c2d_4->Divide(1,1);
	c2d_4->cd(1);
	h2d_inel->Draw("");
	h2d_midp->SetMarkerSize(0.5);	
	h2d_midp->Draw("same");
	c2d_4->Print(Form("%s/ntrklen_chi2pid_%s_inel_midp.eps",fout_path.Data(),str_cut.Data()));
*/

	//[1d distributions]
	//ntrklen
	TCanvas *c1d_0=new TCanvas(Form("c1d0"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	c1d_0->Divide(1,1);
	c1d_0->cd(1);
	//c1d_0->cd(1)->SetLogy();
	//c1d_0->cd(1)->SetGrid();

	//TH2D *f2d_ntrklen=new TH2D("f2d_ntrklen","",100,-0.02,1.2,100,0,5600);
	TH2D *f2d_ntrklen=new TH2D("f2d_ntrklen","",150, 0,150,1000,0,6600);
	f2d_ntrklen->GetXaxis()->SetTitle("Proton Track Length [cm]");
	//f2d_ntrklen->GetXaxis()->SetTitle("#chi^{2} PID");
   	f2d_ntrklen->GetXaxis()->SetNdivisions(-502);
	f2d_ntrklen->Draw();
	hs_ntrklen->Draw("same");
		
	TLegend *leg = new TLegend(0.2,0.6,0.8,0.9);
	leg->SetFillStyle(0);
	//leg->AddEntry(trklen_data, "Data", "ep");
	leg->AddEntry(ntrklen_inel, "Inel","f");
	leg->AddEntry(ntrklen_el, "El","f");
	leg->AddEntry(ntrklen_midcosmic,"misID:cosmic","f");
	leg->AddEntry(ntrklen_midp, "misID:p","f");
	leg->AddEntry(ntrklen_midpi, "misID:#pi","f");
	leg->AddEntry(ntrklen_midmu,"misID:#mu","f");
	leg->AddEntry(ntrklen_mideg, "misID:e/#gamma","f");
	leg->AddEntry(ntrklen_midother, "misID:other","f");
	leg->SetNColumns(3);
	leg->Draw();
	c1d_0->Print(Form("%s/%s_%s.eps", fout_path.Data(), rep.Data(), str_cut.Data()));
	



/*
	TH1D *trklen_inel=(TH1D *)f0->Get(Form("trklen_reco_inel_%s",str_cut.Data()));
	TH1D *trklen_el=(TH1D *)f0->Get(Form("trklen_reco_el_%s",str_cut.Data()));

	TH1D *trklen_midcosmic=(TH1D *)f0->Get(Form("trklen_reco_midcosmic_%s",str_cut.Data()));
	TH1D *trklen_midpi=(TH1D *)f0->Get(Form("trklen_reco_midpi_%s",str_cut.Data()));
	TH1D *trklen_midp=(TH1D *)f0->Get(Form("trklen_reco_midp_%s",str_cut.Data()));

	TH1D *trklen_midmu=(TH1D *)f0->Get(Form("trklen_reco_midmu_%s",str_cut.Data()));
	TH1D *trklen_mideg=(TH1D *)f0->Get(Form("trklen_reco_mideg_%s",str_cut.Data()));
	TH1D *trklen_midother=(TH1D *)f0->Get(Form("trklen_reco_midother_%s",str_cut.Data()));

	trklen_inel->SetFillColor(2); trklen_inel->SetLineColor(2);
	trklen_el->SetFillColor(4); trklen_el->SetLineColor(4);

	trklen_midp->SetFillColor(3); trklen_midp->SetLineColor(3);
	trklen_midcosmic->SetFillColor(5); trklen_midcosmic->SetLineColor(5);
	trklen_midpi->SetFillColor(6); trklen_midpi->SetLineColor(6);

	trklen_midmu->SetFillColor(28); trklen_midmu->SetLineColor(28);
	trklen_mideg->SetFillColor(30); trklen_mideg->SetLineColor(30);
	trklen_midother->SetFillColor(15); trklen_midother->SetLineColor(15);

	int n_mc_inel=trklen_inel->Integral();
	int n_mc_el=trklen_el->Integral();
	int n_mc_midcosmic=trklen_midcosmic->Integral();
	int n_mc_midpi=trklen_midpi->Integral();
	int n_mc_midp=trklen_midpi->Integral();
	int n_mc_midmu=trklen_midmu->Integral();
	int n_mc_mideg=trklen_mideg->Integral();
	int n_mc_midother=trklen_midother->Integral();
	int n_mc=n_mc_inel+n_mc_el+n_mc_midcosmic+n_mc_midpi+n_mc_midp+n_mc_midmu+n_mc_mideg+n_mc_midother;

	double norm_mc=(double)n_data/(double)n_mc;
	trklen_inel->Scale(norm_mc);
	trklen_el->Scale(norm_mc);
	trklen_midcosmic->Scale(norm_mc);
	trklen_midpi->Scale(norm_mc);
	trklen_midp->Scale(norm_mc);
	trklen_midmu->Scale(norm_mc);
	trklen_mideg->Scale(norm_mc);
	trklen_midother->Scale(norm_mc);

	//--------------------------------------------------------------------------------------------//

	trklen_inel->GetXaxis()->SetTitle("Reco Track Length [cm]");

	THStack* hs=new THStack("hs","");
	hs->Add(trklen_inel);
	hs->Add(trklen_el);

	hs->Add(trklen_midcosmic);
	hs->Add(trklen_midp);
	hs->Add(trklen_midpi);

	hs->Add(trklen_midmu);
	hs->Add(trklen_mideg);
	hs->Add(trklen_midother);


	TCanvas *c_=new TCanvas(Form("c"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_->Divide(1,1);
	c_->cd(1);
	TH2D *f2d=new TH2D("f2d",Form("%s",str_cut.Data()),136,-4,132,1100,0,1100);
	f2d->GetXaxis()->SetTitle("Reco Track Length [cm]");
	f2d->Draw();
	//c_->cd(1)->SetLogy();
	hs->Draw("hist same");
	trklen_data->Draw("ep same");

	TLegend *leg = new TLegend(0.2,0.6,0.8,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(trklen_data, "Data", "ep");
	leg->AddEntry(trklen_inel, "Inel","f");
	leg->AddEntry(trklen_el, "El","f");

	leg->AddEntry(trklen_midcosmic,"misID:cosmic","f");
	leg->AddEntry(trklen_midp, "misID:p","f");
	leg->AddEntry(trklen_midpi, "misID:#pi","f");

	leg->AddEntry(trklen_midmu,"misID:#mu","f");
	leg->AddEntry(trklen_mideg, "misID:e/#gamma","f");
	leg->AddEntry(trklen_midother, "misID:other","f");

	leg->SetNColumns(3);
	leg->Draw();

	c_->Print(Form("%s/trklenreco_%s.eps",fout_path.Data(),str_cut.Data()));

*/


}
