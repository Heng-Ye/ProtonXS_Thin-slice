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
	TH1D *h1d_data=(TH1D *)f_data->Get(Form("h1d_%s_%s",rep.Data(), str_cut.Data()));
	int n_data=h1d_data->Integral();

	//TH1D *trklen_data=(TH1D *)f_data->Get(Form("trklen_reco_%s",str_cut.Data()));
	//trklen_data->SetMarkerColor(1);
	//trklen_data->SetLineColor(1);
	
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
	TH1D *h1d_inel=(TH1D *)f_mc->Get(Form("h1d_%s_%s_inel",rep.Data(), str_cut.Data()));
	TH1D *h1d_el=(TH1D *)f_mc->Get(Form("h1d_%s_%s_el", rep.Data(), str_cut.Data()));

	TH1D *h1d_midcosmic=(TH1D *)f_mc->Get(Form("h1d_%s_%s_midcosmic", rep.Data(), str_cut.Data()));
	TH1D *h1d_midpi=(TH1D *)f_mc->Get(Form("h1d_%s_%s_midpi", rep.Data(), str_cut.Data()));
	TH1D *h1d_midp=(TH1D *)f_mc->Get(Form("h1d_%s_%s_midp", rep.Data(), str_cut.Data()));

	TH1D *h1d_midmu=(TH1D *)f_mc->Get(Form("h1d_%s_%s_midmu", rep.Data(), str_cut.Data()));
	TH1D *h1d_mideg=(TH1D *)f_mc->Get(Form("h1d_%s_%s_mideg",rep.Data(), str_cut.Data()));
	TH1D *h1d_midother=(TH1D *)f_mc->Get(Form("h1d_%s_%s_midother", rep.Data(), str_cut.Data()));

	TH1D *h1d_recoinel=(TH1D *)f_mc->Get(Form("h1d_%s_RecoInel", rep.Data()));

	h1d_inel->SetFillColor(2); h1d_inel->SetLineColor(2);
	h1d_el->SetFillColor(4); h1d_el->SetLineColor(4);

	h1d_midp->SetFillColor(3); h1d_midp->SetLineColor(3);
	h1d_midcosmic->SetFillColor(5); h1d_midcosmic->SetLineColor(5);
	h1d_midpi->SetFillColor(6); h1d_midpi->SetLineColor(6);

	h1d_midmu->SetFillColor(28); h1d_midmu->SetLineColor(28);
	h1d_mideg->SetFillColor(30); h1d_mideg->SetLineColor(30);
	h1d_midother->SetFillColor(15); h1d_midother->SetLineColor(15);

	h1d_recoinel->SetFillColor(0); h1d_recoinel->SetLineColor(17);
	h1d_recoinel->SetLineWidth(2); h1d_recoinel->SetLineStyle(7);

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
	h1d_recoinel->Scale(norm_mc);





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
/*
	//[1]
	TCanvas *c2d_0=new TCanvas(Form("c2d0"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c2d_0->Divide(1,1);
	c2d_0->cd(1);
	h2d->Draw("");
	c2d_0->Print(Form("%s/h1d_chi2pid_%s.eps",fout_path.Data(),str_cut.Data()));
	//h2d->Draw("colz");
	//c2d_0->Print(Form("%s/h1d_chi2pid_%s_colz.eps",fout_path.Data(),str_cut.Data()));

	//[2]
	TCanvas *c2d_1=new TCanvas(Form("c2d1"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c2d_1->Divide(1,1);
	c2d_1->cd(1);
	h2d_inel->Draw("");
	c2d_1->Print(Form("%s/h1d_chi2pid_%s_inel.eps",fout_path.Data(),str_cut.Data()));
	//h2d_inel->Draw("colz");
	//c2d_1->Print(Form("%s/h1d_chi2pid_%s_inel_colz.eps",fout_path.Data(),str_cut.Data()));

	//[2]
	TCanvas *c2d_2=new TCanvas(Form("c2d2"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c2d_2->Divide(1,1);
	c2d_2->cd(1);
	h2d_el->Draw("");
	c2d_2->Print(Form("%s/h1d_chi2pid_%s_el.eps",fout_path.Data(),str_cut.Data()));
	//h2d_el->Draw("colz");
	//c2d_2->Print(Form("%s/h1d_chi2pid_%s_el_colz.eps",fout_path.Data(),str_cut.Data()));

	//[3]
	TCanvas *c2d_3=new TCanvas(Form("c2d3"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c2d_3->Divide(1,1);
	c2d_3->cd(1);
	h2d_midp->Draw("");
	c2d_3->Print(Form("%s/h1d_chi2pid_%s_midp.eps",fout_path.Data(),str_cut.Data()));
	//h2d_midp->Draw("colz");
	//c2d_3->Print(Form("%s/h1d_chi2pid_%s_midp_colz.eps",fout_path.Data(),str_cut.Data()));

	//[2]&[3]
	TCanvas *c2d_4=new TCanvas(Form("c2d4"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c2d_4->Divide(1,1);
	c2d_4->cd(1);
	h2d_inel->Draw("");
	h2d_midp->SetMarkerSize(0.5);	
	h2d_midp->Draw("same");
	c2d_4->Print(Form("%s/h1d_chi2pid_%s_inel_midp.eps",fout_path.Data(),str_cut.Data()));
*/

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
	//c1d_0->Print(Form("%s/%s_%s.eps", fout_path.Data(), rep.Data(), str_cut.Data()));
	//c1d_0->Print(Form("%s/%s_%s_zoom_logy.eps", fout_path.Data(), rep.Data(), str_cut.Data()));
	//c1d_0->Print(Form("%s/%s_%s_zoom2_logy.eps", fout_path.Data(), rep.Data(), str_cut.Data()));
	c1d_0->Print(Form("%s/%s_%s_all_logy.eps", fout_path.Data(), rep.Data(), str_cut.Data()));
	



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
