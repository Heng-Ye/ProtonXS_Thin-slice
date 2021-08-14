void plot_ntrklen(TString fin, TString outpath){

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	//gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 

	gStyle->SetOptStat(0);

	TFile *f0 = TFile::Open(fin.Data());
	TH1D *ntrklen=(TH1D *)f0->Get("ntrklen");
	TH1D *ntrklen_PureInel=(TH1D *)f0->Get("ntrklen_PureInel");


	ntrklen->SetLineColor(1);  
	ntrklen_PureInel->SetLineColor(2);  
	ntrklen_PureInel->SetFillColor(2);  

	int n_rb=2;

	TH2D* frame2d=new TH2D("frame2d","", 120, 0, 1.2, 200, 0, 2000); //zend_2d
	frame2d->SetTitle("MC Protons;Track Length/CSDA [a.u.];");


	TCanvas *c0=new TCanvas("c0",""); 
	c0->Divide(1,1);
	c0->cd(1);
	frame2d->Draw();
	ntrklen_PureInel->Draw("same");
	ntrklen->Draw("same");

	TLegend *leg5 = new TLegend(0.15,0.65,0.75,0.85);
	leg5->SetFillStyle(0);
	leg5->AddEntry(ntrklen, "Total", "l");
	leg5->AddEntry(ntrklen_PureInel, "Truth Inel. Protons", "f");
	leg5->Draw();

	c0->Print(Form("%s/ntrklen.eps",outpath.Data()));


}

