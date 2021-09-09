
void plot_bkgfit_recosliceid(TString fin_data, TString fin, TString fout_path, TString fout_filename, TString rep, TString title1) {

	//TString rep="h_recosliceid_allevts_cuts";
	TString title=Form("%s; Reco SliceID; Events", title1.Data());
	//TString title=Form("%s; Reco SliceID; Events", "All Protons");

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	//gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	gStyle->SetOptStat(0);
	gStyle->SetPadLeftMargin(0.13);

	//get data ----------------------------------------------------//
	TFile *f_data = TFile::Open(fin_data.Data());
	TH1D *h1d_data=(TH1D *)fin_data->Get(Form("%s",rep.Data()));
	int n_data=h1d_data->Integral();
	cout<<"n_data:"<<n_data<<endl;

	//get mc -----------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fin.Data());
	TH1D *h1d_inel=(TH1D *)f_mc->Get(Form("%s_inel",rep.Data()));
	TH1D *h1d_el=(TH1D *)f_mc->Get(Form("%s_el", rep.Data()));
	TH1D *h1d_midcosmic=(TH1D *)f_mc->Get(Form("%s_midcosmic", rep.Data()));
	TH1D *h1d_midpi=(TH1D *)f_mc->Get(Form("%s_midpi", rep.Data()));
	TH1D *h1d_midp=(TH1D *)f_mc->Get(Form("%s_midp", rep.Data()));
	TH1D *h1d_midmu=(TH1D *)f_mc->Get(Form("%s_midmu", rep.Data()));
	TH1D *h1d_mideg=(TH1D *)f_mc->Get(Form("%s_mideg",rep.Data()));
	TH1D *h1d_midother=(TH1D *)f_mc->Get(Form("%s_midother", rep.Data()));

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

	//h1d_mc->GetXaxis()->SetTitle("#chi^{2}PID");
	
	cout<<"n_mc:"<<n_mc<<endl;
	//cout<<"h1d_mc->Integral():"<<h1d_mc->Integral()<<endl;

	THStack* hs=new THStack("hs","");
	hs->Add(h1d_inel);
	hs->Add(h1d_el);
	hs->Add(h1d_midp);
	hs->Add(h1d_midcosmic);
	hs->Add(h1d_midpi);
	hs->Add(h1d_midmu);
	hs->Add(h1d_mideg);
	hs->Add(h1d_midother);

	TCanvas *c_ = new TCanvas("c_", "c_", 900,600);
	float xmin=-2;
	float xmax=26;
        float ymax=4100;

	TH2D *f2d=new TH2D("f2d",Form("%s",""),100,xmin,xmax,3000,0,ymax); //logy
	f2d->SetTitle(title.Data());
	f2d->GetYaxis()->SetTitleOffset(1.3);
	f2d->Draw();
	hs->Draw("same hist");

	TLegend *leg = new TLegend(0.2,0.5,0.8,0.9);
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

	c_->Print(Form("%s/%s.eps",fout_path.Data(), fout_filename.Data()));



}
