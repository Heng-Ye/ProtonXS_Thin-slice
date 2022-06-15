
void plotKEff() {
	TString fdata="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_ke.root";
	TString fmc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/mc_ke.root";
	TString outpath="./plots_KE/";

	//TString str_obs="_recoinel";
	//TString str_cut="Reco InEl";
	//TString str_obs="";
	//TString str_cut="All";
	TString str_obs="_recoel";
	TString str_cut="Reco El";

	//plot style --------------------------------------------------------------------//
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	gStyle->SetOptStat(0);

	//read MC -----------------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fmc.Data());
	TH1D *keff=(TH1D*)f_mc->Get(Form("h1d_keff%s",str_obs.Data())); //keff at z=0
	TH1D *keff0=(TH1D*)f_mc->Get(Form("h1d_keff0%s",str_obs.Data())); //keff at tpc entrance
	TH1D *keff2=(TH1D*)f_mc->Get(Form("h1d_keff2%s",str_obs.Data())); //keff=ke_beam-const. E-loss

	TH1D *keff_inel=(TH1D*)f_mc->Get(Form("h1d_keff%s_inel",str_obs.Data())); //keff at z=0
	TH1D *keff_el=(TH1D*)f_mc->Get(Form("h1d_keff%s_el",str_obs.Data())); //keff at z=0
	TH1D *keff_midcosmic=(TH1D*)f_mc->Get(Form("h1d_keff%s_midcosmic",str_obs.Data())); //keff at z=0
	TH1D *keff_midpi=(TH1D*)f_mc->Get(Form("h1d_keff%s_midpi",str_obs.Data())); //keff at z=0
	TH1D *keff_midp=(TH1D*)f_mc->Get(Form("h1d_keff%s_midp",str_obs.Data())); //keff at z=0
	TH1D *keff_midmu=(TH1D*)f_mc->Get(Form("h1d_keff%s_midmu",str_obs.Data())); //keff at z=0
	TH1D *keff_mideg=(TH1D*)f_mc->Get(Form("h1d_keff%s_mideg",str_obs.Data())); //keff at z=0
	TH1D *keff_midother=(TH1D*)f_mc->Get(Form("h1d_keff%s_midother",str_obs.Data())); //keff at z=0

	
	keff_inel->SetFillColor(2); keff_inel->SetLineColor(2);
	keff_el->SetFillColor(7); keff_el->SetLineColor(7);
	keff_midcosmic->SetFillColor(5); keff_midcosmic->SetLineColor(5);
	keff_midpi->SetFillColor(6); keff_midpi->SetLineColor(6);
	keff_midp->SetFillColor(3); keff_midp->SetLineColor(3);
	keff_midmu->SetFillColor(28); keff_midmu->SetLineColor(28);
	keff_mideg->SetFillColor(30); keff_mideg->SetLineColor(30);
	keff_midother->SetFillColor(15); keff_midother->SetLineColor(15);
	
	keff->SetLineColor(16);
	keff->SetLineStyle(2);

	keff0->SetLineColor(1);

	keff2->SetLineColor(4);


	THStack* hs=new THStack("hs","");
	hs->Add(keff_inel);
	hs->Add(keff_el);

	hs->Add(keff_midcosmic);
	hs->Add(keff_midp);
	hs->Add(keff_midpi);

	hs->Add(keff_midmu);
	hs->Add(keff_mideg);
	hs->Add(keff_midother);



	TCanvas *c_=new TCanvas(Form("c"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_->Divide(1,1);
	c_->cd(1);
	TH2D *f2d=new TH2D("f2d",Form("%s",str_cut.Data()), 800,0.0,800,100,0.1,3000000);
	f2d->GetXaxis()->SetTitle("Proton Energy at TPC Front Face [MeV]");
	f2d->Draw();
	c_->cd(1)->SetLogy();
	hs->Draw("hist same");
	keff->Draw("same");
	keff0->Draw("same");
	keff2->Draw("same");



	TLegend *leg = new TLegend(0.2,0.6,0.8,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(keff, "KE_{ff} (z=0)","l");
   	leg->AddEntry((TObject*)0, "", "");	
   	leg->AddEntry((TObject*)0, "", "");	
	leg->AddEntry(keff0, "KE_{ff} (1st point)","l");	
   	leg->AddEntry((TObject*)0, "", "");	
   	leg->AddEntry((TObject*)0, "", "");	
	leg->AddEntry(keff2, "KE_{ff}=KE_{beam}-#DeltaE_{TPC}","l");	
   	leg->AddEntry((TObject*)0, "", "");	
   	leg->AddEntry((TObject*)0, "", "");	
	leg->AddEntry(keff_inel, "Inel","f");
	leg->AddEntry(keff_el, "El","f");

	leg->AddEntry(keff_midcosmic,"misID:cosmic","f");
	leg->AddEntry(keff_midp, "misID:p","f");
	leg->AddEntry(keff_midpi, "misID:#pi","f");

	leg->AddEntry(keff_midmu,"misID:#mu","f");
	leg->AddEntry(keff_mideg, "misID:e/#gamma","f");
	leg->AddEntry(keff_midother, "misID:other","f");

	leg->SetNColumns(3);
	leg->Draw();

	c_->Print(Form("%s/keff%s.eps",outpath.Data(),str_obs.Data()));



}
