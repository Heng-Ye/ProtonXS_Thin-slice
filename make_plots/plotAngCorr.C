void plotAngCorr(TString fin, TString fout_path){

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	//gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 

	gStyle->SetOptStat(0);


	TFile *file = TFile::Open(fin.Data());

	TH1D *trueAngCorr = (TH1D*)file->Get("true_AngCorr");
	TH1D *recoAngCorr = (TH1D*)file->Get("reco_AngCorr");
	//TGraphErrors *gr_trueAngCorr = (TGraphErrors*)file->Get("true_AngCorr");
	//TGraphErrors *gr_recoAngCorr = (TGraphErrors*)file->Get("reco_AngCorr");
	//TGraphErrors *gr_reco_trueAngCorr = (TGraphErrors*)file->Get("gr_reco_trueAngCorr");

	TCanvas *c1 = new TCanvas("c1","c1");
	c1->cd(1)->SetLogy();	
	recoAngCorr->SetTitle("");
	//trueAngCorr->GetXaxis()->SetTitle("Slice ID");
	//trueAngCorr->GetXaxis()->SetRangeUser(-1, 23);
	//trueAngCorr->GetYaxis()->SetRangeUser(0.7, 1);
	//trueAngCorr->GetYaxis()->SetTitle("Angular correction");
	recoAngCorr->SetLineColor(2);
	recoAngCorr->SetMarkerColor(2);
	recoAngCorr->Draw("hist");
	trueAngCorr->Draw("hist same");
	TLegend *leg = new TLegend(0.2,0.3,0.6,0.5);
	leg->SetFillStyle(0);
	leg->AddEntry(trueAngCorr, "True angular correction", "pe");
	leg->AddEntry(recoAngCorr, "Reco angular correction", "pe");
	leg->Draw();

	/*
	   TCanvas *c2 = new TCanvas("c2","c2");
	   gr_reco_trueAngCorr->Draw("ape");
	   gr_reco_trueAngCorr->SetTitle("");
	   gr_reco_trueAngCorr->GetXaxis()->SetTitle("Slice ID");
	   gr_reco_trueAngCorr->GetXaxis()->SetRangeUser(-1, 23);
	   gr_reco_trueAngCorr->GetYaxis()->SetTitle("Reco - True angular correction (MeV)");
	   */

	c1->Print(Form("%s/AngCorr.eps",fout_path.Data()));
	//c2->Print(Form("%s/AngCorrdiff.eps",fout_path.Data()));

}
