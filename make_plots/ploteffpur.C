void ploteffpur(TString fin, TString fout_path){
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);

	gStyle->SetOptStat(0);

	TFile *file = TFile::Open(fin.Data());

	TH1D *eff_num_Int = (TH1D*)file->Get("eff_num_Int");
	TH1D *eff_den_Int = (TH1D*)file->Get("eff_den_Int");
	TH1D *eff_num_Inc = (TH1D*)file->Get("eff_num_Inc");
	TH1D *eff_den_Inc = (TH1D*)file->Get("eff_den_Inc");
	TH1D *pur_num_Int = (TH1D*)file->Get("pur_num_Int");
	TH1D *pur_num_Inc = (TH1D*)file->Get("pur_num_Inc");
	TH1D *pur_den = (TH1D*)file->Get("pur_den");
	TH2D *response_SliceID_Int = (TH2D*)file->Get("response_SliceID_Int");
	TH2D *response_SliceID_Inc = (TH2D*)file->Get("response_SliceID_Inc");

	TCanvas *c1 = new TCanvas("c1","c1");
	TEfficiency *eff_proton = 0;
	if (TEfficiency::CheckConsistency(*eff_num_Inc, *eff_den_Inc)){
		eff_proton = new TEfficiency(*eff_num_Inc, *eff_den_Inc);
	}
	eff_proton->SetTitle("All Protons;True slice ID;Efficiency");
	eff_proton->Draw();

	TCanvas *c2 = new TCanvas("c2","c2");
	TEfficiency *eff_protoninel = 0;
	if (TEfficiency::CheckConsistency(*eff_num_Int, *eff_den_Int)){
		eff_protoninel = new TEfficiency(*eff_num_Int, *eff_den_Int);
	}
	eff_protoninel->SetTitle("Proton Inelastic Scatterings;True slice ID;Efficiency");
	eff_protoninel->Draw();

	TCanvas *c3 = new TCanvas("c3","c3");
	TEfficiency *pur_proton = 0;
	if (TEfficiency::CheckConsistency(*pur_num_Inc, *pur_den)){
		pur_proton = new TEfficiency(*pur_num_Inc, *pur_den);
	}
	pur_proton->SetTitle("All Protons;Reco slice ID;Purity");
	pur_proton->Draw();

	TCanvas *c4 = new TCanvas("c4","c4");
	TEfficiency *pur_protoninel = 0;
	if (TEfficiency::CheckConsistency(*pur_num_Int, *pur_den)){
		pur_protoninel = new TEfficiency(*pur_num_Int, *pur_den);
	}
	pur_protoninel->SetTitle("Proton Inelastic Scatterings;Reco slice ID;Purity");
	pur_protoninel->Draw();

	TCanvas *c5 = new TCanvas("c5","c5");
	response_SliceID_Int->Draw("colz");

	TCanvas *c6 = new TCanvas("c6","c6");
	response_SliceID_Inc->Draw("colz");

	c1->Print(Form("%s/protoneff.eps",fout_path.Data()));
	c2->Print(Form("%s/protonineleff.eps",fout_path.Data()));
	c3->Print(Form("%s/protonpur.eps",fout_path.Data()));
	c4->Print(Form("%s/protoninelpur.eps",fout_path.Data()));
	c5->Print(Form("%s/protonres.eps",fout_path.Data()));
	c6->Print(Form("%s/protoninelres.eps",fout_path.Data()));

	c1->Print(Form("%s/protoneff.eps",fout_path.Data()));
	c2->Print(Form("%s/protonineleff.eps",fout_path.Data()));
	c3->Print(Form("%s/protonpur.eps",fout_path.Data()));
	c4->Print(Form("%s/protoninelpur.eps",fout_path.Data()));
	c5->Print(Form("%s/protonres.eps",fout_path.Data()));
	c6->Print(Form("%s/protoninelres.eps",fout_path.Data()));

}
