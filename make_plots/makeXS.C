#include "/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4_MC_1GeV_reco1_sce_datadriven/SliceParams.h"


void makeXS(TString fin, TString outpath){

	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	gStyle->SetTitleX(0.5);
	//gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23); 
	gStyle->SetOptStat(0);
	gStyle->SetPadLeftMargin(0.13);


	TFile *file = TFile::Open(fin.Data());
	TH1D *eff_num_Int = (TH1D*)file->Get("eff_num_Int");
	TH1D *eff_den_Int = (TH1D*)file->Get("eff_den_Int");
	TH1D *eff_num_Inc = (TH1D*)file->Get("eff_num_Inc");
	TH1D *eff_den_Inc = (TH1D*)file->Get("eff_den_Inc");
	TH1D *pur_num_Int = (TH1D*)file->Get("pur_num_Int");
	TH1D *pur_num_Inc = (TH1D*)file->Get("pur_num_Inc");
	TH1D *pur_den = (TH1D*)file->Get("pur_den");
	TH1D *pur_den_Inc = (TH1D*)file->Get("pur_den_Inc");
	TH1D *pur_den_Int = (TH1D*)file->Get("pur_den_Int");

	TH1D *reco_AngCorr = (TH1D*)file->Get("reco_AngCorr");
	TH1D *true_AngCorr = (TH1D*)file->Get("true_AngCorr");

	TH1D *h_truesliceid_uf = (TH1D*)file->Get("h_truesliceid_uf");
	TH1D *h_truesliceid_inelastic_uf = (TH1D*)file->Get("h_truesliceid_inelastic_uf");

	TH1D *eff_Int = (TH1D*)eff_num_Int->Clone("eff_Int");
	eff_Int->Divide(eff_den_Int);

	TH1D *eff_Inc = (TH1D*)eff_num_Inc->Clone("eff_Inc");
	eff_Inc->Divide(eff_den_Inc);

	TH1D *pur_Int = (TH1D*)pur_num_Int->Clone("pur_Int");
	pur_Int->Divide(pur_den_Int);

	TH1D *pur_Inc = (TH1D*)pur_num_Inc->Clone("pur_Inc");
	pur_Inc->Divide(pur_den_Inc);

	TH1D *h_truesliceid_all = (TH1D*)file->Get("h_truesliceid_all");
	TH1D *h_truesliceid_cuts = (TH1D*)file->Get("h_truesliceid_cuts");
	TH1D *h_truesliceid_inelastic_all = (TH1D*)file->Get("h_truesliceid_inelastic_all");
	TH1D *h_truesliceid_inelastic_cuts = (TH1D*)file->Get("h_truesliceid_inelastic_cuts");
	TH1D *h_recosliceid_allevts_cuts = (TH1D*)file->Get("h_recosliceid_allevts_cuts");
	TH1D *h_recosliceid_cuts = (TH1D*)file->Get("h_recosliceid_cuts");
	TH1D *h_recosliceid_inelastic_cuts = (TH1D*)file->Get("h_recosliceid_inelastic_cuts");
	TH1D *h_recosliceid_recoinelastic_cuts = (TH1D*)file->Get("h_recosliceid_recoinelastic_cuts");

	TH1D *hinc = (TH1D*)h_recosliceid_allevts_cuts->Clone("hinc");
	//TH1D *hint = (TH1D*)h_recosliceid_allevts_cuts->Clone("hint");
	TH1D *hint = (TH1D*)h_recosliceid_recoinelastic_cuts->Clone("hint");
	hinc->Multiply(pur_Inc);



	TCanvas *c1 = new TCanvas("c1","c1");
	h_recosliceid_allevts_cuts->SetLineColor(7);
	h_recosliceid_allevts_cuts->SetMarkerColor(7);
	//h_recosliceid_allevts_cuts->SetTitle("All Protons;Reco SliceID;Events");
	h_recosliceid_allevts_cuts->SetTitle(";Reco SliceID;Events");
	h_recosliceid_allevts_cuts->GetYaxis()->SetRangeUser(0,3000);
	h_recosliceid_allevts_cuts->GetYaxis()->SetTitleOffset(1.3);
	h_recosliceid_allevts_cuts->DrawCopy();
	hinc->DrawCopy("same");
	h_recosliceid_cuts->SetLineColor(2);
	h_recosliceid_cuts->Draw("same hist");
	//TLegend *leg1 = new TLegend(0.5,0.6,0.8,0.9);
	//TLegend *leg1 = new TLegend(0.5-0.3,0.6,0.8-0.3,0.9);
	TLegend *leg1 = new TLegend(0.14,0.55,0.64,0.9);
	leg1->SetFillStyle(0);
	leg1->AddEntry(h_recosliceid_allevts_cuts,"Selected","ple");
	leg1->AddEntry(hinc, "Selected #times purity","ple");
	leg1->AddEntry(h_recosliceid_cuts,"Selected true protons","l");
	leg1->Draw();

	TCanvas *c2 = new TCanvas("c2","c2");
	h_truesliceid_uf->SetLineColor(7);
	h_truesliceid_uf->SetMarkerColor(7);
	//h_truesliceid_uf->SetTitle("All Protons; True SliceID; Events");
	h_truesliceid_uf->SetTitle("; True SliceID; Events");
	h_truesliceid_uf->GetYaxis()->SetRangeUser(0,4200);
	h_truesliceid_uf->GetYaxis()->SetTitleOffset(1.3);
	h_truesliceid_uf->Draw();
	//hinc->SetLineColor(3);
	//hinc->SetMarkerColor(3);
	hinc->DrawCopy("same");
	//  eff_den_Inc->SetLineColor(2);
	//  eff_den_Inc->SetMarkerColor(2);
	//  eff_den_Inc->Draw("same hist");
	h_truesliceid_all->SetLineColor(2);
	h_truesliceid_all->SetLineWidth(4);
	h_truesliceid_all->SetMarkerColor(2);
	h_truesliceid_all->Draw("same hist");
	//TLegend *leg2 = new TLegend(0.5,0.6,0.8,0.9);
	//TLegend *leg2 = new TLegend(0.5,0.6,0.8,0.9);
	TLegend *leg2 = new TLegend(0.25,0.6,0.75,0.9);
	leg2->SetFillStyle(0);
	leg2->AddEntry(hinc, "Selected #times purity","ple");
	leg2->AddEntry(h_truesliceid_uf,"Unfolded protons","ple");
	leg2->AddEntry(h_truesliceid_all,"True protons","l");
	leg2->Draw();



	hint->Multiply(pur_Int);
	float ymin_recoinel=0;
	float ymax_recoinel=1900;
	TCanvas *c3 = new TCanvas("c3","c3");
	//h_recosliceid_allevts_cuts->SetLineColor(4);
	//h_recosliceid_allevts_cuts->SetMarkerColor(4);
	//h_recosliceid_allevts_cuts->SetTitle("Proton Inelastic Scatterings;Reco SliceID;Events");
	//h_recosliceid_allevts_cuts->Draw();

	h_recosliceid_recoinelastic_cuts->SetLineColor(7);
	h_recosliceid_recoinelastic_cuts->SetMarkerColor(7);
	//h_recosliceid_recoinelastic_cuts->SetTitle("Proton Inelastic Scatterings;Reco SliceID;Events");
	h_recosliceid_recoinelastic_cuts->SetTitle(";Reco SliceID;Events");
	h_recosliceid_recoinelastic_cuts->GetYaxis()->SetTitleOffset(1.3);
	h_recosliceid_recoinelastic_cuts->GetYaxis()->SetRangeUser(ymin_recoinel,ymax_recoinel);
	h_recosliceid_recoinelastic_cuts->Draw();

	h_recosliceid_inelastic_cuts->SetLineColor(2);
	h_recosliceid_inelastic_cuts->SetLineWidth(4);
	//h_recosliceid_inelastic_cuts->SetFillColor(2);
	//h_recosliceid_inelastic_cuts->SetFillStyle(3354);
	h_recosliceid_inelastic_cuts->Draw("same hist");
	hint->DrawCopy("same");
	//TLegend *leg3 = new TLegend(0.5,0.6,0.8,0.9);
	//TLegend *leg3 = new TLegend(0.5,0.6,0.8,0.9);
	TLegend *leg3 = new TLegend(0.45,0.55,0.88,0.9);
	leg3->SetFillStyle(0);
	leg3->AddEntry(h_recosliceid_allevts_cuts,"Selected","ple");
	leg3->AddEntry(hint, "Selected #times purity","ple");
	leg3->AddEntry(h_recosliceid_inelastic_cuts,"Selected true protons","l");
	leg3->Draw();

	TCanvas *c4 = new TCanvas("c4","c4");
	h_truesliceid_inelastic_uf->SetLineColor(7);
	h_truesliceid_inelastic_uf->SetMarkerColor(7);
	//h_truesliceid_inelastic_uf->SetTitle("Proton Inelastic Scatterings; True SliceID; Events");
	h_truesliceid_inelastic_uf->SetTitle("; True SliceID; Events");
	h_truesliceid_inelastic_uf->GetYaxis()->SetRangeUser(0,3800);
	h_truesliceid_inelastic_uf->GetYaxis()->SetTitleOffset(1.3);
	h_truesliceid_inelastic_uf->Draw();
	//hint->SetLineColor(3);
	//hint->SetMarkerColor(3);
	h_truesliceid_inelastic_all->SetLineColor(2);
	h_truesliceid_inelastic_all->SetLineWidth(4);
	h_truesliceid_inelastic_all->SetMarkerColor(2);
	h_truesliceid_inelastic_all->Draw("same hist");
	hint->DrawCopy("same");

	//TLegend *leg4 = new TLegend(0.5,0.6,0.8,0.9);
	//TLegend *leg4 = new TLegend(0.5,0.6,0.8,0.9);
	TLegend *leg4 = new TLegend(0.35,0.56,0.85,0.9);
	leg4->SetFillStyle(0);
	leg4->AddEntry(hint, "Selected #times purity","ple");
	leg4->AddEntry(h_truesliceid_inelastic_uf,"Unfolded protons","ple");
	leg4->AddEntry(h_truesliceid_inelastic_all,"True protons","l");
	leg4->Draw();

	double Ninc[nthinslices] = {0};
	double Nint[nthinslices] = {0};
	double err_inc[nthinslices] = {0};
	double err_int[nthinslices] = {0};
	double SliceID[nthinslices] = {0};

	for (int i = 0; i<nthinslices; ++i){
		SliceID[i] = i;
		Nint[i] = h_truesliceid_inelastic_uf->GetBinContent(i+2);
		err_int[i] = h_truesliceid_inelastic_uf->GetBinError(i+2);
		for (int j = i; j<=nthinslices; ++j){
			Ninc[i] += h_truesliceid_uf->GetBinContent(j+2);
			err_inc[i] += pow(h_truesliceid_uf->GetBinError(j+2),2);
			//cout<<i<<" "<<j<<" "<<h_truesliceid_uf->GetBinContent(j+2)<<" "<<Ninc[i]<<endl;
		}
		err_inc[i] = sqrt(err_inc[i]);
	}

	TGraphErrors *gr_inc = new TGraphErrors(nthinslices, SliceID, Ninc, 0, err_inc);
	TGraphErrors *gr_int = new TGraphErrors(nthinslices, SliceID, Nint, 0, err_int);
	TCanvas *c6 = new TCanvas("c6","c6");
	gr_inc->SetTitle("");
	gr_inc->SetLineWidth(2);
	gr_inc->SetLineColor(7);
	gr_inc->SetMarkerColor(7);
	gr_inc->GetXaxis()->SetTitle("Slice ID");
	gr_inc->GetYaxis()->SetTitle("N_{Inc}");
	gr_inc->GetXaxis()->SetRangeUser(-1, 23);
	gr_inc->Draw("ape");

	TCanvas *c7 = new TCanvas("c7","c7");
	gr_int->SetTitle("");
	gr_int->SetLineWidth(2);
	gr_int->SetLineColor(7);
	gr_int->SetMarkerColor(7);
	gr_int->GetXaxis()->SetTitle("Slice ID");
	gr_int->GetYaxis()->SetTitle("N_{Int}");
	gr_int->GetXaxis()->SetRangeUser(-1, 23);
	gr_int->Draw("ape");

	double NA=6.02214076e23;
	double MAr=39.95; //gmol
	double Density = 1.4; // g/cm^3

	double xs[nthinslices] = {0};
	double err_xs[nthinslices] = {0};
	double incE[nthinslices] = {0};

	//TGraphErrors *gr_trueincE = (TGraphErrors*)file->Get("gr_trueincE"); //use truth KE as KE definition
	TGraphErrors *gr_trueincE = (TGraphErrors*)file->Get("gr_recoincE"); //use reco KE as KE definition
	TGraphErrors *gr_truexs = (TGraphErrors*)file->Get("gr_truexs");
	for (int i = 0; i<nthinslices; ++i){
		////err_xs[i] = MAr/(Density*NA*thinslicewidth)*1e27*sqrt(N_int[i]+pow(N_int[i],2)/N_inc[i])/N_incidents[i];
		////incE[i] = gr_trueincE->GetPointY(i);
		//xs[i] = MAr/(Density*NA*thinslicewidth/true_AngCorr->GetMean())*log(Ninc[i]/(Ninc[i]-Nint[i]))*1e27;
		//err_xs[i] = MAr/(Density*NA*thinslicewidth/true_AngCorr->GetMean())*1e27*sqrt(pow(Nint[i]*err_inc[i]/Ninc[i]/(Ninc[i]-Nint[i]),2)+pow(err_int[i]/(Ninc[i]-Nint[i]),2));
		xs[i] = MAr/(Density*NA*thinslicewidth/1.)*log(Ninc[i]/(Ninc[i]-Nint[i]))*1e27;
		err_xs[i] = MAr/(Density*NA*thinslicewidth/1.)*1e27*sqrt(pow(Nint[i]*err_inc[i]/Ninc[i]/(Ninc[i]-Nint[i]),2)+pow(err_int[i]/(Ninc[i]-Nint[i]),2));
		incE[i] = gr_trueincE->GetY()[i];
		//std::cout<<i<<" "<<Ninc[i]<<" "<<Nint[i]<<" "<<xs[i]<<" "<<incE[i]<<std::endl;
	}




	TFile f2("/dune/data2/users/hyliao/GeantReweight/xs_cascade/proton_cross_section.root");
	TGraph *total_inel_KE = (TGraph*)f2.Get("inel_KE");

	TGraphErrors *gr_recoxs = new TGraphErrors(nthinslices, incE, xs, 0, err_xs);
	gStyle->SetPadLeftMargin(0.13);
	TCanvas *c5 = new TCanvas("c5", "c5");
	float ymax=1500;
	float xmax=450;
	float xmin=0;
	TH2D *f2d=new TH2D("f2d","",450,xmin,xmax,ymax,0,ymax);
	f2d->GetXaxis()->SetTitle("Proton Kinetic Energy [MeV]");
	f2d->GetYaxis()->SetTitleOffset(1.3);
	//f2d->GetYaxis()->SetTitle("#sigma_{inelastic} [mb]");
	f2d->GetYaxis()->SetTitle("P-Ar inelastic cross section [mb]");
	f2d->Draw("");

	//gr_recoxs->SetTitle("Proton Inelastic Cross Section");
	//gr_recoxs->GetXaxis()->SetTitle("Proton Kinetic Energy (MeV)");
	//gr_recoxs->GetYaxis()->SetTitle("#sigma_{inelastic} (mb)");
	total_inel_KE->SetLineColor(2);
	total_inel_KE->Draw("c same");
	gr_recoxs->SetLineWidth(2);
	//gr_recoxs->SetEndErrorSize(0);
	gr_recoxs->Draw("p same");
	//gr_recoxs->Draw("ape same");
	gr_truexs->SetMarkerColor(3);
	gr_truexs->SetLineColor(3);
	//gr_truexs->Draw("pe same");
	gr_truexs->Draw("p same");

	TLegend *leg5 = new TLegend(0.3,0.55,0.8,0.85);
	leg5->SetFillStyle(0);
	leg5->AddEntry(gr_recoxs, "MC with reco information", "pe");
	leg5->AddEntry(gr_truexs, "MC with truth information", "pe");
	//leg5->AddEntry(total_inel_KE, "Geant4 v4_10_6_p01c", "l");
	leg5->AddEntry(total_inel_KE, "Geant4", "l");
	leg5->Draw();

   	//pDUNE Logo
        float logo_y=ymax+15;
        TLatex **txt_pdune1=new TLatex*[1];
        txt_pdune1[0]=new TLatex(xmin, logo_y, Form("#bf{DUNE:ProtoDUNE-SP}"));
        txt_pdune1[0]->SetTextColor(1);
        txt_pdune1[0]->Draw();
        //
        //Beam Logo
        TLatex **txt_p1=new TLatex*[1];
        //txt_p1[0]=new TLatex(xmax-6.3, logo_y, Form("Protons (1 GeV/c)")); //x:0-20
        txt_p1[0]=new TLatex(xmax-150, logo_y, Form("Protons (1 GeV/c)")); //x:0-40
        txt_p1[0]->SetTextColor(1);
        txt_p1[0]->SetTextSize(0.05);
        txt_p1[0]->Draw();


	//c1->Print("plots/xs_sliceidinc_reco.eps");
	//c2->Print("plots/xs_sliceidinc_true.eps");
	//c3->Print("plots/xs_sliceidint_reco.eps");
	//c4->Print("plots/xs_sliceidint_true.eps");
	//c5->Print("plots/xs_pi+inel.eps");
	//c6->Print("plots/xs_Ninc.eps");
	//c7->Print("plots/xs_Nint.eps");

	c1->Print(Form("%s/xs_sliceidinc_reco.eps",outpath.Data()));
	c2->Print(Form("%s/xs_sliceidinc_true.eps",outpath.Data()));
	c3->Print(Form("%s/xs_sliceidint_reco.eps",outpath.Data()));
	c4->Print(Form("%s/xs_sliceidint_true.eps",outpath.Data()));
	c5->Print(Form("%s/xs_p_inel.eps",outpath.Data()));
	c6->Print(Form("%s/xs_Ninc.eps",outpath.Data()));
	c7->Print(Form("%s/xs_Nint.eps",outpath.Data()));

}

