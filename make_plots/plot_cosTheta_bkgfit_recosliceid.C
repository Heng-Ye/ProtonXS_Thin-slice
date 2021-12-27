//#include "../headers/BasicAnaFunc.h"
#include <vector>
#include <iostream>
#include "../headers/SliceParams.h"
#include "../headers/TemplateFitter.h"

void plot_cosTheta_bkgfit_recosliceid() {

	//TString fin=Form("../prod4a_bkgstudy_dx4cm_25slcs.root");
	//TString fin_data=Form("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_run5387_prod4reco2_bkgstudy_dx4cm_25slcs.root");
	TString fin=Form("../prod4a_bkgstudy_dx4cm_25slcs_largerbin.root");
	TString fin_data=Form("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_run5387_prod4reco2_bkgstudy_dx4cm_25slcs_largerbin.root");
	TString fout_path=Form("./plots_bkgstudy");

	//TString fin=Form("../prod4a_bkgstudy_dx4cm_25slcs_largerbin_Testsample.root");
	//TString fin_data=Form("../prod4a_bkgstudy_dx4cm_25slcs_largerbin_Validationsample.root");
	//TString fout_path=Form("./plots_bkgstudy_fakedata");
	
	//TString fout_filename=Form("cosTheta");
	//TString rep="";
	//TString title1="";

	//TString rep="h_recosliceid_allevts_cuts";
	//TString title=Form("%s; Reco SliceID; Events", title1.Data());
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

	//get data -------------------------------------------------------------------------------------//
	TFile *f_data = TFile::Open(fin_data.Data());
	TH1D *h_data[nthinslices+2];
	for (int k=0; k<(int)nthinslices+2; ++k) {
		h_data[k]=(TH1D *)f_data->Get(Form("h_cosTheta_Pos_recosliceid_%d",k-1)); 
	}
	TH1D *h_data_all=(TH1D *)f_data->Get(Form("h_cosTheta_Pos_all"));
	TH1D *h_data_all_nosliceidOne=(TH1D *)f_data->Get(Form("h_cosTheta_Pos_all_nosliceidOne"));

	int n_data=h_data_all->Integral();
	cout<<"n_data:"<<n_data<<endl;

	//get mc -----------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fin.Data());

	TH1D *h_mc[nthinslices+2];
	TH1D *h_mc_inel[nthinslices+2];
	TH1D *h_mc_el[nthinslices+2];
	TH1D *h_mc_midcosmic[nthinslices+2];
	TH1D *h_mc_midpi[nthinslices+2];
	TH1D *h_mc_midp[nthinslices+2];
	TH1D *h_mc_midmu[nthinslices+2];
	TH1D *h_mc_mideg[nthinslices+2];
	TH1D *h_mc_midother[nthinslices+2];
	for (int k=0; k<(int)nthinslices+2; ++k) {
		h_mc[k]=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_recosliceid_%d",k-1)); 

		h_mc_inel[k]=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_inel_recosliceid_%d",k-1)); 
		h_mc_el[k]=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_el_recosliceid_%d",k-1)); 
		h_mc_midcosmic[k]=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_midcosmic_recosliceid_%d",k-1)); 
		h_mc_midpi[k]=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_midpi_recosliceid_%d",k-1)); 
		h_mc_midp[k]=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_midp_recosliceid_%d",k-1)); 
		h_mc_midmu[k]=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_midmu_recosliceid_%d",k-1)); 
		h_mc_mideg[k]=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_mideg_recosliceid_%d",k-1)); 
		h_mc_midother[k]=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_midother_recosliceid_%d",k-1)); 

		h_mc_inel[k]->SetFillColor(2); h_mc_inel[k]->SetLineColor(2);
		h_mc_el[k]->SetFillColor(4); h_mc_el[k]->SetLineColor(4);
		h_mc_midp[k]->SetFillColor(3); h_mc_midp[k]->SetLineColor(3);
		h_mc_midcosmic[k]->SetFillColor(5); h_mc_midcosmic[k]->SetLineColor(5);
		h_mc_midpi[k]->SetFillColor(6); h_mc_midpi[k]->SetLineColor(6);
		h_mc_midmu[k]->SetFillColor(28); h_mc_midmu[k]->SetLineColor(28);
		h_mc_mideg[k]->SetFillColor(30); h_mc_mideg[k]->SetLineColor(30);
		h_mc_midother[k]->SetFillColor(15); h_mc_midother[k]->SetLineColor(15);
	}

	TH1D *h_mc_all=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all"));
	TH1D *h_mc_all_inel=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all_inel"));
	TH1D *h_mc_all_el=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all_el"));
	TH1D *h_mc_all_midcosmic=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all_midcosmic"));
	TH1D *h_mc_all_midpi=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all_midpi"));
	TH1D *h_mc_all_midp=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all_midp"));
	TH1D *h_mc_all_midmu=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all_midmu"));
	TH1D *h_mc_all_mideg=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all_mideg"));
	TH1D *h_mc_all_midother=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all_midother"));

	h_mc_all_inel->SetFillColor(2); h_mc_all_inel->SetLineColor(2);
	h_mc_all_el->SetFillColor(4); h_mc_all_el->SetLineColor(4);
	h_mc_all_midp->SetFillColor(3); h_mc_all_midp->SetLineColor(3);
	h_mc_all_midcosmic->SetFillColor(5); h_mc_all_midcosmic->SetLineColor(5);
	h_mc_all_midpi->SetFillColor(6); h_mc_all_midpi->SetLineColor(6);
	h_mc_all_midmu->SetFillColor(28); h_mc_all_midmu->SetLineColor(28);
	h_mc_all_mideg->SetFillColor(30); h_mc_all_mideg->SetLineColor(30);
	h_mc_all_midother->SetFillColor(15); h_mc_all_midother->SetLineColor(15);

	TH1D *h_mc_all_nosliceidOne=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all_nosliceidOne"));
	TH1D *h_mc_all_nosliceidOne_inel=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all_nosliceidOne_inel"));
	TH1D *h_mc_all_nosliceidOne_el=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all_nosliceidOne_el"));
	TH1D *h_mc_all_nosliceidOne_midcosmic=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all_nosliceidOne_midcosmic"));
	TH1D *h_mc_all_nosliceidOne_midpi=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all_nosliceidOne_midpi"));
	TH1D *h_mc_all_nosliceidOne_midp=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all_nosliceidOne_midp"));
	TH1D *h_mc_all_nosliceidOne_midmu=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all_nosliceidOne_midmu"));
	TH1D *h_mc_all_nosliceidOne_mideg=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all_nosliceidOne_mideg"));
	TH1D *h_mc_all_nosliceidOne_midother=(TH1D *)f_mc->Get(Form("h_cosTheta_Pos_all_nosliceidOne_midother"));

	h_mc_all_nosliceidOne_inel->SetFillColor(2); h_mc_all_nosliceidOne_inel->SetLineColor(2);
	h_mc_all_nosliceidOne_el->SetFillColor(4); h_mc_all_nosliceidOne_el->SetLineColor(4);
	h_mc_all_nosliceidOne_midp->SetFillColor(3); h_mc_all_nosliceidOne_midp->SetLineColor(3);
	h_mc_all_nosliceidOne_midcosmic->SetFillColor(5); h_mc_all_nosliceidOne_midcosmic->SetLineColor(5);
	h_mc_all_nosliceidOne_midpi->SetFillColor(6); h_mc_all_nosliceidOne_midpi->SetLineColor(6);
	h_mc_all_nosliceidOne_midmu->SetFillColor(28); h_mc_all_nosliceidOne_midmu->SetLineColor(28);
	h_mc_all_nosliceidOne_mideg->SetFillColor(30); h_mc_all_nosliceidOne_mideg->SetLineColor(30);
	h_mc_all_nosliceidOne_midother->SetFillColor(15); h_mc_all_nosliceidOne_midother->SetLineColor(15);

	int n_mc=h_mc_all->Integral();
	double norm_mc=(double)n_data/(double)n_mc;
	cout<<"n_mc:"<<n_mc<<endl;


	TH1D *MC2_all_nosliceidOne=(TH1D *)h_mc_all_nosliceidOne_midp->Clone();
	TH1D *MC1_all_nosliceidOne=(TH1D *)h_mc_all_nosliceidOne->Clone();
	MC1_all_nosliceidOne->Add(h_mc_all_nosliceidOne_midp,-1);
	MC2_all_nosliceidOne->Scale(norm_mc);
	MC1_all_nosliceidOne->Scale(norm_mc);


	TH1D *MC2_all=(TH1D *)h_mc_all_midp->Clone();
	TH1D *MC1_all=(TH1D *)h_mc_all->Clone();
	MC1_all->Add(h_mc_all_midp,-1);
	MC2_all->Scale(norm_mc);
	MC1_all->Scale(norm_mc);




	//normalization ------------------------------------//
	//slice-by-slice
	int n_bin_roi=1;
	vector<double> sliceid;
	vector<double> data;
	vector<double> err_data;
	vector<double> mc1; //mc(mc except misID:p)
	vector<double> err_mc1; //mc(mc except misID:p)
	vector<double> mc2; //mc misID:p
	vector<double> err_mc2; //mc misID:p

	for (size_t k=0; k<nthinslices+2; ++k) {
		sliceid.push_back((double)k-1);
		double n_data=h_data[k]->GetBinContent(n_bin_roi);
		double err_n_data=h_data[k]->GetBinError(n_bin_roi);
		data.push_back(n_data);
		err_data.push_back(err_n_data);

		double n_mc_inel=h_mc_inel[k]->GetBinContent(n_bin_roi);
		double err_n_mc_inel=h_mc_inel[k]->GetBinError(n_bin_roi);
	
		double n_mc_el=h_mc_el[k]->GetBinContent(n_bin_roi);
		double err_n_mc_el=h_mc_el[k]->GetBinError(n_bin_roi);

		double n_mc_midp=h_mc_midp[k]->GetBinContent(n_bin_roi);
		double err_n_mc_midp=h_mc_midp[k]->GetBinError(n_bin_roi);
	
		double n_mc_midcosmic=h_mc_midcosmic[k]->GetBinContent(n_bin_roi);
		double err_n_mc_midcosmic=h_mc_midcosmic[k]->GetBinError(n_bin_roi);

		double n_mc_midpi=h_mc_midpi[k]->GetBinContent(n_bin_roi);
		double err_n_mc_midpi=h_mc_midpi[k]->GetBinError(n_bin_roi);

		double n_mc_midmu=h_mc_midmu[k]->GetBinContent(n_bin_roi);
		double err_n_mc_midmu=h_mc_midmu[k]->GetBinError(n_bin_roi);

		double n_mc_mideg=h_mc_mideg[k]->GetBinContent(n_bin_roi);
		double err_n_mc_mideg=h_mc_mideg[k]->GetBinError(n_bin_roi);

		double n_mc_midother=h_mc_midother[k]->GetBinContent(n_bin_roi);
		double err_n_mc_midother=h_mc_midother[k]->GetBinError(n_bin_roi);

		mc2.push_back(norm_mc*n_mc_midp);
		err_mc2.push_back(norm_mc*err_n_mc_midp);
		double n_mc1=n_mc_inel+n_mc_el+n_mc_midcosmic+n_mc_midpi+n_mc_midmu+n_mc_mideg+n_mc_midother;
		double err_n_mc1=sqrt(pow(err_n_mc_inel,2)+pow(err_n_mc_el,2)+pow(err_n_mc_midcosmic,2)+pow(err_n_mc_midpi,2)+pow(err_n_mc_midmu,2)+pow(err_n_mc_mideg,2)+pow(err_n_mc_midother,2));
		mc1.push_back(norm_mc*n_mc1);
		err_mc1.push_back(norm_mc*err_n_mc1);

		//mc norm
		h_mc[k]->Scale(norm_mc);
		h_mc_inel[k]->Scale(norm_mc);
		h_mc_el[k]->Scale(norm_mc);
		h_mc_midcosmic[k]->Scale(norm_mc);
		h_mc_midpi[k]->Scale(norm_mc);
		h_mc_midp[k]->Scale(norm_mc);
		h_mc_midmu[k]->Scale(norm_mc);
		h_mc_mideg[k]->Scale(norm_mc);
		h_mc_midother[k]->Scale(norm_mc);
	}
	//all
	h_mc_all->Scale(norm_mc);
	h_mc_all_inel->Scale(norm_mc);
	h_mc_all_el->Scale(norm_mc);
	h_mc_all_midp->Scale(norm_mc);
	h_mc_all_midcosmic->Scale(norm_mc);
	h_mc_all_midpi->Scale(norm_mc);
	h_mc_all_midmu->Scale(norm_mc);
	h_mc_all_mideg->Scale(norm_mc);
	h_mc_all_midother->Scale(norm_mc);

	//all_nosliceidOne
	h_mc_all_nosliceidOne->Scale(norm_mc);
	h_mc_all_nosliceidOne_inel->Scale(norm_mc);
	h_mc_all_nosliceidOne_el->Scale(norm_mc);
	h_mc_all_nosliceidOne_midp->Scale(norm_mc);
	h_mc_all_nosliceidOne_midcosmic->Scale(norm_mc);
	h_mc_all_nosliceidOne_midpi->Scale(norm_mc);
	h_mc_all_nosliceidOne_midmu->Scale(norm_mc);
	h_mc_all_nosliceidOne_mideg->Scale(norm_mc);
	h_mc_all_nosliceidOne_midother->Scale(norm_mc);
	//--------------------------------------------------//


	//cosTheta dists.
	vector<double> r_origin;
	vector<double> err_r_origin;
	vector<double> scal_midp;
	vector<double> err_scal_midp;
	vector<double> r_raw;
	vector<double> err_r_raw;

	for (int k=0; k<(int)nthinslices+2; ++k) {
	
		//stack histogram
		THStack* hs=new THStack("hs","");
		hs->Add(h_mc_inel[k]);
		hs->Add(h_mc_el[k]);
		hs->Add(h_mc_midp[k]);
		hs->Add(h_mc_midcosmic[k]);
		hs->Add(h_mc_midpi[k]);
		hs->Add(h_mc_midmu[k]);
		hs->Add(h_mc_mideg[k]);
		hs->Add(h_mc_midother[k]);

		//data/mc ratio
        	TH1D *R=(TH1D*)h_data[k]->Clone();
		R->Add(h_mc_inel[k], -1);
		R->Add(h_mc_el[k], -1);
		R->Add(h_mc_midcosmic[k], -1);
		R->Add(h_mc_midpi[k], -1);
		R->Add(h_mc_midmu[k], -1);
		R->Add(h_mc_mideg[k], -1);
		R->Add(h_mc_midother[k], -1);
		R->Divide(R, h_mc_midp[k]);

		//data//MC ratio
        	TH1D *R_raw=(TH1D*)h_data[k]->Clone();
        	R_raw->Divide(h_data[k], h_mc[k]);

		r_origin.push_back(R->GetBinContent(n_bin_roi)); //for bin that cosTheta<0.9
		err_r_origin.push_back(R->GetBinError(n_bin_roi));

		r_raw.push_back(R_raw->GetBinContent(n_bin_roi+1)); 
		err_r_raw.push_back(R_raw->GetBinError(n_bin_roi+1));


		//get bin content of cosTheta<0.9
		TH1D *DATA=new TH1D("DATA","", 1, 0, 1);
		TH1D *MC_1=new TH1D("MC_1","", 1, 0, 1); //:mc(mc except misID:p)
		TH1D *MC_2=new TH1D("MC_2","", 1, 0, 1); //mc(misID:p)
		DATA->SetBinContent(1, data.at(k));
		DATA->SetBinError(1, err_data.at(k));

		MC_1->SetBinContent(1, mc1.at(k));
		MC_1->SetBinError(1, err_mc1.at(k));

		MC_2->SetBinContent(1, mc2.at(k));
		MC_2->SetBinError(1, err_mc1.at(k));


		//fitting
  		TemplateFitter fitter;
		
		//Min.(h0-h1-s*h2)
		//h0:data
		//h1:mc(mc except misID:p)
		//h2:mc(misID:p)			
		fitter.SetHistograms(DATA, MC_1, MC_2);
    		//fitter.SetFitRange(4,-1);
    		fitter.Fit();
		scal_midp.push_back(fitter.GetPar());
    		err_scal_midp.push_back(fitter.GetParError());

		//draw fig
		TCanvas *c_ = new TCanvas("c_", "c_", 900,700);
        	TPad *pad1 = new TPad("pad1","pad1",0,0.29,1.,1.);
        	pad1->SetBottomMargin(0.03);
        	pad1->SetGridx();
        	pad1->SetGridy();
        	pad1->SetLogy();
        	pad1->Draw();
        	pad1->cd();

		float xmin=0;
		float xmax=1;
        	//float ymax=90;
        	//float ymax=1000;
        	float ymax=10000;

		TString title=Form("SliceID %d",k-1);

		TH2D *f2d=new TH2D("f2d",Form("%s",""),100,xmin,xmax,3000,0,ymax); //logy
		f2d->SetTitle(title.Data());
		f2d->GetYaxis()->SetTitleOffset(1.);
		f2d->Draw();
		hs->Draw("same hist");
		h_data[k]->Draw("ep same");

		TLegend *leg = new TLegend(0.2,0.65,0.8,0.9);
		leg->SetFillStyle(0);
		//leg->AddEntry(h_data[k], "Data","ep");
		leg->AddEntry(h_data[k], "Fake Data","ep");
		leg->AddEntry(h_mc_inel[k], "Inel","f");
		leg->AddEntry(h_mc_el[k], "El","f");

		leg->AddEntry(h_mc_midcosmic[k],"misID:cosmic","f");
		leg->AddEntry(h_mc_midp[k], "misID:p","f");
		leg->AddEntry(h_mc_midpi[k], "misID:#pi","f");

		leg->AddEntry(h_mc_midmu[k],"misID:#mu","f");
		leg->AddEntry(h_mc_mideg[k], "misID:e/#gamma","f");
		leg->AddEntry(h_mc_midother[k], "misID:other","f");

		leg->SetNColumns(3);
		leg->Draw();

        	c_->cd();
        	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
        	pad2->SetTopMargin(0.02);
        	pad2->SetBottomMargin(0.25);
        	pad2->SetGridx();
        	pad2->SetGridy();
        	pad2->Draw();
        	pad2->cd();

        	//TH2D *f2d2=new TH2D("f2d2",Form("%s","Pos"),100,xmin,xmax,64,-2,6.); //liny
        	TH2D *f2d2=new TH2D("f2d2",Form("%s","Pos"),100,xmin,xmax,64,-1,5.); //liny for fake data
        	//f2d2->SetTitle(";cos#theta;Data/MC");
        	f2d2->SetTitle(";cos#theta; S");
        	f2d2->GetXaxis()->SetLabelSize(0.1);
        	f2d2->GetYaxis()->SetLabelSize(0.1);

        	f2d2->GetXaxis()->SetTitleSize(0.1);
        	f2d2->GetYaxis()->SetTitleSize(0.1);
        	f2d2->GetYaxis()->SetTitleOffset(.5);

        	f2d2->Draw();
        	TLine *line1=new TLine(xmin,1,xmax,1);
        	line1->SetLineColor(7);
        	line1->SetLineWidth(2);
        	line1->Draw("same");
        	R->Draw("ep same");

		c_->Print(Form("%s/cosTheta_sliceID%d.eps",fout_path.Data(), k-1));
		if (k==0) {
		  c_->Print(Form("%s/cosTheta_sliceID.pdf(",fout_path.Data()));
		}
		else {
			if (k==(int)nthinslices+1) {
		  		c_->Print(Form("%s/cosTheta_sliceID.pdf)",fout_path.Data()));
			}
			else {
		  		c_->Print(Form("%s/cosTheta_sliceID.pdf",fout_path.Data()));
			}
		}


		delete hs;
		delete R;
		delete pad1;
		delete f2d;
		delete leg;
		delete pad2;
		delete c_;
		delete DATA;
		delete MC_1;
		delete MC_2;
	} 

	//data/MC ratio after proton correction (cosTheta>0.9) -----------------------------------------------------------------//
	vector<double> r_corr;
	vector<double> err_r_corr;
	vector<double> r_corr_roi;
	vector<double> err_r_corr_roi;
	for (int k=0; k<(int)nthinslices+2; ++k) {
        	TH1D *R=(TH1D*)h_data[k]->Clone();

        	TH1D *MC_corr=(TH1D*)h_mc_inel[k]->Clone();
		MC_corr->Add(h_mc_el[k]);
		MC_corr->Add(h_mc_midcosmic[k]);
		MC_corr->Add(h_mc_midpi[k]);
		MC_corr->Add(h_mc_midmu[k]);
		MC_corr->Add(h_mc_mideg[k]);
		MC_corr->Add(h_mc_midother[k]);

		TH1D *MC_midp_corr=(TH1D*)h_mc_midp[k]->Clone();
		MC_midp_corr->Scale(scal_midp.at(k));
		MC_corr->Add(MC_midp_corr);

		R->Divide(R, MC_corr);

		r_corr.push_back(R->GetBinContent(n_bin_roi)); //cosTheta<=.9
		err_r_corr.push_back(R->GetBinError(n_bin_roi)); 

		r_corr_roi.push_back(R->GetBinContent(n_bin_roi+1)); ////cosTheta>.9
		err_r_corr_roi.push_back(R->GetBinError(n_bin_roi+1));
		
		cout<<"r (cosTheta>0.9):"<<R->GetBinContent(n_bin_roi+1)<<" +- "<<R->GetBinError(n_bin_roi+1)<<endl;
	}	


	cout<<"\n\n\n\nNbin data:"<<h_data[1]->GetNbinsX()<<endl;




	TCanvas *c_r = new TCanvas("c_r", "c_r", 900,1600);
	c_r->Divide(1,4);

	TGraphErrors *gr_r = new TGraphErrors(nthinslices, &sliceid.at(0), &r_origin.at(0), 0, &err_r_origin.at(0));
	TGraphErrors *gr_scal = new TGraphErrors(nthinslices, &sliceid.at(0), &scal_midp.at(0), 0, &err_scal_midp.at(0));
	TGraphErrors *gr_r_corr = new TGraphErrors(nthinslices, &sliceid.at(0), &r_corr.at(0), 0, &err_r_corr.at(0));
	TGraphErrors *gr_r_corr_roi = new TGraphErrors(nthinslices, &sliceid.at(0), &r_corr_roi.at(0), 0, &err_r_corr_roi.at(0));
	TGraphErrors *gr_r_raw = new TGraphErrors(nthinslices, &sliceid.at(0), &r_raw.at(0), 0, &err_r_raw.at(0));

	//gr_scal->SetMarkerColor(2);
	//gr_scal->SetLineColor(2);
	c_r->cd(1);
        gr_r->SetTitle("cos#Theta<0.9; SliceID; S");
	gr_r->Draw("ap");

	c_r->cd(2);
        gr_scal->SetTitle("cos#Theta<0.9; SliceID; Proton Correction (misid:p)");
	gr_scal->Draw("ap");

	c_r->cd(3);
        //gr_r_corr->SetTitle("cos#Theta<0.9; SliceID; Data/MC after Proton Correction");
        gr_r_corr->SetTitle("cos#Theta<0.9; SliceID; Fake Data/MC after Proton Correction");
	gr_r_corr->Draw("ap");

	c_r->cd(4);
	//TH2D *f2d=new TH2D("f2d",Form(""),26,-1, 26);
	gr_r_raw->SetTitle("cos#Theta>0.96; SliceID; Data/MC");
	//gr_r_raw->SetTitle("cos#Theta>=0.9; SliceID; Fake Data/MC");
	gr_r_raw->GetYaxis()->SetRangeUser(0,2.2);
	gr_r_raw->Draw("ap");
	gr_r_corr_roi->SetMarkerColor(2);
	gr_r_corr_roi->SetLineColor(2);
	gr_r_corr_roi->SetMarkerStyle(24);
	//gr_r_corr_roi->SetTitle("cos#Theta>=0.9; SliceID; Data/MC after Proton Correction");
	gr_r_corr_roi->Draw("psame");

	TLine *l=new TLine(-3,1,25,1);
	l->SetLineStyle(2);

	TLegend *leg_r = new TLegend(0.2,0.65,0.8,0.9);
	leg_r->SetFillStyle(0);
	leg_r->AddEntry(gr_r_raw, "Before Proton Correction","ep");
	leg_r->AddEntry(gr_r_corr_roi, "After Proton Correction","ep");
	leg_r->Draw();	
	l->Draw("same");

	c_r->Print(Form("%s/sliceID_Ratio.eps",fout_path.Data()));


	TCanvas *c_r_ = new TCanvas("c_r_", "c_r_", 900, 600);
	c_r_->Divide(1,1);
	gr_r_raw->Draw("ap");
	gr_r_corr_roi->Draw("psame");
	leg_r->Draw();	
	l->Draw("same");
	c_r_->Print(Form("%s/sliceID_Ratio_result.eps",fout_path.Data()));


/*
	//over-all fit ----------------------------------------------------------------------------------------------//
		//h0:data
		//h1:mc(mc except misID:p)
		//h2:mc(misID:p)			
	//fitting A
  	TemplateFitter fitter_all_nosliceidOne;
	fitter_all_nosliceidOne.SetHistograms(h_data_all_nosliceidOne, MC1_all_nosliceidOne, MC2_all_nosliceidOne);
	fitter_all_nosliceidOne.Fit();
	double scal_all_nosliceidOne=fitter_all_nosliceidOne.GetPar();
	double err_scal_all_nosliceidOne=fitter_all_nosliceidOne.GetParError();


	
	cout<<"\n\n\nscal_all_nosliceidOne:"<<scal_all_nosliceidOne<<" +- "<<err_scal_all_nosliceidOne<<"\n\n\n"<<endl;


	//fitting B
  	TemplateFitter fitter_all;
	fitter_all.SetHistograms(h_data_all, MC1_all, MC2_all);
	fitter_all.Fit();
	double scal_all=fitter_all.GetPar();
	double err_scal_all=fitter_all.GetParError();

	cout<<"\n\n\nscal_all:"<<scal_all<<" +- "<<err_scal_all<<endl;
*/


}
