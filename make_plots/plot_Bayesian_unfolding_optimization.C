#include "TVector3.h"
#include <TH2.h>
#include <TH1.h>
#include "TH2D.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TLine.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH3F.h>
#include <TString.h>
#include <TProfile2D.h>
#include <THStack.h>
#include "TGraph.h" 
#include "TGraphSmooth.h" 
#include "TParameter.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TSystem.h"
#include "string"
#include "vector"
#include "TSpline.h"
#include "TH3F.h"
#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>

#include <stdio.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

using namespace std;
using namespace ROOT::Math;

R__LOAD_LIBRARY(libRooUnfold.so) //load share lib

void plot_Bayesian_unfolding_optimization() {

	//load data
	TString outpath="./plots_Bayesian_Lcurve/";
        TString fmc="../prod4areco2_mc_ESliceE_dE20MeV_40slcs_beamxy_nobmrw_kebeamff_v09_39_01_test_only.root";
        TString fdata="../prod4areco2_mc_ESliceE_dE20MeV_40slcs_beamxy_nobmrw_kebeamff_v09_39_01_test_only.root";
        TString fmc_valid="../prod4areco2_mc_ESliceE_dE20MeV_40slcs_beamxy_nobmrw_kebeamff_v09_39_01_valid_only.root";

	
	//reco string pre-fix --------------------------------------//
	TString str_inc=Form("h_recosliceid_allevts_cuts");
	TString str_st_inc=Form("h_reco_st_sliceid_allevts_cuts");
	TString str_int=Form("h_recosliceid_recoinelastic_cuts");

	//true string pre-fix -----------------------------------------//
	TString str_inc_true=Form("h_truesliceid_allevts_cuts");
	TString str_st_inc_true=Form("h_true_st_sliceid_allevts_cuts");
	TString str_int_true=Form("h_truesliceid_recoinelastic_cuts");

        //plot style --------------------------------------------------------------------//
        gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
        gROOT->SetStyle("protoDUNEStyle");
        gROOT->ForceStyle();
        gStyle->SetTitleX(0.5);
        gStyle->SetTitleAlign(23);
        gStyle->SetOptStat(0);

	//read data ------------------------------------------------------------------------------------------//
	TFile *f_data = TFile::Open(fdata.Data());
	TH1D *data_inc=(TH1D*)f_data->Get(str_inc.Data()); //recosliceID after beam quality cuts
	TH1D *data_st_inc=(TH1D*)f_data->Get(str_st_inc.Data()); //reco_st_sliceID after beam quality cuts
	//TH1D *data_int=(TH1D*)f_data->Get("h_recosliceid_inelastic_cuts"); //h_recosliceid_inelastic_cuts
	TH1D *data_int=(TH1D*)f_data->Get(str_int.Data()); //h_recosliceid_inelastic_cuts
	data_inc->SetName("data_inc");	
	data_st_inc->SetName("data_st_inc");	
	data_int->SetName("data_int");	

	int n_data_inc=data_inc->Integral();
	int n_data_st_inc=data_st_inc->Integral();
	int n_data_int=data_int->Integral();

	//read mc [after bmrw] ------------------------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fmc.Data());

	//get truth inc & int
	//TH1D *mc_true_incidents=(TH1D *)f_mc->Get("h_true_incidents");
	//TH1D *mc_true_st_incidents=(TH1D *)f_mc->Get("h_true_st_incidents");
	//TH1D *mc_true_interactions=(TH1D *)f_mc->Get("h_true_interactions");

	//get truesliceID of int &inc
	TH1D *mc_truesliceID_inel=(TH1D *)f_mc->Get("h_truesliceid_inelastic_all"); //IsPureInEL
	TH1D *mc_truesliceID_all=(TH1D *)f_mc->Get("h_truesliceid_all"); //all protons
	TH1D *mc_true_st_sliceID_all=(TH1D *)f_mc->Get("h_true_st_sliceid_all"); //all protons

	//get mc reco slice IDs
	//inc
	TH1D* mc_inc_all=(TH1D*)f_mc->Get(Form("%s",str_inc.Data()));
	TH1D* mc_inc_inel=(TH1D*)f_mc->Get(Form("%s_inel",str_inc.Data()));
	TH1D* mc_inc_el=(TH1D*)f_mc->Get(Form("%s_el",str_inc.Data()));
	TH1D* mc_inc_midcosmic=(TH1D*)f_mc->Get(Form("%s_midcosmic",str_inc.Data()));
	TH1D* mc_inc_midpi=(TH1D*)f_mc->Get(Form("%s_midpi",str_inc.Data()));
	TH1D* mc_inc_midp=(TH1D*)f_mc->Get(Form("%s_midp",str_inc.Data()));
	TH1D* mc_inc_midmu=(TH1D*)f_mc->Get(Form("%s_midmu",str_inc.Data()));
	TH1D* mc_inc_mideg=(TH1D*)f_mc->Get(Form("%s_mideg",str_inc.Data()));
	TH1D* mc_inc_midother=(TH1D*)f_mc->Get(Form("%s_midother",str_inc.Data()));

	//inc_st
	TH1D* mc_st_inc_all=(TH1D*)f_mc->Get(Form("%s",str_st_inc.Data()));
	TH1D* mc_st_inc_inel=(TH1D*)f_mc->Get(Form("%s_inel",str_st_inc.Data()));
	TH1D* mc_st_inc_el=(TH1D*)f_mc->Get(Form("%s_el",str_st_inc.Data()));
	TH1D* mc_st_inc_midcosmic=(TH1D*)f_mc->Get(Form("%s_midcosmic",str_st_inc.Data()));
	TH1D* mc_st_inc_midpi=(TH1D*)f_mc->Get(Form("%s_midpi",str_st_inc.Data()));
	TH1D* mc_st_inc_midp=(TH1D*)f_mc->Get(Form("%s_midp",str_st_inc.Data()));
	TH1D* mc_st_inc_midmu=(TH1D*)f_mc->Get(Form("%s_midmu",str_st_inc.Data()));
	TH1D* mc_st_inc_mideg=(TH1D*)f_mc->Get(Form("%s_mideg",str_st_inc.Data()));
	TH1D* mc_st_inc_midother=(TH1D*)f_mc->Get(Form("%s_midother",str_st_inc.Data()));

	//int
	TH1D* mc_int_all=(TH1D*)f_mc->Get(Form("%s",str_int.Data()));
	TH1D* mc_int_inel=(TH1D*)f_mc->Get(Form("%s_inel",str_int.Data()));
	TH1D* mc_int_el=(TH1D*)f_mc->Get(Form("%s_el",str_int.Data()));
	TH1D* mc_int_midcosmic=(TH1D*)f_mc->Get(Form("%s_midcosmic",str_int.Data()));
	TH1D* mc_int_midpi=(TH1D*)f_mc->Get(Form("%s_midpi",str_int.Data()));
	TH1D* mc_int_midp=(TH1D*)f_mc->Get(Form("%s_midp",str_int.Data()));
	TH1D* mc_int_midmu=(TH1D*)f_mc->Get(Form("%s_midmu",str_int.Data()));
	TH1D* mc_int_mideg=(TH1D*)f_mc->Get(Form("%s_mideg",str_int.Data()));
	TH1D* mc_int_midother=(TH1D*)f_mc->Get(Form("%s_midother",str_int.Data()));

	mc_inc_inel->SetFillColor(2); mc_inc_inel->SetLineColor(2);
	//mc_inc_el->SetFillColor(4); mc_inc_el->SetLineColor(4);
	mc_inc_el->SetFillColor(62); mc_inc_el->SetLineColor(62);
	mc_inc_midp->SetFillColor(3); mc_inc_midp->SetLineColor(3);
	mc_inc_midcosmic->SetFillColor(5); mc_inc_midcosmic->SetLineColor(5);
	mc_inc_midpi->SetFillColor(6); mc_inc_midpi->SetLineColor(6);
	mc_inc_midmu->SetFillColor(28); mc_inc_midmu->SetLineColor(28);
	mc_inc_mideg->SetFillColor(30); mc_inc_mideg->SetLineColor(30);
	mc_inc_midother->SetFillColor(15); mc_inc_midother->SetLineColor(15);

	mc_st_inc_inel->SetFillColor(2); mc_st_inc_inel->SetLineColor(2);
	//mc_st_inc_el->SetFillColor(4); mc_st_inc_el->SetLineColor(4);
	mc_st_inc_el->SetFillColor(62); mc_st_inc_el->SetLineColor(62);
	mc_st_inc_midp->SetFillColor(3); mc_st_inc_midp->SetLineColor(3);
	mc_st_inc_midcosmic->SetFillColor(5); mc_st_inc_midcosmic->SetLineColor(5);
	mc_st_inc_midpi->SetFillColor(6); mc_st_inc_midpi->SetLineColor(6);
	mc_st_inc_midmu->SetFillColor(28); mc_st_inc_midmu->SetLineColor(28);
	mc_st_inc_mideg->SetFillColor(30); mc_st_inc_mideg->SetLineColor(30);
	mc_st_inc_midother->SetFillColor(15); mc_st_inc_midother->SetLineColor(15);

	int n_mc_inc_inel=mc_inc_inel->Integral();
	int n_mc_inc_el=mc_inc_el->Integral();
	int n_mc_inc_midcosmic=mc_inc_midcosmic->Integral();
	int n_mc_inc_midpi=mc_inc_midpi->Integral();
	int n_mc_inc_midp=mc_inc_midp->Integral();
	int n_mc_inc_midmu=mc_inc_midmu->Integral();
	int n_mc_inc_mideg=mc_inc_mideg->Integral();
	int n_mc_inc_midother=mc_inc_midother->Integral();
	int n_mc_inc=n_mc_inc_inel+n_mc_inc_el+n_mc_inc_midcosmic+n_mc_inc_midpi+n_mc_inc_midp+n_mc_inc_midmu+n_mc_inc_mideg+n_mc_inc_midother;
	double norm_mc_inc=(double)n_data_inc/(double)n_mc_inc;
	cout<<"norm_mc_inc="<<norm_mc_inc<<"="<<"n_data_inc:"<<n_data_inc<<"/"<<"n_mc_inc:"<<n_mc_inc<<endl;
	cout<<"INC purity (El+Inel)/(all except misid:p):"<<100.*(float)(n_mc_inc_inel+n_mc_inc_el)/(float)(n_mc_inc)<<endl;
	cout<<"n_mc_inc_midp:"<<n_mc_inc_midp<<endl;
	cout<<"n_mc_inc_inel+n_mc_inc_el="<<n_mc_inc_inel+n_mc_inc_el<<endl;

	int n_mc_st_inc_inel=mc_st_inc_inel->Integral();
	int n_mc_st_inc_el=mc_st_inc_el->Integral();
	int n_mc_st_inc_midcosmic=mc_st_inc_midcosmic->Integral();
	int n_mc_st_inc_midpi=mc_st_inc_midpi->Integral();
	int n_mc_st_inc_midp=mc_st_inc_midp->Integral();
	int n_mc_st_inc_midmu=mc_st_inc_midmu->Integral();
	int n_mc_st_inc_mideg=mc_st_inc_mideg->Integral();
	int n_mc_st_inc_midother=mc_st_inc_midother->Integral();
	int n_mc_st_inc=n_mc_st_inc_inel+n_mc_st_inc_el+n_mc_st_inc_midcosmic+n_mc_st_inc_midpi+n_mc_st_inc_midp+n_mc_st_inc_midmu+n_mc_st_inc_mideg+n_mc_st_inc_midother;
	double norm_mc_st_inc=(double)n_data_st_inc/(double)n_mc_st_inc;
	cout<<"norm_mc_st_inc="<<norm_mc_st_inc<<"="<<"n_data_inc:"<<n_data_st_inc<<"/"<<"n_mc_st_inc:"<<n_mc_st_inc<<endl;
	cout<<"INC purity (El+Inel)/(all except misid:p):"<<100.*(float)(n_mc_st_inc_inel+n_mc_st_inc_el)/(float)(n_mc_st_inc)<<endl;
	cout<<"n_mc_st_inc_midp:"<<n_mc_st_inc_midp<<endl;
	cout<<"n_mc_st_inc_inel+n_mc_st_inc_el="<<n_mc_st_inc_inel+n_mc_st_inc_el<<endl;

	mc_inc_inel->Scale(norm_mc_inc);
	mc_inc_el->Scale(norm_mc_inc);
	mc_inc_midcosmic->Scale(norm_mc_inc);
	mc_inc_midpi->Scale(norm_mc_inc);
	mc_inc_midp->Scale(norm_mc_inc);
	mc_inc_midmu->Scale(norm_mc_inc);
	mc_inc_mideg->Scale(norm_mc_inc);
	mc_inc_midother->Scale(norm_mc_inc);

	mc_st_inc_inel->Scale(norm_mc_st_inc);
	mc_st_inc_el->Scale(norm_mc_st_inc);
	mc_st_inc_midcosmic->Scale(norm_mc_st_inc);
	mc_st_inc_midpi->Scale(norm_mc_st_inc);
	mc_st_inc_midp->Scale(norm_mc_st_inc);
	mc_st_inc_midmu->Scale(norm_mc_st_inc);
	mc_st_inc_mideg->Scale(norm_mc_st_inc);
	mc_st_inc_midother->Scale(norm_mc_st_inc);


	mc_int_inel->SetFillColor(2); mc_int_inel->SetLineColor(2);
	//mc_int_el->SetFillColor(4); mc_int_el->SetLineColor(4);
	mc_int_el->SetFillColor(62); mc_int_el->SetLineColor(62);
	mc_int_midp->SetFillColor(3); mc_int_midp->SetLineColor(3);
	mc_int_midcosmic->SetFillColor(5); mc_int_midcosmic->SetLineColor(5);
	mc_int_midpi->SetFillColor(6); mc_int_midpi->SetLineColor(6);
	mc_int_midmu->SetFillColor(28); mc_int_midmu->SetLineColor(28);
	mc_int_mideg->SetFillColor(30); mc_int_mideg->SetLineColor(30);
	mc_int_midother->SetFillColor(15); mc_int_midother->SetLineColor(15);

	int n_mc_int_inel=mc_int_inel->Integral();
	int n_mc_int_el=mc_int_el->Integral();
	int n_mc_int_midcosmic=mc_int_midcosmic->Integral();
	int n_mc_int_midpi=mc_int_midpi->Integral();
	int n_mc_int_midp=mc_int_midp->Integral();
	int n_mc_int_midmu=mc_int_midmu->Integral();
	int n_mc_int_mideg=mc_int_mideg->Integral();
	int n_mc_int_midother=mc_int_midother->Integral();
	int n_mc_int=n_mc_int_inel+n_mc_int_el+n_mc_int_midcosmic+n_mc_int_midpi+n_mc_int_midp+n_mc_int_midmu+n_mc_int_mideg+n_mc_int_midother;
	double norm_mc_int=(double)n_data_int/(double)n_mc_int;
	cout<<"INT purity (inel/all):"<<100.*n_mc_int_inel/n_mc_int<<endl;
	cout<<"INT purity (el/all):"<<100.*n_mc_int_el/n_mc_int<<endl;
	cout<<"INT purity (misidp/all):"<<100.*n_mc_int_midp/n_mc_int<<endl;

	mc_int_inel->Scale(norm_mc_int);
	mc_int_el->Scale(norm_mc_int);
	mc_int_midcosmic->Scale(norm_mc_int);
	mc_int_midpi->Scale(norm_mc_int);
	mc_int_midp->Scale(norm_mc_int);
	mc_int_midmu->Scale(norm_mc_int);
	mc_int_mideg->Scale(norm_mc_int);
	mc_int_midother->Scale(norm_mc_int);

	
	//bkg subtraction ---------------------------------------------------------------------------------------------------------------//
	//data
	TH1D* data_inc_bkgfree=(TH1D *)data_inc->Clone("data_inc_bkgfree"); data_inc_bkgfree->SetName("data_inc_bkgfree");	
	TH1D* data_st_inc_bkgfree=(TH1D *)data_st_inc->Clone("data_st_inc_bkgfree"); data_st_inc_bkgfree->SetName("data_st_inc_bkgfree");	
	TH1D* data_int_bkgfree=(TH1D *)data_int->Clone("data_int_bkgfree"); data_int_bkgfree->SetName("data_int_bkgfree");	

	double scal_fact_misidp=-1;
	//inc
	//data_inc_bkgfree->Add(mc_inc_el, -1);
	//data_inc_bkgfree->Add(mc_inc_midcosmic, -1);
	//data_inc_bkgfree->Add(mc_inc_midpi, -1);
	//data_inc_bkgfree->Add(mc_inc_midmu, -1);
	//data_inc_bkgfree->Add(mc_inc_mideg, -1);
	//data_inc_bkgfree->Add(mc_inc_midother, -1);
	//data_inc_bkgfree->Multiply(pur_inc);
	//data_inc_bkgfree->Add(mc_inc_midp, scal_fact_misidp);
	data_st_inc_bkgfree->Add(mc_st_inc_midp, scal_fact_misidp);

	//
	//Note: Numerical value of errorbar after subtraction is correct, i.e. (s-b)+-sqrt(s+b)	
	
	//int
	data_int_bkgfree->Add(mc_int_el, scal_fact_misidp);
	data_int_bkgfree->Add(mc_int_midp, scal_fact_misidp);
	//data_int_bkgfree->Add(mc_int_midcosmic, -1);
	//data_int_bkgfree->Add(mc_int_midpi, -1);
	//data_int_bkgfree->Add(mc_int_midmu, -1);
	//data_int_bkgfree->Add(mc_int_mideg, -1);
	//data_int_bkgfree->Add(mc_int_midother, -1);
	



	//response matrix from mc -----------------------------------------------------------------------------------------------------//
	TFile *f_mc_valid = TFile::Open(fmc_valid.Data());
	//Response matrix as a 2D-histogram: (x,y)=(measured,truth)
	RooUnfoldResponse *res_inc=(RooUnfoldResponse*)f_mc_valid->Get("response_SliceID_Inc"); res_inc->SetName("res_inc");
	RooUnfoldResponse *res_st_inc=(RooUnfoldResponse*)f_mc_valid->Get("response_st_SliceID_Inc"); res_st_inc->SetName("res_st_inc");
	RooUnfoldResponse *res_int=(RooUnfoldResponse*)f_mc_valid->Get("response_SliceID_Int"); res_int->SetName("res_int");

	//Bayesian Unfolding
	RooUnfoldBayes uf_inc (res_inc, data_inc_bkgfree, 50); //inc
	RooUnfoldBayes uf_st_inc (res_st_inc, data_st_inc_bkgfree, 50); //st_inc
	RooUnfoldBayes uf_int (res_int, data_int_bkgfree, 50); //int

	//unfolding
  	TH1D *data_inc_uf;
  	TH1D *data_st_inc_uf;
  	TH1D *data_int_uf;
	data_inc_uf=(TH1D* )uf_inc.Hreco();
	data_st_inc_uf=(TH1D* )uf_st_inc.Hreco();
  	data_int_uf=(TH1D *)uf_int.Hreco();

	//apply cut for normalization purpose -----------------//
	mc_truesliceID_all->GetXaxis()->SetRangeUser(1,42);
	data_inc_bkgfree->GetXaxis()->SetRangeUser(1,42);

	data_inc_bkgfree->SetLineColor(1);
	mc_truesliceID_all->SetLineColor(2);
	data_inc_uf->SetLineColor(4);

	
        TCanvas *c_inc = new TCanvas("c_inc", "c_inc", 1400, 900);
	c_inc->Divide(1,1);
	c_inc->cd(1)->SetLogy();
	TH2D *f2d_inc=new TH2D("f2d_inc",Form(""), 31, 10, 41, 600, 0, 60000);
	//TH2D *f2d_inc=new TH2D("f2d_inc",Form(""), 42, -1, 41, 100, 0, 10000);
	f2d_inc->SetTitle("Incident Histogram; SliceID; Counts");
	f2d_inc->Draw();
	data_inc_bkgfree->Draw("ep same");
	mc_truesliceID_all->Draw("hist same");
	data_inc_uf->Scale((double)mc_truesliceID_all->Integral()/(double)data_inc_uf->Integral());
	data_inc_uf->Draw("hist same");
		
	c_inc->Print(Form("%sinc.eps",outpath.Data()));

	data_inc_bkgfree->SetLineColor(1);
	mc_truesliceID_all->SetLineColor(2);
	data_inc_uf->SetLineColor(4);


	//Interaction ---------------------------------------------------------------------------------------//
	//apply cut for normalization purpose -----------------//
	mc_truesliceID_inel->GetXaxis()->SetRangeUser(1,42);
	data_int_bkgfree->GetXaxis()->SetRangeUser(1,42);

	data_int_bkgfree->SetLineColor(1);
	mc_truesliceID_inel->SetLineColor(2);
	data_int_uf->SetLineColor(4);

	
        TCanvas *c_int = new TCanvas("c_int", "c_int", 1400, 900);
	c_int->Divide(1,1);
	c_int->cd(1);
	//TH2D *f2d_inc=new TH2D("f2d_inc",Form(""), 42, -1, 41, 600, 0, 60000);
	TH2D *f2d_int=new TH2D("f2d_int",Form(""), 42, -1, 41, 100, 0, 10000);
	f2d_int->SetTitle("Interaction Histogram; SliceID; Counts");
	f2d_int->Draw();
	data_int_bkgfree->Draw("ep same");
	mc_truesliceID_inel->Draw("hist same");
	//data_int_uf->Scale((double)mc_truesliceID_all->Integral()/(double)data_inc_uf->Integral());
	data_int_uf->Draw("hist same");
		
	c_int->Print(Form("%sint.eps",outpath.Data()));



	//get truth spectra
	//get reco spectra
	//unfolding
	//calculate chi^2
	//plot L-curve: chi^2 vs spiky-ness (# of iterations)





}
