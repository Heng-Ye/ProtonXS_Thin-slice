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

#include "../headers/BetheBloch.h"
#include "../headers/ESliceParams.h"
#include "../headers/BasicParameters.h"
#include "../headers/BasicFunctions.h"


using namespace std;
using namespace ROOT::Math;

R__LOAD_LIBRARY(libRooUnfold.so) //load share lib

	void plot_Bayesian_unfolding_optimization() {

		//load data
		TString outpath="./plots_Bayesian_Lcurve/";
		//TString fmc="../prod4areco2_mc_ESliceE_dE20MeV_40slcs_beamxy_nobmrw_kebeamff_v09_39_01_test_only.root";
		//TString fdata="../prod4areco2_mc_ESliceE_dE20MeV_40slcs_beamxy_nobmrw_kebeamff_v09_39_01_test_only.root";
		//TString fmc_valid="../prod4areco2_mc_ESliceE_dE20MeV_40slcs_beamxy_nobmrw_kebeamff_v09_39_01_valid_only.root";

		//TString fmc="../prod4areco2_mc_ESliceE_dE20MeV_40slcs_beamxy_nobmrw_kebeamff_v09_39_01_All.root";
		//TString fdata="../prod4areco2_mc_ESliceE_dE20MeV_40slcs_beamxy_nobmrw_kebeamff_v09_39_01_All.root";
		//TString fmc_valid="../prod4areco2_mc_ESliceE_dE20MeV_40slcs_beamxy_nobmrw_kebeamff_v09_39_01_All.root";

		TString fmc="../prod4areco2_mc_ESliceE_dE20MeV_30slcs_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold.root";
		TString fdata="../prod4areco2_mc_ESliceE_dE20MeV_30slcs_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold.root";
		TString fmc_valid="../prod4areco2_mc_ESliceE_dE20MeV_30slcs_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold.root";

		//reco string pre-fix --------------------------------------//
		TString str_inc=Form("h_recosliceid_allevts_cuts");
		TString str_st_inc=Form("h_reco_st_sliceid_allevts_cuts");
		TString str_2d_inc=Form("h2d_recosliceid_allevts_cuts");
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

		//read data -------------------------------------------------------------------------------------------------------//
		TFile *f_data = TFile::Open(fdata.Data());
		TH1D *data_inc=(TH1D*)f_data->Get(str_inc.Data()); //recosliceID after beam quality cuts
		TH1D *data_st_inc=(TH1D*)f_data->Get(str_st_inc.Data()); //reco_st_sliceID after beam quality cuts
		TH2D *data2d_inc=(TH2D*)f_data->Get(str_2d_inc.Data()); //reco_st_sliceID vs recosliceID after beam quality cuts
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

		//inc_2d
		TH2D* mc2d_inc_all=(TH2D*)f_mc->Get(Form("%s",str_2d_inc.Data()));
		TH2D* mc2d_inc_inel=(TH2D*)f_mc->Get(Form("%s_inel",str_2d_inc.Data()));
		TH2D* mc2d_inc_el=(TH2D*)f_mc->Get(Form("%s_el",str_2d_inc.Data()));
		TH2D* mc2d_inc_midcosmic=(TH2D*)f_mc->Get(Form("%s_midcosmic",str_2d_inc.Data()));
		TH2D* mc2d_inc_midpi=(TH2D*)f_mc->Get(Form("%s_midpi",str_2d_inc.Data()));
		TH2D* mc2d_inc_midp=(TH2D*)f_mc->Get(Form("%s_midp",str_2d_inc.Data()));
		TH2D* mc2d_inc_midmu=(TH2D*)f_mc->Get(Form("%s_midmu",str_2d_inc.Data()));
		TH2D* mc2d_inc_mideg=(TH2D*)f_mc->Get(Form("%s_mideg",str_2d_inc.Data()));
		TH2D* mc2d_inc_midother=(TH2D*)f_mc->Get(Form("%s_midother",str_2d_inc.Data()));

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

		//2d_inc
		int n_mc2d_inc_inel=mc2d_inc_inel->Integral();
		int n_mc2d_inc_el=mc2d_inc_el->Integral();
		int n_mc2d_inc_midcosmic=mc2d_inc_midcosmic->Integral();
		int n_mc2d_inc_midpi=mc2d_inc_midpi->Integral();
		int n_mc2d_inc_midp=mc2d_inc_midp->Integral();
		int n_mc2d_inc_midmu=mc2d_inc_midmu->Integral();
		int n_mc2d_inc_mideg=mc2d_inc_mideg->Integral();
		int n_mc2d_inc_midother=mc2d_inc_midother->Integral();
		int n_mc2d_inc=n_mc2d_inc_inel+n_mc2d_inc_el+n_mc2d_inc_midcosmic+n_mc2d_inc_midpi+n_mc2d_inc_midp+n_mc2d_inc_midmu+n_mc2d_inc_mideg+n_mc2d_inc_midother;

		double norm_mc2d_inc=(double)data2d_inc->Integral()/(double)n_mc2d_inc;
		//double norm_mc2d_inc=(double)data2d_inc->Integral()/(double)mc2d_inc_all->Integral();
		mc2d_inc_inel->Scale(norm_mc2d_inc);
		mc2d_inc_el->Scale(norm_mc2d_inc);
		mc2d_inc_midcosmic->Scale(norm_mc2d_inc);
		mc2d_inc_midpi->Scale(norm_mc2d_inc);
		mc2d_inc_midp->Scale(norm_mc2d_inc);
		mc2d_inc_midmu->Scale(norm_mc2d_inc);
		mc2d_inc_mideg->Scale(norm_mc2d_inc);
		mc2d_inc_midother->Scale(norm_mc2d_inc);
		mc2d_inc_all->Scale(norm_mc2d_inc);


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
		TH2D* data2d_inc_bkgfree=(TH2D *)data2d_inc->Clone("data2d_inc_bkgfree"); data2d_inc_bkgfree->SetName("data2d_inc_bkgfree");	
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
		data2d_inc_bkgfree->Add(mc2d_inc_midp, scal_fact_misidp);

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
		RooUnfoldResponse *res2d_inc=(RooUnfoldResponse*)f_mc_valid->Get("response_SliceID_2D"); res2d_inc->SetName("res2d_inc");
		RooUnfoldResponse *res_int=(RooUnfoldResponse*)f_mc_valid->Get("response_SliceID_Int"); res_int->SetName("res_int");


		//Some constants ----------------------------------------------//
		double xs_const=MAr/(Density*NA*thinslicewidth)*1e27;
		//[0]KE estimation
		BetheBloch BB(2212);
		double KE[nthinslices] = {0};
		double err_KE[nthinslices] = {0};
		double dEdx[nthinslices] = {0};

		bool KE_selection[nthinslices] = {0};
		double KE_thr=420;
		int NDF=0;
		for (int i = 0; i<nthinslices; ++i) {
			KE[i]=Emax-((double)i+0.5)*thinslicewidth; //av_KE
			//KE[i]=true_incE[i]->GetMean();
			err_KE[i]=(double)thinslicewidth/2.;
			dEdx[i]=BB.meandEdx(KE[i]); // MeV/cm

			cout<<"["<<i<<"] KE="<<KE[i]<<endl;

			if (KE[i]<KE_thr) {
				KE_selection[i]=1;
				NDF++;
			}
			else {
				KE_selection[i]=0;
			}
		}
		//-------------------------------------------------------------//

		//true xs --------------------------------------
		double sliceid[nthinslices] = {0};	
		double true_inc[nthinslices]={0};
		double true_int[nthinslices]={0};
		double err_true_inc[nthinslices]={0};
		double err_true_int[nthinslices]={0};
		double true_xs[nthinslices] = {0};
		double err_true_xs[nthinslices] = {0};
		vector<double> vec_true_xs;
		vector<double> vec_err_true_xs;
		vector<double> vec_ke;
		vector<double> vec_zero;
		for (int i = 0; i<nthinslices; ++i){
			sliceid[i]=i+.5;

			//TRUE ------------------------------------------------------------------------------------------
			//[1a]true inc/int from truth sliceID dists.
			true_int[i] = mc_truesliceID_inel->GetBinContent(i+2);
			err_true_int[i] = mc_truesliceID_inel->GetBinError(i+2);

			//new treatment for INC ----------------------------------------------------
			for (int j=0; j<=i; ++j){
				true_inc[i]+=mc_true_st_sliceID_all->GetBinContent(j+2);
				err_true_inc[i]+=pow(mc_true_st_sliceID_all->GetBinError(j+2),2);
			}

			for (int j=0; j<=i-1; ++j){
				true_inc[i]-=mc_truesliceID_all->GetBinContent(j+2);
				err_true_inc[i]+=pow(mc_truesliceID_all->GetBinError(j+2),2);
			}
			//--------------------------------------------------------------------------

			err_true_inc[i] = sqrt(err_true_inc[i]);

			true_xs[i]=xs_const*dEdx[i]*log(true_inc[i]/(true_inc[i]-true_int[i]));
			err_true_xs[i]=xs_const*dEdx[i]*sqrt(true_int[i]+pow(true_int[i],2)/true_inc[i])/true_inc[i];

			//for chi2 calc -----------------------------------
			if (KE_selection[i]==1) {
				vec_true_xs.push_back(true_xs[i]);
				vec_err_true_xs.push_back(err_true_xs[i]);
				vec_ke.push_back(KE[i]);
				vec_zero.push_back(0);
			}
		}	
		TGraphErrors *gr_truexs = new TGraphErrors(vec_ke.size(), &vec_ke.at(0), &vec_true_xs.at(0), &vec_zero.at(0), &vec_err_true_xs.at(0));

		//Bayesian Unfolding ------------------------------------------------------------------
		//int n_iteration=4;
		int Nmax_iteration=51;
		vector<double> iter;
		vector<double> chi2;

		for (int k=1; k<Nmax_iteration; ++k) {
			int n_iteration=k;
			RooUnfoldBayes uf_inc (res_inc, data_inc_bkgfree, n_iteration); //inc
			RooUnfoldBayes uf_st_inc (res_st_inc, data_st_inc_bkgfree, n_iteration); //st_inc
			RooUnfoldBayes uf2d_inc (res2d_inc, data2d_inc_bkgfree, n_iteration); //inc_2d
			RooUnfoldBayes uf_int (res_int, data_int_bkgfree, n_iteration); //int

			//unfolding
			TH1D *data_inc_uf;
			TH1D *data_st_inc_uf;

			TH2D *data2d_inc_uf; //2d
			TH1D *data1d_inc_uf; //2d->1d
			TH1D *data1d_st_inc_uf; //2d->1d

			TH1D *data_int_uf;

			data_inc_uf=(TH1D* )uf_inc.Hreco();
			data_st_inc_uf=(TH1D* )uf_st_inc.Hreco();
			data2d_inc_uf=(TH2D*)uf2d_inc.Hreco();
			data1d_inc_uf=(TH1D*)data2d_inc_uf->ProjectionY();
			data1d_st_inc_uf=(TH1D*)data2d_inc_uf->ProjectionX();
			data_int_uf=(TH1D *)uf_int.Hreco();

			//chi^2 calculation ------------------------
			double chi2_each=0;

			double reco_inc[nthinslices] = {0};
			double reco_int[nthinslices] = {0};
			double err_reco_inc[nthinslices] = {0};
			double err_reco_int[nthinslices] = {0};
			double reco_xs[nthinslices] = {0};
			double err_reco_xs[nthinslices] = {0};

			vector<double> vec_reco_xs;
			vector<double> vec_err_reco_xs;

			for (int i = 0; i<nthinslices; ++i) {
				//RECO -------------------------------------------------------------------------------------------
				reco_int[i]=data_int_uf->GetBinContent(i+2);
				err_reco_int[i]=data_int_uf->GetBinError(i+2);

				//New way to calculate INC -------------------------------------------------
				for (int j=0; j<=i; ++j) {
					//2d unfolding
					reco_inc[i]+=data1d_st_inc_uf->GetBinContent(j+2); //2d->1d
					err_reco_inc[i]+=pow(data1d_st_inc_uf->GetBinError(j+2),2); //2d->1d
				}

				for (int j=0; j<=i-1; ++j){
					//2d unfolding
					reco_inc[i]-=data1d_inc_uf->GetBinContent(j+2);
					err_reco_inc[i]+=pow(data1d_inc_uf->GetBinError(j+2),2);
				}
				//--------------------------------------------------------------------------

				//do NOT use the old way to calculate INC, bias at high KE! --------------------------------------------------------------------------------------------------
				//std::cout<<"sliceid["<<i<<"]="<<sliceid[i]<<" reco_inc["<<i<<"]="<<reco_inc[i]<<"KE["<<i<<"]="<<KE[i]<<" KE_selection["<<i<<"]="<<KE_selection[i]<<std::endl;
				err_reco_inc[i] = sqrt(err_reco_inc[i]);

				//reco xs
				reco_xs[i]=xs_const*dEdx[i]*log(reco_inc[i]/(reco_inc[i]-reco_int[i]));
				err_reco_xs[i]=xs_const*dEdx[i]*sqrt(reco_int[i]+pow(reco_int[i],2)/reco_inc[i])/reco_inc[i];

				//for chi2 calc -----------------------------------
				if (KE_selection[i]==1) {
					vec_reco_xs.push_back(reco_xs[i]);
					vec_err_reco_xs.push_back(err_reco_xs[i]);
				}

			}

			chi2_each=neyman_chi2_data_mc(vec_reco_xs, vec_err_reco_xs, vec_true_xs, vec_err_true_xs);
			iter.push_back(k);
			chi2.push_back(chi2_each);

			std::cout<<"chi2:"<<chi2_each<<std::endl;

			TGraphErrors *gr_recoxs = new TGraphErrors(vec_ke.size(), &vec_ke.at(0), &vec_reco_xs.at(0), &vec_zero.at(0), &vec_err_reco_xs.at(0));
			gr_recoxs->SetNameTitle("gr_recoxs", "; Energy (MeV); Cross-section [mb]");

			//show results ------------------------------------------------------
			TCanvas *c_all = new TCanvas("c_all", "c_all", 900, 600);
			c_all->Divide(1,1);
			c_all->cd(1);
			TH2D *f2d=new TH2D("f2d",Form("N(iteration)=%d ; Energy (MeV); Cross-section [mb]",k), 420,0,420,900,300,1200);
			f2d->Draw();
			gr_truexs->SetNameTitle(Form("N(iteration)=%d ; Energy (MeV); Cross-section [mb]",k));
			gr_truexs->SetMarkerColor(2);
			gr_truexs->SetLineColor(2);
			gr_recoxs->SetLineColor(1);
			gr_recoxs->SetMarkerColor(1);
			gr_recoxs->SetMarkerStyle(24);

			gr_truexs->Draw("p same");
			gr_recoxs->Draw("p same");

			TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
			leg->SetFillStyle(0);
			leg->AddEntry(gr_truexs, "Truth", "ep");
			leg->AddEntry(gr_recoxs, "Reco", "ep");
			leg->AddEntry((TObject*)0, Form("#chi^{2}/NDF=%.2f/%d",chi2_each,NDF), "");
			leg->Draw();

			TString out_pdf;
			if (k==1) out_pdf=Form("xs_50.pdf(");
			if (k>1&&k<Nmax_iteration-2) out_pdf=Form("xs_50.pdf(");
			if (k==Nmax_iteration-1) out_pdf=Form("xs_50.pdf)");

			c_all->Print(outpath+out_pdf);
		}


		TGraph *gr_iter_chi2 = new TGraph(iter.size(), &iter.at(0), &chi2.at(0));
		gr_iter_chi2->SetName("gr_iter_chi2");

		TFile *f_out_int = new TFile(outpath+Form("iter_chi2_50.root"),"RECREATE");
			gr_iter_chi2->Write();
		f_out_int->Close();




		/*
		//show results ------------------------------------------------------
		TCanvas *c_all = new TCanvas("c_all", "c_all", 900, 1400);
		c_all->Divide(1,3);
		c_all->cd(1);
		mc_truesliceID_inel->SetLineColor(2);
		mc_truesliceID_inel->SetMarkerColor(2);
		data_int_uf->SetMarkerStyle(20);
		data_int_uf->Draw("ep ");
		mc_truesliceID_inel->Draw("hist same");

		c_all->cd(2);
		mc_truesliceID_all->SetLineColor(2);
		mc_truesliceID_all->SetMarkerColor(2);
		data1d_inc_uf->SetMarkerStyle(20);
		data1d_inc_uf->Draw("ep");
		mc_truesliceID_all->Draw("hist same");

		c_all->cd(3);
		mc_true_st_sliceID_all->SetLineColor(2);
		mc_true_st_sliceID_all->SetMarkerColor(2);
		data1d_st_inc_uf->SetMarkerStyle(20);
		data1d_st_inc_uf->Draw("ep");
		mc_true_st_sliceID_all->Draw("hist same");


		cout<<"mc_truesliceID_inel:"<<mc_truesliceID_inel->Integral()<<endl;
		cout<<"mc_true_st_sliceID_all:"<<mc_true_st_sliceID_all->Integral()<<endl;
		cout<<"mc_truesliceID_all:"<<mc_truesliceID_all->Integral()<<endl;

*/



		/*
		   data1d_inc_uf->SetLineColor(1);
		   mc_truesliceID_all->SetLineColor(3);
		   mc_truesliceID_all->SetMarkerColor(3);
		   mc_truesliceID_all->Draw("hist");
		   data1d_inc_uf->Draw("hist same");

*/


		//TH1D *mc_truesliceID_inel=(TH1D *)f_mc->Get("h_truesliceid_inelastic_all"); //IsPureInEL
		//TH1D *mc_truesliceID_all=(TH1D *)f_mc->Get("h_truesliceid_all"); //all protons
		//TH1D *mc_true_st_sliceID_all=(TH1D *)f_mc->Get("h_true_st_sliceid_all"); //all protons











		//c_inc->cd(1)->SetLogy();
		//TH2D *f2d_inc=new TH2D("f2d_inc",Form(""), 31, 10, 41, 600, 0, 60000);
		//TH2D *f2d_inc=new TH2D("f2d_inc",Form(""), 42, -1, 41, 100, 0, 10000);
		//f2d_inc->SetTitle("Incident Histogram; SliceID; Counts");
		//f2d_inc->Draw();
		//data_inc_bkgfree->Draw("ep same");
		//mc_truesliceID_all->Draw("hist same")






		/*
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

*/



	}
