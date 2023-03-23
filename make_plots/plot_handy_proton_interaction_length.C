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

void plot_handy_proton_interaction_length() {


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

		//plot style --------------------------------------------------------------------//
		gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
		gROOT->SetStyle("protoDUNEStyle");
		gROOT->ForceStyle();
		gStyle->SetTitleX(0.5);
		gStyle->SetTitleAlign(23);
		gStyle->SetOptStat(0);

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
		std::cout<<"Density="<<Density<<std::endl;





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
		vector<double> vec_inel_len;
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
				vec_inel_len.push_back(MAr/(true_xs[i]*NA*Density*(1e-27))); //unit:cm
			}
		}	

		TGraph *gr_inellen_ke = new TGraph(vec_ke.size(), &vec_ke.at(0), &vec_inel_len.at(0));


			//show results ------------------------------------------------------
			TCanvas *c_all = new TCanvas("c_all", "c_all", 900, 600);
			c_all->Divide(1,1);
			c_all->cd(1);
			TH2D *f2d=new TH2D("f2d",Form("Proton Inelastic-scattering ; Energy (MeV); Interaction Length [cm]"), 420,0,420,10,0,100);
			f2d->Draw();
			//gr_recoxs->SetMarkerColor(1);
			//gr_recoxs->SetMarkerStyle(24);

			gr_inellen_ke->Draw("p same");


}
