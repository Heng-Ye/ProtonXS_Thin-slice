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
#include <algorithm>

#include "../headers/BasicParameters.h"
#include "../headers/BasicFunctions.h"

void make_SYSdataXS() {
	//read xs files -----------------------------------------------------------------------//
        //TString str_cen="./xs_files/xs_Eslice_dE20MeV_40slcs_bmrw_data_cen.root"; //xs_central

	//Systematic uncertainty due to the shift of KE at TPC FF (shift is due to error propagation of fitting)
        //TString str_up="./xs_files/xs_Eslice_dE20MeV_40slcs_bmrwup_data.root"; //xs_up
        //TString str_dn="./xs_files/xs_Eslice_dE20MeV_40slcs_bmrwdn_data.root"; //xs_dn
        //TString str_output="./xs_files/xs_Eslice_dE20MeV_40slcs_bmrw_data_cen_sysbmrw.root"; //out

	//These upstream E-loss is over-counting the contribution, neglect the contribution
        //TString str_up="./xs_files/xs_Eslice_dE20MeV_40slcs_elossup_data.root"; //xs_up
        //TString str_dn="./xs_files/xs_Eslice_dE20MeV_40slcs_elossdn_data.root"; //xs_dn
        //TString str_output="./xs_files/xs_Eslice_dE20MeV_40slcs_bmrw_data_cen_sysbeloss.root"; //out

        //TString str_up="./xs_files/xs_Eslice_dE20MeV_40slcs_bmrw_data_cen.root"; //xs_up
        //TString str_dn="./xs_files/xs_Eslice_dE20MeV_40slcs_BB_data.root"; //xs_dn
        //TString str_output="./xs_files/xs_Eslice_dE20MeV_40slcs_bmrw_data_cen_sysBB.root"; //out

	//New sliceIDs:(st=ceil, end=floor) ------------------------------------------------------------------------------------------------------------
        TString str_cen="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_data_cen_newslcid.root"; //xs_central
	
	//BB
        //TString str_up="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_data_cen_newslcid.root"; //xs_up
        //TString str_dn="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_data_BBI_newslcid.root"; //xs_dn
        //TString str_output="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_SYS_BB_newslcid.root"; //out

	//BMRW-fit
        //TString str_up="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_elossup_data_newslcid.root"; //xs_up
        //TString str_dn="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_elossdn_data_newslcid.root"; //xs_dn
        //TString str_output="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_SYS_BMRW_newslcid.root"; //out

	//El-scale
        //TString str_up="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_data_ElbkgScalingsubtraction_up_newslcid.root"; //xs_up
        //TString str_dn="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_data_ElbkgScalingsubtraction_dn_newslcid.root"; //xs_dn
        //TString str_output="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_SYS_ElScale_newslcid.root"; //out

	//MisIDP-scale
        //TString str_up="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_data_MisIDPbkgScalingsubtraction_up_newslcid.root"; //xs_up
        //TString str_dn="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_data_MisIDPbkgScalingsubtraction_dn_newslcid.root"; //xs_dn
        //TString str_output="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_SYS_MisIDPScale_newslcid.root"; //out

	//Beamline inst.
        TString str_up="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_data_beamline_up_newslcid.root"; //xs_up
        TString str_dn="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_data_beamline_dn_newslcid.root"; //xs_dn
        TString str_output="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_SYS_beamline_newslcid.root"; //out


        //plot style --------------------------------------------------------------------//
        gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
        gROOT->SetStyle("protoDUNEStyle");
        gROOT->ForceStyle();
        gStyle->SetTitleX(0.5);
        gStyle->SetTitleAlign(23);
        gStyle->SetOptStat(0);

	//read data ---------------------------------------------------------------//
	TFile *f_cen = TFile::Open(str_cen.Data()); //xs_central
		TGraphErrors* tr_xs_cen=(TGraphErrors* )f_cen->Get("gr_recoxs");
	TFile *f_up = TFile::Open(str_up.Data()); //xs_up
		TGraphErrors* tr_xs_up=(TGraphErrors* )f_up->Get("gr_recoxs");
	TFile *f_dn = TFile::Open(str_dn.Data()); //xs_dn
		TGraphErrors* tr_xs_dn=(TGraphErrors* )f_dn->Get("gr_recoxs");

	Double_t *ke=tr_xs_cen->GetX();
	Double_t *cen_xs=tr_xs_cen->GetY();
	Double_t *err_ke=tr_xs_cen->GetEX();
	Double_t *err_cen_xs=tr_xs_cen->GetEY();

	Double_t *up_xs=tr_xs_up->GetY();
	Double_t *dn_xs=tr_xs_dn->GetY();

	vector<double> KE;
	vector<double> err_KE;
	vector<double> XS;
	vector<double> err_XS_CEN;
	vector<double> err_XS_SYS;
	vector<double> err_XS_ALL;

	vector<double> err_SYS_KE;
	vector<double> err_ALL_KE;

	int n=tr_xs_cen->GetN();
	for (int i=0; i<n; i++) { 
	  	//printf("%g %g %g %g\n",ke[i],cen_xs[i],err_ke[i], err_cen_xs[i]);

		KE.push_back(ke[i]);
		err_KE.push_back(err_ke[i]);

		XS.push_back(cen_xs[i]);

		double up_minus_cen_xs=up_xs[i]-cen_xs[i];
		double cen_minus_dn_xs=cen_xs[i]-dn_xs[i];

		//find max among 3 values -------------------
		double up_max_xs=max(up_xs[i], cen_xs[i]);
		up_max_xs=max(up_max_xs, dn_xs[i]);

		//find min among 3 values -------------------
		double dn_min_xs=min(up_xs[i], cen_xs[i]);
		dn_min_xs=min(dn_min_xs, dn_xs[i]);
	
		double err_sys=(up_max_xs-dn_min_xs)/2.;

		cout<<ke[i]<<" cen_xs:"<<cen_xs[i]<<" err_sys:"<<err_sys<<" | up_max_xs:"<<up_max_xs<<" | dn_min_xs:"<<dn_min_xs<<" :: "<<up_xs[i]<<" "<<cen_xs[i]<<" "<<dn_xs[i]<<endl;

		//combine -----------------------------------------------------//
		err_XS_CEN.push_back(err_cen_xs[i]);
		err_XS_SYS.push_back(err_sys);

		double err_all=sqrt(pow(err_sys,2)+pow(err_cen_xs[i],2));
		err_XS_ALL.push_back(err_all);

		//Uncertainties due to fiber position shift --------------------------------------------------------------------------//
		double err_pbeam_fiber_pos_shift=1.3/100.; //in percentage
		double dke_up=1000.*p2ke(Pbeam_sys(ke2p(ke[i]/1000.), err_pbeam_fiber_pos_shift, 1))-ke[i]; //up value
		double dke_dn=ke[i]-1000.*p2ke(Pbeam_sys(ke2p(ke[i]/1000.), err_pbeam_fiber_pos_shift, -1)); //dn value
		double dke=dke_up; if (dke_up<0) dke_up=-1.*dke_up;

		err_SYS_KE.push_back(dke);
		err_ALL_KE.push_back(sqrt(pow(dke,2)+pow(err_ke[i],2)));
	}

	//output the combined xs results -------------------------------------------------------------------------------------------------------------------------------------------//
	TGraphErrors* reco_xs_sysonly=new TGraphErrors(n, &KE.at(0), &XS.at(0), &err_SYS_KE.at(0), &err_XS_SYS.at(0));	
	TGraphErrors* reco_xs_all=new TGraphErrors(n, &KE.at(0), &XS.at(0), &err_ALL_KE.at(0), &err_XS_ALL.at(0));
	//auto gr = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);	

	TFile *fout=new TFile(Form("%s",str_output.Data()),"recreate");	
		reco_xs_sysonly->Write("reco_xs_sysonly");
		reco_xs_all->Write("reco_xs_all");
		tr_xs_cen->Write("reco_xs");
	fout->Close();


}
