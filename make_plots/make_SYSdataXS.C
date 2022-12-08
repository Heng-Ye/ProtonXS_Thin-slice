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

void make_SYSdataXS() {
	//read xs files -----------------------------------------------------------------------//
        TString str_cen="./xs_files/xs_Eslice_dE20MeV_40slcs_bmrw_data_cen.root"; //xs_central
        TString str_up="./xs_files/xs_Eslice_dE20MeV_40slcs_bmrwup_data.root"; //xs_up
        TString str_dn="./xs_files/xs_Eslice_dE20MeV_40slcs_bmrwdn_data.root"; //xs_dn
        TString str_output="./xs_files/xs_Eslice_dE20MeV_40slcs_bmrw_data_cen_sysbmrw.root"; //xs_dn

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
	vector<double> err_XS_SYS_UP;
	vector<double> err_XS_SYS_DN;
	vector<double> err_XS_ALL_UP;
	vector<double> err_XS_ALL_DN;

	int n=tr_xs_cen->GetN();
	for (int i=0; i<n; i++) { 
	  	//printf("%g %g %g %g\n",ke[i],cen_xs[i],err_ke[i], err_cen_xs[i]);

		KE.push_back(ke[i]);
		err_KE.push_back(err_ke[i]);

		XS.push_back(cen_xs[i]);

		double up_minus_cen_xs=up_xs[i]-cen_xs[i];
		double cen_minus_dn_xs=cen_xs[i]-dn_xs[i];

		//up_xs case -----------------------------//
		double err_sys_up1=0;
		double err_sys_dn1=0;
		if (up_minus_cen_xs>0) { //up_xs>cen_xs
			err_sys_up1=up_minus_cen_xs;
		} //up_xs>cen_xs
		else if (up_minus_cen_xs==0) {
			//do nothing
		} 
		else {
			err_sys_dn1=-up_minus_cen_xs;
		}

		//dn xs case -----------------------------//
		double err_sys_dn2=0;
		double err_sys_up2=0;
		if (cen_minus_dn_xs>0) { //cen_xs>dn_xs
			err_sys_dn2=cen_minus_dn_xs;
		} //cen_xs>dn_xs
		else if (cen_minus_dn_xs==0) {
		}
		else {
			err_sys_up2=-cen_minus_dn_xs;
		}

		//combine -----------------------------------------------------//
		double err_sys_up=sqrt(pow(err_sys_up1,2)+pow(err_sys_up2,2));
		double err_sys_dn=sqrt(pow(err_sys_dn1,2)+pow(err_sys_dn2,2));

		err_XS_CEN.push_back(err_cen_xs[i]);
		err_XS_SYS_UP.push_back(err_sys_up);
		err_XS_SYS_DN.push_back(err_sys_dn);

		double err_all_up=sqrt(pow(err_sys_up,2)+pow(err_cen_xs[i],2));
		double err_all_dn=sqrt(pow(err_sys_dn,2)+pow(err_cen_xs[i],2));
		err_XS_ALL_UP.push_back(err_all_up);
		err_XS_ALL_DN.push_back(err_all_dn);

	}

	//output the combined xs results ------------------------------------------------------------------------------------------------------------------------------//
	TGraphAsymmErrors* reco_xs_sysonly=new TGraphAsymmErrors(n, &KE.at(0), &XS.at(0), &err_KE.at(0), &err_KE.at(0), &err_XS_SYS_DN.at(0), &err_XS_SYS_UP.at(0));	
	TGraphAsymmErrors* reco_xs_all=new TGraphAsymmErrors(n, &KE.at(0), &XS.at(0), &err_KE.at(0), &err_KE.at(0), &err_XS_ALL_DN.at(0), &err_XS_ALL_UP.at(0));	

	TFile *fout=new TFile(Form("%s",str_output.Data()),"recreate");	
		reco_xs_sysonly->Write("reco_xs_sysonly");
		reco_xs_all->Write("reco_xs_all");
		tr_xs_cen->Write("reco_xs");
	fout->Close();


}
