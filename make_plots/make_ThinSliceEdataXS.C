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

//#include "../headers/SliceParams.h"
#include "../headers/BetheBloch.h"
#include "../headers/ESliceParams.h"

//constants ---------------------//
double NA=6.02214076e23;
double MAr=39.95; //gmol
double Density = 1.4; // g/cm^3
double true_cosangle = 1.;
//-------------------------------//

using namespace std;
using namespace ROOT::Math;


R__LOAD_LIBRARY(libRooUnfold.so) //load share lib
//#include "../headers/Unfold.h"
//#include "../headers/RooUnfold.h"
//#include "../headers/RooUnfold.h"
#include "../headers/RooUnfoldBayes.h"
#include "../headers/RooUnfoldResponse.h"

void make_ThinSliceEdataXS() {

	//TString outpath="./plots_XS_ThinsliceE/";
	//TString outpath="./plots_XS_ThinsliceE_StSlicdIDCeil/";
	//TString outpath="./plots_XS_ThinsliceE_StSlicdID+0.5/";
	//
	//TString outpath="./plots_DataXS_ThinsliceE_StSlicdID+0.5/";
	//TString outpath="./plots_XS_ThinsliceE_StSlicdIDp0.5_nobmrw/";

	TString outpath="./plots_XS_Eslice/";

	//TString outpath="./plots_XS_ThinsliceE_StSlicdID++/";
	//TString outpath="./plots_XS_ThinsliceE_plus1/";
	//TString outpath="./plots_XS_ThinsliceE_keff30/";

	//data & mc ----------------------------------------------------------------//
        //TString fmc="../prod4areco2_mc_ThinSliceE_dE30MeV_20slcs_nobmrw.root";
        //TString fdata="../prod4areco2_mc_ThinSliceE_dE30MeV_20slcs_nobmrw.root";

        //TString fmc="../prod4areco2_mc_ThinSliceE_dE15MeV_40slcs_nobmrw.root";
        //TString fdata="../prod4areco2_mc_ThinSliceE_dE15MeV_40slcs_nobmrw.root";

        //TString fmc="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw.root";
        //TString fdata="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw.root";

        //TString fmc="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_StIDplus1.root";
        //TString fdata="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_StIDplus1.root";

        //TString fmc="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_StIDminus1.root";
        //TString fdata="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_StIDminus1.root";

        //TString fmc="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_KEffZ0.root";
        //TString fdata="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_KEffZ0.root";

        //TString fmc="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_KEffZ0_stIDplus1.root";
        //TString fdata="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_KEffZ0_stIDplus1.root";

        //TString fmc="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_KEff30_old.root";
        //TString fdata="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_KEff30_old.root";

        //TString fmc="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw.root";
        //TString fdata="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw.root";

        //TString fmc="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_KEffZ0_stIDplus1.root";
        //TString fdata="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_KEffZ0_stIDplus1.root";

        //TString fmc="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_skip1stslc.root";
        //TString fdata="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_skip1stslc.root";

        //TString fmc="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_stslcplus1.root";
        //TString fdata="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_stslcplus1.root";

        //TString fmc="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_stslcplus0.5.root";
        //TString fdata="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_stslcplus0.5.root";

        //TString fmc="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_bmrw_stslcplus0.5.root";
        //TString fdata="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_bmrw_stslcplus0.5.root";

        //TString fmc="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_bmrw_stslcceil.root";
        //TString fdata="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_bmrw_stslcceil.root";

        //TString fmc="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_stslcceil.root";
        //TString fdata="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_stslcceil.root";


        //TString fmc="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_bmrw_stslcplus0.5.root";
        //TString fdata="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_run5387_prod4a_eslice_dx20cm_40slcs_stid+0.5.root";

        //TString fmc="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_stslcplus0.5.root";
        //TString fdata="../prod4areco2_mc_ThinSliceE_dE20MeV_40slcs_nobmrw_stslcplus0.5.root";

        TString fmc="../prod4areco2_mc_ESliceE_dE20MeV_40slcs_beamxy_bmrw_v09_39_01.root";
        TString fdata="/dune/data2/users/hyliao/protonana/v09_39_01/XS/prod4a_Eslice_dE20MeV_40slcs_beamxy_runAll_v09_39_01.root";

	//reco string pre-fix --------------------------------------//
	TString str_inc=Form("h_recosliceid_allevts_cuts");
	TString str_st_inc=Form("h_reco_st_sliceid_allevts_cuts");
	TString str_int=Form("h_recosliceid_recoinelastic_cuts");

        //plot style --------------------------------------------------------------------//
        gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
        gROOT->SetStyle("protoDUNEStyle");
        gROOT->ForceStyle();
        gStyle->SetTitleX(0.5);
        gStyle->SetTitleAlign(23);
        gStyle->SetOptStat(0);

	//read data ----------------------------------------------------------------------------------------------------------------//
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
	TH1D *mc_true_incidents=(TH1D *)f_mc->Get("h_true_incidents");
	TH1D *mc_true_st_incidents=(TH1D *)f_mc->Get("h_true_st_incidents");
	TH1D *mc_true_interactions=(TH1D *)f_mc->Get("h_true_interactions");

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


	THStack* hs_inc=new THStack("hs_inc","");
	hs_inc->Add(mc_inc_inel);
	hs_inc->Add(mc_inc_el);
	hs_inc->Add(mc_inc_midcosmic);
	hs_inc->Add(mc_inc_midpi);
	hs_inc->Add(mc_inc_midp);
	hs_inc->Add(mc_inc_midmu);
	hs_inc->Add(mc_inc_mideg);
	hs_inc->Add(mc_inc_midother);

	THStack* hs_st_inc=new THStack("hs_st_inc","");
	hs_st_inc->Add(mc_st_inc_inel);
	hs_st_inc->Add(mc_st_inc_el);
	hs_st_inc->Add(mc_st_inc_midcosmic);
	hs_st_inc->Add(mc_st_inc_midpi);
	hs_st_inc->Add(mc_st_inc_midp);
	hs_st_inc->Add(mc_st_inc_midmu);
	hs_st_inc->Add(mc_st_inc_mideg);
	hs_st_inc->Add(mc_st_inc_midother);

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
	//Note: Numerical value of errorbar after scaling is correct, i.e. s*sqrt(n)	

	THStack* hs_int=new THStack("hs_int","");
	hs_int->Add(mc_int_inel);
	hs_int->Add(mc_int_el);
	hs_int->Add(mc_int_midcosmic);
	hs_int->Add(mc_int_midpi);
	hs_int->Add(mc_int_midp);
	hs_int->Add(mc_int_midmu);
	hs_int->Add(mc_int_mideg);
	hs_int->Add(mc_int_midother);

	//mc, nobmrw
	//TFile *f_mc_nobmrw = TFile::Open(fmc_nobmrw.Data());
	//TH1D* mc_inc_all_nobmrw=(TH1D*)f_mc_nobmrw->Get(Form("%s",str_inc.Data()));
	//TH1D* mc_int_all_nobmrw=(TH1D*)f_mc_nobmrw->Get(Form("%s",str_int.Data()));
	//mc_inc_all_nobmrw->SetLineColor(9);
	//mc_int_all_nobmrw->SetLineColor(9);

	//mc_inc_all_nobmrw->Scale((double)n_data_inc/(double)mc_inc_all_nobmrw->Integral());
	//mc_int_all_nobmrw->Scale((double)n_data_int/(double)mc_int_all_nobmrw->Integral());
	

	//purity & efficiency 
	TH1D* pur_inc=(TH1D* )f_mc->Get("pur_Inc"); //inc
	TH1D *eff_inc=(TH1D* )f_mc->Get("eff_Inc"); //inc

	TH1D* pur_st_inc=(TH1D* )f_mc->Get("pur_st_Inc"); //inc
	TH1D *eff_st_inc=(TH1D* )f_mc->Get("eff_st_Inc"); //inc

	TH1D* pur_int=(TH1D* )f_mc->Get("pur_Int"); //int
	TH1D *eff_int=(TH1D* )f_mc->Get("eff_Int"); //int

	//bkg subtraction --------------------------------------------------------------------------//
	TH1D* data_inc_bkgfree=(TH1D *)data_inc->Clone("data_inc_bkgfree"); data_inc_bkgfree->SetName("data_inc_bkgfree");	
	TH1D* data_st_inc_bkgfree=(TH1D *)data_st_inc->Clone("data_st_inc_bkgfree"); data_st_inc_bkgfree->SetName("data_st_inc_bkgfree");	
	TH1D* data_int_bkgfree=(TH1D *)data_int->Clone("data_int_bkgfree"); data_int_bkgfree->SetName("data_int_bkgfree");	

	//inc
	//data_inc_bkgfree->Add(mc_inc_el, -1);
	//data_inc_bkgfree->Add(mc_inc_midcosmic, -1);
	//data_inc_bkgfree->Add(mc_inc_midpi, -1);
	data_inc_bkgfree->Add(mc_inc_midp, -1);
	//data_inc_bkgfree->Add(mc_inc_midmu, -1);
	//data_inc_bkgfree->Add(mc_inc_mideg, -1);
	//data_inc_bkgfree->Add(mc_inc_midother, -1);
	//data_inc_bkgfree->Multiply(pur_inc);
	//
	data_st_inc_bkgfree->Add(mc_st_inc_midp, -1);

	//
	//Note: Numerical value of errorbar after subtraction is correct, i.e. (s-b)+-sqrt(s+b)	
	
	//int
	data_int_bkgfree->Add(mc_int_el, -1);
	//data_int_bkgfree->Add(mc_int_midcosmic, -1);
	//data_int_bkgfree->Add(mc_int_midpi, -1);
	data_int_bkgfree->Add(mc_int_midp, -1);
	//data_int_bkgfree->Add(mc_int_midmu, -1);
	//data_int_bkgfree->Add(mc_int_mideg, -1);
	//data_int_bkgfree->Add(mc_int_midother, -1);

	//unfolding ---------------------------------------------------------------------------------------------//
	//response matrix from mc
	//Response matrix as a 2D-histogram: (x,y)=(measured,truth)
	RooUnfoldResponse *res_inc=(RooUnfoldResponse*)f_mc->Get("response_SliceID_Inc"); res_inc->SetName("res_inc");
	RooUnfoldResponse *res_st_inc=(RooUnfoldResponse*)f_mc->Get("response_st_SliceID_Inc"); res_st_inc->SetName("res_st_inc");

	//Response using shape
	//TH1D* res_Inc_reco=(TH1D*)f_mc->Get("res_Inc_reco"); res_Inc_reco->SetName("res_Inc_reco");
	//TH1D* res_Inc_truth=(TH1D*)f_mc->Get("res_Inc_truth"); res_Inc_truth->SetName("res_Inc_truth");
	//RooUnfoldResponse *res_inc=new RooUnfoldResponse(res_Inc_reco, res_Inc_truth, "res_inc","res_inc");
	
	//sansity check on the response matrix
	std::cout<<"res_inc->GetDimensionMeasured():"<<res_inc->GetDimensionMeasured()<<std::endl;
	std::cout<<"res_inc->GetDimensionTruth():"<<res_inc->GetDimensionTruth()<<std::endl;
	std::cout<<"res_inc->GetNbinsMeasured():"<<res_inc->GetNbinsMeasured()<<std::endl;
	std::cout<<"res_inc->GetNbinsTruth():"<<res_inc->GetNbinsTruth()<<std::endl; 

	//Unfolding
	RooUnfoldBayes uf_inc (res_inc, data_inc_bkgfree, 4); //inc
	RooUnfoldBayes uf_st_inc (res_st_inc, data_st_inc_bkgfree, 4); //inc

	//std::cout<<"GetIterations:"<<uf_inc.GetIterations()<<std::endl;
	//std::cout<<"GetSmoothing:"<<uf_inc.GetSmoothing()<<std::endl;
	//uf_inc.Print();

	RooUnfoldResponse *res_int=(RooUnfoldResponse*)f_mc->Get("response_SliceID_Int"); res_int->SetName("res_int"); //int
	//TH1D* res_Int_reco=(TH1D*)f_mc->Get("res_Int_reco"); res_Int_reco->SetName("res_Int_reco");
	//TH1D* res_Int_truth=(TH1D*)f_mc->Get("res_Int_truth"); res_Int_truth->SetName("res_Int_truth");
	//RooUnfoldResponse res_int(res_Int_reco, res_Int_truth, res_int);
	//RooUnfoldResponse *res_int=new RooUnfoldResponse(res_Int_reco, res_Int_truth, "res_int","res_int");
	RooUnfoldBayes uf_int (res_int, data_int_bkgfree, 4);

	//std::cout<<"data_inc_bkgfree x-axis bin:"<<data_inc_bkgfree->GetNbinsX()<<std::endl;
	//std::cout<<"data_int_bkgfree x-axis bin:"<<data_int_bkgfree->GetNbinsX()<<std::endl;
	//std::cout<<"response_SliceID_Inc x-axis bin:"<<(TH2D *)f_mc->Get("response_SliceID_Int")->GetNbinsX()<<std::endl;
	
	//unfolding
  	TH1D *data_inc_uf;
  	TH1D *data_st_inc_uf;
  	TH1D *data_int_uf;
	data_inc_uf=(TH1D* )uf_inc.Hreco();
	data_st_inc_uf=(TH1D* )uf_st_inc.Hreco();
  	data_int_uf=(TH1D *)uf_int.Hreco();
	data_inc_uf->SetNameTitle("data_inc_uf", "Unfolded incident protons; Slice ID; Events");
	data_st_inc_uf->SetNameTitle("data_st_inc_uf", "Unfolded incident protons; Start Slice ID; Events");
	data_int_uf->SetNameTitle("data_int_uf", "Unfolded interacting protons; Slice ID; Events");
	//unfold.IncludeSystematics(2); // Default 1, propagates both statistical+systematics, =2 no errors included.


	//draw slideID to KE conversion map ----------------------------------------//
	TCanvas *c_map=new TCanvas(Form("c_map"),"",900, 600);
        TH2D *f2d_map=new TH2D("f2d_map","",nthinslices,0,nthinslices,500,Emin,Emax);
	f2d_map->SetTitle("; SliceID; Proton Kinetic Energy [MeV]");
	
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.02); 
	gStyle->SetPadLeftMargin(0.12); 

	c_map->Divide(1,1);
	c_map->cd(1);
	//c_map->SetGrid(1,1);
	f2d_map->Draw();

	TBox *box[nthinslices];
	TText *tt[nthinslices];	
        for (int i = 0; i<nthinslices; ++i){
	  double id=i+0.5;
	  double ke=Emax-((double)i+0.5)*thinslicewidth;
	  std::cout<<"("<<id<<","<<ke<<")"<<endl;	

	  box[i]=new TBox(id-0.5,ke-thinslicewidth*.5,id+0.5,ke+thinslicewidth*.5);
	  box[i]->SetLineColor(4);
	  box[i]->SetLineWidth(2);
	  box[i]->SetFillColor(5);
          box[i]->SetFillStyle(0);
	  box[i]->Draw("same");
	
	  tt[i]=new TText(id-0.3, ke-thinslicewidth*.3,Form("%d",i));
	  tt[i]->SetTextColor(kRed);
	  tt[i]->SetTextSize(0.03);
	  tt[i]->Draw("same");
	  //box[i]->SetFillColorAlpha(0,0.5);
	}	
	c_map->Print(Form("%sID2KE_map.eps",outpath.Data()));




	//get KEs from MC --------------------------------------------------------//
	double slcid[nthinslices] = {0};
	TH1D *reco_incE[nthinslices];
	TH1D *true_incE[nthinslices];
	double avg_recoincE[nthinslices] = {0};
	double avg_trueincE[nthinslices] = {0};
	double err_recoincE[nthinslices] = {0};
	double err_trueincE[nthinslices] = {0};
	double rms_trueincE[nthinslices] = {0};
	double rms_recoincE[nthinslices] = {0};

	double reco_trueincE[nthinslices] = {0};
	double err_reco_trueincE[nthinslices] = {0};

	cout<<"nthinslices:"<<nthinslices<<endl;
        for (int i = 0; i<nthinslices; ++i){
                slcid[i]=i+.5;

		reco_incE[i]=(TH1D* )f_mc->Get(Form("reco_incE_%d",i));
		true_incE[i]=(TH1D* )f_mc->Get(Form("true_incE_%d",i));

                avg_trueincE[i] = true_incE[i]->GetMean();
                avg_recoincE[i] = reco_incE[i]->GetMean();
                err_trueincE[i] = true_incE[i]->GetMeanError();
                err_recoincE[i] = reco_incE[i]->GetMeanError();
                rms_trueincE[i] = true_incE[i]->GetRMS();
                rms_recoincE[i] = reco_incE[i]->GetRMS();

		reco_trueincE[i] = avg_recoincE[i] - avg_trueincE[i];
		err_reco_trueincE[i] = sqrt(pow(err_trueincE[i],2)+pow(err_recoincE[i],2));

		//cout<<"avg_recoincE["<<i<<"]:"<<avg_recoincE[i]<<endl;
        }
	
	//Calc XS ----------------------------------------------------------------------------------------------------------------------------//
	double xs_const=MAr/(Density*NA*thinslicewidth)*1e27;
	cout<<"xs_const:"<<xs_const<<" [mb]"<<endl;

	//[0]KE estimation
	BetheBloch BB(2212);
  	double KE[nthinslices] = {0};
	double err_KE[nthinslices] = {0};
  	double dEdx[nthinslices] = {0};
        for (int i = 0; i<nthinslices; ++i) {
    		KE[i]=Emax-((double)i+0.5)*thinslicewidth; //av_KE
		//KE[i]=true_incE[i]->GetMean();
		err_KE[i]=(double)thinslicewidth/2.;
		dEdx[i]=BB.meandEdx(KE[i]); // MeV/cm
        }

	//[1] true xs
	double true_inc[nthinslices]={0};
	double true_int[nthinslices]={0};
	double err_true_inc[nthinslices]={0};
	double err_true_int[nthinslices]={0};
	double true_xs[nthinslices] = {0};
	double err_true_xs[nthinslices] = {0};
  	for (int i = 0; i<nthinslices; ++i){
		//[1a]true inc/int from truth sliceID dists.
    		true_int[i] = mc_truesliceID_inel->GetBinContent(i+2);
    		err_true_int[i] = mc_truesliceID_inel->GetBinError(i+2);

    		for (int j=i; j<=nthinslices; ++j){
      			true_inc[i]+=mc_truesliceID_all->GetBinContent(j+2);
      			err_true_inc[i]+=pow(mc_truesliceID_all->GetBinError(j+2),2);
    		}
    		for (int j=i+1; j<=nthinslices; ++j){
      			true_inc[i]-=mc_true_st_sliceID_all->GetBinContent(j+2);
      			err_true_inc[i]+=pow(mc_true_st_sliceID_all->GetBinError(j+2),2);
    		}
    		err_true_inc[i] = sqrt(err_true_inc[i]);

		//if (true_inc[i]&&true_int[i]) {
			//note that 1/dE has been included in the xs_constant
			true_xs[i]=xs_const*dEdx[i]*log(true_inc[i]/(true_inc[i]-true_int[i]));
			err_true_xs[i]=xs_const*dEdx[i]*sqrt(true_int[i]+pow(true_int[i],2)/true_inc[i])/true_inc[i];
		//}
  	}

	//[2] reco xs 
	//inc & int
	double sliceid[nthinslices] = {0};
	double reco_inc[nthinslices] = {0};
	double reco_int[nthinslices] = {0};
	double err_reco_inc[nthinslices] = {0};
	double err_reco_int[nthinslices] = {0};
	double err_sliceid[nthinslices] = {0};
	double zero[nthinslices] = {0};

        //for (int i = 0; i<nthinslices; ++i) {
		//cout<<"data_inc_uf["<<i+2<<"]:"<<data_inc_uf->GetBinContent(i+2)<<endl;
	//}
	
	double reco_xs[nthinslices] = {0};
	double err_reco_xs[nthinslices] = {0};
        for (int i = 0; i<nthinslices; ++i) {
		sliceid[i]=i+.5;

		//int[inelastic]
		reco_int[i]=data_int_uf->GetBinContent(i+2);
		err_reco_int[i]=data_int_uf->GetBinError(i+2);

		//inc[all]
    		for (int j=i; j<=nthinslices; ++j){
      			reco_inc[i]+=data_inc_uf->GetBinContent(j+2);
      			err_reco_inc[i]+=pow(data_inc_uf->GetBinError(j+2),2);
    		}
    		for (int j=i+1; j<=nthinslices; ++j){
      			reco_inc[i]-=data_st_inc_uf->GetBinContent(j+2);
      			err_reco_inc[i]+=pow(data_st_inc_uf->GetBinError(j+2),2);
    		}
    		err_reco_inc[i] = sqrt(err_reco_inc[i]);

		//reco xs
		reco_xs[i]=xs_const*dEdx[i]*log(reco_inc[i]/(reco_inc[i]-reco_int[i]));
		err_reco_xs[i]=xs_const*dEdx[i]*sqrt(reco_int[i]+pow(reco_int[i],2)/reco_inc[i])/reco_inc[i];
		cout<<"reco_xs["<<i<<"]="<<reco_xs[i]<<endl;
		
	}

	//for (int ij = 0; ij<=true_sliceID; ++ij){
		//if (ij<nthinslices) ++true_incidents[ij];
	//}

        //for (int i = 0; i<nthinslices; ++i) {
		//if (reco_inc[i]&&reco_int[i]) {
    			//reco_xs[i]=xs_const*dEdx[i]*log(reco_inc[i]/(reco_inc[i]-reco_int[i]));
    			//err_reco_xs[i]=xs_const*dEdx[i]*sqrt(pow(reco_int[i]*err_reco_inc[i]/reco_inc[i]/(reco_inc[i]-reco_int[i]),2)+pow(err_reco_int[i]/(reco_inc[i]-reco_int[i]),2));
		//}
        //}
  	TGraphErrors *gr_truexs = new TGraphErrors(nthinslices, &KE[0], &true_xs[0], &err_KE[0], &err_true_xs[0]);
  	TGraphErrors *gr_recoxs = new TGraphErrors(nthinslices, &KE[0], &reco_xs[0], &err_KE[0], &err_reco_xs[0]);

  	TGraphErrors *gr_reco_inc = new TGraphErrors(nthinslices, slcid, reco_inc, err_sliceid, err_reco_inc);
  	TGraphErrors *gr_reco_int = new TGraphErrors(nthinslices, slcid, reco_int, err_sliceid, err_reco_int);

  	TGraphErrors *gr_sliceid_recoincE = new TGraphErrors(nthinslices, slcid, avg_recoincE, err_sliceid, err_recoincE);
  	TGraphErrors *gr_sliceid_trueincE = new TGraphErrors(nthinslices, slcid, avg_trueincE, err_sliceid, err_trueincE);
	TGraphErrors *gr_sliceid_recotrueincE = new TGraphErrors(nthinslices, &(slcid[0]), &(reco_trueincE[0]), 0, &(err_reco_trueincE[0]));

  	TGraphErrors *gr_slcid_trueE = new TGraphErrors(nthinslices, &sliceid[0], &avg_trueincE[0], &zero[0], &rms_trueincE[0]);
  	TGraphErrors *gr_slcid_recoE = new TGraphErrors(nthinslices, &sliceid[0], &avg_recoincE[0], &zero[0], &rms_recoincE[0]);

	gr_reco_inc->SetNameTitle("gr_reco_inc", "; SliceID; Events");
	gr_reco_int->SetNameTitle("gr_reco_int", "; ; ");

  	gr_truexs->SetNameTitle("gr_truexs", "; Energy (MeV); Cross-section [mb]");
  	gr_recoxs->SetNameTitle("gr_recoxs", "; Energy (MeV); Cross-section [mb]");

	gr_sliceid_recoincE->SetNameTitle("gr_sliceid_recoincE; SliceID; Proton Kinetic Energy [MeV]");
	gr_sliceid_trueincE->SetNameTitle("gr_sliceid_trueincE; SliceID; Proton Kinetic Energy [MeV]");

	gr_sliceid_recotrueincE->SetNameTitle("gr_sliceid_recotrueincE; SliceID; Reco KE - True KE [MeV]");

	gr_slcid_trueE->SetNameTitle("gr_slcid_trueE; Slice ID; Proton Kinetic Energy [MeV]");
        gr_slcid_recoE->SetNameTitle("gr_slcid_recoE; Slice ID; Proton Kinetic Energy [MeV]");


	//All protons ---------------------------------------------------------------------//
	//[0]sliceID vs KE	
	TCanvas *c_slcid_ke=new TCanvas(Form("c_slcid_ke"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.02); 
	gStyle->SetPadLeftMargin(0.12); 

	c_slcid_ke->Divide(1,1);
	c_slcid_ke->cd(1);
	float ymax_ke=800;
	//TH2D *f2d_slcid_ke=new TH2D("f2d_slcid_ke",Form("%s",""),20,0,20,600,0,ymax_ke);
	TH2D *f2d_slcid_ke=new TH2D("f2d_slcid_ke",Form("%s",""),nthinslices,0,nthinslices,600,0,ymax_ke);
	f2d_slcid_ke->SetTitle(";SliceID; Proton Kinetic Energy [MeV] (mean #pm RMS)");
	//f2d_inc->GetYaxis()->SetTitleOffset(1.3);
	f2d_slcid_ke->Draw();
	gr_slcid_trueE->SetMarkerColor(2);
	gr_slcid_trueE->SetLineColor(2);
	gr_slcid_trueE->Draw("p same");
	gr_slcid_recoE->Draw("p same");

	TLegend *leg_slcid_ke = new TLegend(0.7,0.6,0.9,0.9);
	leg_slcid_ke->SetFillStyle(0);
	leg_slcid_ke->AddEntry(gr_slcid_trueE, "Truth", "ep");
	leg_slcid_ke->AddEntry(gr_slcid_recoE, "Reco","ep");
        leg_slcid_ke->Draw();
	c_slcid_ke->Print(Form("%sslc_vs_KE.eps",outpath.Data()));



	//[0]data with mc components (start slice) --------------------------------------------------------------------//
	TCanvas *c_st_inc=new TCanvas(Form("c_st_inc"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.02); 
	gStyle->SetPadLeftMargin(0.12); 

	c_st_inc->Divide(1,1);
	c_st_inc->cd(1);
	float ymax_st_evt=20000;
	TH2D *f2d_st_inc=new TH2D("f2d_st_inc",Form("%s","All Protons"),nthinslices+2,-1,nthinslices+1,100,0,ymax_st_evt);
	f2d_st_inc->SetTitle("All Protons; Reco Start SliceID ; Events");
	f2d_st_inc->GetYaxis()->SetTitleOffset(1.3);
	f2d_st_inc->Draw();
	hs_st_inc->Draw("hist same");
	data_st_inc->Draw("ep same");

	TLegend *leg_st_inc = new TLegend(0.2,0.6,0.8,0.9);
	leg_st_inc->SetFillStyle(0);
	leg_st_inc->AddEntry(data_inc, "Fake Data", "ep");
	leg_st_inc->AddEntry(mc_inc_inel, "Inel","f");
	leg_st_inc->AddEntry(mc_inc_el, "El","f");

	leg_st_inc->AddEntry(mc_inc_midcosmic,"misID:cosmic","f");
	leg_st_inc->AddEntry(mc_inc_midp, "misID:p","f");
	leg_st_inc->AddEntry(mc_inc_midpi, "misID:#pi","f");

	leg_st_inc->AddEntry(mc_inc_midmu,"misID:#mu","f");
	leg_st_inc->AddEntry(mc_inc_mideg, "misID:e/#gamma","f");
	leg_st_inc->AddEntry(mc_inc_midother, "misID:other","f");

	leg_st_inc->SetNColumns(3);
	leg_st_inc->Draw();
	c_st_inc->Print(Form("%sdata_reco_st_sliceid_inc.eps",outpath.Data()));



	//[1]data with mc components
	TCanvas *c_inc=new TCanvas(Form("c_inc"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.02); 
	gStyle->SetPadLeftMargin(0.12); 

	c_inc->Divide(1,1);
	c_inc->cd(1);
	float ymax_evt=6000; //mc
	//float ymax_evt=2500; //data, run5387
	//TH2D *f2d_inc=new TH2D("f2d_inc",Form("%s","All Protons"),22,-1,21,80,0,ymax_evt);
	TH2D *f2d_inc=new TH2D("f2d_inc",Form("%s","All Protons"),nthinslices+2,-1,nthinslices+1,80,0,ymax_evt);
	f2d_inc->SetTitle("All Protons; Reco SliceID (End point); Events");
	f2d_inc->GetYaxis()->SetTitleOffset(1.3);
	f2d_inc->Draw();
	hs_inc->Draw("hist same");
	data_inc->Draw("ep same");

	TLegend *leg_inc = new TLegend(0.2,0.6,0.8,0.9);
	leg_inc->SetFillStyle(0);
	leg_inc->AddEntry(data_inc, "Fake Data", "ep");
	leg_inc->AddEntry(mc_inc_inel, "Inel","f");
	leg_inc->AddEntry(mc_inc_el, "El","f");

	leg_inc->AddEntry(mc_inc_midcosmic,"misID:cosmic","f");
	leg_inc->AddEntry(mc_inc_midp, "misID:p","f");
	leg_inc->AddEntry(mc_inc_midpi, "misID:#pi","f");

	leg_inc->AddEntry(mc_inc_midmu,"misID:#mu","f");
	leg_inc->AddEntry(mc_inc_mideg, "misID:e/#gamma","f");
	leg_inc->AddEntry(mc_inc_midother, "misID:other","f");

	leg_inc->SetNColumns(3);
	leg_inc->Draw();
	c_inc->Print(Form("%sdata_recosliceid_inc.eps",outpath.Data()));

	//[2]bkg subtraction
	TCanvas *c_inc2=new TCanvas(Form("c_inc2"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.02); 
	gStyle->SetPadLeftMargin(0.12); 

	c_inc2->Divide(1,1);
	c_inc2->cd(1);
	//TH2D *f2d_inc2=new TH2D("f2d_inc2",Form("%s","INC"),22,-1,21,80,0,ymax_evt);
	TH2D *f2d_inc2=new TH2D("f2d_inc2",Form("%s","INC"),nthinslices+2,-1,nthinslices+1,80,0,ymax_evt);
	f2d_inc2->SetTitle("All protons; Reco SliceID (End point); Events");
	f2d_inc2->GetYaxis()->SetTitleOffset(1.3);
	f2d_inc2->Draw();

	THStack* hs_inc_sg=new THStack("hs_inc_sg","");
	hs_inc_sg->Add(mc_inc_inel);
	hs_inc_sg->Add(mc_inc_el);
	hs_inc_sg->Draw("hist same");

	data_inc->Draw("ep same");
	data_inc_bkgfree->SetMarkerColor(1);
	data_inc_bkgfree->SetLineColor(1);
	data_inc_bkgfree->SetLineWidth(3);
	data_inc_bkgfree->SetMarkerStyle(24);
	data_inc_bkgfree->Draw("ep same");

	TLegend *leg_inc2 = new TLegend(0.2,0.6,0.8,0.9);
	leg_inc2->SetFillStyle(0);
	leg_inc2->AddEntry(data_inc, "Selected", "ep");
	leg_inc2->AddEntry(data_inc_bkgfree, "Selected+BKG Subtraction","ep");
	leg_inc2->AddEntry(mc_inc_inel, "Selected MC Truth [Inel]","f");
	leg_inc2->AddEntry(mc_inc_el, "Selected MC Truth [El]","f");
	//leg_inc2->SetNColumns(3);
	leg_inc2->Draw();
	c_inc2->Print(Form("%sdata_recosliceid_bkgsub_inc.eps",outpath.Data()));


	//[2]bkg subtraction(start sliceID)
	TCanvas *c_st_inc2=new TCanvas(Form("c_st_inc2"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.02); 
	gStyle->SetPadLeftMargin(0.12); 

	c_st_inc2->Divide(1,1);
	c_st_inc2->cd(1);
	TH2D *f2d_st_inc2=new TH2D("f2d_st_inc2",Form("%s","INC_st"),nthinslices+2,-1,nthinslices+1,80,0,ymax_st_evt);
	f2d_st_inc2->SetTitle("All protons; Reco Start SliceID ; Events");
	f2d_st_inc2->GetYaxis()->SetTitleOffset(1.3);
	f2d_st_inc2->Draw();

	THStack* hs_st_inc_sg=new THStack("hs_st_inc_sg","");
	hs_st_inc_sg->Add(mc_st_inc_inel);
	hs_st_inc_sg->Add(mc_st_inc_el);
	hs_st_inc_sg->Draw("hist same");

	data_st_inc->Draw("ep same");
	data_st_inc_bkgfree->SetMarkerColor(1);
	data_st_inc_bkgfree->SetLineColor(1);
	data_st_inc_bkgfree->SetLineWidth(3);
	data_st_inc_bkgfree->SetMarkerStyle(24);
	data_st_inc_bkgfree->Draw("ep same");

	TLegend *leg_st_inc2 = new TLegend(0.2,0.6,0.8,0.9);
	leg_st_inc2->SetFillStyle(0);
	leg_st_inc2->AddEntry(data_st_inc, "Selected", "ep");
	leg_st_inc2->AddEntry(data_st_inc_bkgfree, "Selected+BKG Subtraction","ep");
	leg_st_inc2->AddEntry(mc_st_inc_inel, "Selected MC Truth [Inel]","f");
	leg_st_inc2->AddEntry(mc_st_inc_el, "Selected MC Truth [El]","f");
	//leg_st_inc2->SetNColumns(3);
	leg_st_inc2->Draw();
	c_st_inc2->Print(Form("%sdata_reco_st_sliceid_bkgsub_inc.eps",outpath.Data()));
	
	
	//[3]data unfolding (true sliceID)
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.02); 
	gStyle->SetPadLeftMargin(0.12); 
	TCanvas *c_inc3=new TCanvas(Form("c_inc3"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_inc3->Divide(1,1);
	c_inc3->cd(1);
	//TH2D *f2d_inc3=new TH2D("f2d_inc3",Form("%s","All protons"),22,-1,21,100,0,ymax_evt);
	TH2D *f2d_inc3=new TH2D("f2d_inc3",Form("%s","All protons"),nthinslices+2,-1,nthinslices+1,100,0,ymax_evt);
	f2d_inc3->SetTitle("All Protons; True SliceID (End point); Events");
	f2d_inc3->GetYaxis()->SetTitleOffset(1.3);
	f2d_inc3->Draw();
	data_inc_bkgfree->Draw("ep same");

	data_inc_uf->SetLineColor(4);
	data_inc_uf->SetMarkerColor(4);
        data_inc_uf->SetMarkerStyle(21);
	data_inc_uf->Draw("ep same");

	mc_true_incidents->SetLineColor(2);
	mc_true_incidents->Draw("hist same");

	TLegend *leg_inc3 = new TLegend(0.4,0.6,0.9,0.9);
	leg_inc3->SetFillStyle(0);
	leg_inc3->AddEntry(data_inc_bkgfree, "Selected+BKG Subtraction","ep");
	leg_inc3->AddEntry(data_inc_uf, "Unfolding","ep");
	leg_inc3->AddEntry(mc_true_incidents,"MC Truth","l");
	leg_inc3->Draw();
	c_inc3->Print(Form("%sdata_recosliceid_uf_inc.eps",outpath.Data()));

	//[3]data unfolding (true start sliceID)
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.02); 
	gStyle->SetPadLeftMargin(0.12); 
	TCanvas *c_st_inc3=new TCanvas(Form("c_st_inc3"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_st_inc3->Divide(1,1);
	c_st_inc3->cd(1);
	//TH2D *f2d_inc3=new TH2D("f2d_inc3",Form("%s","All protons"),22,-1,21,100,0,ymax_evt);
	TH2D *f2d_st_inc3=new TH2D("f2d_st_inc3",Form("%s","All protons"),nthinslices+2,-1,nthinslices+1,100,0,ymax_st_evt);
	f2d_st_inc3->SetTitle("All Protons; True Start SliceID; Events");
	f2d_st_inc3->GetYaxis()->SetTitleOffset(1.3);
	f2d_st_inc3->Draw();
	data_st_inc_bkgfree->Draw("ep same");

	data_st_inc_uf->SetLineColor(4);
	data_st_inc_uf->SetMarkerColor(4);
        data_st_inc_uf->SetMarkerStyle(21);
	data_st_inc_uf->Draw("ep same");

	mc_true_st_incidents->SetLineColor(2);
	mc_true_st_incidents->Draw("hist same");

	TLegend *leg_st_inc3 = new TLegend(0.4,0.6,0.9,0.9);
	leg_st_inc3->SetFillStyle(0);
	leg_st_inc3->AddEntry(data_st_inc_bkgfree, "Selected+BKG Subtraction","ep");
	leg_st_inc3->AddEntry(data_st_inc_uf, "Unfolding","ep");
	leg_st_inc3->AddEntry(mc_true_st_incidents,"MC Truth","l");
	leg_st_inc3->Draw();
	c_st_inc3->Print(Form("%sdata_reco_st_sliceid_uf_inc.eps",outpath.Data()));



	//[4] eff & purity
	TCanvas *c_inc4=new TCanvas(Form("c_inc4"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.1); 
	gStyle->SetPadLeftMargin(0.1); 

	c_inc4->Divide(1,1);
	c_inc4->cd(1);
	float ymax_eff_inc=1.2;
	//TH2D *f2d_inc4=new TH2D("f2d_inc4",Form("%s",""),22,-1,21,10,0,ymax_eff_inc);
	TH2D *f2d_inc4=new TH2D("f2d_inc4",Form("%s",""),nthinslices+2,-1,nthinslices+1,10,0,ymax_eff_inc);
	f2d_inc4->SetTitle("All protons; Reco SliceID; Efficiency");
	f2d_inc4->Draw("");
	eff_inc->Draw("same");

	double av_eff_inc=0;
	double err_av_eff_inc=0;
	double cnt_eff_inc=0;
	cout<<"\neff_inc->GetNbinsX():"<<eff_inc->GetNbinsX()<<endl;
	for (int k=0; k<eff_inc->GetNbinsX(); ++k) {
		if (k==0) cout<<"eff 1st bin, slice #:"<<eff_inc->GetBinCenter(k)<<"\n\n"<<endl;
		if (eff_inc->GetBinContent(k+1)) { 
			cnt_eff_inc++;
			av_eff_inc+=eff_inc->GetBinContent(k+1);
			err_av_eff_inc+=pow(eff_inc->GetBinError(k+1),2);
		}
	}
	av_eff_inc/=cnt_eff_inc;
	err_av_eff_inc=sqrt(err_av_eff_inc)/cnt_eff_inc;
	cout<<"av_eff_inc="<<av_eff_inc<<" +- "<<err_av_eff_inc<<endl;

	TLine *line_av_eff_inc=new TLine(-1,av_eff_inc,nthinslices+1,av_eff_inc);
	line_av_eff_inc->SetLineColor(7);
	line_av_eff_inc->SetLineWidth(2);
	line_av_eff_inc->SetLineStyle(2);
	line_av_eff_inc->Draw("same");

	c_inc4->Update();
   	//scale hint1 to the pad coordinates
   	float rightmax_inc = ymax_eff_inc; //display max in y-axis for pur. dist.
   	float scale_inc = gPad->GetUymax()/rightmax_inc;
   	pur_inc->Scale(scale_inc);
	pur_inc->SetMarkerColor(2);
	pur_inc->SetLineColor(2);
	pur_inc->Draw("same");

	double av_pur_inc=0;
	double err_av_pur_inc=0;
	double cnt_pur_inc=0;
	for (int k=0; k<pur_inc->GetNbinsX(); ++k) {
		av_pur_inc+=pur_inc->GetBinContent(k+1);
		err_av_pur_inc+=pow(pur_inc->GetBinError(k+1),2);
		if (pur_inc->GetBinContent(k+1)) cnt_pur_inc++;

	}
	av_pur_inc/=cnt_pur_inc;
	err_av_pur_inc=sqrt(err_av_pur_inc)/cnt_pur_inc;
	cout<<"av_pur_inc="<<av_pur_inc<<" +- "<<err_av_pur_inc<<endl;

 
   	//draw an axis on the right side
   	TGaxis *axis_inc = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
         gPad->GetUxmax(), gPad->GetUymax(),0,rightmax_inc,510,"+L");
   	axis_inc->SetLineColor(2);
   	axis_inc->SetLabelColor(2);
	axis_inc->SetTitle("Purity");
	axis_inc->SetTextColor(2);
   	axis_inc->Draw();

	TLine *line_av_pur_inc=new TLine(-1,av_pur_inc,nthinslices+1,av_pur_inc);
	line_av_pur_inc->SetLineColor(3);
	line_av_pur_inc->SetLineWidth(2);
	line_av_pur_inc->SetLineStyle(2);
	line_av_pur_inc->Draw("same");

	c_inc4->Print(Form("%sdata_eff_pur_inc.eps",outpath.Data()));


	//[4] eff & purity (start sliceID)
	TCanvas *c_st_inc4=new TCanvas(Form("c_st_inc4"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.1); 
	gStyle->SetPadLeftMargin(0.1); 

	c_st_inc4->Divide(1,1);
	c_st_inc4->cd(1);
	float ymax_eff_st_inc=1.2;
	TH2D *f2d_st_inc4=new TH2D("f2d_st_inc4",Form("%s",""),nthinslices+2,-1,nthinslices+1,10,0,ymax_eff_st_inc);
	f2d_st_inc4->SetTitle("All protons; Reco Start SliceID; Efficiency");
	f2d_st_inc4->Draw("");
	eff_st_inc->Draw("same");

	double av_eff_st_inc=0;
	double err_av_eff_st_inc=0;
	double cnt_eff_st_inc=0;
	cout<<"\neff_st_inc->GetNbinsX():"<<eff_st_inc->GetNbinsX()<<endl;
	for (int k=0; k<eff_st_inc->GetNbinsX(); ++k) {
		if (k==0) cout<<"eff 1st bin, slice #:"<<eff_st_inc->GetBinCenter(k)<<"\n\n"<<endl;
		if (eff_st_inc->GetBinContent(k+1)) { 
			cnt_eff_st_inc++;
			av_eff_st_inc+=eff_st_inc->GetBinContent(k+1);
			err_av_eff_st_inc+=pow(eff_st_inc->GetBinError(k+1),2);
		}
	}
	av_eff_st_inc/=cnt_eff_st_inc;
	err_av_eff_st_inc=sqrt(err_av_eff_st_inc)/cnt_eff_st_inc;
	cout<<"av_eff_st_inc="<<av_eff_st_inc<<" +- "<<err_av_eff_st_inc<<endl;

	TLine *line_av_eff_st_inc=new TLine(-1,av_eff_st_inc,nthinslices+1,av_eff_st_inc);
	line_av_eff_st_inc->SetLineColor(7);
	line_av_eff_st_inc->SetLineWidth(2);
	line_av_eff_st_inc->SetLineStyle(2);
	line_av_eff_st_inc->Draw("same");

	c_st_inc4->Update();
   	//scale hint1 to the pad coordinates
   	float rightmax_st_inc = ymax_eff_st_inc; //display max in y-axis for pur. dist.
   	float scale_st_inc = gPad->GetUymax()/rightmax_st_inc;
   	pur_st_inc->Scale(scale_st_inc);
	pur_st_inc->SetMarkerColor(2);
	pur_st_inc->SetLineColor(2);
	pur_st_inc->Draw("same");

	double av_pur_st_inc=0;
	double err_av_pur_st_inc=0;
	double cnt_pur_st_inc=0;
	for (int k=0; k<pur_st_inc->GetNbinsX(); ++k) {
		av_pur_st_inc+=pur_st_inc->GetBinContent(k+1);
		err_av_pur_st_inc+=pow(pur_st_inc->GetBinError(k+1),2);
		if (pur_st_inc->GetBinContent(k+1)) cnt_pur_st_inc++;

	}
	av_pur_st_inc/=cnt_pur_st_inc;
	err_av_pur_st_inc=sqrt(err_av_pur_st_inc)/cnt_pur_st_inc;
	cout<<"av_pur_st_inc="<<av_pur_st_inc<<" +- "<<err_av_pur_st_inc<<endl;

 
   	//draw an axis on the right side
   	TGaxis *axis_st_inc = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
         gPad->GetUxmax(), gPad->GetUymax(),0,rightmax_st_inc,510,"+L");
   	axis_st_inc->SetLineColor(2);
   	axis_st_inc->SetLabelColor(2);
	axis_st_inc->SetTitle("Purity");
	axis_st_inc->SetTextColor(2);
   	axis_st_inc->Draw();

	TLine *line_av_pur_st_inc=new TLine(-1,av_pur_st_inc,nthinslices+1,av_pur_st_inc);
	line_av_pur_st_inc->SetLineColor(3);
	line_av_pur_st_inc->SetLineWidth(2);
	line_av_pur_st_inc->SetLineStyle(2);
	line_av_pur_st_inc->Draw("same");

	c_st_inc4->Print(Form("%sdata_eff_pur_st_inc.eps",outpath.Data()));



	//[5]Ninc distribution
	TCanvas *c_inc5=new TCanvas(Form("c_inc5"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.02); 
	gStyle->SetPadLeftMargin(0.12); 
	c_inc5->Divide(1,1);
	c_inc5->cd(1);
	//TH2D *f2d_inc5=new TH2D("f2d_inc5",Form("%s",""),21,0,21,10,0,150000);
	TH2D *f2d_inc5=new TH2D("f2d_inc5",Form("%s",""),nthinslices+1,0,nthinslices+1,10,0,150000);
	f2d_inc5->SetTitle("All protons; True SliceID; Events");
	f2d_inc5->Draw();
	mc_true_incidents->SetLineColor(2);
	mc_true_incidents->Draw("same");
	gr_reco_inc->Draw("ep same");
	c_inc5->Print(Form("%sincident_inc.eps",outpath.Data()));

	
	//inelastic scattering protons ------------------------------------------------//
	//[1]data with mc components
	TCanvas *c_int=new TCanvas(Form("c_int"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.02); 
	gStyle->SetPadLeftMargin(0.12); 

	c_int->Divide(1,1);
	c_int->cd(1);
	//TH2D *f2d_int=new TH2D("f2d_int",Form("%s","Proton Inelastic Scatterings"),21,0,21,80,0,ymax_evt);
	//ymax_evt=500;
	ymax_evt=6000;
	TH2D *f2d_int=new TH2D("f2d_int",Form("%s","Proton Inelastic Scatterings"),nthinslices+1,0,nthinslices+1,80,0,ymax_evt);
	f2d_int->SetTitle("Proton Inelastic Scatterings; Reco SliceID; Events");
	f2d_int->GetYaxis()->SetTitleOffset(1.3);
	f2d_int->Draw();
	hs_int->Draw("hist same");
	//mc_int_all_nobmrw->Draw("hist same");
	data_int->Draw("ep same");

	TLegend *leg_int = new TLegend(0.2,0.6,0.8,0.9);
	leg_int->SetFillStyle(0);
	leg_int->AddEntry(data_int, "Fake Data", "ep");
	leg_int->AddEntry(mc_int_inel, "Inel","f");
	leg_int->AddEntry(mc_int_el, "El","f");

	leg_int->AddEntry(mc_int_midcosmic,"misID:cosmic","f");
	leg_int->AddEntry(mc_int_midp, "misID:p","f");
	leg_int->AddEntry(mc_int_midpi, "misID:#pi","f");

	leg_int->AddEntry(mc_int_midmu,"misID:#mu","f");
	leg_int->AddEntry(mc_int_mideg, "misID:e/#gamma","f");
	leg_int->AddEntry(mc_int_midother, "misID:other","f");

	leg_int->SetNColumns(3);
	leg_int->Draw();
	c_int->Print(Form("%sfakedata_recosliceid_int.eps",outpath.Data()));


	//[2]data with bkg subtraction
	TCanvas *c_int2=new TCanvas(Form("c_int2"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.02); 
	gStyle->SetPadLeftMargin(0.12); 

	c_int2->Divide(1,1);
	c_int2->cd(1);
	TH2D *f2d_int2=new TH2D("f2d_int2",Form("%s","INT"),nthinslices+1,0,nthinslices+1,80,0,ymax_evt);
	f2d_int2->SetTitle("Proton Inelasic Scatterings; Reco SliceID; Events");
	f2d_int2->GetYaxis()->SetTitleOffset(1.3);
	f2d_int2->Draw();
	mc_int_inel->SetLineColor(2); //mc 
	mc_int_inel->Draw("hist same");
	data_int->Draw("ep same");
	data_int_bkgfree->SetMarkerColor(1);
	data_int_bkgfree->SetLineColor(1);
	data_int_bkgfree->SetLineWidth(3);
	data_int_bkgfree->SetMarkerStyle(24);
	data_int_bkgfree->Draw("ep same");

	TLegend *leg_int2 = new TLegend(0.2,0.6,0.8,0.9);
	leg_int2->SetFillStyle(0);
	leg_int2->AddEntry(data_int, "Selected", "ep");
	leg_int2->AddEntry(data_int_bkgfree, "Selected+BKG Subtraction","ep");
	leg_int2->AddEntry(mc_int_inel, "Selected MC Truth [Inel]","f");
	leg_int2->Draw();

	c_int2->Print(Form("%sdata_recosliceid_bkgsub_int.eps",outpath.Data()));


	//[3]data unfolding (true sliceID)
	TCanvas *c_int3=new TCanvas(Form("c_int3"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.02); 
	gStyle->SetPadLeftMargin(0.12); 

	c_int3->Divide(1,1);
	c_int3->cd(1);
	//TH2D *f2d_int3=new TH2D("f2d_int3",Form("%s","INT"),22,-1,21,100,0,ymax_evt);
	TH2D *f2d_int3=new TH2D("f2d_int3",Form("%s","INT"),nthinslices+1,-1,nthinslices+1,100,0,ymax_evt);
	f2d_int3->SetTitle("Proton Inelasic Scatterings; True SliceID; Events");
	f2d_int3->GetYaxis()->SetTitleOffset(1.3);
	f2d_int3->Draw();
	data_int_bkgfree->Draw("ep same");

	data_int_uf->SetLineColor(4);
	data_int_uf->SetMarkerColor(4);
        data_int_uf->SetMarkerStyle(21);
	data_int_uf->Draw("ep same");

	mc_true_interactions->SetLineColor(2);
	mc_true_interactions->Draw("hist same");

	TLegend *leg_int3 = new TLegend(0.2,0.6,0.8,0.9);
	leg_int3->SetFillStyle(0);
	leg_int3->AddEntry(data_int_bkgfree, "Selected+BKG Subtraction","ep");
	leg_int3->AddEntry(data_int_uf, "Unfolding","ep");
	leg_int3->AddEntry(mc_true_interactions,"MC Truth [Inel]","l");
	//leg_int2->SetNColumns(3);
	leg_int3->Draw();

	c_int3->Print(Form("%sdata_recosliceid_uf_int.eps",outpath.Data()));


	//[4] eff & purity
	TCanvas *c_int4=new TCanvas(Form("c_int4"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 
	gStyle->SetPadRightMargin(0.1); 
	gStyle->SetPadLeftMargin(0.1); 

	c_int4->Divide(1,1);
	c_int4->cd(1);
	float ymax_eff_int=1.2;
	//TH2D *f2d_int4=new TH2D("f2d_int4",Form("%s",""),22,-1,21,10,0,ymax_eff_int);
	TH2D *f2d_int4=new TH2D("f2d_int4",Form("%s",""),nthinslices+2,-1,nthinslices+1,10,0,ymax_eff_int);
	f2d_int4->SetTitle("Proton Inelastic Scatterings; Reco SliceID; Efficiency");
	f2d_int4->Draw("");
	eff_int->Draw("same");

	double av_eff_int=0;
	double err_av_eff_int=0;
	double cnt_eff_int=0;
	for (int k=0; k<eff_int->GetNbinsX(); ++k) {
		if (eff_int->GetBinContent(k+1)) { 
			cnt_eff_int++;
			av_eff_int+=eff_int->GetBinContent(k+1);
			err_av_eff_int+=pow(eff_int->GetBinError(k+1),2);
		}
	}
	av_eff_int/=cnt_eff_int;
	err_av_eff_int=sqrt(err_av_eff_int)/cnt_eff_int;
	cout<<"av_eff_int="<<av_eff_int<<" +- "<<err_av_eff_int<<endl;

	TLine *line_av_eff_int=new TLine(-1,av_eff_int,nthinslices+1,av_eff_int);
	line_av_eff_int->SetLineColor(7);
	line_av_eff_int->SetLineWidth(2);
	line_av_eff_int->SetLineStyle(2);
	line_av_eff_int->Draw("same");


	c_int4->Update();
   	//scale hint1 to the pad coordinates
   	float rightmax_int = ymax_eff_int; //display max in y-axis for pur. dist.
   	float scale_int = gPad->GetUymax()/rightmax_int;
   	pur_int->Scale(scale_int);
	pur_int->SetMarkerColor(2);
	pur_int->SetLineColor(2);
	pur_int->Draw("same");

	double av_pur_int=0;
	double err_av_pur_int=0;
	double cnt_pur_int=0;
	for (int k=0; k<pur_int->GetNbinsX(); ++k) {
		av_pur_int+=pur_int->GetBinContent(k+1);
		err_av_pur_int+=pow(pur_int->GetBinError(k+1),2);
		if (pur_int->GetBinContent(k+1)) cnt_pur_int++;

	}
	av_pur_int/=cnt_pur_int;
	err_av_pur_int=sqrt(err_av_pur_int)/cnt_pur_int;
	cout<<"av_pur_int="<<av_pur_int<<" +- "<<err_av_pur_int<<endl;

 
   	//draw an axis on the right side
   	TGaxis *axis_int = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
         gPad->GetUxmax(), gPad->GetUymax(),0,rightmax_int,510,"+L");
   	axis_int->SetLineColor(2);
   	axis_int->SetLabelColor(2);
	axis_int->SetTitle("Purity");
	axis_int->SetTextColor(2);
   	axis_int->Draw();

	//TLine *line_av_pur_int=new TLine(-1,av_pur_int,21,av_pur_int);
	TLine *line_av_pur_int=new TLine(-1,av_pur_int,nthinslices+1,av_pur_int);
	line_av_pur_int->SetLineColor(3);
	line_av_pur_int->SetLineWidth(2);
	line_av_pur_int->SetLineStyle(2);
	line_av_pur_int->Draw("same");

	c_int4->Print(Form("%sdata_eff_pur_int.eps",outpath.Data()));



	//reco inc and int --------------------------------------------------//
        //gStyle->SetPadLeftMargin(0.13);
        gStyle->SetPadRightMargin(0.0);
 	gStyle->SetPadLeftMargin(0.16);
    	gStyle->SetPadRightMargin(0.04);//0.02
        TCanvas *c_inct = new TCanvas("c_inct", "c_inct", 1400, 900);
	c_inct->Divide(1,1);
	c_inct->cd(1);
	c_inct->cd(1)->SetLogy();

        //TH2D *f2d_inct=new TH2D("f2d_inct","",21,0,21,6000,100,1500000);
        TH2D *f2d_inct=new TH2D("f2d_inct","",nthinslices+1,0,nthinslices+1,6000,100,1500000);
        f2d_inct->SetTitle("Reconstructed Incident/Inelastic Protons; SliceID; Events");
        f2d_inct->GetYaxis()->SetTitleOffset(1.3);
        f2d_inct->Draw("");
        gr_reco_inc->SetLineWidth(3);
        gr_reco_int->SetLineWidth(3);
        gr_reco_inc->SetMarkerColor(1);
        gr_reco_inc->SetMarkerStyle(20);
        gr_reco_int->SetMarkerColor(2);
        gr_reco_int->SetLineColor(2);
        gr_reco_int->SetMarkerStyle(20);
	gr_reco_inc->Draw("p same");
	gr_reco_int->Draw("p same");

        TLegend *leg_inct = new TLegend(0.7,0.7,0.9,0.9);
	//leg_inct->SetLineColor(1);
        leg_inct->SetFillStyle(0);
	//leg_inct->SetHeader("Reconstructed Incident/Interacting Protons");
        leg_inct->AddEntry(gr_reco_inc, "Incident", "ep");
        leg_inct->AddEntry(gr_reco_int, "Inelastic", "ep");
        leg_inct->Draw();
	c_inct->Print(Form("%sreco_inc_int.eps",outpath.Data()));

	//sliceID vs KE
	gStyle->SetEndErrorSize(0);
	TCanvas *c_id_ke=new TCanvas(Form("c_id_ke"),"",900, 600);
        gStyle->SetPadRightMargin(0.0);
 	gStyle->SetPadLeftMargin(0.16);
    	gStyle->SetPadRightMargin(0.04);//0.02
	c_id_ke->Divide(1,1);
	c_id_ke->cd(1);

	gr_sliceid_recoincE->SetMarkerColor(2);
	gr_sliceid_recoincE->SetLineColor(2);
	gr_sliceid_trueincE->SetMarkerColor(1);
	
        //TH2D *f2d_id_ke=new TH2D("f2d_id_ke","",21,0,21,500,50,500);	
        TH2D *f2d_id_ke=new TH2D("f2d_id_ke","",nthinslices+1,0,nthinslices+1,500,50,500);	
	f2d_id_ke->SetTitle("; SliceID; Proton Kinetic Energy [MeV]");
	f2d_id_ke->Draw();
	gr_sliceid_trueincE->Draw("p same");
	gr_sliceid_recoincE->Draw("p same");

        TLegend *leg_idke = new TLegend(0.6,0.6,0.9,0.85);
        leg_idke->SetFillStyle(0);
        leg_idke->AddEntry(gr_sliceid_trueincE, "KE Truth", "pe");
        leg_idke->AddEntry(gr_sliceid_recoincE, "KE Reco", "pe");
        leg_idke->Draw();
	c_id_ke->Print("%strue_recoKE_sliceID.eps",outpath.Data());

	//reco KE - true KE ----------------------------------------------//
        gStyle->SetPadRightMargin(0.0);
 	gStyle->SetPadLeftMargin(0.16);
    	gStyle->SetPadRightMargin(0.04);//0.02
        TCanvas *c_dke = new TCanvas("c_dke", "c_dke", 1400, 900);
	c_dke->Divide(1,1);
	c_dke->cd(1);

        //TH2D *f2d_id_dke=new TH2D("f2d_id_dke","",21,0,21,30,-20,10);	
        TH2D *f2d_id_dke=new TH2D("f2d_id_dke","",nthinslices+1,0,nthinslices+1,30,-20,10);	
	f2d_id_dke->SetTitle("; SliceID; Reco KE - True KE [MeV]");
	f2d_id_dke->Draw();
	gr_sliceid_recotrueincE->SetMarkerColor(2);
	gr_sliceid_recotrueincE->SetLineColor(2);
	gr_sliceid_recotrueincE->Draw("p same");
	c_dke->Print("%sdKE_sliceID.eps",outpath.Data());
	

/*
	//demo KE 
        gStyle->SetPadRightMargin(0.0);
 	gStyle->SetPadLeftMargin(0.16);
    	gStyle->SetPadRightMargin(0.04);//0.02

        TCanvas *c_ke = new TCanvas("c_ke", "c_ke", 1400, 900);
	c_ke->Divide(1,1);
	c_ke->cd(1)->SetLogy();
	c_ke->cd(1);
	
	TH1D *TrueKE_DEMO_0=(TH1D *)f_mc->Get(Form("true_incE_%d",6));
	TrueKE_DEMO_0->SetLineColor(2);
	TrueKE_DEMO_0->SetLineWidth(2);
	TrueKE_DEMO_0->SetFillColor(2);
	TrueKE_DEMO_0->SetFillStyle(3353);
	TrueKE_DEMO_0->GetXaxis()->SetTitle("Proton Kinetic Energy [MeV]");
	TrueKE_DEMO_0->GetYaxis()->SetTitle("Events");
	TrueKE_DEMO_0->Draw("hist");
	cout<<"\n\nTrueKE_DEMO_0:"<<TrueKE_DEMO_0->GetMean()<<" +- "<<TrueKE_DEMO_0->GetMeanError()<<endl;
	cout<<"gPad->GetUymax():"<<gPad->GetUymax()<<endl;
	//TLine lll(TrueKE_DEMO_0->GetMean(),0,TrueKE_DEMO_0->GetMean(),90000);
	TLine lll(TrueKE_DEMO_0->GetMean(),0,TrueKE_DEMO_0->GetMean(),30000);
	lll.SetLineStyle(2);
	lll.SetLineWidth(2);
        lll.Draw("same");

         TLegend *leg_slc = new TLegend(0.5,0.6,0.95,0.9);
         leg_slc->SetFillStyle(0);
         leg_slc->AddEntry((TObject*)0, Form("<KE>=%.2f #pm %.2f MeV", TrueKE_DEMO_0->GetMean(), TrueKE_DEMO_0->GetMeanError()), "");
         leg_slc->Draw();


	c_ke->Print("%strueKE_sliceID6.eps",outpath.Data());



	//include SCE OFF sample
	
	//TFile *f_mc_sceoff = TFile::Open("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/PDSPProd4a_MC_6GeV_gen_datadriven_reco1_sce_off_v1/prod4a_mc_thinslice_dx5cm_20slcs_sceoff.root");


        for (int i = 0; i<nthinslices; ++i){

        	TCanvas *c_slc = new TCanvas("c_slc", "c_slc", 1400, 900);
		c_slc->Divide(1,1);
		c_slc->cd(1)->SetLogy();
		c_slc->cd(1);

		TString str_slcslc;
        	if (i==0) str_slcslc=Form("%strueKE_slc-by-slc.pdf(",outpath.Data()); 
		if (i>0&&i<nthinslices-1) str_slcslc=Form("%strueKE_slc-by-slc.pdf",outpath.Data()); 
		if (i==nthinslices-1) str_slcslc=Form("%strueKE_slc-by-slc.pdf)",outpath.Data()); 

		TH1D *htmp=(TH1D* )f_mc->Get(Form("true_incE_%d",i));
		htmp->GetXaxis()->SetTitle("Proton Kinetic Energy [MeV]");
		htmp->GetYaxis()->SetTitle("Events");
		htmp->SetLineColor(1);
		htmp->SetLineWidth(2);
		htmp->SetFillColor(0);
		htmp->SetFillStyle(3353);

		TH1D *htmp2=(TH1D* )f_mc->Get(Form("reco_incE_%d",i));
		htmp2->SetLineColor(2);
		htmp2->SetLineWidth(2);
		htmp2->Scale((float)htmp->Integral()/(float)htmp2->Integral());

		//TH1D *htmp3=(TH1D* )f_mc_sceoff->Get(Form("true_incE_%d",i));
		//htmp3->SetLineColor(4);
		//htmp3->SetLineWidth(4);
		//htmp3->Scale((float)htmp->Integral()/(float)htmp3->Integral());



		htmp->Draw("hist");
		htmp2->Draw("hist same");
		//htmp3->Draw("hist same");


        	TLegend *leg_slc = new TLegend(0.6,0.6,0.9,0.85);
        	leg_slc->SetFillStyle(0);

		TLegendEntry* text_ke[10];		
		text_ke[0]=leg_slc->AddEntry(htmp, Form("True KE"), "l");
        	text_ke[1]=leg_slc->AddEntry((TObject*)0, Form("True <KE>:%.2f MeV",htmp->GetMean()), ""); 
        	text_ke[2]=leg_slc->AddEntry((TObject*)0, Form("True RMS KE:%.2f MeV",htmp->GetRMS()), ""); 

		text_ke[3]=leg_slc->AddEntry(htmp2, Form("Reco KE"), "l"); text_ke[3]->SetTextColor(2);
        	text_ke[4]=leg_slc->AddEntry((TObject*)0, Form("Reco <KE>:%.2f MeV",htmp2->GetMean()), ""); text_ke[4]->SetTextColor(2);
        	text_ke[5]=leg_slc->AddEntry((TObject*)0, Form("Reco RMS KE:%.2f MeV",htmp2->GetRMS()), ""); text_ke[5]->SetTextColor(2);


		//text_ke[6]=leg_slc->AddEntry(htmp, Form("True KE (SCE OFF)"), "l"); text_ke[6]->SetTextColor(4);
        	//text_ke[7]=leg_slc->AddEntry((TObject*)0, Form("True <KE>:%.2f MeV",htmp3->GetMean()), ""); text_ke[7]->SetTextColor(4);
        	//text_ke[8]=leg_slc->AddEntry((TObject*)0, Form("True RMS KE:%.2f MeV",htmp3->GetRMS()), ""); text_ke[8]->SetTextColor(4);

		leg_slc->Draw();

		c_slc->Print(str_slcslc.Data());
		delete c_slc;
		delete leg_slc;

	}
*/


	//xs result -------------------------------------------------------------------------------------//
        TFile f_xs("/dune/data2/users/hyliao/GeantReweight/xs_cascade/proton_cross_section.root");
        TGraph *total_inel_KE = (TGraph*)f_xs.Get("inel_KE");

        //gStyle->SetPadLeftMargin(0.13);
        gStyle->SetPadRightMargin(0.0);
 	gStyle->SetPadLeftMargin(0.16);
    	gStyle->SetPadRightMargin(0.04);//0.02
        TCanvas *c_xs = new TCanvas("c_xs", "c_xs", 1400, 900);
	c_xs->Divide(1,1);
	c_xs->cd(1);

        float ymax=1200;
        float xmax=460;
        //float xmax=620;
        float xmin=0;

        TH2D *f2d_xs=new TH2D("f2d_xs","",450,xmin,xmax,ymax,0,ymax);
        f2d_xs->GetXaxis()->SetTitle("Proton Kinetic Energy [MeV]");
        f2d_xs->GetYaxis()->SetTitleOffset(1.3);
        f2d_xs->GetYaxis()->SetTitle("P-Ar inelastic cross section [mb]");
        f2d_xs->Draw("");

        total_inel_KE->SetLineColor(2);
        total_inel_KE->Draw("c same");
        gr_truexs->SetLineWidth(2);
        gr_truexs->SetMarkerColor(3);
        gr_truexs->SetLineColor(3);
        gr_truexs->SetMarkerStyle(21);
        gr_truexs->SetMarkerSize(1.5);
        gr_truexs->Draw("p same");

        gr_recoxs->SetLineWidth(2);
        gr_recoxs->SetMarkerColor(1);
        gr_recoxs->SetLineColor(1);
        gr_recoxs->SetMarkerStyle(20);
        gr_recoxs->SetMarkerSize(1.5);
        gr_recoxs->Draw("p same");

	double ke_ff_mean_stop=4.12972e+02;
	double ke_ff_sigma_stop=4.08988e+01;
	TLine l_mean(ke_ff_mean_stop,0,ke_ff_mean_stop,ymax);
	TLine l_left(ke_ff_mean_stop-3.*ke_ff_sigma_stop,0,ke_ff_mean_stop-3.*ke_ff_sigma_stop,ymax);
	TLine l_right(ke_ff_mean_stop+3.*ke_ff_sigma_stop,0,ke_ff_mean_stop+3.*ke_ff_sigma_stop,ymax);
	l_mean.SetLineColor(1);
	l_left.SetLineColor(1);
        l_right.SetLineColor(1);

	l_left.SetLineStyle(2);	
	l_right.SetLineStyle(2);	

	//l_mean.Draw("same");
	//l_left.Draw("same");
	//l_right.Draw("same");

        TLegend *leg_xs = new TLegend(0.6,0.6,0.9,0.85);
        leg_xs->SetFillStyle(0);
        leg_xs->AddEntry(total_inel_KE, "Geant4", "l");
        leg_xs->AddEntry(gr_truexs, "MC Truth", "pe");
        //leg_xs->AddEntry(gr_recoxs, "MC Reco", "pe");
        leg_xs->Draw();


	//c_xs->Print(Form("%sxs_data.eps",outpath.Data()));
	c_xs->Print(Form("%sxs_data_noreco.eps",outpath.Data()));


	//output files ----------------------------------------------------------------------------------------------//
	gr_recoxs->SetName("gr_recoxs");
	gr_truexs->SetName("gr_truexs");


/*
	//TFile *fout_xs = new TFile("./xs_rootfiles/xs_Eslice_dE20MeV_40slcs_nobmrw_stslcplus0.5.root","recreate");
	TFile *fout_xs = new TFile("./xs_rootfiles/xs_Eslice_dE20MeV_40slcs_nobmrw_stslcplus0.5_new.root","recreate");
	gr_recoxs->Write();
	gr_truexs->Write();
	fout_xs->Close();
*/



/*

	//ke vs length ----------------------------------------------------------//
	TH2D *KEbb_truetrklen_all=(TH2D *)f_mc->Get("KEbb_truetrklen_all");
	TH2D *KEbb_truetrklen_inel=(TH2D *)f_mc->Get("KEbb_truetrklen_inel");
	TH2D *KEbb_recotrklen_all=(TH2D *)f_mc->Get("KEbb_recotrklen_all");
	TH2D *KEbb_recotrklen_inel=(TH2D *)f_mc->Get("KEbb_recotrklen_inel");

        gStyle->SetPadRightMargin(0.0);
 	gStyle->SetPadLeftMargin(0.16);
    	gStyle->SetPadRightMargin(0.04);//0.02

        TCanvas *c_kebb_trklen = new TCanvas("c_kebb_trklen", "", 1200, 900);
	c_kebb_trklen->Divide(1,1);
	c_kebb_trklen->cd(1);

        float yymax=800;
        float xxmax=120;
        float xxmin=0;

        TH2D *f2d_ketrklen=new TH2D("f2d_ketrklen","",1200,xxmin,xxmax,yymax,0,yymax);
        f2d_ketrklen->GetYaxis()->SetTitle("Proton Kinetic Energy [MeV]");
        f2d_ketrklen->GetXaxis()->SetTitle("True Track Length [cm]");
        f2d_ketrklen->GetYaxis()->SetTitleOffset(1.3);

        f2d_ketrklen->Draw("");
	KEbb_truetrklen_all->Draw("colz same");
	c_kebb_trklen->Print(Form("%skebb_truetrklen_all.png",outpath.Data()));

	KEbb_truetrklen_inel->Draw("colz same");
	c_kebb_trklen->Print(Form("%skebb_truetrklen_inel.png",outpath.Data()));
	
        f2d_ketrklen->GetXaxis()->SetTitle("Reco Track Length [cm]");	
	KEbb_recotrklen_all->Draw("colz same");
	c_kebb_trklen->Print(Form("%skebb_recotrklen_all.png",outpath.Data()));
	
	KEbb_recotrklen_inel->Draw("colz same");
	c_kebb_trklen->Print(Form("%skebb_recotrklen_inel.png",outpath.Data()));
	
	TBox *box2[nthinslices];
	TText *tt2[nthinslices];	
        for (int i = 0; i<nthinslices; ++i){
	  //double id=i+0.5;
	  double ke=Emax-((double)i+0.5)*thinslicewidth;

	  box2[i]=new TBox(xxmin,ke-thinslicewidth*.5,xxmax,ke+thinslicewidth*.5);
	  box2[i]->SetLineColor(6);
	  box2[i]->SetLineWidth(2);
	  box2[i]->SetFillColor(5);
          box2[i]->SetFillStyle(0);
	  box2[i]->Draw("same");
	
	  tt2[i]=new TText(xxmin+10, ke-thinslicewidth*.3,Form("%d",i));
	  tt2[i]->SetTextColor(6);
	  tt2[i]->SetTextSize(0.03);
	  tt2[i]->Draw("same");
	  //box[i]->SetFillColorAlpha(0,0.5);
	}	
	c_kebb_trklen->Print(Form("%skebb_recotrklen_inel_withID.png",outpath.Data()));
	


	//zoom
        yymax=600;
        xxmax=10;
        xxmin=0;

        TH2D *f2d_ketrklen_zoom=new TH2D("f2d_ketrklen_zoom","",200,xxmin,xxmax,yymax,0,yymax);
        f2d_ketrklen_zoom->GetYaxis()->SetTitle("Proton Kinetic Energy [MeV]");
        f2d_ketrklen_zoom->GetXaxis()->SetTitle("True Track Length [cm]");
        f2d_ketrklen_zoom->GetYaxis()->SetTitleOffset(1.3);

        f2d_ketrklen_zoom->Draw("");
	KEbb_truetrklen_all->Draw("colz same");
	c_kebb_trklen->Print(Form("%skebb_truetrklen_all_zoom.png",outpath.Data()));

	KEbb_truetrklen_inel->Draw("colz same");
	c_kebb_trklen->Print(Form("%skebb_truetrklen_inel_zoom.png",outpath.Data()));
	
        f2d_ketrklen->GetXaxis()->SetTitle("Reco Track Length [cm]");	
	KEbb_recotrklen_all->Draw("colz same");
	c_kebb_trklen->Print(Form("%skebb_recotrklen_all_zoom.png",outpath.Data()));
	
	KEbb_recotrklen_inel->Draw("colz same");
	c_kebb_trklen->Print(Form("%skebb_recotrklen_inel_zoom.png",outpath.Data()));


*/
	


}
