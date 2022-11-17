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
#include "./headers/RooUnfoldBayes.h"
#include "./headers/RooUnfoldResponse.h"

#include "./headers/BasicParameters.h"
#include "./headers/BasicFunctions.h"
#include "./headers/BetheBloch.h"
#include "./headers/ESliceParams.h"

void Wiener_SVD_Wrapper() {

	//Input files ----------------------------------------------------------------------------------------------------------------//
        //TString fmc="./prod4areco2_mc_ESliceE_dE20MeV_40slcs_beamxy_bmrwkebeamff_v09_39_01.root";
        //TString fdata="/dune/data2/users/hyliao/protonana/v09_39_01/XS/prod4a_Eslice_dE20MeV_40slcs_beamxy_runAll_v09_39_01.root";
        TString fmc="./prod4areco2_mc_ESliceE_dE20MeV_40slcs_beamxy_nobmrw_kebeamff_v09_39_01.root";
        TString fdata="./prod4areco2_mc_ESliceE_dE20MeV_40slcs_beamxy_nobmrw_kebeamff_v09_39_01.root";

	//Output files -------------------------------------------//
	TString outpath="./Output_Wiener_SVD/MC_nobmrw/";
        TString fout_inc=outpath+"wiener_svd_inc.root";
        TString fout_inc_st=outpath+"wiener_svd_inc_st.root";
        TString fout_int=outpath+"wiener_svd_int.root";

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


	//get mc true slice IDs -------------------------------------------------------------//
	//inc
	TH1D* mc_inc_true_all=(TH1D*)f_mc->Get(Form("%s",str_inc_true.Data()));
	TH1D* mc_inc_true_inel=(TH1D*)f_mc->Get(Form("%s_inel",str_inc_true.Data()));
	TH1D* mc_inc_true_el=(TH1D*)f_mc->Get(Form("%s_el",str_inc_true.Data()));
	TH1D* mc_inc_true_midcosmic=(TH1D*)f_mc->Get(Form("%s_midcosmic",str_inc_true.Data()));
	TH1D* mc_inc_true_midpi=(TH1D*)f_mc->Get(Form("%s_midpi",str_inc_true.Data()));
	TH1D* mc_inc_true_midp=(TH1D*)f_mc->Get(Form("%s_midp",str_inc_true.Data()));
	TH1D* mc_inc_true_midmu=(TH1D*)f_mc->Get(Form("%s_midmu",str_inc_true.Data()));
	TH1D* mc_inc_true_mideg=(TH1D*)f_mc->Get(Form("%s_mideg",str_inc_true.Data()));
	TH1D* mc_inc_true_midother=(TH1D*)f_mc->Get(Form("%s_midother",str_inc_true.Data()));

	//inc_st
	TH1D* mc_st_inc_true_all=(TH1D*)f_mc->Get(Form("%s",str_st_inc_true.Data()));
	TH1D* mc_st_inc_true_inel=(TH1D*)f_mc->Get(Form("%s_inel",str_st_inc_true.Data()));
	TH1D* mc_st_inc_true_el=(TH1D*)f_mc->Get(Form("%s_el",str_st_inc_true.Data()));
	TH1D* mc_st_inc_true_midcosmic=(TH1D*)f_mc->Get(Form("%s_midcosmic",str_st_inc_true.Data()));
	TH1D* mc_st_inc_true_midpi=(TH1D*)f_mc->Get(Form("%s_midpi",str_st_inc_true.Data()));
	TH1D* mc_st_inc_true_midp=(TH1D*)f_mc->Get(Form("%s_midp",str_st_inc_true.Data()));
	TH1D* mc_st_inc_true_midmu=(TH1D*)f_mc->Get(Form("%s_midmu",str_st_inc_true.Data()));
	TH1D* mc_st_inc_true_mideg=(TH1D*)f_mc->Get(Form("%s_mideg",str_st_inc_true.Data()));
	TH1D* mc_st_inc_true_midother=(TH1D*)f_mc->Get(Form("%s_midother",str_st_inc_true.Data()));

	//int
	TH1D* mc_int_true_all=(TH1D*)f_mc->Get(Form("%s",str_int_true.Data()));
	TH1D* mc_int_true_inel=(TH1D*)f_mc->Get(Form("%s_inel",str_int_true.Data()));
	TH1D* mc_int_true_el=(TH1D*)f_mc->Get(Form("%s_el",str_int_true.Data()));
	TH1D* mc_int_true_midcosmic=(TH1D*)f_mc->Get(Form("%s_midcosmic",str_int_true.Data()));
	TH1D* mc_int_true_midpi=(TH1D*)f_mc->Get(Form("%s_midpi",str_int_true.Data()));
	TH1D* mc_int_true_midp=(TH1D*)f_mc->Get(Form("%s_midp",str_int_true.Data()));
	TH1D* mc_int_true_midmu=(TH1D*)f_mc->Get(Form("%s_midmu",str_int_true.Data()));
	TH1D* mc_int_true_mideg=(TH1D*)f_mc->Get(Form("%s_mideg",str_int_true.Data()));
	TH1D* mc_int_true_midother=(TH1D*)f_mc->Get(Form("%s_midother",str_int_true.Data()));

	int n_mc_inc_true_inel=mc_inc_true_inel->Integral();
	int n_mc_inc_true_el=mc_inc_true_el->Integral();
	int n_mc_inc_true_midcosmic=mc_inc_true_midcosmic->Integral();
	int n_mc_inc_true_midpi=mc_inc_true_midpi->Integral();
	int n_mc_inc_true_midp=mc_inc_true_midp->Integral();
	int n_mc_inc_true_midmu=mc_inc_true_midmu->Integral();
	int n_mc_inc_true_mideg=mc_inc_true_mideg->Integral();
	int n_mc_inc_true_midother=mc_inc_true_midother->Integral();
	int n_mc_inc_true=n_mc_inc_true_inel+n_mc_inc_true_el+n_mc_inc_true_midcosmic+n_mc_inc_true_midpi+n_mc_inc_true_midp+n_mc_inc_true_midmu+n_mc_inc_true_mideg+n_mc_inc_true_midother;
	double norm_mc_inc_true=(double)n_data_inc/(double)n_mc_inc_true;

	int n_mc_st_inc_true_inel=mc_st_inc_true_inel->Integral();
	int n_mc_st_inc_true_el=mc_st_inc_true_el->Integral();
	int n_mc_st_inc_true_midcosmic=mc_st_inc_true_midcosmic->Integral();
	int n_mc_st_inc_true_midpi=mc_st_inc_true_midpi->Integral();
	int n_mc_st_inc_true_midp=mc_st_inc_true_midp->Integral();
	int n_mc_st_inc_true_midmu=mc_st_inc_true_midmu->Integral();
	int n_mc_st_inc_true_mideg=mc_st_inc_true_mideg->Integral();
	int n_mc_st_inc_true_midother=mc_st_inc_true_midother->Integral();
	int n_mc_st_inc_true=n_mc_st_inc_true_inel+n_mc_st_inc_true_el+n_mc_st_inc_true_midcosmic+n_mc_st_inc_true_midpi+n_mc_st_inc_true_midp+n_mc_st_inc_true_midmu+n_mc_st_inc_true_mideg+n_mc_st_inc_true_midother;
	double norm_mc_st_inc_true=(double)n_data_st_inc/(double)n_mc_st_inc_true;

	int n_mc_int_true_inel=mc_int_true_inel->Integral();
	int n_mc_int_true_el=mc_int_true_el->Integral();
	int n_mc_int_true_midcosmic=mc_int_true_midcosmic->Integral();
	int n_mc_int_true_midpi=mc_int_true_midpi->Integral();
	int n_mc_int_true_midp=mc_int_true_midp->Integral();
	int n_mc_int_true_midmu=mc_int_true_midmu->Integral();
	int n_mc_int_true_mideg=mc_int_true_mideg->Integral();
	int n_mc_int_true_midother=mc_int_true_midother->Integral();
	int n_mc_int_true=n_mc_int_true_inel+n_mc_int_true_el+n_mc_int_true_midcosmic+n_mc_int_true_midpi+n_mc_int_true_midp+n_mc_int_true_midmu+n_mc_int_true_mideg+n_mc_int_true_midother;
	double norm_mc_int_true=(double)n_data_int/(double)n_mc_int_true;

	//-----------------------------------------------------------------------------------//

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
	

	mc_inc_true_inel->Scale(norm_mc_inc_true);
	mc_inc_true_el->Scale(norm_mc_inc_true);
	mc_inc_true_midcosmic->Scale(norm_mc_inc_true);
	mc_inc_true_midpi->Scale(norm_mc_inc_true);
	mc_inc_true_midp->Scale(norm_mc_inc_true);
	mc_inc_true_midmu->Scale(norm_mc_inc_true);
	mc_inc_true_mideg->Scale(norm_mc_inc_true);
	mc_inc_true_midother->Scale(norm_mc_inc_true);

	mc_st_inc_true_inel->Scale(norm_mc_st_inc_true);
	mc_st_inc_true_el->Scale(norm_mc_st_inc_true);
	mc_st_inc_true_midcosmic->Scale(norm_mc_st_inc_true);
	mc_st_inc_true_midpi->Scale(norm_mc_st_inc_true);
	mc_st_inc_true_midp->Scale(norm_mc_st_inc_true);
	mc_st_inc_true_midmu->Scale(norm_mc_st_inc_true);
	mc_st_inc_true_mideg->Scale(norm_mc_st_inc_true);
	mc_st_inc_true_midother->Scale(norm_mc_st_inc_true);



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


	mc_int_true_inel->Scale(norm_mc_int_true);
	mc_int_true_el->Scale(norm_mc_int_true);
	mc_int_true_midcosmic->Scale(norm_mc_int_true);
	mc_int_true_midpi->Scale(norm_mc_int_true);
	mc_int_true_midp->Scale(norm_mc_int_true);
	mc_int_true_midmu->Scale(norm_mc_int_true);
	mc_int_true_mideg->Scale(norm_mc_int_true);
	mc_int_true_midother->Scale(norm_mc_int_true);

	//Note: Numerical value of errorbar after scaling is correct, i.e. s*sqrt(n)	

	//MC truth ---------------------------------------------------------------------------------------
	//inc
	double scale_true=(double)n_data_inc/(double)mc_truesliceID_all->Integral();
	mc_truesliceID_all->Scale(scale_true);
	mc_true_st_sliceID_all->Scale(scale_true);

	//int
	mc_truesliceID_inel->Scale(scale_true);
	//------------------------------------------------------------------------------------------------


	
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
	data_inc_bkgfree->Add(mc_inc_midp, scal_fact_misidp);
	//data_inc_bkgfree->Add(mc_inc_midmu, -1);
	//data_inc_bkgfree->Add(mc_inc_mideg, -1);
	//data_inc_bkgfree->Add(mc_inc_midother, -1);
	//data_inc_bkgfree->Multiply(pur_inc);
	//
	data_st_inc_bkgfree->Add(mc_st_inc_midp, scal_fact_misidp);

	//
	//Note: Numerical value of errorbar after subtraction is correct, i.e. (s-b)+-sqrt(s+b)	
	
	//int
	data_int_bkgfree->Add(mc_int_el, scal_fact_misidp);
	//data_int_bkgfree->Add(mc_int_midcosmic, -1);
	//data_int_bkgfree->Add(mc_int_midpi, -1);
	data_int_bkgfree->Add(mc_int_midp, scal_fact_misidp);
	//data_int_bkgfree->Add(mc_int_midmu, -1);
	//data_int_bkgfree->Add(mc_int_mideg, -1);
	//data_int_bkgfree->Add(mc_int_midother, -1);
	
	//MC ------------------------------------------------------------------------------------------------------------------------//
	//reco
	//inc
	TH1D* mc_inc_bkgfree=(TH1D *)mc_inc_inel->Clone("mc_inc_bkgfree"); mc_inc_bkgfree->SetName("mc_inc_bkgfree");
        mc_inc_bkgfree->Add(mc_inc_el, 1);
	mc_inc_bkgfree->Add(mc_inc_midcosmic, 1);
	mc_inc_bkgfree->Add(mc_inc_midpi, 1);
	mc_inc_bkgfree->Add(mc_inc_midmu, 1);
	mc_inc_bkgfree->Add(mc_inc_mideg, 1);
	mc_inc_bkgfree->Add(mc_inc_midother, 1);
	
	//inc_st
	TH1D* mc_st_inc_bkgfree=(TH1D *)mc_st_inc_inel->Clone("mc_st_inc_bkgfree"); mc_st_inc_bkgfree->SetName("mc_st_inc_bkgfree");
        mc_st_inc_bkgfree->Add(mc_st_inc_el, 1);
	mc_st_inc_bkgfree->Add(mc_st_inc_midcosmic, 1);
	mc_st_inc_bkgfree->Add(mc_st_inc_midpi, 1);
	mc_st_inc_bkgfree->Add(mc_st_inc_midmu, 1);
	mc_st_inc_bkgfree->Add(mc_st_inc_mideg, 1);
	mc_st_inc_bkgfree->Add(mc_st_inc_midother, 1);
	
	//int
	TH1D* mc_int_bkgfree=(TH1D *)mc_int_inel->Clone("mc_int_bkgfree"); mc_int_bkgfree->SetName("mc_int_bkgfree");
        mc_int_bkgfree->Add(mc_int_el, 1);
	mc_int_bkgfree->Add(mc_int_midcosmic, 1);
	mc_int_bkgfree->Add(mc_int_midpi, 1);
	mc_int_bkgfree->Add(mc_int_midmu, 1);
	mc_int_bkgfree->Add(mc_int_mideg, 1);
	mc_int_bkgfree->Add(mc_int_midother, 1);

	//truth
	//inc
	TH1D* mc_inc_true_bkgfree=(TH1D *)mc_inc_true_inel->Clone("mc_inc_true_bkgfree"); mc_inc_true_bkgfree->SetName("mc_inc_true_bkgfree");
        mc_inc_true_bkgfree->Add(mc_inc_true_el, 1);
	mc_inc_true_bkgfree->Add(mc_inc_true_midcosmic, 1);
	mc_inc_true_bkgfree->Add(mc_inc_true_midpi, 1);
	mc_inc_true_bkgfree->Add(mc_inc_true_midmu, 1);
	mc_inc_true_bkgfree->Add(mc_inc_true_mideg, 1);
	mc_inc_true_bkgfree->Add(mc_inc_true_midother, 1);
	
	//inc_st
	TH1D* mc_st_inc_true_bkgfree=(TH1D *)mc_st_inc_true_inel->Clone("mc_st_inc_true_bkgfree"); mc_st_inc_true_bkgfree->SetName("mc_st_inc_true_bkgfree");
        mc_st_inc_true_bkgfree->Add(mc_st_inc_true_el, 1);
	mc_st_inc_true_bkgfree->Add(mc_st_inc_true_midcosmic, 1);
	mc_st_inc_true_bkgfree->Add(mc_st_inc_true_midpi, 1);
	mc_st_inc_true_bkgfree->Add(mc_st_inc_true_midmu, 1);
	mc_st_inc_true_bkgfree->Add(mc_st_inc_true_mideg, 1);
	mc_st_inc_true_bkgfree->Add(mc_st_inc_true_midother, 1);
	
	//int
	TH1D* mc_int_true_bkgfree=(TH1D *)mc_int_true_inel->Clone("mc_int_true_bkgfree"); mc_int_true_bkgfree->SetName("mc_int_true_bkgfree");
        mc_int_true_bkgfree->Add(mc_int_true_el, 1);
	mc_int_true_bkgfree->Add(mc_int_true_midcosmic, 1);
	mc_int_true_bkgfree->Add(mc_int_true_midpi, 1);
	mc_int_true_bkgfree->Add(mc_int_true_midmu, 1);
	mc_int_true_bkgfree->Add(mc_int_true_mideg, 1);
	mc_int_true_bkgfree->Add(mc_int_true_midother, 1);


	//response matrix from mc --------------------------------------------------------------------------------------------------//
	//Response matrix as a 2D-histogram: (x,y)=(measured,truth)
	RooUnfoldResponse *res_inc=(RooUnfoldResponse*)f_mc->Get("response_SliceID_Inc"); res_inc->SetName("res_inc");
	RooUnfoldResponse *res_st_inc=(RooUnfoldResponse*)f_mc->Get("response_st_SliceID_Inc"); res_st_inc->SetName("res_st_inc");
	RooUnfoldResponse *res_int=(RooUnfoldResponse*)f_mc->Get("response_SliceID_Int"); res_int->SetName("res_int");
	
	//format conversion
	TH2D *h2d_res_inc=(TH2D*)res_inc->Hresponse();
	TH2D *h2d_res_st_inc=(TH2D*)res_st_inc->Hresponse();
	TH2D *h2d_res_int=(TH2D*)res_int->Hresponse();

	//sansity check on the response matrix
	std::cout<<"res_inc->GetDimensionMeasured():"<<res_inc->GetDimensionMeasured()<<std::endl;
	std::cout<<"res_inc->GetDimensionTruth():"<<res_inc->GetDimensionTruth()<<std::endl;
	std::cout<<"res_int->GetDimensionTruth():"<<res_int->GetDimensionTruth()<<std::endl;

	std::cout<<"res_inc->GetNbinsMeasured():"<<res_inc->GetNbinsMeasured()<<std::endl; //x
	std::cout<<"res_inc->GetNbinsTruth():"<<res_inc->GetNbinsTruth()<<std::endl;  //y

        //for (size_t i=0; i<res_inc->GetNbinsMeasured(); ++i) { //x
        	//for (size_t j=0; j<res_inc->GetNbinsTruth(); ++j) { //y
			//double res_inc_each=h2d_res_inc->GetBinContent(i,j);
			//if (res_inc_each>2000) std::cout<<"res_inc_each["<<i<<"]["<<j<<"]="<<res_inc_each<<std::endl;
		//} //y
	//} //x

	//Construct covariance matrix -----------------------------------------------------------------------------//
	Int_t nb=data_inc->GetNbinsX();
	Double_t xlo=data_inc->GetXaxis()->GetXmin(), xhi= data_inc->GetXaxis()->GetXmax(), xb= (xhi-xlo)/nb;
	TH2D *h2d_cov_inc=new TH2D("h2d_cov_inc", "", nb, xlo, xhi, nb, xlo, xhi);	
	TH2D *h2d_cov_inc_st=new TH2D("h2d_cov_inc_st", "", nb, xlo, xhi, nb, xlo, xhi);	
	TH2D *h2d_cov_int=new TH2D("h2d_cov_int", "", nb, xlo, xhi, nb, xlo, xhi);	
	
	std::cout<<"nb="<<nb<<std::endl;
	std::cout<<"xlo="<<xlo<<std::endl;
	std::cout<<"xhi="<<xhi<<std::endl;
	std::cout<<"xb="<<xb<<std::endl;

	//Construct covariance matrix: x-data recoSliceID measurement; y-MC recoSliceID measurement
        for (int i=0; i<nb; ++i) { //x

		double cov_inc=0;
		double cov_inc_st=0;
		double cov_int=0;

		//Data:recoslice with bkg subtraction -------------------------//
		double m_inc=data_inc_bkgfree->GetBinContent(i); //inc
		double m_st_inc=data_st_inc_bkgfree->GetBinContent(i); //inc_st
		double m_int=data_int_bkgfree->GetBinContent(i); //int

		//MC:recoslice ---------------------------------------------------------//
		double p_inc=mc_inc_bkgfree->GetBinContent(i); //inc
		double p_st_inc=mc_st_inc_bkgfree->GetBinContent(i); //inc_st
		double p_int=mc_int_bkgfree->GetBinContent(i); //int
		
        	for (int j=0; j<nb; ++j) { //y
			if (i==j) { //diagonal element
				cov_inc=cov_cnp(m_inc, p_inc);
				cov_inc_st=cov_cnp(m_st_inc, p_st_inc);
				cov_int=cov_cnp(m_int, p_int);
				
				h2d_cov_inc->Fill(i, j, cov_inc);
				h2d_cov_inc_st->Fill(i, j, cov_inc_st);
				h2d_cov_int->Fill(i, j, cov_int);
			} //diagonal element
			if (i!=j) { //non-diagonal element
				h2d_cov_inc->Fill(i, j, 0);
				h2d_cov_inc_st->Fill(i, j, 0);
				h2d_cov_int->Fill(i, j, 0);
			} //non-diagonal element	
        	} //y

	} //x
	//---------------------------------------------------------------------------------------------------------//

	//Output files --------------------------------------------------------//
	//inc
	TFile *f_out_inc = new TFile(fout_inc.Data(),"RECREATE");
		mc_inc_true_bkgfree->Write("htrue_signal");
		data_inc_bkgfree->Write("hmeas");
		h2d_res_inc->Write("hR");
		h2d_cov_inc->Write("hcov_tot");
	f_out_inc->Close();


	//inc_st
	TFile *f_out_inc_st = new TFile(fout_inc_st.Data(),"RECREATE");
		mc_st_inc_true_bkgfree->Write("htrue_signal");
		data_st_inc_bkgfree->Write("hmeas");
		h2d_res_st_inc->Write("hR");
		h2d_cov_inc_st->Write("hcov_tot");
	f_out_inc_st->Close();

	//int
	TFile *f_out_int = new TFile(fout_int.Data(),"RECREATE");
		mc_int_true_bkgfree->Write("htrue_signal");
		data_int_bkgfree->Write("hmeas");
		h2d_res_int->Write("hR");
		h2d_cov_int->Write("hcov_tot");
	f_out_int->Close();
	//---------------------------------------------------------------------//



}
