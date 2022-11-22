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
	TString outpath="./Wiener_SVD_files/MC_nobmrw/";
        TString fout_inc=outpath+"input_wiener_svd_inc.root";
        TString fout_inc_st=outpath+"input_wiener_svd_inc_st.root";
        TString fout_int=outpath+"input_wiener_svd_int.root";

	//Output files2 -------------------------------------------//
        TString fout2_inc=outpath+"uf_wiener_svd_inc.root";
        TString fout2_inc_st=outpath+"uf_wiener_svd_inc_st.root";
        TString fout2_int=outpath+"uf_wiener_svd_int.root";

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


	//get mc true slice IDs ---------------------------------------------------------------------//
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

	TH2D *h2d_sg_inc=(TH2D*)res_inc->Htruth();
	TH2D *h2d_sg_int=(TH2D*)res_int->Htruth();
	TH2D *h2d_sg_st_inc=(TH2D*)res_st_inc->Htruth();

	TH1D *h1d_sg_inc=(TH1D*)res_inc->Htruth();
	TH1D *h1d_sg_int=(TH1D*)res_int->Htruth();
	TH1D *h1d_sg_st_inc=(TH1D*)res_st_inc->Htruth();

	//sansity check on the response matrix
	std::cout<<"\nres_inc->GetDimensionMeasured():"<<res_inc->GetDimensionMeasured()<<std::endl;
	std::cout<<"res_inc->GetDimensionTruth():"<<res_inc->GetDimensionTruth()<<std::endl;
	std::cout<<"res_int->GetDimensionTruth():"<<res_int->GetDimensionTruth()<<"\n"<<std::endl;

	std::cout<<"res_inc->GetNbinsMeasured():"<<res_inc->GetNbinsMeasured()<<std::endl; //x
	std::cout<<"res_inc->GetNbinsTruth():"<<res_inc->GetNbinsTruth()<<std::endl;  //y

	std::cout<<"\nh2d_res_inc->GetNbinsX():"<<h2d_res_inc->GetNbinsX()<<std::endl;
	std::cout<<"h2d_res_inc->GetNbinsY():"<<h2d_res_inc->GetNbinsY()<<"\n"<<std::endl;

	std::cout<<"h1d_sg_int->GetNbinsX():"<<h1d_sg_int->GetNbinsX()<<std::endl;
	std::cout<<"h1d_sg_int->GetNbinsY():"<<h1d_sg_int->GetNbinsY()<<"\n"<<std::endl;


/*
        TH1D *h1d_sg_inc=new TH1D("h1d_sg_inc","", nb_res, xlo_res, xhi_res);
        TH1D *h1d_sg_st_inc=new TH1D("h1d_sg_st_inc","",nb_res, xlo_res, xhi_res);
        TH1D *h1d_sg_int=new TH1D("h1d_sg_int","",nb_res, xlo_res, xhi_res);

	int key=1;
        for (size_t i=1; i<=res_inc->GetNbinsMeasured(); ++i) { //x
        	for (size_t j=1; j<=res_inc->GetNbinsTruth(); ++j) { //y
			h1d_sg_inc->SetBinContent(key, h2d_sg_inc->GetBinContent(i,j));
			h1d_sg_inc->SetBinError(key, h2d_sg_inc->GetBinError(i,j));

			h1d_sg_st_inc->SetBinContent(key, h2d_sg_st_inc->GetBinContent(i,j));
			h1d_sg_st_inc->SetBinError(key, h2d_sg_st_inc->GetBinError(i,j));

			h1d_sg_int->SetBinContent(key, h2d_sg_int->GetBinContent(i,j));
			h1d_sg_int->SetBinError(key, h2d_sg_int->GetBinError(i,j));

			++key;
		} //y
	} //x
	std::cout<<"\nres_inc->GetNbinsMeasured()="<<res_inc->GetNbinsMeasured()<<std::endl;
	std::cout<<"key:"<<key<<"\n"<<std::endl;
*/

	//Response Matrix for Wiener SVD --------------------------------------------------------------------------------------------------//
	//M=R*S: S:truth; M:reco
	Int_t nb_res=h2d_res_inc->GetNbinsX();
	Double_t xlo_res=h2d_res_inc->GetXaxis()->GetXmin(), xhi_res=h2d_res_inc->GetXaxis()->GetXmax(), xb_res=(xhi_res-xlo_res)/nb_res;
	TH2D *h2d_hR_inc=new TH2D("h2d_hR_inc","",nb_res, xlo_res, xhi_res,nb_res, xlo_res, xhi_res);
	TH2D *h2d_hR_st_inc=new TH2D("h2d_hR_st_inc","",nb_res, xlo_res, xhi_res,nb_res, xlo_res, xhi_res);
	TH2D *h2d_hR_int=new TH2D("h2d_hR_int","",nb_res, xlo_res, xhi_res,nb_res, xlo_res, xhi_res);

        for (size_t i=1; i<=res_inc->GetNbinsMeasured(); ++i) { //x
		//double norm_inc=1;
		double norm_sg_inc=h1d_sg_inc->GetBinContent(i);
		double norm_sg_st_inc=h1d_sg_st_inc->GetBinContent(i);
		double norm_sg_int=h1d_sg_int->GetBinContent(i);

        	for (size_t j=1; j<=res_inc->GetNbinsTruth(); ++j) { //y
			double res_inc=h2d_res_inc->GetBinContent(j,i);
			double res_st_inc=h2d_res_st_inc->GetBinContent(j,i);
			double res_int=h2d_res_int->GetBinContent(j,i);

			//response matrix rotation
			if (norm_sg_inc!=0) h2d_hR_inc->SetBinContent(j,i, res_inc/norm_sg_inc);
			if (norm_sg_st_inc!=0) h2d_hR_st_inc->SetBinContent(j,i, res_st_inc/norm_sg_st_inc);
			if (norm_sg_int!=0)  h2d_hR_int->SetBinContent(j,i, res_int/norm_sg_int);
		} //y
	} //x
	//---------------------------------------------------------------------------------------------------------------------------------//




/*
        for (size_t i=0; i<h2d_hR_inc->GetNbinsX(); ++i) { //X
		double sum_res_inc=0;
		double sum_res_inc_st=0;
		double sum_res_int=0;

        	for (size_t j=0; j<h2d_hR_inc->GetNbinsY(); ++j) { //y
			sum_res_inc+=h2d_hR_inc->GetBinContent(j,i);
			sum_res_inc_st+=h2d_hR_inc_st->GetBinContent(j,i);
			sum_res_int+=h2d_hR_int->GetBinContent(j,i);
		} //y
		
		//normalization
        	for (size_t j=0; j<h2d_hR_inc->GetNbinsY(); ++j) { //y
			if (sum_res_inc!=0) h2d_hR_inc->SetBinContent(j,i, (h2d_hR_inc->GetBinContent(j,i))/sum_res_inc); 
			if (sum_res_inc_st!=0) h2d_hR_inc_st->SetBinContent(j,i, (h2d_hR_inc_st->GetBinContent(j,i))/sum_res_inc_st); 
			if (sum_res_int!=0) h2d_hR_int->SetBinContent(j,i, (h2d_hR_int->GetBinContent(j,i))/sum_res_int); 
		} //y

	} //X
*/



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
        for (int i=1; i<=nb; ++i) { //x

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
		
		cov_inc=cov_cnp(m_inc, p_inc);
		cov_inc_st=cov_cnp(m_st_inc, p_st_inc);
		cov_int=cov_cnp(m_int, p_int);
			
		//CNP approach	
		h2d_cov_inc->Fill(i, i, cov_inc);
		h2d_cov_inc_st->Fill(i, i, cov_inc_st);
		h2d_cov_int->Fill(i, i, cov_int);

		//h2d_cov_inc->Fill(i, i, pow(data_inc_bkgfree->GetBinError(i),2));
		//h2d_cov_inc_st->Fill(i, i, pow(data_st_inc_bkgfree->GetBinError(i),2));
		//h2d_cov_int->Fill(i, i, pow(data_int_bkgfree->GetBinError(i),2));
	} //x
	//---------------------------------------------------------------------------------------------------------//

	//Output files --------------------------------------------------------//
	//TH1D htrue_signal; a model of true signal
	//TH1D hmeas; measured spectrum
	//TH2D hcov_tot; covariance matrix of measured spectrum uncertainty
	//TH2D hR; resposne matrix (X-axis: measured, Y-axis: true)

	//inc
	TFile *f_out_inc = new TFile(fout_inc.Data(),"RECREATE");
		mc_inc_true_bkgfree->Write("htrue_signal");
		data_inc->Write("hmeas_before_bkgsub");
		data_inc_bkgfree->Write("hmeas");
		h2d_hR_inc->Write("hR");
		h2d_cov_inc->Write("hcov_tot");
	f_out_inc->Close();

	//inc_st
	TFile *f_out_inc_st = new TFile(fout_inc_st.Data(),"RECREATE");
		mc_st_inc_true_bkgfree->Write("htrue_signal");
		data_st_inc->Write("hmeas_before_bkgsub");
		data_st_inc_bkgfree->Write("hmeas");
		h2d_hR_st_inc->Write("hR");
		h2d_cov_inc_st->Write("hcov_tot");
	f_out_inc_st->Close();

	//int
	TFile *f_out_int = new TFile(fout_int.Data(),"RECREATE");
		mc_int_true_bkgfree->Write("htrue_signal");
		data_int_bkgfree->Write("hmeas");
		h2d_hR_int->Write("hR");
		h2d_cov_int->Write("hcov_tot");
		data_int->Write("hmeas_before_bkgsub");
	f_out_int->Close();
	//---------------------------------------------------------------------//

	//Output files2 -------------------------------------------------------------------------------------------------//
	TString str_uf_inc=Form("./Wiener-SVD-Unfolding/Example %s %s 2 0", fout_inc.Data(), fout2_inc.Data());
	TString str_uf_inc_st=Form("./Wiener-SVD-Unfolding/Example %s %s 2 0", fout_inc_st.Data(), fout2_inc_st.Data());
	TString str_uf_int=Form("./Wiener-SVD-Unfolding/Example %s %s 2 0", fout_int.Data(), fout2_int.Data());

	std::cout<<"\n"<<std::endl;
	std::cout<<str_uf_inc<<std::endl;
	gSystem->Exec(str_uf_inc.Data());

	std::cout<<str_uf_inc_st<<std::endl;
	gSystem->Exec(str_uf_inc_st.Data());

	std::cout<<str_uf_int<<std::endl;
	gSystem->Exec(str_uf_int.Data());
	//----------------------------------------------------------------------------------------------------------------//

/*
INSTRUCTION

Input:
    true/expected signal spectrum (models);
    measured spectrum;
    covariance matrix of the uncertainty of the measured spectrum;
    response matrix from signal (corresponding to matrix column, TH2D Y-axis) to measurement (corresponding to matrix row, TH2D X-axis), M = R*S;
    Choice of addional matrix for smoothness, e.g. 2nd derivative matrix.

Output:
    Unfolded spectrum
    Additional smearing matrix
    Wiener filter
    Covariance matrix of the unfolded spectrum

An example of application can be found in Example.C
Usage:
make clean;
make;
Example input.root output.root 2 0 (2 means: 2nd derivative matrix; 0 means the number of measured events in each bin normalized by m(i)^0)
*/


}
