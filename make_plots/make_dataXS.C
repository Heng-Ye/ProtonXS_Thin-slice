#include "TGraphErrors.h"
#include "TVector3.h"
//#include "RooUnfoldBayes.h"

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

#include "../headers/SliceParams.h"




//constants ---------------------//
double NA=6.02214076e23;
double MAr=39.95; //gmol
double Density = 1.4; // g/cm^3
double true_cosangle = 1.;
//-------------------------------//

using namespace std;
using namespace ROOT::Math;


//R__LOAD_LIBRARY(/dune/app/users/hyliao/WORK/larsoft_mydev/RooUnfold/build/libRooUnfold.so) //load share lib
//#include "../headers/Unfold.h"
//#include "../headers/RooUnfold.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

void make_dataXS() {

	TString outpath="./plots_XS/";
        TString fmc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/prod4areco2_data_thinslice_dx5cm_20slcs.root";
        TString fdata="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_run5387_prod4a_thinslice_dx5cm_20slcs.root";

	TString str_inc=Form("h_recosliceid_allevts_cuts");
	TString str_int=Form("h_recosliceid_inelastic_cuts");




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
	TH1D *data_int=(TH1D*)f_data->Get(str_int.Data()); //h_recosliceid_inelastic_cuts
	data_inc->SetName("data_inc");	
	data_int->SetName("data_int");	

	//read mc [after bmrw] ------------------------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fmc.Data());

	//response matrix
	RooUnfoldResponse *response_SliceID_Inc=(RooUnfoldResponse*)f_mc->Get("response_SliceID_Inc"); //inc
	RooUnfoldResponse *response_SliceID_Int=(RooUnfoldResponse*)f_mc->Get("response_SliceID_Int"); //int

	//get mc reco slice IDs
	//inc
	TH1D* mc_inc_all=(TH1D*)f_mc->Get(Form("%s",str_inc.Data()));
	TH1D* mc_inc_inel=(TH1D*)f_mc->Get(Form("%s_inel",str_inc.Data()));
	TH1D* mc_inc_el=(TH1D*)f_mc->Get(Form("%s_el",str_inc.Data()));
	TH1D* mc_inc_midcosmic=(TH1D*)f_mc->Get(Form("%s_midcosmic",str_inc.Data()));
	
h_recosliceid_allevts_cuts_midcosmic
h_recosliceid_allevts_cuts_midpi
h_recosliceid_allevts_cuts_midp
h_recosliceid_allevts_cuts_midmu
h_recosliceid_allevts_cuts_mideg
h_recosliceid_allevts_cuts_midother


h_recosliceid_recoinelastic_cuts_inel
h_recosliceid_recoinelastic_cuts_el
h_recosliceid_recoinelastic_cuts_midcosmic
h_recosliceid_recoinelastic_cuts_midpi
h_recosliceid_recoinelastic_cuts_midp
h_recosliceid_recoinelastic_cuts_midmu
h_recosliceid_recoinelastic_cuts_mideg
h_recosliceid_recoinelastic_cuts_midother

	//purity & efficiency
	TH1D* pur_inc=(TH1D* )f_mc->Get("pur_Inc"); //inc
	TH1D *eff_inc=(TH1D* )f_mc->Get("eff_Inc"); //inc

	TH1D* pur_int=(TH1D* )f_mc->Get("pur_Int"); //int
	TH1D *eff_int=(TH1D* )f_mc->Get("eff_Int"); //int

	//get KEs
	double slcid[nthinslices] = {0};
	TH1D *reco_incE[nthinslices];
	TH1D *true_incE[nthinslices];
	double avg_recoincE[nthinslices] = {0};
	double avg_trueincE[nthinslices] = {0};
	double err_recoincE[nthinslices] = {0};
	double err_trueincE[nthinslices] = {0};
        for (int i = 0; i<nthinslices; ++i){
                slcid[i] = i;
		reco_incE[i]=(TH1D* )f_mc->Get(Form("reco_incE_%d",i));
		true_incE[i]=(TH1D* )f_mc->Get(Form("true_incE_%d",i));

                avg_trueincE[i] = true_incE[i]->GetMean();
                err_trueincE[i] = true_incE[i]->GetMeanError();
                avg_recoincE[i] = reco_incE[i]->GetMean();
                err_recoincE[i] = reco_incE[i]->GetMeanError();
        }
	


/*
	//recosliceID
	//inc
	TH1D* mc_inc=(TH1D *)f_mc->Get(Form("%s", str_inc.Data()));


	TH1D* mc_inc_inel=(TH1D *)f_mc->Get(Form("%s", str_inc.Data()));

h_recosliceid_allevts_cuts
h_recosliceid_allevts_cuts_inel
h_recosliceid_allevts_cuts_el
h_recosliceid_allevts_cuts_midcosmic
h_recosliceid_allevts_cuts_midpi
h_recosliceid_allevts_cuts_midp
h_recosliceid_allevts_cuts_midmu
h_recosliceid_allevts_cuts_mideg
h_recosliceid_allevts_cuts_midother

	//int
h_recosliceid_recoinelastic_cuts
h_recosliceid_recoinelastic_cuts_inel
h_recosliceid_recoinelastic_cuts_el
h_recosliceid_recoinelastic_cuts_midcosmic
h_recosliceid_recoinelastic_cuts_midpi
h_recosliceid_recoinelastic_cuts_midp
h_recosliceid_recoinelastic_cuts_midmu
h_recosliceid_recoinelastic_cuts_mideg
h_recosliceid_recoinelastic_cuts_midother


TH2D *h2d=(TH2D *)f_mc->Get(Form("h2d_trklen_ke2%s",str_obs.Data()));
*/	


	//inc and int after purity correction
	TH1D *data_inc_pur=(TH1D *)data_inc->Clone(); 
	TH1D *data_int_pur=(TH1D *)data_int->Clone(); 
	data_inc_pur->SetName("data_inc_pur");	
	data_int_pur->SetName("data_int_pur");	

        data_inc_pur->Multiply(pur_inc);
        data_int_pur->Multiply(pur_int);


	//recosliceID after unfolding
	//RooUnfoldResponse response_SliceID_Inc;

	//virtual RooUnfoldResponse& Setup (const TH1* measured, const TH1* truth, const TH2* response);  // set up from already-filled histograms
	//RooUnfoldBayes   inc_uf (&response_SliceID_Inc, data_inc_pur, 4);    // OR
	//Unfold uf(nthinslices+2, -1, nthinslices+1);
	//uf.response_SliceID_Inc.Clone("response_SliceID_Inc");
	//uf.response_SliceID_Inc.Hresponse()=response_SliceID_Inc;


  	//RooUnfoldBayes   unfold (&response, hMeas, 4);    // OR
  	//TH1D* hReco= (TH1D*) unfold.Hreco();



	//TCanvas *c_ = new TCanvas("c_", "c_", 900, 600);
	//response_SliceID_Inc->Draw("colz");
	//uf.response_SliceID_Inc.Draw("colz");
	//uf.response_SliceID_Inc=response_SliceID_Inc;
	//uf.response_SliceID_Int=response_SliceID_Int;
	//data unfolding
        //RooUnfoldBayes   unfold_Inc (response_SliceID_Inc, data_inc_pur, 4);
        //RooUnfoldBayes   unfold_Int (&response_SliceID_Int, data_int_pur, 4);


	//h_truesliceid_uf = (TH1D*) unfold_Inc.Hreco();

//TH2D *hint = (TH2D*)response_SliceID_Int.Hresponse()

/*
	//xs
        for (int i = 0; i<nthinslices; ++i){
                //std::cout<<i<<" "<<avg_trueincE[i]<<std::endl;
                if (true_incidents[i] && true_interactions[i]){
                        //truexs[i] = MAr/(Density*NA*thinslicewidth/true_AngCorr->GetMean())*log(true_incidents[i]/(true_incidents[i]-true_interactions[i]))*1e27;
                        //err_truexs[i] = MAr/(Density*NA*thinslicewidth/true_AngCorr->GetMean())*1e27*sqrt(true_interactions[i]+pow(true_interactions[i],2)/true_incidents[i])/true_incidents[i];
                        truexs[i] = MAr/(Density*NA*thinslicewidth/true_cosangle)*log(true_incidents[i]/(true_incidents[i]-true_interactions[i]))*1e27;
                        err_truexs[i] = MAr/(Density*NA*thinslicewidth/true_cosangle)*1e27*sqrt(true_interactions[i]+pow(true_interactions[i],2)/true_incidents[i])/true_incidents[i];
                }
        }

        TGraphErrors *gr_trueincE = new TGraphErrors(nthinslices, &(slcid[0]), &(avg_trueincE[0]), 0, &(err_trueincE[0]));
        TGraphErrors *gr_recoincE = new TGraphErrors(nthinslices, &(slcid[0]), &(avg_recoincE[0]), 0, &(err_recoincE[0]));
        TGraphErrors *gr_reco_trueincE = new TGraphErrors(nthinslices, &(slcid[0]), &(reco_trueincE[0]), 0, &(err_reco_trueincE[0]));

        gr_trueincE->Write("gr_trueincE");
        gr_recoincE->Write("gr_recoincE");
        gr_reco_trueincE->Write("gr_reco_trueincE");




        TH1D *hinc = (TH1D*)h_recosliceid_allevts_cuts->Clone("hinc");
        //TH1D *hint = (TH1D*)h_recosliceid_allevts_cuts->Clone("hint");
        TH1D *hint = (TH1D*)h_recosliceid_recoinelastic_cuts->Clone("hint");
        hinc->Multiply(uf.pur_Inc);
        hint->Multiply(uf.pur_Int);

        //  RooUnfoldBayes   unfold_Inc (&uf.response_SliceID_Inc, uf.pur_num_Inc, 4);
        //  RooUnfoldBayes   unfold_Int (&uf.response_SliceID_Int, uf.pur_num_Int, 4);

        RooUnfoldBayes   unfold_Inc (&uf.response_SliceID_Inc, hinc, 4);
        RooUnfoldBayes   unfold_Int (&uf.response_SliceID_Int, hint, 4);

        //RooUnfoldBayes   unfold_Inc (&uf.response_SliceID_Inc, hinc, 12);
        //RooUnfoldBayes   unfold_Int (&uf.response_SliceID_Int, hint, 12);

        //RooUnfoldSvd     unfold_Inc (&uf.response_SliceID_Inc, hinc, 20);   // OR
        //RooUnfoldSvd     unfold_Int (&uf.response_SliceID_Int, hint, 20);   // OR

        //  RooUnfoldSvd     unfold_Inc (&uf.response_SliceID_Inc, uf.pur_num_Inc, 20);   // OR
        //  RooUnfoldSvd     unfold_Int (&uf.response_SliceID_Int, uf.pur_num_Int, 20);   // OR

        h_truesliceid_uf = (TH1D*) unfold_Inc.Hreco();
        h_truesliceid_inelastic_uf = (TH1D*) unfold_Int.Hreco();



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



*/











	//plot slice ID with mc thuth labels
	//inc
	//int
	
	//xs
	


}
