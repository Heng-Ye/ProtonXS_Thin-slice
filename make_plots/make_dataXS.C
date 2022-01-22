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

#include "../headers/SliceParams.h"

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

void make_dataXS() {

	TString outpath="./plots_XS/";

	//data & mc
        //TString fmc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/prod4a_thinslice_dx5cm_20slcs_new.root"; //mc, .5*test+.5*valid, no_bmrw
        //TString fmc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/prod4areco2_data_thinslice_dx5cm_20slcs.root";
        //TString fmc_nobmrw="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/prod4areco2_data_thinslice_dx5cm_20slcs_nobmrw.root";
        //TString fdata="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_run5387_prod4a_thinslice_dx5cm_20slcs.root";

	//fake data test
        //TString fmc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/prod4areco2_mc_thinslice_dx5cm_20slcs.root";
        //TString fdata="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/prod4areco2_mc_thinslice_dx5cm_20slcs.root";

        TString fmc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/prod4areco2_mc_thinslice_dx5cm_20slcs_nobmrw.root";
        TString fdata="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/prod4areco2_mc_thinslice_dx5cm_20slcs_nobmrw.root";

        //TString fmc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/prod4areco2_data_thinslice_dx4cm_25slcs.root";
        //TString fdata="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_run5387_prod4a_thinslice_dx4cm_25slcs.root";

	TString str_inc=Form("h_recosliceid_allevts_cuts");
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
	//TH1D *data_int=(TH1D*)f_data->Get("h_recosliceid_inelastic_cuts"); //h_recosliceid_inelastic_cuts
	TH1D *data_int=(TH1D*)f_data->Get(str_int.Data()); //h_recosliceid_inelastic_cuts
	data_inc->SetName("data_inc");	
	data_int->SetName("data_int");	

	int n_data_inc=data_inc->Integral();
	int n_data_int=data_int->Integral();


	//read mc [after bmrw] ------------------------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fmc.Data());

	//get truth inc & int
	TH1D *mc_true_incidents=(TH1D *)f_mc->Get("h_true_incidents");
	TH1D *mc_true_interactions=(TH1D *)f_mc->Get("h_true_interactions");

	//get truesliceID of int &inc
	TH1D *mc_truesliceID_inel=(TH1D *)f_mc->Get("h_truesliceid_inelastic_all"); //IsPureInEL
	TH1D *mc_truesliceID_all=(TH1D *)f_mc->Get("h_truesliceid_all"); //all protons


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
	
	mc_inc_inel->SetFillColor(2); mc_inc_inel->SetLineColor(2);
	mc_inc_el->SetFillColor(4); mc_inc_el->SetLineColor(4);
	mc_inc_midp->SetFillColor(3); mc_inc_midp->SetLineColor(3);
	mc_inc_midcosmic->SetFillColor(5); mc_inc_midcosmic->SetLineColor(5);
	mc_inc_midpi->SetFillColor(6); mc_inc_midpi->SetLineColor(6);
	mc_inc_midmu->SetFillColor(28); mc_inc_midmu->SetLineColor(28);
	mc_inc_mideg->SetFillColor(30); mc_inc_mideg->SetLineColor(30);
	mc_inc_midother->SetFillColor(15); mc_inc_midother->SetLineColor(15);

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

	mc_inc_inel->Scale(norm_mc_inc);
	mc_inc_el->Scale(norm_mc_inc);
	mc_inc_midcosmic->Scale(norm_mc_inc);
	mc_inc_midpi->Scale(norm_mc_inc);
	mc_inc_midp->Scale(norm_mc_inc);
	mc_inc_midmu->Scale(norm_mc_inc);
	mc_inc_mideg->Scale(norm_mc_inc);
	mc_inc_midother->Scale(norm_mc_inc);

	THStack* hs_inc=new THStack("hs_inc","");
	hs_inc->Add(mc_inc_inel);
	hs_inc->Add(mc_inc_el);
	hs_inc->Add(mc_inc_midcosmic);
	hs_inc->Add(mc_inc_midpi);
	hs_inc->Add(mc_inc_midp);
	hs_inc->Add(mc_inc_midmu);
	hs_inc->Add(mc_inc_mideg);
	hs_inc->Add(mc_inc_midother);

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
	mc_int_el->SetFillColor(4); mc_int_el->SetLineColor(4);
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

	TH1D* pur_int=(TH1D* )f_mc->Get("pur_Int"); //int
	TH1D *eff_int=(TH1D* )f_mc->Get("eff_Int"); //int

	//bkg subtraction --------------------------------------------------------------------------//
	TH1D* data_inc_bkgfree=(TH1D *)data_inc->Clone("data_inc_bkgfree"); data_inc_bkgfree->SetName("data_inc_bkgfree");	
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
	RooUnfoldResponse *res_inc=(RooUnfoldResponse*)f_mc->Get("response_SliceID_Inc"); res_inc->SetName("res_inc"); //inc
	RooUnfoldBayes uf_inc (res_inc, data_inc_bkgfree, 4); //inc

	//sansity check on the response matrix
	std::cout<<"res_inc->GetDimensionMeasured():"<<res_inc->GetDimensionMeasured()<<std::endl;
	std::cout<<"res_inc->GetDimensionTruth():"<<res_inc->GetDimensionTruth()<<std::endl;
	std::cout<<"res_inc->GetNbinsMeasured():"<<res_inc->GetNbinsMeasured()<<std::endl;
	std::cout<<"res_inc->GetNbinsTruth():"<<res_inc->GetNbinsTruth()<<std::endl; 

	//std::cout<<"GetIterations:"<<uf_inc.GetIterations()<<std::endl;
	//std::cout<<"GetSmoothing:"<<uf_inc.GetSmoothing()<<std::endl;
	//uf_inc.Print();

	RooUnfoldResponse *res_int=(RooUnfoldResponse*)f_mc->Get("response_SliceID_Int"); res_int->SetName("res_int"); //int
	RooUnfoldBayes uf_int (res_int, data_int_bkgfree, 4);

	//std::cout<<"data_inc_bkgfree x-axis bin:"<<data_inc_bkgfree->GetNbinsX()<<std::endl;
	//std::cout<<"data_int_bkgfree x-axis bin:"<<data_int_bkgfree->GetNbinsX()<<std::endl;
	//std::cout<<"response_SliceID_Inc x-axis bin:"<<(TH2D *)f_mc->Get("response_SliceID_Int")->GetNbinsX()<<std::endl;
	
	//unfolding
  	TH1D *data_inc_uf;
  	TH1D *data_int_uf;
	data_inc_uf=(TH1D* )uf_inc.Hreco();
  	data_int_uf=(TH1D *)uf_int.Hreco();
	data_inc_uf->SetNameTitle("data_inc_uf", "Unfolded incident protons; Slice ID; Events");
	data_int_uf->SetNameTitle("data_int_uf", "Unfolded interacting protons; Slice ID; Events");
	//unfold.IncludeSystematics(2); // Default 1, propagates both statistical+systematics, =2 no errors included.


	//get KEs from MC --------------------------------------------------------//
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
		//cout<<"avg_recoincE["<<i<<"]:"<<avg_recoincE[i]<<endl;
        }
	
	//Calc XS ----------------------------------------------------------------------------------------------------------------------------//
	double xs_const=MAr/(Density*NA*thinslicewidth)*1e27;
	//[1] true xs
	double true_xs[nthinslices] = {0};
	double err_true_xs[nthinslices] = {0};
	double true_inc[nthinslices] = {0};
	double true_int[nthinslices] = {0};
        for (int i = 0; i<nthinslices; ++i) {
		true_inc[i]=mc_true_incidents->GetBinContent(i+1);
		true_int[i]=mc_true_interactions->GetBinContent(i+1);
		if (true_inc[i]&&true_int[i]) {
			true_xs[i]=xs_const*log(true_inc[i]/(true_inc[i]-true_int[i]));
			err_true_xs[i]=xs_const*sqrt(true_int[i]+pow(true_int[i],2)/true_inc[i])/true_inc[i];
		}
        }

	//[2] reco xs 
	//inc & int
	double sliceid[nthinslices] = {0};
	double reco_inc[nthinslices] = {0};
	double reco_int[nthinslices] = {0};
	double err_reco_inc[nthinslices] = {0};
	double err_reco_int[nthinslices] = {0};
	double err_sliceid[nthinslices] = {0};

        for (int i = 0; i<nthinslices; ++i) {
		sliceid[i]=i+.5;

		//int[inelastic]
		reco_int[i]=data_int_uf->GetBinContent(i+2);
		err_reco_int[i]=data_int_uf->GetBinError(i+2);

		//inc[all]
    		for (int j=i; j<=nthinslices; ++j){
      			reco_inc[i] += data_inc_uf->GetBinContent(j+2);
      			err_reco_inc[i] += pow(data_inc_uf->GetBinError(j+2),2);
    		}
		err_reco_inc[i]=sqrt(err_reco_inc[i]);
	}

	//for (int ij = 0; ij<=true_sliceID; ++ij){
		//if (ij<nthinslices) ++true_incidents[ij];
	//}

	double reco_xs[nthinslices] = {0};
	double err_reco_xs[nthinslices] = {0};
        for (int i = 0; i<nthinslices; ++i) {
		if (reco_inc[i]&&reco_int[i]) {
    			reco_xs[i]=xs_const*log(reco_inc[i]/(reco_inc[i]-reco_int[i]));
    			err_reco_xs[i]=xs_const*sqrt(pow(reco_int[i]*err_reco_inc[i]/reco_inc[i]/(reco_inc[i]-reco_int[i]),2)+pow(err_reco_int[i]/(reco_inc[i]-reco_int[i]),2));
		}
        }
  	TGraphErrors *gr_truexs = new TGraphErrors(nthinslices, avg_trueincE, true_xs, err_trueincE, err_true_xs);
  	TGraphErrors *gr_recoxs = new TGraphErrors(nthinslices, avg_trueincE, reco_xs, err_trueincE, err_reco_xs);

  	TGraphErrors *gr_reco_inc = new TGraphErrors(nthinslices, sliceid, reco_inc, err_sliceid, err_reco_inc);
	gr_reco_inc->SetNameTitle("gr_reco_inc", "; SliceID; Events");

  	gr_truexs->SetNameTitle("gr_truexs", "; Energy (MeV); Cross-section [mb]");
  	gr_recoxs->SetNameTitle("gr_recoxs", "; Energy (MeV); Cross-section [mb]");

	//All protons ---------------------------------------------------------------------//
	//[1]data with mc components
	TCanvas *c_inc=new TCanvas(Form("c_inc"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_inc->Divide(1,1);
	c_inc->cd(1);
	float ymax_evt=6000;
	TH2D *f2d_inc=new TH2D("f2d_inc",Form("%s","All Protons"),21,0,21,80,0,ymax_evt);
	f2d_inc->GetXaxis()->SetTitle("Reco SliceID (End point)");
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
	c_inc->Print(Form("./plot_dataXS/data_recosliceid_inc.eps"));

	//[2]bkg subtraction
	TCanvas *c_inc2=new TCanvas(Form("c_inc2"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_inc2->Divide(1,1);
	c_inc2->cd(1);
	TH2D *f2d_inc2=new TH2D("f2d_inc2",Form("%s","INC"),21,0,21,80,0,ymax_evt);
	f2d_inc2->SetTitle("All protons; Reco SliceID (End point); Events");
	f2d_inc2->Draw();

	THStack* hs_inc_sg=new THStack("hs_inc_sg","");
	hs_inc_sg->Add(mc_inc_inel);
	hs_inc_sg->Add(mc_inc_el);
	hs_inc_sg->Draw("hist same");

	data_inc->Draw("ep same");
	data_inc_bkgfree->SetMarkerColor(7);
	data_inc_bkgfree->SetLineColor(7);
	data_inc_bkgfree->SetMarkerStyle(24);
	data_inc_bkgfree->Draw("ep same");

	TLegend *leg_inc2 = new TLegend(0.2,0.6,0.8,0.9);
	leg_inc2->SetFillStyle(0);
	leg_inc2->AddEntry(data_inc, "Selected", "ep");
	leg_inc2->AddEntry(data_inc_bkgfree, "Selected+BKG Subtraction","ep");
	leg_inc2->AddEntry(mc_inc_inel, "Selected MC Truth [Inel.]","f");
	leg_inc2->AddEntry(mc_inc_el, "Selected MC Truth [El.]","f");
	//leg_inc2->SetNColumns(3);
	leg_inc2->Draw();
	c_inc2->Print(Form("./plot_dataXS/data_recosliceid_bkgsub_inc.eps"));
	
	
	//[3]data unfolding (true sliceID)
	TCanvas *c_inc3=new TCanvas(Form("c_inc3"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_inc3->Divide(1,1);
	c_inc3->cd(1);
	TH2D *f2d_inc3=new TH2D("f2d_inc3",Form("%s","All protons"),22,-1,21,100,0,10000);
	f2d_inc3->SetTitle("All protons; True SliceID (End point); Events");
	f2d_inc3->Draw();
	data_inc_bkgfree->Draw("ep same");

	data_inc_uf->SetLineColor(4);
	data_inc_uf->SetMarkerColor(4);
        data_inc_uf->SetMarkerStyle(21);
	data_inc_uf->Draw("ep same");

	mc_truesliceID_all->SetLineColor(2);
	mc_truesliceID_all->Draw("hist same");

	//TH1D* data_inc_bkgfree_eff=(TH1D *)data_inc_bkgfree->Clone();
	//data_inc_bkgfree_eff->Divide(eff_inc);
	//data_inc_bkgfree_eff->Draw("hist same");

	TLegend *leg_inc3 = new TLegend(0.2,0.6,0.8,0.9);
	leg_inc3->SetFillStyle(0);
	leg_inc3->AddEntry(data_inc_bkgfree, "Selected Protons after BKG Subtraction","ep");
	leg_inc3->AddEntry(data_inc_uf, "Unfolding","ep");
	leg_inc3->AddEntry(mc_truesliceID_all,"True Protons","l");
	leg_inc3->Draw();
	c_inc3->Print(Form("./plot_dataXS/data_recosliceid_uf_inc.eps"));


	//[4] eff & purity
	TCanvas *c_inc4=new TCanvas(Form("c_inc4"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_inc4->Divide(1,1);
	c_inc4->cd(1);
	TH2D *f2d_inc4=new TH2D("f2d_inc4",Form("%s",""),22,-1,21,10,0,1.2);
	f2d_inc4->SetTitle("All protons; Reco SliceID; Efficiency");
	f2d_inc4->Draw();
	eff_inc->Draw("same");
	pur_inc->SetMarkerColor(2);
	pur_inc->SetLineColor(2);
	pur_inc->Draw("same");
	c_inc4->Print(Form("./plot_dataXS/data_eff_pur_inc.eps"));


	//[5]Ninc distribution
	TCanvas *c_inc5=new TCanvas(Form("c_inc5"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_inc5->Divide(1,1);
	c_inc5->cd(1);
	TH2D *f2d_inc5=new TH2D("f2d_inc5",Form("%s",""),22,-1,21,10,0,150000);
	f2d_inc5->SetTitle("All protons; True SliceID; Events");
	f2d_inc5->Draw();
	mc_true_incidents->SetLineColor(2);
	mc_true_incidents->Draw("same");
	gr_reco_inc->Draw("ep same");
	c_inc5->Print(Form("./plot_dataXS/incident_inc.eps"));

	
	//inelastic scattering protons ------------------------------------------------//
	//[1]data with mc components
	TCanvas *c_int=new TCanvas(Form("c_int"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_int->Divide(1,1);
	c_int->cd(1);
	TH2D *f2d_int=new TH2D("f2d_int",Form("%s","Proton Inelastic Scatterings"),21,0,21,80,0,ymax_evt);
	f2d_int->GetXaxis()->SetTitle("Reco SliceID");
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
	c_int->Print(Form("./plot_dataXS/fakedata_recosliceid_int.eps"));


	//[2]data with bkg subtraction
	TCanvas *c_int2=new TCanvas(Form("c_int2"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_int2->Divide(1,1);
	c_int2->cd(1);
	TH2D *f2d_int2=new TH2D("f2d_int2",Form("%s","INT"),21,0,21,80,0,ymax_evt);
	f2d_int2->SetTitle("Proton Inelasic Scatterings; Reco SliceID; Events");
	f2d_int2->Draw();
	mc_int_inel->SetLineColor(2); //mc 
	mc_int_inel->Draw("hist same");
	data_int->Draw("ep same");
	data_int_bkgfree->SetMarkerColor(3);
	data_int_bkgfree->SetLineColor(3);
	data_int_bkgfree->SetMarkerStyle(24);
	data_int_bkgfree->Draw("ep same");

	TLegend *leg_int2 = new TLegend(0.2,0.6,0.8,0.9);
	leg_int2->SetFillStyle(0);
	leg_int2->AddEntry(data_int, "Selected", "ep");
	leg_int2->AddEntry(data_int_bkgfree, "Selected+BKG Subtraction","ep");
	leg_int2->AddEntry(mc_int_inel, "Selected MC Truth","f");
	leg_int2->Draw();

	c_int2->Print(Form("./plot_dataXS/data_recosliceid_bkgsub_int.eps"));


	//[3]data unfolding (true sliceID)
	TCanvas *c_int3=new TCanvas(Form("c_int3"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_int3->Divide(1,1);
	c_int3->cd(1);
	TH2D *f2d_int3=new TH2D("f2d_int3",Form("%s","INT"),22,-1,21,100,0,10000);
	f2d_int3->SetTitle("Proton Inelasic Scatterings; True SliceID; Events");
	f2d_int3->Draw();
	data_int_bkgfree->Draw("ep same");

	data_int_uf->SetLineColor(4);
	data_int_uf->SetMarkerColor(4);
        data_int_uf->SetMarkerStyle(21);
	data_int_uf->Draw("ep same");

	mc_truesliceID_inel->SetLineColor(2);
	mc_truesliceID_inel->Draw("hist same");

	//mc_true_interactions->Draw("hist same");

	//TH1D* data_int_bkgfree_eff=(TH1D *)data_int_bkgfree->Clone();
	//data_int_bkgfree_eff->Divide(eff_int);
	//data_int_bkgfree_eff->Draw("hist same");


	TLegend *leg_int3 = new TLegend(0.2,0.6,0.8,0.9);
	leg_int3->SetFillStyle(0);
	leg_int3->AddEntry(data_int_bkgfree, "Selected Inelastic Protons after BKG Subtraction","ep");
	leg_int3->AddEntry(data_int_uf, "Unfolding","ep");
	leg_int3->AddEntry(mc_truesliceID_inel,"True Inelastic Scattering Protons","l");
	//leg_int2->SetNColumns(3);
	leg_int3->Draw();

	c_int3->Print(Form("./plot_dataXS/data_recosliceid_uf_int.eps"));


	//[4] eff & purity
	TCanvas *c_int4=new TCanvas(Form("c_int4"),"",900, 600);
	gStyle->SetTitleX(0.5); 
	gStyle->SetTitleAlign(23); 

	c_int4->Divide(1,1);
	c_int4->cd(1);
	TH2D *f2d_int4=new TH2D("f2d_int4",Form("%s","INT"),22,-1,21,10,0,1.2);
	f2d_int4->SetTitle("Proton Inelasic Scatterings; Reco SliceID; Efficiency");
	f2d_int4->Draw();
	eff_int->Draw("same");
	c_int4->Update();

	pur_int->SetStats(kFALSE);
  	Float_t rightmax_pur_int=1.2*pur_int->GetMaximum();
	//Double_t min_pur_int=pur_int->GetMinimum();
  	Float_t scale_pur_int=gPad->GetUymax()/rightmax_pur_int;
  	pur_int->Scale(scale_pur_int);
	pur_int->SetMarkerColor(2);
	pur_int->SetLineColor(2);
	pur_int->Draw("same");

  	TGaxis *ax_pur_int = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
                   gPad->GetUxmax(), gPad->GetUymax(),
                   0, rightmax_pur_int, 510, "+L");
	ax_pur_int->SetLineColor(kRed);
   	ax_pur_int->SetLabelColor(kRed);
  	ax_pur_int->Draw("");
  	//
	c_int4->Print(Form("./plot_dataXS/data_eff_pur_int.eps"));











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

        float ymax=1500;
        float xmax=450;
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

        TLegend *leg_xs = new TLegend(0.6,0.6,0.9,0.85);
        leg_xs->SetFillStyle(0);
        leg_xs->AddEntry(total_inel_KE, "Geant4", "l");
        leg_xs->AddEntry(gr_truexs, "MC Truth", "pe");
        leg_xs->AddEntry(gr_recoxs, "MC Reco", "pe");
        leg_xs->Draw();
	c_xs->Print(Form("./plot_dataXS/xs_data.eps"));






	


}
