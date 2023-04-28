#define ProtonESliceDataAll_cxx
#include "ProtonESliceDataAll.h"

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

#include "./cali/dedx_function_35ms.h"
#include "./headers/BasicParameters.h"
#include "./headers/BasicFunctions.h"
#include "./headers/ESliceParams.h"
#include "./headers/util.h"
#include "./headers/ESlice.h"
#include "./headers/BetheBloch.h"

#include <chrono>
#include <ratio>
#include <thread>

using namespace std;
using namespace ROOT::Math;

/////////////////////////////////
// define the parametric line equation
void line(double t, const double *p, double &x, double &y, double &z) {
	// a parametric line is define from 6 parameters but 4 are independent
	// x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
	// can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
	x = p[0] + p[1]*t;
	y = p[2] + p[3]*t;
	z = t;
}

bool first = true;

// function Object to be minimized
struct SumDistance2 {
	// the TGraph is a data member of the object
	TGraph2D *fGraph;

	SumDistance2(TGraph2D *g) : fGraph(g) {}

	// calculate distance line-point
	double distance2(double x,double y,double z, const double *p) {
		// distance line point is D= | (xp-x0) cross  ux |
		// where ux is direction of line and x0 is a point in the line (like t = 0)
		XYZVector xp(x,y,z);
		XYZVector x0(p[0], p[2], 0. );
		XYZVector x1(p[0] + p[1], p[2] + p[3], 1. );
		XYZVector u = (x1-x0).Unit();
		double d2 = ((xp-x0).Cross(u)).Mag2();
		return d2;
	}

	// implementation of the function to be minimized
	double operator() (const double *par) {
		assert(fGraph != 0);
		double * x = fGraph->GetX();
		double * y = fGraph->GetY();
		double * z = fGraph->GetZ();
		int npoints = fGraph->GetN();
		double sum = 0;
		for (int i  = 0; i < npoints; ++i) {
			double d = distance2(x[i],y[i],z[i],par);
			sum += d;
		}
		if (first) {
			//std::cout << "Total Initial distance square = " << sum << std::endl;
		}
		first = false;
		return sum;
	}

};
/////////////////////////////////


void ProtonESliceDataAll::Loop() {
	if (fChain == 0) return;

	//time_start
	auto t_st = std::chrono::high_resolution_clock::now();	

	//MC beam momentum -----------------------//
	double mm1=1007.1482; //MC prod4a [spec]
	double ss1=60.703307; //MC prod4a [spec]
	double mu_min=mm1-3.*ss1;
	double mu_max=mm1+3.*ss1;

	//Beam momentum reweighting -----------------------------//
	//MC KE beam Gaussian 
	//double m1=3.89270e+02; //mc, keff with const E-loss
	//double s1=4.49638e+01; //mc, keff with const E-loss
	//double a1=7.06341e+02; //mc, keff with const E-loss

	//double m2=3.93027e+02; //data, keff with const E-loss
	//double s2=5.18623e+01; //data, keff with const E-loss
	//double a2=6.09665e+02; //data, keff with const E-loss

	//double xmin=0.; //emin [MeV]
	//double xmax=1000.; //emax [MeV]

	//double mu_min=m1-3.*s1;
	//double mu_max=m1+3.*s1;

	//TF1 *gng=new TF1(Form("gng"),agovg,xmin,xmax,6);
	//gng->SetParameter(0,m1);
	//gng->SetParameter(1,s1);
	//gng->SetParameter(2,a1);

	//gng->SetParameter(3,m2);
	//gng->SetParameter(4,s2);
	//gng->SetParameter(5,a2);
	//-----------------------------------------------------------------//

	//[MC] Energy loss using stopping protons
	double mean_Eloss_upstream=19.3073;
	double err_mean_Eloss_upstream=0.187143;
	double sigma_Eloss_upstream=18.7378;
	double err_sigma_Eloss_upstream=0.140183;

	//const E-loss using stopping protons ---------
	double Eloss_mc_hy_stop=19.542/0.998495;
	//p[0]:19.542;   //err_p[0]:0.126113
	//p[1]:0.998495; //err_p[1]:0.00549534

	double R_fit_hy=1.0008142352819318;
	double er_R_fit_hy=0.04629667706788889;

	//const. E-loss assumption ---------------------------------------------------------------------------------------	
	double const_eloss_mc=47.0058/1.00097; //const E-loss from fit (calo)
	//#p[0]:47.0058 err_p[0]:0.372157 p[1]:-1.00097 err_p[1]:0.00787403

	//New weighting func (using KEff_fit_stop at TPC FF as a reference) ------//
	//double mu_denom_data=411.06602388610895; //old
	//double sg_denom_data=47.075678784947826; //old

	double mu_denom_data=411.05145837595467; //new with event-by-event R corr (const E-loss using stopping protons)
	double sg_denom_data=47.48714821962207; //new with event-by-event R corr (const E-loss using stopping protons)

	//double mu_denom_data=411.06645442311424; //new with event-by-event (KEHY(fit) using stopping protons)
	//double sg_denom_data=47.076122305960645; //new with event-by-event (KEHY(fit) using stopping protons)

	//Data //
	//i= 0  m= 411.05145837595467 s= 47.48714821962207 [kebeam-dE]*R (R~1)
	//i= 1  m= 411.06645442311424 s= 47.076122305960645 [KE(Fit)]

	//double mu_nom_data=390.81237292943916; //for data (now use event-by-event correction)
	//double sg_nom_data=47.52091718691363; //for data (now use event-by-event correction)

	//double mu_nom_mc=388.560260293186; //for mc(KEbeam-const_from_calo)
	//double sg_nom_mc=43.13168235197187; //formc

	double mu_nom_mc=416.224743039812; //for mc(KEbeam-const) with R=1 (R=Ratio of KEff(Fit)/(KEbeam-dE))
	double sg_nom_mc=42.786018962508784; //

	//double mu_nom_mc=416.1620092367158; //for mc [KE(Fit)]
	//double sg_nom_mc=40.48356757740762; //

	//i= 0  m= 416.2247430398121 s= 42.786018962508784 [kebeam-dE]*1 (R=1)
	//i= 1  m= 416.1620092367158 s= 40.48356757740762 [KE(Fit)] 

	//reference of error of perv. calc
	//p0           4.16231e+02   2.81061e-01  -3.22787e-05   3.06053e-05
	//p1           4.34283e+01   1.97065e-01   1.96200e-04  -1.35452e-04

	//E-dept E-loss --------------------------------//
	//mc
	double p0_edept_stop=1.97726e+01;
	double p1_edept_stop=-1.37827e-01;
	double p2_edept_stop=3.16264e-04;

	double err_p0_edept_stop=3.30157e+01;
	double err_p1_edept_stop=1.63940e-01;
	double err_p2_edept_stop=2.00009e-04;

	//weighting func. (ke) ------------------------
	TF1 *kerw=new TF1(Form("kerw"),govg,0,800,4);
	kerw->SetParameter(0, mu_nom_mc);
	kerw->SetParameter(1, sg_nom_mc);
	kerw->SetParameter(2, mu_denom_data);
	kerw->SetParameter(3, sg_denom_data);

	//ke cut range	
	//double mu_kemin=mu_nom_mc-3.*sg_nom_mc;
	//double mu_kemax=mu_nom_mc+3.*sg_nom_mc;
	double mu_kemin=mu_nom_mc-5.*sg_nom_mc;
	double mu_kemax=mu_nom_mc+5.*sg_nom_mc;
	//double mu_kemin=mu_nom_mc-6.*sg_nom_mc;
	//double mu_kemax=mu_nom_mc+6.*sg_nom_mc;
	//------------------------------------------------------------------------//

	//weighting func from E-dept E-loss upstream ----------------------------
	//weighting func: bmrw between reco and data
	TF1 *kerw_edept=new TF1(Form("kerw_edept"),govg,0,800,4);
	kerw_edept->SetParameter(0, 4.15248e+02); //mu of keffbeam_stop_mc
	kerw_edept->SetParameter(1, 3.72018e+01); //sigma of keffbeam_stop_mc
	kerw_edept->SetParameter(2, 4.09837e+02); //mu of keffbeam_stop_data
	kerw_edept->SetParameter(3, 4.28448e+01); //sigma of keffbeam_stop_data

	//keffbeam_stop_data
	//1  p0           4.09837e+02   2.93549e-01  -6.91597e-05   2.09147e-05
	//2  p1           4.28448e+01   2.02807e-01   1.95520e-04  -1.15037e-04
	//3  p2           4.11054e+02   3.35511e+00   3.35511e+00   4.19084e-05
	//keff_stop_mc
	//1  p0           4.15916e+02   2.41791e-01  -2.65188e-06   1.04777e-06
	//2  p1           3.99959e+01   1.55280e-01   9.25748e-05  -5.88385e-05
	//3  p2           4.38820e+02   3.09842e+00   3.09842e+00   3.72791e-05
	//fit_keffbeam_stop_mc
	//1  p0           4.15248e+02   2.24060e-01   1.38021e-05  -4.12830e-06
	//2  p1           3.72018e+01   1.58229e-01   1.60618e-04  -1.49393e-04
	//3  p2           4.75150e+02   3.44224e+00   3.44224e+00   3.04656e-05

	//MisID:P Modeling (good fit for KE<300 MeV) -------------------------------------------------------------------------------------------------
	double p0_data_misidp=14489;       double err_p0_data_misidp=381.785;
	double p1_data_misidp=10.053;      double err_p1_data_misidp=0.00129971;
	double p2_data_misidp=0.00319447;  double err_p2_data_misidp=7.24379e-05;
	double p3_data_misidp=-22889.9;    double err_p3_data_misidp=30.1843;

	double p0_mc_misidp=15109.6;       double err_p0_mc_misidp=502.325;
	double p1_mc_misidp=9.99314;       double err_p1_mc_misidp=0.00785002;
	double p2_mc_misidp=0.00393017;    double err_p2_mc_misidp=0.000127054;
	double p3_mc_misidp=-21553;        double err_p3_mc_misidp=171.736;
	TF1 * kerw_misidp = new TF1("kerw_misidp","([0]*ROOT::Math::lognormal_pdf(x,[1],[2],[3]))/([4]*ROOT::Math::lognormal_pdf(x,[5],[6],[7]))",0,800);
	kerw_misidp->SetParameter(0, p0_data_misidp);
	kerw_misidp->SetParameter(1, p1_data_misidp);
	kerw_misidp->SetParameter(2, p2_data_misidp);
	kerw_misidp->SetParameter(3, p3_data_misidp);

	kerw_misidp->SetParameter(4, p0_mc_misidp);
	kerw_misidp->SetParameter(5, p1_mc_misidp);
	kerw_misidp->SetParameter(6, p2_mc_misidp);
	kerw_misidp->SetParameter(7, p3_mc_misidp);

	//unfolding config. ---------------------------------------------------------------------------------------------------//
	//Unfold uf(nthinslices+2, -1, nthinslices+1);

	//2d unfolding

	TH2D* h2d_reco = new TH2D("h2d_reco","h2d_reco", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	TH2D* h2d_true = new TH2D("h2d_true","h2d_true", p::true_nbins, p::true_bins, p::true_nbins, p::true_bins);
	TH2D* h2d_reco2 = new TH2D("h2d_reco2","h2d_reco2", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	TH2D* h2d_true2 = new TH2D("h2d_true2","h2d_true2", p::true_nbins, p::true_bins, p::true_nbins, p::true_bins);	

	//TH2D* h2d_reco = new TH2D("h2d_reco","h2d_reco", nthinslices+2, -1, nthinslices+1, nthinslices+2, -1, nthinslices+1);
	//TH2D* h2d_true = new TH2D("h2d_true","h2d_true", nthinslices+2, -1, nthinslices+1, nthinslices+2, -1, nthinslices+1);
	//TH2D* h2d_reco2 = new TH2D("h2d_reco2","h2d_reco2", nthinslices+2, -1, nthinslices+1, nthinslices+2, -1, nthinslices+1);
	//TH2D* h2d_true2 = new TH2D("h2d_true2","h2d_true2", nthinslices+2, -1, nthinslices+1, nthinslices+2, -1, nthinslices+1);
	Unfold uf(h2d_reco, h2d_true, h2d_reco2, h2d_true2);
	//Unfold uf(nthinslices+2, -1, nthinslices+1, h2d_reco, h2d_true, h2d_reco2, h2d_true2);
	//Unfold uf(nthinslices+2, -1, nthinslices+1, h2d_reco, h2d_true);
	//---------------------------------------------------------------------------------------------------------------------//

	//ThinSlice config. --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
	//SetOutputFileName(Form("prod4areco2_mc_ESliceE_dE%dMeV_%dslcs_beamxy_bmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold.root", name_thinslicewidth, nthinslices)); //output file name
	//SetOutputFileName(Form("prod4areco2_mc_ESliceE_dE%dMeV_%dslcs_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold.root", name_thinslicewidth, nthinslices)); //output file name
	//SetOutputFileName(Form("prod4areco2_mc_ESliceE_dE%dMeV_%dslcs_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2x2dunfold.root", name_thinslicewidth, nthinslices)); //output file name
	//SetOutputFileName(Form("prod4areco2_mc_ESliceE_dE%dMeV_%dslcs_beamxy_bmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold_constEloss.root", name_thinslicewidth, nthinslices)); //output file name
	//SetOutputFileName(Form("prod4areco2_mc_ESliceE_dE%dMeV_%dslcs_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold_constEloss.root", name_thinslicewidth, nthinslices)); //output file name
	//SetOutputFileName(Form("prod4areco2_mc_ESliceE_dE%dMeV_%dslcs_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold_constEloss_keff1sthitstudy.root", name_thinslicewidth, nthinslices)); //output file name
	//SetOutputFileName(Form("prod4areco2_mc_ESliceE_dE%dMeV_%dslcs_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold_constEloss_keff1sthitstudy_rmshorttrk.root", name_thinslicewidth, nthinslices)); //output file name
	//SetOutputFileName(Form("prod4areco2_mc_ESliceE_dE%dMeV_%dslcs_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold_constEloss_keff1sthitstudy.root", name_thinslicewidth, nthinslices)); //output file name
	//SetOutputFileName(Form("prod4areco2_mc_ESliceE_dE%dMeV_%dslcs_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold_constEloss_keff2ndhitstudy.root", name_thinslicewidth, nthinslices)); //output file name
	//SetOutputFileName(Form("prod4areco2_mc_ESliceE_dE%dMeV_%dslcs_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold_constEloss_modify_unphysicalbinassign.root", name_thinslicewidth, nthinslices)); //output file name
	//SetOutputFileName(Form("prod4areco2_mc_ESliceE_dE%dMeV_%dslcs_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold_constEloss_modify_upstreaminttreatment.root", name_thinslicewidth, nthinslices)); //output file name
	//SetOutputFileName(Form("prod4areco2_mc_ESliceE_dE%dMeV_%dslcs_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold_constEloss_modify_stendatthesamebin.root", name_thinslicewidth, nthinslices)); //output file name
	//SetOutputFileName(Form("prod4areco2_mc_ESliceE_dE%dMeV_%dslcs_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold_constEloss.root", name_thinslicewidth, nthinslices)); //output file name
	//SetOutputFileName(Form("prod4areco2_mc_ESliceE_dE%dMeV_%dslcs_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold_constEloss_exp_shorttruetrk.root", name_thinslicewidth, nthinslices)); //output file name
	//SetOutputFileName(Form("prod4areco2_mc_ESliceE_dE%dMeV_%dslcs_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold_constEloss_skip_shorttruetrk.root", name_thinslicewidth, nthinslices)); //output file name
	//SetOutputFileName(Form("prod4areco2_mc_ESliceE_dE%dMeV_%dslcs_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold_constEloss_0.1MC.root", name_thinslicewidth, nthinslices)); //output file name
	SetOutputFileName(Form("prod4areco2_mc_ESliceE_dynamicbin_beamxy_nobmrw_kebeamff_v09_39_01_ceil_ignoreincompleteSlice_2dunfold_constEloss.root")); //output file name
	//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

	//Basic configure ------//
	BetheBloch BB;
	BB.SetPdgCode(pdg);

	//book histograms --//
	BookHistograms();

	//dump txt files
	//ofstream myfile;
	//myfile.open ("true_ke_gt_800MeV.txt");

	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	Long64_t nbytes = 0, nb = 0;
	//bool isTestSample=true;
	int true_sliceID = -1, reco_sliceID = -1;
	int true_st_sliceID = -1, reco_st_sliceID = -1;
	int true_int_sliceID = -1, reco_int_sliceID = -1; //for 3D unfolding
	for (Long64_t jentry=0; jentry<nentries;jentry++) { //main entry loop
	//for (Long64_t jentry=0; jentry<0.1*nentries;jentry++) { //main entry loop
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		//isTestSample = true;
		//isTestSample = false; //validation only
		//isTestSample = true; //test only
		//if (ientry%2 == 0) isTestSample = false; //Divide MC sample by 2 parts: test+ufold
		//if (isTestSample) continue; //only validate sample
		//if (!isTestSample) continue; //only test sample

		true_sliceID = -1;
		reco_sliceID = -1;

		true_st_sliceID = -1;
		reco_st_sliceID = -1;

		true_int_sliceID = -1; //3d unfolding
		reco_int_sliceID = -1; //3d unfolding

		//only select protons	
		//if (primary_truth_Pdg!=pdg) continue; //only interested in protons
		if (beamtrackPdg!=pdg) continue; //only interested in protons
		//std::cout<<"beamtrackPdg:"<<beamtrackPdg<<std::endl;
		//n_tot++;

		//Event Selection Cut -- Part 1 ----------------------------------//
		bool IsBeamMatch=false; //if recostructed the right track (recoID=truthID)
		bool IsPandoraSlice=false; //pandora slice cut (can pandora reconstruct this track)
		bool IsCaloSize=false; //if calo size not empty
		bool IsIntersection=false; //if any track intersect with our reco track		
		if (primary_truth_Isbeammatched==1) IsBeamMatch=true;
		if (isprimarytrack==1&&isprimaryshower==0) IsPandoraSlice=true; 
		if (!primtrk_hitz->empty()) IsCaloSize=true;
		if (timeintersection->size()) IsIntersection=true;
		//----------------------------------------------------------------//

		//cout<<"\n"<<endl;
		//cout<<"run/subrun:event:"<<run<<" "<<subrun<<" "<<event<<endl;
		//cout<<"primaryID:"<<primaryID<<endl;
		//cout<<"IsPandoraSlice:"<<IsPandoraSlice<<" | isprimarytrack:"<<isprimarytrack<<" isprimaryshower:"<<isprimaryshower<<endl;
		//cout<<"IsCaloSize:"<<IsCaloSize<<endl;
		//cout<<"primary_truth_EndProcess:"<<primary_truth_EndProcess->c_str()<<endl;
		//cout<<"Isendpoint_outsidetpc:"<<Isendpoint_outsidetpc<<endl;
		//cout<<"IsBeamMatch:"<<IsBeamMatch<<endl;
		//cout<<"primary_truth_byE_origin="<<primary_truth_byE_origin<<""<<endl;
		//cout<<"primary_truth_byE_PDG="<<primary_truth_byE_PDG<<""<<endl;
		//cout<<"primary_truth_Pdg:"<<primary_truth_Pdg<<endl;
		//cout<<"beamtrackPdg:"<<beamtrackPdg<<endl;

		//Truth label of Primarytrack_End ------------------------------------------------------------------------------------------------//
		bool IsPureInEL=false; //inel
		bool IsPureEL=false; //el
		//bool IsPureMCS=false; //no hadron scattering

		//if (primary_truth_EndProcess->c_str()!=NULL) n_true_end++;
		if (strcmp(primary_truth_EndProcess->c_str(),"protonInelastic")==0) {
			IsPureInEL=true;
		}
		else { //hIoni
			IsPureEL=true;
			/*
			   if (interactionProcesslist->size()) { //size of interactionProcesslist >=0
			   cout<<"interactionProcesslist->size():"<<interactionProcesslist->size()<<endl;	
			   for(size_t iiii=0; iiii<interactionProcesslist->size(); iiii++) { //loop over all true interaction hits in this track
			   try {
			   double intx=interactionX->at(iiii);
			   double inty=interactionY->at(iiii);
			   double intz=interactionZ->at(iiii); 
			   cout<<"["<<iiii<<"] process:"<<interactionProcesslist->at(iiii)<<" z:"<<intz<<endl;

			   if(strcmp(interactionProcesslist->at(iiii).c_str(),"hadElastic")==0) {
			   IsPureEL=1;
			   }
			   }
			   catch (const std::out_of_range & ex) {
			   std::cout << "out_of_range Exception Caught :: interactionProcesslist" << ex.what() << std::endl;
			   n_processmap_error++;
			   }
			//if (intz<0) { //if interaction outside tpc
			//if(strcmp(interactionProcesslist->at(iiii).c_str(),"Transportation")!=0) {
			} //loop over all true interaction hits in this track 
			} //size of interactionProcesslist >=0
			*/
		} //hIoni

		//if (IsPureInEL==0&&IsPureEL==0) {
		//	//IsPureMCS=1;
		//}

		////if (strcmp(primary_truth_EndProcess->c_str(),"hIoni")==0) {
		////IsPureEL=true;
		////}
		////if (strcmp(primary_truth_EndProcess->c_str(),"CoulombScat")==0) {
		////IsPureMCS=true;
		////}
		//--------------------------------------------------------------------------------------------------------------------------------//

		//for (size_t j=0; j<beamtrk_z->size(); ++j) { //MCParticle loop
		//cout<<"beamtrk_z["<<j<<"]"<<beamtrk_z->at(j)<<" beamtrk_Eng["<<"]"<<beamtrk_Eng->at(j)<<endl;
		//} //MCParticle loop

		//Get true start/end point -----------------------------------------------------------------------//
		//double true_endz=-99; if (beamtrk_z->size()>1) true_endz=beamtrk_z->at(-1+beamtrk_z->size()); 
		//double true_endy=-99; if (beamtrk_y->size()>1) true_endy=beamtrk_y->at(-1+beamtrk_y->size()); 
		//double true_endx=-99; if (beamtrk_x->size()>1) true_endx=beamtrk_x->at(-1+beamtrk_x->size()); 
		double true_endz=primary_truth_EndPosition_MC[2]; 
		double true_endy=primary_truth_EndPosition_MC[1]; 
		double true_endx=primary_truth_EndPosition_MC[0];

		double true_stz=primary_truth_StartPosition_MC[2];
		double true_sty=primary_truth_StartPosition_MC[1];
		double true_stx=primary_truth_StartPosition_MC[0];

		bool IsTrueEndOutside=false;
		if (true_endz<0.) {
			IsTrueEndOutside=true;
		}
		//cout<<"trueEnd z/y/x:"<<true_endz<<"/"<<true_endy<<"/"<<true_endx<<endl;
		//cout<<"trueSt z/y/x:"<<true_stz<<"/"<<true_sty<<"/"<<true_stx<<endl;
		//cout<<"InEL EL MCS:"<<IsPureInEL<<" "<<IsPureEL<<" "<<IsPureMCS<<endl;
		//cout<<"InEL EL:"<<IsPureInEL<<" "<<IsPureEL<<" "<<endl;
		//cout<<"IsTrueEndOutside:"<<IsTrueEndOutside<<endl;
		//if (IsPureInEL==1) cout<<"Summary(TrueEnd, Endoutside, Bm, Orig, EPDG):("<<1<<", "<<IsTrueEndOutside<<", "<<IsBeamMatch<<", "<<primary_truth_byE_origin<<", "<<primary_truth_byE_PDG<<")"<<endl;	
		//if (IsPureEL==1) cout<<"Summary(TrueEnd, Endoutside, Bm, Orig, EPDG):("<<2<<", "<<IsTrueEndOutside<<", "<<IsBeamMatch<<", "<<primary_truth_byE_origin<<", "<<primary_truth_byE_PDG<<")"<<endl;	
		//if (IsPureMCS==1) cout<<"Summary(TrueEnd, Endoutside, Bm, Orig, EPDG):("<<3<<", "<<IsTrueEndOutside<<", "<<IsBeamMatch<<", "<<primary_truth_byE_origin<<", "<<primary_truth_byE_PDG<<")"<<endl;	

		//First point of MCParticle entering TPC ------------------------------------------------------------------------//
		bool is_beam_at_ff=false; //if the beam reach tpc
		int key_reach_tpc=-99;
		if (beamtrk_z->size()){
			for (size_t kk=0; kk<beamtrk_z->size(); ++kk) {  //loop over all beam hits
				double zpos_beam=beamtrk_z->at(kk);
				//if (1000.*beamtrk_Eng->at(kk)>800.) myfile<<"beamtrk_Eng["<<kk<<"]="<<1000.*beamtrk_Eng->at(kk)<<"\n";	
				if (zpos_beam>=0) {
					key_reach_tpc=(int)kk;
					break;
				}
			} //loop over all beam hits
			//Fill1DHist(KEtrue_Beam, (1000.*beamtrk_Eng->at(0)));
		} 
		if (key_reach_tpc!=-99) { is_beam_at_ff=true; }
		//cout<<"key_reach_tpc:"<<key_reach_tpc<<endl;	
		//cout<<"is_beam_at_ff:"<<is_beam_at_ff<<endl;

		//Get true trklen ---------------------------------------------------------------------------------------//
		int key_st = 0;
		double tmp_z = 9999;
		vector<double> true_trklen_accum;
		true_trklen_accum.reserve(beamtrk_z->size()); // initialize true_trklen_accum
		for (int iz=0; iz<(int)beamtrk_z->size(); iz++) {
			if (abs(beamtrk_z->at(iz)) < tmp_z){
				tmp_z = abs(beamtrk_z->at(iz));
				key_st = iz; // find the point where the beam enters the TPC (find the smallest abs(Z))
			}
			if (is_beam_at_ff) true_trklen_accum[iz] = 0.; // initialize true_trklen_accum [beam at ff]
			if (!is_beam_at_ff) true_trklen_accum[iz] = -1; // initialize true_trklen_accum [beam not at ff]
		}

		//fix on the truth length by adding distance between 1st tpc hit to front face ------------------------------------------------------//
		//[1] 3D projection on TPC front face
		double zproj_beam=0; //set beam z at ff
		double yproj_beam=0; //ini. value
		double xproj_beam=0; //ini. value
		int n_fit=3; //num of points used for fitting
		if (beamtrk_z->size()) {

			int key_fit_st=0;
			int key_fit_ed=-1+(int)beamtrk_z->size();
			if (key_reach_tpc!=-99) { //true track reach TPC FF
				key_fit_st=key_reach_tpc-1;
				key_fit_ed=key_reach_tpc+1;
			} //true track reach TPC FF

			//B.C.
			if (key_fit_st<0) key_fit_st=0;
			if (key_fit_ed>(-1+(int)beamtrk_z->size())) key_fit_ed=-1+(int)beamtrk_z->size();	

			//cout<<"beamtrk_z->size():"<<beamtrk_z->size()<<endl;
			//cout<<"key_reach_tpc:"<<key_reach_tpc<<endl;
			//std::cout<<"key_fit_st-ed:"<<key_fit_st<<"-"<<key_fit_ed<<std::endl;

			//start 3D line fit
			TGraph2D *gr=new TGraph2D();
			//cout<<"ck0"<<endl;
			//for (int N=key_fit_st; N<key_fit_ed; N++) {
			int nsize_fit=n_fit;
			if ((1+(key_fit_ed-key_fit_st))<n_fit) nsize_fit=1+(key_fit_ed-key_fit_st);
			//if (key_fit_ed>-1+(int)beamtrk_z->size()) n_fit=((int)beamtrk_z->size())-key_fit_st; //in case really short track
			for (int N=0; N<nsize_fit; N++) {
				gr->SetPoint(N, beamtrk_x->at(N+key_fit_st), beamtrk_y->at(N+key_fit_st), beamtrk_z->at(N+key_fit_st));
			}
			//cout<<"ck1"<<endl;
			//Initialization of parameters
			//int N=(int)Z_RECO.size();
			double ini_p1=(beamtrk_x->at(key_fit_ed)-beamtrk_x->at(key_fit_st))/(beamtrk_z->at(key_fit_ed)-beamtrk_z->at(key_fit_st));
			double ini_p0=beamtrk_x->at(key_fit_st)-ini_p1*beamtrk_z->at(key_fit_st);
			double ini_p3=beamtrk_y->at(key_fit_ed)-beamtrk_y->at(key_fit_st);
			double ini_p2=beamtrk_y->at(key_fit_st)-ini_p3*beamtrk_z->at(key_fit_st);
			//cout<<"ck2"<<endl;

			ROOT::Fit::Fitter  fitter;
			// make the functor objet
			SumDistance2 sdist(gr);
			ROOT::Math::Functor fcn(sdist,4);

			// set the function and the initial parameter values
			double pStart[4]={ini_p0, ini_p1, ini_p2, ini_p3};   
			fitter.SetFCN(fcn,pStart);
			//cout<<"ck3"<<endl;

			// set step sizes different than default ones (0.3 times parameter values)
			for (int ik = 0; ik < 4; ++ik) fitter.Config().ParSettings(ik).SetStepSize(0.01);
			//cout<<"ck4"<<endl;

			bool ok = fitter.FitFCN();
			//if (!ok) {
			//Error("line3Dfit","Line3D Fit failed");
			//return 1;
			//}
			//cout<<"ck5"<<endl;

			const ROOT::Fit::FitResult & result = fitter.Result();
			//std::cout << "Total final distance square " << result.MinFcnValue() << std::endl;
			//result.Print(std::cout);
			//cout<<"ck6"<<endl;

			// get fit parameters
			const double * parFit = result.GetParams();
			yproj_beam=result.Parameter(2)+result.Parameter(3)*zproj_beam;
			xproj_beam=result.Parameter(0)+result.Parameter(1)*zproj_beam;
			//cout<<"ck7"<<endl;

			delete gr;
		}

		//[2] Range compensation ----------------------------------------------------------//
		double range_true_patch=0;
		if (is_beam_at_ff) { //is beam at ff
			//calculate distance 1st hit and pojected point at TPC front face
			range_true_patch = sqrt( pow(beamtrk_x->at(key_reach_tpc)-xproj_beam, 2)+
					pow(beamtrk_y->at(key_reach_tpc)-yproj_beam, 2)+	
					pow(beamtrk_z->at(key_reach_tpc)-zproj_beam, 2) );
			//range_true_patch=0; //no fix on true len
		} //if entering tpc

		//true_trklen_accum
		double range_true=-9999;
		if (is_beam_at_ff) { //is beam at ff
			for (int iz=key_reach_tpc+1; iz<(int)beamtrk_z->size(); iz++) {
				if (iz == key_reach_tpc+1) range_true = range_true_patch;
				range_true += sqrt( pow(beamtrk_x->at(iz)-beamtrk_x->at(iz-1), 2)+
						pow(beamtrk_y->at(iz)-beamtrk_y->at(iz-1), 2)+	
						pow(beamtrk_z->at(iz)-beamtrk_z->at(iz-1), 2) );						    	
				true_trklen_accum[iz] = range_true;
			}
		} //is beam at ff						    	

		//fix on the truth length by adding distance between 1st tpc hit to front face ---------------------------------------------------------------------------//
		//cout<<"range_true:"<<range_true<<endl;
		//cout<<"key_st:"<<key_st<<endl;
		//for (size_t j=0; j<beamtrk_z->size(); ++j) { //MCParticle loop
		//cout<<"beamtrk_z["<<j<<"]:"<<beamtrk_z->at(j)<<" beamtrk_Eng["<<j<<"]:"<<beamtrk_Eng->at(j)<<" true_trklen_accum["<<j<<"]:"<<true_trklen_accum[j]<<endl;
		//} //MCParticle loop
		//Get reco info ----------------------------------------------------------------------------------//

		//Evt Classification =====================================================================//
		//signal -----------------------------//
		bool kinel=false;
		bool kel=false;
		//bool kmcs=false;
		if (IsBeamMatch) { //beam-match
			if (IsPureInEL) kinel=true;
			if (IsPureEL) kel=true;
			//if (IsPureMCS) kmcs=true;
		} //beam-match

		//background ------------------------------------------------------------------------//
		bool kMIDcosmic=false; //beam or cosmic
		bool kMIDpi=false; //+-pi
		bool kMIDp=false; //p
		bool kMIDmu=false; //mu
		bool kMIDeg=false; //e/gamma
		bool kMIDother=false; //other
		if (!IsBeamMatch) { //!beam-match
			if (primary_truth_byE_origin==2) { 
				kMIDcosmic=true;
			}
			else if (std::abs(primary_truth_byE_PDG)==211) {
				kMIDpi=true;
			}
			else if (primary_truth_byE_PDG==2212) {
				kMIDp=true;
			}
			else if (std::abs(primary_truth_byE_PDG)==13) {
				kMIDmu=true;
			}
			else if (std::abs(primary_truth_byE_PDG)==11 || primary_truth_byE_PDG==22) {
				kMIDeg=true;
			}
			else {
				kMIDother=true;
			}
		} //!beam-match	
		//cout<<"kMIDcosmic:"<<kMIDcosmic<<endl;
		//Evt Classification =====================================================================//

		//reco pos info & cut -----------------------------------------------------------------//
		double reco_stx=-99, reco_sty=-99, reco_stz=-99;
		double reco_endx=-99, reco_endy=-99, reco_endz=-99;
		bool IsPos=false;
		if (IsCaloSize) {
			reco_stx=primtrk_hitx->at(0); 
			reco_sty=primtrk_hity->at(0);
			reco_stz=primtrk_hitz->at(0);

			reco_endx=primtrk_hitx->at(primtrk_dedx->size()-1);	
			reco_endy=primtrk_hity->at(primtrk_dedx->size()-1);
			reco_endz=primtrk_hitz->at(primtrk_dedx->size()-1);

			//reco_startX_sce->Fill(reco_stx);
			//reco_startY_sce->Fill(reco_sty);
			//reco_startZ_sce->Fill(reco_stz);

			//Fill1DHist(reco_startX_sce, reco_stx);
			//Fill1DHist(reco_startY_sce, reco_sty);
			//Fill1DHist(reco_startZ_sce, reco_stz);

			double beam_dx=(reco_stx-mean_StartX)/sigma_StartX;
			double beam_dy=(reco_sty-mean_StartY)/sigma_StartY;
			double beam_dz=(reco_stz-mean_StartZ)/sigma_StartZ;
			double beam_dxy=sqrt(pow(beam_dx,2)+pow(beam_dy,2));	

			//hdeltaX->Fill(beam_dx);
			//hdeltaY->Fill(beam_dy);
			//hdeltaZ->Fill(beam_dz);
			//hdeltaXY->Fill(beam_dxy);

			if (beam_dx>=dx_min&&beam_dx<=dx_max) { //dx
				if (beam_dy>=dy_min&&beam_dy<=dy_max) { //dy
					if (beam_dz>=dz_min&&beam_dz<=dz_max) { //dz
						if (beam_dxy>=dxy_min&&beam_dxy<=dxy_max) { //dxy
							IsPos=true;
						} //dxy
					} //dz
				} //dy
			} //dx

		}

		//cosine_theta/cut ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
		bool IsCosine=false;
		double cosine_beam_spec_primtrk=-99; 
		//cosine_beam_spec_primtrk=beamDirx_spec->at(0)*primaryStartDirection[0]+beamDiry_spec->at(0)*primaryStartDirection[1]+beamDirz_spec->at(0)*primaryStartDirection[2]; //cosine between beam_spec and primary trk direction(trk before SCE corr)

		TVector3 dir;
		if (IsCaloSize) { //calosize	
			//trk direction after SCE corr.
			TVector3 pt0(primtrk_hitx->at(0), primtrk_hity->at(0), primtrk_hitz->at(0));
			TVector3 pt1(primtrk_hitx->at(-1+primtrk_hitx->size()), primtrk_hity->at(-1+primtrk_hity->size()), primtrk_hitz->at(-1+primtrk_hitz->size()));
			//TVector3 dir = pt1 - pt0;
			dir = pt1 - pt0;
			dir = dir.Unit();

			//beam direction
			//TVector3 beamdir(cos(beam_angleX_mc*TMath::Pi()/180), cos(beam_angleY_mc*TMath::Pi()/180), cos(beam_angleZ_mc*TMath::Pi()/180));
			TVector3 beamdir(beamDirx_spec->at(0),beamDiry_spec->at(0),beamDirz_spec->at(0));
			beamdir = beamdir.Unit();
			//beam_costh = dir.Dot(beamdir);
			cosine_beam_spec_primtrk=dir.Dot(beamdir);

			if (cosine_beam_spec_primtrk<0) { cosine_beam_spec_primtrk=-1.*cosine_beam_spec_primtrk; }
			//if (cosine_beam_spec_primtrk>cosine_beam_primtrk_min) { IsCosine=true; }
			if (cosine_beam_spec_primtrk>costh_min&&cosine_beam_spec_primtrk<costh_max) { IsCosine=true; }

		} //calosize

		//xy-cut (has been merged in the BQ cut in the new version)
		//bool IsXY=false;		
		//double x0_tmp=0, y0_tmp=0, z0_tmp=0; //start-pos, before sce
		//if (primaryEndPosition[2]>primaryStartPosition[2]) { //check if Pandora flip the sign
		//x0_tmp=primaryStartPosition[0];
		//y0_tmp=primaryStartPosition[1];
		//z0_tmp=primaryStartPosition[2];
		//} //check if Pandora flip the sign
		//else {
		//x0_tmp=primaryEndPosition[0];
		//y0_tmp=primaryEndPosition[1];
		//z0_tmp=primaryEndPosition[2];
		//}
		//if ((pow(((x0_tmp-mean_x)/dev_x),2)+pow(((y0_tmp-mean_y)/dev_y),2))<=1.) IsXY=true;

		bool IsMisidpRich=false;
		if (IsPos&&IsCaloSize&&IsPandoraSlice) {
			if (cosine_beam_spec_primtrk<=0.9) IsMisidpRich=true;
		}

		//reco trklen patch -------------------------------------------------------------------------------------------------------------------------------//
		//[0]key of 1st hit that z>=0, usually this key=0
		int key_reach_tpc_reco=-99;
		if (IsCaloSize) { //if calo size not empty
			for (size_t kk=0; kk<primtrk_hitz->size(); ++kk) {  //loop over all reco hits
				double zpos_reco=primtrk_hitz->at(kk);
				if (zpos_reco>=-1.) { //gt 0cm
					key_reach_tpc_reco=(int)kk;
					break;
				} //gt 0cm
			} //loop over all reco hits
		} //if calo size not empty

		//fix on the reco length by adding distance between 1st tpc hit to front face ------------------------------------------------------//
		/*
		//[1] 3D projection on TPC front face
		double zproj_reco=0; //set beam z at ff
		double yproj_reco=0; //ini. value
		double xproj_reco=0; //ini. value
		int n_fit_reco=3; //num of points used for fitting
		if (IsCaloSize&&key_reach_tpc_reco>-99) { //if calo size not empty
		if (primtrk_hitz->at(-1+primtrk_hitz->size())>=0) { //end_point inside tpc
		int key_fit_st=0;
		int key_fit_ed=2;
		if (key_reach_tpc_reco>=0) {
		key_fit_st=key_reach_tpc_reco;
		key_fit_ed=key_fit_st+2;
		}

		//B.C.
		if (key_fit_st<0) key_fit_st=0;
		if (key_fit_ed>(-1+(int)primtrk_hitz->size())) key_fit_ed=-1+(int)primtrk_hitz->size();	

		//if (primtrk_hitz->size()<=3) { 
		//key_fit_st=0;
		//key_fit_ed=-1+(int)primtrk_hitz->size(); //in case really short track
		//}
		//cout<<"\n(reco)key_fit_st-key_fit_ed:"<<key_fit_st<<"-"<<key_fit_ed<<endl;
		//cout<<"primtrk_hitz->size():"<<primtrk_hitz->size()<<endl;
		//cout<<"primtrk_hitz->at(max-1):"<<primtrk_hitz->at(-1+primtrk_hitz->size())<<endl;
		//cout<<"primtrk_hitz->at(FF):"<<primtrk_hitz->at(key_reach_tpc_reco)<<endl;

		//start 3D line fit
		TGraph2D *gr=new TGraph2D();
		int nsize_fit=n_fit_reco;
		if ((key_fit_ed-key_fit_st+1)<nsize_fit) nsize_fit=key_fit_ed-key_fit_st+1;
		//cout<<"check point0_0"<<endl;
		//cout<<"nsize_fit:"<<nsize_fit<<endl;
		for (int N=0; N<nsize_fit; N++) {
		gr->SetPoint(N, primtrk_hitx->at(N+key_fit_st), primtrk_hity->at(N+key_fit_st), primtrk_hitz->at(N+key_fit_st));
		//cout<<"(x,y,z):("<<primtrk_hitx->at(N+key_fit_st)<<","<<primtrk_hity->at(N+key_fit_st)<<","<<primtrk_hitz->at(N+key_fit_st)<<")"<<endl;
		//cout<<"check point0_x"<<endl;
		}
		//cout<<"check point0_1"<<endl;
		//cout<<"key_fit_st-key_fit_ed:"<<key_fit_st<<"-"<<key_fit_ed<<endl;
		double ini_p1=(primtrk_hitx->at(key_fit_ed)-primtrk_hitx->at(key_fit_st))/(primtrk_hitz->at(key_fit_ed)-primtrk_hitz->at(key_fit_st));
		double ini_p0=primtrk_hitx->at(key_fit_st)-ini_p1*primtrk_hitz->at(key_fit_st);
		double ini_p3=primtrk_hity->at(key_fit_ed)-primtrk_hity->at(key_fit_st);
		double ini_p2=primtrk_hity->at(key_fit_st)-ini_p3*primtrk_hitz->at(key_fit_st);
		//cout<<"check point0"<<endl;

		ROOT::Fit::Fitter  fitter;
		// make the functor objet
		SumDistance2 sdist(gr);
		ROOT::Math::Functor fcn(sdist,4);

		// set the function and the initial parameter values
		double pStart[4]={ini_p0, ini_p1, ini_p2, ini_p3};   
		fitter.SetFCN(fcn,pStart);

		// set step sizes different than default ones (0.3 times parameter values)
		for (int ik = 0; ik < 4; ++ik) fitter.Config().ParSettings(ik).SetStepSize(0.01);

		bool ok = fitter.FitFCN();
		//cout<<"ok:"<<ok<<endl;
		//if (!ok) {
		//Error("line3Dfit","Line3D Fit failed");
		//return 1;
		//}
		//cout<<"ck5"<<endl;

		const ROOT::Fit::FitResult & result = fitter.Result();
		//std::cout << "Total final distance square " << result.MinFcnValue() << std::endl;
		//result.Print(std::cout);

		//cout<<"check point1"<<endl;
		// get fit parameters
		const double * parFit = result.GetParams();
		yproj_reco=result.Parameter(2)+result.Parameter(3)*zproj_reco;
		xproj_reco=result.Parameter(0)+result.Parameter(1)*zproj_reco;
		//cout<<"(xproj_reco,yproj_reco,zproj_reco)=("<<xproj_reco<<","<<yproj_reco<<","<<zproj_reco<<")"<<endl;

		delete gr;
		} //end_point inside tpc
		} //if calo size not empty
		//cout<<"check point"<<endl;
		*/
			//-------------------------------------------------------------------------------------------------------------------------------------------------//


			//reco calorimetry ---------------------------------------------------------------------------//
			double range_reco_patch=0;
		//cout<<"\nkey_reach_tpc_reco:"<<key_reach_tpc_reco<<endl;
		//cout<<"primtrk_hitz->size()="<<primtrk_hitz->size()<<endl;
		//cout<<"(xproj_reco,yproj_reco,zproj_reco)=("<<xproj_reco<<","<<yproj_reco<<","<<zproj_reco<<")"<<endl;
		/*
		   if (IsCaloSize) { //if calo size not empty
		   if (key_reach_tpc_reco>-99) { //end_point inside tpc
		//cout<<"primtrk_hitz->at(TPC_FF)="<<primtrk_hitz->at(key_reach_tpc_reco)<<endl;
		//calculate distance 1st hit and pojected point at TPC front face
		range_reco_patch = sqrt( pow(primtrk_hitx->at(0)-xproj_reco, 2)+
		pow(primtrk_hity->at(0)-yproj_reco, 2)+	
		pow(primtrk_hitz->at(0)-zproj_reco, 2) );
		//range_true_reco=0; //no fix on reco len
		} //end_point inside tpc
		if (primtrk_hitz->at(0)<0) range_reco_patch=-range_reco_patch;
		} //if calo size not empty
		*/
		//cout<<"key_reach_tpc_reco="<<key_reach_tpc_reco<<endl;
		//cout<<"range_reco_patch="<<range_reco_patch<<endl;
		//cout<<"check pointYYYYYYYYYYYYYY"<<endl;

		int index_reco_endz=0;
		double wid_reco_max=-9999;
		double range_reco=0;
		vector<double> reco_trklen_accum;
		double reco_calo_MeV=0;
		double kereco_range=0;
		double kereco_range2=0;
		vector<double> EDept;
		vector<double> DEDX;
		vector<double> DX;
		double pid=-99; 
		vector<double> trkdedx;
		vector<double> trkres;
		if (IsCaloSize) { //if calo size not empty
			for (size_t h=0; h<primtrk_dedx->size(); ++h) { //loop over reco hits of a given track
				double hitx_reco=primtrk_hitx->at(h);
				double hity_reco=primtrk_hity->at(h);
				double hitz_reco=primtrk_hitz->at(h);
				double resrange_reco=primtrk_resrange->at(h);

				double dqdx=primtrk_dqdx->at(h);
				double pitch=primtrk_pitch->at(h);

				int wid_reco=primtrk_wid->at(-1+primtrk_wid->size()-h);
				double pt_reco=primtrk_pt->at(-1+primtrk_wid->size()-h);

				//if (wid_reco==-9999) continue; //outside TPC
				if (wid_reco>wid_reco_max) { 
					wid_reco_max=wid_reco;
					index_reco_endz=(int)-1+primtrk_wid->size()-h;
				}

				double cali_dedx=0.;
				cali_dedx=dedx_function_35ms(dqdx, hitx_reco, hity_reco, hitz_reco);

				EDept.push_back(cali_dedx*pitch);
				DEDX.push_back(cali_dedx);
				DX.push_back(pitch);

				//if (IsPureInEL) rangereco_dedxreco_TrueInEL->Fill(range_reco-resrange_reco, cali_dedx);
				//if (IsPureEL) rangereco_dedxreco_TrueEL->Fill(range_reco-resrange_reco, cali_dedx);
				//if (IsPureMCS) rangereco_dedxreco_TrueMCS->Fill(range_reco-resrange_reco, cali_dedx);

				//if (h==0) range_reco=range_reco_patch; //with range_patch correction: if z0>=0, add length; if z0<0, subtract length
				//if (h==0&&key_reach_tpc_reco<0) { //1st hit outside tpc
				//range_reco=-sqrt( pow(primtrk_hitx->at(h)-primtrk_hitx->at(key_reach_tpc_reco), 2)+
				//pow(primtrk_hity->at(h)-primtrk_hity->at(key_reach_tpc_reco), 2)+
				//pow(primtrk_hitz->at(h)-primtrk_hitz->at(key_reach_tpc_reco), 2) );
				//} //1st hit outside tpc
				if (h>=1) {
					range_reco += sqrt( pow(primtrk_hitx->at(h)-primtrk_hitx->at(h-1), 2)+
							pow(primtrk_hity->at(h)-primtrk_hity->at(h-1), 2)+
							pow(primtrk_hitz->at(h)-primtrk_hitz->at(h-1), 2) );
					//reco_trklen_accum[h] = range_reco;
					reco_trklen_accum.push_back(range_reco);
				}

				reco_calo_MeV+=cali_dedx*pitch;
				kereco_range+=pitch*dedx_predict(resrange_reco);
				kereco_range2+=pitch*(double)gr_predict_dedx_resrange->Eval(resrange_reco);

				/*
				   if (kinel) rangereco_dedxreco_TrueInEL->Fill(range_reco, cali_dedx);
				   if (kel) { 
				   rangereco_dedxreco_TrueEL->Fill(range_reco, cali_dedx);
				   rr_dedx_truestop->Fill(resrange_reco, cali_dedx);
				   }
				   */

				trkdedx.push_back(cali_dedx);
				trkres.push_back(resrange_reco);

			} //loop over reco hits of a given track

			pid=chi2pid(trkdedx,trkres); //pid using stopping proton hypothesis

		} //if calo size not empty
		//cout<<"check pointZZZZZZZZZZZZZZZZZZZ"<<endl;


		//if (reco_trklen_accum.size()) { //if calo size not empty
		//cout<<"primtrk_hitz->at(0):"<<primtrk_hitz->at(0)<<endl;
		//cout<<"primtrk_hitz->at(1):"<<primtrk_hitz->at(1)<<endl;
		//cout<<"reco_trklen_accum.at(0)="<<reco_trklen_accum.at(0)<<endl;
		//cout<<"key_reach_tpc_reco:"<<key_reach_tpc_reco<<endl;
		//for (int uuu=0; uuu<=key_reach_tpc_reco; ++uuu) {
		//cout<<"<primtrk_hitz->at("<<uuu<<")="<<primtrk_hitz->at(uuu)<<endl;
		//}
		//} //if calo size not empty




		//Reco stopping/Inel p cut ---------------------------------------------------------------------------------------------------------//
		bool IsRecoStop=false;
		bool IsRecoInEL=false;
		bool IsRecoEL=false;
		double mom_beam_spec=-99; mom_beam_spec=beamMomentum_spec->at(0);
		double bx_spec=beamPosx_spec->at(0);
		double by_spec=beamPosy_spec->at(0);
		//double range_reco=-99; if (!primtrk_range->empty()) range_reco=primtrk_range->at(0); //reco primary trklen

		//mc p0
		double btrk_px=-99; btrk_px=beamtrk_Px->at(0);
		double btrk_py=-99; btrk_py=beamtrk_Py->at(0);
		double btrk_pz=-99; btrk_pz=beamtrk_Pz->at(0);
		double btrk_p=sqrt(btrk_px*btrk_px+btrk_py*btrk_py+btrk_pz*btrk_pz);
		double mom_beam=btrk_p; //GeV/c
		double mom_beam_MeV=1000.*mom_beam;

		//beam XY cut to remove E-loss events upstream
		bool IsBeamXY=false;
		if ((pow(((bx_spec-meanX_mc)/(1.5*rmsX_mc)),2)+pow(((by_spec-meanY_mc)/(1.5*rmsY_mc)),2))<=1.) IsBeamXY=true;

		//beam-mom cut (within 3-sigma)
		bool IsBeamMom=false;
		if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) IsBeamMom=true;

		//beam quality cut --------------------------------------------------------------------------------------------------------------//
		bool IsBQ=false;
		//if (IsCosine&&IsPos) IsBQ=true;
		if (IsBeamXY&&IsBeamMom&&IsCosine&&IsPos) IsBQ=true;

		double csda_val_spec=csda_range_vs_mom_sm->Eval(mom_beam_spec);

		if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) IsRecoStop=true; //old cut
		//if ((range_reco/csda_val_spec)<min_norm_trklen_csda) IsRecoInEL=true; //old cut

		if ((range_reco/csda_val_spec)<min_norm_trklen_csda) { //inel region
			if (pid>pid_1) IsRecoInEL=true; 
			if (pid<=pid_1) IsRecoEL=true; 
		} //inel region
		if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) { //stopping p region
			if (pid>pid_2) IsRecoInEL=true; 
			if (pid<=pid_2) IsRecoEL=true;
		} //stopping p region

		//kinetic energies -------------------------------------------------------------------//
		//double ke_beam=1000.*p2ke(mom_beam); //ke_beam
		double ke_beam_spec=p2ke(mom_beam_spec); //ke_beam_spec [GeV]
		double ke_beam_spec_MeV=1000.*ke_beam_spec; //ke_beam_spec [MeV]
		double ke_trklen=1000.*ke_vs_csda_range_sm->Eval(range_reco); //[unit: MeV]
		double p_trklen=ke2p(ke_trklen);
		//double ke_simide=0;
		//for (int hk=0; hk<(int)primtrk_true_edept->size(); ++hk) { //loop over simIDE points
		//ke_simide+=primtrk_true_edept->at(hk);
		//} //loop over simIDE points

		double KE_ff=0;
		double KE_1st=0;
		//double KE_ff30=0;
		//if (is_beam_at_ff) KE_ff=1000.*beamtrk_Eng->at(key_reach_tpc); //unit:MeV
		if (is_beam_at_ff) { 		
			KE_ff=ke_ff; //use KE exactly at z=0
			KE_1st=1000.*beamtrk_Eng->at(key_reach_tpc);
			//double KE_1st_predict=BB.KEAtLength(KE_ff, range_true_patch);
			//KE_ff30=BB.KEAtLength(KE_ff, 30.);

			//h2d_R_kE1st->Fill(KE_1st, KE_1st_predict/KE_1st);
		}

		//KEff (reco) with const E-loss assumption ---------------------------------------------------------------------------------------------------------------------//
		//double mean_Elosscalo_stop=(4.95958e+01)/(1.00489e+00); //using fit [no bmrw]
		//double KE_ff_reco=ke_beam_spec_MeV-mean_Elosscalo_stop;
		//double KEend_reco=0;
		//KEend_reco=KE_ff_reco-reco_calo_MeV;

		double slope_kff_keffbeamR=0.859605;
		double ke_ffbeam_MeV=(ke_beam_spec_MeV-Eloss_mc_hy_stop)*R_fit_hy; //const E-loss with correction (R_fit_hy is essentially one)
		//double ke_ffbeam_MeV=ke_beam_spec_MeV-(p0_edept_stop+p1_edept_stop*ke_beam_spec_MeV+p2_edept_stop*pow(ke_beam_spec_MeV,2)); //E-dept E-loss (latest version!)
		//double ke_ffbeam_MeV=ke_ff; //KEff(truth)

		//Hypothetical length -------------------------------------------------------------------------------------//
		double fitted_length=-1; 
		double tmp_fitted_length=BB.Fit_dEdx_Residual_Length(trkdedx, trkres, pdg, false);
		if (tmp_fitted_length>0) fitted_length=tmp_fitted_length;
		double fitted_KE=-50; 
		if (fitted_length>0) fitted_KE=BB.KEFromRangeSpline(fitted_length);

		//Bethe-Bloch ---------------------------------------------------------//
		double kebb=-9999.; kebb=BB.KEAtLength(ke_ffbeam_MeV, range_reco);

		double KE_ff_reco=ke_ffbeam_MeV; //KE_ff_reco exactly at TPC FF [default]
		//cout<<"\ncheck0"<<endl;
		//if (reco_trklen_accum.size()>=2) { 
		//if (key_reach_tpc_reco!=-99) { 
		//KE_ff_reco=BB.KEAtLength(ke_ffbeam_MeV, reco_trklen_accum.at(1));
		//cout<<"reco_trklen_accum.at(2)="<<reco_trklen_accum.at(1)<<endl;
		//cout<<"KE_ff_reco="<<KE_ff_reco<<endl;
		//kebb=BB.KEAtLength(ke_ffbeam_MeV, range_reco);

		//KE_ff_reco=BB.KEAtLength(ke_ffbeam_MeV, range_reco_patch);

		//}
		//if (key_reach_tpc_reco!=-99&&primtrk_hitz->at(key_reach_tpc_reco)<0) KE_ff_reco=BB.KEAtLength(ke_ffbeam_MeV, range_reco_patch); //KE_ff_reco for the 1st hit [if size of hits !=0]
		//}
		//cout<<"check1\n"<<endl;
		double KEend_reco=kebb;

		//trklen cut ----------------------------------//
		bool IsRecoLONG_Trk=true;
		if (range_reco<=30.) IsRecoLONG_Trk=false;

		//bmrw ------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
		double mom_rw_minchi2=1; //weight for beam-momentum-reweight
		//if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) { //beam-mom (within 3-sigma)
		//mom_rw_minchi2=bmrw_func->Eval(mom_beam_spec*1000.); //bmrw, set weight if beam mom. within 3-sigma
		//} //beam-mom (within 3-sigma)
		//bmrw based on keff [const E-loss]
		//if ((ke_beam_spec_MeV-mean_Elosscalo_stop)>=mu_min&&(ke_beam_spec_MeV-mean_Elosscalo_stop)<=mu_max) mom_rw_minchi2=gng->Eval(ke_beam_spec_MeV-mean_Elosscalo_stop); //bmrw

		//if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) mom_rw_minchi2=kerw->Eval(ke_ffbeam_MeV); //new bmrw (const E-loss) [use me]
		//if (ke_ff>=mu_kemin&&ke_ff<=mu_kemax) mom_rw_minchi2=kerw->Eval(ke_ff); //new bmrw (using truth)
		//if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) mom_rw_minchi2=kerw_edept->Eval(ke_ffbeam_MeV); //new bmrw (edept E-loss) [use me::latest version]

		//misid:p rw ----------------------------------------//
		double misidp_rw=1;
		//if (kebb<=280.) misidp_rw=kerw_misidp->Eval(kebb);

		//KEend ---------------------------------------------------------------------------//
		//double KEend_true=0;
		double KEend_true=-1;
		if (beamtrk_Eng->size()) KEend_true=1000.*(beamtrk_Eng->at(-2+beamtrk_Eng->size()));

		//KEs ---------------------------------------------------------------------------------------//
		//double Eloss_upstream=0; 
		//if (KE_ff>0) Eloss_upstream=
		//double dEbb_true=0; if (range_true>=0&&KE_ff>0) dEbb_true=BB.KEAtLength(KE_ff, range_true);
		//double dEbb_reco=0; if (range_reco>=0&&KE_ff>0) dEbb_reco=BB.KEAtLength(KE_ff, range_reco);

		//double KEbb_true=-1; KEbb_true=KE_ff-BB.KEAtLength(KE_ff, range_true);
		//double KEbb_reco=-1; KEbb_reco=KE_ff-BB.KEAtLength(KE_ff, range_reco);
		//---------------------------------------------------------------------------------------------------------------//

		//countings -------------------------------------------------------//
		/*
		   if (IsPandoraSlice) n_pan_tot++;
		   if (IsPandoraSlice&&IsCaloSize) n_calsz_tot++;
		   if (IsPandoraSlice&&IsCaloSize&&IsBQ) n_bq_tot++;
		   if (IsPandoraSlice&&IsCaloSize&&IsBQ&&IsRecoInEL) n_recoinel_tot++;
		   if (IsPandoraSlice&&IsCaloSize&&IsBQ&&IsRecoEL) n_recoel_tot++;

		//[0]pure inel
		if (kinel) { //pure inel
		n_inel++;
		if (IsPandoraSlice) { //pandora
		n_inel_pan++;
		if (IsCaloSize) { //calosz
		n_inel_calsz++;
		if (IsBQ) { //bq
		n_inel_bq++;
		if (IsRecoInEL) { //reco inel
		n_inel_recoinel++;
		} //reco inel

		if (IsRecoEL) { //reco el
		n_inel_recoel++;
		} //reco el

		} //bq
		} //calosz
		} //pandora
		} //pure inel		

		//[1]pure el
		if (kel) { //pure el
		n_el++;
		if (IsPandoraSlice) { //pandora
		n_el_pan++;
		if (IsCaloSize) { //calosz
		n_el_calsz++;
		if (IsBQ) { //bq
		n_el_bq++;
		if (IsRecoInEL) { //reco inel
		n_el_recoinel++;
		} //reco inel

		if (IsRecoEL) { //reco el
		n_el_recoel++;
		} //reco el

		} //bq
		} //calosz
		} //pandora
		} //pure el		

		//[3]MID:Cosmic
		if (kMIDcosmic) { //mid:cosmic
		n_midcosmic++;
		if (IsPandoraSlice) { //pandora
		n_midcosmic_pan++;
		if (IsCaloSize) { //calosz
		n_midcosmic_calsz++;
		if (IsBQ) { //bq
		n_midcosmic_bq++;
		if (IsRecoInEL) { //reco inel
		n_midcosmic_recoinel++;
		} //reco inel
		if (IsRecoEL) { //reco el
		n_midcosmic_recoel++;
		} //reco el

		} //bq
		} //calosz
		} //pandora
		} //mid:cosmic

		//[4]MID:midpi
		if (kMIDpi) { //mid:pi
			n_midpi++;
			if (IsPandoraSlice) { //pandora
				n_midpi_pan++;
				if (IsCaloSize) { //calosz
					n_midpi_calsz++;
					if (IsBQ) { //bq
						n_midpi_bq++;
						if (IsRecoInEL) { //reco inel
							n_midpi_recoinel++;
						} //reco inel

						if (IsRecoEL) { //reco el
							n_midpi_recoel++;
						} //reco el

					} //bq
				} //calosz
			} //pandora
		} //mid:pi

		//[5]MID:midp
		if (kMIDp) { //mid:p
			n_midp++;
			if (IsPandoraSlice) { //pandora
				n_midp_pan++;
				if (IsCaloSize) { //calosz
					n_midp_calsz++;
					if (IsBQ) { //bq
						n_midp_bq++;
						if (IsRecoInEL) { //reco inel
							n_midp_recoinel++;
						} //reco inel

						if (IsRecoEL) { //reco el
							n_midp_recoel++;
						} //reco el

					} //bq
				} //calosz
			} //pandora
		} //mid:p

		//[6]MID:midmu
		if (kMIDmu) { //mid:mu
			n_midmu++;
			if (IsPandoraSlice) { //pandora
				n_midmu_pan++;
				if (IsCaloSize) { //calosz
					n_midmu_calsz++;
					if (IsBQ) { //bq
						n_midmu_bq++;
						if (IsRecoInEL) { //reco inel
							n_midmu_recoinel++;
						} //reco inel

						if (IsRecoEL) { //reco el
							n_midmu_recoel++;
						} //reco el

					} //bq
				} //calosz
			} //pandora
		} //mid:mu

		//[7]MID:mideg
		if (kMIDeg) { //mid:eg
			n_mideg++;
			if (IsPandoraSlice) { //pandora
				n_mideg_pan++;
				if (IsCaloSize) { //calosz
					n_mideg_calsz++;
					if (IsBQ) { //bq
						n_mideg_bq++;
						if (IsRecoInEL) { //reco inel
							n_mideg_recoinel++;
						} //reco inel

						if (IsRecoEL) { //reco el
							n_mideg_recoel++;
						} //reco el

					} //bq
				} //calosz
			} //pandora
		} //mid:eg

		//[8]MID:midother
		if (kMIDother) { //mid:other
			n_midother++;
			if (IsPandoraSlice) { //pandora
				n_midother_pan++;
				if (IsCaloSize) { //calosz
					n_midother_calsz++;
					if (IsBQ) { //bq
						n_midother_bq++;
						if (IsRecoInEL) { //reco inel
							n_midother_recoinel++;
						} //reco inel

						if (IsRecoEL) { //reco el
							n_midother_recoel++;
						} //reco el

					} //bq
				} //calosz
			} //pandora
		} //mid:other
		*/
			//countings -------------------------------------------------------//

			//ntrklen --------------------------------------------------------------------//
			//if (IsPandoraSlice&&IsCaloSize&&IsBQ) { } 
			//bkg-rich
			//if (IsMisidpRich) { }

			//true trklen vs ke ------------------------------------------------------------------//
			/*
			   if (!beamtrk_Eng->empty()) { //if true container not empty
			   for (int hk=0; hk<(int)beamtrk_Eng->size()-1; ++hk) { //loop over true hits
			   double thisLen=true_trklen_accum[hk];
			   double this_incE = 1000.*beamtrk_Eng->at(hk); //MeV
			   } //loop over true hits
			   } //if true container not empty
			   */

		//E-slice method ------------------------------------------------------------------------------------------------------------//
		//true start slice ID
		//true_st_sliceID=int((Emax-KE_ff)/thinslicewidth+0.5); //old sliceID assignment
		//true_st_sliceID=int(ceil((Emax-KE_ff)/thinslicewidth)); //ignore the incomplete slideID
		//true_st_sliceID=(Emax-KE_ff)/thinslicewidth; //thin-slice

      		for (true_st_sliceID=0; true_st_sliceID<p::true_nbins-2; ++true_st_sliceID) {
        		if (KE_ff > p::true_KE[true_st_sliceID]) break;
      		}

		//if (true_st_sliceID<0) true_st_sliceID=0; //KE higher than Emax, put it to over-flow bin
		//if (true_st_sliceID >= nthinslices) true_st_sliceID = nthinslices;

		//double KE_true=BB.KEAtLength(KE_ff, range_true);
		//true_sliceID = int((Emax-KE_true)/thinslicewidth);
		double KE_true=KEend_true;
		//true_sliceID = int(floor((Emax-KE_true)/thinslicewidth));
		//true_sliceID = (Emax-KE_true)/thinslicewidth;
		//if (true_sliceID >= nthinslices||KE_true<0) true_sliceID = nthinslices;
		//if (true_sliceID < 0) true_sliceID = 0; //KE higher than Emax, put it to over-flow bin
		//if (true_sliceID >= nthinslices) true_sliceID = nthinslices;
      		for (true_sliceID=0; true_sliceID<p::true_nbins-2; ++true_sliceID) {
        		if (KE_true > p::true_KE[true_sliceID+1]) break;
      		}

		//sansity check -----------------//
		//place where xs rise
		//if (KE_ff>480.&&KE_ff<520.) {
		//if (true_st_sliceID>true_sliceID) {
		//std::cout<<"\n[Me!] true_stz:"<<true_stz<<" true_endz:"<<true_endz<<" KE_ff:"<<KE_ff<<" KE_1st:"<<KE_1st<<" KEend_true:"<<KEend_true<<" | true_st_sliceID:"<<true_st_sliceID<<" true_sliceID:"<<true_sliceID<<", range_true:"<<range_true<<", is_beam_at_ff:"<<is_beam_at_ff<<std::endl;
		//}
		//else {
		//std::cout<<"\n[Check!] true_stz:"<<true_stz<<" true_endz:"<<true_endz<<" KE_ff:"<<KE_ff<<" KE_1st:"<<KE_1st<<" KEend_true:"<<KEend_true<<" | true_st_sliceID:"<<true_st_sliceID<<" true_sliceID:"<<true_sliceID<<", range_true:"<<range_true<<", is_beam_at_ff:"<<is_beam_at_ff<<std::endl;
		//}
		//}
		//-------------------------------//

		//if (true_endz < 0) { //Upstream-interaction
		//true_st_sliceID=-1;
		//true_sliceID = -1;
		//Note. XS measurement starts from ID=0, ID=-1 will not be considered in the XS measurement
		//p.s. ID=-1 serves as over-flow bin and un-physical bin
		//}


		//if (true_st_sliceID==true_sliceID) {
		//std::cout<<"\n[same ID!] true_stz:"<<true_stz<<" true_endz:"<<true_endz<<" KE_ff:"<<KE_ff<<" KE_1st:"<<KE_1st<<" KEend_true:"<<KEend_true<<" | true_st_sliceID:"<<true_st_sliceID<<" true_sliceID:"<<true_sliceID<<std::endl;
		//}
		if (true_st_sliceID>true_sliceID) { //get rid of super short track
			//The case happens when the track length is very short, if true_st_sliceID>true_sliceID, treat this case as unphysical, discard this type of measurement
			//if true_st_sliceID==true_sliceID for short tracks, physical and will be considered in the XS measurement
			//std::cout<<"\ntrue_stz:"<<true_stz<<" true_endz:"<<true_endz<<" KE_ff:"<<KE_ff<<" KE_1st:"<<KE_1st<<" KEend_true:"<<KEend_true<<" | true_st_sliceID:"<<true_st_sliceID<<" true_sliceID:"<<true_sliceID<<std::endl;
			true_st_sliceID=-1;
			true_sliceID=-1;

		} //get rid of super short track

		//for 3D unfolding ----------------//
		true_int_sliceID=true_sliceID;
		if (!IsPureInEL) {
			true_int_sliceID=-1;
		}
		//---------------------------------//

		//evt selection cuts
		bool PassCuts_INT=false; //all bq cut+reco inel cut
		bool PassCuts_INC=false; //all bq cut
		if (IsPandoraSlice&&IsCaloSize) {
			//reco_sliceID = int(primtrk_range->at(0)/thinslicewidth);
			//reco_sliceID = int(reco_endz/thinslicewidth);
			//reco_sliceID = int(range_reco/thinslicewidth);
			//cout<<"reco_endz:"<<reco_endz<<"  reco_sliceID:"<<reco_sliceID<<endl;

			//double keff_reco=ke_beam_spec_MeV-mean_Eloss_upstream;
			//double keff_reco=KE_ff_reco;
			//double KE_reco30=BB.KEAtLength(keff_reco, 30.);
			//reco_st_sliceID=int((Emax-KE_ff_reco)/thinslicewidth+0.5); //old slideID assignment
			//if (range_reco>=30.) reco_st_sliceID=int((Emax-KE_reco30)/thinslicewidth);
			//if (range_reco<30) reco_st_sliceID=-1;
			//reco_st_sliceID=int(ceil((Emax-KE_ff_reco)/thinslicewidth)); //ignore incomplete slideID assignment
			//reco_st_sliceID=(Emax-KE_ff_reco)/thinslicewidth;
			//if (reco_st_sliceID<0) reco_st_sliceID=0; //KE higher than Emax
			//if (reco_st_sliceID > nthinslices) reco_st_sliceID = nthinslices;
			//if (reco_endz < 0) reco_st_sliceID = -1; //un-physical

      			for (reco_st_sliceID=0; reco_st_sliceID<p::reco_nbins-2; ++reco_st_sliceID) {
        			if (KE_ff_reco > p::reco_KE[reco_st_sliceID]) break;
      			}

			//double KE_reco=BB.KEAtLength(keff_reco, range_reco);
			//reco_sliceID = int((Emax-KEend_reco)/thinslicewidth);
			//reco_sliceID = int(floor((Emax-KEend_reco)/thinslicewidth));
			//reco_sliceID = (Emax-KEend_reco)/thinslicewidth;
			//if (reco_sliceID >= nthinslices||KEend_reco<0) reco_sliceID = nthinslices;
			//if (reco_sliceID >= nthinslices) reco_sliceID = nthinslices;
			//if (reco_sliceID < 0) reco_sliceID = 0;
			//if (reco_endz<0) reco_sliceID = -1;

      			for (reco_sliceID=0; reco_sliceID<p::reco_nbins-2; ++reco_sliceID) {
        			if (KEend_reco > p::reco_KE[reco_sliceID+1]) break;
      			}


			if (reco_st_sliceID>reco_sliceID) { //incomplete charge deposition
				//std::cout<<"\nreco_stz:"<<reco_stz<<" reco_endz:"<<reco_endz<<" KE_ff_reco:"<<KE_ff_reco<<" KEend_reco:"<<KEend_reco<<" reco_st_sliceID:"<<reco_st_sliceID<<" reco_sliceID:"<<reco_sliceID<<std::endl;
				reco_st_sliceID=-1;
				reco_sliceID=-1;
			} //incomplete charge deposition

			//for 3d/2x2d unfolding ---------//
			reco_int_sliceID=reco_sliceID;
			//------------------------------//

			if (IsBQ&&IsRecoInEL) {
				//if (IsRecoLONG_Trk&&IsBQ&&IsRecoInEL) {
				PassCuts_INT=true; //for INT 
			}
			if (IsBQ) {
				//if (IsRecoLONG_Trk&&IsBQ) {
				PassCuts_INC=true; //for INC
				if (!IsRecoInEL) { //if NOT pass RecoInel selection
					reco_int_sliceID=-1; //for 3d/2x2d unfolding
				} //if NOT pass RecoInel selection
			}

			//double this_calo_MeV=0;
			//double this_KE=KE_ff_reco;
			//for (size_t ih=0; ih<primtrk_dedx->size(); ++ih) {
			//double thisLen=reco_trklen_accum[ih];
			//double thisKE=BB.KEAtLength(keff_reco, thisLen); //len2ke conversion
			//KEbb_recotrklen_all->Fill(thisLen, thisKE); 

			//this_KE-=EDept.at(ih);	
			//KEcalo_recotrklen_all->Fill(thisLen, this_KE); 
			//}

			}

			//if (isTestSample) { //if test sample
			cout<<"\ncheck_0"<<endl;
			h_truesliceid_all->Fill(true_sliceID, mom_rw_minchi2);
			h_true_st_sliceid_all->Fill(true_st_sliceID, mom_rw_minchi2);
			//} //if test sample
			//else { //if NOT test sample
			cout<<"\ncheck_0_1"<<endl;
			uf.eff_den_Inc->Fill(true_sliceID, mom_rw_minchi2);
			uf.eff_den_st_Inc->Fill(true_st_sliceID, mom_rw_minchi2);
			cout<<"\ncheck_0_2"<<endl;
			//for (int ij=true_st_sliceID; ij<=true_sliceID; ++ij){
			//if (true_sliceID < nthinslices && true_sliceID>=0){
			//for (int ij=0; ij<=true_sliceID+1; ++ij){
			//for (int ij=0; ij<=nthinslices+1; ++ij) {
			//if (ij<nthinslices) ++true_incidents[ij+1];
			//if (ij==(true_sliceID+1)) ++true_incidents[ij];
			//if (ij==(true_st_sliceID+1)) ++true_st_incidents[ij];
			//}
			//} //if NOT test sample
			if (PassCuts_INC&&IsBeamMatch) { //if passing all basic cuts
				//if (isTestSample) { //if test sample
				h_recosliceid_cuts->Fill(reco_sliceID, mom_rw_minchi2); 
				h_truesliceid_cuts->Fill(true_sliceID, mom_rw_minchi2);

				h_reco_st_sliceid_cuts->Fill(reco_st_sliceID, mom_rw_minchi2); 
				h_true_st_sliceid_cuts->Fill(true_st_sliceID, mom_rw_minchi2);
				//} //if test sample
				//else{ //if NOT test sample
				uf.eff_num_Inc->Fill(true_sliceID, mom_rw_minchi2);
				uf.eff_num_st_Inc->Fill(true_st_sliceID, mom_rw_minchi2);

				uf.pur_num_Inc->Fill(reco_sliceID, mom_rw_minchi2);
				uf.pur_num_st_Inc->Fill(reco_st_sliceID, mom_rw_minchi2);

				uf.response_SliceID_Inc.Fill(reco_sliceID, true_sliceID, mom_rw_minchi2);
				uf.response_st_SliceID_Inc.Fill(reco_st_sliceID, true_st_sliceID, mom_rw_minchi2);
				uf.response_SliceID_2D.Fill(reco_st_sliceID, reco_sliceID, true_st_sliceID, true_sliceID, mom_rw_minchi2);
				uf.response_SliceID_Int_2D.Fill(reco_sliceID, reco_int_sliceID, true_sliceID, true_int_sliceID, mom_rw_minchi2);

				uf.res_Inc_reco->Fill(reco_sliceID, mom_rw_minchi2);
				uf.res_st_Inc_reco->Fill(reco_st_sliceID, mom_rw_minchi2);

				uf.res_Inc_truth->Fill(true_sliceID, mom_rw_minchi2);
				uf.res_st_Inc_truth->Fill(true_st_sliceID, mom_rw_minchi2);
				//} //if NOT test sample
			} //if passing all basic cuts
			else { //if NOT passing all cuts
				//if (!isTestSample){
				uf.response_SliceID_Inc.Miss(true_sliceID, mom_rw_minchi2);
				uf.response_st_SliceID_Inc.Miss(true_st_sliceID, mom_rw_minchi2);
				uf.response_SliceID_2D.Miss(true_st_sliceID, true_sliceID, mom_rw_minchi2);
				uf.response_SliceID_Int_2D.Miss(true_sliceID, true_int_sliceID, mom_rw_minchi2);

				uf.res_Inc_truth->Fill(true_sliceID, mom_rw_minchi2);
				uf.res_st_Inc_truth->Fill(true_st_sliceID, mom_rw_minchi2);
				//std::cout<<true_sliceID<<std::endl;
				//}
			} //if NOT passing all cuts

			//INT histograms ------------------------------------------------------------//
			if (IsPureInEL) { //pure inel
				//n_test_recoinel++;
				//if (IsBeamMatch) n_kinel++;

				//if (isTestSample){
				h_truesliceid_inelastic_all->Fill(true_sliceID, mom_rw_minchi2);
				//if (true_sliceID<=0) {
				//std::cout<<"WRONG WRONG!!:: true_sliceID of IsPureInEL="<<true_sliceID<<" | KEend_true="<<KEend_true<<" | is_beam_at_ff="<<is_beam_at_ff<<" | true_endz="<<true_endz<<std::endl;
				//}
				//}
				//else{ //NOT test sample for unfolding
				uf.eff_den_Int->Fill(true_sliceID, mom_rw_minchi2);
				//n_test_recoinel_sample++;
				//if (IsBeamMatch) n_kinel2++;
				//} //NOT test sample for unfolding

				if (PassCuts_INT&&IsBeamMatch) { //if pass reco inel cuts
					//if (isTestSample){
					h_recosliceid_inelastic_cuts->Fill(reco_sliceID, mom_rw_minchi2);
					h_truesliceid_inelastic_cuts->Fill(true_sliceID, mom_rw_minchi2);
					//}
					//else{
					uf.eff_num_Int->Fill(true_sliceID, mom_rw_minchi2);
					uf.pur_num_Int->Fill(reco_sliceID, mom_rw_minchi2);
					uf.response_SliceID_Int.Fill(reco_sliceID, true_sliceID, mom_rw_minchi2);

					uf.res_Int_reco->Fill(reco_sliceID, mom_rw_minchi2);
					uf.res_Int_truth->Fill(true_sliceID, mom_rw_minchi2);
					//}
				} //if pass reco inel cuts
				else{ //if NOT pass all basic cuts
					//if (!isTestSample) { 
					uf.response_SliceID_Int.Miss(true_sliceID, mom_rw_minchi2);
					uf.res_Int_truth->Fill(true_sliceID, mom_rw_minchi2);
					//uf.response_SliceID_Int.Miss(true_sliceID, mom_rw_minchi2);
					//}
				} //if not pass all basic cuts
			} //pure inel

			if (PassCuts_INC) { //if pass all cuts
				//if (isTestSample){ //if test sample
				//h_recosliceid_allevts_cuts->Fill(reco_sliceID, mom_rw_minchi2);
				//h_reco_st_sliceid_allevts_cuts->Fill(reco_st_sliceID, mom_rw_minchi2);
				h2d_recosliceid_allevts_cuts->Fill(reco_st_sliceID, reco_sliceID, mom_rw_minchi2);
				h2d_recosliceid_recoinelastic_cuts->Fill(reco_sliceID, reco_int_sliceID, mom_rw_minchi2);

				//h_truesliceid_allevts_cuts->Fill(true_sliceID, mom_rw_minchi2);
				//h_true_st_sliceid_allevts_cuts->Fill(true_st_sliceID, mom_rw_minchi2);
				h2d_truesliceid_allevts_cuts->Fill(true_st_sliceID, true_sliceID, mom_rw_minchi2);

				if (kinel) { 
					//h_recosliceid_allevts_cuts_inel->Fill(reco_sliceID, mom_rw_minchi2);
					//h_reco_st_sliceid_allevts_cuts_inel->Fill(reco_st_sliceID, mom_rw_minchi2);
					h2d_recosliceid_allevts_cuts_inel->Fill(reco_st_sliceID, reco_sliceID, mom_rw_minchi2);
					h2d_recosliceid_recoinelastic_cuts_inel->Fill(reco_sliceID, reco_int_sliceID, mom_rw_minchi2);

					//h_truesliceid_allevts_cuts_inel->Fill(true_sliceID, mom_rw_minchi2);
					//h_true_st_sliceid_allevts_cuts_inel->Fill(true_st_sliceID, mom_rw_minchi2);
					h2d_truesliceid_allevts_cuts_inel->Fill(true_st_sliceID, true_sliceID, mom_rw_minchi2);
				}
				if (kel) { 
					//h_recosliceid_allevts_cuts_el->Fill(reco_sliceID, mom_rw_minchi2);
					//h_reco_st_sliceid_allevts_cuts_el->Fill(reco_st_sliceID, mom_rw_minchi2);
					h2d_recosliceid_allevts_cuts_el->Fill(reco_st_sliceID, reco_sliceID, mom_rw_minchi2);
					h2d_recosliceid_recoinelastic_cuts_el->Fill(reco_sliceID, reco_int_sliceID, mom_rw_minchi2);

					//h_truesliceid_allevts_cuts_el->Fill(true_sliceID, mom_rw_minchi2);
					//h_true_st_sliceid_allevts_cuts_el->Fill(true_st_sliceID, mom_rw_minchi2);
					h2d_truesliceid_allevts_cuts_el->Fill(true_st_sliceID, true_sliceID, mom_rw_minchi2);
				}
				if (kMIDcosmic) { 
					//h_recosliceid_allevts_cuts_midcosmic->Fill(reco_sliceID, mom_rw_minchi2);
					//h_reco_st_sliceid_allevts_cuts_midcosmic->Fill(reco_st_sliceID, mom_rw_minchi2);
					h2d_recosliceid_allevts_cuts_midcosmic->Fill(reco_st_sliceID, reco_sliceID, mom_rw_minchi2);
					h2d_recosliceid_recoinelastic_cuts_midcosmic->Fill(reco_sliceID, reco_int_sliceID, mom_rw_minchi2);

					//h_truesliceid_allevts_cuts_midcosmic->Fill(true_sliceID, mom_rw_minchi2);
					//h_true_st_sliceid_allevts_cuts_midcosmic->Fill(true_st_sliceID, mom_rw_minchi2);
					h2d_truesliceid_allevts_cuts_midcosmic->Fill(true_st_sliceID, true_sliceID, mom_rw_minchi2);
				}
				if (kMIDpi) { 
					//h_recosliceid_allevts_cuts_midpi->Fill(reco_sliceID, mom_rw_minchi2);
					//h_reco_st_sliceid_allevts_cuts_midpi->Fill(reco_st_sliceID, mom_rw_minchi2);
					h2d_recosliceid_allevts_cuts_midpi->Fill(reco_st_sliceID, reco_sliceID, mom_rw_minchi2);
					h2d_recosliceid_recoinelastic_cuts_midpi->Fill(reco_sliceID, reco_int_sliceID, mom_rw_minchi2);

					//h_truesliceid_allevts_cuts_midpi->Fill(true_sliceID, mom_rw_minchi2);
					//h_true_st_sliceid_allevts_cuts_midpi->Fill(true_st_sliceID, mom_rw_minchi2);
					h2d_truesliceid_allevts_cuts_midpi->Fill(true_st_sliceID, true_sliceID, mom_rw_minchi2);
				}
				if (kMIDp) { 
					//h_recosliceid_allevts_cuts_midp->Fill(reco_sliceID, mom_rw_minchi2*misidp_rw);
					//h_reco_st_sliceid_allevts_cuts_midp->Fill(reco_st_sliceID, mom_rw_minchi2);
					h2d_recosliceid_allevts_cuts_midp->Fill(reco_st_sliceID, reco_sliceID, mom_rw_minchi2);
					h2d_recosliceid_recoinelastic_cuts_midp->Fill(reco_sliceID, reco_int_sliceID, mom_rw_minchi2);

					//h_truesliceid_allevts_cuts_midp->Fill(true_sliceID, mom_rw_minchi2*misidp_rw);
					//h_true_st_sliceid_allevts_cuts_midp->Fill(true_st_sliceID, mom_rw_minchi2);
					h2d_truesliceid_allevts_cuts_midp->Fill(true_st_sliceID, true_sliceID, mom_rw_minchi2);
				}
				if (kMIDmu) { 
					//h_recosliceid_allevts_cuts_midmu->Fill(reco_sliceID, mom_rw_minchi2);
					//h_reco_st_sliceid_allevts_cuts_midmu->Fill(reco_st_sliceID, mom_rw_minchi2);
					h2d_recosliceid_allevts_cuts_midmu->Fill(reco_st_sliceID, reco_sliceID, mom_rw_minchi2);
					h2d_recosliceid_recoinelastic_cuts_midmu->Fill(reco_sliceID, reco_int_sliceID, mom_rw_minchi2);

					//h_truesliceid_allevts_cuts_midmu->Fill(true_sliceID, mom_rw_minchi2);
					//h_true_st_sliceid_allevts_cuts_midmu->Fill(true_st_sliceID, mom_rw_minchi2);
					h2d_truesliceid_allevts_cuts_midmu->Fill(true_st_sliceID, true_sliceID, mom_rw_minchi2);
				}
				if (kMIDeg) { 
					//h_recosliceid_allevts_cuts_mideg->Fill(reco_sliceID, mom_rw_minchi2);
					//h_reco_st_sliceid_allevts_cuts_mideg->Fill(reco_st_sliceID, mom_rw_minchi2);
					h2d_recosliceid_allevts_cuts_mideg->Fill(reco_st_sliceID, reco_sliceID, mom_rw_minchi2);
					h2d_recosliceid_recoinelastic_cuts_mideg->Fill(reco_sliceID, reco_int_sliceID, mom_rw_minchi2);

					//h_truesliceid_allevts_cuts_mideg->Fill(true_sliceID, mom_rw_minchi2);
					//h_true_st_sliceid_allevts_cuts_mideg->Fill(true_st_sliceID, mom_rw_minchi2);
					h2d_truesliceid_allevts_cuts_mideg->Fill(true_st_sliceID, true_sliceID, mom_rw_minchi2);
				}
				if (kMIDother) { 
					//h_recosliceid_allevts_cuts_midother->Fill(reco_sliceID, mom_rw_minchi2);
					//h_reco_st_sliceid_allevts_cuts_midother->Fill(reco_st_sliceID, mom_rw_minchi2);
					h2d_recosliceid_allevts_cuts_midother->Fill(reco_st_sliceID, reco_sliceID, mom_rw_minchi2);
					h2d_recosliceid_recoinelastic_cuts_midother->Fill(reco_sliceID, reco_int_sliceID, mom_rw_minchi2);

					//h_truesliceid_allevts_cuts_midother->Fill(reco_sliceID, mom_rw_minchi2);
					//h_true_st_sliceid_allevts_cuts_midother->Fill(reco_st_sliceID, mom_rw_minchi2);
					h2d_truesliceid_allevts_cuts_midother->Fill(true_st_sliceID, true_sliceID, mom_rw_minchi2);
				}
				//} //if test sample
				//else { //if NOT test sample
				uf.pur_den->Fill(reco_sliceID, mom_rw_minchi2);
				uf.pur_den_Inc->Fill(reco_sliceID, mom_rw_minchi2);
				uf.pur_den_st_Inc->Fill(reco_st_sliceID, mom_rw_minchi2);
				//} //if NOT test sample
			} //if pass all cuts
			cout<<"check_1"<<endl;

			if (PassCuts_INT) { //if pass reco inel cut
				//if (isTestSample){ //if test sample
				h_recosliceid_recoinelastic_cuts->Fill(reco_sliceID, mom_rw_minchi2);
				h_truesliceid_recoinelastic_cuts->Fill(true_sliceID, mom_rw_minchi2);
				if (kinel) { 
					h_recosliceid_recoinelastic_cuts_inel->Fill(reco_sliceID, mom_rw_minchi2);
					h_truesliceid_recoinelastic_cuts_inel->Fill(true_sliceID, mom_rw_minchi2);
				}
				if (kel) { 
					h_recosliceid_recoinelastic_cuts_el->Fill(reco_sliceID, mom_rw_minchi2);
					h_truesliceid_recoinelastic_cuts_el->Fill(true_sliceID, mom_rw_minchi2);
				}
				if (kMIDcosmic) { 
					h_recosliceid_recoinelastic_cuts_midcosmic->Fill(reco_sliceID, mom_rw_minchi2);
					h_truesliceid_recoinelastic_cuts_midcosmic->Fill(true_sliceID, mom_rw_minchi2);
				}
				if (kMIDpi) { 
					h_recosliceid_recoinelastic_cuts_midpi->Fill(reco_sliceID, mom_rw_minchi2);
					h_truesliceid_recoinelastic_cuts_midpi->Fill(true_sliceID, mom_rw_minchi2);
				}
				if (kMIDp) { 
					h_recosliceid_recoinelastic_cuts_midp->Fill(reco_sliceID, mom_rw_minchi2*misidp_rw);
					h_truesliceid_recoinelastic_cuts_midp->Fill(true_sliceID, mom_rw_minchi2*misidp_rw);
				}
				if (kMIDmu) { 
					h_recosliceid_recoinelastic_cuts_midmu->Fill(reco_sliceID, mom_rw_minchi2);
					h_truesliceid_recoinelastic_cuts_midmu->Fill(true_sliceID, mom_rw_minchi2);
				}
				if (kMIDeg) { 
					h_recosliceid_recoinelastic_cuts_mideg->Fill(reco_sliceID, mom_rw_minchi2);
					h_truesliceid_recoinelastic_cuts_mideg->Fill(true_sliceID, mom_rw_minchi2);
				}
				if (kMIDother) { 
					h_recosliceid_recoinelastic_cuts_midother->Fill(reco_sliceID, mom_rw_minchi2);
					h_truesliceid_recoinelastic_cuts_midother->Fill(true_sliceID, mom_rw_minchi2);
				}
				//} //if test sample
				//else { //if NOT test sample
				uf.pur_den_Int->Fill(reco_sliceID, mom_rw_minchi2);
				//} //if NOT test sample
			} //if pass reco inel cut

			//reco/truth KEs
			if (IsPureInEL) { //is pure inelastic
				//if (!isTestSample){ //if validation sample
				//if (true_sliceID < nthinslices && true_sliceID>=0){
				//++true_interactions[true_sliceID+1];
				//}
				//} //if validation sample

				//reco ke
				/*
				   if (IsCaloSize==true&&IsBeamMatch==true) { //calo & beam_match
				   std::vector<std::vector<double>> vincE(nthinslices);
				   double thisKE=KE_ff_reco;
				//for (size_t ih=0; ih<zreco_rawindex.size(); ++ih) {
				for (size_t ih=0; ih<primtrk_dedx->size(); ++ih) {
				//double this_calo_z=zreco_widreco[ih].second;
				//int this_sliceID = int(this_calo_z/thinslicewidth);
				//double this_calo_len=zreco_lenreco[ih].second;
				//double this_dE=zreco_de[ih].second;
				double thisZ=primtrk_hitz->at(ih);
				double thisLen=reco_trklen_accum[ih];

				//int this_sliceID = int(thisLen/thinslicewidth);
				//ke_reco-=this_dE;

				int this_sliceID=-1;
				thisKE-=EDept.at(ih);
				//double thiskeff_reco=ke_beam_spec_MeV-mean_Eloss_upstream;
				//double thisKE=BB.KEAtLength(thiskeff_reco, thisLen); //len2ke conversion
				this_sliceID=int((Emax-thisKE)/thinslicewidth);
				if (this_sliceID < 0) this_sliceID = -1;
				if (thisZ < 0) this_sliceID = -1; //up-stram INT, set trklen_accum to -1 by default
				if (this_sliceID >= nthinslices) this_sliceID = nthinslices;
				//KEbb_recotrklen_inel->Fill(thisLen,thisKE);
				//KEcalo_recotrklen_inel->Fill(thisLen, thisKE);

				if (this_sliceID>=nthinslices) continue;
				if (this_sliceID<0) continue;

				double this_incE = thisKE;
				vincE[this_sliceID].push_back(this_incE);


				KEreco_z_inel->Fill(thisZ, this_incE);
				KEreco_range_inel->Fill(thisLen, this_incE);
				reco_z_range_inel->Fill(thisZ, thisLen);

				dEdx_range_inel->Fill(thisLen, DEDX.at(ih));
				dx_range_inel->Fill(thisLen, DX.at(ih));
				dE_range_inel->Fill(thisLen, EDept.at(ih));

				}

				//reco_AngCorr->Fill(dir.Z());
				} //calo & beam_match
				*/

				//truth KEs
				/*
				   if (!beamtrk_Eng->empty()) { //if true container not empty
				//if (!beamtrk_Eng->empty()&&(1000.*beamtrk_Eng->at(0)<600.)) { //if true container not empty
				std::vector<std::vector<double>> vincE_true(nthinslices);
				for (int hk=0; hk<(int)beamtrk_Eng->size()-1; ++hk) { //loop over true hits
				double thisZ=beamtrk_z->at(hk);
				//int this_sliceID = int(thisZ/thinslicewidth);
				double thisLen=true_trklen_accum[hk];
				//int this_sliceID = int(thisLen/thinslicewidth);
				double this_incE = 1000.*beamtrk_Eng->at(hk); //MeV
				//double thisKE=BB.KEAtLength(KE_ff,thisLen); //len2ke conversion
				int this_sliceID=-1;
				this_sliceID=int((Emax-this_incE)/thinslicewidth);

				if (this_sliceID < 0) this_sliceID = -1;
				if (thisZ < 0) this_sliceID = -1; //up-stram INT, set trklen_accum to -1 by default
				if (this_sliceID >= nthinslices) this_sliceID = nthinslices;
				//if (thisLen < 0) true_sliceID = -1; //up-stram INT, set trklen_accum to -1 by default

				//KEbb_truetrklen_inel->Fill(thisLen, this_incE);
				//KEcalo_truetrklen_inel->Fill(thisLen, this_incE);

				if (this_sliceID>=nthinslices) continue;
				if (this_sliceID<0) continue;
				vincE_true[this_sliceID].push_back(this_incE);

				//test
				//double fX1=20.4585; //X of 1st point
				//double fY1=591.912; //Y of 1st point
				//double fX2=106.418; //X of 2nd point
				//double fY2=250.525; //Y of 2nd point
				//double m=(fY2-fY1)/(fX2-fX1);
				//double b=fY1-m*fX1;

				//if (this_incE>(b+m*thisLen)) {
				//KEtrue_z_inel->Fill(thisZ, this_incE);
				//KEtrue_range_inel->Fill(thisLen, this_incE);
				//true_z_range_inel->Fill(thisZ, thisLen);
				//true_z_range_KEtrue_inel->Fill(thisZ, thisLen, this_incE);
				//}	
				}//loop over true hits (last point always has KE = 0)

				//KEtrue_Beam_inel->Fill(1000.*beamtrk_Eng->at(0));
				//KEtrue_ff_inel->Fill(ke_ff);

				} //if true container not empty
				*/
			} //if pure inelastic


			} //main entry loop

			//save true inc & int arrays to histograms -------------------------------------//
			//for (int iii = 0; iii<nthinslices+2; ++iii){
			//h_true_incidents->SetBinContent(iii+1, true_incidents[iii]);
			//h_true_st_incidents->SetBinContent(iii+1, true_st_incidents[iii]);
			//h_true_interactions->SetBinContent(iii+1, true_interactions[iii]);
			//}

			//save results -------//
			uf.SaveHistograms();
			SaveHistograms();


			//myfile.close();

			//end_time ---------------------------------------------------------------------------//
			auto t_end = std::chrono::high_resolution_clock::now();

			// floating-point duration: no duration_cast needed
			std::chrono::duration<double, std::milli> fp_ms = t_st - t_end;

			// integral duration: requires duration_cast
			auto int_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_st);

			// converting integral duration to integral duration of shorter divisible time unit:
			// no duration_cast needed
			std::chrono::duration<long, std::micro> int_usec = int_ms;

			std::cout << "\n\nf() took " << fp_ms.count() << " ms, "
				<< "or " << int_ms.count() << " whole milliseconds "
				<< "(which is " << int_usec.count() << " whole microseconds)" << std::endl;
			std::cout<<"\n Execution time:"<<(int_ms.count()/1000.)/60.<<" (min.)"<<std::endl;
			//------------------------------------------------------------------------------------//



			}
