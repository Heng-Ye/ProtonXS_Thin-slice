#define ProtonBetheBlochKE_cxx
#include "ProtonBetheBlochKE.h"

#include <TH2.h>
#include <TH1.h>
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
#include "TVector3.h"

#include <stdio.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "./cali/dedx_function_35ms.h"
#include "./headers/BasicParameters.h"
#include "./headers/BasicFunctions.h"
#include "./headers/util.h"
//#include "./headers/SliceParams.h"
//#include "./headers/ThinSlice.h"
#include "./headers/BetheBloch.h"
#include "./headers/sg_smmoth.h"

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
			std::cout << "Total Initial distance square = " << sum << std::endl;
		}
		first = false;
		return sum;
	}

};
/////////////////////////////////


void ProtonBetheBlochKE::Loop() {
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	Long64_t nbytes = 0, nb = 0;
	bool isTestSample=true;
	//int true_sliceID = -1, reco_sliceID = -1;

	//mean and sigma of KEbeam_spec
	double mean_kebeam=4.38618e+02;
	double err_mean_kebeam=1.73088e-01;
	double sigma_kebeam=4.49648e+01;
	double err_sigma_kebeam=1.34969e-01;

	//energy loss
	double mean_mc_kebeam_stop=438.783;
	double err_mean_mc_kebeam_stop=0.167262;
	double sigma_mc_kebeam_stop=44.763;
	double err_sigma_mc_kebeam_stop=0.120428;

	double mean_mc_kerange_stop=405.356;
	double err_mean_mc_kerange_stop=0.225953;
	double sigma_mc_kerange_stop=41.8589;
	double err_sigma_mc_kerange_stop=0.153558;

	double Eloss_mc=33.4266;
	double err_Eloss_mc=0.281126;

	//Basic configure ------//
	//BetheBloch BB;
	//BB.SetPdgCode(pdg);

	//Kinetic energies -------------------------------------------------------------------------------------------------------//
	//int nke=320;
	int nke=160;
	double kemin=-800;
	double kemax=800;

	//ff
	TH1D *h1d_kebeam=new TH1D("h1d_kebeam","", nke, kemin, kemax); //ke from beamline inst.
	TH1D *h1d_ke0=new TH1D("h1d_ke0","", nke, kemin, kemax); //ke from ke_truth0
	
	TH1D *h1d_keff_truth=new TH1D("h1d_keff_truth","", nke, kemin, kemax); //keff_truth
	TH1D *h1d_keff_truth_el=new TH1D("h1d_keff_truth_el","", nke, kemin, kemax); //keff_truth
	TH1D *h1d_keff_truth_inel=new TH1D("h1d_keff_truth_inel","", nke, kemin, kemax); //keff_truth

	TH1D *h1d_keff_reco=new TH1D("h1d_keff_reco","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_keff_reco_el=new TH1D("h1d_keff_reco_el","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_keff_reco_inel=new TH1D("h1d_keff_reco_inel","", nke, kemin, kemax); //keff_reco

	TH1D *h1d_keff_reco_constErange=new TH1D("h1d_keff_reco_constErange","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_keff_reco_constErange_el=new TH1D("h1d_keff_reco_constErange_el","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_keff_reco_constErange_inel=new TH1D("h1d_keff_reco_constErange_inel","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_keff_reco_constErange_midcosmic=new TH1D("h1d_keff_reco_constErange_midcosmic","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_keff_reco_constErange_midpi=new TH1D("h1d_keff_reco_constErange_midpi","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_keff_reco_constErange_midp=new TH1D("h1d_keff_reco_constErange_midp","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_keff_reco_constErange_midmu=new TH1D("h1d_keff_reco_constErange_midmu","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_keff_reco_constErange_mideg=new TH1D("h1d_keff_reco_constErange_mideg","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_keff_reco_constErange_midother=new TH1D("h1d_keff_reco_constErange_midother","", nke, kemin, kemax); //keff_reco


	TH1D *h1d_Edept_reco=new TH1D("h1d_Edept_reco","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_Edept_reco_el=new TH1D("h1d_Edept_reco_el","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_Edept_reco_inel=new TH1D("h1d_Edept_reco_inel","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_Edept_reco_midcosmic=new TH1D("h1d_Edept_reco_midcosmic","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_Edept_reco_midpi=new TH1D("h1d_Edept_reco_midpi","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_Edept_reco_midp=new TH1D("h1d_Edept_reco_midp","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_Edept_reco_midmu=new TH1D("h1d_Edept_reco_midmu","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_Edept_reco_mideg=new TH1D("h1d_Edept_reco_mideg","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_Edept_reco_midother=new TH1D("h1d_Edept_reco_midother","", nke, kemin, kemax); //keff_reco


	TH1D *h1d_keff_reco_Edept_range_el=new TH1D("h1d_keff_reco_Edept_range_el","", nke, kemin, kemax); //
	TH1D *h1d_keff_reco_Edept_calo_el=new TH1D("h1d_keff_reco_Edept_calo_el","", nke, kemin, kemax); //

	TH1D *h1d_keff_reco_Edept_range_inel=new TH1D("h1d_keff_reco_Edept_range_inel","", nke, kemin, kemax); //
	TH1D *h1d_keff_reco_Edept_calo_inel=new TH1D("h1d_keff_reco_Edept_calo_inel","", nke, kemin, kemax); //

	//normalized trklen
	TH1D *h1d_ntrklen=new TH1D("h1d_ntrklen","", 140, 0, 1.4); //
	TH1D *h1d_ntrklen_recostop=new TH1D("h1d_ntrklen_recostop","", 140, 0, 1.4); //
	TH1D *h1d_ntrklen_el=new TH1D("h1d_ntrklen_el","", 140, 0, 1.4); //
	TH1D *h1d_ntrklen_inel=new TH1D("h1d_ntrklen_inel","", 140, 0, 1.4); //

	//some surveys
	//TH2D *h2d_keff_r0ratio_inel=new TH2D("h2d_keff_r0ratio_inel","", nke, kemin, kemax,1000,-50,50); //
	//TH2D *h2d_ntrklen_r0ratio_inel=new TH2D("h2d_ntrklen_r0ratio_inel","", 140, 0, 1.4, 100 ,-50,50); //
	

	//end-truth
	TH1D *h1d_keend_truth_el=new TH1D("h1d_keend_truth_el","", nke, kemin, kemax);
	TH1D *h1d_keend_truth_inel=new TH1D("h1d_keend_truth_inel","", nke, kemin, kemax);

	//TH1D *h1d_kebb_truth_el=new TH1D("h1d_kebb_truth_el","", nke, kemin, kemax); //keff_truth+range_truth
	//TH1D *h1d_kebb_truth_inel=new TH1D("h1d_kebb_truth_inel","", nke, kemin, kemax); 
	
	//end-reco
	//kebb
	//TH1D *h1d_kebb_reco_constErange_el=new TH1D("h1d_kebb_reco_constErange_el","", nke, kemin, kemax); //keff_reco
	//TH1D *h1d_kebb_reco_constErange_inel=new TH1D("h1d_kebb_reco_constErange_inel","", nke, kemin, kemax); //keff_reco

	//TH1D *h1d_kebb_reco_Edept_range_el=new TH1D("h1d_kebb_reco_Edept_range_el","", nke, kemin, kemax); //
	//TH1D *h1d_kebb_reco_Edept_range_inel=new TH1D("h1d_kebb_reco_Edept_range_inel","", nke, kemin, kemax); //

	//kebb
	TH1D *h1d_kebb_reco_constErange=new TH1D("h1d_kebb_reco_constErange","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_reco_constErange_el=new TH1D("h1d_kebb_reco_constErange_el","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_reco_constErange_inel=new TH1D("h1d_kebb_reco_constErange_inel","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_reco_constErange_midcosmic=new TH1D("h1d_kebb_reco_constErange_midcosmic","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_reco_constErange_midpi=new TH1D("h1d_kebb_reco_constErange_midpi","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_reco_constErange_midp=new TH1D("h1d_kebb_reco_constErange_midp","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_reco_constErange_midmu=new TH1D("h1d_kebb_reco_constErange_midmu","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_reco_constErange_mideg=new TH1D("h1d_kebb_reco_constErange_mideg","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_reco_constErange_midother=new TH1D("h1d_kebb_reco_constErange_midother","", nke, kemin, kemax); //keff_reco


	TH1D *h1d_kebb_recoInel_constEcalo=new TH1D("h1d_kebb_recoInel_constEcalo","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_recoInel_constEcalo_el=new TH1D("h1d_kebb_recoInel_constEcalo_el","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_recoInel_constEcalo_inel=new TH1D("h1d_kebb_recoInel_constEcalo_inel","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_recoInel_constEcalo_midcosmic=new TH1D("h1d_kebb_recoInel_constEcalo_midcosmic","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_recoInel_constEcalo_midpi=new TH1D("h1d_kebb_recoInel_constEcalo_midpi","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_recoInel_constEcalo_midp=new TH1D("h1d_kebb_recoInel_constEcalo_midp","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_recoInel_constEcalo_midmu=new TH1D("h1d_kebb_recoInel_constEcalo_midmu","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_recoInel_constEcalo_mideg=new TH1D("h1d_kebb_recoInel_constEcalo_mideg","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_recoInel_constEcalo_midother=new TH1D("h1d_kebb_recoInel_constEcalo_midother","", nke, kemin, kemax); //keff_reco


	TH1D *h1d_kebb_recoEl_constEcalo=new TH1D("h1d_kebb_recoEl_constEcalo","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_recoEl_constEcalo_el=new TH1D("h1d_kebb_recoEl_constEcalo_el","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_recoEl_constEcalo_inel=new TH1D("h1d_kebb_recoEl_constEcalo_inel","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_recoEl_constEcalo_midcosmic=new TH1D("h1d_kebb_recoEl_constEcalo_midcosmic","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_recoEl_constEcalo_midpi=new TH1D("h1d_kebb_recoEl_constEcalo_midpi","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_recoEl_constEcalo_midp=new TH1D("h1d_kebb_recoEl_constEcalo_midp","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_recoEl_constEcalo_midmu=new TH1D("h1d_kebb_recoEl_constEcalo_midmu","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_recoEl_constEcalo_mideg=new TH1D("h1d_kebb_recoEl_constEcalo_mideg","", nke, kemin, kemax); //keff_reco
	TH1D *h1d_kebb_recoEl_constEcalo_midother=new TH1D("h1d_kebb_recoEl_constEcalo_midother","", nke, kemin, kemax); //keff_reco

	
	TH1D *h1d_kebb_reco_Edept_calo_el=new TH1D("h1d_kebb_reco_Edept_calo_el","", nke, kemin, kemax); //keff_reco
	//TH1D *h1d_kebb_reco_Edept_calo_inel=new TH1D("h1d_kebb_reco_Edept_calo_inel","", nke, kemin, kemax); //keff_reco

	//KEend(reco)
	//all
/*
	TH1D *h1d_kebb_reco=new TH1D("h1d_kebb_reco","", nke, kemin, kemax);
	TH1D *h1d_kebb_reco_inel=new TH1D("h1d_kebb_reco_inel","", nke, kemin, kemax);
	TH1D *h1d_kebb_reco_el=new TH1D("h1d_kebb_reco_el","", nke, kemin, kemax);
	TH1D *h1d_kebb_reco_midcosmic=new TH1D("h1d_kebb_reco_midcosmic","", nke, kemin, kemax);
	TH1D *h1d_kebb_reco_midpi=new TH1D("h1d_kebb_reco_midpi","", nke, kemin, kemax);
	TH1D *h1d_kebb_reco_midp=new TH1D("h1d_kebb_reco_midp","", nke, kemin, kemax);
	TH1D *h1d_kebb_reco_midmu=new TH1D("h1d_kebb_reco_midmu","", nke, kemin, kemax);
	TH1D *h1d_kebb_reco_mideg=new TH1D("h1d_kebb_reco_mideg","", nke, kemin, kemax);
	TH1D *h1d_kebb_reco_midother=new TH1D("h1d_kebb_reco_midother","", nke, kemin, kemax);
*/
	TH2D *h2d_trklen_KEbb=new TH2D("h2d_trklen_KEbb","",140, 0, 140, nke, kemin, kemax);

	//RecoInel
/*
	TH1D *h1d_kebb_recoInel=new TH1D("h1d_kebb_recoInel","", nke, kemin, kemax);
	TH1D *h1d_kebb_recoInel_inel=new TH1D("h1d_kebb_recoInel_inel","", nke, kemin, kemax);
	TH1D *h1d_kebb_recoInel_el=new TH1D("h1d_kebb_recoInel_el","", nke, kemin, kemax);
	TH1D *h1d_kebb_recoInel_midcosmic=new TH1D("h1d_kebb_recoInel_midcosmic","", nke, kemin, kemax);
	TH1D *h1d_kebb_recoInel_midpi=new TH1D("h1d_kebb_recoInel_midpi","", nke, kemin, kemax);
	TH1D *h1d_kebb_recoInel_midp=new TH1D("h1d_kebb_recoInel_midp","", nke, kemin, kemax);
	TH1D *h1d_kebb_recoInel_midmu=new TH1D("h1d_kebb_recoInel_midmu","", nke, kemin, kemax);
	TH1D *h1d_kebb_recoInel_mideg=new TH1D("h1d_kebb_recoInel_mideg","", nke, kemin, kemax);
	TH1D *h1d_kebb_recoInel_midother=new TH1D("h1d_kebb_recoInel_midother","", nke, kemin, kemax);
*/
	TH2D *h2d_trklen_KEbb_RecoInel=new TH2D("h2d_trklen_KEbb_RecoInel","",140, 0, 140, nke, kemin, kemax);

	//RecoEl
/*
	TH1D *h1d_kebb_recoEl=new TH1D("h1d_kebb_recoEl","", nke, kemin, kemax);
	TH1D *h1d_kebb_recoEl_inel=new TH1D("h1d_kebb_recoEl_inel","", nke, kemin, kemax);
	TH1D *h1d_kebb_recoEl_el=new TH1D("h1d_kebb_recoEl_el","", nke, kemin, kemax);
	TH1D *h1d_kebb_recoEl_midcosmic=new TH1D("h1d_kebb_recoEl_midcosmic","", nke, kemin, kemax);
	TH1D *h1d_kebb_recoEl_midpi=new TH1D("h1d_kebb_recoEl_midpi","", nke, kemin, kemax);
	TH1D *h1d_kebb_recoEl_midp=new TH1D("h1d_kebb_recoEl_midp","", nke, kemin, kemax);
	TH1D *h1d_kebb_recoEl_midmu=new TH1D("h1d_kebb_recoEl_midmu","", nke, kemin, kemax);
	TH1D *h1d_kebb_recoEl_mideg=new TH1D("h1d_kebb_recoEl_mideg","", nke, kemin, kemax);
	TH1D *h1d_kebb_recoEl_midother=new TH1D("h1d_kebb_recoEl_midother","", nke, kemin, kemax);
*/
	TH2D *h2d_trklen_KEbb_RecoEl=new TH2D("h2d_trklen_KEbb_RecoEl","",140, 0, 140, nke, kemin, kemax);

	//E-loss scan
	const int n_eloss=40;
	vector<double> Eloss;
	TH1D *h1d_KEend_tune_bmrw[n_eloss];
	double Eloss_start=37.;
	double dEloss=0.5;
	for (int k=0; k<n_eloss; ++k) {
		double eloss_step=Eloss_start+(double)k*dEloss;
		Eloss.push_back(eloss_step);

		h1d_KEend_tune_bmrw[k]=new TH1D(Form("h1d_KEend_tune_bmrw_%d",k),"", nke, kemin, kemax);
	}

	//------------------------------------------------------------------------------------------------------------------------//

	//Beam momentum reweighting ----------------------------------------------------------------------------------------------//
	//MC KE beam Gaussian 
	double m1=3.89270e+02; //mc, keff with const E-loss
	double s1=4.49638e+01; //mc, keff with const E-loss
	double a1=7.06341e+02; //mc, keff with const E-loss

	double m2=3.93027e+02; //data, keff with const E-loss
	double s2=5.18623e+01; //data, keff with const E-loss
	double a2=6.09665e+02; //data, keff with const E-loss
		
	double xmin=0.; //pmin [MeV]
	double xmax=1000.; //pmax [MeV]

	double mu_min=m1-3.*s1;
	double mu_max=m1+3.*s1;

	TF1 *agng=new TF1(Form("agng"),agovg,xmin,xmax,6);
	agng->SetParameter(0,m1);
	agng->SetParameter(1,s1);
	agng->SetParameter(2,a1);

	agng->SetParameter(3,m2);
	agng->SetParameter(4,s2);
	agng->SetParameter(5,a2);


	//MC beam momentum -----------------------//
	double mm1=1007.1482; //MC prod4a [spec]
	double ss1=60.703307; //MC prod4a [spec]
	double mmu_min=mm1-3.*ss1;
	double mmu_max=mm1+3.*ss1;

	double xxmin=0.; //pmin [MeV/c]
	double xxmax=2000.; //pmax [MeV/c]
	TF1 *g1=new TF1("g1",fitg,xxmin,xxmax,2);
	g1->SetName("g1");
	g1->SetParameter(0,mm1);
	g1->SetParameter(1,ss1);

	//mu range
	double dmu=0.0005;
	double mu_st=1.01;
	int nmu=71;

	double dsigma=0.002;
	//double sigma_st=1.5;
	//int nsigma=250;
	double sigma_st=1.6;
	int nsigma=350;

	//mu x sigma
	const int n_mu_sigma=(const int)nmu*nsigma;
	int n_1d=nmu*nsigma; 
	TF1 **gn=new TF1*[n_mu_sigma];
	TF1 **gng=new TF1*[n_mu_sigma];

	//use trklen as an observable for reweighting
	//TH1D *h1d_trklen_rw[n_mu_sigma];

	int cnt_array=0;
	int index_original=0;
	//int index_minchi2=13331; //index of minchi2(this index is wrong)
	//int index_minchi2=17537; //index of minchi2
	int index_minchi2=11623; //index of minchi2(with beamXY cut)
	for (int imu=0; imu<nmu; ++imu){ //mu loop
		double frac_mu=mu_st-(double)imu*dmu;
		double mu=mm1*frac_mu;
		for (int isigma=0; isigma<nsigma; ++isigma){ //sigma loop
			double frac_sigma=sigma_st-(double)isigma*dsigma;
			double sigma=ss1*frac_sigma;

			//if (mu==m1&&sigma==s1) { //no rw
			if (std::abs(mu-mm1)<0.0001&&std::abs(sigma-ss1)<0.0001) { //no rw
				index_original=cnt_array;
				mu=mm1;
				sigma=ss1;
			} //no rw

			//Gaussian with changed mean and sigma
			gn[cnt_array]=new TF1(Form("gn_%d",cnt_array),fitg,xxmin,xxmax,2);
			gn[cnt_array]->SetParameter(0,mu);
			gn[cnt_array]->SetParameter(1,sigma);

			//weighting func. (beam mom)
			gng[cnt_array]=new TF1(Form("gng_%d",cnt_array),govg,xxmin,xxmax,4);
			gng[cnt_array]->SetParameter(0,mm1);
			gng[cnt_array]->SetParameter(1,ss1);
			gng[cnt_array]->SetParameter(2,mu);
			gng[cnt_array]->SetParameter(3,sigma);

			//prepare rw histograms
			//h1d_trklen_rw[cnt_array]=new TH1D(Form("h1d_trklen_rw_%d",cnt_array),Form("f_{#mu}:%.2f f_{#sigma}:%.2f #oplus RecoStop Cut",frac_mu,frac_sigma),n_b,b_min,b_max);
			//h1d_trklen_rw[cnt_array]->GetXaxis()->SetTitle("Track Length [cm]");

			cnt_array++;
			} //sigma loop
	} //mu loop

	//------------------------------------------------------------------------------------------------------------------------//

	for (Long64_t jentry=0; jentry<nentries;jentry++) { //main entry loop
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
                if (jentry%10000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;

		isTestSample = true;
		if (ientry%2 == 0) isTestSample = false; //Divide MC sample by 2 parts: test+ufold

		//only select protons	
		if (beamtrackPdg!=pdg) continue; //only interested in protons

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

		//Truth label of Primarytrack_End ------------------------------------------------------------------------------------------------//
		bool IsPureInEL=false; //inel
		bool IsPureEL=false; //el

		if (strcmp(primary_truth_EndProcess->c_str(),"protonInelastic")==0) { IsPureInEL=true; }
		else { IsPureEL=true; }
		//--------------------------------------------------------------------------------------------------------------------------------//

		//Get true start/end point -----------------------------------------------------------------------//
		double true_endz=primary_truth_EndPosition_MC[2]; 
		double true_endy=primary_truth_EndPosition_MC[1]; 
		double true_endx=primary_truth_EndPosition_MC[0];

		double true_stz=primary_truth_StartPosition_MC[2];
		double true_sty=primary_truth_StartPosition_MC[1];
		double true_stx=primary_truth_StartPosition_MC[0];

		bool IsTrueEndOutside=false;
		if (true_endz<0.) IsTrueEndOutside=true;

		//Get reco info ----------------------------------------------------------------------------------//
		//Evt Classification -----------------------------------------------------------------------------//
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
		//Evt Classification -----------------------------------------------------------------------------//

		//reco pos info & cut
		double reco_stx=-99, reco_sty=-99, reco_stz=-99;
		double reco_endx=-99, reco_endy=-99, reco_endz=-99;
		//double dx_reco_stx__bposx_ff=-999, dy_reco_sty__bposy_ff=-999, dz_reco_stz__bposz_ff=-999;
		bool IsPos=false;
		if (IsCaloSize) {
			reco_stx=primtrk_hitx->at(0); 
			reco_sty=primtrk_hity->at(0);
			reco_stz=primtrk_hitz->at(0);

			reco_endx=primtrk_hitx->at(primtrk_dedx->size()-1);	
			reco_endy=primtrk_hity->at(primtrk_dedx->size()-1);
			reco_endz=primtrk_hitz->at(primtrk_dedx->size()-1);

			double beam_dx=(reco_stx-mean_StartX)/sigma_StartX;
			double beam_dy=(reco_sty-mean_StartY)/sigma_StartY;
			double beam_dz=(reco_stz-mean_StartZ)/sigma_StartZ;
			double beam_dxy=sqrt(pow(beam_dx,2)+pow(beam_dy,2));	

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

		//cosine_theta/cut
		bool IsCosine=false;
		double cosine_beam_spec_primtrk=-999; 
		//cosine_beam_spec_primtrk=beamDirx_spec->at(0)*primaryStartDirection[0]+beamDiry_spec->at(0)*primaryStartDirection[1]+beamDirz_spec->at(0)*primaryStartDirection[2]; //cosine between beam_spec and primary trk direction(no SCE corr.)
		TVector3 dir;
		if (IsCaloSize) {	
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
		}

		if (cosine_beam_spec_primtrk<0) { cosine_beam_spec_primtrk=-1.*cosine_beam_spec_primtrk; }
		if (cosine_beam_spec_primtrk>cosine_beam_primtrk_min) { IsCosine=true; }

		//xy-cut
		bool IsXY=false;		
		double reco_stx_noSCE=0, reco_sty_noSCE=0, reco_stz_noSCE=0; //start-pos, before sce
		if (primaryEndPosition[2]>primaryStartPosition[2]) { //check if Pandora flip the sign
			reco_stx_noSCE=primaryStartPosition[0];
			reco_sty_noSCE=primaryStartPosition[1];
			reco_stz_noSCE=primaryStartPosition[2];
		} //check if Pandora flip the sign
		else {
			reco_stx_noSCE=primaryEndPosition[0];
			reco_sty_noSCE=primaryEndPosition[1];
			reco_stz_noSCE=primaryEndPosition[2];
		}
		if ((pow(((reco_stx_noSCE-mean_x)/dev_x),2)+pow(((reco_sty_noSCE-mean_y)/dev_y),2))<=1.) IsXY=true;

		//Intersection cut
		//bool IsIntersection=false;		
		//if (timeintersection->size()) IsIntersection=true; //over-lapping track cut

		//beam quality cut
		bool IsBQ=false;
		if (IsCosine&&IsPos) IsBQ=true;

		int index_reco_endz=0;
		double wid_reco_max=-9999;
		double range_reco=-999;
		vector<double> reco_trklen_accum;
  		reco_trklen_accum.reserve(primtrk_hitz->size());
		double kereco_calo=0;
		double kereco_range=0;
		double kereco_range2=0;
		vector<double> EDept;
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

			if (wid_reco>wid_reco_max) { 
				wid_reco_max=wid_reco;
				index_reco_endz=(int)-1+primtrk_wid->size()-h;
			}

			double cali_dedx=0.;
			cali_dedx=dedx_function_35ms(dqdx, hitx_reco, hity_reco, hitz_reco);
			EDept.push_back(cali_dedx*pitch);

			if (h==1) range_reco=0;
			if (h>=1) {
    					range_reco += sqrt( pow(primtrk_hitx->at(h)-primtrk_hitx->at(h-1), 2)+
					    		    pow(primtrk_hity->at(h)-primtrk_hity->at(h-1), 2)+
					    		    pow(primtrk_hitz->at(h)-primtrk_hitz->at(h-1), 2) );
					reco_trklen_accum[h] = range_reco;
							//TString sav_reco=Form("trklen_reco.push_back(%.4f);\n", range_reco);
							//myfile_reco<<sav_reco.Data();
			}

			kereco_calo+=cali_dedx*pitch;
			kereco_range+=pitch*dedx_predict(resrange_reco);
			kereco_range2+=pitch*(double)gr_predict_dedx_resrange->Eval(resrange_reco);

		  } //loop over reco hits of a given track


		  //range_reco=primtrk_range->at(0);
		} //if calo size not empty

		//Reco stopping/Inel p cut
		bool IsRecoStop=false;
		bool IsRecoInEL=false;
		bool IsRecoEL=false;
		double mom_beam_spec=-99; mom_beam_spec=beamMomentum_spec->at(0);
		double bx_spec=beamPosx_spec->at(0);
		double by_spec=beamPosy_spec->at(0);

		//double range_reco=-99; if (!primtrk_range->empty()) range_reco=primtrk_range->at(0); //reco primary trklen
		double csda_val_spec=csda_range_vs_mom_sm->Eval(mom_beam_spec);
		double ntrklen_reco=range_reco/csda_val_spec;

		if (ntrklen_reco>=min_norm_trklen_csda&&ntrklen_reco<max_norm_trklen_csda) IsRecoStop=true;
		//if ((range_reco/csda_val_spec)<min_norm_trklen_csda) IsRecoInEL=true;
		//
		


                //PID parameter & cut
                double pid=-99;
                double median_dedx=-99;
                //if (IsBQ&&IsCaloSize&&IsPandoraSlice) { //beam quality cut
                if (IsCaloSize&&IsPandoraSlice) { //beam quality cut
                        //calo
                        vector<double> trkdedx;
                        vector<double> trkres;
                        for (size_t h=0; h<primtrk_hitz->size(); ++h) { //loop over reco hits of a given track
                                double hitx_reco=primtrk_hitx->at(h);
                                double hity_reco=primtrk_hity->at(h);
                                double hitz_reco=primtrk_hitz->at(h);
                                double rr_reco=primtrk_resrange->at(h);

                                double dqdx=primtrk_dqdx->at(h);
                                double pitch=primtrk_pitch->at(h);

                                int wid_reco=primtrk_wid->at(-1+primtrk_wid->size()-h);
                                double pt_reco=primtrk_pt->at(-1+primtrk_wid->size()-h);

                                double dedx_reco=0.;
                                dedx_reco=dedx_function_35ms(dqdx, hitx_reco, hity_reco, hitz_reco);
                                //ke_calo_MeV+=cali_dedx*pitch; //prod4a, pandoracalinoxyzt

                                trkdedx.push_back(dedx_reco);
                                trkres.push_back(rr_reco);

				double r0=pow(dedx_reco,0.42)/17.-rr_reco;
				double dedx_pred=17.*pow(rr_reco,-0.42);

				//cout<<"z:"<<hitz_reco<<" dedx_reco:"<<dedx_reco<<" dedx_pred:"<<dedx_pred<<" rr_reco:"<<rr_reco<<" r0:"<<r0<<" len="<<range_reco-rr_reco<<endl;
				//myfile_reco<<"dedx_reco.push_back("<<dedx_reco<<"); rr_reco.push_back("<<rr_reco<<"); \n";
				

                        } //loop over reco hits of a given track
                        pid=chi2pid(trkdedx,trkres); //pid using stopping proton hypothesis

                        median_dedx=TMath::Median(trkdedx.size(), &trkdedx.at(0));
                } //beam quality cut

		if ((range_reco/csda_val_spec)<min_norm_trklen_csda) { //inel region
			if (pid>pid_1) IsRecoInEL=true; 
			if (pid<=pid_1) IsRecoEL=true; 
		} //inel region
		if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) { //stopping p region
			if (pid>pid_2) IsRecoInEL=true; 
			if (pid<=pid_2) IsRecoEL=true;
		} //stopping p region


		//kinetic energies ---------------------------------------------------------------//
		//double ke_beam=1000.*p2ke(mom_beam); //ke_beam
		double ke_beam_spec=p2ke(mom_beam_spec); //ke_beam_spec [GeV]
		double ke_beam_spec_MeV=1000.*ke_beam_spec; //ke_beam_spec [MeV]
		double ke_trklen=ke_vs_csda_range_sm->Eval(range_reco); //[unit: GeV]
		double ke_trklen_MeV=1000.*ke_trklen; //[unit: MeV]
		double ke0=-999;

		bool IsBeamMom=false; //apply 3-sigma cut and remove tail events
		if (ke_beam_spec_MeV>=(mean_kebeam-3.*sigma_kebeam)&&ke_beam_spec_MeV<=(mean_kebeam+3.*sigma_kebeam)) IsBeamMom=true;
		if (IsBeamMom==false) continue; //if not within 3-sigma, skip the event

		//First point of MCParticle entering TPC ------------------------------------------------------------------//
		bool is_beam_at_ff=false;
		int key_reach_tpc=-99;
		if (beamtrk_z->size()){
			ke0=1000.*(beamtrk_Eng->at(0));
			for (size_t kk=0; kk<beamtrk_z->size(); ++kk) {  //loop over all beam hits
				double zpos_beam=beamtrk_z->at(kk);
				if (zpos_beam>=0) {
					key_reach_tpc=(int)kk;
					break;
				}
			} //loop over all beam hits

			//for (size_t kk=0; kk<beamtrk_z->size(); ++kk) {  //loop over all beam hits
			//cout<<"["<<kk<<"] beamtrk_z:"<<beamtrk_z->at(kk) <<" beamtrk_Eng:"<<beamtrk_Eng->at(kk)<<endl;
			//} //loop over all beam hits
		} 
		if (key_reach_tpc!=-99) { is_beam_at_ff=true; }
		cout<<"key_reach_tpc:"<<key_reach_tpc<<endl;	
		cout<<"is_beam_at_ff:"<<is_beam_at_ff<<endl;

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


		double KE_ff_tpc=-999;
		double KE_ff=0;

		if (is_beam_at_ff) {
			KE_ff=ke_ff; //use KE exactly at z=0 
			KE_ff_tpc=1000.*(beamtrk_Eng->at(key_reach_tpc)); //MeV
		}

		//double p_trklen=ke2p(ke_trklen);
		//double ke_simide=0;
		//for (int hk=0; hk<(int)primtrk_true_edept->size(); ++hk) { //loop over simIDE points
		//ke_simide+=primtrk_true_edept->at(hk);
		//} //loop over simIDE points

		//fix on the truth length by adding distance between 1st tpc hit to front face ------------------------------------------------------//
		//[1] 3D projection on TPC front face
		double zproj_beam=0; //set beam z at ff
		double yproj_beam=0; //ini. value
		double xproj_beam=0; //ini. value
		int n_fit=3; //num of points used for fitting
                if (beamtrk_z->size()) {

			int key_fit_st=0;
			int key_fit_ed=-1+(int)beamtrk_z->size();
			if (key_reach_tpc!=-99) {
				key_fit_st=key_reach_tpc-1;
				key_fit_ed=key_reach_tpc+1;
			}
			if (key_fit_st<0) key_fit_st=0;
			if (key_fit_ed>(-1+(int)beamtrk_z->size())) key_fit_ed=-1+(int)beamtrk_z->size();	

			cout<<"beamtrk_z->size():"<<beamtrk_z->size()<<endl;
			cout<<"key_reach_tpc:"<<key_reach_tpc<<endl;
			std::cout<<"key_fit_st-ed:"<<key_fit_st<<"-"<<key_fit_ed<<std::endl;

			//start 3D line fit
			TGraph2D *gr=new TGraph2D();
			//cout<<"ck0"<<endl;
   		  	//for (int N=key_fit_st; N<key_fit_ed; N++) {
   		  	int nsize_fit=n_fit;
			if ((1+(key_fit_ed-key_fit_st))<n_fit) nsize_fit=1+(key_fit_ed-key_fit_st);
			if ((int)beamtrk_z->size()<=n_fit) nsize_fit=(int)beamtrk_z->size(); //in case really short track
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
			if (!ok) {
				Error("line3Dfit","Line3D Fit failed");
				//return 1;
			}
			//cout<<"ck5"<<endl;
				
			const ROOT::Fit::FitResult & result = fitter.Result();
			std::cout << "Total final distance square " << result.MinFcnValue() << std::endl;
			result.Print(std::cout);
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

				
							//double len_true=true_trklen_accum.at(mm);
							//TString sav_true=Form("trklen_true.push_back(%.4f);\n",range_true);
							//myfile_true<<sav_true.Data();
							//myfile_reco<<sav_true.Data();

		  }


		} //is beam at ff
		double ke_truetrklen=ke_vs_csda_range_sm->Eval(range_true); //[unit: GeV]
		double ke_truetrklen_MeV=1000.*ke_truetrklen; //[unit: MeV]
						    	
		//fix on the truth length by adding distance between 1st tpc hit to front face ------------------------------------------------------//

		//KEs ---------------------------------------------------------------------------------------------------//
		//energy loss using stopping protons
		double mean_Eloss_upstream=19.3073; //unit:MeV
		double err_mean_Eloss_upstream=0.187143;
		double sigma_Eloss_upstream=18.7378;
		double err_sigma_Eloss_upstream=0.140183;

		double mean_Elossrange_stop=433.441-405.371; //KEstop using range, unit:MeV
		//double mean_Elosscalo_stop=433.441-379.074; //
		//double mean_Elosscalo_stop=(4.77927e+01)/(1.00480e+00); //using fit [bmrw]
		double mean_Elosscalo_stop=(4.95958e+01)/(1.00489e+00); //using fit [no bmrw]
		//fit result
   		//p0           4.77927e+01   3.62629e-01   1.40469e-08   0.00000e+00
   		//p1          -1.00480e+00   7.70658e-03   7.70658e-03  -4.73181e-08

   		//p0           4.95958e+01   2.69311e-01   3.68850e-08   5.23619e-11
   		//p1          -1.00489e+00   5.73159e-03   5.73159e-03   1.03334e-07


		double mean_ntrklen=9.00928e-01;
		double sigma_ntrklen=7.61209e-02;



		double Eloss_upstream=0; 
		double Eloss_upstream_reco=0;
		double dKE_recotruth_1=0;
		double dKE_recotruth_2=0;
		if (is_beam_at_ff) { //if beam reach tpc
			Eloss_upstream=ke_beam_spec_MeV-KE_ff;
			Eloss_upstream_reco=ke_beam_spec_MeV-mean_Eloss_upstream;
		} //if beam reach tpc
		else { //if beam NOT reach tpc
		} //if beam NOT reach TPC

		//double ke_trklen=ke_vs_csda_range_sm->Eval(range_reco); //[unit: GeV]
		//double ke_trklen_MeV=1000.*ke_trklen; //[unit: MeV]

		//double dEbb_true=0; if (range_true>=0&&is_beam_at_ff) dEbb_true=BB.KEFromRangeSpline(range_true);
		//double dEbb_reco=0; if (range_reco>=0&&is_beam_at_ff) dEbb_reco=BB.KEFromRangeSpline(range_reco);
		//double dEbb=0; dEbb=dEbb_true-dEbb_reco;

		//mean and sigma of e-loss [mc]
		//double mean_Eloss_upstream=2.04949e+01;
		//double err_mean_Eloss_upstream=1.47946e-01;
		//double sigma_Eloss_upstream=1.90634e+01;
		//double err_sigma_Eloss_upstream=1.23169e-01;

		
		//double keff_reco=ke_beam_spec_MeV-mean_Eloss_upstream;


		double ke_calo_MeV=0;
		//double av_dedx=0;     int n_dedx=0;
		//double av_dedx_all=0; int n_dedx_all=0;
		vector<double> trkdedx; 
		vector<double> trkres;
		//vector<double> R0;
		double av_r0=0;
		int n_r0=0;
		if (IsPandoraSlice&&IsBQ&&IsCaloSize) { //if calo size not empty
			//calo
			for (size_t h=0; h<primtrk_hitz->size(); ++h) { //loop over reco hits of a given track
				double hitx_reco=primtrk_hitx->at(h);
				double hity_reco=primtrk_hity->at(h);
				double hitz_reco=primtrk_hitz->at(h);
				double resrange_reco=primtrk_resrange->at(h);

				double dqdx=primtrk_dqdx->at(h);
				double pitch=primtrk_pitch->at(h);

				int wid_reco=primtrk_wid->at(-1+primtrk_wid->size()-h);
				double pt_reco=primtrk_pt->at(-1+primtrk_wid->size()-h);

				double cali_dedx=0.;
				cali_dedx=dedx_function_35ms(dqdx, hitx_reco, hity_reco, hitz_reco);
				ke_calo_MeV+=cali_dedx*pitch; //prod4a, pandoracalinoxyzt

				trkdedx.push_back(cali_dedx);
				trkres.push_back(resrange_reco);

				double tmp_r0=pow((17./cali_dedx),1./0.42)-resrange_reco;
				//R0.push_back(tmp_r0);
				av_r0+=tmp_r0;
				n_r0++;

				//average dE/dx ----------------------------------//
				//if (cali_dedx<30.&&cali_dedx>1.) { 
				//if (h==-1+(primtrk_hitz->size())) { 
					//av_dedx+=cali_dedx; n_dedx++;
				//}
				//av_dedx_all+=cali_dedx; n_dedx_all++;
				//average dE/dx ----------------------------------//

				//if (IsRecoStop) {
					//h2d_rr_dedx_stop->Fill(resrange_reco, cali_dedx);
					//h2d_dedx_rr_stop->Fill(cali_dedx, resrange_reco);
				//}

				//double len_true=true_trklen_accum.at(mm);
			} //loop over reco hits of a given track

		} //if calo size not empty
		av_r0/=(double)n_r0;

		//KEff for el
        	double pff_truth=ke2p(ke_ff/1000.);
        	double csda_truth=csda_range_vs_mom_sm->Eval(pff_truth);
		double KEcsda=1000.*ke_vs_csda_range_sm->Eval(range_reco);

		//keff for inel --------------------------------------------------//
/*
		vector<double> r0_sg;
		vector<double> res_r0;
		double av_res_r0=0;
		double av2_res_r0=0;
		double rms_res_r0=0;
		double up_rms_res_r0=1000;
		double dn_rms_res_r0=-1000;
		double r0_predict=0;
		int n_predict=0;
		if (IsPandoraSlice&&IsBQ&&IsCaloSize) { //if calo size not empty
			//low-pass filter on r0
  			r0_sg=sg_smooth(R0, 3, 3); //window size, pol order

			for (size_t jj=0; jj<R0.size(); ++jj) {
				double tmp_res_r0=R0.at(jj)-r0_sg.at(jj);
				av_res_r0+=tmp_res_r0;
				av2_res_r0+=tmp_res_r0*tmp_res_r0;
				res_r0.push_back(tmp_res_r0);
			}
			if (R0.size()) { 
				av_res_r0/=(double)R0.size();
				av2_res_r0/=(double)R0.size();
			}
			rms_res_r0=sqrt(av2_res_r0-av_res_r0*av_res_r0);
			up_rms_res_r0=av_res_r0+rms_res_r0;
			dn_rms_res_r0=av_res_r0-rms_res_r0;

			for (size_t jj=0; jj<R0.size(); ++jj) {
				double tmp_res_r0=R0.at(jj)-r0_sg.at(jj);
				if (tmp_res_r0>dn_rms_res_r0&&tmp_res_r0<up_rms_res_r0) {
					r0_predict+=R0.at(jj);
					n_predict++;
				}
			}
			if (n_predict) r0_predict/=n_predict;

		} //if calo size not empty

		double KEcsda_inel=ke_beam_spec_MeV-mean_Elossrange_stop;
		if (n_predict) {
			KEcsda_inel=1000.*ke_vs_csda_range_sm->Eval(r0_predict+range_reco);
		}
*/
		//double r0_reco=(mean_ntrklen-ntrklen_reco)*csda_val_spec;
		//double KEcsda_inel=1000.*ke_vs_csda_range_sm->Eval(r0_reco+range_reco);
		//if (ntrklen_reco>mean_ntrklen) KEcsda_inel=ke_beam_spec_MeV-mean_Elossrange_stop;
		//
		//double KEcsda_inel=ke_beam_spec_MeV-mean_Elossrange_stop;
		//if (av_r0>0) KEcsda_inel=1000.*ke_vs_csda_range_sm->Eval(av_r0+range_reco);
		//keff for inel --------------------------------------------------//
		
		//double KEcsda_inel=1000.*ke_vs_csda_range_sm->Eval(csda_truth); //given precise KEff, cossesponding CSDA range
		//double Edept_range=KEcsda_inel;


		//KEend -------------------------------------------------------------------------------------------------------------------//
		double KEend_true=1000.*(beamtrk_Eng->at(-2+beamtrk_Eng->size()));
		//double KEbb_true=-1; KEbb_true=BB.KEAtLength(ke_ff, range_true);

		double KEff_reco=ke_beam_spec_MeV-mean_Elossrange_stop;
		//double KEbb_reco_constErange=-1; KEbb_reco_constErange=BB.KEAtLength(KEff_reco, range_reco);
		//double KEbb_reco_Edept_range=-1; KEbb_reco_Edept_range=BB.KEAtLength(KEcsda, range_reco);
		//double KEbb_reco_Edept_range_inel=-1; KEbb_reco_Edept_range_inel=BB.KEAtLength(KEcsda_inel, range_reco);

		//double KEcalo_reco_constEcalo=-1; KEcalo_reco_constEcalo=ke_beam_spec_MeV-mean_Elosscalo_stop-ke_calo_MeV;
		//double KEcalo_reco_Edept_el=-1; KEcalo_reco_Edept_el=KEcsda-ke_calo_MeV;
		//double KEcalo_reco_Edept_inel=-1; KEcalo_reco_Edept_inel=KEcsda_inel-ke_calo_MeV;

/*
		//av dE/dx ----------------------------------------------------------------//
		if (n_dedx>0) av_dedx/=(double)n_dedx;
		//if (n_dedx_all>0) av_dedx_all/=(double)n_dedx_all;
		double av_rr_reco=0;
		if (n_dedx!=0) av_rr_reco=pow((17./av_dedx),1./0.42);
		//else av_rr_reco=pow((17./av_dedx_all),1./0.42);

		double Aint=0.58;
		double Edept_range=(17./Aint)*(pow(av_rr_reco, Aint));
		//double keff_reco5=ke_calo_MeV+Edept_range;
		//double KEbb_reco5=-1; KEbb_reco5=BB.KEAtLength(keff_reco5, range_reco);
		//av dE/dx ----------------------------------------------------------------//
*/



		if (IsPandoraSlice&&IsBQ&&IsCaloSize) { //basic cuts
			Fill1DHist(h1d_kebeam, ke_beam_spec_MeV);
			Fill1DHist(h1d_ke0, ke0);
			Fill1DHist(h1d_keff_truth, ke_ff);
			Fill1DHist(h1d_keff_reco_constErange, KEff_reco);
			Fill1DHist(h1d_ntrklen, ntrklen_reco);

			double mom_rw_minchi2=1;
			if ((mom_beam_spec*1000.)>=mmu_min&&(mom_beam_spec*1000.)<=mmu_max) mom_rw_minchi2=gng[index_minchi2]->Eval(mom_beam_spec*1000.); //bmrw(use me!)
			//if ((ke_beam_spec_MeV-mean_Elosscalo_stop)>=mu_min&&(ke_beam_spec_MeV-mean_Elosscalo_stop)<=mu_max) mom_rw_minchi2=gng->Eval(ke_beam_spec_MeV-mean_Elosscalo_stop); //bmrw

			//h1d_keff_reco_constErange->Fill(KEff_reco, mom_rw_minchi2);
			//h2d_trklen_KEbb->Fill(range_reco, KEbb_reco_constErange, mom_rw_minchi2);
			//h1d_kebb_reco_constErange->Fill(KEbb_reco_constErange, mom_rw_minchi2);
			//h1d_Edept_reco->Fill(KEcsda, mom_rw_minchi2);

			for (int k=0; k<n_eloss; ++k) {
				h1d_KEend_tune_bmrw[k]->Fill((ke_beam_spec_MeV-Eloss.at(k)-ke_calo_MeV), mom_rw_minchi2);
			}

			if (IsRecoInEL)	{ //IsRecoInEL
				//h1d_kebb_reco->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				//h1d_kebb_recoInel->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				//h2d_trklen_KEbb_RecoInel->Fill(range_reco, KEbb_reco_constErange, mom_rw_minchi2);

/*
				if (kinel) {
					h1d_kebb_reco_inel->Fill(KEbb_reco_constErange, mom_rw_minchi2);
					h1d_kebb_recoInel_inel->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				}
				if (kel) {
					h1d_kebb_reco_el->Fill(KEbb_reco_constErange, mom_rw_minchi2);
					h1d_kebb_recoInel_el->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				}
				if (kMIDcosmic) {
					h1d_kebb_reco_midcosmic->Fill(KEbb_reco_constErange, mom_rw_minchi2);
					h1d_kebb_recoInel_midcosmic->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				}
				if (kMIDpi) {
					h1d_kebb_reco_midpi->Fill(KEbb_reco_constErange, mom_rw_minchi2);
					h1d_kebb_recoInel_midpi->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				}
				if (kMIDp) {
					h1d_kebb_reco_midp->Fill(KEbb_reco_constErange, mom_rw_minchi2);
					h1d_kebb_recoInel_midp->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				}
				if (kMIDmu) {
					h1d_kebb_reco_midmu->Fill(KEbb_reco_constErange, mom_rw_minchi2);
					h1d_kebb_recoInel_midmu->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				}
				if (kMIDeg) {
					h1d_kebb_reco_mideg->Fill(KEbb_reco_constErange, mom_rw_minchi2);
					h1d_kebb_recoInel_mideg->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				}
				if (kMIDother) {
					h1d_kebb_reco_midother->Fill(KEbb_reco_constErange, mom_rw_minchi2);
					h1d_kebb_recoInel_midother->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				}
*/

/*

				h1d_kebb_recoInel_constEcalo->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				if (kinel) {
					h1d_kebb_recoInel_constEcalo_inel->Fill(KEbb_reco_constErange, mom_rw_minchi2);
					h1d_Edept_reco_inel->Fill(KEcsda, mom_rw_minchi2);
				}
				if (kel) {
					h1d_kebb_recoInel_constEcalo_el->Fill(KEbb_reco_constErange, mom_rw_minchi2);
					h1d_Edept_reco_el->Fill(KEcsda, mom_rw_minchi2);
				}
				if (kMIDcosmic) {
					h1d_kebb_recoInel_constEcalo_midcosmic->Fill(KEbb_reco_constErange, mom_rw_minchi2);
					h1d_Edept_reco_midcosmic->Fill(KEcsda, mom_rw_minchi2);
				}
				if (kMIDpi) {
					h1d_kebb_recoInel_constEcalo_midpi->Fill(KEbb_reco_constErange, mom_rw_minchi2);
					h1d_Edept_reco_midpi->Fill(KEcsda, mom_rw_minchi2);
				}
				if (kMIDp) {
					h1d_kebb_recoInel_constEcalo_midp->Fill(KEbb_reco_constErange, mom_rw_minchi2);
					h1d_Edept_reco_midp->Fill(KEcsda, mom_rw_minchi2);
				}
				if (kMIDmu) {
					h1d_kebb_recoInel_constEcalo_midmu->Fill(KEbb_reco_constErange, mom_rw_minchi2);
					h1d_Edept_reco_midmu->Fill(KEcsda, mom_rw_minchi2);
				}
				if (kMIDeg) {
					h1d_kebb_recoInel_constEcalo_mideg->Fill(KEbb_reco_constErange, mom_rw_minchi2);
					h1d_Edept_reco_mideg->Fill(KEcsda, mom_rw_minchi2);
				}
				if (kMIDother) {
					h1d_kebb_recoInel_constEcalo_midother->Fill(KEbb_reco_constErange, mom_rw_minchi2);
					h1d_Edept_reco_midother->Fill(KEcsda, mom_rw_minchi2);
				}
*/


			} //IsRecoInEL

			if (IsRecoEL) { //IsRecoEL
				//h1d_kebb_reco->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
				//h1d_kebb_recoEl->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
				//h2d_trklen_KEbb_RecoEl->Fill(range_reco, KEbb_reco_constErange, mom_rw_minchi2);

/*
				if (kinel) {
					h1d_kebb_reco_inel->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
					h1d_kebb_recoEl_inel->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
				}
				if (kel) {
					h1d_kebb_reco_el->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
					h1d_kebb_recoEl_el->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
				}
				if (kMIDcosmic) {
					h1d_kebb_reco_midcosmic->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
					h1d_kebb_recoEl_midcosmic->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
				}
				if (kMIDpi) {
					h1d_kebb_reco_midpi->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
					h1d_kebb_recoEl_midpi->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
				}
				if (kMIDp) {
					h1d_kebb_reco_midp->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
					h1d_kebb_recoEl_midp->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
				}
				if (kMIDmu) {
					h1d_kebb_reco_midmu->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
					h1d_kebb_recoEl_midmu->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
				}
				if (kMIDeg) {
					h1d_kebb_reco_mideg->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
					h1d_kebb_recoEl_mideg->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
				}
				if (kMIDother) {
					h1d_kebb_reco_midother->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
					h1d_kebb_recoEl_midother->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
				}
*/

/*

				h1d_kebb_recoEl_constEcalo->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				if (kinel) {
					h1d_kebb_recoEl_constEcalo_inel->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				}
				if (kel) {
					h1d_kebb_recoEl_constEcalo_el->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				}
				if (kMIDcosmic) {
					h1d_kebb_recoEl_constEcalo_midcosmic->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				}
				if (kMIDpi) {
					h1d_kebb_recoEl_constEcalo_midpi->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				}
				if (kMIDp) {
					h1d_kebb_recoEl_constEcalo_midp->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				}
				if (kMIDmu) {
					h1d_kebb_recoEl_constEcalo_midmu->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				}
				if (kMIDeg) {
					h1d_kebb_recoEl_constEcalo_mideg->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				}
				if (kMIDother) {
					h1d_kebb_recoEl_constEcalo_midother->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				}

*/

			} //IsRecoEL




/*

			if (kel) { 

				Fill1DHist(h1d_keff_truth_el, ke_ff);
				Fill1DHist(h1d_keff_reco_constErange_el, ke_beam_spec_MeV-mean_Elossrange_stop);
				h1d_keff_reco_constErange_el->Fill(ke_beam_spec_MeV-mean_Elosscalo_stop, mom_rw_minchi2);

				Fill1DHist(h1d_keff_reco_Edept_range_el, KEcsda); 
				Fill1DHist(h1d_keff_reco_Edept_calo_el, KEcsda);

				Fill1DHist(h1d_ntrklen_el, ntrklen_reco);

				Fill1DHist(h1d_keend_truth_el, KEend_true);
				//Fill1DHist(h1d_kebb_truth_el, KEbb_true);

				//Fill1DHist(h1d_kebb_reco_constErange_el, KEbb_reco_constErange);
				//Fill1DHist(h1d_kebb_reco_Edept_range_el, KEbb_reco_Edept_range);

				h1d_kebb_reco_constErange_el->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				Fill1DHist(h1d_kebb_reco_Edept_calo_el, KEbb_reco_Edept_el); 
				
	
			}


			if (kinel) { 
				Fill1DHist(h1d_keff_truth_inel, ke_ff);
				Fill1DHist(h1d_keff_reco_constErange_inel, ke_beam_spec_MeV-mean_Elossrange_stop);
				h1d_keff_reco_constErange_inel->Fill(ke_beam_spec_MeV-mean_Elosscalo_stop, mom_rw_minchi2);

				Fill1DHist(h1d_keff_reco_Edept_range_inel, KEcsda_inel); 
				//Fill1DHist(h1d_keff_reco_Edept_calo_inel, ke_calo_MeV+Edept_range);

				Fill1DHist(h1d_ntrklen_inel, ntrklen_reco);

				Fill1DHist(h1d_keend_truth_inel, KEend_true);
				//Fill1DHist(h1d_kebb_truth_inel, KEbb_true);

				//Fill1DHist(h1d_kebb_reco_constErange_inel, KEbb_reco_constErange);
				//Fill1DHist(h1d_kebb_reco_Edept_range_inel, KEbb_reco_Edept_range_inel);

				h1d_kebb_reco_constErange_inel->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				//Fill1DHist(h1d_kebb_reco_Edept_calo_inel, KEbb_reco_Edept_inel); 

 
				//h2d_keff_r0ratio_inel->Fill(ke_ff, (av_r0/csda_truth));
				//h2d_ntrklen_r0ratio_inel->Fill(ntrklen_reco, (av_r0/csda_truth));

			}

			if (kMIDcosmic) {
				h1d_kebb_reco_constErange_midcosmic->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				h1d_keff_reco_constErange_midcosmic->Fill(ke_beam_spec_MeV-mean_Elosscalo_stop, mom_rw_minchi2);
			}
			if (kMIDpi) {
				h1d_kebb_reco_constErange_midpi->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				h1d_keff_reco_constErange_midpi->Fill(ke_beam_spec_MeV-mean_Elosscalo_stop, mom_rw_minchi2);
			}
			if (kMIDp) {
				h1d_kebb_reco_constErange_midp->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				h1d_keff_reco_constErange_midp->Fill(ke_beam_spec_MeV-mean_Elosscalo_stop, mom_rw_minchi2);
			}
			if (kMIDmu) {
				h1d_kebb_reco_constErange_midmu->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				h1d_keff_reco_constErange_midmu->Fill(ke_beam_spec_MeV-mean_Elosscalo_stop, mom_rw_minchi2);
			}
			if (kMIDeg) {
				h1d_kebb_reco_constErange_mideg->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				h1d_keff_reco_constErange_mideg->Fill(ke_beam_spec_MeV-mean_Elosscalo_stop, mom_rw_minchi2);
			}
			if (kMIDother) {
				h1d_kebb_reco_constErange_midother->Fill(KEbb_reco_constErange, mom_rw_minchi2);
				h1d_keff_reco_constErange_midother->Fill(ke_beam_spec_MeV-mean_Elosscalo_stop, mom_rw_minchi2);
			}


			if (IsRecoStop) {
				Fill1DHist(h1d_ntrklen_recostop, ntrklen_reco);
			}
*/

		} //basic cuts

	} //main entry loop


	//save results...
   	//TFile *fout = new TFile("mc_ke_advancedKEff_nobmrw.root","RECREATE");
   	//TFile *fout = new TFile("mc_kebb_bmrw.root","RECREATE");
   	//TFile *fout = new TFile("mc_kebb_nobmrw_ElossTune.root","RECREATE");
   	//TFile *fout = new TFile("mc_kebb_nobmrw.root","RECREATE");
   	//TFile *fout = new TFile("mc_kebb_bmrw.root","RECREATE");
   	//TFile *fout = new TFile("mc_kebb_nobmrw.root","RECREATE");
   	//TFile *fout = new TFile("mc_kecalo_bmrw_ElossTune.root","RECREATE");
   	TFile *fout = new TFile("mc_kecalo_bmrw_beamxy_ElossTune.root","RECREATE");

/*
		h1d_kebb_truth_el->Write();
		h1d_kebb_truth_inel->Write();

		h1d_kebb_reco_constErange_el->Write();
		h1d_kebb_reco_constErange_inel->Write();

		h1d_kebb_reco_Edept_range_el->Write();
		h1d_kebb_reco_Edept_range_inel->Write();
*/


/*
		h1d_kebeam->Write();
		h1d_ke0->Write();

		h1d_keff_truth->Write();
		h1d_keff_truth_el->Write();
		h1d_keff_truth_inel->Write();

		h1d_keff_reco_constErange->Write();
		h1d_keff_reco_constErange_el->Write();
		h1d_keff_reco_constErange_inel->Write();

		h1d_keff_reco_constErange->Write();
		h1d_keff_reco_constErange_el->Write();
		h1d_keff_reco_constErange_inel->Write();
		h1d_keff_reco_constErange_midcosmic->Write();
		h1d_keff_reco_constErange_midpi->Write();
		h1d_keff_reco_constErange_midp->Write();
		h1d_keff_reco_constErange_midmu->Write();
		h1d_keff_reco_constErange_mideg->Write();
		h1d_keff_reco_constErange_midother->Write();

		h1d_keff_reco_Edept_range_el->Write();
		h1d_keff_reco_Edept_calo_el->Write();

		h1d_keff_reco_Edept_range_inel->Write();
		//h1d_keff_reco_Edept_calo_inel->Write();


		h1d_ntrklen->Write();
		h1d_ntrklen_recostop->Write();
		h1d_ntrklen_el->Write();
		h1d_ntrklen_inel->Write();

		//h2d_keff_r0ratio_inel->Write();
		//h2d_ntrklen_r0ratio_inel->Write();

		h1d_keend_truth_el->Write();
		h1d_keend_truth_inel->Write();


		h1d_kebb_reco_constErange->Write();
		h1d_kebb_reco_constErange_el->Write();
		h1d_kebb_reco_constErange_inel->Write();
		h1d_kebb_reco_constErange_midcosmic->Write();
		h1d_kebb_reco_constErange_midpi->Write();
		h1d_kebb_reco_constErange_midp->Write();
		h1d_kebb_reco_constErange_midmu->Write();
		h1d_kebb_reco_constErange_mideg->Write();
		h1d_kebb_reco_constErange_midother->Write();



		h1d_kebb_reco_Edept_calo_el->Write();
		//h1d_kebb_reco_Edept_calo_inel->Write();
*/

/*
		h1d_kebb_reco->Write();
		h1d_kebb_reco_inel->Write();
		h1d_kebb_reco_el->Write();
		h1d_kebb_reco_midcosmic->Write();
		h1d_kebb_reco_midpi->Write();
		h1d_kebb_reco_midp->Write();
		h1d_kebb_reco_midmu->Write();
		h1d_kebb_reco_mideg->Write();
		h1d_kebb_reco_midother->Write();

		h1d_kebb_recoInel->Write();
		h1d_kebb_recoInel_inel->Write();
		h1d_kebb_recoInel_el->Write();
		h1d_kebb_recoInel_midcosmic->Write();
		h1d_kebb_recoInel_midpi->Write();
		h1d_kebb_recoInel_midp->Write();
		h1d_kebb_recoInel_midmu->Write();
		h1d_kebb_recoInel_mideg->Write();
		h1d_kebb_recoInel_midother->Write();

		h1d_kebb_recoEl->Write();
		h1d_kebb_recoEl_inel->Write();
		h1d_kebb_recoEl_el->Write();
		h1d_kebb_recoEl_midcosmic->Write();
		h1d_kebb_recoEl_midpi->Write();
		h1d_kebb_recoEl_midp->Write();
		h1d_kebb_recoEl_midmu->Write();
		h1d_kebb_recoEl_mideg->Write();
		h1d_kebb_recoEl_midother->Write();
*/

/*
		h2d_trklen_KEbb->Write();
		h2d_trklen_KEbb_RecoInel->Write();
		h2d_trklen_KEbb_RecoEl->Write();


		h1d_kebb_recoInel_constEcalo->Write();
		h1d_kebb_recoInel_constEcalo_inel->Write();
		h1d_kebb_recoInel_constEcalo_el->Write();
		h1d_kebb_recoInel_constEcalo_midcosmic->Write();
		h1d_kebb_recoInel_constEcalo_midpi->Write();
		h1d_kebb_recoInel_constEcalo_midp->Write();
		h1d_kebb_recoInel_constEcalo_midmu->Write();
		h1d_kebb_recoInel_constEcalo_mideg->Write();
		h1d_kebb_recoInel_constEcalo_midother->Write();


		h1d_kebb_recoEl_constEcalo->Write();
		h1d_kebb_recoEl_constEcalo_inel->Write();
		h1d_kebb_recoEl_constEcalo_el->Write();
		h1d_kebb_recoEl_constEcalo_midcosmic->Write();
		h1d_kebb_recoEl_constEcalo_midpi->Write();
		h1d_kebb_recoEl_constEcalo_midp->Write();
		h1d_kebb_recoEl_constEcalo_midmu->Write();
		h1d_kebb_recoEl_constEcalo_mideg->Write();
		h1d_kebb_recoEl_constEcalo_midother->Write();
*/

		for (int k=0; k<n_eloss; ++k) {
			h1d_KEend_tune_bmrw[k]->Write();
		}

/*
		h1d_Edept_reco->Write();
		h1d_Edept_reco_el->Write();
		h1d_Edept_reco_inel->Write();
		h1d_Edept_reco_midcosmic->Write();
		h1d_Edept_reco_midpi->Write();
		h1d_Edept_reco_midp->Write();
		h1d_Edept_reco_midmu->Write();
		h1d_Edept_reco_mideg->Write();
		h1d_Edept_reco_midother->Write();
*/

	fout->Close();



}
