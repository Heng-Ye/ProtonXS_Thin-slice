#define ProtonKE_cxx
#include "ProtonKE.h"

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


void ProtonKE::Loop() {
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	Long64_t nbytes = 0, nb = 0;
	bool isTestSample=true;
	//int true_sliceID = -1, reco_sliceID = -1;

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
	BetheBloch BB;
	BB.SetPdgCode(pdg);

	//Kinetic energies -------------------------------------------------------------------------------------------------------//
	int nke=160;
	double kemin=0;
	double kemax=800;

	//beam
	TH1D *h1d_ke0=new TH1D("h1d_ke0","", nke, kemin, kemax); //ke from ke_truth0
	TH1D *h1d_kebeam=new TH1D("h1d_kebeam","", nke, kemin, kemax); //ke from beamline inst.
	TH1D *h1d_kebeam_bmrw=new TH1D("h1d_kebeam_bmrw","", nke, kemin, kemax); //ke from beamline inst.
	TH1D *h1d_kebeam_stop=new TH1D("h1d_kebeam_stop","", nke, kemin, kemax); //ke from beamline inst. (stopping protons)
	TH1D *h1d_kebeam_stop_bmrw=new TH1D("h1d_kebeam_stop_bmrw","", nke, kemin, kemax); //ke from beamline inst. (stopping protons)

	//ff [truth info] 
	TH1D *h1d_keff=new TH1D("h1d_keff","", nke, kemin, kemax); //keff at z=0
	TH1D *h1d_keff2=new TH1D("h1d_keff2","", nke, kemin, kemax); //keff=ke_beam-const. E-loss
	TH1D *h1d_keff0=new TH1D("h1d_keff0","", nke, kemin, kemax); //keff at tpc entrance

	TH1D *h1d_keff_inel=new TH1D("h1d_keff_inel","", nke, kemin, kemax); 
	TH1D *h1d_keff_el=new TH1D("h1d_keff_el","", nke, kemin, kemax); 
	TH1D *h1d_keff_midcosmic=new TH1D("h1d_keff_midcosmic","", nke, kemin, kemax); 
	TH1D *h1d_keff_midpi=new TH1D("h1d_keff_midpi","", nke, kemin, kemax); 
	TH1D *h1d_keff_midp=new TH1D("h1d_keff_midp","", nke, kemin, kemax); 
	TH1D *h1d_keff_midmu=new TH1D("h1d_keff_midmu","", nke, kemin, kemax); 
	TH1D *h1d_keff_mideg=new TH1D("h1d_keff_mideg","", nke, kemin, kemax); 
	TH1D *h1d_keff_midother=new TH1D("h1d_keff_midother","", nke, kemin, kemax); 

	TH1D *h1d_keff_reco=new TH1D("h1d_keff_reco","", nke, kemin, kemax); //keff at z=0


	TH1D *h1d_keff_stop=new TH1D("h1d_keff_stop","", nke, kemin, kemax); //ke at ff (stopping protons)
	TH1D *h1d_keff2_stop=new TH1D("h1d_keff2_stop","", nke, kemin, kemax); //ke at ff (stopping protons)
	TH1D *h1d_keff0_stop=new TH1D("h1d_keff0_stop","", nke, kemin, kemax); //ke at ff (stopping protons)


	//ke of stopping protons
	TH1D *h1d_kerange_stop=new TH1D("h1d_kerange_stop","", nke, kemin, kemax); //range-based calc. (stopping protons)
	TH1D *h1d_kerange_stop_bmrw=new TH1D("h1d_kerange_stop_bmrw","", nke, kemin, kemax); //range-based calc. (stopping protons)
	TH1D *h1d_kecalo=new TH1D("h1d_kecalo","", nke, kemin, kemax); //calorimetric-based calc. (all reco. protons)
	TH1D *h1d_kecalo_stop=new TH1D("h1d_kecalo_stop","", nke, kemin, kemax); //calorimetric-based calc. (stopping protons)
	TH1D *h1d_kecalo_stop_bmrw=new TH1D("h1d_kecalo_stop_bmrw","", nke, kemin, kemax); //calorimetric-based calc. (stopping protons)

	TH1D *h1d_kecalo_recoinel=new TH1D("h1d_kecalo_recoinel","", nke, kemin, kemax); //calorimetric-based calc. (reco. inel. protons)
	TH1D *h1d_kecalo_recoinel_bmrw=new TH1D("h1d_kecalo_recoinel_bmrw","", nke, kemin, kemax); //calorimetric-based calc. (reco. inel. protons)

	TH1D *h1d_kecalo_recoel=new TH1D("h1d_kecalo_recoel","", nke, kemin, kemax); //calorimetric-based calc. (reco. inel. protons)

	//recoinel_with_truth_labels
	TH1D *h1d_kecalo_recoinel_bmrw_inel=new TH1D("h1d_kecalo_recoinel_bmrw_inel","", nke, kemin, kemax); //calorimetric-based calc. (reco. inel. protons)
	TH1D *h1d_kecalo_recoinel_bmrw_el=new TH1D("h1d_kecalo_recoinel_bmrw_el","", nke, kemin, kemax); //calorimetric-based calc. (reco. inel. protons)
	TH1D *h1d_kecalo_recoinel_bmrw_midcosmic=new TH1D("h1d_kecalo_recoinel_bmrw_midcosmic","", nke, kemin, kemax); //calorimetric-based calc. (reco. inel. protons)
	TH1D *h1d_kecalo_recoinel_bmrw_midpi=new TH1D("h1d_kecalo_recoinel_bmrw_midpi","", nke, kemin, kemax); //calorimetric-based calc. (reco. inel. protons)
	TH1D *h1d_kecalo_recoinel_bmrw_midp=new TH1D("h1d_kecalo_recoinel_bmrw_midp","", nke, kemin, kemax); //calorimetric-based calc. (reco. inel. protons)
	TH1D *h1d_kecalo_recoinel_bmrw_midmu=new TH1D("h1d_kecalo_recoinel_bmrw_midmu","", nke, kemin, kemax); //calorimetric-based calc. (reco. inel. protons)
	TH1D *h1d_kecalo_recoinel_bmrw_mideg=new TH1D("h1d_kecalo_recoinel_bmrw_mideg","", nke, kemin, kemax); //calorimetric-based calc. (reco. inel. protons)
	TH1D *h1d_kecalo_recoinel_bmrw_midother=new TH1D("h1d_kecalo_recoinel_bmrw_midother","", nke, kemin, kemax); //calorimetric-based calc. (reco. inel. protons)

	//keff of recoinel
	TH1D *h1d_keff_recoinel=new TH1D("h1d_keff_recoinel","", nke, kemin, kemax); //keff at z=0
	TH1D *h1d_keff_recoinel_inel=new TH1D("h1d_keff_recoinel_inel","", nke, kemin, kemax); 
	TH1D *h1d_keff_recoinel_el=new TH1D("h1d_keff_recoinel_el","", nke, kemin, kemax); 
	TH1D *h1d_keff_recoinel_midcosmic=new TH1D("h1d_keff_recoinel_midcosmic","", nke, kemin, kemax); 
	TH1D *h1d_keff_recoinel_midpi=new TH1D("h1d_keff_recoinel_midpi","", nke, kemin, kemax); 
	TH1D *h1d_keff_recoinel_midp=new TH1D("h1d_keff_recoinel_midp","", nke, kemin, kemax); 
	TH1D *h1d_keff_recoinel_midmu=new TH1D("h1d_keff_recoinel_midmu","", nke, kemin, kemax); 
	TH1D *h1d_keff_recoinel_mideg=new TH1D("h1d_keff_recoinel_mideg","", nke, kemin, kemax); 
	TH1D *h1d_keff_recoinel_midother=new TH1D("h1d_keff_recoinel_midother","", nke, kemin, kemax); 

	TH1D *h1d_keff0_recoinel=new TH1D("h1d_keff0_recoinel","", nke, kemin, kemax); //keff at tpc entrance
	TH1D *h1d_keff2_recoinel=new TH1D("h1d_keff2_recoinel","", nke, kemin, kemax); //keff=ke_beam-const. E-loss

	//keff of recoel
	TH1D *h1d_keff_recoel=new TH1D("h1d_keff_recoel","", nke, kemin, kemax); //keff at z=0
	TH1D *h1d_keff_recoel_inel=new TH1D("h1d_keff_recoel_inel","", nke, kemin, kemax); 
	TH1D *h1d_keff_recoel_el=new TH1D("h1d_keff_recoel_el","", nke, kemin, kemax); 
	TH1D *h1d_keff_recoel_midcosmic=new TH1D("h1d_keff_recoel_midcosmic","", nke, kemin, kemax); 
	TH1D *h1d_keff_recoel_midpi=new TH1D("h1d_keff_recoel_midpi","", nke, kemin, kemax); 
	TH1D *h1d_keff_recoel_midp=new TH1D("h1d_keff_recoel_midp","", nke, kemin, kemax); 
	TH1D *h1d_keff_recoel_midmu=new TH1D("h1d_keff_recoel_midmu","", nke, kemin, kemax); 
	TH1D *h1d_keff_recoel_mideg=new TH1D("h1d_keff_recoel_mideg","", nke, kemin, kemax); 
	TH1D *h1d_keff_recoel_midother=new TH1D("h1d_keff_recoel_midother","", nke, kemin, kemax); 

	TH1D *h1d_keff0_recoel=new TH1D("h1d_keff0_recoel","", nke, kemin, kemax); //keff at tpc entrance
	TH1D *h1d_keff2_recoel=new TH1D("h1d_keff2_recoel","", nke, kemin, kemax); //keff=ke_beam-const. E-loss

	//keff bias
	int ndkeff=320;
	double dkeff_min=-800;
	double dkeff_max=800;

	TH1D *h1d_dkeff_recoinel=new TH1D("h1d_dkeff_recoinel","", ndkeff, dkeff_min, dkeff_max); 
	TH1D *h1d_dkeff_recoel=new TH1D("h1d_dkeff_recoel","", ndkeff, dkeff_min, dkeff_max); 
	TH1D *h1d_dkeff_stop=new TH1D("h1d_dkeff_stop","", ndkeff, dkeff_min, dkeff_max); 
	TH1D *h1d_dkeff_all=new TH1D("h1d_dkeff_all","", ndkeff, dkeff_min, dkeff_max); 
	TH1D *h1d_dkeff2_all=new TH1D("h1d_dkeff2_all","", ndkeff, dkeff_min, dkeff_max); 

	//ke vs trklen
	int n2d_trklen=140;
	double n2d_trklen_min=0;
	double n2d_trklen_max=140;
	TH2D *h2d_trklen_ke_recoinel=new TH2D("h2d_trklen_ke_recoinel","", n2d_trklen, n2d_trklen_min, n2d_trklen_max, nke, kemin, kemax); 
	TH2D *h2d_trklen_ke2_recoinel=new TH2D("h2d_trklen_ke2_recoinel","", n2d_trklen, n2d_trklen_min, n2d_trklen_max, nke, kemin, kemax);
	TH2D *h2d_trklen_dke_recoinel=new TH2D("h2d_trklen_dke_recoinel","", n2d_trklen, n2d_trklen_min, n2d_trklen_max, ndkeff, dkeff_min, dkeff_max);
	
	TH2D *h2d_trklen_ke_recoel=new TH2D("h2d_trklen_ke_recoel","", n2d_trklen, n2d_trklen_min, n2d_trklen_max, nke, kemin, kemax); 
	TH2D *h2d_trklen_ke2_recoel=new TH2D("h2d_trklen_ke2_recoel","", n2d_trklen, n2d_trklen_min, n2d_trklen_max, nke, kemin, kemax);
	TH2D *h2d_trklen_dke_recoel=new TH2D("h2d_trklen_dke_recoel","", n2d_trklen, n2d_trklen_min, n2d_trklen_max, ndkeff, dkeff_min, dkeff_max);

	TH2D *h2d_trklen_ke_stop=new TH2D("h2d_trklen_ke_stop","", n2d_trklen, n2d_trklen_min, n2d_trklen_max, nke, kemin, kemax); 
	TH2D *h2d_trklen_ke2_stop=new TH2D("h2d_trklen_ke2_stop","", n2d_trklen, n2d_trklen_min, n2d_trklen_max, nke, kemin, kemax);
	TH2D *h2d_trklen_dke_stop=new TH2D("h2d_trklen_dke_stop","", n2d_trklen, n2d_trklen_min, n2d_trklen_max, ndkeff, dkeff_min, dkeff_max);
	
	TH2D *h2d_trklen_ke_all=new TH2D("h2d_trklen_ke_all","", n2d_trklen, n2d_trklen_min, n2d_trklen_max, nke, kemin, kemax); 
	TH2D *h2d_trklen_ke2_all=new TH2D("h2d_trklen_ke2_all","", n2d_trklen, n2d_trklen_min, n2d_trklen_max, nke, kemin, kemax);
	TH2D *h2d_trklen_dke_all=new TH2D("h2d_trklen_dke_all","", n2d_trklen, n2d_trklen_min, n2d_trklen_max, ndkeff, dkeff_min, dkeff_max);







	//upstream energy loss
	int ndke=300;
	double dkemin=0;
	double dkemax=600;

	TH1D *h1d_dEbb=new TH1D("h1d_dEbb","", 1000, -500,500); //edept_bb_true-edept_bb_reco
	TH1D *h1d_dEbb_stop=new TH1D("h1d_dEbb_stop","", 1000, -500,500); //edept_bb_true-edept_bb_reco
	TH1D *h1d_dEbb_recoinel=new TH1D("h1d_dEbb_recoinel","", 1000, -500,500); //edept_bb_true-edept_bb_reco
	TH1D *h1d_dEbb_recoel=new TH1D("h1d_dEbb_recoel","", 1000, -500,500); //edept_bb_true-edept_bb_reco

	TH1D *h1d_dKEbb=new TH1D("h1d_dKEbb","", 1000, -500,500); //KEbb_true-KEbb_reco
	TH2D *h2d_dKEbb_KEtrue=new TH2D("h2d_dKEbb_KEtrue","", 800,0,800, 1000, -500, 500); 
	TH2D *h2d_dKEbb_KEreco=new TH2D("h2d_dKEbb_KEreco","", 800,0,800, 1000, -500, 500); 


	TH2D *h2d_dE1_KEtrue=new TH2D("h2d_dE1_KEtrue","", 800,0,800, 1000, -500, 500); 
	TH2D *h2d_dE1_KEreco=new TH2D("h2d_dE1_KEreco","", 800,0,800, 1000, -500, 500); 
	TH2D *h2d_dE2_KEtrue=new TH2D("h2d_dE2_KEtrue","", 800,0,800, 1000, -500, 500); 
	TH2D *h2d_dE2_KEreco=new TH2D("h2d_dE2_KEreco","", 800,0,800, 1000, -500, 500); 



	//reco_label
	TH1D *h1d_dke=new TH1D("h1d_dke","", ndke, dkemin, dkemax); //ke_beam-ke_ff
	TH1D *h1d_dke0=new TH1D("h1d_dke0","", ndke, dkemin, dkemax); //ke0-ke_ff
	TH1D *h1d_dke_stop=new TH1D("h1d_dke_stop","", ndke, dkemin, dkemax); //ke_beam-ke_ff
	TH1D *h1d_dke_recoinel=new TH1D("h1d_dke_recoinel","", ndke, dkemin, dkemax); //ke_beam-ke_ff
	TH1D *h1d_dke_recoel=new TH1D("h1d_dke_recoel","", ndke, dkemin, dkemax); //ke_beam-ke_ff

	//reco_label
	TH1D *h1d_dKE=new TH1D("h1d_dKE","", ndke, dkemin, dkemax); //ke_beam-E_dept
	TH1D *h1d_dKE_stop=new TH1D("h1d_dKE_stop","", ndke, dkemin, dkemax); //ke_beam-E_dept
	TH1D *h1d_dKE_recoinel=new TH1D("h1d_dKE_recoinel","", ndke, dkemin, dkemax); //ke_beam-E_dept
	TH1D *h1d_dKE_recoel=new TH1D("h1d_dKE_recoel","", ndke, dkemin, dkemax); //ke_beam-E_dept

	//true_label
	TH1D *h1d_dke_el=new TH1D("h1d_dke_el","", ndke, dkemin, dkemax); //ke_beam-ke_ff
	TH1D *h1d_dke_inel=new TH1D("h1d_dke_inel","", ndke, dkemin, dkemax); //ke_beam-ke_ff
	TH1D *h1d_dke_misidp=new TH1D("h1d_dke_misidp","", ndke, dkemin, dkemax); //ke_beam-ke_ff

	//KEbb [all protons]
	TH1D *h1d_KEbb_reco=new TH1D("h1d_KEbb_reco","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_true=new TH1D("h1d_KEbb_true","", ndke, dkemin, dkemax);
	
	TH1D *h1d_KEbb_reco_inel=new TH1D("h1d_KEbb_reco_inel","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_reco_el=new TH1D("h1d_KEbb_reco_el","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_reco_midcosmic=new TH1D("h1d_KEbb_reco_midcosmic","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_reco_midpi=new TH1D("h1d_KEbb_reco_midpi","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_reco_midp=new TH1D("h1d_KEbb_reco_midp","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_reco_midmu=new TH1D("h1d_KEbb_reco_midmu","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_reco_mideg=new TH1D("h1d_KEbb_reco_mideg","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_reco_midother=new TH1D("h1d_KEbb_reco_midother","", ndke, dkemin, dkemax);

	//KEbb [reco inel]
	TH1D *h1d_KEbb_reco_RecoInel=new TH1D("h1d_KEbb_reco_RecoInel","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_true_RecoInel=new TH1D("h1d_KEbb_true_RecoInel","", ndke, dkemin, dkemax);
	
	TH1D *h1d_KEbb_reco_inel_RecoInel=new TH1D("h1d_KEbb_reco_inel_RecoInel","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_reco_el_RecoInel=new TH1D("h1d_KEbb_reco_el_RecoInel","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_reco_midcosmic_RecoInel=new TH1D("h1d_KEbb_reco_midcosmic_RecoInel","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_reco_midpi_RecoInel=new TH1D("h1d_KEbb_reco_midpi_RecoInel","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_reco_midp_RecoInel=new TH1D("h1d_KEbb_reco_midp_RecoInel","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_reco_midmu_RecoInel=new TH1D("h1d_KEbb_reco_midmu_RecoInel","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_reco_mideg_RecoInel=new TH1D("h1d_KEbb_reco_mideg_RecoInel","", ndke, dkemin, dkemax);
	TH1D *h1d_KEbb_reco_midother_RecoInel=new TH1D("h1d_KEbb_reco_midother_RecoInel","", ndke, dkemin, dkemax);


	TH2D *h2d_dKEbb_KEtrue_RecoInel=new TH2D("h2d_dKEbb_KEtrue_RecoInel","", 800,0,800, 1000, -500, 500); 
	TH2D *h2d_dKEbb_KEreco_RecoInel=new TH2D("h2d_dKEbb_KEreco_RecoInel","", 800,0,800, 1000, -500, 500); 

	//------------------------------------------------------------------------------------------------------------------------//

	//Beam momentum reweighting ----------------------------------------------------------------------------------------------//
	//MC Beam Mom Gaussian 
	double m1=1007.1482; //MC prod4a [spec]
	double s1=60.703307; //MC prod4a [spec]

	//momentum cut range	
	double mu_min=m1-3.*s1;
	double mu_max=m1+3.*s1;

	//default gaussian
	int nx=250;	
	double xmin=0.; //pmin [MeV/c]
	double xmax=2000.; //pmax [MeV/c]
	TF1 *g1=new TF1("g1",fitg,xmin,xmax,2);
	g1->SetName("g1");
	g1->SetParameter(0,m1);
	g1->SetParameter(1,s1);

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
	int index_minchi2=13331; //index of minchi2
	for (int imu=0; imu<nmu; ++imu){ //mu loop
		double frac_mu=mu_st-(double)imu*dmu;
		double mu=m1*frac_mu;
		for (int isigma=0; isigma<nsigma; ++isigma){ //sigma loop
			double frac_sigma=sigma_st-(double)isigma*dsigma;
			double sigma=s1*frac_sigma;

			//if (mu==m1&&sigma==s1) { //no rw
			if (std::abs(mu-m1)<0.0001&&std::abs(sigma-s1)<0.0001) { //no rw
				index_original=cnt_array;
				mu=m1;
				sigma=s1;
			} //no rw

			//Gaussian with changed mean and sigma
			gn[cnt_array]=new TF1(Form("gn_%d",cnt_array),fitg,xmin,xmax,2);
			gn[cnt_array]->SetParameter(0,mu);
			gn[cnt_array]->SetParameter(1,sigma);

			//weighting func. (beam mom)
			gng[cnt_array]=new TF1(Form("gng_%d",cnt_array),govg,xmin,xmax,4);
			gng[cnt_array]->SetParameter(0,m1);
			gng[cnt_array]->SetParameter(1,s1);
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
			}

			kereco_calo+=cali_dedx*pitch;
			kereco_range+=pitch*dedx_predict(resrange_reco);
			kereco_range2+=pitch*(double)gr_predict_dedx_resrange->Eval(resrange_reco);

			//if (kinel) rangereco_dedxreco_TrueInEL->Fill(range_reco, cali_dedx);
			//if (kel) { 
						//rangereco_dedxreco_TrueEL->Fill(range_reco, cali_dedx);
						//rr_dedx_truestop->Fill(resrange_reco, cali_dedx);
			//}
		  } //loop over reco hits of a given track
		  //range_reco=primtrk_range->at(0);
		} //if calo size not empty

		//Reco stopping/Inel p cut
		bool IsRecoStop=false;
		bool IsRecoInEL=false;
		bool IsRecoEL=false;
		double mom_beam_spec=-99; mom_beam_spec=beamMomentum_spec->at(0);
		//double range_reco=-99; if (!primtrk_range->empty()) range_reco=primtrk_range->at(0); //reco primary trklen
		double csda_val_spec=csda_range_vs_mom_sm->Eval(mom_beam_spec);

		if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) IsRecoStop=true;
		//if ((range_reco/csda_val_spec)<min_norm_trklen_csda) IsRecoInEL=true;

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
				//true_trklen_accum[iz] = range_true;
		  }
		} //is beam at ff						    	
		//fix on the truth length by adding distance between 1st tpc hit to front face ------------------------------------------------------//

		//KEs ---------------------------------------------------------------------------------------------------//
		//energy loss using stopping protons
		double mean_Eloss_upstream=19.3073;
		double err_mean_Eloss_upstream=0.187143;
		double sigma_Eloss_upstream=18.7378;
		double err_sigma_Eloss_upstream=0.140183;


		double Eloss_upstream=0; 
		double Eloss_upstream_reco=0;
		double dKE_recotruth_1=0;
		double dKE_recotruth_2=0;
		if (is_beam_at_ff) { //if beam reach tpc
			Eloss_upstream=ke_beam_spec_MeV-KE_ff;
			Eloss_upstream_reco=ke_beam_spec_MeV-mean_Eloss_upstream;
		} //if beam reach tpc

		//double ke_trklen=ke_vs_csda_range_sm->Eval(range_reco); //[unit: GeV]
		//double ke_trklen_MeV=1000.*ke_trklen; //[unit: MeV]

		double dEbb_true=0; if (range_true>=0&&is_beam_at_ff) dEbb_true=BB.KEFromRangeSpline(range_true);
		double dEbb_reco=0; if (range_reco>=0&&is_beam_at_ff) dEbb_reco=BB.KEFromRangeSpline(range_reco);
		double dEbb=0; dEbb=dEbb_true-dEbb_reco;

		//mean and sigma of e-loss [mc]
		//double mean_Eloss_upstream=2.04949e+01;
		//double err_mean_Eloss_upstream=1.47946e-01;
		//double sigma_Eloss_upstream=1.90634e+01;
		//double err_sigma_Eloss_upstream=1.23169e-01;

		
		double keff_reco=ke_beam_spec_MeV-mean_Eloss_upstream;
		//double keff_reco=(1000.*beamtrk_Eng->at(0))-mean_Eloss_upstream;

		double KEbb_true=-1; KEbb_true=BB.KEAtLength(KE_ff, range_true);
		//double KEbb_reco=-1; KEbb_reco=BB.KEAtLength(KE_ff, range_reco);
		double KEbb_reco=-1; KEbb_reco=BB.KEAtLength((keff_reco), range_reco);
		double dKEbb=0; dKEbb=KEbb_reco-KEbb_true;


		double ke_calo_MeV=0;
		if (IsPandoraSlice&&IsBQ&&IsCaloSize) { //if calo size not empty
			//no-sce
			//h2d_xy_noSCE->Fill(reco_stx_noSCE, reco_sty_noSCE);
			//h1d_zst_noSCE->Fill(reco_stz_noSCE);
			//after sce
			//h2d_xy_SCE->Fill(reco_stx, reco_sty);
			//h1d_zst_SCE->Fill(reco_stz);
			
			//calo
			vector<double> trkdedx; 
			vector<double> trkres;
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

				//if (IsRecoStop) {
					//h2d_rr_dedx_recoSTOP->Fill(resrange_reco, cali_dedx);
				//}

			} //loop over reco hits of a given track

			//if (IsRecoStop) chi2pid_recostop->Fill(chi2pid(trkdedx,trkres));
			//if (IsRecoInEL) chi2pid_recoinel->Fill(chi2pid(trkdedx,trkres));
			//if (IsPureMCS) chi2pid_truestop->Fill(chi2pid(trkdedx,trkres));
			//if (kel) chi2pid_trueel->Fill(chi2pid(trkdedx,trkres));
			//if (kinel) chi2pid_trueinel->Fill(chi2pid(trkdedx,trkres)); 

		} //if calo size not empty

		if (IsPandoraSlice&&IsBQ&&IsCaloSize) { //basic cuts
			h1d_kebeam->Fill(ke_beam_spec_MeV);
			h1d_ke0->Fill(ke0);
			Fill1DHist(h1d_keff, ke_ff);
			Fill1DHist(h1d_keff_reco, keff_reco);

			h1d_keff2->Fill(ke_beam_spec_MeV-Eloss_mc);
			h1d_keff0->Fill(KE_ff_tpc);
			//Fill1DHist(h1d_dkeff_all, ke_beam_spec_MeV-Eloss_mc-ke_ff);
			Fill1DHist(h1d_dkeff_all, keff_reco-ke_ff);
			h2d_trklen_ke_all->Fill(range_reco, (ke_ff-ke_calo_MeV));
			h2d_trklen_ke2_all->Fill(range_reco, (ke_beam_spec_MeV-Eloss_mc-ke_calo_MeV));
			h2d_trklen_dke_all->Fill(range_reco, (ke_beam_spec_MeV-Eloss_mc-ke_ff));

			if (kinel) { 
				h1d_keff_inel->Fill(ke_ff);
			}
			if (kel) { 
				h1d_keff_el->Fill(ke_ff);
			}
			if (kMIDcosmic) { 
				h1d_keff_midcosmic->Fill(ke_ff);
			}
			if (kMIDpi) { 
				h1d_keff_midpi->Fill(ke_ff);
			}
			if (kMIDp) { 
				h1d_keff_midp->Fill(ke_ff);
			}
			if (kMIDmu) { 
				h1d_keff_midmu->Fill(ke_ff);
			}
			if (kMIDeg) { 
				h1d_keff_mideg->Fill(ke_ff);
			}
			if (kMIDother) { 
				h1d_keff_midother->Fill(ke_ff);					
			}


			//h1d_trklen->Fill(range_reco);	
			//h1d_zend->Fill(reco_endz);

			double mom_rw_minchi2=1.;
			if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) mom_rw_minchi2=gng[index_minchi2]->Eval(mom_beam_spec*1000.); //bmrw

			h1d_kebeam_bmrw->Fill(ke_beam_spec_MeV, mom_rw_minchi2);
			//h1d_trklen_bmrw->Fill(range_reco, mom_rw_minchi2); 
			//h1d_zend_bmrw->Fill(reco_endz, mom_rw_minchi2);	

			//if (IsXY) { //xy
				//h1d_trklen_XY->Fill(range_reco);
				//h1d_trklen_bmrw_XY->Fill(range_reco, mom_rw_minchi2);

				//h1d_zend_XY->Fill(reco_endz);
				//h1d_zend_bmrw_XY->Fill(reco_endz, mom_rw_minchi2);
			//} //xy

			if (is_beam_at_ff) {
				h1d_dke->Fill(Eloss_upstream);
				h1d_dke0->Fill(ke0-KE_ff);
			}


			Fill1DHist(h1d_dEbb, dEbb);
			Fill1DHist(h1d_dKEbb, dKEbb);
			h2d_dKEbb_KEtrue->Fill(KEbb_true, dKEbb);	
			h2d_dKEbb_KEreco->Fill(KEbb_reco, dKEbb);

			h2d_dE1_KEtrue->Fill(KEbb_true, dEbb);
			h2d_dE1_KEreco->Fill(KEbb_reco, dEbb);

			h2d_dE2_KEtrue->Fill(KEbb_true, keff_reco-KE_ff);
			h2d_dE2_KEreco->Fill(KEbb_reco, keff_reco-KE_ff);


			Fill1DHist(h1d_KEbb_true, KEbb_true);
			Fill1DHist(h1d_KEbb_reco, KEbb_reco);
		
			h1d_dKE->Fill(ke_beam_spec_MeV-ke_calo_MeV);
			h1d_kecalo->Fill(ke_calo_MeV);

			if (kinel) { //inel
				h1d_dke_inel->Fill(Eloss_upstream);
				Fill1DHist(h1d_KEbb_reco_inel, KEbb_reco);
			} //inel

			if (kel) { //el
				h1d_dke_el->Fill(Eloss_upstream);
				Fill1DHist(h1d_KEbb_reco_el, KEbb_reco);
			} //el

			if (kMIDp) { //misid:p
				h1d_dke_misidp->Fill(Eloss_upstream);
				Fill1DHist(h1d_KEbb_reco_midp, KEbb_reco);
			} //misid:p

			if (kMIDcosmic) { 
				Fill1DHist(h1d_KEbb_reco_midcosmic, KEbb_reco);
			}

			if (kMIDpi) { 
				Fill1DHist(h1d_KEbb_reco_midpi, KEbb_reco);
			}

			if (kMIDmu) { 
				Fill1DHist(h1d_KEbb_reco_midmu, KEbb_reco);
			}

			if (kMIDeg) { 
				Fill1DHist(h1d_KEbb_reco_mideg, KEbb_reco);
			}

			if (kMIDother) { 
				Fill1DHist(h1d_KEbb_reco_midother, KEbb_reco);
			}




			if (IsRecoInEL) {
				Fill1DHist(h1d_dEbb_recoinel, dEbb);


				Fill1DHist(h1d_KEbb_true_RecoInel, KEbb_true);
				Fill1DHist(h1d_KEbb_reco_RecoInel, KEbb_reco);


				h2d_dKEbb_KEtrue_RecoInel->Fill(KEbb_true, dKEbb);	
				h2d_dKEbb_KEreco_RecoInel->Fill(KEbb_reco, dKEbb);



				h1d_dke_recoinel->Fill(Eloss_upstream);
				h1d_kecalo_recoinel->Fill(ke_calo_MeV);
				h1d_keff_recoinel->Fill(ke_ff);
				h1d_keff0_recoinel->Fill(KE_ff_tpc);
				h1d_keff2_recoinel->Fill(ke_beam_spec_MeV-Eloss_mc);
				h1d_kecalo_recoinel_bmrw->Fill(ke_calo_MeV,mom_rw_minchi2);
				h1d_dKE_recoinel->Fill(ke_beam_spec_MeV-ke_calo_MeV);

				Fill1DHist(h1d_dkeff_recoinel, ke_beam_spec_MeV-Eloss_mc-ke_ff);
				h2d_trklen_ke_recoinel->Fill(range_reco, (ke_ff-ke_calo_MeV));
				h2d_trklen_ke2_recoinel->Fill(range_reco, (ke_beam_spec_MeV-Eloss_mc-ke_calo_MeV));
				h2d_trklen_dke_recoinel->Fill(range_reco, (ke_beam_spec_MeV-Eloss_mc-ke_ff));

				if (kinel) { 
					h1d_kecalo_recoinel_bmrw_inel->Fill(ke_calo_MeV, mom_rw_minchi2);
					h1d_keff_recoinel_inel->Fill(ke_ff);

					Fill1DHist(h1d_KEbb_reco_inel_RecoInel, KEbb_reco);

				}
				if (kel) { 
					h1d_kecalo_recoinel_bmrw_el->Fill(ke_calo_MeV, mom_rw_minchi2);
					h1d_keff_recoinel_el->Fill(ke_ff);
					Fill1DHist(h1d_KEbb_reco_el_RecoInel, KEbb_reco);

				}
				if (kMIDcosmic) { 
					h1d_kecalo_recoinel_bmrw_midcosmic->Fill(ke_calo_MeV, mom_rw_minchi2);
					h1d_keff_recoinel_midcosmic->Fill(ke_ff);
					Fill1DHist(h1d_KEbb_reco_midcosmic_RecoInel, KEbb_reco);
				}
				if (kMIDpi) { 
					h1d_kecalo_recoinel_bmrw_midpi->Fill(ke_calo_MeV, mom_rw_minchi2);
					h1d_keff_recoinel_midpi->Fill(ke_ff);
					Fill1DHist(h1d_KEbb_reco_midpi_RecoInel, KEbb_reco);
				}
				if (kMIDp) { 
					h1d_kecalo_recoinel_bmrw_midp->Fill(ke_calo_MeV, mom_rw_minchi2);
					h1d_keff_recoinel_midp->Fill(ke_ff);
					Fill1DHist(h1d_KEbb_reco_midp_RecoInel, KEbb_reco);
				}
				if (kMIDmu) { 
					h1d_kecalo_recoinel_bmrw_midmu->Fill(ke_calo_MeV, mom_rw_minchi2);
					h1d_keff_recoinel_midmu->Fill(ke_ff);
					Fill1DHist(h1d_KEbb_reco_midmu_RecoInel, KEbb_reco);
				}
				if (kMIDeg) { 
					h1d_kecalo_recoinel_bmrw_mideg->Fill(ke_calo_MeV, mom_rw_minchi2);
					h1d_keff_recoinel_mideg->Fill(ke_ff);
					Fill1DHist(h1d_KEbb_reco_mideg_RecoInel, KEbb_reco);
				}
				if (kMIDother) { 
					h1d_kecalo_recoinel_bmrw_midother->Fill(ke_calo_MeV, mom_rw_minchi2);
					h1d_keff_recoinel_midother->Fill(ke_ff);
					Fill1DHist(h1d_KEbb_reco_midother_RecoInel, KEbb_reco);
				}
			}
			if (IsRecoEL) {
				Fill1DHist(h1d_dEbb_recoel, dEbb);

				h1d_dke_recoel->Fill(Eloss_upstream);
				h1d_kecalo_recoel->Fill(ke_calo_MeV);
				h1d_dKE_recoel->Fill(ke_beam_spec_MeV-ke_calo_MeV);

				h1d_keff_recoel->Fill(ke_ff);
				h1d_keff0_recoel->Fill(KE_ff_tpc);
				h1d_keff2_recoel->Fill(ke_beam_spec_MeV-Eloss_mc);

				Fill1DHist(h1d_dkeff_recoel, ke_beam_spec_MeV-Eloss_mc-ke_ff);
				h2d_trklen_ke_recoel->Fill(range_reco, (ke_ff-ke_calo_MeV));
				h2d_trklen_ke2_recoel->Fill(range_reco, (ke_beam_spec_MeV-Eloss_mc-ke_calo_MeV));
				h2d_trklen_dke_recoel->Fill(range_reco, (ke_beam_spec_MeV-Eloss_mc-ke_ff));

				if (kinel) { 
					h1d_keff_recoel_inel->Fill(ke_ff);
				}
				if (kel) { 
					h1d_keff_recoel_el->Fill(ke_ff);
				}
				if (kMIDcosmic) { 
					h1d_keff_recoel_midcosmic->Fill(ke_ff);
				}
				if (kMIDpi) { 
					h1d_keff_recoel_midpi->Fill(ke_ff);
				}
				if (kMIDp) { 
					h1d_keff_recoel_midp->Fill(ke_ff);
				}
				if (kMIDmu) { 
					h1d_keff_recoel_midmu->Fill(ke_ff);
				}
				if (kMIDeg) { 
					h1d_keff_recoel_mideg->Fill(ke_ff);
				}
				if (kMIDother) { 
					h1d_keff_recoel_midother->Fill(ke_ff);					
				}

			}

			if (IsRecoStop) { //reco stop
				Fill1DHist(h1d_dEbb_stop, dEbb);
				h1d_kebeam_stop->Fill(ke_beam_spec_MeV);
				h1d_keff_stop->Fill(ke_ff);
				h1d_keff0_stop->Fill(KE_ff_tpc);
				h1d_keff2_stop->Fill(ke_beam_spec_MeV-Eloss_mc);

				Fill1DHist(h1d_dkeff_stop, ke_beam_spec_MeV-Eloss_mc-ke_ff);
				h2d_trklen_ke_stop->Fill(range_reco, (ke_ff-ke_calo_MeV));
				h2d_trklen_ke2_stop->Fill(range_reco, (ke_beam_spec_MeV-Eloss_mc-ke_calo_MeV));
				h2d_trklen_dke_stop->Fill(range_reco, (ke_beam_spec_MeV-Eloss_mc-ke_ff));

				h1d_kecalo_stop->Fill(ke_calo_MeV);
				h1d_kecalo_stop_bmrw->Fill(ke_calo_MeV,mom_rw_minchi2);
				h1d_kerange_stop->Fill(ke_trklen_MeV);
				
				h1d_dke_stop->Fill(Eloss_upstream);
				h1d_dKE_stop->Fill(ke_beam_spec_MeV-ke_calo_MeV);

				double mom_rw_minchi2=1.;
				if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) mom_rw_minchi2=gng[index_minchi2]->Eval(mom_beam_spec*1000.); //bmrw

				h1d_kebeam_stop_bmrw->Fill(ke_beam_spec_MeV, mom_rw_minchi2);
				h1d_kerange_stop_bmrw->Fill(ke_trklen_MeV, mom_rw_minchi2);
			} //reco stop
		} //basic cuts

	} //main entry loop

	//Fit Gaussians on momenta ...
	TF1* kebeam_fit; kebeam_fit=VFit(h1d_kebeam, 2);
	kebeam_fit->SetName("kebeam_fit");

	TF1* kebeam_bmrw_fit; kebeam_bmrw_fit=VFit(h1d_kebeam_bmrw, 2);
	kebeam_bmrw_fit->SetName("kebeam_bmrw_fit");

	TF1* kebeam_stop_fit; kebeam_stop_fit=VFit(h1d_kebeam_stop, 2);
	kebeam_stop_fit->SetName("kebeam_stop_fit");

	TF1* kecalo_stop_fit; kecalo_stop_fit=VFit(h1d_kecalo_stop, 2);
	kecalo_stop_fit->SetName("kecalo_stop_fit");

	TF1* kerange_stop_fit; kerange_stop_fit=VFit(h1d_kerange_stop, 2);
	kerange_stop_fit->SetName("kerange_stop_fit");

	TF1* keff_stop_fit; keff_stop_fit=VFit(h1d_keff_stop, 2);
	keff_stop_fit->SetName("keff_stop_fit");

	TF1* kerange_stop_bmrw_fit; kerange_stop_bmrw_fit=VFit(h1d_kerange_stop_bmrw, 2);
	kerange_stop_bmrw_fit->SetName("kerange_stop_bmrw_fit");

	TF1* kebeam_stop_bmrw_fit; kebeam_stop_bmrw_fit=VFit(h1d_kebeam_stop_bmrw, 2);
	kebeam_stop_bmrw_fit->SetName("kebeam_stop_bmrw_fit");



	//save results...
   	//TFile *fout = new TFile("mc_ke.root","RECREATE");
   	TFile *fout = new TFile("mc_ke0truth.root","RECREATE");
   	//TFile *fout = new TFile("mc_ke_crosscheck_KEfftrue.root","RECREATE");
		h1d_kebeam->Write();
		h1d_kebeam_bmrw->Write();
		h1d_ke0->Write();

		h1d_kebeam_stop->Write();
		h1d_kebeam_stop_bmrw->Write();

		h1d_keff->Write();
		h1d_keff2->Write();
		h1d_keff0->Write();

		h1d_keff_reco->Write();

		h1d_keff_inel->Write();
		h1d_keff_el->Write();
		h1d_keff_midcosmic->Write();
		h1d_keff_midpi->Write();
		h1d_keff_midp->Write();
		h1d_keff_midmu->Write();
		h1d_keff_mideg->Write();
		h1d_keff_midother->Write();

		h1d_keff_stop->Write();
		h1d_keff0_stop->Write();
		h1d_keff2_stop->Write();

		kebeam_fit->Write();
		kebeam_bmrw_fit->Write();
		kebeam_stop_fit->Write();
		kebeam_stop_bmrw_fit->Write();

		kecalo_stop_fit->Write();
		kerange_stop_fit->Write();
		kerange_stop_bmrw_fit->Write();
		keff_stop_fit->Write();

		h1d_kecalo->Write();
		h1d_kecalo_stop->Write();
		h1d_kecalo_stop_bmrw->Write();

		h1d_kecalo_recoinel->Write();

		h1d_keff_recoinel->Write();
		h1d_keff_recoinel_inel->Write();
		h1d_keff_recoinel_el->Write();
		h1d_keff_recoinel_midcosmic->Write();
		h1d_keff_recoinel_midpi->Write();
		h1d_keff_recoinel_midp->Write();
		h1d_keff_recoinel_midmu->Write();
		h1d_keff_recoinel_mideg->Write();
		h1d_keff_recoinel_midother->Write();

		h1d_keff0_recoinel->Write();
		h1d_keff2_recoinel->Write();



		h1d_keff_recoel->Write();
		h1d_keff_recoel_inel->Write();
		h1d_keff_recoel_el->Write();
		h1d_keff_recoel_midcosmic->Write();
		h1d_keff_recoel_midpi->Write();
		h1d_keff_recoel_midp->Write();
		h1d_keff_recoel_midmu->Write();
		h1d_keff_recoel_mideg->Write();
		h1d_keff_recoel_midother->Write();

		h1d_keff0_recoel->Write();
		h1d_keff2_recoel->Write();


		//Edept (recoinel) after reweighting --------
		h1d_kecalo_recoinel_bmrw->Write();
		h1d_kecalo_recoinel_bmrw_inel->Write();
		h1d_kecalo_recoinel_bmrw_el->Write();
		h1d_kecalo_recoinel_bmrw_midcosmic->Write();
		h1d_kecalo_recoinel_bmrw_midpi->Write();
		h1d_kecalo_recoinel_bmrw_midp->Write();
		h1d_kecalo_recoinel_bmrw_midmu->Write();
		h1d_kecalo_recoinel_bmrw_mideg->Write();
		h1d_kecalo_recoinel_bmrw_midother->Write();

		h1d_kecalo_recoel->Write();
		h1d_kerange_stop->Write();
		h1d_kerange_stop_bmrw->Write();

		h1d_dke0->Write();
		h1d_dke->Write();
		h1d_dke_stop->Write();
		h1d_dke_recoinel->Write();
		h1d_dke_recoel->Write();

		h1d_dke_el->Write();
		h1d_dke_inel->Write();
		h1d_dke_misidp->Write();

		h1d_dKE->Write();
		h1d_dKE_stop->Write();
		h1d_dKE_recoinel->Write();
		h1d_dKE_recoel->Write();



		h1d_dkeff_recoinel->Write();
		h1d_dkeff_recoel->Write();
		h1d_dkeff_stop->Write();
		h1d_dkeff_all->Write();

		h2d_trklen_ke_recoinel->Write();
		h2d_trklen_ke2_recoinel->Write();
		h2d_trklen_dke_recoinel->Write();

		h2d_trklen_ke_recoel->Write();
		h2d_trklen_ke2_recoel->Write();
		h2d_trklen_dke_recoel->Write();

		h2d_trklen_ke_stop->Write();
		h2d_trklen_ke2_stop->Write();
		h2d_trklen_dke_stop->Write();

		h2d_trklen_ke_all->Write();
		h2d_trklen_ke2_all->Write();
		h2d_trklen_dke_all->Write();


		h1d_dEbb->Write();
		h1d_dEbb_stop->Write();
		h1d_dEbb_recoinel->Write();
		h1d_dEbb_recoel->Write();
		h1d_dKEbb->Write();
		h2d_dKEbb_KEtrue->Write();
		h2d_dKEbb_KEreco->Write();


		h1d_KEbb_true->Write();
		h1d_KEbb_reco->Write();

		h1d_KEbb_reco_inel->Write();
		h1d_KEbb_reco_el->Write();
		h1d_KEbb_reco_midcosmic->Write();
		h1d_KEbb_reco_midpi->Write();
		h1d_KEbb_reco_midp->Write();
		h1d_KEbb_reco_midmu->Write();
		h1d_KEbb_reco_mideg->Write();
		h1d_KEbb_reco_midother->Write();


		h2d_dKEbb_KEtrue_RecoInel->Write();
		h2d_dKEbb_KEreco_RecoInel->Write();
		h1d_KEbb_reco_RecoInel->Write();
		h1d_KEbb_true_RecoInel->Write();

		h1d_KEbb_reco_inel_RecoInel->Write();
		h1d_KEbb_reco_el_RecoInel->Write();
		h1d_KEbb_reco_midcosmic_RecoInel->Write();
		h1d_KEbb_reco_midpi_RecoInel->Write();
		h1d_KEbb_reco_midp_RecoInel->Write();
		h1d_KEbb_reco_midmu_RecoInel->Write();
		h1d_KEbb_reco_mideg_RecoInel->Write();
		h1d_KEbb_reco_midother_RecoInel->Write();

		h2d_dE1_KEtrue->Write();
		h2d_dE1_KEreco->Write();
		h2d_dE2_KEtrue->Write();
		h2d_dE2_KEreco->Write();

	fout->Close();



}
