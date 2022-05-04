#define ProtonDataDrivenBKGMeas_BetheBloch_cxx
#include "ProtonDataDrivenBKGMeas_BetheBloch.h"

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
#include "./headers/ESliceParams.h"
#include "./headers/util.h"
//#include "./headers/ESlice.h"
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

//[MC] Energy loss using stopping protons
double mean_Eloss_upstream=19.3073;
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

void ProtonDataDrivenBKGMeas_BetheBloch::Loop() {
	if (fChain == 0) return;

	//MC beam momentum -----------------------//
	double mm1=1007.1482; //MC prod4a [spec]
	double ss1=60.703307; //MC prod4a [spec]
	double mmu_min=mm1-3.*ss1;
	double mmu_max=mm1+3.*ss1;

	//Booking histograms -------------------------------------------------------------------------------//
	int nke=160;
	//int nke=80;
	double kemin=-800;
	double kemax=800;

	//reco ------------------------------------------------------------------------------------------
	//end-point//--
	//inc_end
	TH1D *ke_reco_All=new TH1D("ke_reco_All","", nke, kemin, kemax);
	TH1D *ke_reco_All_inel=new TH1D("ke_reco_All_inel","", nke, kemin, kemax);
	TH1D *ke_reco_All_el=new TH1D("ke_reco_All_el","", nke, kemin, kemax);
	TH1D *ke_reco_All_midcosmic=new TH1D("ke_reco_All_midcosmic","", nke, kemin, kemax);
	TH1D *ke_reco_All_midpi=new TH1D("ke_reco_All_midpi","", nke, kemin, kemax);
	TH1D *ke_reco_All_midp=new TH1D("ke_reco_All_midp","", nke, kemin, kemax);
	TH1D *ke_reco_All_midmu=new TH1D("ke_reco_All_midmu","", nke, kemin, kemax);
	TH1D *ke_reco_All_mideg=new TH1D("ke_reco_All_mideg","", nke, kemin, kemax);
	TH1D *ke_reco_All_midother=new TH1D("ke_reco_All_midother","", nke, kemin, kemax);

	TH1D *ke_reco_All_RecoStop=new TH1D("ke_reco_All_RecoStop","", nke, kemin, kemax);
	TH2D *trklen_ke_reco_All=new TH2D("trklen_ke_reco_All","", 140, 0, 140, nke, kemin, kemax);


	//reco_inel
	TH1D *ke_reco_RecoInEl=new TH1D("ke_reco_RecoInEl","", nke, kemin, kemax);
	TH1D *ke_reco_RecoInEl_inel=new TH1D("ke_reco_RecoInEl_inel","", nke, kemin, kemax);
	TH1D *ke_reco_RecoInEl_el=new TH1D("ke_reco_RecoInEl_el","", nke, kemin, kemax);
	TH1D *ke_reco_RecoInEl_midcosmic=new TH1D("ke_reco_RecoInEl_midcosmic","", nke, kemin, kemax);
	TH1D *ke_reco_RecoInEl_midpi=new TH1D("ke_reco_RecoInEl_midpi","", nke, kemin, kemax);
	TH1D *ke_reco_RecoInEl_midp=new TH1D("ke_reco_RecoInEl_midp","", nke, kemin, kemax);
	TH1D *ke_reco_RecoInEl_midmu=new TH1D("ke_reco_RecoInEl_midmu","", nke, kemin, kemax);
	TH1D *ke_reco_RecoInEl_mideg=new TH1D("ke_reco_RecoInEl_mideg","", nke, kemin, kemax);
	TH1D *ke_reco_RecoInEl_midother=new TH1D("ke_reco_RecoInEl_midother","", nke, kemin, kemax);

	//reco_el
	TH1D *ke_reco_RecoEl=new TH1D("ke_reco_RecoEl","", nke, kemin, kemax);
	TH1D *ke_reco_RecoEl_inel=new TH1D("ke_reco_RecoEl_inel","", nke, kemin, kemax);
	TH1D *ke_reco_RecoEl_el=new TH1D("ke_reco_RecoEl_el","", nke, kemin, kemax);
	TH1D *ke_reco_RecoEl_midcosmic=new TH1D("ke_reco_RecoEl_midcosmic","", nke, kemin, kemax);
	TH1D *ke_reco_RecoEl_midpi=new TH1D("ke_reco_RecoEl_midpi","", nke, kemin, kemax);
	TH1D *ke_reco_RecoEl_midp=new TH1D("ke_reco_RecoEl_midp","", nke, kemin, kemax);
	TH1D *ke_reco_RecoEl_midmu=new TH1D("ke_reco_RecoEl_midmu","", nke, kemin, kemax);
	TH1D *ke_reco_RecoEl_mideg=new TH1D("ke_reco_RecoEl_mideg","", nke, kemin, kemax);
	TH1D *ke_reco_RecoEl_midother=new TH1D("ke_reco_RecoEl_midother","", nke, kemin, kemax);

	TH1D *keHD_reco_RecoEl=new TH1D("keHD_reco_RecoEl","", 20000, -100, 700);

	//misidp-rich
	TH1D *ke_reco_MidP=new TH1D("ke_reco_MidP","", nke, kemin, kemax);
	TH1D *ke_reco_MidP_inel=new TH1D("ke_reco_MidP_inel","", nke, kemin, kemax);
	TH1D *ke_reco_MidP_el=new TH1D("ke_reco_MidP_el","", nke, kemin, kemax);
	TH1D *ke_reco_MidP_midcosmic=new TH1D("ke_reco_MidP_midcosmic","", nke, kemin, kemax);
	TH1D *ke_reco_MidP_midpi=new TH1D("ke_reco_MidP_midpi","", nke, kemin, kemax);
	TH1D *ke_reco_MidP_midp=new TH1D("ke_reco_MidP_midp","", nke, kemin, kemax);
	TH1D *ke_reco_MidP_midmu=new TH1D("ke_reco_MidP_midmu","", nke, kemin, kemax);
	TH1D *ke_reco_MidP_mideg=new TH1D("ke_reco_MidP_mideg","", nke, kemin, kemax);
	TH1D *ke_reco_MidP_midother=new TH1D("ke_reco_MidP_midother","", nke, kemin, kemax);

	//start//--
	//inc_end
	TH1D *keff_reco_All=new TH1D("keff_reco_All","", nke, kemin, kemax);
	TH1D *keff_reco_All_inel=new TH1D("keff_reco_All_inel","", nke, kemin, kemax);
	TH1D *keff_reco_All_el=new TH1D("keff_reco_All_el","", nke, kemin, kemax);
	TH1D *keff_reco_All_midcosmic=new TH1D("keff_reco_All_midcosmic","", nke, kemin, kemax);
	TH1D *keff_reco_All_midpi=new TH1D("keff_reco_All_midpi","", nke, kemin, kemax);
	TH1D *keff_reco_All_midp=new TH1D("keff_reco_All_midp","", nke, kemin, kemax);
	TH1D *keff_reco_All_midmu=new TH1D("keff_reco_All_midmu","", nke, kemin, kemax);
	TH1D *keff_reco_All_mideg=new TH1D("keff_reco_All_mideg","", nke, kemin, kemax);
	TH1D *keff_reco_All_midother=new TH1D("keff_reco_All_midother","", nke, kemin, kemax);

	//reco_inel
	TH1D *keff_reco_RecoInEl=new TH1D("keff_reco_RecoInEl","", nke, kemin, kemax);
	TH1D *keff_reco_RecoInEl_inel=new TH1D("keff_reco_RecoInEl_inel","", nke, kemin, kemax);
	TH1D *keff_reco_RecoInEl_el=new TH1D("keff_reco_RecoInEl_el","", nke, kemin, kemax);
	TH1D *keff_reco_RecoInEl_midcosmic=new TH1D("keff_reco_RecoInEl_midcosmic","", nke, kemin, kemax);
	TH1D *keff_reco_RecoInEl_midpi=new TH1D("keff_reco_RecoInEl_midpi","", nke, kemin, kemax);
	TH1D *keff_reco_RecoInEl_midp=new TH1D("keff_reco_RecoInEl_midp","", nke, kemin, kemax);
	TH1D *keff_reco_RecoInEl_midmu=new TH1D("keff_reco_RecoInEl_midmu","", nke, kemin, kemax);
	TH1D *keff_reco_RecoInEl_mideg=new TH1D("keff_reco_RecoInEl_mideg","", nke, kemin, kemax);
	TH1D *keff_reco_RecoInEl_midother=new TH1D("keff_reco_RecoInEl_midother","", nke, kemin, kemax);

	//reco_el
	TH1D *keff_reco_RecoEl=new TH1D("keff_reco_RecoEl","", nke, kemin, kemax);
	TH1D *keff_reco_RecoEl_inel=new TH1D("keff_reco_RecoEl_inel","", nke, kemin, kemax);
	TH1D *keff_reco_RecoEl_el=new TH1D("keff_reco_RecoEl_el","", nke, kemin, kemax);
	TH1D *keff_reco_RecoEl_midcosmic=new TH1D("keff_reco_RecoEl_midcosmic","", nke, kemin, kemax);
	TH1D *keff_reco_RecoEl_midpi=new TH1D("keff_reco_RecoEl_midpi","", nke, kemin, kemax);
	TH1D *keff_reco_RecoEl_midp=new TH1D("keff_reco_RecoEl_midp","", nke, kemin, kemax);
	TH1D *keff_reco_RecoEl_midmu=new TH1D("keff_reco_RecoEl_midmu","", nke, kemin, kemax);
	TH1D *keff_reco_RecoEl_mideg=new TH1D("keff_reco_RecoEl_mideg","", nke, kemin, kemax);
	TH1D *keff_reco_RecoEl_midother=new TH1D("keff_reco_RecoEl_midother","", nke, kemin, kemax);

	//misidp-rich
	TH1D *keff_reco_MidP=new TH1D("keff_reco_MidP","", nke, kemin, kemax);
	TH1D *keff_reco_MidP_inel=new TH1D("keff_reco_MidP_inel","", nke, kemin, kemax);
	TH1D *keff_reco_MidP_el=new TH1D("keff_reco_MidP_el","", nke, kemin, kemax);
	TH1D *keff_reco_MidP_midcosmic=new TH1D("keff_reco_MidP_midcosmic","", nke, kemin, kemax);
	TH1D *keff_reco_MidP_midpi=new TH1D("keff_reco_MidP_midpi","", nke, kemin, kemax);
	TH1D *keff_reco_MidP_midp=new TH1D("keff_reco_MidP_midp","", nke, kemin, kemax);
	TH1D *keff_reco_MidP_midmu=new TH1D("keff_reco_MidP_midmu","", nke, kemin, kemax);
	TH1D *keff_reco_MidP_mideg=new TH1D("keff_reco_MidP_mideg","", nke, kemin, kemax);
	TH1D *keff_reco_MidP_midother=new TH1D("keff_reco_MidP_midother","", nke, kemin, kemax);


	TH2D *ke_kebeam_reco_RecoInEl=new TH2D("ke_kebeam_reco_RecoInEl","", nke, kemin, kemax, nke, kemin, kemax);
	TH2D *ke_kebeam_reco_RecoEl=new TH2D("ke_kebeam_reco_RecoEl","", nke, kemin, kemax, nke, kemin, kemax);



	//true ------------------------------------------------------------------------------------------
	//end-point//--
	//inc_end
	TH1D *ke_true_All=new TH1D("ke_true_All","", nke, kemin, kemax);
	TH1D *ke_true_All_inel=new TH1D("ke_true_All_inel","", nke, kemin, kemax);
	TH1D *ke_true_All_el=new TH1D("ke_true_All_el","", nke, kemin, kemax);
	TH1D *ke_true_All_midcosmic=new TH1D("ke_true_All_midcosmic","", nke, kemin, kemax);
	TH1D *ke_true_All_midpi=new TH1D("ke_true_All_midpi","", nke, kemin, kemax);
	TH1D *ke_true_All_midp=new TH1D("ke_true_All_midp","", nke, kemin, kemax);
	TH1D *ke_true_All_midmu=new TH1D("ke_true_All_midmu","", nke, kemin, kemax);
	TH1D *ke_true_All_mideg=new TH1D("ke_true_All_mideg","", nke, kemin, kemax);
	TH1D *ke_true_All_midother=new TH1D("ke_true_All_midother","", nke, kemin, kemax);

	//reco_inel
	TH1D *ke_true_RecoInEl=new TH1D("ke_true_RecoInEl","", nke, kemin, kemax);
	TH1D *ke_true_RecoInEl_inel=new TH1D("ke_true_RecoInEl_inel","", nke, kemin, kemax);
	TH1D *ke_true_RecoInEl_el=new TH1D("ke_true_RecoInEl_el","", nke, kemin, kemax);
	TH1D *ke_true_RecoInEl_midcosmic=new TH1D("ke_true_RecoInEl_midcosmic","", nke, kemin, kemax);
	TH1D *ke_true_RecoInEl_midpi=new TH1D("ke_true_RecoInEl_midpi","", nke, kemin, kemax);
	TH1D *ke_true_RecoInEl_midp=new TH1D("ke_true_RecoInEl_midp","", nke, kemin, kemax);
	TH1D *ke_true_RecoInEl_midmu=new TH1D("ke_true_RecoInEl_midmu","", nke, kemin, kemax);
	TH1D *ke_true_RecoInEl_mideg=new TH1D("ke_true_RecoInEl_mideg","", nke, kemin, kemax);
	TH1D *ke_true_RecoInEl_midother=new TH1D("ke_true_RecoInEl_midother","", nke, kemin, kemax);

	//reco_el
	TH1D *ke_true_RecoEl=new TH1D("ke_true_RecoEl","", nke, kemin, kemax);
	TH1D *ke_true_RecoEl_inel=new TH1D("ke_true_RecoEl_inel","", nke, kemin, kemax);
	TH1D *ke_true_RecoEl_el=new TH1D("ke_true_RecoEl_el","", nke, kemin, kemax);
	TH1D *ke_true_RecoEl_midcosmic=new TH1D("ke_true_RecoEl_midcosmic","", nke, kemin, kemax);
	TH1D *ke_true_RecoEl_midpi=new TH1D("ke_true_RecoEl_midpi","", nke, kemin, kemax);
	TH1D *ke_true_RecoEl_midp=new TH1D("ke_true_RecoEl_midp","", nke, kemin, kemax);
	TH1D *ke_true_RecoEl_midmu=new TH1D("ke_true_RecoEl_midmu","", nke, kemin, kemax);
	TH1D *ke_true_RecoEl_mideg=new TH1D("ke_true_RecoEl_mideg","", nke, kemin, kemax);
	TH1D *ke_true_RecoEl_midother=new TH1D("ke_true_RecoEl_midother","", nke, kemin, kemax);

	//misidp-rich
	TH1D *ke_true_MidP=new TH1D("ke_true_MidP","", nke, kemin, kemax);
	TH1D *ke_true_MidP_inel=new TH1D("ke_true_MidP_inel","", nke, kemin, kemax);
	TH1D *ke_true_MidP_el=new TH1D("ke_true_MidP_el","", nke, kemin, kemax);
	TH1D *ke_true_MidP_midcosmic=new TH1D("ke_true_MidP_midcosmic","", nke, kemin, kemax);
	TH1D *ke_true_MidP_midpi=new TH1D("ke_true_MidP_midpi","", nke, kemin, kemax);
	TH1D *ke_true_MidP_midp=new TH1D("ke_true_MidP_midp","", nke, kemin, kemax);
	TH1D *ke_true_MidP_midmu=new TH1D("ke_true_MidP_midmu","", nke, kemin, kemax);
	TH1D *ke_true_MidP_mideg=new TH1D("ke_true_MidP_mideg","", nke, kemin, kemax);
	TH1D *ke_true_MidP_midother=new TH1D("ke_true_MidP_midother","", nke, kemin, kemax);

	//start//--
	//inc_end
	TH1D *keff_true_All=new TH1D("keff_true_All","", nke, kemin, kemax);
	TH1D *keff_true_All_inel=new TH1D("keff_true_All_inel","", nke, kemin, kemax);
	TH1D *keff_true_All_el=new TH1D("keff_true_All_el","", nke, kemin, kemax);
	TH1D *keff_true_All_midcosmic=new TH1D("keff_true_All_midcosmic","", nke, kemin, kemax);
	TH1D *keff_true_All_midpi=new TH1D("keff_true_All_midpi","", nke, kemin, kemax);
	TH1D *keff_true_All_midp=new TH1D("keff_true_All_midp","", nke, kemin, kemax);
	TH1D *keff_true_All_midmu=new TH1D("keff_true_All_midmu","", nke, kemin, kemax);
	TH1D *keff_true_All_mideg=new TH1D("keff_true_All_mideg","", nke, kemin, kemax);
	TH1D *keff_true_All_midother=new TH1D("keff_true_All_midother","", nke, kemin, kemax);

	//reco_inel
	TH1D *keff_true_RecoInEl=new TH1D("keff_true_RecoInEl","", nke, kemin, kemax);
	TH1D *keff_true_RecoInEl_inel=new TH1D("keff_true_RecoInEl_inel","", nke, kemin, kemax);
	TH1D *keff_true_RecoInEl_el=new TH1D("keff_true_RecoInEl_el","", nke, kemin, kemax);
	TH1D *keff_true_RecoInEl_midcosmic=new TH1D("keff_true_RecoInEl_midcosmic","", nke, kemin, kemax);
	TH1D *keff_true_RecoInEl_midpi=new TH1D("keff_true_RecoInEl_midpi","", nke, kemin, kemax);
	TH1D *keff_true_RecoInEl_midp=new TH1D("keff_true_RecoInEl_midp","", nke, kemin, kemax);
	TH1D *keff_true_RecoInEl_midmu=new TH1D("keff_true_RecoInEl_midmu","", nke, kemin, kemax);
	TH1D *keff_true_RecoInEl_mideg=new TH1D("keff_true_RecoInEl_mideg","", nke, kemin, kemax);
	TH1D *keff_true_RecoInEl_midother=new TH1D("keff_true_RecoInEl_midother","", nke, kemin, kemax);

	//reco_el
	TH1D *keff_true_RecoEl=new TH1D("keff_true_RecoEl","", nke, kemin, kemax);
	TH1D *keff_true_RecoEl_inel=new TH1D("keff_true_RecoEl_inel","", nke, kemin, kemax);
	TH1D *keff_true_RecoEl_el=new TH1D("keff_true_RecoEl_el","", nke, kemin, kemax);
	TH1D *keff_true_RecoEl_midcosmic=new TH1D("keff_true_RecoEl_midcosmic","", nke, kemin, kemax);
	TH1D *keff_true_RecoEl_midpi=new TH1D("keff_true_RecoEl_midpi","", nke, kemin, kemax);
	TH1D *keff_true_RecoEl_midp=new TH1D("keff_true_RecoEl_midp","", nke, kemin, kemax);
	TH1D *keff_true_RecoEl_midmu=new TH1D("keff_true_RecoEl_midmu","", nke, kemin, kemax);
	TH1D *keff_true_RecoEl_mideg=new TH1D("keff_true_RecoEl_mideg","", nke, kemin, kemax);
	TH1D *keff_true_RecoEl_midother=new TH1D("keff_true_RecoEl_midother","", nke, kemin, kemax);

	//misidp-rich
	TH1D *keff_true_MidP=new TH1D("keff_true_MidP","", nke, kemin, kemax);
	TH1D *keff_true_MidP_inel=new TH1D("keff_true_MidP_inel","", nke, kemin, kemax);
	TH1D *keff_true_MidP_el=new TH1D("keff_true_MidP_el","", nke, kemin, kemax);
	TH1D *keff_true_MidP_midcosmic=new TH1D("keff_true_MidP_midcosmic","", nke, kemin, kemax);
	TH1D *keff_true_MidP_midpi=new TH1D("keff_true_MidP_midpi","", nke, kemin, kemax);
	TH1D *keff_true_MidP_midp=new TH1D("keff_true_MidP_midp","", nke, kemin, kemax);
	TH1D *keff_true_MidP_midmu=new TH1D("keff_true_MidP_midmu","", nke, kemin, kemax);
	TH1D *keff_true_MidP_mideg=new TH1D("keff_true_MidP_mideg","", nke, kemin, kemax);
	TH1D *keff_true_MidP_midother=new TH1D("keff_true_MidP_midother","", nke, kemin, kemax);

	//Booking histograms -------------------------------------------------------------------------------//

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
	//-----------------------------------------------------------------//

	//int nxx=250;	
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
	int index_minchi2=13331; //index of minchi2
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







	//Name of output file ------------------------------------------------------------------------------------------------------------//
	//TString str_out=Form("mc_kebbbkg_nobmrw.root"); //allow ke<0 and set ke=-700 if under-estimation of keff
	//TString str_out=Form("mc_kebbbkg_bmrw.root"); //allow ke<0 and set ke=-700 if under-estimation of keff
	//TString str_out=Form("mc_kebbbkg_nobmrw_HD.root"); //allow ke<0 and set ke=-700 if under-estimation of keff
	TString str_out=Form("mc_kebbbkg_bmrw_HD.root"); //allow ke<0 and set ke=-700 if under-estimation of keff

	//Basic configure ------//
	BetheBloch BB;
	BB.SetPdgCode(pdg);

	//book histograms --//
	//BookHistograms();

	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	Long64_t nbytes = 0, nb = 0;
	bool isTestSample=true;
	int true_sliceID = -1, reco_sliceID = -1;
	//int true_st_sliceID = -1, reco_st_sliceID = -1;
	for (Long64_t jentry=0; jentry<nentries;jentry++) { //main entry loop
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		isTestSample = true;
		if (ientry%2 == 0) isTestSample = false; //Divide MC sample by 2 parts: test+ufold
		//if (isTestSample) continue; //only validate sample
		//if (!isTestSample) continue; //only test sample

		//true_sliceID = -1;
		//reco_sliceID = -1;

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
			if (key_reach_tpc!=-99) {
				key_fit_st=key_reach_tpc-1;
				key_fit_ed=key_reach_tpc+1;
			}
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
			}
		} //is beam at ff						    	

		//fix on the truth length by adding distance between 1st tpc hit to front face ------------------------------------------------------//
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
			//reco_cosineTheta->Fill(cosine_beam_spec_primtrk);
			//Fill1DHist(reco_cosineTheta, cosine_beam_spec_primtrk);
			//if (kinel) reco_cosineTheta_inel->Fill(cosine_beam_spec_primtrk);
			//if (kel) reco_cosineTheta_el->Fill(cosine_beam_spec_primtrk);
			//if (kMIDcosmic) reco_cosineTheta_midcosmic->Fill(cosine_beam_spec_primtrk);
			//if (kMIDpi) reco_cosineTheta_midpi->Fill(cosine_beam_spec_primtrk);
			//if (kMIDp) reco_cosineTheta_midp->Fill(cosine_beam_spec_primtrk);
			//if (kMIDmu) reco_cosineTheta_midmu->Fill(cosine_beam_spec_primtrk);
			//if (kMIDeg) reco_cosineTheta_mideg->Fill(cosine_beam_spec_primtrk);
			//if (kMIDother) reco_cosineTheta_midother->Fill(cosine_beam_spec_primtrk);

			//if (kinel) Fill1DHist(reco_cosineTheta_inel, cosine_beam_spec_primtrk);
			//if (kel) Fill1DHist(reco_cosineTheta_el, cosine_beam_spec_primtrk);
			//if (kMIDcosmic) Fill1DHist(reco_cosineTheta_midcosmic, cosine_beam_spec_primtrk);
			//if (kMIDpi) Fill1DHist(reco_cosineTheta_midpi, cosine_beam_spec_primtrk);
			//if (kMIDp) Fill1DHist(reco_cosineTheta_midp, cosine_beam_spec_primtrk);
			//if (kMIDmu) Fill1DHist(reco_cosineTheta_midmu, cosine_beam_spec_primtrk);
			//if (kMIDeg) Fill1DHist(reco_cosineTheta_mideg, cosine_beam_spec_primtrk);
			//if (kMIDother) Fill1DHist(reco_cosineTheta_midother, cosine_beam_spec_primtrk);
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

		//beam quality cut --------------//
		bool IsBQ=false;
		if (IsCosine&&IsPos) IsBQ=true;

		bool IsMisidpRich=false;
		if (IsPos&&IsCaloSize&&IsPandoraSlice) {
			if (cosine_beam_spec_primtrk<=0.9) IsMisidpRich=true;
		}


		//reco calorimetry ---------------------------------------------------------------------------//
		int index_reco_endz=0;
		double wid_reco_max=-9999;
		double range_reco=-999;
		vector<double> reco_trklen_accum;
		reco_trklen_accum.reserve(primtrk_hitz->size());
		double reco_calo_MeV=0;
		//double kereco_range=0;
		//double kereco_range2=0;
		vector<double> EDept;
		vector<double> DEDX;
		vector<double> DX;
		double pid=-99; 
		if (IsCaloSize) { //if calo size not empty
			vector<double> trkdedx;
			vector<double> trkres;
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

				if (h==1) range_reco=0;
				if (h>=1) {
					range_reco += sqrt( pow(primtrk_hitx->at(h)-primtrk_hitx->at(h-1), 2)+
							pow(primtrk_hity->at(h)-primtrk_hity->at(h-1), 2)+
							pow(primtrk_hitz->at(h)-primtrk_hitz->at(h-1), 2) );
					reco_trklen_accum[h] = range_reco;
				}

				reco_calo_MeV+=cali_dedx*pitch;
				//kereco_range+=pitch*dedx_predict(resrange_reco);
				//kereco_range2+=pitch*(double)gr_predict_dedx_resrange->Eval(resrange_reco);

				//if (kinel) rangereco_dedxreco_TrueInEL->Fill(range_reco, cali_dedx);
				//if (kel) { 
				//rangereco_dedxreco_TrueEL->Fill(range_reco, cali_dedx);
				//rr_dedx_truestop->Fill(resrange_reco, cali_dedx);
				//}

				trkdedx.push_back(cali_dedx);
				trkres.push_back(resrange_reco);

			} //loop over reco hits of a given track

			pid=chi2pid(trkdedx,trkres); //pid using stopping proton hypothesis

		} //if calo size not empty

		//Reco stopping/Inel p cut ---------------------------------------------------------------------------------------------------------//
		bool IsRecoStop=false;
		bool IsRecoInEL=false;
		bool IsRecoEL=false;
		double mom_beam_spec=-99; mom_beam_spec=beamMomentum_spec->at(0);
		//double range_reco=-99; if (!primtrk_range->empty()) range_reco=primtrk_range->at(0); //reco primary trklen

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
		//if (pid>pid_2) IsRecoInEL=true; 
		//if (pid<=pid_2) IsRecoEL=true;

		//kinetic energies -------------------------------------------------------------------//
		//double ke_beam=1000.*p2ke(mom_beam); //ke_beam
		double ke_beam_spec=p2ke(mom_beam_spec); //ke_beam_spec [GeV]
		double ke_beam_spec_MeV=1000.*ke_beam_spec; //ke_beam_spec [MeV]
		double ke_trklen=1000.*ke_vs_csda_range_sm->Eval(range_reco); //[unit: MeV]
		double p_trklen=ke2p(ke_trklen);
		double ke_simide=0;
		for (int hk=0; hk<(int)primtrk_true_edept->size(); ++hk) { //loop over simIDE points
			ke_simide+=primtrk_true_edept->at(hk);
		} //loop over simIDE points

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
		double KE_ff_true=KE_ff;

		//KEff (reco) with const E-loss assumption -----------------------------------//
		//double mean_Elosscalo_stop=(4.95958e+01)/(1.00489e+00); //using fit [no bmrw]
		//double KE_ff_reco=ke_beam_spec_MeV-mean_Elosscalo_stop;
		double KE_ff_reco=ke_beam_spec_MeV-mean_Elossrange_stop; //kebb
		double KEend_reco=0;
		//KEend_reco=KE_ff_reco-reco_calo_MeV;		
		KEend_reco=BB.KEAtLength(KE_ff_reco, range_reco);		

		//KEend ---------------------------------------------------------------------------//
		double KEend_true=0;
		if (beamtrk_Eng->size()) KEend_true=1000.*(beamtrk_Eng->at(-2+beamtrk_Eng->size()));

		//KEs ---------------------------------------------------------------------------------------//
		//double Eloss_upstream=0; 
		//if (KE_ff>0) Eloss_upstream=
		//double dEbb_true=0; if (range_true>=0&&KE_ff>0) dEbb_true=BB.KEAtLength(KE_ff, range_true);
		//double dEbb_reco=0; if (range_reco>=0&&KE_ff>0) dEbb_reco=BB.KEAtLength(KE_ff, range_reco);

		//double KEbb_true=-1; KEbb_true=KE_ff-BB.KEAtLength(KE_ff, range_true);
		//double KEbb_reco=-1; KEbb_reco=KE_ff-BB.KEAtLength(KE_ff, range_reco);
		//---------------------------------------------------------------------------------------------------------------//

		//bmrw -------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
		double mom_rw_minchi2=1; //weight for beam-momentum-reweight
		//if ((ke_beam_spec_MeV-mean_Elosscalo_stop)>=mu_min&&(ke_beam_spec_MeV-mean_Elosscalo_stop)<=mu_max) mom_rw_minchi2=agng->Eval(ke_beam_spec_MeV-mean_Elosscalo_stop); //bmrw
		if ((mom_beam_spec*1000.)>=mmu_min&&(mom_beam_spec*1000.)<=mmu_max) mom_rw_minchi2=gng[index_minchi2]->Eval(mom_beam_spec*1000.); //bmrw

		//Fill histograms -------------------------------------------------------------------------------------------//
		if (IsPandoraSlice&&IsCaloSize&&IsBQ) {  //basic cuts
			//inc
			Fill1DWHist(ke_reco_All, KEend_reco, mom_rw_minchi2);
			Fill1DWHist(ke_true_All, KEend_true, mom_rw_minchi2);
			trklen_ke_reco_All->Fill(range_reco, KEend_reco, mom_rw_minchi2);

			Fill1DWHist(ke_reco_All_RecoStop, ke_trklen, mom_rw_minchi2);

			if (kinel) { 
				Fill1DWHist(ke_reco_All_inel, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_All_inel, KEend_true, mom_rw_minchi2);
			}
			if (kel) { 
				Fill1DWHist(ke_reco_All_el, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_All_el, KEend_true, mom_rw_minchi2);
			}
			if (kMIDcosmic) { 
				Fill1DWHist(ke_reco_All_midcosmic, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_All_midcosmic, KEend_true, mom_rw_minchi2);
			}
			if (kMIDpi) { 
				Fill1DWHist(ke_reco_All_midpi, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_All_midpi, KEend_true, mom_rw_minchi2);
			}
			if (kMIDp) { 
				Fill1DWHist(ke_reco_All_midp, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_All_midp, KEend_true, mom_rw_minchi2);
			}
			if (kMIDmu) { 
				Fill1DWHist(ke_reco_All_midmu, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_All_midmu, KEend_true, mom_rw_minchi2);
			}
			if (kMIDeg) { 
				Fill1DWHist(ke_reco_All_mideg, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_All_mideg, KEend_true, mom_rw_minchi2);
			}
			if (kMIDother) { 
				Fill1DWHist(ke_reco_All_midother, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_All_midother, KEend_true, mom_rw_minchi2);
			}



			if (IsRecoInEL) { //reco inel
				Fill1DWHist(ke_reco_RecoInEl, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_RecoInEl, KEend_true, mom_rw_minchi2);
				ke_kebeam_reco_RecoInEl->Fill(KEend_reco, ke_beam_spec_MeV);

				if (kinel) {
					Fill1DWHist(ke_reco_RecoInEl_inel, KEend_reco, mom_rw_minchi2);
					Fill1DWHist(ke_true_RecoInEl_inel, KEend_true, mom_rw_minchi2);
				}
				if (kel) { 
					Fill1DWHist(ke_reco_RecoInEl_el, KEend_reco, mom_rw_minchi2);
					Fill1DWHist(ke_true_RecoInEl_el, KEend_true, mom_rw_minchi2);
				}
				if (kMIDcosmic) { 
					Fill1DWHist(ke_reco_RecoInEl_midcosmic, KEend_reco, mom_rw_minchi2);
					Fill1DWHist(ke_true_RecoInEl_midcosmic, KEend_true, mom_rw_minchi2);
				}
				if (kMIDpi) { 
					Fill1DWHist(ke_reco_RecoInEl_midpi, KEend_reco, mom_rw_minchi2);
					Fill1DWHist(ke_true_RecoInEl_midpi, KEend_true, mom_rw_minchi2);
				}
				if (kMIDp) { 
					Fill1DWHist(ke_reco_RecoInEl_midp, KEend_reco, mom_rw_minchi2);
					Fill1DWHist(ke_true_RecoInEl_midp, KEend_true, mom_rw_minchi2);
				}
				if (kMIDmu) { 
					Fill1DWHist(ke_reco_RecoInEl_midmu, KEend_reco, mom_rw_minchi2);
					Fill1DWHist(ke_true_RecoInEl_midmu, KEend_true, mom_rw_minchi2);
				}
				if (kMIDeg) { 
					Fill1DWHist(ke_reco_RecoInEl_mideg, KEend_reco, mom_rw_minchi2);
					Fill1DWHist(ke_true_RecoInEl_mideg, KEend_true, mom_rw_minchi2);
				}
				if (kMIDother) { 
					Fill1DWHist(ke_reco_RecoInEl_midother, KEend_reco, mom_rw_minchi2);
					Fill1DWHist(ke_true_RecoInEl_midother, KEend_true, mom_rw_minchi2);
				}

			} //reco inel

			if (IsRecoEL) { //reco el
				Fill1DWHist(ke_reco_RecoEl, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(keHD_reco_RecoEl, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_RecoEl, KEend_true, mom_rw_minchi2);
				ke_kebeam_reco_RecoEl->Fill(KEend_reco, ke_beam_spec_MeV);
				if (kinel) {
					Fill1DWHist(ke_reco_RecoEl_inel, KEend_reco, mom_rw_minchi2);
					Fill1DWHist(ke_true_RecoEl_inel, KEend_true, mom_rw_minchi2);
				}
				if (kel) { 
					Fill1DWHist(ke_reco_RecoEl_el, KEend_reco, mom_rw_minchi2);
					Fill1DWHist(ke_true_RecoEl_el, KEend_true, mom_rw_minchi2);
				}
				if (kMIDcosmic) { 
					Fill1DWHist(ke_reco_RecoEl_midcosmic, KEend_reco, mom_rw_minchi2);
					Fill1DWHist(ke_true_RecoEl_midcosmic, KEend_true, mom_rw_minchi2);
				}
				if (kMIDpi) { 
					Fill1DWHist(ke_reco_RecoEl_midpi, KEend_reco, mom_rw_minchi2);
					Fill1DWHist(ke_true_RecoEl_midpi, KEend_true, mom_rw_minchi2);
				}
				if (kMIDp) { 
					Fill1DWHist(ke_reco_RecoEl_midp, KEend_reco, mom_rw_minchi2);
					Fill1DWHist(ke_true_RecoEl_midp, KEend_true, mom_rw_minchi2);
				}
				if (kMIDmu) { 
					Fill1DWHist(ke_reco_RecoEl_midmu, KEend_reco, mom_rw_minchi2);
					Fill1DWHist(ke_true_RecoEl_midmu, KEend_true, mom_rw_minchi2);
				}
				if (kMIDeg) { 
					Fill1DWHist(ke_reco_RecoEl_mideg, KEend_reco, mom_rw_minchi2);
					Fill1DWHist(ke_true_RecoEl_mideg, KEend_true, mom_rw_minchi2);
				}
				if (kMIDother) { 
					Fill1DWHist(ke_reco_RecoEl_midother, KEend_reco, mom_rw_minchi2);
					Fill1DWHist(ke_true_RecoEl_midother, KEend_true, mom_rw_minchi2);
				}

			} //reco el
		} //basic cuts


		if (IsMisidpRich) { //misidp-rich
			Fill1DWHist(ke_reco_MidP, KEend_reco, mom_rw_minchi2);
			Fill1DWHist(ke_true_MidP, KEend_true, mom_rw_minchi2);
			if (kinel) {
				Fill1DWHist(ke_reco_MidP_inel, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_MidP_inel, KEend_true, mom_rw_minchi2);
			}
			if (kel) { 
				Fill1DWHist(ke_reco_MidP_el, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_MidP_el, KEend_true, mom_rw_minchi2);
			}
			if (kMIDcosmic) { 
				Fill1DWHist(ke_reco_MidP_midcosmic, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_MidP_midcosmic, KEend_true, mom_rw_minchi2);
			}
			if (kMIDpi) { 
				Fill1DWHist(ke_reco_MidP_midpi, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_MidP_midpi, KEend_true, mom_rw_minchi2);
			}
			if (kMIDp) { 
				Fill1DWHist(ke_reco_MidP_midp, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_MidP_midp, KEend_true, mom_rw_minchi2);
			}
			if (kMIDmu) { 
				Fill1DWHist(ke_reco_MidP_midmu, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_MidP_midmu, KEend_true, mom_rw_minchi2);
			}
			if (kMIDeg) { 
				Fill1DWHist(ke_reco_MidP_mideg, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_MidP_mideg, KEend_true, mom_rw_minchi2);
			}
			if (kMIDother) { 
				Fill1DWHist(ke_reco_MidP_midother, KEend_reco, mom_rw_minchi2);
				Fill1DWHist(ke_true_MidP_midother, KEend_true, mom_rw_minchi2);
			}
		} //misidp-rich


		//ff
		if (IsPandoraSlice&&IsCaloSize&&IsBQ) {  //basic cuts
			//inc
			Fill1DWHist(keff_reco_All, KE_ff_reco, mom_rw_minchi2);
			Fill1DWHist(keff_true_All, KE_ff_true, mom_rw_minchi2);

			if (kinel) { 
				Fill1DWHist(keff_reco_All_inel, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_All_inel, KE_ff_true, mom_rw_minchi2);
			}
			if (kel) { 
				Fill1DWHist(keff_reco_All_el, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_All_el, KE_ff_true, mom_rw_minchi2);
			}
			if (kMIDcosmic) { 
				Fill1DWHist(keff_reco_All_midcosmic, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_All_midcosmic, KE_ff_true, mom_rw_minchi2);
			}
			if (kMIDpi) { 
				Fill1DWHist(keff_reco_All_midpi, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_All_midpi, KE_ff_true, mom_rw_minchi2);
			}
			if (kMIDp) { 
				Fill1DWHist(keff_reco_All_midp, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_All_midp, KE_ff_true, mom_rw_minchi2);
			}
			if (kMIDmu) { 
				Fill1DWHist(keff_reco_All_midmu, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_All_midmu, KE_ff_true, mom_rw_minchi2);
			}
			if (kMIDeg) { 
				Fill1DWHist(keff_reco_All_mideg, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_All_mideg, KE_ff_true, mom_rw_minchi2);
			}
			if (kMIDother) { 
				Fill1DWHist(keff_reco_All_midother, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_All_midother, KE_ff_true, mom_rw_minchi2);
			}

			if (IsRecoInEL) { //reco inel
				Fill1DWHist(keff_reco_RecoInEl, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_RecoInEl, KE_ff_true, mom_rw_minchi2);
				if (kinel) {
					Fill1DWHist(keff_reco_RecoInEl_inel, KE_ff_reco, mom_rw_minchi2);
					Fill1DWHist(keff_true_RecoInEl_inel, KE_ff_true, mom_rw_minchi2);
				}
				if (kel) { 
					Fill1DWHist(keff_reco_RecoInEl_el, KE_ff_reco, mom_rw_minchi2);
					Fill1DWHist(keff_true_RecoInEl_el, KE_ff_true, mom_rw_minchi2);
				}
				if (kMIDcosmic) { 
					Fill1DWHist(keff_reco_RecoInEl_midcosmic, KE_ff_reco, mom_rw_minchi2);
					Fill1DWHist(keff_true_RecoInEl_midcosmic, KE_ff_true, mom_rw_minchi2);
				}
				if (kMIDpi) { 
					Fill1DWHist(keff_reco_RecoInEl_midpi, KE_ff_reco, mom_rw_minchi2);
					Fill1DWHist(keff_true_RecoInEl_midpi, KE_ff_true, mom_rw_minchi2);
				}
				if (kMIDp) { 
					Fill1DWHist(keff_reco_RecoInEl_midp, KE_ff_reco, mom_rw_minchi2);
					Fill1DWHist(keff_true_RecoInEl_midp, KE_ff_true, mom_rw_minchi2);
				}
				if (kMIDmu) { 
					Fill1DWHist(keff_reco_RecoInEl_midmu, KE_ff_reco, mom_rw_minchi2);
					Fill1DWHist(keff_true_RecoInEl_midmu, KE_ff_true, mom_rw_minchi2);
				}
				if (kMIDeg) { 
					Fill1DWHist(keff_reco_RecoInEl_mideg, KE_ff_reco, mom_rw_minchi2);
					Fill1DWHist(keff_true_RecoInEl_mideg, KE_ff_true, mom_rw_minchi2);
				}
				if (kMIDother) { 
					Fill1DWHist(keff_reco_RecoInEl_midother, KE_ff_reco, mom_rw_minchi2);
					Fill1DWHist(keff_true_RecoInEl_midother, KE_ff_true, mom_rw_minchi2);
				}

			} //reco inel

			if (IsRecoEL) { //reco el
				Fill1DWHist(keff_reco_RecoEl, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_RecoEl, KE_ff_true, mom_rw_minchi2);
				if (kinel) {
					Fill1DWHist(keff_reco_RecoEl_inel, KE_ff_reco, mom_rw_minchi2);
					Fill1DWHist(keff_true_RecoEl_inel, KE_ff_true, mom_rw_minchi2);
				}
				if (kel) { 
					Fill1DWHist(keff_reco_RecoEl_el, KE_ff_reco, mom_rw_minchi2);
					Fill1DWHist(keff_true_RecoEl_el, KE_ff_true, mom_rw_minchi2);
				}
				if (kMIDcosmic) { 
					Fill1DWHist(keff_reco_RecoEl_midcosmic, KE_ff_reco, mom_rw_minchi2);
					Fill1DWHist(keff_true_RecoEl_midcosmic, KE_ff_true, mom_rw_minchi2);
				}
				if (kMIDpi) { 
					Fill1DWHist(keff_reco_RecoEl_midpi, KE_ff_reco, mom_rw_minchi2);
					Fill1DWHist(keff_true_RecoEl_midpi, KE_ff_true, mom_rw_minchi2);
				}
				if (kMIDp) { 
					Fill1DWHist(keff_reco_RecoEl_midp, KE_ff_reco, mom_rw_minchi2);
					Fill1DWHist(keff_true_RecoEl_midp, KE_ff_true, mom_rw_minchi2);
				}
				if (kMIDmu) { 
					Fill1DWHist(keff_reco_RecoEl_midmu, KE_ff_reco, mom_rw_minchi2);
					Fill1DWHist(keff_true_RecoEl_midmu, KE_ff_true, mom_rw_minchi2);
				}
				if (kMIDeg) { 
					Fill1DWHist(keff_reco_RecoEl_mideg, KE_ff_reco, mom_rw_minchi2);
					Fill1DWHist(keff_true_RecoEl_mideg, KE_ff_true, mom_rw_minchi2);
				}
				if (kMIDother) { 
					Fill1DWHist(keff_reco_RecoEl_midother, KE_ff_reco, mom_rw_minchi2);
					Fill1DWHist(keff_true_RecoEl_midother, KE_ff_true, mom_rw_minchi2);
				}
			} //reco el
		} //basic cuts


		if (IsMisidpRich) { //misidp-rich
			Fill1DWHist(keff_reco_MidP, KE_ff_reco, mom_rw_minchi2);
			Fill1DWHist(keff_true_MidP, KE_ff_true, mom_rw_minchi2);
			if (kinel) {
				Fill1DWHist(keff_reco_MidP_inel, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_MidP_inel, KE_ff_true, mom_rw_minchi2);
			}
			if (kel) { 
				Fill1DWHist(keff_reco_MidP_el, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_MidP_el, KE_ff_true, mom_rw_minchi2);
			}
			if (kMIDcosmic) { 
				Fill1DWHist(keff_reco_MidP_midcosmic, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_MidP_midcosmic, KE_ff_true, mom_rw_minchi2);
			}
			if (kMIDpi) { 
				Fill1DWHist(keff_reco_MidP_midpi, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_MidP_midpi, KE_ff_true, mom_rw_minchi2);
			}
			if (kMIDp) { 
				Fill1DWHist(keff_reco_MidP_midp, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_MidP_midp, KE_ff_true, mom_rw_minchi2);
			}
			if (kMIDmu) { 
				Fill1DWHist(keff_reco_MidP_midmu, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_MidP_midmu, KE_ff_true, mom_rw_minchi2);
			}
			if (kMIDeg) { 
				Fill1DWHist(keff_reco_MidP_mideg, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_MidP_mideg, KE_ff_true, mom_rw_minchi2);
			}
			if (kMIDother) { 
				Fill1DWHist(keff_reco_MidP_midother, KE_ff_reco, mom_rw_minchi2);
				Fill1DWHist(keff_true_MidP_midother, KE_ff_true, mom_rw_minchi2);
			}
		} //misidp-rich




		} //main entry loop


		//save results -------//
		TFile *fout = new TFile(Form("%s",str_out.Data()),"RECREATE");
		trklen_ke_reco_All->Write();

		ke_true_All->Write();
		ke_true_All_inel->Write();
		ke_true_All_el->Write();
		ke_true_All_midcosmic->Write();
		ke_true_All_midpi->Write();
		ke_true_All_midp->Write();
		ke_true_All_midmu->Write();
		ke_true_All_mideg->Write();
		ke_true_All_midother->Write();

		ke_true_RecoInEl->Write();
		ke_true_RecoInEl_inel->Write();
		ke_true_RecoInEl_el->Write();
		ke_true_RecoInEl_midcosmic->Write();
		ke_true_RecoInEl_midpi->Write();
		ke_true_RecoInEl_midp->Write();
		ke_true_RecoInEl_midmu->Write();
		ke_true_RecoInEl_mideg->Write();
		ke_true_RecoInEl_midother->Write();

		ke_true_RecoEl->Write();
		ke_true_RecoEl_inel->Write();
		ke_true_RecoEl_el->Write();
		ke_true_RecoEl_midcosmic->Write();
		ke_true_RecoEl_midpi->Write();
		ke_true_RecoEl_midp->Write();
		ke_true_RecoEl_midmu->Write();
		ke_true_RecoEl_mideg->Write();
		ke_true_RecoEl_midother->Write();

		ke_true_MidP->Write();
		ke_true_MidP_inel->Write();
		ke_true_MidP_el->Write();
		ke_true_MidP_midcosmic->Write();
		ke_true_MidP_midpi->Write();
		ke_true_MidP_midp->Write();
		ke_true_MidP_midmu->Write();
		ke_true_MidP_mideg->Write();
		ke_true_MidP_midother->Write();

		keff_true_All->Write();
		keff_true_All_inel->Write();
		keff_true_All_el->Write();
		keff_true_All_midcosmic->Write();
		keff_true_All_midpi->Write();
		keff_true_All_midp->Write();
		keff_true_All_midmu->Write();
		keff_true_All_mideg->Write();
		keff_true_All_midother->Write();

		keff_true_RecoInEl->Write();
		keff_true_RecoInEl_inel->Write();
		keff_true_RecoInEl_el->Write();
		keff_true_RecoInEl_midcosmic->Write();
		keff_true_RecoInEl_midpi->Write();
		keff_true_RecoInEl_midp->Write();
		keff_true_RecoInEl_midmu->Write();
		keff_true_RecoInEl_mideg->Write();
		keff_true_RecoInEl_midother->Write();

		keff_true_RecoEl->Write();
		keff_true_RecoEl_inel->Write();
		keff_true_RecoEl_el->Write();
		keff_true_RecoEl_midcosmic->Write();
		keff_true_RecoEl_midpi->Write();
		keff_true_RecoEl_midp->Write();
		keff_true_RecoEl_midmu->Write();
		keff_true_RecoEl_mideg->Write();
		keff_true_RecoEl_midother->Write();

		keff_true_MidP->Write();
		keff_true_MidP_inel->Write();
		keff_true_MidP_el->Write();
		keff_true_MidP_midcosmic->Write();
		keff_true_MidP_midpi->Write();
		keff_true_MidP_midp->Write();
		keff_true_MidP_midmu->Write();
		keff_true_MidP_mideg->Write();
		keff_true_MidP_midother->Write();





		ke_reco_All->Write();
		ke_reco_All_inel->Write();
		ke_reco_All_el->Write();
		ke_reco_All_midcosmic->Write();
		ke_reco_All_midpi->Write();
		ke_reco_All_midp->Write();
		ke_reco_All_midmu->Write();
		ke_reco_All_mideg->Write();
		ke_reco_All_midother->Write();

		ke_reco_All_RecoStop->Write();

		ke_reco_RecoInEl->Write();
		ke_reco_RecoInEl_inel->Write();
		ke_reco_RecoInEl_el->Write();
		ke_reco_RecoInEl_midcosmic->Write();
		ke_reco_RecoInEl_midpi->Write();
		ke_reco_RecoInEl_midp->Write();
		ke_reco_RecoInEl_midmu->Write();
		ke_reco_RecoInEl_mideg->Write();
		ke_reco_RecoInEl_midother->Write();

		keHD_reco_RecoEl->Write();
		ke_reco_RecoEl->Write();
		ke_reco_RecoEl_inel->Write();
		ke_reco_RecoEl_el->Write();
		ke_reco_RecoEl_midcosmic->Write();
		ke_reco_RecoEl_midpi->Write();
		ke_reco_RecoEl_midp->Write();
		ke_reco_RecoEl_midmu->Write();
		ke_reco_RecoEl_mideg->Write();
		ke_reco_RecoEl_midother->Write();

		ke_reco_MidP->Write();
		ke_reco_MidP_inel->Write();
		ke_reco_MidP_el->Write();
		ke_reco_MidP_midcosmic->Write();
		ke_reco_MidP_midpi->Write();
		ke_reco_MidP_midp->Write();
		ke_reco_MidP_midmu->Write();
		ke_reco_MidP_mideg->Write();
		ke_reco_MidP_midother->Write();

		keff_reco_All->Write();
		keff_reco_All_inel->Write();
		keff_reco_All_el->Write();
		keff_reco_All_midcosmic->Write();
		keff_reco_All_midpi->Write();
		keff_reco_All_midp->Write();
		keff_reco_All_midmu->Write();
		keff_reco_All_mideg->Write();
		keff_reco_All_midother->Write();

		keff_reco_RecoInEl->Write();
		keff_reco_RecoInEl_inel->Write();
		keff_reco_RecoInEl_el->Write();
		keff_reco_RecoInEl_midcosmic->Write();
		keff_reco_RecoInEl_midpi->Write();
		keff_reco_RecoInEl_midp->Write();
		keff_reco_RecoInEl_midmu->Write();
		keff_reco_RecoInEl_mideg->Write();
		keff_reco_RecoInEl_midother->Write();

		keff_reco_RecoEl->Write();
		keff_reco_RecoEl_inel->Write();
		keff_reco_RecoEl_el->Write();
		keff_reco_RecoEl_midcosmic->Write();
		keff_reco_RecoEl_midpi->Write();
		keff_reco_RecoEl_midp->Write();
		keff_reco_RecoEl_midmu->Write();
		keff_reco_RecoEl_mideg->Write();
		keff_reco_RecoEl_midother->Write();

		keff_reco_MidP->Write();
		keff_reco_MidP_inel->Write();
		keff_reco_MidP_el->Write();
		keff_reco_MidP_midcosmic->Write();
		keff_reco_MidP_midpi->Write();
		keff_reco_MidP_midp->Write();
		keff_reco_MidP_midmu->Write();
		keff_reco_MidP_mideg->Write();
		keff_reco_MidP_midother->Write();

		ke_kebeam_reco_RecoInEl->Write();
		ke_kebeam_reco_RecoEl->Write();
		fout->Close();




		}
