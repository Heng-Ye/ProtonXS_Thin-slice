#define ProtonApplyMomentumReweight_cxx
#include "ProtonApplyMomentumReweight.h"

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



void ProtonApplyMomentumReweight::Loop() {
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	Long64_t nbytes = 0, nb = 0;
	bool isTestSample=true;
	//int true_sliceID = -1, reco_sliceID = -1;

	//Weighted Gaussians & Weighting function ----------------------------------------------------------------------------------------------------------------------------------------------//
	int n_b=30;
	double b_min=0;
	double b_max=150;


	int nx=250;	
	double xmin=0.; //pmin [MeV/c]
	double xmax=2000.; //pmax [MeV/c]
	TH1D *h1d_p0=new TH1D("h1d_p0","",nx,xmin,xmax);
	TH1D *h1d_p0_stop=new TH1D("h1d_p0_stop","",nx,xmin,xmax);
	TH1D *h1d_pbeam=new TH1D("h1d_pbeam","",nx,xmin,xmax);
	TH1D *h1d_pbeam_stop=new TH1D("h1d_pbeam_stop","",nx,xmin,xmax);

	TH1D *h1d_pff=new TH1D("h1d_pff","",nx,xmin,xmax);
	TH1D *h1d_pff_stop=new TH1D("h1d_pff_stop","",nx,xmin,xmax);

	TH1D *h1d_prange_stop=new TH1D("h1d_prange_stop","",nx,xmin,xmax);
	TH1D *h1d_pcalo_stop=new TH1D("h1d_pcalo_stop","",nx,xmin,xmax);

	h1d_p0->Sumw2();
	h1d_p0_stop->Sumw2();
	h1d_pbeam->Sumw2();
	h1d_pff->Sumw2();
	h1d_pff_stop->Sumw2();
	h1d_prange_stop->Sumw2();
	h1d_pcalo_stop->Sumw2();

	//MC Beam Mom Gaussian 
	//double m1=1007.1482; //MC prod4a [spec]
	//double s1=60.703307; //MC prod4a [spec]
	double m1=997.969; //MC prod4a [truth]
	double s1=54.4602; //MC prod4a [truth]


	//const E-loss using stopping protons ---------
	double Eloss_mc_hy_stop=19.542/0.998495;
	//p[0]:19.542;   //err_p[0]:0.126113
	//p[1]:0.998495; //err_p[1]:0.00549534

	double R_fit_hy=1.0008142352819318;
	double er_R_fit_hy=0.04629667706788889;


	//const. E-loss assumption	
	double const_eloss_mc=47.0058/1.00097; //const E-loss from fit (calo)
	//#p[0]:47.0058 err_p[0]:0.372157 p[1]:-1.00097 err_p[1]:0.00787403

	//New weighting func (using KEff_fit_stop at TPC FF as a reference) ------//
	//double mu_denom_data=411.06602388610895; //old
	//double sg_denom_data=47.075678784947826; //old

	double mu_denom_data=411.05145837595467; //new with event-by-event R corr
	double sg_denom_data=47.48714821962207; //new with event-by-event R corr

	//double mu_nom_data=390.81237292943916; //for data (now use event-by-event correction)
	//double sg_nom_data=47.52091718691363; //for data (now use event-by-event correction)

	//double mu_nom_mc=388.560260293186; //for mc(KEbeam-const_from_calo)
	//double sg_nom_mc=43.13168235197187; //formc

	double mu_nom_mc=416.224743039812; //for mc(KEbeam-const) with R=1
	double sg_nom_mc=42.786018962508784; //

	//i= 0  m= 416.2247430398121 s= 42.786018962508784 [kebeam-dE]*1 (R=1)
	//i= 1  m= 416.1620092367158 s= 40.48356757740762 [KE(Fit)] 

	//weighting func. (ke)
	TF1 *kerw=new TF1(Form("kerw"),govg,0,800,4);
	kerw->SetParameter(0, mu_nom_mc);
	kerw->SetParameter(1, sg_nom_mc);
	kerw->SetParameter(2, mu_denom_data);
	kerw->SetParameter(3, sg_denom_data);

	//ke cut range	
	double mu_kemin=mu_nom_mc-3.*sg_nom_mc;
	double mu_kemax=mu_nom_mc+3.*sg_nom_mc;
	//------------------------------------------------------------------------//

	//momentum cut range	
	double mu_min=1007.1482-3.*60.703307;
	double mu_max=1007.1482+3.*60.703307;

	//default gaussian
	TF1 *g1=new TF1("g1",fitg,xmin,xmax,2);
	g1->SetName("g1");
	g1->SetParameter(0,m1);
	g1->SetParameter(1,s1);


	//mu range
	double dmu=0.0005;
	double mu_st=1.01;
	//int nmu=71;
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
	//int index_minchi2=13331; //index of minchi2(use spec)
	//int index_minchi2=15494; //index of minchi2 (use MC truth)
	//int index_minchi2=3407; //index of minchi2 (use calo)
	int index_minchi2=1270; //index of minchi2 (use calo+rm track)
	
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


	//trklen
	TH1D *h1d_trklen_stop=new TH1D(Form("h1d_trklen_stop"),Form("reco stop"),n_b,b_min,b_max);
	TH1D *h1d_trklen_stop_XY=new TH1D(Form("h1d_trklen_stop_XY"),Form("reco stop with xy cut"),n_b,b_min,b_max);
	
	TH1D *h1d_trklen=new TH1D(Form("h1d_trklen"),Form("reco+BQ"),n_b,b_min,b_max);
	TH1D *h1d_trklen_XY=new TH1D(Form("h1d_trklen_XY"),Form("reco+BQ+XY"),n_b,b_min,b_max);
	TH1D *h1d_trklen_bmrw=new TH1D(Form("h1d_trklen_bmrw"),Form("reco_bmrw+BQ"),n_b,b_min,b_max);
	TH1D *h1d_trklen_bmrw_XY=new TH1D(Form("h1d_trklen_bmrw_XY"),Form("reco_bmrw+BQ+XY"),n_b,b_min,b_max);

	//zend
        int dz=350;
        float z_st=-100;
        float z_end=150;

	TH1D *h1d_zend = new TH1D("h1d_zend", "reco+BQ", dz, z_st, z_end);
	TH1D *h1d_zend_XY = new TH1D("h1d_zend_XY", "reco+BQ+XY", dz, z_st, z_end);
	TH1D *h1d_zend_bmrw = new TH1D("h1d_zend_bmrw", "reco_bmrw+BQ", dz, z_st, z_end);
	TH1D *h1d_zend_bmrw_XY = new TH1D("h1d_zend_bmrw_XY", "reco_bmrw+BQ+XY", dz, z_st, z_end);
	//Weighted Gaussians & Weighting function ----------------------------------------------------------------------------------------------------------------------------------------------//


	//KEs ---------------------------------------------------------------------------------------------------------------------------------//
	int nx_trklen=150;
	double xmin_trklen=0;
	double xmax_trklen=150;
	int ny_edept=500;
	double ymin_edept=-200;
	double ymax_edept=800;

	TH1D *h1d_ke0=new TH1D("h1d_ke0","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_ke0_stop=new TH1D("h1d_ke0_stop","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_ke0_el=new TH1D("h1d_ke0_el","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_ke0_inel=new TH1D("h1d_ke0_inel","",ny_edept,ymin_edept,ymax_edept);
	h1d_ke0->Sumw2();
	h1d_ke0_stop->Sumw2();
	h1d_ke0_el->Sumw2();
	h1d_ke0_inel->Sumw2();

	TH1D *h1d_kebeam=new TH1D("h1d_kebeam","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_kebeam_stop=new TH1D("h1d_kebeam_stop","",ny_edept,ymin_edept,ymax_edept);
	h1d_kebeam->Sumw2();
	h1d_kebeam_stop->Sumw2();

	TH1D *h1d_keff=new TH1D("h1d_keff","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_keff_stop=new TH1D("h1d_keff_stop","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_keff_el=new TH1D("h1d_keff_el","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_keff_inel=new TH1D("h1d_keff_inel","",ny_edept,ymin_edept,ymax_edept);
	h1d_keff->Sumw2();
	h1d_keff_stop->Sumw2();
	h1d_keff_el->Sumw2();
	h1d_keff_inel->Sumw2();

	//TH1D *h1d_kehy=new TH1D("h1d_kehy","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_kehy_stop=new TH1D("h1d_kehy_stop","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_kehy_inel=new TH1D("h1d_kehy_inel","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_kehy_el=new TH1D("h1d_kehy_el","",ny_edept,ymin_edept,ymax_edept);
	h1d_kehy_stop->Sumw2();
	h1d_kehy_inel->Sumw2();
	h1d_kehy_el->Sumw2();

	TH1D *h1d_dke_keffbeam_keff_el=new TH1D("h1d_dke_keffbeam_keff_el","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_dke_keffbeam_keff_inel=new TH1D("h1d_dke_keffbeam_keff_inel","",ny_edept,ymin_edept,ymax_edept);

	TH1D *h1d_keffbeam_stop=new TH1D("h1d_keffbeam_stop","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_keffbeam_el_noxy=new TH1D("h1d_keffbeam_el_noxy","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_keffbeam_el=new TH1D("h1d_keffbeam_el","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_keffbeam_el_inel=new TH1D("h1d_keffbeam_el_inel","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_keffbeam_el_el=new TH1D("h1d_keffbeam_el_el","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_keffbeam_el_midcosmic=new TH1D("h1d_keffbeam_el_midcosmic","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_keffbeam_el_midpi=new TH1D("h1d_keffbeam_el_midpi","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_keffbeam_el_midp=new TH1D("h1d_keffbeam_el_midp","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_keffbeam_el_midmu=new TH1D("h1d_keffbeam_el_midmu","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_keffbeam_el_mideg=new TH1D("h1d_keffbeam_el_mideg","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_keffbeam_el_midother=new TH1D("h1d_keffbeam_el_midother","",ny_edept,ymin_edept,ymax_edept); 
	h1d_keffbeam_stop->Sumw2();
	h1d_keffbeam_el->Sumw2();
	h1d_keffbeam_el_inel->Sumw2();
	h1d_keffbeam_el_el->Sumw2();
	h1d_keffbeam_el_midcosmic->Sumw2();
	h1d_keffbeam_el_midpi->Sumw2();
	h1d_keffbeam_el_midp->Sumw2();
	h1d_keffbeam_el_midmu->Sumw2();
	h1d_keffbeam_el_mideg->Sumw2();
	h1d_keffbeam_el_midother->Sumw2();

	TH1D *h1d_keffbeam_inel=new TH1D("h1d_keffbeam_inel","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_keffbeam_inel_inel=new TH1D("h1d_keffbeam_inel_inel","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_keffbeam_inel_el=new TH1D("h1d_keffbeam_inel_el","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_keffbeam_inel_midcosmic=new TH1D("h1d_keffbeam_inel_midcosmic","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_keffbeam_inel_midpi=new TH1D("h1d_keffbeam_inel_midpi","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_keffbeam_inel_midp=new TH1D("h1d_keffbeam_inel_midp","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_keffbeam_inel_midmu=new TH1D("h1d_keffbeam_inel_midmu","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_keffbeam_inel_mideg=new TH1D("h1d_keffbeam_inel_mideg","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_keffbeam_inel_midother=new TH1D("h1d_keffbeam_inel_midother","",ny_edept,ymin_edept,ymax_edept); 
	h1d_keffbeam_inel->Sumw2();
	h1d_keffbeam_inel_inel->Sumw2();
	h1d_keffbeam_inel_el->Sumw2();
	h1d_keffbeam_inel_midcosmic->Sumw2();
	h1d_keffbeam_inel_midpi->Sumw2();
	h1d_keffbeam_inel_midp->Sumw2();
	h1d_keffbeam_inel_midmu->Sumw2();
	h1d_keffbeam_inel_mideg->Sumw2();
	h1d_keffbeam_inel_midother->Sumw2();

	TH1D *h1d_kerange_stop=new TH1D("h1d_kerange_stop","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kecalo_stop=new TH1D("h1d_kecalo_stop","", ny_edept, ymin_edept, ymax_edept);
	h1d_kerange_stop->Sumw2();
	h1d_kecalo_stop->Sumw2();

	TH1D *h1d_kend_calo_stop=new TH1D("h1d_kend_calo_stop","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kend_calo_el=new TH1D("h1d_kend_calo_el","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kend_calo_el_inel=new TH1D("h1d_kend_calo_el_inel","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_calo_el_el=new TH1D("h1d_kend_calo_el_el","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_calo_el_midcosmic=new TH1D("h1d_kend_calo_el_midcosmic","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_calo_el_midpi=new TH1D("h1d_kend_calo_el_midpi","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_calo_el_midp=new TH1D("h1d_kend_calo_el_midp","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_calo_el_midmu=new TH1D("h1d_kend_calo_el_midmu","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_calo_el_mideg=new TH1D("h1d_kend_calo_el_mideg","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_calo_el_midother=new TH1D("h1d_kend_calo_el_midother","",ny_edept,ymin_edept,ymax_edept); 
	h1d_kend_calo_stop->Sumw2();
	h1d_kend_calo_el->Sumw2();
	h1d_kend_calo_el_inel->Sumw2();
	h1d_kend_calo_el_el->Sumw2();
	h1d_kend_calo_el_midcosmic->Sumw2();
	h1d_kend_calo_el_midpi->Sumw2();
	h1d_kend_calo_el_midp->Sumw2();
	h1d_kend_calo_el_midmu->Sumw2();
	h1d_kend_calo_el_mideg->Sumw2();
	h1d_kend_calo_el_midother->Sumw2();

	TH1D *h1d_kend_calo_inel=new TH1D("h1d_kend_calo_inel","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kend_calo_inel_inel=new TH1D("h1d_kend_calo_inel_inel","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_calo_inel_el=new TH1D("h1d_kend_calo_inel_el","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_calo_inel_midcosmic=new TH1D("h1d_kend_calo_inel_midcosmic","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_calo_inel_midpi=new TH1D("h1d_kend_calo_inel_midpi","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_calo_inel_midp=new TH1D("h1d_kend_calo_inel_midp","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_calo_inel_midmu=new TH1D("h1d_kend_calo_inel_midmu","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_calo_inel_mideg=new TH1D("h1d_kend_calo_inel_mideg","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_calo_inel_midother=new TH1D("h1d_kend_calo_inel_midother","",ny_edept,ymin_edept,ymax_edept); 
	h1d_kend_calo_inel->Sumw2();
	h1d_kend_calo_inel_inel->Sumw2();
	h1d_kend_calo_inel_el->Sumw2();
	h1d_kend_calo_inel_midcosmic->Sumw2();
	h1d_kend_calo_inel_midpi->Sumw2();
	h1d_kend_calo_inel_midp->Sumw2();
	h1d_kend_calo_inel_midmu->Sumw2();
	h1d_kend_calo_inel_mideg->Sumw2();
	h1d_kend_calo_inel_midother->Sumw2();

	TH1D *h1d_kend_bb_stop=new TH1D("h1d_kend_bb_stop","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kend_bb_el=new TH1D("h1d_kend_bb_el","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kend_bb_el_inel=new TH1D("h1d_kend_bb_el_inel","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_bb_el_el=new TH1D("h1d_kend_bb_el_el","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_bb_el_midcosmic=new TH1D("h1d_kend_bb_el_midcosmic","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_bb_el_midpi=new TH1D("h1d_kend_bb_el_midpi","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_bb_el_midp=new TH1D("h1d_kend_bb_el_midp","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_bb_el_midmu=new TH1D("h1d_kend_bb_el_midmu","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_bb_el_mideg=new TH1D("h1d_kend_bb_el_mideg","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_bb_el_midother=new TH1D("h1d_kend_bb_el_midother","",ny_edept,ymin_edept,ymax_edept); 
	h1d_kend_bb_stop->Sumw2();
	h1d_kend_bb_el->Sumw2();
	h1d_kend_bb_el_inel->Sumw2();
	h1d_kend_bb_el_el->Sumw2();
	h1d_kend_bb_el_midcosmic->Sumw2();
	h1d_kend_bb_el_midpi->Sumw2();
	h1d_kend_bb_el_midp->Sumw2();
	h1d_kend_bb_el_midmu->Sumw2();
	h1d_kend_bb_el_mideg->Sumw2();
	h1d_kend_bb_el_midother->Sumw2();

	TH1D *h1d_kend_bb_inel=new TH1D("h1d_kend_bb_inel","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kend_bb_inel_inel=new TH1D("h1d_kend_bb_inel_inel","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_bb_inel_el=new TH1D("h1d_kend_bb_inel_el","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_bb_inel_midcosmic=new TH1D("h1d_kend_bb_inel_midcosmic","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_bb_inel_midpi=new TH1D("h1d_kend_bb_inel_midpi","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_bb_inel_midp=new TH1D("h1d_kend_bb_inel_midp","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_bb_inel_midmu=new TH1D("h1d_kend_bb_inel_midmu","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_bb_inel_mideg=new TH1D("h1d_kend_bb_inel_mideg","",ny_edept,ymin_edept,ymax_edept); 
	TH1D *h1d_kend_bb_inel_midother=new TH1D("h1d_kend_bb_inel_midother","",ny_edept,ymin_edept,ymax_edept); 
	h1d_kend_bb_inel->Sumw2();
	h1d_kend_bb_inel_inel->Sumw2();
	h1d_kend_bb_inel_el->Sumw2();
	h1d_kend_bb_inel_midcosmic->Sumw2();
	h1d_kend_bb_inel_midpi->Sumw2();
	h1d_kend_bb_inel_midp->Sumw2();
	h1d_kend_bb_inel_midmu->Sumw2();
	h1d_kend_bb_inel_mideg->Sumw2();
	h1d_kend_bb_inel_midother->Sumw2();

	TH1D *h1d_kend_true_stop=new TH1D("h1d_kend_true_stop","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kend_true_el=new TH1D("h1d_kend_true_el","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kend_true_inel=new TH1D("h1d_kend_true_inel","", ny_edept, ymin_edept, ymax_edept);
	h1d_kend_true_stop->Sumw2();
	h1d_kend_true_el->Sumw2();
	h1d_kend_true_inel->Sumw2();

	TH1D *h1d_kend_bbtrue_stop=new TH1D("h1d_kend_bbtrue_stop","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kend_bbtrue_el=new TH1D("h1d_kend_bbtrue_el","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kend_bbtrue_inel=new TH1D("h1d_kend_bbtrue_inel","", ny_edept, ymin_edept, ymax_edept);
	h1d_kend_bbtrue_stop->Sumw2();
	h1d_kend_bbtrue_el->Sumw2();
	h1d_kend_bbtrue_inel->Sumw2();

	TH1D *h1d_kend_bbtruelength_stop=new TH1D("h1d_kend_bbtruelength_stop","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kend_bbtruelength_el=new TH1D("h1d_kend_bbtruelength_el","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kend_bbtruelength_inel=new TH1D("h1d_kend_bbtruelength_inel","", ny_edept, ymin_edept, ymax_edept);
	h1d_kend_bbtruelength_stop->Sumw2();
	h1d_kend_bbtruelength_el->Sumw2();
	h1d_kend_bbtruelength_inel->Sumw2();

	TH1D *h1d_kend_bbcorr_el=new TH1D("h1d_kend_bbcorr_el","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kend_bbcorr_inel=new TH1D("h1d_kend_bbcorr_inel","", ny_edept, ymin_edept, ymax_edept);

	//dedx_rr ---------------------------------------------------------------------------//
	TH2D *h2d_rr_dedx_recoSTOP=new TH2D("h2d_rr_dedx_recoSTOP","",240,0,120,90,0,30);

	//chi2_pid ------------------------------------------------------------//
	TH1D *chi2pid_recostop=new TH1D("chi2pid_recostop","",500,0,100);
	TH1D *chi2pid_truestop=new TH1D("chi2pid_truestop","",500,0,100);
	TH1D *chi2pid_trueel=new TH1D("chi2pid_trueel","",500,0,100);
	TH1D *chi2pid_trueinel=new TH1D("chi2pid_trueinel","",500,0,100);
	TH1D *chi2pid_recoinel=new TH1D("chi2pid_recoinel","",500,0,100);

	//xy-dist
	TH2D *h2d_xy_noSCE=new TH2D("h2d_xy_noSCE","", 70,-60,10,60,390,450); //nosce
	TH2D *h2d_xy_SCE=new TH2D("h2d_xy_SCE","", 70,-60,10,60,390,450); //after sce
	
	//z-st
	TH1D *h1d_zst_noSCE=new TH1D("h1d_zst_noSCE","",110,-10,100);
	TH1D *h1d_zst_SCE=new TH1D("h1d_zst_SCE","",110,-10,100);

	//Ratio plot -----------------------------------------------------------------------------------------------------------------//
	TH2D *h2d_trklen_ratio_KEbbfit_KEffbeam_stop=new TH2D("h2d_trklen_ratio_KEbbfit_KEffbeam_stop","", 14000,0,140,4000,-20,20);
	TH2D *h2d_trklen_ratio_KEbbfit_KEffbeam_el=new TH2D("h2d_trklen_ratio_KEbbfit_KEffbeam_el","", 14000,0,140,4000,-20,20);
	TH2D *h2d_trklen_ratio_KEbbfit_KEffbeam_inel=new TH2D("h2d_trklen_ratio_KEbbfit_KEffbeam_inel","", 14000,0,140,4000,-20,20);
	
	TH2D *h2d_trklen_ratio_KEbbtruth_KEffbeam_stop=new TH2D("h2d_trklen_ratio_KEbbtruth_KEffbeam_stop","", 14000,0,140,4000,-20,20);
	TH2D *h2d_trklen_ratio_KEbbtruth_KEffbeam_el=new TH2D("h2d_trklen_ratio_KEbbtruth_KEffbeam_el","", 14000,0,140,4000,-20,20);
	TH2D *h2d_trklen_ratio_KEbbtruth_KEffbeam_inel=new TH2D("h2d_trklen_ratio_KEbbtruth_KEffbeam_inel","", 14000,0,140,4000,-20,20);

	//ff	
	TH2D *h2d_keffbeam_keff_stop=new TH2D("h2d_keffbeam_keff_stop","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);
	TH2D *h2d_keffbeam_keff_el=new TH2D("h2d_keffbeam_keff_el","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);
	TH2D *h2d_keffbeam_keff_inel=new TH2D("h2d_keffbeam_keff_inel","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);
	
	TH2D *h2d_kehy_keff_stop=new TH2D("h2d_kehy_keff_stop","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);
	TH2D *h2d_kehy_keff_el=new TH2D("h2d_kehy_keff_el","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);
	TH2D *h2d_kehy_keff_inel=new TH2D("h2d_kehy_keff_inel","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);

	//end
	TH2D *h2d_kebb_keendtruth_stop=new TH2D("h2d_kebb_keendtruth_stop","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);
	TH2D *h2d_kebb_keendtruth_el=new TH2D("h2d_kebb_keendtruth_el","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);
	TH2D *h2d_kebb_keendtruth_inel=new TH2D("h2d_kebb_keendtruth_inel","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);

	TH2D *h2d_kebbtruelength_keendtruth_stop=new TH2D("h2d_kebbtruelength_keendtruth_stop","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);
	TH2D *h2d_kebbtruelength_keendtruth_el=new TH2D("h2d_kebbtruelength_keendtruth_el","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);
	TH2D *h2d_kebbtruelength_keendtruth_inel=new TH2D("h2d_kebbtruelength_keendtruth_inel","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);
	TH2D *h2d_kebbtruelength_keendtruth_trueel=new TH2D("h2d_kebbtruelength_keendtruth_trueel","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);


	
	TH2D *h2d_kecalo_keendtruth_stop=new TH2D("h2d_kecalo_keendtruth_stop","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);
	TH2D *h2d_kecalo_keendtruth_el=new TH2D("h2d_kecalo_keendtruth_el","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);
	TH2D *h2d_kecalo_keendtruth_inel=new TH2D("h2d_kecalo_keendtruth_inel","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);

	TH2D *h2d_kebb_kebbtruth_stop=new TH2D("h2d_kebb_kebbtruth_stop","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);
	TH2D *h2d_kebb_kebbtruth_el=new TH2D("h2d_kebb_kebbtruth_el","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);
	TH2D *h2d_kebb_kebbtruth_inel=new TH2D("h2d_kebb_kebbtruth_inel","", ny_edept, ymin_edept, ymax_edept, ny_edept, ymin_edept, ymax_edept);

	TH1D *h1d_ratio_rangetrue_rangereco_stop=new TH1D("h1d_ratio_rangetrue_rangereco_stop","",4000,-20,20);
	TH1D *h1d_ratio_rangetrue_rangereco_el=new TH1D("h1d_ratio_rangetrue_rangereco_el","",4000,-20,20);
	TH1D *h1d_ratio_rangetrue_rangereco_inel=new TH1D("h1d_ratio_rangetrue_rangereco_inel","",4000,-20,20);

	TH2D *h2d_recorange_truerange_stop=new TH2D("h2d_recorange_truerange_stop","", 1400, 0, 140, 1400, 0, 140);
	TH2D *h2d_recorange_truerange_el=new TH2D("h2d_recorange_truerange_el","", 1400, 0, 140, 1400, 0, 140);
	TH2D *h2d_recorange_truerange_inel=new TH2D("h2d_recorange_truerange_inel","", 1400, 0, 140, 1400, 0, 140);

	TH1D *h1d_ratio_rangehy_rangereco_stop=new TH1D("h1d_ratio_rangehy_rangereco_stop","",4000,-20,20);
	TH1D *h1d_ratio_rangehy_rangereco_el=new TH1D("h1d_ratio_rangehy_rangereco_el","",4000,-20,20);

	TH2D *h2d_recorange_rangehy_el=new TH2D("h2d_recorange_rangehy_el","", 1400, 0, 140, 1400, 0, 140);
	TH2D *h2d_recorange_rangehy_stop=new TH2D("h2d_recorange_rangehy_stop","", 1400, 0, 140, 1400, 0, 140);


	//Efficiency vs trklen ------------------------------------------------------------------------------------------//
	//2D histograms
	int nx__trklen=1400;
	double xmin__trklen=0;
	double xmax__trklen=140;
	int ny_eff=2000;
	double ymin_eff=-20;
	double ymax_eff=20;
	TH2D *h2d_trklen_eff_KEconst_KEff_all=new TH2D("h2d_trklen_eff_KEconst_KEff_all","", nx__trklen, xmin__trklen, xmax__trklen, ny_eff, ymin_eff, ymax_eff);
	TH2D *h2d_trklen_eff_KEconst_KEff_stop=new TH2D("h2d_trklen_eff_KEconst_KEff_stop","", nx__trklen, xmin__trklen, xmax__trklen, ny_eff, ymin_eff, ymax_eff);
	TH2D *h2d_trklen_eff_KEconst_KEff_el=new TH2D("h2d_trklen_eff_KEconst_KEff_el","", nx__trklen, xmin__trklen, xmax__trklen, ny_eff, ymin_eff, ymax_eff);
	TH2D *h2d_trklen_eff_KEconst_KEff_inel=new TH2D("h2d_trklen_eff_KEconst_KEff_inel","", nx__trklen, xmin__trklen, xmax__trklen, ny_eff, ymin_eff, ymax_eff);

	TH2D *h2d_trklen_eff_KEhy_KEff_all=new TH2D("h2d_trklen_eff_KEhy_KEff_all","", nx__trklen, xmin__trklen, xmax__trklen, ny_eff, ymin_eff, ymax_eff);
	TH2D *h2d_trklen_eff_KEhy_KEff_stop=new TH2D("h2d_trklen_eff_KEhy_KEff_stop","", nx__trklen, xmin__trklen, xmax__trklen, ny_eff, ymin_eff, ymax_eff);
	TH2D *h2d_trklen_eff_KEhy_KEff_el=new TH2D("h2d_trklen_eff_KEhy_KEff_el","", nx__trklen, xmin__trklen, xmax__trklen, ny_eff, ymin_eff, ymax_eff);
	TH2D *h2d_trklen_eff_KEhy_KEff_inel=new TH2D("h2d_trklen_eff_KEhy_KEff_inel","", nx__trklen, xmin__trklen, xmax__trklen, ny_eff, ymin_eff, ymax_eff);

	TH2D *h2d_trklen_eff_KEhy_KEconst_all=new TH2D("h2d_trklen_eff_KEhy_KEconst_all","", nx__trklen, xmin__trklen, xmax__trklen, ny_eff, ymin_eff, ymax_eff);
	TH2D *h2d_trklen_eff_KEhy_KEconst_stop=new TH2D("h2d_trklen_eff_KEhy_KEconst_stop","", nx__trklen, xmin__trklen, xmax__trklen, ny_eff, ymin_eff, ymax_eff);
	TH2D *h2d_trklen_eff_KEhy_KEconst_el=new TH2D("h2d_trklen_eff_KEhy_KEconst_el","", nx__trklen, xmin__trklen, xmax__trklen, ny_eff, ymin_eff, ymax_eff);
	TH2D *h2d_trklen_eff_KEhy_KEconst_inel=new TH2D("h2d_trklen_eff_KEhy_KEconst_inel","", nx__trklen, xmin__trklen, xmax__trklen, ny_eff, ymin_eff, ymax_eff);

	//1D histograms
	//TH1D *h1d_KEconst_all=new TH1D("h1d_KEconst_all","",nx_trklen, xmin_trklen, xmax_trklen);
	//TH1D *h1d_KEconst_el=new TH1D("h1d_KEconst_el","",nx_trklen, xmin_trklen, xmax_trklen);
	//TH1D *h1d_KEconst_inel=new TH1D("h1d_KEconst_inel","",nx_trklen, xmin_trklen, xmax_trklen);
	
	//TH1D *h1d_KEhy_all=new TH1D("h1d_KEhy_all","",nx_trklen, xmin_trklen, xmax_trklen);
	//TH1D *h1d_KEhy_el=new TH1D("h1d_KEhy_el","",nx_trklen, xmin_trklen, xmax_trklen);
	//TH1D *h1d_KEhy_inel=new TH1D("h1d_KEhy_inel","",nx_trklen, xmin_trklen, xmax_trklen);
	
	//TH1D *h1d_KEff_all=new TH1D("h1d_KEff_all","",nx_trklen, xmin_trklen, xmax_trklen);
	//TH1D *h1d_KEff_el=new TH1D("h1d_KEff_el","",nx_trklen, xmin_trklen, xmax_trklen);
	//TH1D *h1d_KEff_inel=new TH1D("h1d_KEff_inel","",nx_trklen, xmin_trklen, xmax_trklen);


	//Basic configure ------//
	BetheBloch BB;
	BB.SetPdgCode(pdg);
	//----------------------//	

	for (Long64_t jentry=0; jentry<nentries;jentry++) { //main entry loop
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		isTestSample = true;
		if (ientry%2 == 0) isTestSample = false; //Divide MC sample by 2 parts: test+ufold
                if (jentry%1000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;

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
		} 
		if (key_reach_tpc!=-99) { is_beam_at_ff=true; }
		//cout<<"key_reach_tpc:"<<key_reach_tpc<<endl;	
		//cout<<"is_beam_at_ff:"<<is_beam_at_ff<<endl;

		//Get true trklen ---------------------------------------------------------------------------------------//
		double range_true=-999;
		int key_st = 0;
		double tmp_z = 9999;
		vector<double> true_trklen_accum;
		//cout<<"beamtrk_z->size():"<<beamtrk_z->size()<<endl;
		true_trklen_accum.reserve(beamtrk_z->size()); // initialize true_trklen_accum
		//cout<<"ck0"<<endl;
		if (is_beam_at_ff) {
			for (int iz=0; iz<(int)beamtrk_z->size(); iz++) {
				if (abs(beamtrk_z->at(iz)) < tmp_z){
					tmp_z = abs(beamtrk_z->at(iz));
					key_st = iz; // find the point where the beam enters the TPC (find the smallest abs(Z))
				}
				true_trklen_accum[iz] = 0.; // initialize true_trklen_accum
			}
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

			//bool ok = fitter.FitFCN();
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

		//[2] Range compensation ---------------------------------------------------------------------
		double range_true_patch=0;
		if (is_beam_at_ff) {
			//calculate distance 1st hit and pojected point at TPC front face
			range_true_patch = sqrt( pow(beamtrk_x->at(key_reach_tpc)-xproj_beam, 2)+
					pow(beamtrk_y->at(key_reach_tpc)-yproj_beam, 2)+	
					pow(beamtrk_z->at(key_reach_tpc)-zproj_beam, 2) );
			//range_true_patch=0; //no fix on true len

			//true_trklen_accum
			for (int iz=key_reach_tpc+1; iz<(int)beamtrk_z->size(); iz++) {
				if (iz == key_reach_tpc+1) range_true = range_true_patch;
					range_true += sqrt( pow(beamtrk_x->at(iz)-beamtrk_x->at(iz-1), 2)+
						pow(beamtrk_y->at(iz)-beamtrk_y->at(iz-1), 2)+	
						pow(beamtrk_z->at(iz)-beamtrk_z->at(iz-1), 2) );						    	
					true_trklen_accum[iz] = range_true;
			}						    	
		}













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
		double bx_spec=beamPosx_spec->at(0);
		double by_spec=beamPosy_spec->at(0);
		
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

		//double range_reco=-99; if (!primtrk_range->empty()) range_reco=primtrk_range->at(0); //reco primary trklen
		double csda_val_spec=csda_range_vs_mom_sm->Eval(mom_beam_spec);

		if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) IsRecoStop=true;
		//if ((range_reco/csda_val_spec)<min_norm_trklen_csda) IsRecoInEL=true;

		//kinetic energies
		double ke_beam_MeV=1000.*p2ke(mom_beam); //ke_beam [MeV]
		double ke_beam_spec=p2ke(mom_beam_spec); //ke_beam_spec [GeV]
		double ke_beam_spec_MeV=1000.*ke_beam_spec; //ke_beam_spec [MeV]
		double ke_trklen=ke_vs_csda_range_sm->Eval(range_reco); //[unit: GeV]
		double ke_trklen_MeV=1000.*ke_trklen; //[unit: MeV]

		//double p_trklen=ke2p(ke_trklen);
		//double ke_simide=0;
		//for (int hk=0; hk<(int)primtrk_true_edept->size(); ++hk) { //loop over simIDE points
		//ke_simide+=primtrk_true_edept->at(hk);
		//} //loop over simIDE points

		double ke_calo_MeV=0;
		double pid=-99;
		//if (IsPandoraSlice&&IsBQ&&IsCaloSize) { //if calo size not empty
		vector<double> trkdedx; 
		vector<double> trkres;
		if (IsCaloSize) { //if calo size not empty
			//no-sce
			h2d_xy_noSCE->Fill(reco_stx_noSCE, reco_sty_noSCE);
			h1d_zst_noSCE->Fill(reco_stz_noSCE);
			//after sce
			h2d_xy_SCE->Fill(reco_stx, reco_sty);
			h1d_zst_SCE->Fill(reco_stz);
			
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

				if (IsRecoStop) {
					h2d_rr_dedx_recoSTOP->Fill(resrange_reco, cali_dedx);
				}

			} //loop over reco hits of a given track
			pid=chi2pid(trkdedx,trkres); //pid using stopping proton hypothesis

			if (IsRecoStop) chi2pid_recostop->Fill(pid);
			if (IsRecoInEL) chi2pid_recoinel->Fill(pid);
			//if (IsPureMCS) chi2pid_truestop->Fill(pid);
			if (kel) chi2pid_trueel->Fill(pid);
			if (kinel) chi2pid_trueinel->Fill(pid); 

		} //if calo size not empty

		if ((range_reco/csda_val_spec)<min_norm_trklen_csda) { //inel region
			if (pid>pid_1) IsRecoInEL=true;
			if (pid<=pid_1) IsRecoEL=true;
		} //inel region
		if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) { //stopping p region
			if (pid>pid_2) IsRecoInEL=true;
			if (pid<=pid_2) IsRecoEL=true;
		}
		
		//Const E-loss ------------------------------------------------------------------------------------------//
		//double ke_ffbeam_MeV=ke_beam_spec_MeV-const_eloss_mc; //const E-loss (our assumption of KE at TPC FF)
		double slope_kff_keffbeamR=0.859605;
		double ke_ffbeam_MeV=(ke_beam_spec_MeV-Eloss_mc_hy_stop); //const E-loss with correction
		//double ke_ffbeam_MeV=(ke_beam_spec_MeV-Eloss_mc_hy_stop)*R_fit_hy; //const E-loss with correction
		//double ke_ffbeam_MeV=(ke_beam_spec_MeV-Eloss_mc_hy_stop)*R_fit_hy*slope_kff_keffbeamR; //const E-loss with correction

		//hypothetical length -------------------------------------------------------------------------------------//
		double fitted_length=-1; 
		double tmp_fitted_length=BB.Fit_dEdx_Residual_Length(trkdedx, trkres, pdg, false);
		//double tmp_fitted_length=BB.Fit_Proton_Residual_Length_Likelihood(trkdedx, trkres, pdg, false);
		if (tmp_fitted_length>0) fitted_length=tmp_fitted_length;
		double fitted_KE=-50; 
		if (fitted_length>0) fitted_KE=BB.KEFromRangeSpline(fitted_length);
	
		//ke at end point ---------------------------------------------------------------------//
		//double kebb=-50; if (fitted_KE>0) kebb=BB.KEAtLength(ke_ffbeam_MeV, range_reco);

		double ratio_range_reco_stop=1.0071824773690985;

		//double kebb=-50; kebb=BB.KEAtLength(ke_ffbeam_MeV, range_reco);
		double kebb=-1000; kebb=BB.KEAtLength(ke_ffbeam_MeV, range_reco);
		double kebb_corr=-1000; kebb_corr=BB.KEAtLength(ke_ffbeam_MeV, range_reco*ratio_range_reco_stop);

		double kebb_fit=-1000; kebb_fit=BB.KEAtLength(fitted_KE, range_reco);
		//double kebb_truth=-50; kebb_truth=BB.KEAtLength(ke_ff, range_reco);
		double kebb_truth=-1000; kebb_truth=BB.KEAtLength(ke_ff, range_true);

		double r_keconst_keff=kebb/kebb_truth;
		double r_kefit_keff=kebb_fit/kebb_truth;
		double r_kefit_keconst=kebb_fit/kebb;

		double kebb_truerange=-1000; kebb_truerange=BB.KEAtLength(ke_ffbeam_MeV, range_true);

		double kecalo=-1000; kecalo=ke_ffbeam_MeV-ke_calo_MeV;
		double kend=-1000; kend=1000.*(beamtrk_Eng->at(-2+beamtrk_Eng->size()));

		//range_ratio -------------------------------------------------//
		double r_range=range_true/range_reco;
		double r_hy_range=fitted_length/range_reco;


		//if (IsBeamMom&&IsPandoraSlice&&IsBQ&&IsCaloSize) { //basic cuts
			//if (IsRecoEL) { //reco el
			//}
		//}

		//if (IsPandoraSlice&&IsBQ&&IsCaloSize) { //basic cuts
		//if (IsBeamXY&&IsPandoraSlice&&IsBQ&&IsCaloSize) { //basic cuts
		//if (IsBeamMom&&IsBeamXY&&IsPandoraSlice&&IsBQ&&IsCaloSize) { //basic cuts
		if (IsBeamMom&&IsPandoraSlice&&IsBQ&&IsCaloSize) { //basic cuts
			double mom_rw_minchi2=1.;
			//if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) mom_rw_minchi2=gng[index_minchi2]->Eval(mom_beam_spec*1000.); //old bmrw
			if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) mom_rw_minchi2=kerw->Eval(ke_ffbeam_MeV); //new bmrw
			//std::cout<<"ke_beam_spec_MeV:"<<ke_beam_spec_MeV<<" ke_ffbeam_MeV:"<<ke_ffbeam_MeV<<" mom_rw_minchi2:"<<mom_rw_minchi2<<std::endl;

			h1d_keffbeam_el_noxy->Fill(ke_ffbeam_MeV, mom_rw_minchi2);

			if (IsBeamXY) { //beam xy

			h1d_ke0->Fill(ke_beam_MeV, mom_rw_minchi2);
			h1d_p0->Fill(mom_beam_MeV, mom_rw_minchi2);
			h1d_kebeam->Fill(ke_beam_spec_MeV, mom_rw_minchi2);
			h1d_pbeam->Fill(1000.*mom_beam_spec, mom_rw_minchi2);

			h1d_keff->Fill(ke_ff, mom_rw_minchi2);
			h1d_pff->Fill(1000.*ke2p(ke_ff/1000.), mom_rw_minchi2);

			//h1d_kehy->Fill(fitted_KE, mom_rw_minchi2);

			h1d_trklen->Fill(range_reco, mom_rw_minchi2);	
			h1d_zend->Fill(reco_endz, mom_rw_minchi2);

			h1d_trklen_bmrw->Fill(range_reco, mom_rw_minchi2); 
			h1d_zend_bmrw->Fill(reco_endz, mom_rw_minchi2);	

			h2d_trklen_eff_KEconst_KEff_all->Fill(range_reco, r_keconst_keff);
			h2d_trklen_eff_KEhy_KEff_all->Fill(range_reco, r_kefit_keff);
			h2d_trklen_eff_KEhy_KEconst_all->Fill(range_reco, r_kefit_keconst);

			//if (IsXY) { //xy
				//h1d_trklen_XY->Fill(range_reco);
				//h1d_trklen_bmrw_XY->Fill(range_reco, mom_rw_minchi2);

				//h1d_zend_XY->Fill(reco_endz);
				//h1d_zend_bmrw_XY->Fill(reco_endz, mom_rw_minchi2);
			//} //xy

			if (IsRecoEL) { //reco el
				h1d_ke0_el->Fill(ke_beam_MeV, mom_rw_minchi2);
				h1d_kehy_el->Fill(fitted_KE, mom_rw_minchi2);
				h1d_keffbeam_el->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
				h1d_kend_calo_el->Fill(kecalo, mom_rw_minchi2);
				h1d_kend_bb_el->Fill(kebb, mom_rw_minchi2);
				h1d_kend_bbcorr_el->Fill(kebb_corr, mom_rw_minchi2);
				h1d_kend_true_el->Fill(kend, mom_rw_minchi2);
				h1d_kend_bbtrue_el->Fill(kebb_truth, mom_rw_minchi2);
				h1d_kend_bbtruelength_el->Fill(kebb_truerange, mom_rw_minchi2);

				h1d_keff_el->Fill(ke_ff, mom_rw_minchi2);
				h2d_trklen_ratio_KEbbfit_KEffbeam_el->Fill(range_reco, (kebb_fit/kebb));
				h2d_trklen_ratio_KEbbtruth_KEffbeam_el->Fill(range_reco, (kebb_truth/kebb));
				h2d_keffbeam_keff_el->Fill(ke_ffbeam_MeV, ke_ff);
				h2d_kehy_keff_el->Fill(fitted_KE, ke_ff);
				h2d_kebb_keendtruth_el->Fill(kebb, kend);
				h2d_kebbtruelength_keendtruth_el->Fill(kebb_truerange, kend);

				h2d_kecalo_keendtruth_el->Fill(kecalo, kend);
				h2d_kebb_kebbtruth_el->Fill(kebb, kebb_truth);

				h1d_ratio_rangetrue_rangereco_el->Fill(r_range);
				h1d_ratio_rangehy_rangereco_el->Fill(r_hy_range);
				h2d_recorange_truerange_el->Fill(range_reco,range_true);
				h2d_recorange_rangehy_el->Fill(range_reco, fitted_length);

				h1d_dke_keffbeam_keff_el->Fill(ke_ffbeam_MeV-ke_ff);

				h2d_trklen_eff_KEconst_KEff_el->Fill(range_reco, r_keconst_keff);
				h2d_trklen_eff_KEhy_KEff_el->Fill(range_reco, r_kefit_keff);
				h2d_trklen_eff_KEhy_KEconst_el->Fill(range_reco, r_kefit_keconst);

				if (kinel) { //inel
					h1d_keffbeam_el_inel->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
					h1d_kend_calo_el_inel->Fill(kecalo, mom_rw_minchi2);
					h1d_kend_bb_el_inel->Fill(kebb, mom_rw_minchi2);
				} //inel
				if (kel) { 
					h1d_keffbeam_el_el->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
					h1d_kend_calo_el_el->Fill(kecalo, mom_rw_minchi2);
					h1d_kend_bb_el_el->Fill(kebb, mom_rw_minchi2);
					h2d_kebbtruelength_keendtruth_trueel->Fill(kebb_truerange, kend);
				}
				if (kMIDcosmic) {
					h1d_keffbeam_el_midcosmic->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
					h1d_kend_calo_el_midcosmic->Fill(kecalo, mom_rw_minchi2); 
					h1d_kend_bb_el_midcosmic->Fill(kebb, mom_rw_minchi2);
				}
				if (kMIDpi) { 
					h1d_keffbeam_el_midpi->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
					h1d_kend_calo_el_midpi->Fill(kecalo, mom_rw_minchi2);
					h1d_kend_bb_el_midpi->Fill(kebb, mom_rw_minchi2);
				}
				if (kMIDp) {
					h1d_keffbeam_el_midp->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
					h1d_kend_calo_el_midp->Fill(kecalo, mom_rw_minchi2);
					h1d_kend_bb_el_midp->Fill(kebb, mom_rw_minchi2);
				}
				if (kMIDmu) { 
					h1d_keffbeam_el_midmu->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
					h1d_kend_calo_el_midmu->Fill(kecalo, mom_rw_minchi2);
					h1d_kend_bb_el_midmu->Fill(kebb, mom_rw_minchi2);
				}
				if (kMIDeg) { 
					h1d_keffbeam_el_mideg->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
					h1d_kend_calo_el_mideg->Fill(kecalo, mom_rw_minchi2);
					h1d_kend_bb_el_mideg->Fill(kebb, mom_rw_minchi2);
				}
				if (kMIDother) { 
					h1d_keffbeam_el_midother->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
					h1d_kend_calo_el_midother->Fill(kecalo, mom_rw_minchi2);
					h1d_kend_bb_el_midother->Fill(kebb, mom_rw_minchi2);
				}

			} //reco el			

			if (IsRecoStop) { //reco stop 
				h1d_ke0_stop->Fill(ke_beam_MeV, mom_rw_minchi2);           h1d_p0_stop->Fill(mom_beam_MeV, mom_rw_minchi2);
				h1d_kebeam_stop->Fill(ke_beam_spec_MeV, mom_rw_minchi2);   h1d_pbeam_stop->Fill(1000.*mom_beam_spec, mom_rw_minchi2);
				h1d_keff_stop->Fill(ke_ff, mom_rw_minchi2);      	   h1d_pff_stop->Fill(1000.*ke2p(ke_ff/1000.), mom_rw_minchi2);

				h1d_kehy_stop->Fill(fitted_KE, mom_rw_minchi2);
				h1d_keffbeam_stop->Fill(ke_ffbeam_MeV, mom_rw_minchi2);

				h1d_trklen_stop->Fill(range_reco, mom_rw_minchi2);
				h1d_kerange_stop->Fill(ke_trklen_MeV, mom_rw_minchi2);     h1d_prange_stop->Fill(1000.*ke2p(ke_trklen), mom_rw_minchi2);
				h1d_kecalo_stop->Fill(ke_calo_MeV, mom_rw_minchi2);	   h1d_pcalo_stop->Fill(1000.*ke2p(ke_calo_MeV/1000.), mom_rw_minchi2);	

				h1d_kend_calo_stop->Fill(kecalo, mom_rw_minchi2);
				h1d_kend_bb_stop->Fill(kebb, mom_rw_minchi2);
				h1d_kend_true_stop->Fill(kend, mom_rw_minchi2);
				h1d_kend_bbtrue_stop->Fill(kebb_truth, mom_rw_minchi2);
				h1d_kend_bbtruelength_stop->Fill(kebb_truerange, mom_rw_minchi2);

				h2d_trklen_ratio_KEbbfit_KEffbeam_stop->Fill(range_reco, (kebb_fit/kebb));
				h2d_trklen_ratio_KEbbtruth_KEffbeam_stop->Fill(range_reco, (kebb_truth/kebb));

				h2d_keffbeam_keff_stop->Fill(ke_ffbeam_MeV, ke_ff);
				h2d_kehy_keff_stop->Fill(fitted_KE, ke_ff);
				h2d_kebb_keendtruth_stop->Fill(kebb, kend);
				h2d_kebbtruelength_keendtruth_stop->Fill(kebb_truerange, kend);

				h2d_kecalo_keendtruth_stop->Fill(kecalo, kend);
				h2d_kebb_kebbtruth_stop->Fill(kebb, kebb_truth);

				h1d_ratio_rangetrue_rangereco_stop->Fill(r_range);
				h1d_ratio_rangehy_rangereco_stop->Fill(r_hy_range);
				h2d_recorange_truerange_stop->Fill(range_reco,range_true);
				h2d_recorange_rangehy_stop->Fill(range_reco, fitted_length);

				h2d_trklen_eff_KEconst_KEff_stop->Fill(range_reco, r_keconst_keff);
				h2d_trklen_eff_KEhy_KEff_stop->Fill(range_reco, r_kefit_keff);
				h2d_trklen_eff_KEhy_KEconst_stop->Fill(range_reco, r_kefit_keconst);
				//if (IsXY) { //xy-cut
					//h1d_trklen_stop_XY->Fill(range_reco);
					//double mom_rw=1.;
					//if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) { //beam-mom cut (within 3-sigma)
						//for (int ig = 0; ig < n_1d; ++ig) { //rw loop
							//mom_rw=gng[ig]->Eval(mom_beam_spec*1000.);
							//h1d_trklen_rw[ig]->Fill(range_reco,mom_rw); //beam-mom rw using stopping protons
						//} //rw loop
					//} //beam-mom cut (within 3-sigma)
					//else { //tail of beam
					//if ((mom_beam_spec*1000.)<mu_min||(mom_beam_spec*1000.)>mu_max) { //tail of the beam
						//for (int ig = 0; ig < n_1d; ++ig) { //rw loop
							//h1d_trklen_rw[ig]->Fill(range_reco); //beam-mom rw 
						//} //rw loop
					//} //tail of the beam	
				//} //xy-cut
			} //reco stop

			if (IsRecoInEL) { //reco inel	
				h1d_ke0_inel->Fill(ke_beam_MeV, mom_rw_minchi2);
				h1d_keff_inel->Fill(ke_ff, mom_rw_minchi2);
				h1d_kehy_inel->Fill(fitted_KE, mom_rw_minchi2);
				h1d_keffbeam_inel->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
				h1d_kend_calo_inel->Fill(kecalo, mom_rw_minchi2);
				h1d_kend_bb_inel->Fill(kebb, mom_rw_minchi2);
				h1d_kend_bbcorr_inel->Fill(kebb_corr, mom_rw_minchi2);
				h1d_kend_true_inel->Fill(kend, mom_rw_minchi2);
				h1d_kend_bbtrue_inel->Fill(kebb_truth, mom_rw_minchi2);
				h1d_kend_bbtruelength_inel->Fill(kebb_truerange, mom_rw_minchi2);

				h2d_trklen_ratio_KEbbfit_KEffbeam_inel->Fill(range_reco, (kebb_fit/kebb));
				h2d_trklen_ratio_KEbbtruth_KEffbeam_inel->Fill(range_reco, (kebb_truth/kebb));

				h2d_keffbeam_keff_inel->Fill(ke_ffbeam_MeV, ke_ff);
				h2d_kehy_keff_inel->Fill(fitted_KE, ke_ff);

				h2d_kebb_keendtruth_inel->Fill(kebb, kend);
				h2d_kebbtruelength_keendtruth_inel->Fill(kebb_truerange, kend);
				h2d_kecalo_keendtruth_inel->Fill(kecalo, kend);
				h2d_kebb_kebbtruth_inel->Fill(kebb, kebb_truth);

				h1d_ratio_rangetrue_rangereco_inel->Fill(r_range);
				h2d_recorange_truerange_inel->Fill(range_reco,range_true);

				h1d_dke_keffbeam_keff_inel->Fill(ke_ffbeam_MeV-ke_ff);

				h2d_trklen_eff_KEconst_KEff_inel->Fill(range_reco, r_keconst_keff);
				h2d_trklen_eff_KEhy_KEff_inel->Fill(range_reco, r_kefit_keff);
				h2d_trklen_eff_KEhy_KEconst_inel->Fill(range_reco, r_kefit_keconst);

				if (kinel) { //inel
					h1d_keffbeam_inel_inel->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
					h1d_kend_calo_inel_inel->Fill(kecalo, mom_rw_minchi2);
					h1d_kend_bb_inel_inel->Fill(kebb, mom_rw_minchi2);
				} //inel
				if (kel) { 
					h1d_keffbeam_inel_el->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
					h1d_kend_calo_inel_el->Fill(kecalo, mom_rw_minchi2);
					h1d_kend_bb_inel_el->Fill(kebb, mom_rw_minchi2);
				}
				if (kMIDcosmic) {
					h1d_keffbeam_inel_midcosmic->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
					h1d_kend_calo_inel_midcosmic->Fill(kecalo, mom_rw_minchi2);
					h1d_kend_bb_inel_midcosmic->Fill(kebb, mom_rw_minchi2);
				}
				if (kMIDpi) { 
					h1d_keffbeam_inel_midpi->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
					h1d_kend_calo_inel_midpi->Fill(kecalo, mom_rw_minchi2);
					h1d_kend_bb_inel_midpi->Fill(kebb, mom_rw_minchi2);
				}
				if (kMIDp) {
					h1d_keffbeam_inel_midp->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
					h1d_kend_calo_inel_midp->Fill(kecalo, mom_rw_minchi2);
					h1d_kend_bb_inel_midp->Fill(kebb, mom_rw_minchi2);
				}
				if (kMIDmu) { 
					h1d_keffbeam_inel_midmu->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
					h1d_kend_calo_inel_midmu->Fill(kecalo, mom_rw_minchi2);
					h1d_kend_bb_inel_midmu->Fill(kebb, mom_rw_minchi2);
				}
				if (kMIDeg) { 
					h1d_keffbeam_inel_mideg->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
					h1d_kend_calo_inel_mideg->Fill(kecalo, mom_rw_minchi2);
					h1d_kend_bb_inel_mideg->Fill(kebb, mom_rw_minchi2);

				}
				if (kMIDother) { 
					h1d_keffbeam_inel_midother->Fill(ke_ffbeam_MeV, mom_rw_minchi2);
					h1d_kend_calo_inel_midother->Fill(kecalo, mom_rw_minchi2);
					h1d_kend_bb_inel_midother->Fill(kebb, mom_rw_minchi2);
				}

			} //reco inel
			} //beam xy
		} //basic cuts

	} //main entry loop

	//Fit Gaussians on momenta ...
	TF1* kebeam_fit; kebeam_fit=VFit(h1d_kebeam, 2);
	kebeam_fit->SetName("kebeam_fit");
	TF1* pbeam_fit; pbeam_fit=VFit(h1d_pbeam, 2);
	pbeam_fit->SetName("pbeam_fit");

	TF1* pbeam_stop_fit; pbeam_stop_fit=VFit(h1d_pbeam_stop, 2);
	pbeam_stop_fit->SetName("pbeam_stop_fit");

	TF1* pff_stop_fit; pff_stop_fit=VFit(h1d_pff_stop, 2);
	pff_stop_fit->SetName("pff_stop_fit");

	//basic bmrw info
	TParameter<Int_t>* bm_nmu=new TParameter<Int_t>("bm_nmu",0.);
	bm_nmu->SetVal(nmu);
	TParameter<Double_t>* bm_dmu=new TParameter<Double_t>("bm_dmu",0.);
	bm_dmu->SetVal(dmu);
	TParameter<Double_t>* bm_mu_st=new TParameter<Double_t>("bm_mu_st",0.);
	bm_mu_st->SetVal(mu_st);

	TParameter<Int_t>* bm_nsigma=new TParameter<Int_t>("bm_nsigma",0.);
	bm_nsigma->SetVal(nsigma);
	TParameter<Double_t>* bm_dsigma=new TParameter<Double_t>("bm_dsigma",0.);
	bm_dsigma->SetVal(dsigma);
	TParameter<Double_t>* bm_sigma_st=new TParameter<Double_t>("bm_sigma_st",0.);
	bm_sigma_st->SetVal(sigma_st);


	TParameter<Double_t>* mu_spec=new TParameter<Double_t>("mu_spec",0.);
	mu_spec->SetVal(m1);

	TParameter<Double_t>* sigma_spec=new TParameter<Double_t>("sigma_spec",0.);
	sigma_spec->SetVal(m1);


	//save results...
   	//TFile *fout = new TFile("mc_proton_beamxy_beammom_afterbmrw.root","RECREATE");
   	//TFile *fout = new TFile("mc_proton_beamxy_beammom_afterbmrw_usespecasinput.root","RECREATE");
   	//TFile *fout = new TFile("mc_proton_beamxy_beammom_calo_afterbmrw.root","RECREATE");
   	//TFile *fout = new TFile("mc_proton_beamxy_beammom_calo_rmxtrack_afterbmrw.root","RECREATE");
   	//TFile *fout = new TFile("mc_proton_beamxy_calo_afterkerw.root","RECREATE");
   	//TFile *fout = new TFile("mc_proton_beamxy_beammom_calo_afterkerw.root","RECREATE");
   	//TFile *fout = new TFile("mc_proton_beamxy_beammom_nobmrw.root","RECREATE");
   	//TFile *fout = new TFile("mc_proton_beamxy_beammom_nobmrw_new.root","RECREATE");
   	//TFile *fout = new TFile("mc_proton_beamxy_beammom_nobmrw_new2.root","RECREATE");
   	//TFile *fout = new TFile("mc_proton_beamxy_beammom_nobmrw_effstudy.root","RECREATE");
   	//TFile *fout = new TFile("mc_proton_beamxy_beammom_nobmrw_effstudy_withLikelihoodFit.root","RECREATE");
   	//TFile *fout = new TFile("mc_proton_beamxy_beammom_nobmrw.root","RECREATE");
   	//TFile *fout = new TFile("mc_proton_beamxy_beammom_nobmrw_noRcorr.root","RECREATE");
   	//TFile *fout = new TFile("mc_proton_beamxy_beammom_bmrw_new.root","RECREATE");
   	TFile *fout = new TFile("mc_proton_beamxy_beammom_bmrw_by_kebeamff.root","RECREATE");
		bm_nmu->Write();
		bm_dmu->Write();
		bm_mu_st->Write();

		bm_nsigma->Write();
		bm_dsigma->Write();
		bm_sigma_st->Write();
			
		h2d_rr_dedx_recoSTOP->Write();
		gr_predict_dedx_resrange->Write();

		h1d_ke0->Write();
		h1d_ke0_el->Write();
		h1d_ke0_inel->Write();
		h1d_ke0_stop->Write();
		h1d_p0->Write();
		h1d_p0_stop->Write();

		h1d_pbeam->Write();
		h1d_pbeam_stop->Write();

		h1d_keff->Write();
		h1d_keff_stop->Write();
		h1d_keff_el->Write();
		h1d_keff_inel->Write();
		h1d_pff->Write();
		h1d_pff_stop->Write();

		h1d_prange_stop->Write();
		h1d_pcalo_stop->Write();

		h1d_kebeam->Write();
		h1d_kebeam_stop->Write();
		kebeam_fit->Write();
		pbeam_fit->Write();
		pbeam_stop_fit->Write();
		pff_stop_fit->Write();

		h1d_kerange_stop->Write();
		h1d_kecalo_stop->Write();

		h2d_xy_noSCE->Write();
		h2d_xy_SCE->Write();
		
		h1d_zst_noSCE->Write();
		h1d_zst_SCE->Write();

		chi2pid_recostop->Write();
		chi2pid_recoinel->Write();
		//chi2pid_truestop->Write();
		chi2pid_trueel->Write();
		chi2pid_trueinel->Write();

		h1d_trklen_stop->Write();
		//h1d_trklen_stop_XY->Write();

		//for (int ig = 0; ig < n_1d; ++ig) { //rw loop
			//h1d_trklen_rw[ig]->Write();
		//} //rw loop


		h1d_trklen->Write();
		h1d_trklen_bmrw->Write();
		//h1d_trklen_XY->Write();
		//h1d_trklen_bmrw_XY->Write();
	
		h1d_zend->Write();
		h1d_zend_bmrw->Write();
		//h1d_zend_XY->Write();
		//h1d_zend_bmrw_XY->Write();


		h1d_keffbeam_stop->Write();
		h1d_keffbeam_inel->Write();
		h1d_keffbeam_inel_inel->Write();
		h1d_keffbeam_inel_el->Write();
		h1d_keffbeam_inel_midcosmic->Write();
		h1d_keffbeam_inel_midpi->Write();
		h1d_keffbeam_inel_midp->Write();
		h1d_keffbeam_inel_midmu->Write();
		h1d_keffbeam_inel_mideg->Write();
		h1d_keffbeam_inel_midother->Write();

		h1d_keffbeam_el_noxy->Write();
		h1d_keffbeam_el->Write();
		h1d_keffbeam_el_inel->Write();
		h1d_keffbeam_el_el->Write();
		h1d_keffbeam_el_midcosmic->Write();
		h1d_keffbeam_el_midpi->Write();
		h1d_keffbeam_el_midp->Write();
		h1d_keffbeam_el_midmu->Write();
		h1d_keffbeam_el_mideg->Write();
		h1d_keffbeam_el_midother->Write();


		h1d_kehy_stop->Write();
		h1d_kehy_inel->Write();
		h1d_kehy_el->Write();

		h1d_kend_calo_el->Write();
		h1d_kend_calo_el_inel->Write();
		h1d_kend_calo_el_el->Write();
		h1d_kend_calo_el_midcosmic->Write();
		h1d_kend_calo_el_midpi->Write();
		h1d_kend_calo_el_midp->Write();
		h1d_kend_calo_el_midmu->Write();
		h1d_kend_calo_el_mideg->Write();
		h1d_kend_calo_el_midother->Write();


		h1d_kend_calo_stop->Write();
		h1d_kend_calo_inel->Write();
		h1d_kend_calo_inel_inel->Write();
		h1d_kend_calo_inel_el->Write();
		h1d_kend_calo_inel_midcosmic->Write();
		h1d_kend_calo_inel_midpi->Write();
		h1d_kend_calo_inel_midp->Write();
		h1d_kend_calo_inel_midmu->Write();
		h1d_kend_calo_inel_mideg->Write();
		h1d_kend_calo_inel_midother->Write();


		h1d_kend_bb_el->Write();
		h1d_kend_bb_el_inel->Write();
		h1d_kend_bb_el_el->Write();
		h1d_kend_bb_el_midcosmic->Write();
		h1d_kend_bb_el_midpi->Write();
		h1d_kend_bb_el_midp->Write();
		h1d_kend_bb_el_midmu->Write();
		h1d_kend_bb_el_mideg->Write();
		h1d_kend_bb_el_midother->Write();



		h1d_kend_bb_stop->Write();
		h1d_kend_bb_inel->Write();
		h1d_kend_bb_inel->Write();
		h1d_kend_bb_inel_inel->Write();
		h1d_kend_bb_inel_el->Write();
		h1d_kend_bb_inel_midcosmic->Write();
		h1d_kend_bb_inel_midpi->Write();
		h1d_kend_bb_inel_midp->Write();
		h1d_kend_bb_inel_midmu->Write();
		h1d_kend_bb_inel_mideg->Write();
		h1d_kend_bb_inel_midother->Write();


		h1d_kend_true_el->Write();
		h1d_kend_true_stop->Write();
		h1d_kend_true_inel->Write();

		h2d_trklen_ratio_KEbbfit_KEffbeam_stop->Write();
		h2d_trklen_ratio_KEbbfit_KEffbeam_el->Write();
		h2d_trklen_ratio_KEbbfit_KEffbeam_inel->Write();

		h2d_trklen_ratio_KEbbtruth_KEffbeam_stop->Write();
		h2d_trklen_ratio_KEbbtruth_KEffbeam_el->Write();
		h2d_trklen_ratio_KEbbtruth_KEffbeam_inel->Write();

		h2d_keffbeam_keff_stop->Write();
		h2d_keffbeam_keff_el->Write();
		h2d_keffbeam_keff_inel->Write();

		h2d_kehy_keff_stop->Write();
		h2d_kehy_keff_el->Write();
		h2d_kehy_keff_inel->Write();

		h1d_kend_bbtrue_stop->Write();
		h1d_kend_bbtrue_el->Write();
		h1d_kend_bbtrue_inel->Write();


		h2d_kebb_keendtruth_stop->Write();
		h2d_kebb_keendtruth_el->Write();
		h2d_kebb_keendtruth_inel->Write();

		h2d_kebb_kebbtruth_stop->Write();
		h2d_kebb_kebbtruth_el->Write();
		h2d_kebb_kebbtruth_inel->Write();

		h2d_kecalo_keendtruth_stop->Write();
		h2d_kecalo_keendtruth_el->Write();
		h2d_kecalo_keendtruth_inel->Write();
	
		h2d_kebbtruelength_keendtruth_stop->Write();
		h2d_kebbtruelength_keendtruth_el->Write();
		h2d_kebbtruelength_keendtruth_inel->Write();
		h2d_kebbtruelength_keendtruth_trueel->Write();

		h1d_kend_bbtruelength_stop->Write();
		h1d_kend_bbtruelength_el->Write();
		h1d_kend_bbtruelength_inel->Write();

		h1d_ratio_rangetrue_rangereco_stop->Write();
		h1d_ratio_rangetrue_rangereco_el->Write();
		h1d_ratio_rangetrue_rangereco_inel->Write();

		h1d_ratio_rangehy_rangereco_stop->Write();
		h1d_ratio_rangehy_rangereco_el->Write();

		h2d_recorange_truerange_stop->Write();
		h2d_recorange_truerange_el->Write();
		h2d_recorange_truerange_inel->Write();

		h2d_recorange_rangehy_el->Write();
		h2d_recorange_rangehy_stop->Write();

		h1d_kend_bbcorr_inel->Write();
		h1d_kend_bbcorr_el->Write();

		h1d_dke_keffbeam_keff_el->Write();
		h1d_dke_keffbeam_keff_inel->Write();

		h2d_trklen_eff_KEconst_KEff_all->Write();
		h2d_trklen_eff_KEhy_KEff_all->Write();

		h2d_trklen_eff_KEconst_KEff_stop->Write();
		h2d_trklen_eff_KEhy_KEff_stop->Write();

		h2d_trklen_eff_KEconst_KEff_el->Write();
		h2d_trklen_eff_KEhy_KEff_el->Write();

		h2d_trklen_eff_KEconst_KEff_inel->Write();
		h2d_trklen_eff_KEhy_KEff_inel->Write();

		h2d_trklen_eff_KEhy_KEconst_all->Write();
		h2d_trklen_eff_KEhy_KEconst_stop->Write();
		h2d_trklen_eff_KEhy_KEconst_el->Write();
		h2d_trklen_eff_KEhy_KEconst_inel->Write();

	fout->Close();



}
