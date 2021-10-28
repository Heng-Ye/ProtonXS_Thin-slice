#define ProtonTrueLen_cxx
#include "ProtonTrueLen.h"

#include <TH2.h>
#include <TH1.h>
#include <TH3.h>
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
//#include <TH3D.h>
//#include <TH3F.h>
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
#include "TGraphAsymmErrors.h"

#include "./cali/dedx_function_35ms.h"
#include "./headers/BasicParameters.h"
#include "./headers/BasicFunctions.h"
#include "./headers/util.h"
//#include "./headers/SliceParams.h"
//#include "./headers/ThinSlice.h"
//#include "./headers/sce_map.h"

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

void ProtonTrueLen::Loop() {
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	Long64_t nbytes = 0, nb = 0;
	bool isTestSample=true;
	//int true_sliceID = -1, reco_sliceID = -1;


	//book histograms --------------------------------------------------------------------------------------------//
	int n_b=150;
	double b_min=0;
	double b_max=150;

	int n_dr=300;
	double dr_min=0;
	double dr_max=150;

	int n_dcos=100;
	double dcos_min=0;
	double dcos_max=1;

	TH1D *h1d_truetrklen_NoCut_inel=new TH1D(Form("h1d_truetrklen_NoCut_inel"), Form(""), n_b, b_min, b_max);	
	TH1D *h1d_truetrklen_Pan_inel=new TH1D(Form("h1d_truetrklen_Pan_inel"), Form(""), n_b, b_min, b_max);	
	TH1D *h1d_truetrklen_CaloSz_inel=new TH1D(Form("h1d_truetrklen_CaloSz_inel"), Form(""), n_b, b_min, b_max);	
	TH1D *h1d_truetrklen_Pos_inel=new TH1D(Form("h1d_truetrklen_Pos_inel"), Form(""), n_b, b_min, b_max);	
	TH1D *h1d_truetrklen_BQ_inel=new TH1D(Form("h1d_truetrklen_BQ_inel"), Form(""), n_b, b_min, b_max);	
	TH1D *h1d_truetrklen_RecoInel_inel=new TH1D(Form("h1d_truetrklen_RecoInel_inel"), Form(""), n_b, b_min, b_max);	
	TH1D *h1d_truetrklen_BeamMatch_inel=new TH1D(Form("h1d_truetrklen_BeamMatch_inel"), Form(""), n_b, b_min, b_max);
	h1d_truetrklen_NoCut_inel->Sumw2();
	h1d_truetrklen_Pan_inel->Sumw2();
	h1d_truetrklen_CaloSz_inel->Sumw2();
	h1d_truetrklen_Pos_inel->Sumw2();
	h1d_truetrklen_BQ_inel->Sumw2();
	h1d_truetrklen_RecoInel_inel->Sumw2();
	h1d_truetrklen_BeamMatch_inel->Sumw2();


	int n_endz=200;
	double endz_min=-50;
	double endz_max=150;
	TH2D *h2d_truetrklen_trueendz_CaloSz_inel=new TH2D(Form("h2d_truetrklen_trueendz_CaloSz_inel"), Form(""), n_b, b_min, b_max, n_endz, endz_min, endz_max);	
	TH1D *h1d_trueendz_CaloSz_inel=new TH1D(Form("h1d_trueendz_CaloSz_inel"), Form(""), n_endz, endz_min, endz_max);	


	TH2D *h2d_truetrklen_ketrue_NoCut_inel=new TH2D(Form("h2d_truetrklen_ketrue_NoCut_inel"), Form(""), n_b, b_min, b_max,500,0,500);	
	TH2D *h2d_truetrklen_ketrue_Pan_inel=new TH2D(Form("h2d_truetrklen_ketrue_Pan_inel"), Form(""), n_b, b_min, b_max,500,0,500);	
	TH2D *h2d_truetrklen_ketrue_CaloSz_inel=new TH2D(Form("h2d_truetrklen_ketrue_CaloSz_inel"), Form(""), n_b, b_min, b_max,500,0,500);	

	TH2D *h2d_truetrklen_keff_NoCut_inel=new TH2D(Form("h2d_truetrklen_keff_NoCut_inel"), Form(""), n_b, b_min, b_max,500,0,500);	

	TH2D *h2d_trueEndZ_ketrue_NoCut_inel=new TH2D(Form("h2d_trueEndZ_ketrue_NoCut_inel"), Form(""), 350, -200, 150,500,0,500);	
	TH2D *h2d_trueZ_ketrue_NoCut_inel=new TH2D(Form("h2d_trueZ_ketrue_NoCut_inel"), Form(""), 3500, -200, 150, 800, 0, 800);	

	int n_x=110*5;
	double x_min=-60;
	double x_max=50;

	int n_y=110*5;
	double y_min=390;
	double y_max=500;

	const int n_scan=25; 
	TH2D *h2d_trueXY_NoCut_inel[n_scan];
	//double truelen_st=-10; //scan from truelen=0
	double truelen_st=-2; //scan from truelen=0
	//double truelen_st=-60; //scan from truelen=0
	double d_truelen=2; //every 2 cm
	vector<double> seg_true_len;
	for (int j=0; j<n_scan; ++j) {
		double seg_st=truelen_st+(double)j*d_truelen;
		seg_true_len.push_back(seg_st);
		h2d_trueXY_NoCut_inel[j]=new TH2D(Form("h2d_trueXY_NoCut_inel_%d",j), Form("EndZ:%.1f-%.1f cm",seg_st,seg_st-d_truelen), n_x, x_min, x_max, n_y, y_min, y_max);
	}


	//TH2D *h2d_trueXY_truelen05_NoCut_inel=new TH2D(Form("h2d_trueXY_truelen05_NoCut_inel"), Form(""), n_x, x_min, x_max, n_y, y_min, y_max);	
	//TH2D *h2d_trueXY_truelen510_NoCut_inel=new TH2D(Form("h2d_trueXY_truelen510_NoCut_inel"), Form(""), n_x, x_min, x_max, n_y, y_min, y_max);	
	//TH2D *h2d_trueXY_truelen1015_NoCut_inel=new TH2D(Form("h2d_trueXY_truelen1015_NoCut_inel"), Form(""), n_x, x_min, x_max, n_y, y_min, y_max);	
	//TH2D *h2d_trueXY_truelen1520_NoCut_inel=new TH2D(Form("h2d_trueXY_truelen1520_NoCut_inel"), Form(""), n_x, x_min, x_max, n_y, y_min, y_max);	
	//TH2D *h2d_trueXY_truelen2025_NoCut_inel=new TH2D(Form("h2d_trueXY_truelen2025_NoCut_inel"), Form(""), n_x, x_min, x_max, n_y, y_min, y_max);	
	//TH2D *h2d_trueXY_truelen2530_NoCut_inel=new TH2D(Form("h2d_trueXY_truelen2530_NoCut_inel"), Form(""), n_x, x_min, x_max, n_y, y_min, y_max);	
	//TH2D *h2d_trueXY_truelen3035_NoCut_inel=new TH2D(Form("h2d_trueXY_truelen3035_NoCut_inel"), Form(""), n_x, x_min, x_max, n_y, y_min, y_max);	

	//TH3D *h3d_trueXYZ_NoCut_inel=new TH3D(Form("h3d_trueXYZ_NoCut_inel"), Form(""), n_x, x_min, x_max, n_y, y_min, y_max,1200,0,120);	



	//theta_yz vs theta_xz
	int n_thetaYZ=3600;
	double thetaYZ_min=-180;
	double thetaYZ_max=180;
	int n_thetaXZ=3600;
	double thetaXZ_min=-180;
	double thetaXZ_max=180;
	TProfile2D* h2d_thetaXZ_thetaYZ_BQ_inel=new TProfile2D("h2d_thetaXZ_thetaYZ_BQ_inel","",n_thetaXZ,thetaXZ_min,thetaXZ_max,n_thetaYZ,thetaYZ_min,thetaYZ_max);
	TProfile2D* h2d_thetaXZ_thetaYZ_BQ_misidp=new TProfile2D("h2d_thetaXZ_thetaYZ_BQ_misidp","",n_thetaXZ,thetaXZ_min,thetaXZ_max,n_thetaYZ,thetaYZ_min,thetaYZ_max);
	TProfile2D* h2d_thetaXZ_thetaYZ_BQ_el=new TProfile2D("h2d_thetaXZ_thetaYZ_BQ_el","",n_thetaXZ,thetaXZ_min,thetaXZ_max,n_thetaYZ,thetaYZ_min,thetaYZ_max);

	TProfile2D* h2d_thetaXZ_thetaYZ_Pos_inel=new TProfile2D("h2d_thetaXZ_thetaYZ_Pos_inel","",n_thetaXZ,thetaXZ_min,thetaXZ_max,n_thetaYZ,thetaYZ_min,thetaYZ_max);
	TProfile2D* h2d_thetaXZ_thetaYZ_Pos_misidp=new TProfile2D("h2d_thetaXZ_thetaYZ_Pos_misidp","",n_thetaXZ,thetaXZ_min,thetaXZ_max,n_thetaYZ,thetaYZ_min,thetaYZ_max);
	TProfile2D* h2d_thetaXZ_thetaYZ_Pos_el=new TProfile2D("h2d_thetaXZ_thetaYZ_Pos_el","",n_thetaXZ,thetaXZ_min,thetaXZ_max,n_thetaYZ,thetaYZ_min,thetaYZ_max);
	//------------------------------------------------------------------------------------------------------------//

        //MC Beam Mom Gaussian 
        double m1=1007.1482; //MC prod4a [spec]
        double s1=60.703307; //MC prod4a [spec]

        //momentum cut range    
        double mu_min=m1-3.*s1;
        double mu_max=m1+3.*s1;

	for (Long64_t jentry=0; jentry<nentries;jentry++) { //main entry loop
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

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

			if (reco_stz>reco_endz) {
				reco_endx=primtrk_hitx->at(0); 
				reco_endy=primtrk_hity->at(0);
				reco_endz=primtrk_hitz->at(0);

				reco_stx=primtrk_hitx->at(primtrk_dedx->size()-1);	
				reco_sty=primtrk_hity->at(primtrk_dedx->size()-1);
				reco_stz=primtrk_hitz->at(primtrk_dedx->size()-1);
			}

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
		double thetaYZ_deg=-999;
		double thetaXZ_deg=-999;
		if (IsCaloSize) {	
			//trk direction after SCE corr.
			TVector3 pt0(primtrk_hitx->at(0), primtrk_hity->at(0), primtrk_hitz->at(0));
			TVector3 pt1(primtrk_hitx->at(-1+primtrk_hitx->size()), primtrk_hity->at(-1+primtrk_hity->size()), primtrk_hitz->at(-1+primtrk_hitz->size()));
			dir = pt1 - pt0;
			dir = dir.Unit();

			//beam direction
      			//TVector3 beamdir(cos(beam_angleX_mc*TMath::Pi()/180), cos(beam_angleY_mc*TMath::Pi()/180), cos(beam_angleZ_mc*TMath::Pi()/180));
      			TVector3 beamdir(beamDirx_spec->at(0),beamDiry_spec->at(0),beamDirz_spec->at(0));
      			beamdir = beamdir.Unit();
      			//beam_costh = dir.Dot(beamdir);
      			cosine_beam_spec_primtrk=dir.Dot(beamdir);

			//thetayz ---------------------------------------------------------------------------------------------------//
			TVector3 dir_yz;
			TVector3 pt0_yz(0, primtrk_hity->at(0), primtrk_hitz->at(0));
			TVector3 pt1_yz(0, primtrk_hity->at(-1+primtrk_hity->size()), primtrk_hitz->at(-1+primtrk_hitz->size()));
			dir_yz = pt1_yz - pt0_yz;
			dir_yz = dir_yz.Unit();
      			TVector3 beamdir_yz(0, beamDiry_spec->at(0),beamDirz_spec->at(0));
      			beamdir_yz = beamdir_yz.Unit();
			thetaYZ_deg=(180./TMath::Pi())*TMath::ACos(dir_yz.Dot(beamdir_yz));
			
			//thetaxz ---------------------------------------------------------------------------------------------------//
			TVector3 dir_xz;
			TVector3 pt0_xz(primtrk_hitx->at(0), 0, primtrk_hitz->at(0));
			TVector3 pt1_xz(primtrk_hitx->at(-1+primtrk_hitx->size()), 0, primtrk_hitz->at(-1+primtrk_hitz->size()));
			dir_xz = pt1_xz - pt0_xz;
			dir_xz = dir_xz.Unit();
      			TVector3 beamdir_xz(beamDirx_spec->at(0), 0, beamDirz_spec->at(0));
      			beamdir_xz = beamdir_xz.Unit();
			thetaXZ_deg=(180./TMath::Pi())*TMath::ACos(dir_xz.Dot(beamdir_xz));
		}

		if (cosine_beam_spec_primtrk<0) { cosine_beam_spec_primtrk=-1.*cosine_beam_spec_primtrk; }
		if (cosine_beam_spec_primtrk>cosine_beam_primtrk_min) { IsCosine=true; }

		//thetayz, thetaxz


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

			trkdedx.push_back(cali_dedx);
			trkres.push_back(resrange_reco);

		  } //loop over reco hits of a given track
		  //range_reco=primtrk_range->at(0);

		  pid=chi2pid(trkdedx,trkres); //pid using stopping proton hypothesis
		} //if calo size not empty


		//Reco stopping/Inel p cut
		bool IsRecoStop=false;
		bool IsRecoInEL=false;
		double mom_beam_spec=-99; mom_beam_spec=beamMomentum_spec->at(0);
		//double range_reco=-99; if (!primtrk_range->empty()) range_reco=primtrk_range->at(0); //reco primary trklen
		double csda_val_spec=csda_range_vs_mom_sm->Eval(mom_beam_spec);

		//if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) IsRecoStop=true;
		//if ((range_reco/csda_val_spec)<min_norm_trklen_csda) IsRecoInEL=true;

		if ((range_reco/csda_val_spec)<min_norm_trklen_csda) { //inel region
			if (pid>pid_1) IsRecoInEL=true; 
			if (pid<=pid_1) IsRecoStop=true; 
		} //inel region
		if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) { //stopping p region
			if (pid>pid_2) IsRecoInEL=true; 
			if (pid<=pid_2) IsRecoStop=true;
		} //stopping p region


		//kinetic energies
		double ke_beam_spec=p2ke(mom_beam_spec); //ke_beam_spec [GeV]
		double ke_beam_spec_MeV=1000.*ke_beam_spec; //ke_beam_spec [MeV]
		double ke_trklen=ke_vs_csda_range_sm->Eval(range_reco); //[unit: GeV]
		double ke_trklen_MeV=1000.*ke_trklen; //[unit: MeV]
		double ke_calo_MeV=0;

		//First point of MCParticle entering TPC ------------------------------------------------------------------------//
		bool is_beam_at_ff=false;
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

		//Get true trklen ---------------------------------------------------------------------------------------//
		double range_true=-999;
		int key_st = 0;
		double tmp_z = 9999;
		vector<double> true_trklen_accum;
		//cout<<"beamtrk_z->size():"<<beamtrk_z->size()<<endl;
		true_trklen_accum.reserve(beamtrk_z->size()); // initialize true_trklen_accum
		//cout<<"ck0"<<endl;
		for (int iz=0; iz<(int)beamtrk_z->size(); iz++) {
			if (abs(beamtrk_z->at(iz)) < tmp_z){
				tmp_z = abs(beamtrk_z->at(iz));
				key_st = iz; // find the point where the beam enters the TPC (find the smallest abs(Z))
			}
			//cout<<"ck0/"<<endl;
			true_trklen_accum[iz] = 0.; // initialize true_trklen_accum
			//cout<<"ck0///"<<endl;
			//if (iz==-1+(int)beamtrk_z->size()) 
		}
		//cout<<"key_st:"<<key_st<<"\n"<<endl;
		//cout<<"ck1"<<endl;
		//for (int iz=key_st+1; iz<(int)beamtrk_z->size(); iz++){

		if (is_beam_at_ff) { //if primary protons entering tpc
			for (int iz=key_reach_tpc+1; iz<(int)beamtrk_z->size(); iz++){
				if (iz == key_reach_tpc+1) range_true = 0;
					range_true += sqrt( pow(beamtrk_x->at(iz)-beamtrk_x->at(iz-1), 2)+
					pow(beamtrk_y->at(iz)-beamtrk_y->at(iz-1), 2)+	
					pow(beamtrk_z->at(iz)-beamtrk_z->at(iz-1), 2) );						    	
					true_trklen_accum[iz] = range_true;


			}
		} //if primary protons entering tpc

		double zproj_beam=0; //set beam z at ff
		double yproj_beam=0; //ini. value
		double xproj_beam=0; //ini. value
		int n_fit=3; //num of points used for fitting
                if (beamtrk_z->size()) {

			int key_fit_st=key_reach_tpc-1;
			int key_fit_ed=key_reach_tpc+1;

			if (key_fit_st<0) key_fit_st=0;
			if (key_fit_ed>(-1+(int)beamtrk_z->size())) key_fit_ed=-1+(int)beamtrk_z->size();	
			//std::cout<<"key_fit_st-ed:"<<key_fit_st<<"-"<<key_fit_ed<<std::endl;

			//start 3D line fit
			TGraph2D *gr=new TGraph2D();
			//cout<<"ck0"<<endl;
   		  	//for (int N=key_fit_st; N<key_fit_ed; N++) {
   		  	int nsize_fit=n_fit;
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
			cout<<"ck7"<<endl;

			delete gr;
		}
	
		//range compensation ------------------------------------------------------------//
		if (is_beam_at_ff) {
			range_true += sqrt( pow(beamtrk_x->at(key_reach_tpc)-xproj_beam, 2)+
					pow(beamtrk_y->at(key_reach_tpc)-yproj_beam, 2)+	
					pow(beamtrk_z->at(key_reach_tpc)-zproj_beam, 2) );						    	
					//true_trklen_accum[iz] = range_true; //FIXME
		}


		//if (run==18800167&&subrun==3&&event==400&&primaryID==64) {
		//if (IsPureInEL&&range_true>4&&range_true<5.5) {
			//cout<<"\n\n key_reach_tpc:"<<key_reach_tpc<<endl;
			//cout<<"range_true:"<<range_true<<endl;
			//cout<<"true_endz:"<<true_endz<<endl;
			//cout<<"run:"<<run<<" subrun:"<<subrun<<" event:"<<event<<" primaryID:"<<primaryID<<endl;
			//for (int iz=0; iz<(int)beamtrk_z->size(); iz++){
				//cout<<"beamtrk_z["<<iz<<"]:"<<beamtrk_z->at(iz)<<endl;
			//}
		//}
		
		//dr=true(x,y,z) - reco(x,y,z)
		bool IsGood_Inel=false;
		//cout<<"true_endz/y/x:"<<true_endz<<"/"<<true_endy<<"/"<<true_endx<<endl;
		//cout<<"beamtrk_z/y/x:"<<beamtrk_z->at(-1+beamtrk_z->size())<<"/"<<beamtrk_y->at(-1+beamtrk_y->size())<<"/"<<beamtrk_x->at(-1+beamtrk_x->size())<<endl;
		double dr=sqrt(pow(reco_endx-true_endx,2)+pow(reco_endy-true_endy,2)+pow(reco_endz-true_endz,2));
		//if (dr<=dr_cut_inel) IsGood_Inel=true;

		//KE_true
		//cout<<"\nsize of beamtrk_Eng:"<<beamtrk_Eng->size()<<endl;
		//cout<<"primtrk_edept_true->size():"<<primtrk_edept_true->size()<<endl;

		double KE_ff=-999;
		if (is_beam_at_ff) { 
			KE_ff=1000.*(beamtrk_Eng->at(key_reach_tpc)); //MeV

			//for (int kk=0; kk<(int)beamtrk_Eng->size(); ++kk) {
				//cout<<"beamtrk_Eng["<<kk<<"]:"<<1000.*beamtrk_Eng->at(kk)<<" | beamtrk_z["<<kk<<"]:"<<beamtrk_z->at(kk)<<endl;
			//}
		}
		//cout<<"KE_ff:"<<KE_ff<<endl;
		double KE_true=KE_ff;
		
		if (primtrk_edept_true->size()) {
			for (int kk=0; kk<(int)primtrk_edept_true->size(); ++kk) {
				KE_true-=primtrk_edept_true->at(kk);
				//cout<<"primtrk_edept_true["<<kk<<"]:"<<primtrk_edept_true->at(kk)<<" | primtrk_hitz_true["<<kk<<"]:"<<primtrk_hitz_true->at(kk)<<endl; //simide
			}
		
		}

		//key_reach_tpc
		//ke_ff=1000.*(beamtrk_Eng->at(key_reach_tpc)); //MeV

	   	//for (size_t h=0; h<primtrk_edept_true->at(k).size(); ++h) { //loop over true hits of a given track
	     	//if ((primtrk_trkid_true->at(k)[h])!=-1) { //remove shower	
	        //ke_true-=primtrk_edept_true->at(k)[h]; //ke_true of each hit


		if (IsPureInEL) { //true inel
			Fill1DHist(h1d_truetrklen_NoCut_inel, range_true);
			double tmp_range_true=range_true;
			if (tmp_range_true<0) tmp_range_true=0;

			h2d_truetrklen_ketrue_NoCut_inel->Fill(tmp_range_true, KE_true);
			h2d_truetrklen_keff_NoCut_inel->Fill(tmp_range_true, KE_ff);
			h2d_trueEndZ_ketrue_NoCut_inel->Fill(beamtrk_z->at(-1+beamtrk_z->size()), KE_true);

			//if (is_beam_at_ff) { //if primary protons entering tpc
			if (beamtrk_z->size()) {
				//for (int iz=key_reach_tpc+1; iz<(int)beamtrk_z->size(); iz++){
				for (int iz=0; iz<(int)beamtrk_z->size(); iz++){
					h2d_trueZ_ketrue_NoCut_inel->Fill(beamtrk_z->at(iz), 1000.*(beamtrk_Eng->at(iz)));

					//if (tmp_range_true>=0&&tmp_range_true<5) h2d_trueXY_truelen05_NoCut_inel->Fill(beamtrk_x->at(iz), beamtrk_y->at(iz));
					//if (tmp_range_true>=5&&tmp_range_true<10) h2d_trueXY_truelen510_NoCut_inel->Fill(beamtrk_x->at(iz), beamtrk_y->at(iz));
					//if (tmp_range_true>=10&&tmp_range_true<15) h2d_trueXY_truelen1015_NoCut_inel->Fill(beamtrk_x->at(iz), beamtrk_y->at(iz));
					//if (tmp_range_true>=15&&tmp_range_true<20) h2d_trueXY_truelen1520_NoCut_inel->Fill(beamtrk_x->at(iz), beamtrk_y->at(iz));
					//if (tmp_range_true>=20&&tmp_range_true<25) h2d_trueXY_truelen2025_NoCut_inel->Fill(beamtrk_x->at(iz), beamtrk_y->at(iz));
					//if (tmp_range_true>=25&&tmp_range_true<30) h2d_trueXY_truelen2530_NoCut_inel->Fill(beamtrk_x->at(iz), beamtrk_y->at(iz));
					//if (tmp_range_true>=30&&tmp_range_true<35) h2d_trueXY_truelen3035_NoCut_inel->Fill(beamtrk_x->at(iz), beamtrk_y->at(iz));
					//h3d_trueXYZ_NoCut_inel->Fill(beamtrk_x->at(iz), beamtrk_y->at(iz), beamtrk_z->at(iz));
					//for (int j=0; j<n_scan; ++j) {
						//double min=seg_true_len.at(j)-d_truelen;
						//double max=seg_true_len.at(j);
						//if (tmp_range_true>=min&&tmp_range_true<max) h2d_trueXY_NoCut_inel[j]->Fill(beamtrk_x->at(iz), beamtrk_y->at(iz));
						//if (true_endz>=min&&true_endz<max) h2d_trueXY_NoCut_inel[j]->Fill(beamtrk_x->at(iz), beamtrk_y->at(iz));
						//if (beamtrk_z->at(iz)>=min&&beamtrk_z->at(iz)<max) h2d_trueXY_NoCut_inel[j]->Fill(beamtrk_x->at(iz), beamtrk_y->at(iz));
					//}
				}
			}
			//} //if primary protons entering tpc




			if (IsPandoraSlice) {
				Fill1DHist(h1d_truetrklen_Pan_inel, range_true);
				h2d_truetrklen_ketrue_Pan_inel->Fill(range_true, KE_true);
				if (IsCaloSize) {
					Fill1DHist(h1d_truetrklen_CaloSz_inel, range_true);
					h2d_truetrklen_ketrue_CaloSz_inel->Fill(range_true, KE_true);

					//double tmp_trueendz=beamtrk_z->at(key_reach_tpc+1);
					double tmp_trueendz=true_endz;
					if (tmp_trueendz<=endz_min) tmp_trueendz=endz_min;
					h2d_truetrklen_trueendz_CaloSz_inel->Fill(range_true, tmp_trueendz);
					Fill1DHist(h1d_trueendz_CaloSz_inel, true_endz);
					if (IsPos) {
						Fill1DHist(h1d_truetrklen_Pos_inel, range_true);
						if (IsBQ) {
							Fill1DHist(h1d_truetrklen_BQ_inel, range_true);
							if (IsRecoInEL) {
								Fill1DHist(h1d_truetrklen_RecoInel_inel, range_true);
								if (IsBeamMatch) {
									Fill1DHist(h1d_truetrklen_BeamMatch_inel, range_true);
								}
							}
						}
				}	}
			}
		} //true inel


		//theta_xz vs theta_yz
		if (IsBQ&&IsCaloSize&&IsPandoraSlice) { //BQ
			if (kinel) {
				h2d_thetaXZ_thetaYZ_BQ_inel->Fill(thetaXZ_deg, thetaYZ_deg, cosine_beam_spec_primtrk);
			}
			if (kel) {
				h2d_thetaXZ_thetaYZ_BQ_el->Fill(thetaXZ_deg, thetaYZ_deg, cosine_beam_spec_primtrk);
			}
			if (kMIDp) {
				h2d_thetaXZ_thetaYZ_BQ_misidp->Fill(thetaXZ_deg, thetaYZ_deg, cosine_beam_spec_primtrk);
			}
		} //BQ
		
		if (IsPos&&IsCaloSize&&IsPandoraSlice) { //Pos
			if (kinel) {
				h2d_thetaXZ_thetaYZ_Pos_inel->Fill(thetaXZ_deg, thetaYZ_deg, cosine_beam_spec_primtrk);
			}
			if (kel) {
				h2d_thetaXZ_thetaYZ_Pos_el->Fill(thetaXZ_deg, thetaYZ_deg, cosine_beam_spec_primtrk);
			}
			if (kMIDp) {
				h2d_thetaXZ_thetaYZ_Pos_misidp->Fill(thetaXZ_deg, thetaYZ_deg, cosine_beam_spec_primtrk);
			}
		} //Pos





	} //main entry loop


	//eff.
	TH1D *h1d_eff_NoCut_inel=(TH1D*)h1d_truetrklen_NoCut_inel->Clone();
	h1d_eff_NoCut_inel->Divide(h1d_eff_NoCut_inel);
	h1d_eff_NoCut_inel->SetName("h1d_eff_NoCut_inel");

	TH1D *h1d_eff_Pan_inel=(TH1D*)h1d_truetrklen_Pan_inel->Clone();
	h1d_eff_Pan_inel->Divide(h1d_truetrklen_NoCut_inel);
	h1d_eff_Pan_inel->SetName("h1d_eff_Pan_inel");

	TH1D *h1d_eff_CaloSz_inel=(TH1D*)h1d_truetrklen_CaloSz_inel->Clone();
	h1d_eff_CaloSz_inel->Divide(h1d_truetrklen_NoCut_inel);
	h1d_eff_CaloSz_inel->SetName("h1d_eff_CaloSz_inel");

	TH1D *h1d_eff_Pos_inel=(TH1D*)h1d_truetrklen_Pos_inel->Clone();
	h1d_eff_Pos_inel->Divide(h1d_truetrklen_NoCut_inel);
	h1d_eff_Pos_inel->SetName("h1d_eff_Pos_inel");

	TH1D *h1d_eff_BQ_inel=(TH1D*)h1d_truetrklen_BQ_inel->Clone();
	h1d_eff_BQ_inel->Divide(h1d_truetrklen_NoCut_inel);
	h1d_eff_BQ_inel->SetName("h1d_eff_BQ_inel");

	TH1D *h1d_eff_RecoInel_inel=(TH1D*)h1d_truetrklen_RecoInel_inel->Clone();
	h1d_eff_RecoInel_inel->Divide(h1d_truetrklen_NoCut_inel);
	h1d_eff_RecoInel_inel->SetName("h1d_eff_RecoInel_inel");

	TH1D *h1d_eff_BeamMatch_inel=(TH1D*)h1d_truetrklen_BeamMatch_inel->Clone();
	h1d_eff_BeamMatch_inel->Divide(h1d_truetrklen_NoCut_inel);
	h1d_eff_BeamMatch_inel->SetName("h1d_eff_BeamMatch_inel");

	//Eff with asym error bars
	vector<double> eff_NoCut_inel;
	vector<double> er_p_NoCut_inel;
	vector<double> er_m_NoCut_inel;

	vector<double> eff_Pan_inel;
	vector<double> er_p_Pan_inel;
	vector<double> er_m_Pan_inel;

	vector<double> x; //true len
	vector<double> errx;

	vector<double> eff_CaloSz_inel;
	vector<double> er_p_CaloSz_inel;
	vector<double> er_m_CaloSz_inel;

	vector<double> eff_Pos_inel;
	vector<double> er_p_Pos_inel;
	vector<double> er_m_Pos_inel;

	vector<double> eff_BQ_inel;
	vector<double> er_p_BQ_inel;
	vector<double> er_m_BQ_inel;

	vector<double> eff_RecoInel_inel;
	vector<double> er_p_RecoInel_inel;
	vector<double> er_m_RecoInel_inel;

	vector<double> eff_BeamMatch_inel;
	vector<double> er_p_BeamMatch_inel;
	vector<double> er_m_BeamMatch_inel;
	
	for(int j=1; j<=h1d_truetrklen_NoCut_inel->GetNbinsX(); ++j) {
		double denom=h1d_truetrklen_NoCut_inel->GetBinContent(j);
		if(denom==0) continue; //no zero in denominator
		
		x.push_back(h1d_truetrklen_NoCut_inel->GetBinCenter(j));
		errx.push_back(0);

		//No cut
		eff_NoCut_inel.push_back(1);
		er_p_NoCut_inel.push_back(0);
		er_m_NoCut_inel.push_back(0);

		//PanS
		double nom_Pan_inel=h1d_truetrklen_Pan_inel->GetBinContent(j);
		eff_Pan_inel.push_back(nom_Pan_inel/denom);
		er_p_Pan_inel.push_back(err_p(nom_Pan_inel,denom));
		er_m_Pan_inel.push_back(err_m(nom_Pan_inel,denom));

		//CaloSz
		double nom_CaloSz_inel=h1d_truetrklen_CaloSz_inel->GetBinContent(j);
		eff_CaloSz_inel.push_back(nom_CaloSz_inel/denom);
		er_p_CaloSz_inel.push_back(err_p(nom_CaloSz_inel,denom));
		er_m_CaloSz_inel.push_back(err_m(nom_CaloSz_inel,denom));

		//Pos
		double nom_Pos_inel=h1d_truetrklen_Pos_inel->GetBinContent(j);
		eff_Pos_inel.push_back(nom_Pos_inel/denom);
		er_p_Pos_inel.push_back(err_p(nom_Pos_inel,denom));
		er_m_Pos_inel.push_back(err_m(nom_Pos_inel,denom));

		//BQ
		double nom_BQ_inel=h1d_truetrklen_BQ_inel->GetBinContent(j);
		eff_BQ_inel.push_back(nom_BQ_inel/denom);
		er_p_BQ_inel.push_back(err_p(nom_BQ_inel,denom));
		er_m_BQ_inel.push_back(err_m(nom_BQ_inel,denom));

		//RecoInel
		double nom_RecoInel_inel=h1d_truetrklen_RecoInel_inel->GetBinContent(j);
		eff_RecoInel_inel.push_back(nom_RecoInel_inel/denom);
		er_p_RecoInel_inel.push_back(err_p(nom_RecoInel_inel,denom));
		er_m_RecoInel_inel.push_back(err_m(nom_RecoInel_inel,denom));

		//Beam Match
		double nom_BeamMatch_inel=h1d_truetrklen_BeamMatch_inel->GetBinContent(j);
		eff_BeamMatch_inel.push_back(nom_BeamMatch_inel/denom);
		er_p_BeamMatch_inel.push_back(err_p(nom_BeamMatch_inel,denom));
		er_m_BeamMatch_inel.push_back(err_m(nom_BeamMatch_inel,denom));
	}

	TGraphAsymmErrors *Eff_NoCut_inel = new TGraphAsymmErrors(eff_NoCut_inel.size(), &x[0], &eff_NoCut_inel[0], &errx[0], &errx[0], &er_m_NoCut_inel[0], &er_m_NoCut_inel[0]);  Eff_NoCut_inel->SetName("Eff_NoCut_inel");
	TGraphAsymmErrors *Eff_Pan_inel = new TGraphAsymmErrors(eff_Pan_inel.size(), &x[0], &eff_Pan_inel[0], &errx[0], &errx[0], &er_m_Pan_inel[0], &er_m_Pan_inel[0]);  Eff_Pan_inel->SetName("Eff_Pan_inel");
	TGraphAsymmErrors *Eff_CaloSz_inel = new TGraphAsymmErrors(eff_CaloSz_inel.size(), &x[0], &eff_CaloSz_inel[0], &errx[0], &errx[0], &er_m_CaloSz_inel[0], &er_m_CaloSz_inel[0]);  Eff_CaloSz_inel->SetName("Eff_CaloSz_inel");
	TGraphAsymmErrors *Eff_Pos_inel = new TGraphAsymmErrors(eff_Pos_inel.size(), &x[0], &eff_Pos_inel[0], &errx[0], &errx[0], &er_m_Pos_inel[0], &er_m_Pos_inel[0]);  Eff_Pos_inel->SetName("Eff_Pos_inel");
	TGraphAsymmErrors *Eff_BQ_inel = new TGraphAsymmErrors(eff_BQ_inel.size(), &x[0], &eff_BQ_inel[0], &errx[0], &errx[0], &er_m_BQ_inel[0], &er_m_BQ_inel[0]);  Eff_BQ_inel->SetName("Eff_BQ_inel");
	TGraphAsymmErrors *Eff_RecoInel_inel = new TGraphAsymmErrors(eff_RecoInel_inel.size(), &x[0], &eff_RecoInel_inel[0], &errx[0], &errx[0], &er_m_RecoInel_inel[0], &er_m_RecoInel_inel[0]);  Eff_RecoInel_inel->SetName("Eff_RecoInel_inel");
	TGraphAsymmErrors *Eff_BeamMatch_inel = new TGraphAsymmErrors(eff_BeamMatch_inel.size(), &x[0], &eff_BeamMatch_inel[0], &errx[0], &errx[0], &er_m_BeamMatch_inel[0], &er_m_BeamMatch_inel[0]);  Eff_BeamMatch_inel->SetName("Eff_BeamMatch_inel");



	//save results ---------------------------------------------------------//
   	//TFile *fout = new TFile("mc_truelen.root","RECREATE");
   	//TFile *fout = new TFile("mc_truelen_around50.root","RECREATE");
   	//TFile *fout = new TFile("mc_truelen_truelen20cm.root","RECREATE");
   	//TFile *fout = new TFile("mc_truelen_keff.root","RECREATE");
   	TFile *fout = new TFile("mc_truelen_new_fixtruelen.root","RECREATE");
		h1d_truetrklen_NoCut_inel->Write();
		h1d_truetrklen_Pan_inel->Write();
		h1d_truetrklen_CaloSz_inel->Write();
		h1d_truetrklen_Pos_inel->Write();
		h1d_truetrklen_BQ_inel->Write();
		h1d_truetrklen_RecoInel_inel->Write();
		h1d_truetrklen_BeamMatch_inel->Write();

		h1d_eff_NoCut_inel->Write();
		h1d_eff_Pan_inel->Write();
		h1d_eff_CaloSz_inel->Write();
		h1d_eff_Pos_inel->Write();
		h1d_eff_BQ_inel->Write();
		h1d_eff_RecoInel_inel->Write();
		h1d_eff_BeamMatch_inel->Write();

		Eff_NoCut_inel->Write();
		Eff_Pan_inel->Write();
		Eff_CaloSz_inel->Write();
		Eff_Pos_inel->Write();
		Eff_BQ_inel->Write();
		Eff_RecoInel_inel->Write();
		Eff_BeamMatch_inel->Write();

		h2d_thetaXZ_thetaYZ_BQ_inel->Write();
		h2d_thetaXZ_thetaYZ_BQ_misidp->Write();
		h2d_thetaXZ_thetaYZ_BQ_el->Write();

		h2d_thetaXZ_thetaYZ_Pos_inel->Write();
		h2d_thetaXZ_thetaYZ_Pos_misidp->Write();
		h2d_thetaXZ_thetaYZ_Pos_el->Write();

		h2d_truetrklen_trueendz_CaloSz_inel->Write();
		h1d_trueendz_CaloSz_inel->Write();



		h2d_truetrklen_ketrue_NoCut_inel->Write();
		h2d_truetrklen_keff_NoCut_inel->Write();

		h2d_truetrklen_ketrue_Pan_inel->Write();
		h2d_truetrklen_ketrue_CaloSz_inel->Write();


		h2d_trueEndZ_ketrue_NoCut_inel->Write();
		h2d_trueZ_ketrue_NoCut_inel->Write();

		//h2d_trueXY_truelen05_NoCut_inel->Write();
		//h2d_trueXY_truelen510_NoCut_inel->Write();
		//h2d_trueXY_truelen1015_NoCut_inel->Write();
		//h2d_trueXY_truelen1520_NoCut_inel->Write();
		//h2d_trueXY_truelen2025_NoCut_inel->Write();
		//h2d_trueXY_truelen2530_NoCut_inel->Write();
		//h2d_trueXY_truelen3035_NoCut_inel->Write();
		//h3d_trueXYZ_NoCut_inel->Write();

		//for (int j=0; j<n_scan; ++j) { 
			//h2d_trueXY_NoCut_inel[j]->Write();
		//}
	fout->Close();


}
