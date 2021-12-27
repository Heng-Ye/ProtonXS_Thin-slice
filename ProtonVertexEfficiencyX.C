#define ProtonVertexEfficiency_cxx
#include "ProtonVertexEfficiency.h"

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

void ProtonVertexEfficiency::Loop() {
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	Long64_t nbytes = 0, nb = 0;
	bool isTestSample=true;
	//int true_sliceID = -1, reco_sliceID = -1;


	//book histograms --------------------------------------------------------------------------------------------//
	int n_b=300;
	double b_min=0;
	double b_max=150;

	int n_dr=600;
	double dr_min=-150;
	double dr_max=150;

	int n_dcos=100;
	double dcos_min=0;
	double dcos_max=1;

	

	//pans
	TH2D *h2d_trklen_dr_PanS=new TH2D(Form("h2d_trklen_dr_PanS"), Form("h2d_trklen_dr_PanS"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_PanS_el=new TH2D(Form("h2d_trklen_dr_PanS_el"), Form("h2d_trklen_dr_PanS_el"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_PanS_inel=new TH2D(Form("h2d_trklen_dr_PanS_inel"), Form("h2d_trklen_dr_PanS_inel"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_PanS_misidp=new TH2D(Form("h2d_trklen_dr_PanS_misidp"), Form("h2d_trklen_dr_PanS_misidp"), n_b, b_min, b_max, n_dr, dr_min, dr_max);

	//calosz
	TH2D *h2d_trklen_dr_CaloSz=new TH2D(Form("h2d_trklen_dr_CaloSz"), Form("h2d_trklen_dr_CaloSz"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_CaloSz_el=new TH2D(Form("h2d_trklen_dr_CaloSz_el"), Form("h2d_trklen_dr_CaloSz_el"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_CaloSz_inel=new TH2D(Form("h2d_trklen_dr_CaloSz_inel"), Form("h2d_trklen_dr_CaloSz_inel"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_CaloSz_misidp=new TH2D(Form("h2d_trklen_dr_CaloSz_misidp"), Form("h2d_trklen_dr_CaloSz_misidp"), n_b, b_min, b_max, n_dr, dr_min, dr_max);

	//pos
	TH2D *h2d_trklen_dr_Pos=new TH2D(Form("h2d_trklen_dr_Pos"), Form("h2d_trklen_dr_Pos"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_Pos_el=new TH2D(Form("h2d_trklen_dr_Pos_el"), Form("h2d_trklen_dr_Pos_el"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_Pos_inel=new TH2D(Form("h2d_trklen_dr_Pos_inel"), Form("h2d_trklen_dr_Pos_inel"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_Pos_misidp=new TH2D(Form("h2d_trklen_dr_Pos_misidp"), Form("h2d_trklen_dr_Pos_misidp"), n_b, b_min, b_max, n_dr, dr_min, dr_max);

	//bq
	TH2D *h2d_trklen_dr_BQ=new TH2D(Form("h2d_trklen_dr_BQ"), Form("h2d_trklen_dr_BQ"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_BQ_el=new TH2D(Form("h2d_trklen_dr_BQ_el"), Form("h2d_trklen_dr_BQ_el"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_BQ_inel=new TH2D(Form("h2d_trklen_dr_BQ_inel"), Form("h2d_trklen_dr_BQ_inel"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_BQ_misidp=new TH2D(Form("h2d_trklen_dr_BQ_misidp"), Form("h2d_trklen_dr_BQ_misidp"), n_b, b_min, b_max, n_dr, dr_min, dr_max);

	//reco inel
	TH2D *h2d_trklen_dr_RecoInel=new TH2D(Form("h2d_trklen_dr_RecoInel"), Form("h2d_trklen_dr_RecoInel"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_RecoInel_el=new TH2D(Form("h2d_trklen_dr_RecoInel_el"), Form("h2d_trklen_dr_RecoInel_el"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_RecoInel_inel=new TH2D(Form("h2d_trklen_dr_RecoInel_inel"), Form("h2d_trklen_dr_RecoInel_inel"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_RecoInel_2060_inel=new TH2D(Form("h2d_trklen_dr_RecoInel_2060_inel"), Form("h2d_trklen_dr_RecoInel_2060_inel"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_RecoInel_0010_inel=new TH2D(Form("h2d_trklen_dr_RecoInel_0010_inel"), Form("h2d_trklen_dr_RecoInel_0010_inel"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_RecoInel_misidp=new TH2D(Form("h2d_trklen_dr_RecoInel_misidp"), Form("h2d_trklen_dr_RecoInel_misidp"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_truetrklen_dr_RecoInel_inel=new TH2D(Form("h2d_truetrklen_dr_RecoInel_inel"), Form("h2d_truetrklen_dr_RecoInel_inel"), n_b, b_min, b_max, n_dr, dr_min, dr_max);

	TH2D *h2d_trklen_dr_RecoInel_GOOD=new TH2D(Form("h2d_trklen_dr_RecoInel_GOOD"), Form("h2d_trklen_dr_RecoInel_GOOD"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_RecoInel_BAD=new TH2D(Form("h2d_trklen_dr_RecoInel_BAD"), Form("h2d_trklen_dr_RecoInel_BAD"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_RecoInel_UGLY=new TH2D(Form("h2d_trklen_dr_RecoInel_UGLY"), Form("h2d_trklen_dr_RecoInel_UGLY"), n_b, b_min, b_max, n_dr, dr_min, dr_max);

	//reco el
	TH2D *h2d_trklen_dr_RecoEl=new TH2D(Form("h2d_trklen_dr_RecoEl"), Form("h2d_trklen_dr_RecoEl"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_RecoEl_el=new TH2D(Form("h2d_trklen_dr_RecoEl_el"), Form("h2d_trklen_dr_RecoEl_el"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_RecoEl_inel=new TH2D(Form("h2d_trklen_dr_RecoEl_inel"), Form("h2d_trklen_dr_RecoEl_inel"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_RecoEl_misidp=new TH2D(Form("h2d_trklen_dr_RecoEl_misidp"), Form("h2d_trklen_dr_RecoEl_misidp"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_truetrklen_dr_RecoEl_el=new TH2D(Form("h2d_truetrklen_dr_RecoEl_el"), Form("h2d_truetrklen_dr_RecoEl_el"), n_b, b_min, b_max, n_dr, dr_min, dr_max);


	//RecoInel+MisID:P-rich
	TH2D *h2d_trklen_dr_RecoInelMIDP=new TH2D(Form("h2d_trklen_dr_RecoInelMIDP"), Form("h2d_trklen_dr_RecoInelMIDP"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_RecoInelMIDP_el=new TH2D(Form("h2d_trklen_dr_RecoInelMIDP_el"), Form("h2d_trklen_dr_RecoInelMIDP_el"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_RecoInelMIDP_inel=new TH2D(Form("h2d_trklen_dr_RecoInelMIDP_inel"), Form("h2d_trklen_dr_RecoInelMIDP_inel"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_trklen_dr_RecoInelMIDP_misidp=new TH2D(Form("h2d_trklen_dr_RecoInelMIDP_misidp"), Form("h2d_trklen_dr_RecoInelMIDP_misidp"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_truetrklen_dr_RecoInelMIDP_inel=new TH2D(Form("h2d_truetrklen_dr_RecoInelMIDP_inel"), Form("h2d_truetrklen_dr_RecoInelMIDP_inel"), n_b, b_min, b_max, n_dr, dr_min, dr_max);
	TH2D *h2d_truetrklen_dr_RecoInelMIDP_misidp=new TH2D(Form("h2d_truetrklen_dr_RecoInelMIDP_misidp"), Form("h2d_truetrklen_dr_RecoInelMIDP_misidp"), n_b, b_min, b_max, n_dr, dr_min, dr_max);



	//cos_vs_dr pos
	int n_cos=100;
	double cos_min=0.;
	double cos_max=1.;
	TH2D *h2d_cos_dr_Pos=new TH2D(Form("h2d_cos_dr_Pos"), Form("h2d_cos_dr_Pos"), n_cos, cos_min, cos_max, n_dr, dr_min, dr_max);
	TH2D *h2d_cos_dr_Pos_el=new TH2D(Form("h2d_cos_dr_Pos_el"), Form("h2d_cos_dr_Pos_el"), n_cos, cos_min, cos_max, n_dr, dr_min, dr_max);
	TH2D *h2d_cos_dr_Pos_inel=new TH2D(Form("h2d_cos_dr_Pos_inel"), Form("h2d_cos_dr_Pos_inel"), n_cos, cos_min, cos_max, n_dr, dr_min, dr_max);
	TH2D *h2d_cos_dr_Pos_misidp=new TH2D(Form("h2d_cos_dr_Pos_misidp"), Form("h2d_cos_dr_Pos_misidp"), n_cos, cos_min, cos_max, n_dr, dr_min, dr_max);

	TH2D *h2d_cos_dr_PosRecoInel=new TH2D(Form("h2d_cos_dr_PosRecoInel"), Form("h2d_cos_dr_PosRecoInel"), n_cos, cos_min, cos_max, n_dr, dr_min, dr_max);
	TH2D *h2d_cos_dr_PosRecoInel_el=new TH2D(Form("h2d_cos_dr_PosRecoInel_el"), Form("h2d_cos_dr_PosRecoInel_el"), n_cos, cos_min, cos_max, n_dr, dr_min, dr_max);
	TH2D *h2d_cos_dr_PosRecoInel_inel=new TH2D(Form("h2d_cos_dr_PosRecoInel_inel"), Form("h2d_cos_dr_PosRecoInel_inel"), n_cos, cos_min, cos_max, n_dr, dr_min, dr_max);
	TH2D *h2d_cos_dr_PosRecoInel_misidp=new TH2D(Form("h2d_cos_dr_PosRecoInel_misidp"), Form("h2d_cos_dr_PosRecoInel_misidp"), n_cos, cos_min, cos_max, n_dr, dr_min, dr_max);

	TH2D *h2d_trklen_cos_PosRecoInel=new TH2D(Form("h2d_trklen_cos_PosRecoInel"), Form("h2d_trklen_cos_PosRecoInel"), n_b, b_min, b_max, n_cos, cos_min, cos_max);
	TH2D *h2d_trklen_cos_PosRecoInel_el=new TH2D(Form("h2d_trklen_cos_PosRecoInel_el"), Form("h2d_trklen_cos_PosRecoInel_el"), n_b, b_min, b_max, n_cos, cos_min, cos_max);
	TH2D *h2d_trklen_cos_PosRecoInel_inel=new TH2D(Form("h2d_trklen_cos_PosRecoInel_inel"), Form("h2d_trklen_cos_PosRecoInel_inel"), n_b, b_min, b_max, n_cos, cos_min, cos_max);
	TH2D *h2d_trklen_cos_PosRecoInel_misidp=new TH2D(Form("h2d_trklen_cos_PosRecoInel_misidp"), Form("h2d_trklen_cos_PosRecoInel_misidp"), n_b, b_min, b_max, n_cos, cos_min, cos_max);

	//Edept_per_recolen vs ntrklen
	int n_avedept=120;
	double avedept_min=0;
	double avedept_max=12;
	int n_ntrklen=150.;
	double ntrklen_min=0.;
	double ntrklen_max=1.5;
	//bq
	TH2D *h2d_ntrklen_avedept_BQ=new TH2D(Form("h2d_ntrklen_avedept_BQ"), Form("h2d_ntrklen_avedept_BQ"), n_ntrklen, ntrklen_min, ntrklen_max, n_avedept, avedept_min, avedept_max);
	TH2D *h2d_ntrklen_avedept_BQ_inel=new TH2D(Form("h2d_ntrklen_avedept_BQ_inel"), Form("h2d_ntrklen_avedept_BQ_inel"), n_ntrklen, ntrklen_min, ntrklen_max, n_avedept, avedept_min, avedept_max);
	TH2D *h2d_ntrklen_avedept_BQ_el=new TH2D(Form("h2d_ntrklen_avedept_BQ_el"), Form("h2d_ntrklen_avedept_BQ_el"), n_ntrklen, ntrklen_min, ntrklen_max, n_avedept, avedept_min, avedept_max);
	TH2D *h2d_ntrklen_avedept_BQ_misidp=new TH2D(Form("h2d_ntrklen_avedept_BQ_misidp"), Form("h2d_ntrklen_avedept_BQ_misidp"), n_ntrklen, ntrklen_min, ntrklen_max, n_avedept, avedept_min, avedept_max);
	//pos
	TH2D *h2d_ntrklen_avedept_Pos=new TH2D(Form("h2d_ntrklen_avedept_Pos"), Form("h2d_ntrklen_avedept_Pos"), n_ntrklen, ntrklen_min, ntrklen_max, n_avedept, avedept_min, avedept_max);
	TH2D *h2d_ntrklen_avedept_Pos_inel=new TH2D(Form("h2d_ntrklen_avedept_Pos_inel"), Form("h2d_ntrklen_avedept_Pos_inel"), n_ntrklen, ntrklen_min, ntrklen_max, n_avedept, avedept_min, avedept_max);
	TH2D *h2d_ntrklen_avedept_Pos_el=new TH2D(Form("h2d_ntrklen_avedept_Pos_el"), Form("h2d_ntrklen_avedept_Pos_el"), n_ntrklen, ntrklen_min, ntrklen_max, n_avedept, avedept_min, avedept_max);
	TH2D *h2d_ntrklen_avedept_Pos_misidp=new TH2D(Form("h2d_ntrklen_avedept_Pos_misidp"), Form("h2d_ntrklen_avedept_Pos_misidp"), n_ntrklen, ntrklen_min, ntrklen_max, n_avedept, avedept_min, avedept_max);


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

		bool IsMIDP=false;
		if (cosine_beam_spec_primtrk<0.9) IsMIDP=true;
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
		double range_reco=99999;
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
		//cout<<"haha"<<endl;
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
		//cout<<"haha1"<<endl;

		if (is_beam_at_ff) { //if primary protons entering tpc
			for (int iz=key_reach_tpc+1; iz<(int)beamtrk_z->size(); iz++){
				if (iz == key_reach_tpc+1) range_true = 0;
					range_true += sqrt( pow(beamtrk_x->at(iz)-beamtrk_x->at(iz-1), 2)+
					pow(beamtrk_y->at(iz)-beamtrk_y->at(iz-1), 2)+	
					pow(beamtrk_z->at(iz)-beamtrk_z->at(iz-1), 2) );						    	
					true_trklen_accum[iz] = range_true;
			}
		} //if primary protons entering tpc

		//cout<<"haha2"<<endl;
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
	
		//range compensation ------------------------------------------------------------//
		//cout<<"haha3"<<endl;
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
		///double dr=sqrt(pow(reco_endx-true_endx,2)+pow(reco_endy-true_endy,2)+pow(reco_endz-true_endz,2));

		double dr=range_reco-range_true;
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

		double tmp_dr=dr;
		if (tmp_dr>dr_max) tmp_dr=dr_max;
		if (tmp_dr<dr_min) tmp_dr=dr_min;	

		//double tmp_trklen=range_true;
		double tmp_trklen=range_reco;
		if (tmp_trklen<b_min) tmp_trklen=b_min;
		if (tmp_trklen>b_max) tmp_trklen=b_max;

		double tmp_cos=cosine_beam_spec_primtrk;
		if (tmp_cos<cos_min) tmp_cos=cos_min;
		if (tmp_cos>cos_max) tmp_cos=cos_max;

		double tmp_ntrklen=range_reco/csda_val_spec;
		if (tmp_ntrklen<ntrklen_min) tmp_ntrklen=ntrklen_min;
		if (tmp_ntrklen>ntrklen_max) tmp_ntrklen=ntrklen_max;
		
		double tmp_avedept=kereco_calo/range_reco;
		if (tmp_avedept<avedept_min) tmp_avedept=avedept_min;
		if (tmp_avedept>avedept_max) tmp_avedept=avedept_max;




		if (IsPandoraSlice) { //PanS
			h2d_trklen_dr_PanS->Fill(tmp_trklen, tmp_dr);
			if (kinel) {
				h2d_trklen_dr_PanS_inel->Fill(tmp_trklen, tmp_dr);
			}
			if (kel) {
				h2d_trklen_dr_PanS_el->Fill(tmp_trklen, tmp_dr);
			}
			if (kMIDp) {
				h2d_trklen_dr_PanS_misidp->Fill(tmp_trklen, tmp_dr);
			}
		} //PanS

		if (IsCaloSize&&IsPandoraSlice) { //CaloSz
			h2d_trklen_dr_CaloSz->Fill(tmp_trklen, tmp_dr);
			if (kinel) {
				h2d_trklen_dr_CaloSz_inel->Fill(tmp_trklen, tmp_dr);
			}
			if (kel) {
				h2d_trklen_dr_CaloSz_el->Fill(tmp_trklen, tmp_dr);
			}
			if (kMIDp) {
				h2d_trklen_dr_CaloSz_misidp->Fill(tmp_trklen, tmp_dr);
			}
		} //CaloSz
	
		if (IsPos&&IsCaloSize&&IsPandoraSlice) { //Pos
			h2d_trklen_dr_Pos->Fill(tmp_trklen, tmp_dr);
			h2d_cos_dr_Pos->Fill(tmp_cos, tmp_dr);
			h2d_ntrklen_avedept_Pos->Fill(tmp_ntrklen, tmp_avedept);
			if (kinel) {
				h2d_trklen_dr_Pos_inel->Fill(tmp_trklen, tmp_dr);
				h2d_cos_dr_Pos_inel->Fill(tmp_cos, tmp_dr);
				h2d_ntrklen_avedept_Pos_inel->Fill(tmp_ntrklen, tmp_avedept);
			}
			if (kel) {
				h2d_trklen_dr_Pos_el->Fill(tmp_trklen, tmp_dr);
				h2d_cos_dr_Pos_el->Fill(tmp_cos, tmp_dr);
				h2d_ntrklen_avedept_Pos_el->Fill(tmp_ntrklen, tmp_avedept);
			}
			if (kMIDp) {
				h2d_trklen_dr_Pos_misidp->Fill(tmp_trklen, tmp_dr);
				h2d_cos_dr_Pos_misidp->Fill(tmp_cos, tmp_dr);
				h2d_ntrklen_avedept_Pos_misidp->Fill(tmp_ntrklen, tmp_avedept);
			}
		} //Pos


		if (IsBQ&&IsCaloSize&&IsPandoraSlice) { //BQ
			h2d_trklen_dr_BQ->Fill(tmp_trklen, tmp_dr);
			h2d_ntrklen_avedept_BQ->Fill(tmp_ntrklen, tmp_avedept);
			if (kinel) {
				h2d_trklen_dr_BQ_inel->Fill(tmp_trklen, tmp_dr);
				h2d_ntrklen_avedept_BQ_inel->Fill(tmp_ntrklen, tmp_avedept);
			}
			if (kel) {
				h2d_trklen_dr_BQ_el->Fill(tmp_trklen, tmp_dr);
				h2d_ntrklen_avedept_BQ_el->Fill(tmp_ntrklen, tmp_avedept);
			}
			if (kMIDp) {
				h2d_trklen_dr_BQ_misidp->Fill(tmp_trklen, tmp_dr);
				h2d_ntrklen_avedept_BQ_misidp->Fill(tmp_ntrklen, tmp_avedept);
			}
		} //BQ

		if (IsRecoInEL&&IsBQ&&IsCaloSize&&IsPandoraSlice) { //RecoInel
			h2d_trklen_dr_RecoInel->Fill(tmp_trklen, tmp_dr);

			if (kinel) {
				h2d_trklen_dr_RecoInel_inel->Fill(tmp_trklen, tmp_dr);
				h2d_truetrklen_dr_RecoInel_inel->Fill(range_true, tmp_dr);
				if (tmp_trklen>=20&&tmp_trklen<=60) h2d_trklen_dr_RecoInel_2060_inel->Fill(tmp_trklen, tmp_dr);
				if (tmp_trklen>=0&&tmp_trklen<=10) h2d_trklen_dr_RecoInel_0010_inel->Fill(tmp_trklen, tmp_dr);

				if (tmp_dr>=-3&&tmp_dr<=3) h2d_trklen_dr_RecoInel_GOOD->Fill(tmp_trklen, tmp_dr);
				if (tmp_dr<-3) h2d_trklen_dr_RecoInel_BAD->Fill(tmp_trklen, tmp_dr);
				if (tmp_dr>3) h2d_trklen_dr_RecoInel_UGLY->Fill(tmp_trklen, tmp_dr);
			}
			if (kel) {
				h2d_trklen_dr_RecoInel_el->Fill(tmp_trklen, tmp_dr);
			}
			if (kMIDp) {
				h2d_trklen_dr_RecoInel_misidp->Fill(tmp_trklen, tmp_dr);
			}
		} //RecoInel


		if (IsRecoStop&&IsBQ&&IsCaloSize&&IsPandoraSlice) { //RecoEl
			h2d_trklen_dr_RecoEl->Fill(tmp_trklen, tmp_dr);

			if (kinel) {
				h2d_trklen_dr_RecoEl_inel->Fill(tmp_trklen, tmp_dr);
				h2d_truetrklen_dr_RecoEl_el->Fill(range_true, tmp_dr);
			}
			if (kel) {
				h2d_trklen_dr_RecoEl_el->Fill(tmp_trklen, tmp_dr);
			}
			if (kMIDp) {
				h2d_trklen_dr_RecoEl_misidp->Fill(tmp_trklen, tmp_dr);
			}
		} //RecoEl

		if (IsRecoInEL&&IsPos&&IsCaloSize&&IsPandoraSlice) { //Pos+RecoInel
			h2d_cos_dr_PosRecoInel->Fill(tmp_cos, tmp_dr);
			h2d_trklen_cos_PosRecoInel->Fill(tmp_trklen, tmp_cos);
			if (kinel) {
				h2d_cos_dr_PosRecoInel_inel->Fill(tmp_cos, tmp_dr);
				h2d_trklen_cos_PosRecoInel_inel->Fill(tmp_trklen, tmp_cos);
			}
			if (kel) {
				h2d_cos_dr_PosRecoInel_el->Fill(tmp_cos, tmp_dr);
				h2d_trklen_cos_PosRecoInel_el->Fill(tmp_trklen, tmp_cos);
			}
			if (kMIDp) {
				h2d_cos_dr_PosRecoInel_misidp->Fill(tmp_cos, tmp_dr);
				h2d_trklen_cos_PosRecoInel_misidp->Fill(tmp_trklen, tmp_cos);
			}
		} //Pos+RecoInel


		if (IsRecoInEL&&IsMIDP&&IsPos&&IsCaloSize&&IsPandoraSlice) { //Pos+RecoInel+MidP
			h2d_trklen_dr_RecoInelMIDP->Fill(tmp_trklen, tmp_dr);
			if (kinel) {
				h2d_trklen_dr_RecoInelMIDP_inel->Fill(tmp_trklen, tmp_dr);
				h2d_truetrklen_dr_RecoInelMIDP_inel->Fill(range_true, tmp_dr);
			}
			if (kel) {
				h2d_trklen_dr_RecoInelMIDP_el->Fill(tmp_trklen, tmp_dr);
			}
			if (kMIDp) {
				h2d_trklen_dr_RecoInelMIDP_misidp->Fill(tmp_trklen, tmp_dr);
				h2d_truetrklen_dr_RecoInelMIDP_misidp->Fill(range_true, tmp_dr);
			}
		} //Pos+RecoInel+MidP






	} //main entry loop


	//save results ---------------------------------------------------------//
   	//TFile *fout = new TFile("mc_vtx.root","RECREATE");
   	TFile *fout = new TFile("mc_vtx_rangereco.root","RECREATE");

		h2d_trklen_dr_PanS->Write();
		h2d_trklen_dr_PanS_el->Write();
		h2d_trklen_dr_PanS_inel->Write();
		h2d_trklen_dr_PanS_misidp->Write();

		h2d_trklen_dr_CaloSz->Write();
		h2d_trklen_dr_CaloSz_el->Write();
		h2d_trklen_dr_CaloSz_inel->Write();
		h2d_trklen_dr_CaloSz_misidp->Write();

		h2d_trklen_dr_Pos->Write();
		h2d_trklen_dr_Pos_el->Write();
		h2d_trklen_dr_Pos_inel->Write();
		h2d_trklen_dr_Pos_misidp->Write();

		h2d_trklen_dr_BQ->Write();
		h2d_trklen_dr_BQ_el->Write();
		h2d_trklen_dr_BQ_inel->Write();
		h2d_trklen_dr_BQ_misidp->Write();

		h2d_trklen_dr_RecoInel->Write();
		h2d_trklen_dr_RecoInel_el->Write();
		h2d_trklen_dr_RecoInel_inel->Write();
		h2d_trklen_dr_RecoInel_misidp->Write();

		h2d_trklen_dr_RecoInel_2060_inel->Write();
		h2d_trklen_dr_RecoInel_0010_inel->Write();

		h2d_trklen_dr_RecoInel_BAD->Write();
		h2d_trklen_dr_RecoInel_UGLY->Write();
		h2d_trklen_dr_RecoInel_GOOD->Write();

		h2d_truetrklen_dr_RecoInel_inel->Write();


		h2d_trklen_dr_RecoEl->Write();
		h2d_trklen_dr_RecoEl_el->Write();
		h2d_trklen_dr_RecoEl_inel->Write();
		h2d_trklen_dr_RecoEl_misidp->Write();
		h2d_truetrklen_dr_RecoEl_el->Write();



		h2d_cos_dr_Pos->Write();
		h2d_cos_dr_Pos_el->Write();
		h2d_cos_dr_Pos_inel->Write();
		h2d_cos_dr_Pos_misidp->Write();

		h2d_ntrklen_avedept_Pos->Write();
		h2d_ntrklen_avedept_Pos_inel->Write();
		h2d_ntrklen_avedept_Pos_el->Write();
		h2d_ntrklen_avedept_Pos_misidp->Write();

		h2d_ntrklen_avedept_BQ->Write();
		h2d_ntrklen_avedept_BQ_inel->Write();
		h2d_ntrklen_avedept_BQ_el->Write();
		h2d_ntrklen_avedept_BQ_misidp->Write();

		h2d_cos_dr_PosRecoInel->Write();
		h2d_cos_dr_PosRecoInel_inel->Write();	
		h2d_cos_dr_PosRecoInel_el->Write();	
		h2d_cos_dr_PosRecoInel_misidp->Write();

		h2d_trklen_cos_PosRecoInel->Write();
		h2d_trklen_cos_PosRecoInel_inel->Write();
		h2d_trklen_cos_PosRecoInel_el->Write();
		h2d_trklen_cos_PosRecoInel_misidp->Write();


		h2d_trklen_dr_RecoInelMIDP->Write();
		h2d_trklen_dr_RecoInelMIDP_el->Write();
		h2d_trklen_dr_RecoInelMIDP_inel->Write();
		h2d_trklen_dr_RecoInelMIDP_misidp->Write();
		h2d_truetrklen_dr_RecoInelMIDP_inel->Write();
		h2d_truetrklen_dr_RecoInelMIDP_misidp->Write();

	fout->Close();


}
