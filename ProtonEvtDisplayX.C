#define ProtonEvtDisplay_cxx
#include "ProtonEvtDisplay.h"

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
#include "TVector3.h"

#include "./cali/dedx_function_35ms.h"
#include "./headers/BasicParameters.h"
#include "./headers/BasicFunctions.h"
#include "./headers/SliceParams.h"
#include "./headers/util.h"
//#include "./headers/ThinSlice.h"
//#include "./headers/BackgroundStudy.h"

//#include "./headers/ThinSlice.h"

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


void ProtonEvtDisplay::Loop() {
	if (fChain == 0) return;

	//various counters ---------------------------------------------------------------------------------------------------------------------------------------------------//
	int n_processmap_error=0; //sansity check: processmap
	int n_true_end=0; //sansity check: true-end label

	int n_test_recoinel=0, n_test_recoinel_sample=0, n_kinel=0, n_kinel2=0, n_kel=0;

	int n_tot=0, n_sg=0, n_bkg=0; //no cut
	int n_el=0, n_inel=0, n_midcosmic=0, n_midpi=0, n_midp=0, n_midmu=0, n_mideg=0, n_midother=0;

	int n_pan_tot=0; //pandora cut
	int n_el_pan=0, n_inel_pan=0, n_midcosmic_pan=0, n_midpi_pan=0, n_midp_pan=0, n_midmu_pan=0, n_mideg_pan=0, n_midother_pan=0;

	int n_calsz_tot=0; //calosz cut
	int n_el_calsz=0, n_inel_calsz=0, n_midcosmic_calsz=0, n_midpi_calsz=0, n_midp_calsz=0, n_midmu_calsz=0, n_mideg_calsz=0, n_midother_calsz=0;

	int n_bq_tot=0; //bq cut
	int n_el_bq=0, n_inel_bq=0, n_midcosmic_bq=0, n_midpi_bq=0, n_midp_bq=0, n_midmu_bq=0, n_mideg_bq=0, n_midother_bq=0;

	int n_recoinel_tot=0; //recoinel cut
	int n_el_recoinel=0, n_inel_recoinel=0, n_midcosmic_recoinel=0, n_midpi_recoinel=0, n_midp_recoinel=0, n_midmu_recoinel=0, n_mideg_recoinel=0, n_midother_recoinel=0;

	//unfolding config. -------------------------//
	//Unfold uf(nthinslices+2, -1, nthinslices+1);

	//ThinSlice config. ---------------------------------------------------------------------------------------------------//
	//SetOutputFileName(Form("prod4a_bkgstudy_dx%dcm_%dslcs_largerbin.root", name_thinslicewidth, nthinslices)); //output file name
	//SetOutputFileName(Form("prod4a_bkgstudy_dx%dcm_%dslcs.root", name_thinslicewidth, nthinslices)); //output file name
	//SetOutputFileName(Form("evt_Demp_%dcm_%dslcs.root", name_thinslicewidth, nthinslices)); //output file name

	//book histograms --//
	//BookHistograms();

	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	Long64_t nbytes = 0, nb = 0;
	bool isTestSample=true;
	int true_sliceID = -1, reco_sliceID = -1;
	for (Long64_t jentry=0; jentry<nentries;jentry++) { //main entry loop
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		isTestSample = true;
		if (ientry%2 == 0) isTestSample = false; //Divide MC sample by 2 parts: test+ufold

		true_sliceID = -1;
		reco_sliceID = -1;

		//only select protons	
		//if (primary_truth_Pdg!=pdg) continue; //only interested in protons
		if (beamtrackPdg!=pdg) continue; //only interested in protons
		//std::cout<<"beamtrackPdg:"<<beamtrackPdg<<std::endl;
		n_tot++;

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
		//bool IsPureMCS=false; //no hadron scattering

		if (primary_truth_EndProcess->c_str()!=NULL) n_true_end++;
		if (strcmp(primary_truth_EndProcess->c_str(),"protonInelastic")==0) {
			IsPureInEL=true;
		}
		else { //hIoni
			IsPureEL=true;
		} //hIoni
		//--------------------------------------------------------------------------------------------------------------------------------//

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

		double KE_ff=0;
		if (is_beam_at_ff) KE_ff=1000.*beamtrk_Eng->at(key_reach_tpc); //unit:MeV
		//KE_ff=ke_ff; //use KE exactly at z=0
		//---------------------------------------------------------------------------------------------------------------//

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


		//[2] Range compensation -------------------------------------------------------
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

			//old position cut
			//if (reco_stz>min1_z&&reco_stz<min2_z) {
			//if (reco_sty>min1_y&&reco_sty<min2_y) {
			//if (reco_stx>min1_x&&reco_stx<min2_x) {
			//IsPos=true;
			//}
			//}
			//}
		}

		//cosine_theta/cut ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
		bool IsCosine=false;
		bool IsFlip=false;
		bool IsMisIDP=false;
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

			if (cosine_beam_spec_primtrk<0) { 
				cosine_beam_spec_primtrk=-1.*cosine_beam_spec_primtrk; 
				IsFlip=true;	
			}
			//if (cosine_beam_spec_primtrk>cosine_beam_primtrk_min) { IsCosine=true; }
			if (cosine_beam_spec_primtrk>costh_min&&cosine_beam_spec_primtrk<costh_max) { IsCosine=true; }
			if (cosine_beam_spec_primtrk<0.9) { IsMisIDP=true; }
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

		//if (IsPos&&IsCaloSize&&IsPandoraSlice) {
			//Fill1DHist(reco_cosineTheta_Pos, cosine_beam_spec_primtrk);
			//if (kinel) Fill1DHist(reco_cosineTheta_Pos_inel, cosine_beam_spec_primtrk);
			//if (kel) Fill1DHist(reco_cosineTheta_Pos_el, cosine_beam_spec_primtrk);
			//if (kMIDcosmic) Fill1DHist(reco_cosineTheta_Pos_midcosmic, cosine_beam_spec_primtrk);
			//if (kMIDpi) Fill1DHist(reco_cosineTheta_Pos_midpi, cosine_beam_spec_primtrk);
			//if (kMIDp) Fill1DHist(reco_cosineTheta_Pos_midp, cosine_beam_spec_primtrk);
			//if (kMIDmu) Fill1DHist(reco_cosineTheta_Pos_midmu, cosine_beam_spec_primtrk);
			//if (kMIDeg) Fill1DHist(reco_cosineTheta_Pos_mideg, cosine_beam_spec_primtrk);
			//if (kMIDother) Fill1DHist(reco_cosineTheta_Pos_midother, cosine_beam_spec_primtrk);
		//}


		//reco calorimetry ---------------------------------------------------------------------------//
		//vector< pair<double,int > > zreco_rawindex; //z, original_index
		//vector< pair<double,double > > zreco_rr; //z, rr
		//vector< pair<double,double > > zreco_de; //z, wid_reco
		//vector< pair<double,double > > zreco_dedx; //z, wid_reco
		//vector< pair<double,double > > zreco_dx; //z, wid_reco
		//vector< pair<double,double > > zreco_xreco; //z, wid_reco
		//vector< pair<double,double > > zreco_yreco; //z, wid_reco
		//vector< pair<double,int > > zreco_widreco; //z, wid_reco
		//vector< pair<double,double > > zreco_lenreco; //z, len_reco

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

				//if (wid_reco==-9999) continue; //outside TPC
				if (wid_reco>wid_reco_max) { 
					wid_reco_max=wid_reco;
					index_reco_endz=(int)-1+primtrk_wid->size()-h;
				}

				double cali_dedx=0.;
				cali_dedx=dedx_function_35ms(dqdx, hitx_reco, hity_reco, hitz_reco);

				EDept.push_back(cali_dedx*pitch);

				//zreco_rawindex.push_back(make_pair(wid_reco, h));
				//zreco_de.push_back(make_pair(wid_reco, cali_dedx*pitch));
				//zreco_dedx.push_back(make_pair(wid_reco, cali_dedx));
				//zreco_dx.push_back(make_pair(wid_reco, pitch));
				//zreco_xreco.push_back(make_pair(wid_reco, hitx_reco));
				//zreco_yreco.push_back(make_pair(wid_reco, hity_reco));
				//zreco_widreco.push_back(make_pair(wid_reco, hitz_reco));
				//zreco_rr.push_back(make_pair(wid_reco, resrange_reco));

				//if (IsPureInEL) rangereco_dedxreco_TrueInEL->Fill(range_reco-resrange_reco, cali_dedx);
				//if (IsPureEL) rangereco_dedxreco_TrueEL->Fill(range_reco-resrange_reco, cali_dedx);
				//if (IsPureMCS) rangereco_dedxreco_TrueMCS->Fill(range_reco-resrange_reco, cali_dedx);

				//zreco_ke.push_back(make_pair(wid_reco, ke_reco));

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

				trkdedx.push_back(cali_dedx);
				trkres.push_back(resrange_reco);

			} //loop over reco hits of a given track

			pid=chi2pid(trkdedx,trkres); //pid using stopping proton hypothesis

		} //if calo size not empty

		//use wid for sorting
		//sort(zreco_rawindex.begin(),zreco_rawindex.end(),myComparison); //sorting based on the first column
		//sort(zreco_rr.begin(),zreco_rr.end(),myComparison); //sorting based on the first column
		//sort(zreco_de.begin(),zreco_de.end(),myComparison); //sorting based on the first column
		//sort(zreco_dedx.begin(),zreco_dedx.end(),myComparison); //sorting based on the first column
		//sort(zreco_dx.begin(),zreco_dx.end(),myComparison); //sorting based on the first column
		//sort(zreco_xreco.begin(),zreco_xreco.end(),myComparison); //sorting based on the first column
		//sort(zreco_yreco.begin(),zreco_yreco.end(),myComparison); //sorting based on the first column
		//sort(zreco_widreco.begin(),zreco_widreco.end(),myComparison); //sorting based on the first column
		//sort(zreco_ke.begin(),zreco_ke.end(),myComparison); //sorting based on the first column

		//reco trklen ---------------------------------------------------------------------------//
		//double range_reco=-999;
		//vector<double> reco_trklen_accum;
		//reco_trklen_accum.reserve(primtrk_hitz->size());
		//for (int iz=1; iz<(int)zreco_rawindex.size(); iz++) {
		//if (iz==1) range_reco=0; 	
		//range_reco += sqrt( pow(zreco_xreco[iz].second-zreco_xreco[iz-1].second, 2)+
		//pow(zreco_yreco[iz].second-zreco_yreco[iz-1].second, 2)+
		//pow(zreco_widreco[iz].second-zreco_widreco[iz-1].second, 2) );
		//reco_trklen_accum[iz] = range_reco;
		//zreco_lenreco.push_back(make_pair(zreco_widreco[iz].first, range_reco));
		//}
		//cout<<"range_reco:"<<range_reco<<endl;

		//Reco stopping/Inel p cut ---------------------------------------------------------------------------------------------------------//
		bool IsRecoStop=false;
		bool IsRecoInEL=false;
		bool IsRecoEL=false;
		double mom_beam_spec=-99; mom_beam_spec=beamMomentum_spec->at(0);
		//double range_reco=-99; if (!primtrk_range->empty()) range_reco=primtrk_range->at(0); //reco primary trklen
		double csda_val_spec=csda_range_vs_mom_sm->Eval(mom_beam_spec);

		//if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) IsRecoStop=true; //old cut
		//if ((range_reco/csda_val_spec)<min_norm_trklen_csda) IsRecoInEL=true; //old cut
		
		if ((range_reco/csda_val_spec)<min_norm_trklen_csda) { //inel region
			if (pid>pid_1) IsRecoInEL=true; 
			if (pid<=pid_1) IsRecoStop=true; 
		} //inel region
		if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) { //stopping p region
			if (pid>pid_2) IsRecoInEL=true; 
			if (pid<=pid_2) IsRecoStop=true;
		} //stopping p region


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


		//Slice ID definition ---------------------------------------------------------------------------------------------------------//
		//true slice ID
		true_sliceID = int(range_true/thinslicewidth); //HY:Make sure size of true_trk_len vector is !0, otherwise lost truth info
		if (true_sliceID < 0) true_sliceID = -1;
		if (true_endz < 0) true_sliceID = -1; //hy added
		if (true_sliceID >= nthinslices) true_sliceID = nthinslices;

		//reco slice ID
		if (IsPandoraSlice&&IsCaloSize) {
			reco_sliceID = int(range_reco/thinslicewidth);
			if (reco_sliceID < 0) reco_sliceID = -1;
			if (reco_endz<0) reco_sliceID = -1;
			if (reco_sliceID >= nthinslices) reco_sliceID = nthinslices;
		}
		//-----------------------------------------------------------------------------------------------------------------------------//

		//some labels with numbers
		//define cuts 
		//define true labels
		//(MC/data) vs reco Slice ID using bkg-rich sample
		//[1]misID:p rich sample
		//if (IsPos&&IsCaloSize&&IsPandoraSlice) { //if Pos
		//if (IsBQ&&IsCaloSize&&IsPandoraSlice) { //if Pos
		//if (!IsBQ&&IsCaloSize&&IsPandoraSlice&&range_reco<=10) { //if not BQ
		//if (IsBQ&&IsCaloSize&&IsPandoraSlice&&range_reco<=20&&range_reco>10) { //if BQ
		//if (beamtrk_z->size()) { //if BQ
		//if (!beamtrk_Eng->empty()) { //if true container not empty
		//if (!beamtrk_Eng->empty()&&IsCaloSize&&IsPandoraSlice) { //if true container not empty
		if (!beamtrk_Eng->empty()) { //if true container not empty
			//if (kinel&&reco_sliceID<=1) {
			//if (reco_sliceID<=5) {
			//if (kMIDp) {
			//if (IsPureEL&&range_reco>=100) {
			//if (IsPureInEL&&range_true>=26&&range_true<28) {
			//if (IsPureInEL&&range_true>12&&range_true<13) {
			//if (IsPureInEL&&range_true>4&&range_true<5.5) {
			//if (IsPureInEL&&range_true>=11&&range_true<15) {
			//if (IsMisIDP&&range_reco<=5.&&IsPos&&IsCaloSize&&IsPandoraSlice&&IsRecoInEL) {
			//hunting for a specific event ------------------------------------------------------------//
			bool IsHunt=false;
				for (int hk=0; hk<(int)beamtrk_Eng->size()-1; ++hk) { //loop over true hits
					double thisZ=beamtrk_z->at(hk);
					//int this_sliceID = int(thisZ/thinslicewidth);
					double thisLen=true_trklen_accum[hk];
					int this_sliceID = int(thisLen/thinslicewidth);
					double this_incE = 1000.*beamtrk_Eng->at(hk); //MeV

					//if (this_sliceID>=nthinslices) continue;
					//if (this_sliceID<0) continue;
					//vincE_true[this_sliceID].push_back(this_incE);

					if (thisLen<2&&thisZ<3&&key_reach_tpc>=0) IsHunt=true;
				}//loop over true hits (last point always has KE = 0)

			//hunting for a specific event ------------------------------------------------------------//
			//
			//if (IsPureInEL&&true_endz<4.&&range_true<2&&key_reach_tpc>=0) { //is pure inelastic
			if (IsHunt) { //is pure inelastic

				TString tit;
				tit=Form("run%d subrun:%d evt:%d trackID:%d",run,subrun,event,primaryID); 

				TString fig_out;
				fig_out=Form("/dune/data2/users/hyliao/evt_display/short_trk_inel/run%d_event%d_trackID%d_pdg%d_trklen%.2fcm_costheta%.2f.png",run,event,primaryID,beamtrackPdg,range_reco,cosine_beam_spec_primtrk);

				gStyle->SetOptStat(0); //no statistics box
				gStyle->SetFrameBorderSize(3);

				TCanvas *c1 = new TCanvas(Form("c1"),tit.Data(),1600,800);
				//cnt_canvas++;
				//c1->SetName(Form("c_%d",cnt_canvas));
				c1->SetTitle(tit.Data());
				c1->Divide(2,1);

				//draw track in y-z view
				//after entering TPC
				TGraph *gr_trk_yz;
				TGraph *gr_trk_xz;
				TGraph *gr_trkend_yz;
				TGraph *gr_trkend_xz;

				if (primtrk_hitz->size()) { //if reco container not empty
					gr_trk_yz=new TGraph(primtrk_hitz->size(), &primtrk_hitz->at(0), &primtrk_hity->at(0));
					gr_trk_xz=new TGraph(primtrk_hitz->size(), &primtrk_hitz->at(0), &primtrk_hitx->at(0));
					gr_trk_yz->SetMarkerStyle(20);
					gr_trk_xz->SetMarkerStyle(20);
					gr_trk_yz->SetMarkerColor(2);
					gr_trk_xz->SetMarkerColor(2);

					//after entering TPC, trkend
					gr_trkend_yz=new TGraph(1, &primtrk_hitz->at(-1+primtrk_hitz->size()), &primtrk_hity->at(-1+primtrk_hitz->size()));
					gr_trkend_xz=new TGraph(1, &primtrk_hitz->at(-1+primtrk_hitz->size()), &primtrk_hitx->at(-1+primtrk_hitz->size()));
					gr_trkend_yz->SetMarkerStyle(20);
					gr_trkend_xz->SetMarkerStyle(20);
					gr_trkend_yz->SetMarkerColor(3);
					gr_trkend_xz->SetMarkerColor(3);
				} //if reco container not empty
				else {
					float fake[1]={-999};
					gr_trk_yz=new TGraph(1, &fake[0], &fake[0]);
					gr_trk_xz=new TGraph(1, &fake[0], &fake[0]);
				}


				//before tpc
				//TGraph *gr_btrk_yz=new TGraph(key_reach_tpc+1, &beamtrk_z->at(0), &beamtrk_y->at(0));
				//TGraph *gr_btrk_xz=new TGraph(key_reach_tpc+1, &beamtrk_z->at(0), &beamtrk_x->at(0));
				TGraph *gr_btrk_yz=new TGraph(beamtrk_z->size(), &beamtrk_z->at(0), &beamtrk_y->at(0));
				TGraph *gr_btrk_xz=new TGraph(beamtrk_z->size(), &beamtrk_z->at(0), &beamtrk_x->at(0));
				gr_btrk_yz->SetMarkerStyle(24);
				gr_btrk_xz->SetMarkerStyle(24);
				gr_btrk_yz->SetMarkerColor(4);
				gr_btrk_xz->SetMarkerColor(4);

				//TH2D *h2d_yz=new TH2D("h2d_yz","",170,-50,120,120,350,450); //bkg
				TH2D *h2d_yz=new TH2D("h2d_yz","",25, -5, 20, 120, 350, 450); //bkg
				c1->cd(1);
				h2d_yz->GetXaxis()->SetTitle("Z [cm]");
				h2d_yz->GetYaxis()->SetTitle("Y [cm]");
				h2d_yz->SetTitle(tit.Data());
				h2d_yz->Draw();
				if (primtrk_hitz->size()) { //if reco container not empty
					gr_trk_yz->Draw("p same");
				 	//gr_trkend_yz->Draw("p same");
				} //if reco container not empty
				gr_btrk_yz->Draw("p same");

				//true vertex
				//if (interactionProcesslist->size()){
				TGraph **mgr_vtx_yz=new TGraph*[interactionProcesslist->size()]; 
				for (size_t k=0; k<interactionProcesslist->size(); ++k) {
				  try {
					mgr_vtx_yz[k]=new TGraph(1, &interactionZ->at(k), &interactionY->at(k));
					mgr_vtx_yz[k]->SetMarkerSize(2);
					if(interactionProcesslist->at(k)=="hadElastic") {
						mgr_vtx_yz[k]->SetMarkerStyle(4);
						mgr_vtx_yz[k]->SetMarkerColor(6);
					}
					if(interactionProcesslist->at(k)=="protonInelastic") {
						mgr_vtx_yz[k]->SetMarkerStyle(5);
						mgr_vtx_yz[k]->SetMarkerColor(6);
					}
					c1->cd(1); 
					mgr_vtx_yz[k]->Draw("p same");
				  }
				  catch (const std::out_of_range & ex) {
			   	    std::cout << "out_of_range Exception Caught :: interactionProcesslist" << ex.what() << std::endl;
			   	  }
				}

				cout<<"ck2"<<""<<endl;

				//TH2D *h2d_xz=new TH2D("h2d_xz","",170,-50,120,120,-120,-10);
				TH2D *h2d_xz=new TH2D("h2d_xz","",25, -5, 20,120,-120,-10);
				c1->cd(2);
				h2d_xz->GetXaxis()->SetTitle("Z [cm]");
				h2d_xz->GetYaxis()->SetTitle("X [cm]");
				TString str_right;
				//if (IsFlip==false) { 
					//str_right=Form("reco_len:%.1f cm|true_len:%.1f cm|SID:%d|cos#Theta:%.2f|#chi^{2}:%.2f|ntrklen:%.2f", range_reco, range_true, reco_sliceID, cosine_beam_spec_primtrk,pid,range_reco/csda_val_spec);
					//str_right=Form("true_len:%.1f cm|true_endz:%.1f|reco_len:%.1f cm|cos#Theta:%.2f|#chi^{2}:%.2f|ntrklen:%.2f", range_true, true_endz, range_reco, cosine_beam_spec_primtrk,pid,range_reco/csda_val_spec);
					//str_right=Form("true_len:%.1f cm|true_endz:%.1f|reco_len:%.1f cm| KE_ff:%.1f MeV| cos#Theta:%.2f", range_true, true_endz, range_reco, KE_ff, cosine_beam_spec_primtrk);
					if (key_reach_tpc>=0) str_right=Form("true_len:%.1f cm|true_endz:%.1f| KE_ff:%.1f MeV| Z_st:%.1f cm| range_patch:%.2f cm", range_true, true_endz, KE_ff, beamtrk_z->at(key_reach_tpc),range_true_patch);
					else str_right=Form("true_len:%.1f cm|true_endz:%.1f| KE_ff:%.1f MeV| beam NOT reach TPC", range_true, true_endz, KE_ff);

				//}
				//if (IsFlip==true) { 
					//str_right=Form("reco_len:%.1f cm | true_trklen:%.1f cm | SliceID:%d | cos#theta:-%.2f", range_reco, range_true, reco_sliceID, cosine_beam_spec_primtrk);
				//}
				h2d_xz->SetTitle(str_right.Data());
				h2d_xz->Draw();
				gr_trk_xz->Draw("p same");
				gr_btrk_xz->Draw("p same");
				//gr_trkend_xz->Draw("p same");

				cout<<"ck3"<<"\n"<<endl;

				TGraph **mgr_vtx_xz=new TGraph*[interactionProcesslist->size()]; 
				for (size_t k=0; k<interactionProcesslist->size(); ++k) {
				  try {
					mgr_vtx_xz[k]=new TGraph(1, &interactionZ->at(k), &interactionX->at(k));
					mgr_vtx_xz[k]->SetMarkerSize(2);
					if(interactionProcesslist->at(k)=="hadElastic") {
						mgr_vtx_xz[k]->SetMarkerStyle(4);
						mgr_vtx_xz[k]->SetMarkerColor(6);
					}
					if(interactionProcesslist->at(k)=="protonInelastic") {
						mgr_vtx_xz[k]->SetMarkerStyle(5);
						mgr_vtx_xz[k]->SetMarkerColor(6);
					} 
					c1->cd(2);
					mgr_vtx_xz[k]->Draw("p same");
				  }
				  catch (const std::out_of_range & ex) {
			   	    std::cout << "out_of_range Exception Caught :: interactionProcesslist" << ex.what() << std::endl;
			   	  }
				}

				c1->Modified(); c1->Update();
				c1->Print(fig_out.Data());


				//for (size_t k=0; k<interactionProcesslist->size(); ++k) {
					//delete mgr_vtx_yz[k];
					//delete mgr_vtx_xz[k];
				//}

				delete c1;
				delete gr_trk_yz;
				delete gr_trk_xz;
				//delete gr_trkend_yz;
				//delete gr_trkend_xz;
				delete gr_btrk_yz;
				delete gr_btrk_xz;

				delete h2d_yz;
				delete h2d_xz;

			}



		} //if Pos

		//some ke calc. -------------------------------------------------------------------------------------//






		} //main entry loop

		//save results -------//
		//SaveHistograms();



}
