#define ProtonKEReweight_cxx
#include "ProtonKEReweight.h"

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


void ProtonKEReweight::Loop() {
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
	//BetheBloch BB;
	//BB.SetPdgCode(pdg);

	//Kinetic energies -------------------------------------------------------------------------------------------------------//
	//int nke=320;
	int nke=160;
	double kemin=-800;
	double kemax=800;
	//------------------------------------------------------------------------------------------------------------------------//

	//Beam momentum reweighting ----------------------------------------------------------------------------------------------//
	//MC Beam Mom Gaussian 
	//double m1=1007.1482; //MC prod4a [spec]
	//double s1=60.703307; //MC prod4a [spec]

	//MC Beam KE Gaussian 
	double m1=3.89270e+02; //MC prod4a [spec]
	double s1=4.49638e+01; //MC prod4a [spec]


	//momentum cut range	
	double mu_min=m1-3.*s1;
	double mu_max=m1+3.*s1;

	//default gaussian
	double xmin=0.; //pmin [MeV/c]
	//double xmax=2000.; //pmax [MeV/c]
	double xmax=1000.; //pmax [MeV/c]
	TF1 *g1=new TF1("g1",fitg,xmin,xmax,2);
	g1->SetName("g1");
	g1->SetParameter(0,m1);
	g1->SetParameter(1,s1);

	//mu range
	double dmu=0.0005;
	double mu_st=1.01;
	//int nmu=71;
	int nmu=141;

	double dsigma=0.002;
	//double sigma_st=1.5;
	//int nsigma=250;
	double sigma_st=1.2;
	int nsigma=350;

	//mu x sigma
	const int n_mu_sigma=(const int)nmu*nsigma;
	int n_1d=nmu*nsigma; 
	TF1 **gn=new TF1*[n_mu_sigma];
	TF1 **gng=new TF1*[n_mu_sigma];
	//use KEcalo-constE as an observable for reweighting
	TH1D *h1d_KEcalo_RecoEl_rw[n_mu_sigma];

	int cnt_array=0;
	int index_original=0;
	//int index_minchi2=13331; //index of minchi2
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
			h1d_KEcalo_RecoEl_rw[cnt_array]=new TH1D(Form("h1d_KEcalo_RecoEl_rw_%d",cnt_array),Form("f_{#mu}:%.2f f_{#sigma}:%.2f #oplus RecoStop Cut",frac_mu,frac_sigma),nke, kemin, kemax);

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

				//double r0=pow(dedx_reco,0.42)/17.-rr_reco;
				//double dedx_pred=17.*pow(rr_reco,-0.42);

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
		//vector<double> true_trklen_accum;
		//true_trklen_accum.reserve(beamtrk_z->size()); // initialize true_trklen_accum
		for (int iz=0; iz<(int)beamtrk_z->size(); iz++) {
			if (abs(beamtrk_z->at(iz)) < tmp_z){
				tmp_z = abs(beamtrk_z->at(iz));
				key_st = iz; // find the point where the beam enters the TPC (find the smallest abs(Z))
			}
			//if (is_beam_at_ff) true_trklen_accum[iz] = 0.; // initialize true_trklen_accum [beam at ff]
			//if (!is_beam_at_ff) true_trklen_accum[iz] = -1; // initialize true_trklen_accum [beam not at ff]
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
		//double av_r0=0;
		//int n_r0=0;
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

				//double tmp_r0=pow((17./cali_dedx),1./0.42)-resrange_reco;
				//R0.push_back(tmp_r0);
				//av_r0+=tmp_r0;
				//n_r0++;

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
		//av_r0/=(double)n_r0;

		//KEff for el
        	//double pff_truth=ke2p(ke_ff/1000.);
        	//double csda_truth=csda_range_vs_mom_sm->Eval(pff_truth);
		//double KEcsda=1000.*ke_vs_csda_range_sm->Eval(range_reco);

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
		//double KEend_true=1000.*(beamtrk_Eng->at(-2+beamtrk_Eng->size()));
		//double KEbb_true=-1; KEbb_true=BB.KEAtLength(ke_ff, range_true);

		//double KEbb_reco_constErange=-1; KEbb_reco_constErange=BB.KEAtLength(ke_beam_spec_MeV-mean_Elossrange_stop, range_reco);
		//double KEbb_reco_Edept_range=-1; KEbb_reco_Edept_range=BB.KEAtLength(KEcsda, range_reco);
		//double KEbb_reco_Edept_range_inel=-1; KEbb_reco_Edept_range_inel=BB.KEAtLength(KEcsda_inel, range_reco);

		double KEcalo_reco_constEcalo=-1; KEcalo_reco_constEcalo=ke_beam_spec_MeV-mean_Elosscalo_stop-ke_calo_MeV;
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


			//double mom_rw_minchi2=1;
			//if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) mom_rw_minchi2=gng[index_minchi2]->Eval(mom_beam_spec*1000.); //bmrw

			//h2d_trklen_KEcalo->Fill(range_reco, KEcalo_reco_constEcalo, mom_rw_minchi2);
			//h1d_kecalo_reco_constEcalo->Fill(KEcalo_reco_constEcalo, mom_rw_minchi2);
			//h1d_Edept_reco->Fill(ke_calo_MeV, mom_rw_minchi2);



			if (IsRecoEL) { //IsRecoEL
				//h1d_kebb_reco->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
				//h1d_kebb_recoEl->Fill(KEbb_reco_Edept_range, mom_rw_minchi2);
				//h2d_trklen_KEcalo_RecoEl->Fill(range_reco, KEcalo_reco_constEcalo, mom_rw_minchi2);


					double mom_rw=1.;
					if (KEcalo_reco_constEcalo>=mu_min&&KEcalo_reco_constEcalo<=mu_max) { //ke cut (within 3-sigma)
						for (int ig = 0; ig < n_1d; ++ig) { //rw loop
							mom_rw=gng[ig]->Eval(KEcalo_reco_constEcalo);
							h1d_KEcalo_RecoEl_rw[ig]->Fill(KEcalo_reco_constEcalo,mom_rw); //beam-mom rw using stopping protons
						} //rw loop
					} //beam-mom cut (within 3-sigma)
					else { //tail of beam
					//if ((mom_beam_spec*1000.)<mu_min||(mom_beam_spec*1000.)>mu_max) { //tail of the beam
						for (int ig = 0; ig < n_1d; ++ig) { //rw loop
							h1d_KEcalo_RecoEl_rw[ig]->Fill(KEcalo_reco_constEcalo); //beam-mom rw 
						} //rw loop
					} //tail of the beam	


			} //IsRecoEL





		} //basic cuts

	} //main entry loop


	//save results...
   	//TFile *fout = new TFile("mc_ke_advancedKEff_nobmrw.root","RECREATE");
   	//TFile *fout = new TFile("mc_kecalo_bmrw.root","RECREATE");
   	//TFile *fout = new TFile("mc_kecalo_nobmrw_ElossTune.root","RECREATE");
   	TFile *fout = new TFile("mc_kerw_tunning.root","RECREATE");
		for (int ig = 0; ig < n_1d; ++ig) { //rw loop
			h1d_KEcalo_RecoEl_rw[ig]->Write();
		} //rw loop
	fout->Close();



}
