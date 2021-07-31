#define ProtonMomentumReweight_cxx
#include "ProtonMomentumReweight.h"

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
//#include "./headers/SliceParams.h"
//#include "./headers/ThinSlice.h"

using namespace std;
using namespace ROOT::Math;

void ProtonMomentumReweight::Loop() {
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


	int nx=1000;	
	double xmin=0.; //pmin [MeV/c]
	double xmax=2000.; //pmax [MeV/c]
	TH1D *h1d_pbeam=new TH1D("h1d_pbeam","",nx,xmin,xmax);
	TH1D *h1d_pbeam_stop=new TH1D("h1d_pbeam_stop","",nx,xmin,xmax);
	TH1D *h1d_pff_stop=new TH1D("h1d_pff_stop","",nx,xmin,xmax);

	//MC Beam Mom Gaussian 
	double m1=1007.1482; //MC prod4a [spec]
	double s1=60.703307; //MC prod4a [spec]

	//momentum cut range	
	double mu_min=m1-3.*s1;
	double mu_max=m1+3.*s1;

	//default gaussian
	TF1 *g1=new TF1("g1",fitg,xmin,xmax,2);
	g1->SetName("g1");
	g1->SetParameter(0,m1);
	g1->SetParameter(1,s1);


	//mu range
	double dmu=0.0005;
	double mu_st=1.01;
	int nmu=71;

	double dsigma=0.002;
	double sigma_st=1.5;
	int nsigma=250;

	//mu x sigma
	const int n_mu_sigma=(const int)nmu*nsigma;
	int n_1d=nmu*nsigma; 
	TF1 **gn=new TF1*[n_mu_sigma];
	TF1 **gng=new TF1*[n_mu_sigma];

	//use trklen as an observable for reweighting
	TH1D *h1d_trklen_rw[n_mu_sigma];

	int cnt_array=0;
	int index_original=0;
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
			h1d_trklen_rw[cnt_array]=new TH1D(Form("h1d_trklen_rw_%d",cnt_array),Form("f_{#mu}:%.2f f_{#sigma}:%.2f #oplus RecoStop Cut",frac_mu,frac_sigma),n_b,b_min,b_max);
			h1d_trklen_rw[cnt_array]->GetXaxis()->SetTitle("Track Length [cm]");

			cnt_array++;
			} //sigma loop
	} //mu loop


	//trklen
	TH1D *h1d_trklen_stop=new TH1D(Form("h1d_trklen_stop"),Form("reco stop"),n_b,b_min,b_max);
	TH1D *h1d_trklen_stop_XY=new TH1D(Form("h1d_trklen_stop_XY"),Form("reco stop with xy cut"),n_b,b_min,b_max);
	
	//Weighted Gaussians & Weighting function ----------------------------------------------------------------------------------------------------------------------------------------------//


	//KEs ---------------------------------------------------------------------------------------------------------------------------------//
	int nx_trklen=150;
	double xmin_trklen=0;
	double xmax_trklen=150;
	int ny_edept=400;
	double ymin_edept=0;
	double ymax_edept=800;

	TH1D *h1d_kebeam=new TH1D("h1d_kebeam","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_keff=new TH1D("h1d_keff","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_keff_stop=new TH1D("h1d_keff_stop","",ny_edept,ymin_edept,ymax_edept);
	TH1D *h1d_kerange_stop=new TH1D("h1d_kerange_stop","", ny_edept, ymin_edept, ymax_edept);
	TH1D *h1d_kecalo_stop=new TH1D("h1d_kecalo_stop","", ny_edept, ymin_edept, ymax_edept);

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
		bool IsPureMCS=false; //no hadron scattering

		//if (primary_truth_EndProcess->c_str()!=NULL) n_true_end++;
		if (strcmp(primary_truth_EndProcess->c_str(),"protonInelastic")==0) {
			IsPureInEL=true;
		}
		else { //if primarytrue_end!=InEL
			if (interactionProcesslist->size()) { //size of interactionProcesslist >=0
				for(size_t iiii=0; iiii<interactionProcesslist->size(); iiii++) { //loop over all true interaction hits in this track
					try {
						double intx=interactionX->at(iiii);
						double inty=interactionY->at(iiii);
						double intz=interactionZ->at(iiii); 

						if(strcmp(interactionProcesslist->at(iiii).c_str(),"hadElastic")==0) {
							IsPureEL=1;
						}
					}
					catch (const std::out_of_range & ex) {
						//std::cout << "out_of_range Exception Caught :: interactionProcesslist" << ex.what() << std::endl;
					}
					//if (intz<0) { //if interaction outside tpc
					//if(strcmp(interactionProcesslist->at(iiii).c_str(),"Transportation")!=0) {
				} //loop over all true interaction hits in this track 
			} //size of interactionProcesslist >=0
		} //if primarytrue_end!=InEL

		if (IsPureInEL==0&&IsPureEL==0) IsPureMCS=1;
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
		bool kmcs=false;
		if (IsBeamMatch) { //beam-match
			if (IsPureInEL) kinel=true;
			if (IsPureEL) kel=true;
			if (IsPureMCS) kmcs=true;
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

			//dx_reco_stx__bposx_ff=reco_stx-bposx_ff;
			//dy_reco_sty__bposy_ff=reco_sty-bposy_ff;
			//dz_reco_stz__bposz_ff=reco_stz-bposz_ff;

			if (reco_stz>min1_z&&reco_stz<min2_z) {
				if (reco_sty>min1_y&&reco_sty<min2_y) {
					if (reco_stx>min1_x&&reco_stx<min2_x) {
						IsPos=true;
					}
				}
			}
		}

		//cosine_theta/cut
		bool IsCosine=false;
		double cosine_beam_spec_primtrk=-999; 
		cosine_beam_spec_primtrk=beamDirx_spec->at(0)*primaryStartDirection[0]+beamDiry_spec->at(0)*primaryStartDirection[1]+beamDirz_spec->at(0)*primaryStartDirection[2]; //cosine between beam_spec and primary trk direction
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

		//Reco stopping/Inel p cut
		bool IsRecoStop=false;
		bool IsRecoInEL=false;
		double mom_beam_spec=-99; mom_beam_spec=beamMomentum_spec->at(0);
		double range_reco=-99; if (!primtrk_range->empty()) range_reco=primtrk_range->at(0); //reco primary trklen
		double csda_val_spec=csda_range_vs_mom_sm->Eval(mom_beam_spec);

		if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) IsRecoStop=true;
		if ((range_reco/csda_val_spec)<min_norm_trklen_csda) IsRecoInEL=true;

		//kinetic energies
		//double ke_beam=1000.*p2ke(mom_beam); //ke_beam
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
		if (IsPandoraSlice&&IsBQ&&IsCaloSize) { //if calo size not empty
			//no-sce
			h2d_xy_noSCE->Fill(reco_stx_noSCE, reco_sty_noSCE);
			h1d_zst_noSCE->Fill(reco_stz_noSCE);
			//after sce
			h2d_xy_SCE->Fill(reco_stx, reco_sty);
			h1d_zst_SCE->Fill(reco_stz);
			
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

				if (IsRecoStop) {
					h2d_rr_dedx_recoSTOP->Fill(resrange_reco, cali_dedx);
				}

			} //loop over reco hits of a given track

			if (IsRecoStop) chi2pid_recostop->Fill(chi2pid(trkdedx,trkres));
			if (IsRecoInEL) chi2pid_recoinel->Fill(chi2pid(trkdedx,trkres));
			if (IsPureMCS) chi2pid_truestop->Fill(chi2pid(trkdedx,trkres));
			if (IsPureEL) chi2pid_trueel->Fill(chi2pid(trkdedx,trkres));
			if (IsPureInEL) chi2pid_trueinel->Fill(chi2pid(trkdedx,trkres)); 

		} //if calo size not empty

		if (IsPandoraSlice&&IsBQ&&IsCaloSize) { //basic cuts
			h1d_kebeam->Fill(ke_beam_spec_MeV);
			h1d_pbeam->Fill(1000.*mom_beam_spec);
			h1d_keff->Fill(ke_ff);
			if (IsRecoStop) { //reco stop 
				h1d_keff_stop->Fill(ke_ff);
				h1d_kerange_stop->Fill(ke_trklen_MeV);
				h1d_kecalo_stop->Fill(ke_calo_MeV);
				h1d_trklen_stop->Fill(range_reco);
				h1d_pbeam_stop->Fill(1000.*ke2p(ke_trklen));
				h1d_pff_stop->Fill(1000.*ke2p(ke_ff/1000.));

				if (IsXY) { //xy-cut
					h1d_trklen_stop_XY->Fill(range_reco);
					double mom_rw=1.;
					if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) { //beam-mom cut (within 3-sigma)
						for (int ig = 0; ig < n_1d; ++ig) { //rw loop
							mom_rw=gng[ig]->Eval(mom_beam_spec*1000.);
							h1d_trklen_rw[ig]->Fill(range_reco,mom_rw); //beam-mom rw using stopping protons
						} //rw loop
					} //beam-mom cut (within 3-sigma)
					else { //tail of beam
					//if ((mom_beam_spec*1000.)<mu_min||(mom_beam_spec*1000.)>mu_max) { //tail of the beam
						for (int ig = 0; ig < n_1d; ++ig) { //rw loop
							h1d_trklen_rw[ig]->Fill(range_reco); //beam-mom rw 
						} //rw loop
					} //tail of the beam	
				} //xy-cut
			} //reco stop
		} //basic cuts

	} //main entry loop

	//Fit Beam ...
	TF1* kebeam_fit; kebeam_fit=VFit(h1d_kebeam, 2);
	kebeam_fit->SetName("kebeam_fit");
	TF1* pbeam_fit; pbeam_fit=VFit(h1d_pbeam, 2);
	pbeam_fit->SetName("pbeam_fit");

	TF1* pbeam_stop_fit; pbeam_stop_fit=VFit(h1d_pbeam_stop, 2);
	pbeam_stop_fit->SetName("pbeam_stop_fit");

	TF1* pff_stop_fit; pff_stop_fit=VFit(h1d_pff_stop, 2);
	pff_stop_fit->SetName("pff_stop_fit");


	//save results...
   	TFile *fout = new TFile("mc_proton_bmrw.root","RECREATE");
				
		h2d_rr_dedx_recoSTOP->Write();
		gr_predict_dedx_resrange->Write();

		h1d_pbeam->Write();
		h1d_kebeam->Write();
		h1d_pbeam_stop->Write();
		h1d_pff_stop->Write();

		kebeam_fit->Write();
		pbeam_fit->Write();
		pbeam_stop_fit->Write();
		pff_stop_fit->Write();

		h1d_keff->Write();
		h1d_keff_stop->Write();
		
		h1d_kerange_stop->Write();
		h1d_kecalo_stop->Write();



		h2d_xy_noSCE->Write();
		h2d_xy_SCE->Write();
		
		h1d_zst_noSCE->Write();
		h1d_zst_SCE->Write();

		chi2pid_recostop->Write();
		chi2pid_recoinel->Write();
		chi2pid_truestop->Write();
		chi2pid_trueel->Write();
		chi2pid_trueinel->Write();

		h1d_trklen_stop->Write();
		h1d_trklen_stop_XY->Write();

		for (int ig = 0; ig < n_1d; ++ig) { //rw loop
			h1d_trklen_rw[ig]->Write();
		} //rw loop

	fout->Close();



}
