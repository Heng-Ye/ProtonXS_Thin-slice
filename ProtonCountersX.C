#define ProtonCounters_cxx
#include "ProtonCounters.h"

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
#include "./headers/SliceParams.h"
#include "./headers/ThinSlice.h"

using namespace std;
using namespace ROOT::Math;

void ProtonCounters::Loop() {
	if (fChain == 0) return;

	//various counters --------------------------------------------------------------------------------------------------------------------------------------------------------------------//
	int n_processmap_error=0; //sansity check: processmap
	int n_true_end=0; //sansity check: true-end label

	int n_tot=0, n_sg=0, n_bkg=0; //no cut
	int n_el=0, n_inel=0, n_mcs=0, n_midcosmic=0, n_midpi=0, n_midp=0, n_midmu=0, n_mideg=0, n_midother=0;

	int n_pan_tot=0; //pandora cut
	int n_el_pan=0, n_inel_pan=0, n_mcs_pan=0, n_midcosmic_pan=0, n_midpi_pan=0, n_midp_pan=0, n_midmu_pan=0, n_mideg_pan=0, n_midother_pan=0;

	int n_calsz_tot=0; //calosz cut
	int n_el_calsz=0, n_inel_calsz=0, n_mcs_calsz=0, n_midcosmic_calsz=0, n_midpi_calsz=0, n_midp_calsz=0, n_midmu_calsz=0, n_mideg_calsz=0, n_midother_calsz=0;

	int n_bq_tot=0; //bq cut
	int n_el_bq=0, n_inel_bq=0, n_mcs_bq=0, n_midcosmic_bq=0, n_midpi_bq=0, n_midp_bq=0, n_midmu_bq=0, n_mideg_bq=0, n_midother_bq=0;

	int n_recoinel_tot=0; //recoinel cut
	int n_el_recoinel=0, n_inel_recoinel=0, n_mcs_recoinel=0, n_midcosmic_recoinel=0, n_midpi_recoinel=0, n_midp_recoinel=0, n_midmu_recoinel=0, n_mideg_recoinel=0, n_midother_recoinel=0;


	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	Long64_t nbytes = 0, nb = 0;
	//bool isTestSample=true;
	//int true_sliceID = -1, reco_sliceID = -1;
	for (Long64_t jentry=0; jentry<nentries;jentry++) { //main entry loop
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		//only select protons	
		//if (primary_truth_Pdg!=pdg) continue; //only interested in protons
		if (beamtrackPdg!=pdg) continue; //only interested in protons
		//std::cout<<"beamtrackPdg:"<<beamtrackPdg<<std::endl;
		n_tot++;

		//MC beam momentum -----------------------//
		double mm1=1007.1482; //MC prod4a [spec]
		double ss1=60.703307; //MC prod4a [spec]
		double mu_min=mm1-3.*ss1;
		double mu_max=mm1+3.*ss1;

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
				} //loop over all true interaction hits in this track 
			} //size of interactionProcesslist >=0
			*/
		} //hIoni

		//if (IsPureInEL==0&&IsPureEL==0) {
			//IsPureMCS=1;
		//}
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

				if (h==0) range_reco=0;
				if (h>=1) {
					range_reco += sqrt( pow(primtrk_hitx->at(h)-primtrk_hitx->at(h-1), 2)+
							pow(primtrk_hity->at(h)-primtrk_hity->at(h-1), 2)+
							pow(primtrk_hitz->at(h)-primtrk_hitz->at(h-1), 2) );
					reco_trklen_accum.push_back(range_reco);
				}

				reco_calo_MeV+=cali_dedx*pitch;
				kereco_range+=pitch*dedx_predict(resrange_reco);
				kereco_range2+=pitch*(double)gr_predict_dedx_resrange->Eval(resrange_reco);

				trkdedx.push_back(cali_dedx);
				trkres.push_back(resrange_reco);

			} //loop over reco hits of a given track

			pid=chi2pid(trkdedx,trkres); //pid using stopping proton hypothesis

		} //if calo size not empty
	
		//xy-cut
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

		//beam XY cut to remove E-loss events upstream
		bool IsBeamXY=false;
		double bx_spec=beamPosx_spec->at(0);
		double by_spec=beamPosy_spec->at(0);
		double mom_beam_spec=-99; mom_beam_spec=beamMomentum_spec->at(0);

		if ((pow(((bx_spec-meanX_mc)/(1.5*rmsX_mc)),2)+pow(((by_spec-meanY_mc)/(1.5*rmsY_mc)),2))<=1.) IsBeamXY=true;

		//beam-mom cut (within 3-sigma)
		bool IsBeamMom=false;
		if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) IsBeamMom=true;

		//Intersection cut
		//bool IsIntersection=false;		
		//if (timeintersection->size()) IsIntersection=true; //over-lapping track cut

		//beam quality cut
		bool IsBQ=false;
		//if (IsCosine&&IsPos) IsBQ=true;
		if (IsBeamXY&&IsBeamMom&&IsCosine&&IsPos) IsBQ=true;
	
		//Reco stopping/Inel p cut
		bool IsRecoStop=false;
		bool IsRecoInEL=false;
		bool IsRecoEL=false;
		//double range_reco=-99; if (!primtrk_range->empty()) range_reco=primtrk_range->at(0); //reco primary trklen
		double csda_val_spec=csda_range_vs_mom_sm->Eval(mom_beam_spec);

		if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) IsRecoStop=true;
		//if ((range_reco/csda_val_spec)<min_norm_trklen_csda) IsRecoInEL=true;

		if ((range_reco/csda_val_spec)<min_norm_trklen_csda) { //inel region
			if (pid>pid_1) IsRecoInEL=true; 
			if (pid<=pid_1) IsRecoEL=true; 
		} //inel region
		if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) { //stopping p region
			if (pid>pid_2) IsRecoInEL=true; 
			if (pid<=pid_2) IsRecoEL=true;
		} //stopping p region


		//kinetic energies
		//double ke_beam=1000.*p2ke(mom_beam); //ke_beam
		double ke_beam_spec=p2ke(mom_beam_spec); //ke_beam_spec [GeV]
		double ke_beam_spec_MeV=1000.*ke_beam_spec; //ke_beam_spec [MeV]
		double ke_trklen=1000.*ke_vs_csda_range_sm->Eval(range_reco); //[unit: MeV]
		double p_trklen=ke2p(ke_trklen);

		//countings -------------------------------------------------------//
		if (IsPandoraSlice) n_pan_tot++;
		if (IsPandoraSlice&&IsCaloSize) n_calsz_tot++;
		if (IsPandoraSlice&&IsCaloSize&&IsBQ) n_bq_tot++;
		if (IsPandoraSlice&&IsCaloSize&&IsBQ&&IsRecoInEL) n_recoinel_tot++;

		//if (IsBeamMatch) { //beam-match
			//if (IsPureInEL) kinel=true;
			//if (IsPureEL) kel=true;
			//if (IsPureMCS) kmcs=true;
		//} //beam-match

		//[0]pure inel
		//if (IsBeamMatch&&IsPureInEL) { //pure inel
		//if (kinel) { //pure inel+beam match(construct the right particle)
		if (IsPureInEL) { //pure inel
			n_inel++;
			if (IsBeamMatch&&IsPandoraSlice) { //pandora
				n_inel_pan++;
				if (IsCaloSize) { //calosz
					n_inel_calsz++;
					if (IsBQ) { //bq
						n_inel_bq++;
						if (IsRecoInEL) { //reco inel
							n_inel_recoinel++;
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //pure inel		

		//[1]pure el
		//if (IsBeamMatch&&IsPureEL) { //pure el
		//if (kel) { //pure el+beam match
		if (IsPureEL) { //pure el
			n_el++;
			if (IsBeamMatch&&IsPandoraSlice) { //pandora
				n_el_pan++;
				if (IsCaloSize) { //calosz
					n_el_calsz++;
					if (IsBQ) { //bq
						n_el_bq++;
						if (IsRecoInEL) { //reco inel
							n_el_recoinel++;
						} //reco inel
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
					} //bq
				} //calosz
			} //pandora
		} //mid:other

		//countings -------------------------------------------------------//




	
	} //main entry loop

	//save results -------//
         
        //counting -- summary -----------------------------------------------------//
	cout<<"\nn_tot:"<<n_tot<<endl;
	cout<<"n_inel:"<<n_inel<<endl;
	cout<<"n_el:"<<n_el<<endl;
	cout<<"n_midcosmic:"<<n_midcosmic<<endl;
	cout<<"n_midpi:"<<n_midpi<<endl;
	cout<<"n_midp:"<<n_midp<<endl;
	cout<<"n_midmu:"<<n_midmu<<endl;
	cout<<"n_mideg:"<<n_mideg<<endl;
	cout<<"n_midother:"<<n_midother<<endl;
	cout<<"n_diff:"<<n_tot-(n_el+n_inel+n_midcosmic+n_midpi+n_midp+n_midmu+n_mideg+n_midother)<<endl;

	cout<<"\nn_pan_tot:"<<n_pan_tot<<endl;
	cout<<"n_inel_pan:"<<n_inel_pan<<endl;
	cout<<"n_el_pan:"<<n_el_pan<<endl;
	cout<<"n_midcosmic_pan:"<<n_midcosmic_pan<<endl;
	cout<<"n_midpi_pan:"<<n_midpi_pan<<endl;
	cout<<"n_midp_pan:"<<n_midp_pan<<endl;
	cout<<"n_midmu_pan:"<<n_midmu_pan<<endl;
	cout<<"n_mideg_pan:"<<n_mideg_pan<<endl;
	cout<<"n_midother_pan"<<n_midother_pan<<endl;
	cout<<"n_diff_pan:"<<n_pan_tot-(n_el_pan+n_inel_pan+n_midcosmic_pan+n_midpi_pan+n_midp_pan+n_midmu_pan+n_mideg_pan+n_midother_pan)<<endl;

	cout<<"\nn_calsz_tot:"<<n_calsz_tot<<endl;
	cout<<"n_inel_calsz:"<<n_inel_calsz<<endl;
	cout<<"n_el_calsz:"<<n_el_calsz<<endl;
	cout<<"n_midcosmic_calsz:"<<n_midcosmic_calsz<<endl;
	cout<<"n_midpi_calsz:"<<n_midpi_calsz<<endl;
	cout<<"n_midp_calsz:"<<n_midp_calsz<<endl;
	cout<<"n_midmu_calsz:"<<n_midmu_calsz<<endl;
	cout<<"n_mideg_calsz:"<<n_mideg_calsz<<endl;
	cout<<"n_midother_calsz"<<n_midother_calsz<<endl;
	cout<<"n_diff_calsz:"<<n_calsz_tot-(n_el_calsz+n_inel_calsz+n_midcosmic_calsz+n_midpi_calsz+n_midp_calsz+n_midmu_calsz+n_mideg_calsz+n_midother_calsz)<<endl;

	cout<<"\nn_bq_tot:"<<n_bq_tot<<endl;
	cout<<"n_inel_bq:"<<n_inel_bq<<endl;
	cout<<"n_el_bq:"<<n_el_bq<<endl;
	cout<<"n_midcosmic_bq:"<<n_midcosmic_bq<<endl;
	cout<<"n_midpi_bq:"<<n_midpi_bq<<endl;
	cout<<"n_midp_bq:"<<n_midp_bq<<endl;
	cout<<"n_midmu_bq:"<<n_midmu_bq<<endl;
	cout<<"n_mideg_bq:"<<n_mideg_bq<<endl;
	cout<<"n_midother_bq"<<n_midother_bq<<endl;
	cout<<"n_diff_bq:"<<n_bq_tot-(n_el_bq+n_inel_bq+n_midcosmic_bq+n_midpi_bq+n_midp_bq+n_midmu_bq+n_mideg_bq+n_midother_bq)<<endl;

	cout<<"\nn_recoinel_tot:"<<n_recoinel_tot<<endl;
	cout<<"n_inel_recoinel:"<<n_inel_recoinel<<endl;
	cout<<"n_el_recoinel:"<<n_el_recoinel<<endl;
	cout<<"n_midcosmic_recoinel:"<<n_midcosmic_recoinel<<endl;
	cout<<"n_midpi_recoinel:"<<n_midpi_recoinel<<endl;
	cout<<"n_midp_recoinel:"<<n_midp_recoinel<<endl;
	cout<<"n_midmu_recoinel:"<<n_midmu_recoinel<<endl;
	cout<<"n_mideg_recoinel:"<<n_mideg_recoinel<<endl;
	cout<<"n_midother_recoinel"<<n_midother_recoinel<<endl;
	cout<<"n_diff_recoinel:"<<n_recoinel_tot-(n_el_recoinel+n_inel_recoinel+n_midcosmic_recoinel+n_midpi_recoinel+n_midp_recoinel+n_midmu_recoinel+n_mideg_recoinel+n_midother_recoinel)<<endl;



        vector<int> cnt_x;
        vector<int> cnt_y_el;
        vector<int> cnt_y_inel;
        vector<int> cnt_y_midcosmic;
        vector<int> cnt_y_midpi;
        vector<int> cnt_y_midp;
        vector<int> cnt_y_midmu;
        vector<int> cnt_y_mideg;
        vector<int> cnt_y_midother;


        for (int ci=1; ci<6; ci++) {
                cnt_x.push_back(ci);
        }
        cnt_y_el.push_back(n_el);
        cnt_y_el.push_back(n_el_pan);
        cnt_y_el.push_back(n_el_calsz);
        cnt_y_el.push_back(n_el_bq);
        cnt_y_el.push_back(n_el_recoinel);

        cnt_y_inel.push_back(n_inel);
        cnt_y_inel.push_back(n_inel_pan);
        cnt_y_inel.push_back(n_inel_calsz);
        cnt_y_inel.push_back(n_inel_bq);
        cnt_y_inel.push_back(n_inel_recoinel);

        cnt_y_midcosmic.push_back(n_midcosmic);
        cnt_y_midcosmic.push_back(n_midcosmic_pan);
        cnt_y_midcosmic.push_back(n_midcosmic_calsz);
        cnt_y_midcosmic.push_back(n_midcosmic_bq);
        cnt_y_midcosmic.push_back(n_midcosmic_recoinel);

        cnt_y_midpi.push_back(n_midpi);
        cnt_y_midpi.push_back(n_midpi_pan);
        cnt_y_midpi.push_back(n_midpi_calsz);
        cnt_y_midpi.push_back(n_midpi_bq);
        cnt_y_midpi.push_back(n_midpi_recoinel);

        cnt_y_midp.push_back(n_midp);
        cnt_y_midp.push_back(n_midp_pan);
        cnt_y_midp.push_back(n_midp_calsz);
        cnt_y_midp.push_back(n_midp_bq);
        cnt_y_midp.push_back(n_midp_recoinel);

        cnt_y_midmu.push_back(n_midmu);
        cnt_y_midmu.push_back(n_midmu_pan);
        cnt_y_midmu.push_back(n_midmu_calsz);
        cnt_y_midmu.push_back(n_midmu_bq);
        cnt_y_midmu.push_back(n_midmu_recoinel);

        cnt_y_mideg.push_back(n_mideg);
        cnt_y_mideg.push_back(n_mideg_pan);
        cnt_y_mideg.push_back(n_mideg_calsz);
        cnt_y_mideg.push_back(n_mideg_bq);
        cnt_y_mideg.push_back(n_mideg_recoinel);

        cnt_y_midother.push_back(n_midother);
        cnt_y_midother.push_back(n_midother_pan);
        cnt_y_midother.push_back(n_midother_calsz);
        cnt_y_midother.push_back(n_midother_bq);
        cnt_y_midother.push_back(n_midother_recoinel);



        TGraph *Counters_el = new TGraph(cnt_x.size(), &cnt_x.at(0), &cnt_y_el.at(0));
        TGraph *Counters_inel = new TGraph(cnt_x.size(), &cnt_x.at(0), &cnt_y_inel.at(0));
        TGraph *Counters_midcosmic = new TGraph(cnt_x.size(), &cnt_x.at(0), &cnt_y_midcosmic.at(0));
        TGraph *Counters_midpi = new TGraph(cnt_x.size(), &cnt_x.at(0), &cnt_y_midpi.at(0));
        TGraph *Counters_midp = new TGraph(cnt_x.size(), &cnt_x.at(0), &cnt_y_midp.at(0));
        TGraph *Counters_midmu = new TGraph(cnt_x.size(), &cnt_x.at(0), &cnt_y_midmu.at(0));
        TGraph *Counters_mideg = new TGraph(cnt_x.size(), &cnt_x.at(0), &cnt_y_mideg.at(0));
        TGraph *Counters_midother = new TGraph(cnt_x.size(), &cnt_x.at(0), &cnt_y_midother.at(0));


        Counters_el->SetName("Counters_el");
        Counters_inel->SetName("Counters_inel");
        Counters_midcosmic->SetName("Counters_midcosmic");
        Counters_midpi->SetName("Counters_midpi");
        Counters_midp->SetName("Counters_midp");
        Counters_midmu->SetName("Counters_midmu");
        Counters_mideg->SetName("Counters_mideg");
        Counters_midother->SetName("Counters_midother");


        TFile *fout = new TFile(Form("Counters_MC.root"),"RECREATE");
		Counters_el->Write();
		Counters_inel->Write();
		Counters_midcosmic->Write();
		Counters_midpi->Write();
		Counters_midp->Write();
		Counters_midmu->Write();
		Counters_mideg->Write();
		Counters_midother->Write();
        fout->Close();







}
