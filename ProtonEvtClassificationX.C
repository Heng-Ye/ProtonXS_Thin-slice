#define ProtonEvtClassification_cxx
#include "ProtonEvtClassification.h"

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

void ProtonEvtClassification::Loop() {
	if (fChain == 0) return;

	//various counters --------------------------------------------------------------------------------------------------------------------//
	int n_tot=0; //no cut
	int n_processmap_error=0; //sansity check: processmap
	int n_true_end=0; //sansity check: true-end label

	vector<string> intTypeName{
				     "Inel",
				     "El",
				     "NoHS",
				     "Up:Inel",
				     "Up:El",
				     "Up:NoHS",
				     "Mis:Inel",
				     "Mis:El",
				     "Mis:NoHS"
		       };
	vector<string> IntType{
				  "IsPureInEL",
				  "IsPureEL",
				  "IsPureMCS",
				  "IsPureInEL",
				  "IsPureEL",
				  "IsPureMCS",
				  "IsPureInEL",
				  "IsPureEL",
				  "IsPureMCS"
		       };	
	vector<string> EvtLabel{
				  "IsBeamMatch",
				  "IsBeamMatch",
				  "IsBeamMatch",
				  "!IsBeamMatch",
				  "!IsBeamMatch",
				  "!IsBeamMatch",
				  "!IsBeamMatch",
				  "!IsBeamMatch",
				  "!IsBeamMatch"
		       };
	vector<string> IntCut;
	IntCut=IntType;
	for (size_t k=0; k<IntType.size(); ++k) {
		IntCut.at(k)+="&&"+EvtLabel.at(k);
	}
	for (size_t k=0; k<IntType.size(); ++k) {
		cout<<IntCut.at(k).c_str()<<endl;
	}
	
	vector<string> Cuts{
				"",
				"IsPandoraSlice",
				"IsPandoraSlice&&IsCaloSize",
				"IsPandoraSlice&&IsCaloSize&&IsBQ",
				"IsPandoraSlice&&IsCaloSize&&IsBQ&&IsRecoInEL"
	};
	const int n_class=9;
	const int n_cut=5;
	int count[n_class][n_cut]={0};


//str_rw_out[cnt_array].c_str()

	//int n_sg=0, n_el=0, n_inel=0, n_mcs=0;
	//int n_misid=0, n_misid_outside=0, n_misid_inside=0;
	//int n_misid_outside_el=0, n_misid_outside_inel=0, n_misid_outside_mcs=0;
	//int n_misid_inside_el=0, n_misid_inside_inel=0, n_misid_inside_mcs=0;

	//Pandora-slice cut
	//int n_pandora_sg=0, n_pandora_el=0, n_pandora_inel=0, n_pandora_mcs=0;
	//int n_pandora_misid=0, n_pandora_misid_outside=0, n_pandora_misid_inside=0;
	//int n_pandora_misid_outside_el=0, n_pandora_misid_outside_inel=0, n_pandora_misid_outside_mcs=0;
	//int n_pandora_misid_inside_el=0, n_pandora_misid_inside_inel=0, n_pandora_misid_inside_mcs=0;


	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	Long64_t nbytes = 0, nb = 0;
	bool isTestSample=true;
	for (Long64_t jentry=0; jentry<nentries;jentry++) { //main entry loop
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		isTestSample = true;
		if (ientry%2 == 0) isTestSample = false; //Divide MC sample by 2 parts: test+ufold

		//only select protons	
		if (beamtrackPdg!=pdg) continue; //only interested in protons
		//std::cout<<"beamtrackPdg:"<<beamtrackPdg<<std::endl;
		n_tot++;

		//Event Selection Cut [PartI] ----------------------------------------------------//
		bool IsBeamMatch=false; //if recostructed the right track (recoID=truthID)
		bool IsPandoraSlice=false; //pandora slice cut (can pandora reconstruct this track)
		bool IsCaloSize=false; //if calo size not empty
		bool IsCosmic=false; //beam or cosmic
		bool IsIntersection=false; //if any track intersect with our reco track		
		if (isprimarytrack==1&&isprimaryshower==0) IsPandoraSlice=true; 
		if (!primtrk_hitz->empty()) IsCaloSize=true;
		if (primary_truth_Isbeammatched==1) IsBeamMatch=true;
		if (primary_truth_byE_origin==2) IsCosmic=true;
		if (timeintersection->size()) IsIntersection=true;
		//------------------------------------------------------------------------------//

		/*
		   cout<<"\n"<<endl;
		   cout<<"run/subrun:event:"<<run<<" "<<subrun<<" "<<event<<endl;
		   cout<<"primaryID:"<<primaryID<<endl;
		   cout<<"IsPandoraSlice:"<<IsPandoraSlice<<" | isprimarytrack:"<<isprimarytrack<<" isprimaryshower:"<<isprimaryshower<<endl;
		   cout<<"IsCaloSize:"<<IsCaloSize<<endl;
		   cout<<"primary_truth_EndProcess:"<<primary_truth_EndProcess->c_str()<<endl;
		   cout<<"Isendpoint_outsidetpc:"<<Isendpoint_outsidetpc<<endl;
		   cout<<"IsBeamMatch:"<<IsBeamMatch<<endl;
		   cout<<"IsCosmic:"<<IsCosmic<<" (primary_truth_byE_origin="<<primary_truth_byE_origin<<")"<<endl;	
		   */

		//Truth label of Primarytrack_End ------------------------------------------------------------------------------------------------//
		bool IsPureInEL=false; //inel
		bool IsPureEL=false; //el
		bool IsPureMCS=false; //no hadron scattering

		if (primary_truth_EndProcess->c_str()!=NULL) n_true_end++;
		if (strcmp(primary_truth_EndProcess->c_str(),"protonInelastic")==0) {
			IsPureInEL=true;
		}
		else { //if primarytrue_end!=InEL
			if (interactionProcesslist->size()) { //size of interactionProcesslist >=0
				//cout<<"interactionProcesslist->size():"<<interactionProcesslist->size()<<endl;	
				for(size_t iiii=0; iiii<interactionProcesslist->size(); iiii++) { //loop over all true interaction hits in this track
					try {
						double intx=interactionX->at(iiii);
						double inty=interactionY->at(iiii);
						double intz=interactionZ->at(iiii); 
						//cout<<"["<<iiii<<"] process:"<<interactionProcesslist->at(iiii)<<" z:"<<intz<<endl;

						if(strcmp(interactionProcesslist->at(iiii).c_str(),"hadElastic")==0) {
							IsPureEL=1;
						}
					}
					catch (const std::out_of_range & ex) {
						//std::cout << "out_of_range Exception Caught :: interactionProcesslist" << ex.what() << std::endl;
						n_processmap_error++;
					}
				} //loop over all true interaction hits in this track 
			} //size of interactionProcesslist >=0
		} //if primarytrue_end!=InEL

		if (IsPureInEL==0&&IsPureEL==0) IsPureMCS=1;
		//cout<<"InEL EL MCS:"<<IsPureInEL<<" "<<IsPureEL<<" "<<IsPureMCS<<endl;
		//--------------------------------------------------------------------------------------------------------------------------------//

		//Get true start/end point -----------------------------------------------------------------------//
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
		//Get reco info ----------------------------------------------------------------------------------//

		//reco pos info & cut
		double reco_stx=-999, reco_sty=-999, reco_stz=-999;
		double reco_endx=-999, reco_endy=-999, reco_endz=-999;
		bool IsPos=false;
		if (IsCaloSize) {
			reco_stx=primtrk_hitx->at(0); 
			reco_sty=primtrk_hity->at(0);
			reco_stz=primtrk_hitz->at(0);

			reco_endx=primtrk_hitx->at(primtrk_dedx->size()-1);	
			reco_endy=primtrk_hity->at(primtrk_dedx->size()-1);
			reco_endz=primtrk_hitz->at(primtrk_dedx->size()-1);

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
		double x0_tmp=0, y0_tmp=0, z0_tmp=0; //start-pos, before sce
		if (primaryEndPosition[2]>primaryStartPosition[2]) { //check if Pandora flip the sign
			x0_tmp=primaryStartPosition[0];
			y0_tmp=primaryStartPosition[1];
			z0_tmp=primaryStartPosition[2];
		} //check if Pandora flip the sign
		else {
			x0_tmp=primaryEndPosition[0];
			y0_tmp=primaryEndPosition[1];
			z0_tmp=primaryEndPosition[2];
		}
		if ((pow(((x0_tmp-mean_x)/dev_x),2)+pow(((y0_tmp-mean_y)/dev_y),2))<=1.) IsXY=true;


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

		int cc=0; //no-cut
		for (int cj=0; cj<n_class; ++cj) { //loop over evt class
		} //loop over class





		for (int kt=0; kt<n_cut; ++kt) { //loop over cuts
			bool BM_cut=false;
			if (kt<=2) {


			


			for (int cj=0; cj<n_class; ++cj) { //loop over evt class
				//true label
				TString mycut = IntCut.at(cj);
				if (mycut.Data()) {
					count[cj][kj]++;
				}
			} //loop overr evt class
		} //loop over cuts


		/*
		if (IsBeamMatch) { //beam-matched track
			n_sg++; //pure signal 	
			n_pandora_sg++; //pure signal after Pandora slice cut	
			if (IsPureInEL) { //pure inel
				n_inel++; //no-cut
				if (IsPandoraSlice) { //pandora slice cut
					n_pandora_inel++; 
				} //pandora slice cut
			} //pure inel
			if (IsPureEL) { //pure el
				n_el++; //no-cut
				if (IsPandoraSlice) { //pandora slice cut
					n_pandora_el++;
				} //pandora slice cut
			} //pure el
			if (IsPureMCS) {
				n_mcs++; //no-cut
				if (IsPandoraSlice) { //pandora slice cut
					n_pandora_mcs++;
				} //pandora slice cut
			}
		} //beam-matched track
		else { //beam-not-matched track
			n_misid++;
			if (IsTrueEndOutside) { //if true_end outside tpc
				n_misid_outside++;
				n_pandora_misid_outside++;
				if (IsPureInEL) { //pure inel
					n_misid_outside_inel++;
					if (IsPandoraSlice) { //pandora slice cut
						n_pandora_misid_outside_inel++;
					} //pandora slice cut
				} //pure inel
				if (IsPureEL) { //pure el
					n_misid_outside_el++;
					if (IsPandoraSlice) { //pandora slice cut
						n_pandora_misid_outside_el++;
					} //pandora slice cut
				} //pure el
				if (IsPureMCS) { //pure mcs
					n_misid_outside_mcs++;
					if (IsPandoraSlice) { //pandora slice cut
						n_pandora_misid_outside_mcs++;
					} //pandora slice cut
				} //pure mcs
			} //if true_end outside tpc
			else { //if true_end inside tpc
				n_misid_inside++;
				if (IsPureInEL) { //pure inel
					n_misid_inside_inel++;
					if (IsPandoraSlice) { //pandora slice cut
						n_misid_inside_inel++;
					} //pandora slice cut
				} //pure inel
				if (IsPureEL) { //pure el
					n_misid_inside_el++;
				} //pure el
				if (IsPureMCS) { //pure mcs
					n_misid_inside_mcs++;
				} //pure mcs
			} //if true_End inside tpc

			
			cout<<"\nBeam-NOT-match Evt"<<endl;
			cout<<"run/subrun:event:"<<run<<" "<<subrun<<" "<<event<<endl;
			cout<<"primaryID:"<<primaryID<<endl;
			cout<<"beamtrackID:"<<beamtrackID<<endl;
			cout<<"IsPandoraSlice:"<<IsPandoraSlice<<" | isprimarytrack:"<<isprimarytrack<<" isprimaryshower:"<<isprimaryshower<<endl;
			cout<<"IsCaloSize:"<<IsCaloSize<<endl;
			cout<<"Primary_truth_EndProcess:"<<primary_truth_EndProcess->c_str()<<endl;
			cout<<"IsEndpoint_Outsidetpc:"<<Isendpoint_outsidetpc<<endl;
			cout<<"IsBeamMatch:"<<IsBeamMatch<<endl;
			cout<<"IsCosmic:"<<IsCosmic<<" (primary_truth_byE_origin="<<primary_truth_byE_origin<<")"<<endl;
			cout<<"InEL EL MCS:"<<IsPureInEL<<" "<<IsPureEL<<" "<<IsPureMCS<<endl;
			cout<<"Primary true_endz:"<<true_endz<<endl;	
			cout<<"Primary reco_endz:"<<reco_endz<<endl;	
			cout<<"NDAUGHTERS:"<<NDAUGHTERS<<endl;
			if (NDAUGHTERS>0) {
				for (size_t jj=0; jj<(size_t)NDAUGHTERS; ++jj) {
					cout<<"   daughter_truth_TrackId["<<jj<<"]: "<<daughter_truth_TrackId[jj]<<endl;
					cout<<"   daughter_truth_Pdg["<<jj<<"]: "<<daughter_truth_Pdg[jj]<<endl;
					cout<<"   daughter_truth_StartPositionZ["<<jj<<"]: "<<daughter_truth_StartPosition[jj][2]<<endl;
				}
			}

		} //beam-not-matched track
			*/

} //main entry loop

//counting -- summary -----------------------------------------------------//
/*
cout<<"\nAll proton counting..."<<endl;
cout<<"n_tot:"<<n_tot<<endl;
cout<<"  n_true_end:"<<n_true_end<<endl;
cout<<"  n_processmap_error:"<<n_processmap_error<<endl;
cout<<"  n_sg:"<<n_sg<<endl;
cout<<"  n_inel:"<<n_inel<<endl;
cout<<"  n_el:"<<n_el<<endl;
cout<<"  n_mcs:"<<n_mcs<<endl;
cout<<"  n_inel+n_el+n_mcs:"<<n_inel+n_el+n_mcs<<endl;
cout<<"  n_misid:"<<n_misid<<endl;
cout<<"     n_misid_inside:"<<n_misid_inside<<endl;
cout<<"         n_misid_inside_inel:"<<n_misid_inside_inel<<endl;
cout<<"         n_misid_inside_el:"<<n_misid_inside_el<<endl;
cout<<"         n_misid_inside_mcs:"<<n_misid_inside_mcs<<endl;
cout<<"         n_misid_inside_inel+n_misid_inside_el+n_misid_inside_mcs:"<<n_misid_inside_inel+n_misid_inside_el+n_misid_inside_mcs<<endl;
cout<<"     n_misid_outside:"<<n_misid_outside<<endl;
cout<<"         n_misid_outside_inel:"<<n_misid_outside_inel<<endl;
cout<<"         n_misid_outside_el:"<<n_misid_outside_el<<endl;
cout<<"         n_misid_outside_mcs:"<<n_misid_outside_mcs<<endl;
cout<<"         n_misid_outside_inel+n_misid_outside_el+n_misid_outside_mcs:"<<n_misid_outside_inel+n_misid_outside_el+n_misid_outside_mcs<<endl;
cout<<"     n_misid_inside+n_misid_outside:"<<n_misid_inside+n_misid_outside<<endl;	
cout<<"  n_sg+n_misid:"<<n_sg+n_misid<<endl;
*/


		for (int c=0; c<n_cut; ++c) { //cut loop
			for (int j=0; j<n_class; ++j) { //class loop
			//int c=0;
				//string evt_selection=IntCut.at(j).c_str()+"&&"+Cuts.at(c).c_str();
				//cout<<evt_selection.c_str()<<" | cnt:"<<count[j][c]<<endl;
				cout<<" cnt:"<<count[j][c]<<endl;
			} //cut loop
			cout<<"\n"<<endl;
		} //class loop	



}
