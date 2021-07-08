#define ProtonThinSlice_cxx
#include "ProtonThinSlice.h"

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

void ProtonThinSlice::Loop() {
	if (fChain == 0) return;

	//various counters --------------------------------------------------------------------------------------------------------------------//
	int n_processmap_error=0; //sansity check: processmap
	int n_true_end=0; //sansity check: true-end label

	int n_tot=0; //no cut
	int n_sg=0, n_bkg=0, n_up=0, n_mis=0;
	int n_el=0, n_inel=0, n_mcs=0, n_upel=0, n_upinel=0, n_upmcs=0, n_misel=0, n_misinel=0, n_mismcs=0;
	int n_pan=0; //pandora cut
	int n_el_pan=0, n_inel_pan=0, n_mcs_pan=0, n_upel_pan=0, n_upinel_pan=0, n_upmcs_pan=0, n_misel_pan=0, n_misinel_pan=0, n_mismcs_pan=0;
	int n_calsz=0; //calosz cut
	int n_el_calsz=0, n_inel_calsz=0, n_mcs_calsz=0, n_upel_calsz=0, n_upinel_calsz=0, n_upmcs_calsz=0, n_misel_calsz=0, n_misinel_calsz=0, n_mismcs_calsz=0;
	int n_bq=0; //bq cut
	int n_el_bq=0, n_inel_bq=0, n_mcs_bq=0, n_upel_bq=0, n_upinel_bq=0, n_upmcs_bq=0, n_misel_bq=0, n_misinel_bq=0, n_mismcs_bq=0;
	int n_recoinel=0; //reco inel cut
	int n_el_recoinel=0, n_inel_recoinel=0, n_mcs_recoinel=0, n_upel_recoinel=0, n_upinel_recoinel=0, n_upmcs_recoinel=0, n_misel_recoinel=0, n_misinel_recoinel=0, n_mismcs_recoinel=0;

	//int n_misid=0, n_misid_intoutside=0, n_misid_inttpc=0, n_pandora_misid=0, n_calosz_misid=0, n_bq_misid=0, n_recoinel_misid=0;
	//int n_pandora=0, n_pandora_el=0, n_pandora_inel=0, n_pandora_mcs=0, n_pandora_mid=0; //pandora slice
	//int n_calosz=0, n_calosz_el=0, n_calosz_inel=0, n_calosz_mcs=0, n_calosz_mid=0; //calo size
	//int n_bq=0, n_bq_el=0, n_bq_inel=0, n_bq_mcs=0, n_bq_mid=0; //calo size
	//int n_recoinel=0, n_recoinel_el=0, n_recoinel_inel=0, n_recoinel_mcs=0, n_recoinel_mid=0; //recoinel
	
	//unfolding config. -------------------------//
	Unfold uf(nthinslices+2, -1, nthinslices+1);

	//ThinSlice config. ---------------------------------------------------------------------------------------------------//
	SetOutputFileName(Form("prod4a_thinslice_dx%dcm_%dslcs.root", name_thinslicewidth, nthinslices)); //output file name

	//book histograms --//
	BookHistograms();


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
		bool IsCosmic=false; //beam or cosmic
		bool IsIntersection=false; //if any track intersect with our reco track		
		if (primary_truth_Isbeammatched==1) IsBeamMatch=true;
		if (isprimarytrack==1&&isprimaryshower==0) IsPandoraSlice=true; 
		if (!primtrk_hitz->empty()) IsCaloSize=true;
		if (primary_truth_byE_origin==2) IsCosmic=true;
		if (timeintersection->size()) IsIntersection=true;
		//----------------------------------------------------------------//
	
		cout<<"\n"<<endl;
		cout<<"run/subrun:event:"<<run<<" "<<subrun<<" "<<event<<endl;
		cout<<"primaryID:"<<primaryID<<endl;
		cout<<"IsPandoraSlice:"<<IsPandoraSlice<<" | isprimarytrack:"<<isprimarytrack<<" isprimaryshower:"<<isprimaryshower<<endl;
		cout<<"IsCaloSize:"<<IsCaloSize<<endl;
		cout<<"primary_truth_EndProcess:"<<primary_truth_EndProcess->c_str()<<endl;
		cout<<"Isendpoint_outsidetpc:"<<Isendpoint_outsidetpc<<endl;
		cout<<"IsBeamMatch:"<<IsBeamMatch<<endl;
		cout<<"IsCosmic:"<<IsCosmic<<" (primary_truth_byE_origin="<<primary_truth_byE_origin<<")"<<endl;
		cout<<"primary_truth_Pdg:"<<primary_truth_Pdg<<endl;
		cout<<"beamtrackPdg:"<<beamtrackPdg<<endl;	

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

		} //if primarytrue_end!=InEL

		if (IsPureInEL==0&&IsPureEL==0) {
			IsPureMCS=1;
		}

		//if (strcmp(primary_truth_EndProcess->c_str(),"hIoni")==0) {
			//IsPureEL=true;
		//}
		//if (strcmp(primary_truth_EndProcess->c_str(),"CoulombScat")==0) {
			//IsPureMCS=true;
		//}
		cout<<"InEL EL MCS:"<<IsPureInEL<<" "<<IsPureEL<<" "<<IsPureMCS<<endl;
		//--------------------------------------------------------------------------------------------------------------------------------//

		for (size_t j=0; j<beamtrk_z->size(); ++j) { //MCParticle loop
			cout<<"beamtrk_z["<<j<<"]"<<beamtrk_z->at(j)<<" beamtrk_Eng["<<"]"<<beamtrk_Eng->at(j)<<endl;
		} //MCParticle loop

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
		cout<<"trueEnd z/y/x:"<<true_endz<<"/"<<true_endy<<"/"<<true_endx<<endl;
		cout<<"trueSt z/y/x:"<<true_stz<<"/"<<true_sty<<"/"<<true_stx<<endl;
		cout<<"IsTrueEndOutside:"<<IsTrueEndOutside<<endl;
		//Get reco info ----------------------------------------------------------------------------------//

		//reco pos info & cut
		double reco_stx=-999, reco_sty=-999, reco_stz=-999;
		double reco_endx=-999, reco_endy=-999, reco_endz=-999;
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
		double ke_trklen=1000.*ke_vs_csda_range_sm->Eval(range_reco); //[unit: MeV]
		double p_trklen=ke2p(ke_trklen);
		double ke_simide=0;
		for (int hk=0; hk<(int)primtrk_true_edept->size(); ++hk) { //loop over simIDE points
			ke_simide+=primtrk_true_edept->at(hk);
		} //loop over simIDE points

		//reco calorimetry
		vector< pair<double,int > > zreco_rawindex; //z, original_index
		vector< pair<double,double > > zreco_rr; //z, rr
		vector< pair<double,double > > zreco_de; //z, wid_reco
		vector< pair<double,double > > zreco_dedx; //z, wid_reco
		vector< pair<double,double > > zreco_dx; //z, wid_reco
		//vector< pair<double,double > > zreco_ke; //z, wid_reco
		vector< pair<double,double > > zreco_xreco; //z, wid_reco
		vector< pair<double,double > > zreco_yreco; //z, wid_reco
		vector< pair<double,int > > zreco_widreco; //z, wid_reco

		int index_reco_endz=0;
		double wid_reco_max=-9999;
		double kereco_calo=0;
		double kereco_range=0;
		double kereco_range2=0;
		if (IsCaloSize) { //if calo size not empty
		  for (size_t h=0; h<primtrk_dedx->size(); ++h) { //loop over reco hits of a given track
			double hitx_reco=primtrk_hitx->at(h);
			double hity_reco=primtrk_hity->at(h);
			double hitz_reco=primtrk_hitz->at(h);
			double resrange_reco=primtrk_resrange->at(h);

			double dqdx=primtrk_dqdx->at(h);
			//double resrange=primtrk_resrange->at(h);
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
			kereco_calo+=cali_dedx*pitch;
			//ke_reco-=cali_dedx*pitch;

			//use dedx from rr
			kereco_range+=pitch*dedx_predict(resrange_reco);
			kereco_range2+=pitch*(double)gr_predict_dedx_resrange->Eval(resrange_reco);

			if (IsRecoStop) {
				//double theory1_dedx=(double)gr_predict_dedx_resrange->Eval(resrange_reco);
				//double theory2_dedx=dedx_predict(resrange_reco);

				rr_dedx_recostop->Fill(resrange_reco,cali_dedx);

				//rr_ddedx1_recostop->Fill(resrange_reco,
			}

			if (IsPureMCS) {
				rr_dedx_truestop->Fill(resrange_reco,cali_dedx);
			}

			//if (IsPureInELProton) cout<<"["<<h<<"] z:"<<hitz_reco<<" dx:"<<pitch<<" wid:"<<wid_reco<<" dedx:"<<cali_dedx<<endl;
			//if (run==39279896&&event==1165) cout<<"["<<h<<"] z:"<<hitz_reco<<" dx:"<<pitch<<" wid:"<<wid_reco<<" dedx:"<<cali_dedx<<endl;
			//if (run==39279896&&event==1165) cout<<"["<<h<<"] z:"<<hitz_reco<<" dx:"<<pitch<<" wid:"<<wid_reco<<" pt:"<<pt_reco<<" dedx:"<<cali_dedx<<" rr:"<<resrange_reco<<endl;
			//if(IsPureMCS) cout<<"["<<h<<"] z:"<<hitz_reco<<" dx:"<<pitch<<" wid:"<<wid_reco<<" dedx:"<<cali_dedx<<endl;

			zreco_rawindex.push_back(make_pair(wid_reco, h));
			zreco_de.push_back(make_pair(wid_reco, cali_dedx*pitch));
			zreco_dedx.push_back(make_pair(wid_reco, cali_dedx));
			zreco_dx.push_back(make_pair(wid_reco, pitch));
			zreco_xreco.push_back(make_pair(wid_reco, hitx_reco));
			zreco_yreco.push_back(make_pair(wid_reco, hity_reco));
			zreco_widreco.push_back(make_pair(wid_reco, hitz_reco));
			zreco_rr.push_back(make_pair(wid_reco, resrange_reco));

			if (IsPureInEL) rangereco_dedxreco_TrueInEL->Fill(range_reco-resrange_reco, cali_dedx);
			if (IsPureEL) rangereco_dedxreco_TrueEL->Fill(range_reco-resrange_reco, cali_dedx);
			if (IsPureMCS) rangereco_dedxreco_TrueMCS->Fill(range_reco-resrange_reco, cali_dedx);

			//zreco_ke.push_back(make_pair(wid_reco, ke_reco));
		  } //loop over reco hits of a given track
		} //if calo size not empty

		if (IsRecoStop) { 
			KE_ff_recostop->Fill(ke_ff);
			KE_calo_recostop->Fill(kereco_calo);
			KE_rrange_recostop->Fill(kereco_range);
			KE_rrange2_recostop->Fill(kereco_range2);
			KE_range_recostop->Fill(ke_trklen);
			KE_simide_recostop->Fill(ke_simide);

			dKE_range_ff_recostop->Fill(ke_trklen-ke_ff);
			dKE_calo_ff_recostop->Fill(kereco_calo-ke_ff);
			dKE_rrange_ff_recostop->Fill(kereco_range-ke_ff);
			dKE_rrange2_ff_recostop->Fill(kereco_range2-ke_ff);

			KE_range_ff_recostop->Fill(kereco_range, ke_ff);
			KE_range_calo_recostop->Fill(kereco_range, kereco_calo);

		}
		//cout<<"\n"<<endl;

		//sort using wid
		sort(zreco_rawindex.begin(),zreco_rawindex.end(),myComparison); //sorting based on the first column
		sort(zreco_rr.begin(),zreco_rr.end(),myComparison); //sorting based on the first column
		sort(zreco_de.begin(),zreco_de.end(),myComparison); //sorting based on the first column
		sort(zreco_dedx.begin(),zreco_dedx.end(),myComparison); //sorting based on the first column
		sort(zreco_dx.begin(),zreco_dx.end(),myComparison); //sorting based on the first column
		sort(zreco_xreco.begin(),zreco_xreco.end(),myComparison); //sorting based on the first column
		sort(zreco_yreco.begin(),zreco_yreco.end(),myComparison); //sorting based on the first column
		sort(zreco_widreco.begin(),zreco_widreco.end(),myComparison); //sorting based on the first column
		//sort(zreco_ke.begin(),zreco_ke.end(),myComparison); //sorting based on the first column

		//countings -------------------------------------------------------//
		if (IsPandoraSlice) n_pan++;
		if (IsPandoraSlice&&IsCaloSize) n_calsz++;
		if (IsPandoraSlice&&IsCaloSize&&IsBQ) n_bq++;
		if (IsPandoraSlice&&IsCaloSize&&IsBQ&&IsRecoInEL) n_recoinel++;

		if (IsBeamMatch) n_sg++;
		if (!IsBeamMatch) n_bkg++;
		if (!IsBeamMatch&&IsTrueEndOutside) n_up++;
		if (!IsBeamMatch&&!IsTrueEndOutside) n_mis++;


		//[0]pure inel
		if (IsBeamMatch&&IsPureInEL) { //pure inel
			n_inel++;
			zend_true_inel_NoCut->Fill(true_endz);
			zend_reco_inel_NoCut->Fill(reco_endz);
			dzend_inel_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_inel_pan++;
				zend_true_inel_PanS->Fill(true_endz);
				zend_reco_inel_PanS->Fill(reco_endz);
				dzend_inel_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_inel_calsz++;
					zend_true_inel_CaloSz->Fill(true_endz);
					zend_reco_inel_CaloSz->Fill(reco_endz);
					dzend_inel_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_inel_bq++;
						zend_true_inel_BQ->Fill(true_endz);
						zend_reco_inel_BQ->Fill(reco_endz);
						dzend_inel_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_inel_recoinel++;
							zend_true_inel_RecoInel->Fill(true_endz);
							zend_reco_inel_RecoInel->Fill(reco_endz);
							dzend_inel_RecoInel->Fill(reco_endz-true_endz);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //pure inel		

		//[1]pure el
		if (IsBeamMatch&&IsPureEL) { //pure el
			n_el++;
			zend_true_el_NoCut->Fill(true_endz);
			zend_reco_el_NoCut->Fill(reco_endz);
			dzend_el_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_el_pan++;
				zend_true_el_PanS->Fill(true_endz);
				zend_reco_el_PanS->Fill(reco_endz);
				dzend_el_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_el_calsz++;
					zend_true_el_CaloSz->Fill(true_endz);
					zend_reco_el_CaloSz->Fill(reco_endz);
					dzend_el_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_el_bq++;
						zend_true_el_BQ->Fill(true_endz);
						zend_reco_el_BQ->Fill(reco_endz);
						dzend_el_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_el_recoinel++;
							zend_true_el_RecoInel->Fill(true_endz);
							zend_reco_el_RecoInel->Fill(reco_endz);
							dzend_el_RecoInel->Fill(reco_endz-true_endz);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //pure el		

		//[2]pure mcs
		if (IsBeamMatch&&IsPureMCS) { //pure mcs
			n_mcs++;
			zend_true_mcs_NoCut->Fill(true_endz);
			zend_reco_mcs_NoCut->Fill(reco_endz);
			dzend_mcs_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_mcs_pan++;
				zend_true_mcs_PanS->Fill(true_endz);
				zend_reco_mcs_PanS->Fill(reco_endz);
				dzend_mcs_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_mcs_calsz++;
					zend_true_mcs_CaloSz->Fill(true_endz);
					zend_reco_mcs_CaloSz->Fill(reco_endz);
					dzend_mcs_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_mcs_bq++;
						zend_true_mcs_BQ->Fill(true_endz);
						zend_reco_mcs_BQ->Fill(reco_endz);
						dzend_mcs_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_mcs_recoinel++;
							zend_true_mcs_RecoInel->Fill(true_endz);
							zend_reco_mcs_RecoInel->Fill(reco_endz);
							dzend_mcs_RecoInel->Fill(reco_endz-true_endz);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //pure mcs

		//[3]pure inel:upstream-int
		if (!IsBeamMatch&&IsTrueEndOutside&&IsPureInEL) { //pure inel:upstream-int
			n_upinel++;
			zend_true_upinel_NoCut->Fill(true_endz);
			zend_reco_upinel_NoCut->Fill(reco_endz);
			dzend_upinel_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_upinel_pan++;
				zend_true_upinel_PanS->Fill(true_endz);
				zend_reco_upinel_PanS->Fill(reco_endz);
				dzend_upinel_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_upinel_calsz++;
					zend_true_upinel_CaloSz->Fill(true_endz);
					zend_reco_upinel_CaloSz->Fill(reco_endz);
					dzend_upinel_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_upinel_bq++;
						zend_true_upinel_BQ->Fill(true_endz);
						zend_reco_upinel_BQ->Fill(reco_endz);
						dzend_upinel_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_upinel_recoinel++;
							zend_true_upinel_RecoInel->Fill(true_endz);
							zend_reco_upinel_RecoInel->Fill(reco_endz);
							dzend_upinel_RecoInel->Fill(reco_endz-true_endz);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //pure inel:upstream-int		

		//[4]pure el:upstream-int
		if (!IsBeamMatch&&IsTrueEndOutside&&IsPureEL) { //pure el:upstream-int
			n_upel++;
			zend_true_upel_NoCut->Fill(true_endz);
			zend_reco_upel_NoCut->Fill(reco_endz);
			dzend_upel_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_upel_pan++;
				zend_true_upel_PanS->Fill(true_endz);
				zend_reco_upel_PanS->Fill(reco_endz);
				dzend_upel_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_upel_calsz++;
					zend_true_upel_CaloSz->Fill(true_endz);
					zend_reco_upel_CaloSz->Fill(reco_endz);
					dzend_upel_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_upel_bq++;
						zend_true_upel_BQ->Fill(true_endz);
						zend_reco_upel_BQ->Fill(reco_endz);
						dzend_upel_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_upel_recoinel++;
							zend_true_upel_RecoInel->Fill(true_endz);
							zend_reco_upel_RecoInel->Fill(reco_endz);
							dzend_upel_RecoInel->Fill(reco_endz-true_endz);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //pure el:upstream-int	

		//[5]pure mcs:upstream-int
		if (!IsBeamMatch&&IsTrueEndOutside&&IsPureMCS) { //pure mcs:upstream-int
			n_upmcs++;
			zend_true_upmcs_NoCut->Fill(true_endz);
			zend_reco_upmcs_NoCut->Fill(reco_endz);
			dzend_upmcs_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_upmcs_pan++;
				zend_true_upmcs_PanS->Fill(true_endz);
				zend_reco_upmcs_PanS->Fill(reco_endz);
				dzend_upmcs_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_upmcs_calsz++;
					zend_true_upmcs_CaloSz->Fill(true_endz);
					zend_reco_upmcs_CaloSz->Fill(reco_endz);
					dzend_upmcs_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_upmcs_bq++;
						zend_true_upmcs_BQ->Fill(true_endz);
						zend_reco_upmcs_BQ->Fill(reco_endz);
						dzend_upmcs_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_upmcs_recoinel++;
							zend_true_upmcs_RecoInel->Fill(true_endz);
							zend_reco_upmcs_RecoInel->Fill(reco_endz);
							dzend_upmcs_RecoInel->Fill(reco_endz-true_endz);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //pure mcs:upstream-int	

		//[6]pure inel+int. inside tpc
		if (!IsBeamMatch&&!IsTrueEndOutside&&IsPureInEL) { //pure inel
			n_misinel++;
			zend_true_misinel_NoCut->Fill(true_endz);
			zend_reco_misinel_NoCut->Fill(reco_endz);
			dzend_misinel_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_misinel_pan++;
				zend_true_misinel_PanS->Fill(true_endz);
				zend_reco_misinel_PanS->Fill(reco_endz);
				dzend_misinel_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_misinel_calsz++;
					zend_true_misinel_CaloSz->Fill(true_endz);
					zend_reco_misinel_CaloSz->Fill(reco_endz);
					dzend_misinel_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_misinel_bq++;
						zend_true_misinel_BQ->Fill(true_endz);
						zend_reco_misinel_BQ->Fill(reco_endz);
						dzend_misinel_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_misinel_recoinel++;
							zend_true_misinel_RecoInel->Fill(true_endz);
							zend_reco_misinel_RecoInel->Fill(reco_endz);
							dzend_misinel_RecoInel->Fill(reco_endz-true_endz);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //pure inel		

		//[7]pure el+int. inside tpc+misid
		if (!IsBeamMatch&&!IsTrueEndOutside&&IsPureEL) { //pure el
			n_misel++;
			zend_true_misel_NoCut->Fill(true_endz);
			zend_reco_misel_NoCut->Fill(reco_endz);
			dzend_misel_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_misel_pan++;
				zend_true_misel_PanS->Fill(true_endz);
				zend_reco_misel_PanS->Fill(reco_endz);
				dzend_misel_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_misel_calsz++;
					zend_true_misel_CaloSz->Fill(true_endz);
					zend_reco_misel_CaloSz->Fill(reco_endz);
					dzend_misel_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_misel_bq++;
						zend_true_misel_BQ->Fill(true_endz);
						zend_reco_misel_BQ->Fill(reco_endz);
						dzend_misel_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_misel_recoinel++;
							zend_true_misel_RecoInel->Fill(true_endz);
							zend_reco_misel_RecoInel->Fill(reco_endz);
							dzend_misel_RecoInel->Fill(reco_endz-true_endz);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //pure el		

		//[8]pure mcs+int. inside tpc+misid
		if (!IsBeamMatch&&!IsTrueEndOutside&&IsPureMCS) { //pure mcs
			n_mismcs++;
			zend_true_mismcs_NoCut->Fill(true_endz);
			zend_reco_mismcs_NoCut->Fill(reco_endz);
			dzend_mismcs_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_mismcs_pan++;
				zend_true_mismcs_PanS->Fill(true_endz);
				zend_reco_mismcs_PanS->Fill(reco_endz);
				dzend_mismcs_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_mismcs_calsz++;
					zend_true_mismcs_CaloSz->Fill(true_endz);
					zend_reco_mismcs_CaloSz->Fill(reco_endz);
					dzend_mismcs_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_mismcs_bq++;
						zend_true_mismcs_BQ->Fill(true_endz);
						zend_reco_mismcs_BQ->Fill(reco_endz);
						dzend_mismcs_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_mismcs_recoinel++;
							zend_true_mismcs_RecoInel->Fill(true_endz);
							zend_reco_mismcs_RecoInel->Fill(reco_endz);
							dzend_mismcs_RecoInel->Fill(reco_endz-true_endz);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //pure mcs
		//countings -------------------------------------------------------//


		//thin-slice method ----------------------------------------------------//
		//true slice ID
		true_sliceID = int(true_endz/thinslicewidth);
		if (true_sliceID < 0) true_sliceID = -1;
		if (true_endz < 0) true_sliceID = -1; //hy added
		if (true_sliceID >= nthinslices) true_sliceID = nthinslices;
		cout<<"true_endz:"<<true_endz<<"  true_sliceID:"<<true_sliceID<<endl;

		//evt selection cuts
		bool PassCuts_INT=false; //all bq cut+reco inel cut
		bool PassCuts_INC=false; //all bq cut
		if (IsPandoraSlice&&IsCaloSize) {
			//reco slice ID
			reco_sliceID = int(reco_endz/thinslicewidth);
			if (reco_sliceID < 0) reco_sliceID = -1;
			if (reco_endz<0) reco_sliceID = -1;
			if (reco_sliceID >= nthinslices) reco_sliceID = nthinslices;
			cout<<"reco_endz:"<<reco_endz<<"  reco_sliceID:"<<reco_sliceID<<endl;

			if (IsBQ&&IsRecoInEL) {
				PassCuts_INT=true; //for INT 
			}
			if (IsBQ) {
				PassCuts_INC=true; //for INC
			}
		}

		//INC histograms ------------------------------------------------------------//
		for (int ij = 0; ij<=true_sliceID; ++ij){
			if (ij<nthinslices) ++true_incidents[ij];
		}

		if (isTestSample) { //if test sample 
			h_truesliceid_all->Fill(true_sliceID);
		} //if test sample
		else { //if NOT test sample
			uf.eff_den_Inc->Fill(true_sliceID);
		} //if NOT test sample
		if (PassCuts_INC&&IsBeamMatch) { //if passing all basic cuts
			if (isTestSample) { //if test sample
				h_recosliceid_cuts->Fill(reco_sliceID);
				h_truesliceid_cuts->Fill(true_sliceID);
			} //if test sample
			else{ //if NOT test sample
				uf.eff_num_Inc->Fill(true_sliceID);
				uf.pur_num_Inc->Fill(reco_sliceID);
				uf.response_SliceID_Inc.Fill(reco_sliceID, true_sliceID);
			} //if NOT test sample
		} //if passing all basic cuts
		else { //if NOT passing all cuts
			if (!isTestSample){
				uf.response_SliceID_Inc.Miss(true_sliceID);
				//std::cout<<true_sliceID<<std::endl;
			}
		} //if NOT passing all cuts

		//INT histograms ------------------------------------------------------------//
		if (IsPureInEL) { //pure inel
			if (isTestSample){
				h_truesliceid_inelastic_all->Fill(true_sliceID);
			}
			else{ //NOT test sample for unfolding
				uf.eff_den_Int->Fill(true_sliceID);
			} //NOT test sample for unfolding

			if (PassCuts_INT&&IsBeamMatch) { //if pass all basic cuts
				if (isTestSample){
					h_recosliceid_inelastic_cuts->Fill(reco_sliceID);
					h_truesliceid_inelastic_cuts->Fill(true_sliceID);
				}
				else{
					uf.eff_num_Int->Fill(true_sliceID);
					uf.pur_num_Int->Fill(reco_sliceID);
					uf.response_SliceID_Int.Fill(reco_sliceID, true_sliceID);
				}
			} //if pass all basic cuts
			else{ //if NOT pass all basic cuts
				if (!isTestSample) uf.response_SliceID_Int.Miss(true_sliceID);
			} //if not pass all basic cuts
		} //pure inel
		
		if (PassCuts_INC) { //if pass all cuts
			if (isTestSample){ //if test sample
				h_recosliceid_allevts_cuts->Fill(reco_sliceID);
			} //if test sample
			else { //if NOT test sample
				uf.pur_den->Fill(reco_sliceID);
			} //if NOT test sample
		} //if pass all cuts

		//reco/truth KEs
		if (IsPureInEL) { //is pure inelastic
			if (true_sliceID < nthinslices && true_sliceID>=0){
				++true_interactions[true_sliceID];
			}
			//reco ke
			if (IsCaloSize==true&&IsBeamMatch==true) { //calo & beam_match
				std::vector<std::vector<double>> vincE(nthinslices);
				double ke_reco=ke_ff;
				for (size_t ih=0; ih<zreco_rawindex.size(); ++ih) {
					double this_calo_z=zreco_widreco[ih].second;
					int this_sliceID = int(this_calo_z/thinslicewidth);
					double this_dE=zreco_de[ih].second;
					ke_reco-=this_dE;

					if (this_sliceID>=nthinslices) continue;
					if (this_sliceID<0) continue;

					double this_incE = ke_reco;
					vincE[this_sliceID].push_back(this_incE);
				}

				for (size_t i = 0; i<vincE.size(); ++i){
					if (!vincE[i].empty()){
						double sum_incE = 0;
						for (size_t j = 0; j<vincE[i].size(); ++j){
							sum_incE += vincE[i][j];
						}
						reco_incE[i]->Fill(sum_incE/vincE[i].size());
					}
				}

				TVector3 pt0(zreco_xreco[0].second,
						zreco_yreco[0].second,
						zreco_xreco[0].first); //ST
				TVector3 pt1(zreco_xreco[-1+zreco_rawindex.size()].second,
						zreco_yreco[-1+zreco_rawindex.size()].second,
						zreco_xreco[-1+zreco_rawindex.size()].first);
				TVector3 dir = pt1 - pt0;
				dir = dir.Unit();
				reco_AngCorr->Fill(dir.Z());
			} //calo & beam_match

			//truth kes
			if (!beamtrk_Eng->empty()) { //if true container not empty
				std::vector<std::vector<double>> vincE_true(nthinslices);
				for (int hk=0; hk<(int)beamtrk_Eng->size()-1; ++hk) { //loop over true hits
					double thisZ=beamtrk_z->at(hk);
					int this_sliceID = int(thisZ/thinslicewidth);
					double this_incE = 1000.*beamtrk_Eng->at(hk); //MeV

					if (this_sliceID>=nthinslices) continue;
					if (this_sliceID<0) continue;
					vincE_true[this_sliceID].push_back(this_incE);
				}//loop over true hits (last point always has KE = 0)

				for (size_t i = 0; i<vincE_true.size(); ++i){
					if (!vincE_true[i].empty()){
						double sum_incE_true = 0;
						for (size_t j = 0; j<vincE_true[i].size(); ++j){
							sum_incE_true += vincE_true[i][j];
						}
						true_incE[i]->Fill(sum_incE_true/vincE_true[i].size());
					}
				}

				TVector3 pt0(true_stx, true_sty, true_stz);
				TVector3 pt1(true_endx, true_endy, true_endz);
				TVector3 dir = pt1 - pt0;
				dir = dir.Unit();
				true_AngCorr->Fill(dir.Z());
			} //if true container not empty
		} //if pure inelastic


		if (!beamtrk_Eng->empty()) { //if true container not empty
			for (int hk=0; hk<(int)beamtrk_Eng->size()-1; ++hk) { //loop over true hits
				double thisZ=beamtrk_z->at(hk);
				double this_incE = 1000.*beamtrk_Eng->at(hk); //MeV

				if (IsPureMCS) ztrue_ketrue_TrueMCS->Fill(thisZ, this_incE);
				if (IsPureEL) ztrue_ketrue_TrueEL->Fill(thisZ, this_incE);
				if (IsPureInEL) ztrue_ketrue_TrueInEL->Fill(thisZ, this_incE);
			} //loop over true hits
		} //if true container not empty

		if (IsPandoraSlice&&IsCaloSize&&IsBQ) {  
			ntrklen->Fill(range_reco/csda_val_spec);
			if (IsPureInEL) ntrklen_PureInel->Fill(range_reco/csda_val_spec); 
			if (IsPureEL) ntrklen_PureEl->Fill(range_reco/csda_val_spec); 
			if (IsPureMCS) ntrklen_PureMCS->Fill(range_reco/csda_val_spec); 
		}


	
	} //main entry loop

	//save results -------//
	uf.SaveHistograms();
	CalcXS(uf);
	SaveHistograms();
         
        //counting -- summary -----------------------------------------------------//
	cout<<"\nn_tot:"<<n_tot<<endl;
	cout<<"n_sg:"<<n_sg<<endl;
	cout<<"n_el:"<<n_el<<endl;
	cout<<"n_inel:"<<n_inel<<endl;
	cout<<"n_mcs:"<<n_mcs<<endl;
	cout<<"n_bkg:"<<n_bkg<<endl;
	cout<<"n_up:"<<n_up<<endl;
	cout<<"n_upel:"<<n_upel<<endl;
	cout<<"n_upinel:"<<n_upinel<<endl;
	cout<<"n_upmcs:"<<n_upmcs<<endl;
	cout<<"n_mis:"<<n_mis<<endl;
	cout<<"n_misel:"<<n_misel<<endl;
	cout<<"n_misinel:"<<n_misinel<<endl;
	cout<<"n_mismcs:"<<n_mismcs<<endl;
	cout<<"diff:"<<n_tot-n_el-n_inel-n_mcs-n_upel-n_upinel-n_upmcs-n_misel-n_misinel-n_mismcs<<endl;
	cout<<"\n"<<endl;

	
	cout<<"n_pan:"<<n_pan<<endl;
	cout<<"n_el_pan:"<<n_el_pan<<endl;
	cout<<"n_inel_pan:"<<n_inel_pan<<endl;
	cout<<"n_mcs_pan:"<<n_mcs_pan<<endl;
	cout<<"n_upel_pan:"<<n_upel_pan<<endl;
	cout<<"n_upinel_pan:"<<n_upinel_pan<<endl;
	cout<<"n_upmcs_pan:"<<n_upmcs_pan<<endl;
	cout<<"n_misel_pan:"<<n_misel_pan<<endl;
	cout<<"n_misinel_pan:"<<n_misinel_pan<<endl;
	cout<<"n_mismcs_pan:"<<n_mismcs_pan<<endl;
	cout<<"diff:"<<n_pan-n_el_pan-n_inel_pan-n_mcs_pan-n_upel_pan-n_upinel_pan-n_upmcs_pan-n_misel_pan-n_misinel_pan-n_mismcs_pan<<endl;
	cout<<"\n"<<endl;

	cout<<"n_calsz:"<<n_calsz<<endl;
	cout<<"n_el_calsz:"<<n_el_calsz<<endl;
	cout<<"n_inel_calsz:"<<n_inel_calsz<<endl;
	cout<<"n_mcs_calsz:"<<n_mcs_calsz<<endl;
	cout<<"n_upel_calsz:"<<n_upel_calsz<<endl;
	cout<<"n_upinel_calsz:"<<n_upinel_calsz<<endl;
	cout<<"n_upmcs_calsz:"<<n_upmcs_calsz<<endl;
	cout<<"n_misel_calsz:"<<n_misel_calsz<<endl;
	cout<<"n_misinel_calsz:"<<n_misinel_calsz<<endl;
	cout<<"n_mismcs_calsz:"<<n_mismcs_calsz<<endl;
	cout<<"diff:"<<n_calsz-n_el_calsz-n_inel_calsz-n_mcs_calsz-n_upel_calsz-n_upinel_calsz-n_upmcs_calsz-n_misel_calsz-n_misinel_calsz-n_mismcs_calsz<<endl;
	cout<<"\n"<<endl;

	cout<<"n_bq:"<<n_bq<<endl;
	cout<<"n_el_bq:"<<n_el_bq<<endl;
	cout<<"n_inel_bq:"<<n_inel_bq<<endl;
	cout<<"n_mcs_bq:"<<n_mcs_bq<<endl;
	cout<<"n_upel_bq:"<<n_upel_bq<<endl;
	cout<<"n_upinel_bq:"<<n_upinel_bq<<endl;
	cout<<"n_upmcs_bq:"<<n_upmcs_bq<<endl;
	cout<<"n_misel_bq:"<<n_misel_bq<<endl;
	cout<<"n_misinel_bq:"<<n_misinel_bq<<endl;
	cout<<"n_mismcs_bq:"<<n_mismcs_bq<<endl;
	cout<<"diff:"<<n_bq-n_el_bq-n_inel_bq-n_mcs_bq-n_upel_bq-n_upinel_bq-n_upmcs_bq-n_misel_bq-n_misinel_bq-n_mismcs_bq<<endl;
	cout<<"\n"<<endl;

	cout<<"n_recoinel:"<<n_recoinel<<endl;
	cout<<"n_el_recoinel:"<<n_el_recoinel<<endl;
	cout<<"n_inel_recoinel:"<<n_inel_recoinel<<endl;
	cout<<"n_mcs_recoinel:"<<n_mcs_recoinel<<endl;
	cout<<"n_upel_recoinel:"<<n_upel_recoinel<<endl;
	cout<<"n_upinel_recoinel:"<<n_upinel_recoinel<<endl;
	cout<<"n_upmcs_recoinel:"<<n_upmcs_recoinel<<endl;
	cout<<"n_misel_recoinel:"<<n_misel_recoinel<<endl;
	cout<<"n_misinel_recoinel:"<<n_misinel_recoinel<<endl;
	cout<<"n_mismcs_recoinel:"<<n_mismcs_recoinel<<endl;
	cout<<"diff:"<<n_recoinel-n_el_recoinel-n_inel_recoinel-n_mcs_recoinel-n_upel_recoinel-n_upinel_recoinel-n_upmcs_recoinel-n_misel_recoinel-n_misinel_recoinel-n_mismcs_recoinel<<endl;



/*
	cout<<"n_true_end:"<<n_true_end<<endl;
	cout<<"n_processmap_error:"<<n_processmap_error<<endl;
	cout<<"  n_inel:"<<n_inel<<endl;
	cout<<"  n_el:"<<n_el<<endl;
	cout<<"  n_mcs:"<<n_mcs<<endl;
	cout<<"  n_misid:"<<n_misid<<endl;
	cout<<"    n_misid_inttpc:"<<n_misid_inttpc<<endl;
	cout<<"    n_misid_intoutside:"<<n_misid_intoutside<<endl;
	cout<<"  --> n_inel+n_el+n_mcs+n_misid="<<n_inel+n_el+n_mcs+n_misid<<endl;
	cout<<"\n"<<endl;	
	cout<<"n_pandora:"<<n_pandora<<endl;
	cout<<"  n_pandora_inel:"<<n_pandora_inel<<endl;
	cout<<"  n_pandora_el:"<<n_pandora_el<<endl;
	cout<<"  n_pandora_mcs:"<<n_pandora_mcs<<endl;
	cout<<"  n_pandora_misid:"<<n_pandora_misid<<endl;
	cout<<"  --> n_pandora_inel+n_pandora_el+n_pandora_mcs+n_pandora_misid="<<n_pandora_inel+n_pandora_el+n_pandora_mcs+n_pandora_misid<<endl;
	cout<<"\n"<<endl;	

	cout<<"n_calosz:"<<n_calosz<<endl;
	cout<<"  n_calosz_inel:"<<n_calosz_inel<<endl;
	cout<<"  n_calosz_el:"<<n_calosz_el<<endl;
	cout<<"  n_calosz_mcs:"<<n_calosz_mcs<<endl;
	cout<<"  n_calosz_misid:"<<n_calosz_misid<<endl;
	cout<<"  --> n_calosz_inel+n_calosz_el+n_calosz_mcs+n_calosz_misid="<<n_calosz_inel+n_calosz_el+n_calosz_mcs+n_calosz_misid<<endl;
	cout<<"\n"<<endl;

	cout<<"n_bq:"<<n_bq<<endl;
	cout<<"  n_bq_inel:"<<n_bq_inel<<endl;
	cout<<"  n_bq_el:"<<n_bq_el<<endl;
	cout<<"  n_bq_mcs:"<<n_bq_mcs<<endl;
	cout<<"  n_bq_misid:"<<n_bq_misid<<endl;
	cout<<"  -->n_bq_inel+n_bq_el+n_bq_mcs+n_bq_misid="<<n_bq_inel+n_bq_el+n_bq_mcs+n_bq_misid<<endl;
	cout<<"\n"<<endl;

	cout<<"n_recoinel:"<<n_recoinel<<endl;
	cout<<"  n_recoinel_inel:"<<n_recoinel_inel<<endl;
	cout<<"  n_recoinel_el:"<<n_recoinel_el<<endl;
	cout<<"  n_recoinel_mcs:"<<n_recoinel_mcs<<endl;
	cout<<"  n_recoinel_misid:"<<n_recoinel_misid<<endl;
	cout<<"  -->n_recoinel_inel+n_recoinel_el+n_recoinel_mcs+n_recoinel_misid="<<n_recoinel_inel+n_recoinel_el+n_recoinel_mcs+n_recoinel_misid<<endl;
	cout<<"\n"<<endl;
*/


}
