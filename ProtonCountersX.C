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

	//various counters --------------------------------------------------------------------------------------------------------------------//
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

		//isTestSample = true;
		//if (ientry%2 == 0) isTestSample = false; //Divide MC sample by 2 parts: test+ufold

		//true_sliceID = -1;
		//reco_sliceID = -1;
	
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
		//cout<<"IsTrueEndOutside:"<<IsTrueEndOutside<<endl;
		//if (IsPureInEL==1) cout<<"Summary(TrueEnd, Endoutside, Bm, Orig, EPDG):("<<1<<", "<<IsTrueEndOutside<<", "<<IsBeamMatch<<", "<<primary_truth_byE_origin<<", "<<primary_truth_byE_PDG<<")"<<endl;	
		//if (IsPureEL==1) cout<<"Summary(TrueEnd, Endoutside, Bm, Orig, EPDG):("<<2<<", "<<IsTrueEndOutside<<", "<<IsBeamMatch<<", "<<primary_truth_byE_origin<<", "<<primary_truth_byE_PDG<<")"<<endl;	
		//if (IsPureMCS==1) cout<<"Summary(TrueEnd, Endoutside, Bm, Orig, EPDG):("<<3<<", "<<IsTrueEndOutside<<", "<<IsBeamMatch<<", "<<primary_truth_byE_origin<<", "<<primary_truth_byE_PDG<<")"<<endl;	
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
		//vector< pair<double,int > > zreco_rawindex; //z, original_index
		//vector< pair<double,double > > zreco_rr; //z, rr
		//vector< pair<double,double > > zreco_de; //z, wid_reco
		//vector< pair<double,double > > zreco_dedx; //z, wid_reco
		//vector< pair<double,double > > zreco_dx; //z, wid_reco
		//vector< pair<double,double > > zreco_ke; //z, wid_reco
		//vector< pair<double,double > > zreco_xreco; //z, wid_reco
		//vector< pair<double,double > > zreco_yreco; //z, wid_reco
		//vector< pair<double,int > > zreco_widreco; //z, wid_reco

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

			//if (IsRecoStop) {
				//double theory1_dedx=(double)gr_predict_dedx_resrange->Eval(resrange_reco);
				//double theory2_dedx=dedx_predict(resrange_reco);

				//rr_dedx_recostop->Fill(resrange_reco,cali_dedx);

				//rr_ddedx1_recostop->Fill(resrange_reco,
			//}

			//if (IsPureMCS) {
				//rr_dedx_truestop->Fill(resrange_reco,cali_dedx);
			//}

			//if (IsPureInELProton) cout<<"["<<h<<"] z:"<<hitz_reco<<" dx:"<<pitch<<" wid:"<<wid_reco<<" dedx:"<<cali_dedx<<endl;
			//if (run==39279896&&event==1165) cout<<"["<<h<<"] z:"<<hitz_reco<<" dx:"<<pitch<<" wid:"<<wid_reco<<" dedx:"<<cali_dedx<<endl;
			//if (run==39279896&&event==1165) cout<<"["<<h<<"] z:"<<hitz_reco<<" dx:"<<pitch<<" wid:"<<wid_reco<<" pt:"<<pt_reco<<" dedx:"<<cali_dedx<<" rr:"<<resrange_reco<<endl;
			//if(IsPureMCS) cout<<"["<<h<<"] z:"<<hitz_reco<<" dx:"<<pitch<<" wid:"<<wid_reco<<" dedx:"<<cali_dedx<<endl;

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
		  } //loop over reco hits of a given track
		} //if calo size not empty

		//if (IsRecoStop) { 
			//KE_ff_recostop->Fill(ke_ff);
			//KE_calo_recostop->Fill(kereco_calo);
			//KE_rrange_recostop->Fill(kereco_range);
			//KE_rrange2_recostop->Fill(kereco_range2);
			//KE_range_recostop->Fill(ke_trklen);
			//KE_simide_recostop->Fill(ke_simide);

			//dKE_range_ff_recostop->Fill(ke_trklen-ke_ff);
			//dKE_calo_ff_recostop->Fill(kereco_calo-ke_ff);
			//dKE_rrange_ff_recostop->Fill(kereco_range-ke_ff);
			//dKE_rrange2_ff_recostop->Fill(kereco_range2-ke_ff);

			//KE_range_ff_recostop->Fill(kereco_range, ke_ff);
			//KE_range_calo_recostop->Fill(kereco_range, kereco_calo);
		//}
		//cout<<"\n"<<endl;

		//sort using wid
		//sort(zreco_rawindex.begin(),zreco_rawindex.end(),myComparison); //sorting based on the first column
		//sort(zreco_rr.begin(),zreco_rr.end(),myComparison); //sorting based on the first column
		//sort(zreco_de.begin(),zreco_de.end(),myComparison); //sorting based on the first column
		//sort(zreco_dedx.begin(),zreco_dedx.end(),myComparison); //sorting based on the first column
		//sort(zreco_dx.begin(),zreco_dx.end(),myComparison); //sorting based on the first column
		//sort(zreco_xreco.begin(),zreco_xreco.end(),myComparison); //sorting based on the first column
		//sort(zreco_yreco.begin(),zreco_yreco.end(),myComparison); //sorting based on the first column
		//sort(zreco_widreco.begin(),zreco_widreco.end(),myComparison); //sorting based on the first column
		//sort(zreco_ke.begin(),zreco_ke.end(),myComparison); //sorting based on the first column

		//countings -------------------------------------------------------//
		if (IsPandoraSlice) n_pan_tot++;
		if (IsPandoraSlice&&IsCaloSize) n_calsz_tot++;
		if (IsPandoraSlice&&IsCaloSize&&IsBQ) n_bq_tot++;
		if (IsPandoraSlice&&IsCaloSize&&IsBQ&&IsRecoInEL) n_recoinel_tot++;

		//if (IsBeamMatch) n_sg++;
		//if (!IsBeamMatch) n_bkg++;
		//if (!IsBeamMatch&&IsTrueEndOutside) n_up++;
		//if (!IsBeamMatch&&!IsTrueEndOutside) n_mis++;


		//[0]pure inel
		//if (IsBeamMatch&&IsPureInEL) { //pure inel
		if (kinel) { //pure inel
			n_inel++;
			//zend_true_inel_NoCut->Fill(true_endz);
			//zend_reco_inel_NoCut->Fill(reco_endz);
			//dzend_inel_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_inel_pan++;
				//zend_true_inel_PanS->Fill(true_endz);
				//zend_reco_inel_PanS->Fill(reco_endz);
				//dzend_inel_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_inel_calsz++;
					//zend_true_inel_CaloSz->Fill(true_endz);
					//zend_reco_inel_CaloSz->Fill(reco_endz);
					//dzend_inel_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_inel_bq++;
						//zend_true_inel_BQ->Fill(true_endz);
						//zend_reco_inel_BQ->Fill(reco_endz);
						//dzend_inel_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_inel_recoinel++;
							//zend_true_inel_RecoInel->Fill(true_endz);
							//zend_reco_inel_RecoInel->Fill(reco_endz);
							//dzend_inel_RecoInel->Fill(reco_endz-true_endz);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //pure inel		

		//[1]pure el
		//if (IsBeamMatch&&IsPureEL) { //pure el
		if (kel) { //pure el
			n_el++;
			//zend_true_el_NoCut->Fill(true_endz);
			//zend_reco_el_NoCut->Fill(reco_endz);
			//dzend_el_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_el_pan++;
				//zend_true_el_PanS->Fill(true_endz);
				//zend_reco_el_PanS->Fill(reco_endz);
				//dzend_el_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_el_calsz++;
					//zend_true_el_CaloSz->Fill(true_endz);
					//zend_reco_el_CaloSz->Fill(reco_endz);
					//dzend_el_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_el_bq++;
						//zend_true_el_BQ->Fill(true_endz);
						//zend_reco_el_BQ->Fill(reco_endz);
						//dzend_el_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_el_recoinel++;
							//zend_true_el_RecoInel->Fill(true_endz);
							//zend_reco_el_RecoInel->Fill(reco_endz);
							//dzend_el_RecoInel->Fill(reco_endz-true_endz);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //pure el		

		//[2]pure mcs
		//if (IsBeamMatch&&IsPureMCS) { //pure mcs
		if (kmcs) { //pure mcs
			n_mcs++;
			//zend_true_mcs_NoCut->Fill(true_endz);
			//zend_reco_mcs_NoCut->Fill(reco_endz);
			//dzend_mcs_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_mcs_pan++;
				//zend_true_mcs_PanS->Fill(true_endz);
				//zend_reco_mcs_PanS->Fill(reco_endz);
				//dzend_mcs_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_mcs_calsz++;
					//zend_true_mcs_CaloSz->Fill(true_endz);
					//zend_reco_mcs_CaloSz->Fill(reco_endz);
					//dzend_mcs_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_mcs_bq++;
						//zend_true_mcs_BQ->Fill(true_endz);
						//zend_reco_mcs_BQ->Fill(reco_endz);
						//dzend_mcs_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_mcs_recoinel++;
							//zend_true_mcs_RecoInel->Fill(true_endz);
							//zend_reco_mcs_RecoInel->Fill(reco_endz);
							//dzend_mcs_RecoInel->Fill(reco_endz-true_endz);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //pure mcs

		//[3]MID:Cosmic
		if (kMIDcosmic) { //mid:cosmic
			n_midcosmic++;
			//zend_true_midcosmic_NoCut->Fill(true_endz);
			//zend_reco_midcosmic_NoCut->Fill(reco_endz);
			//dzend_midcosmic_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_midcosmic_pan++;
				//zend_true_midcosmic_PanS->Fill(true_endz);
				//zend_reco_midcosmic_PanS->Fill(reco_endz);
				//dzend_midcosmic_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_midcosmic_calsz++;
					//zend_true_midcosmic_CaloSz->Fill(true_endz);
					//zend_reco_midcosmic_CaloSz->Fill(reco_endz);
					//dzend_midcosmic_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_midcosmic_bq++;
						//zend_true_midcosmic_BQ->Fill(true_endz);
						//zend_reco_midcosmic_BQ->Fill(reco_endz);
						//dzend_midcosmic_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_midcosmic_recoinel++;
							//zend_true_midcosmic_RecoInel->Fill(true_endz);
							//zend_reco_midcosmic_RecoInel->Fill(reco_endz);
							//dzend_midcosmic_RecoInel->Fill(reco_endz-true_endz);
						} //reco inel
					} //bq
				} //calosz
			} //pandora

			
		} //mid:cosmic

		//[4]MID:midpi
		if (kMIDpi) { //mid:pi
			n_midpi++;
			//zend_true_midpi_NoCut->Fill(true_endz);
			//zend_reco_midpi_NoCut->Fill(reco_endz);
			//dzend_midpi_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_midpi_pan++;
				//zend_true_midpi_PanS->Fill(true_endz);
				//zend_reco_midpi_PanS->Fill(reco_endz);
				//dzend_midpi_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_midpi_calsz++;
					//zend_true_midpi_CaloSz->Fill(true_endz);
					//zend_reco_midpi_CaloSz->Fill(reco_endz);
					//dzend_midpi_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_midpi_bq++;
						//zend_true_midpi_BQ->Fill(true_endz);
						//zend_reco_midpi_BQ->Fill(reco_endz);
						//dzend_midpi_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_midpi_recoinel++;
							//zend_true_midpi_RecoInel->Fill(true_endz);
							//zend_reco_midpi_RecoInel->Fill(reco_endz);
							//dzend_midpi_RecoInel->Fill(reco_endz-true_endz);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //mid:pi

		//[5]MID:midp
		if (kMIDp) { //mid:p
			n_midp++;
			//zend_true_midp_NoCut->Fill(true_endz);
			//zend_reco_midp_NoCut->Fill(reco_endz);
			//dzend_midp_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_midp_pan++;
				//zend_true_midp_PanS->Fill(true_endz);
				//zend_reco_midp_PanS->Fill(reco_endz);
				//dzend_midp_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_midp_calsz++;
					//zend_true_midp_CaloSz->Fill(true_endz);
					//zend_reco_midp_CaloSz->Fill(reco_endz);
					//dzend_midp_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_midp_bq++;
						//zend_true_midp_BQ->Fill(true_endz);
						//zend_reco_midp_BQ->Fill(reco_endz);
						//dzend_midp_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_midp_recoinel++;
							//zend_true_midp_RecoInel->Fill(true_endz);
							//zend_reco_midp_RecoInel->Fill(reco_endz);
							//dzend_midp_RecoInel->Fill(reco_endz-true_endz);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //mid:p


		//[6]MID:midmu
		if (kMIDmu) { //mid:mu
			n_midmu++;
			//zend_true_midmu_NoCut->Fill(true_endz);
			//zend_reco_midmu_NoCut->Fill(reco_endz);
			//dzend_midmu_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_midmu_pan++;
				//zend_true_midmu_PanS->Fill(true_endz);
				//zend_reco_midmu_PanS->Fill(reco_endz);
				//dzend_midmu_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_midmu_calsz++;
					//zend_true_midmu_CaloSz->Fill(true_endz);
					//zend_reco_midmu_CaloSz->Fill(reco_endz);
					//dzend_midmu_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_midmu_bq++;
						//zend_true_midmu_BQ->Fill(true_endz);
						//zend_reco_midmu_BQ->Fill(reco_endz);
						//dzend_midmu_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_midmu_recoinel++;
							//zend_true_midmu_RecoInel->Fill(true_endz);
							//zend_reco_midmu_RecoInel->Fill(reco_endz);
							//dzend_midmu_RecoInel->Fill(reco_endz-true_endz);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //mid:mu


		//[7]MID:mideg
		if (kMIDeg) { //mid:eg
			n_mideg++;
			//zend_true_mideg_NoCut->Fill(true_endz);
			//zend_reco_mideg_NoCut->Fill(reco_endz);
			//dzend_mideg_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_mideg_pan++;
				//zend_true_mideg_PanS->Fill(true_endz);
				//zend_reco_mideg_PanS->Fill(reco_endz);
				//dzend_mideg_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_mideg_calsz++;
					//zend_true_mideg_CaloSz->Fill(true_endz);
					//zend_reco_mideg_CaloSz->Fill(reco_endz);
					//dzend_mideg_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_mideg_bq++;
						//zend_true_mideg_BQ->Fill(true_endz);
						//zend_reco_mideg_BQ->Fill(reco_endz);
						//dzend_mideg_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_mideg_recoinel++;
							//zend_true_mideg_RecoInel->Fill(true_endz);
							//zend_reco_mideg_RecoInel->Fill(reco_endz);
							//dzend_mideg_RecoInel->Fill(reco_endz-true_endz);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //mid:eg


		//[8]MID:midother
		if (kMIDother) { //mid:other
			n_midother++;
			//zend_true_midother_NoCut->Fill(true_endz);
			//zend_reco_midother_NoCut->Fill(reco_endz);
			//dzend_midother_NoCut->Fill(reco_endz-true_endz);
			if (IsPandoraSlice) { //pandora
				n_midother_pan++;
				//zend_true_midother_PanS->Fill(true_endz);
				//zend_reco_midother_PanS->Fill(reco_endz);
				//dzend_midother_PanS->Fill(reco_endz-true_endz);
				if (IsCaloSize) { //calosz
					n_midother_calsz++;
					//zend_true_midother_CaloSz->Fill(true_endz);
					//zend_reco_midother_CaloSz->Fill(reco_endz);
					//dzend_midother_CaloSz->Fill(reco_endz-true_endz);
					if (IsBQ) { //bq
						n_midother_bq++;
						//zend_true_midother_BQ->Fill(true_endz);
						//zend_reco_midother_BQ->Fill(reco_endz);
						//dzend_midother_BQ->Fill(reco_endz-true_endz);
						if (IsRecoInEL) { //reco inel
							n_midother_recoinel++;
							//zend_true_midother_RecoInel->Fill(true_endz);
							//zend_reco_midother_RecoInel->Fill(reco_endz);
							//dzend_midother_RecoInel->Fill(reco_endz-true_endz);
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
	cout<<"n_el:"<<n_el<<endl;
	cout<<"n_inel:"<<n_inel<<endl;
	cout<<"n_mcs:"<<n_mcs<<endl;
	cout<<"n_midcosmic:"<<n_midcosmic<<endl;
	cout<<"n_midpi:"<<n_midpi<<endl;
	cout<<"n_midp:"<<n_midp<<endl;
	cout<<"n_midmu:"<<n_midmu<<endl;
	cout<<"n_mideg:"<<n_mideg<<endl;
	cout<<"n_midother"<<n_midother<<endl;
	cout<<"n_diff:"<<n_tot-(n_el+n_inel+n_mcs+n_midcosmic+n_midpi+n_midp+n_midmu+n_mideg+n_midother)<<endl;

	cout<<"\nn_pan_tot:"<<n_pan_tot<<endl;
	cout<<"n_el_pan:"<<n_el_pan<<endl;
	cout<<"n_inel_pan:"<<n_inel_pan<<endl;
	cout<<"n_mcs_pan:"<<n_mcs_pan<<endl;
	cout<<"n_midcosmic_pan:"<<n_midcosmic_pan<<endl;
	cout<<"n_midpi_pan:"<<n_midpi_pan<<endl;
	cout<<"n_midp_pan:"<<n_midp_pan<<endl;
	cout<<"n_midmu_pan:"<<n_midmu_pan<<endl;
	cout<<"n_mideg_pan:"<<n_mideg_pan<<endl;
	cout<<"n_midother_pan"<<n_midother_pan<<endl;
	cout<<"n_diff_pan:"<<n_pan_tot-(n_el_pan+n_inel_pan+n_mcs_pan+n_midcosmic_pan+n_midpi_pan+n_midp_pan+n_midmu_pan+n_mideg_pan+n_midother_pan)<<endl;

	cout<<"\nn_calsz_tot:"<<n_calsz_tot<<endl;
	cout<<"n_el_calsz:"<<n_el_calsz<<endl;
	cout<<"n_inel_calsz:"<<n_inel_calsz<<endl;
	cout<<"n_mcs_calsz:"<<n_mcs_calsz<<endl;
	cout<<"n_midcosmic_calsz:"<<n_midcosmic_calsz<<endl;
	cout<<"n_midpi_calsz:"<<n_midpi_calsz<<endl;
	cout<<"n_midp_calsz:"<<n_midp_calsz<<endl;
	cout<<"n_midmu_calsz:"<<n_midmu_calsz<<endl;
	cout<<"n_mideg_calsz:"<<n_mideg_calsz<<endl;
	cout<<"n_midother_calsz"<<n_midother_calsz<<endl;
	cout<<"n_diff_calsz:"<<n_calsz_tot-(n_el_calsz+n_inel_calsz+n_mcs_calsz+n_midcosmic_calsz+n_midpi_calsz+n_midp_calsz+n_midmu_calsz+n_mideg_calsz+n_midother_calsz)<<endl;

	cout<<"\nn_bq_tot:"<<n_bq_tot<<endl;
	cout<<"n_el_bq:"<<n_el_bq<<endl;
	cout<<"n_inel_bq:"<<n_inel_bq<<endl;
	cout<<"n_mcs_bq:"<<n_mcs_bq<<endl;
	cout<<"n_midcosmic_bq:"<<n_midcosmic_bq<<endl;
	cout<<"n_midpi_bq:"<<n_midpi_bq<<endl;
	cout<<"n_midp_bq:"<<n_midp_bq<<endl;
	cout<<"n_midmu_bq:"<<n_midmu_bq<<endl;
	cout<<"n_mideg_bq:"<<n_mideg_bq<<endl;
	cout<<"n_midother_bq"<<n_midother_bq<<endl;
	cout<<"n_diff_bq:"<<n_bq_tot-(n_el_bq+n_inel_bq+n_mcs_bq+n_midcosmic_bq+n_midpi_bq+n_midp_bq+n_midmu_bq+n_mideg_bq+n_midother_bq)<<endl;

	cout<<"\nn_recoinel_tot:"<<n_recoinel_tot<<endl;
	cout<<"n_el_recoinel:"<<n_el_recoinel<<endl;
	cout<<"n_inel_recoinel:"<<n_inel_recoinel<<endl;
	cout<<"n_mcs_recoinel:"<<n_mcs_recoinel<<endl;
	cout<<"n_midcosmic_recoinel:"<<n_midcosmic_recoinel<<endl;
	cout<<"n_midpi_recoinel:"<<n_midpi_recoinel<<endl;
	cout<<"n_midp_recoinel:"<<n_midp_recoinel<<endl;
	cout<<"n_midmu_recoinel:"<<n_midmu_recoinel<<endl;
	cout<<"n_mideg_recoinel:"<<n_mideg_recoinel<<endl;
	cout<<"n_midother_recoinel"<<n_midother_recoinel<<endl;
	cout<<"n_diff_recoinel:"<<n_recoinel_tot-(n_el_recoinel+n_inel_recoinel+n_mcs_recoinel+n_midcosmic_recoinel+n_midpi_recoinel+n_midp_recoinel+n_midmu_recoinel+n_mideg_recoinel+n_midother_recoinel)<<endl;



        vector<int> cnt_x;
        vector<int> cnt_y_el;
        vector<int> cnt_y_inel;
        vector<int> cnt_y_mcs;
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

        cnt_y_mcs.push_back(n_mcs);
        cnt_y_mcs.push_back(n_mcs_pan);
        cnt_y_mcs.push_back(n_mcs_calsz);
        cnt_y_mcs.push_back(n_mcs_bq);
        cnt_y_mcs.push_back(n_mcs_recoinel);

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
        TGraph *Counters_mcs = new TGraph(cnt_x.size(), &cnt_x.at(0), &cnt_y_mcs.at(0));
        TGraph *Counters_midcosmic = new TGraph(cnt_x.size(), &cnt_x.at(0), &cnt_y_midcosmic.at(0));
        TGraph *Counters_midpi = new TGraph(cnt_x.size(), &cnt_x.at(0), &cnt_y_midpi.at(0));
        TGraph *Counters_midp = new TGraph(cnt_x.size(), &cnt_x.at(0), &cnt_y_midp.at(0));
        TGraph *Counters_midmu = new TGraph(cnt_x.size(), &cnt_x.at(0), &cnt_y_midmu.at(0));
        TGraph *Counters_mideg = new TGraph(cnt_x.size(), &cnt_x.at(0), &cnt_y_mideg.at(0));
        TGraph *Counters_midother = new TGraph(cnt_x.size(), &cnt_x.at(0), &cnt_y_midother.at(0));


        Counters_el->SetName("Counters_el");
        Counters_inel->SetName("Counters_inel");
        Counters_mcs->SetName("Counters_mcs");
        Counters_midcosmic->SetName("Counters_midcosmic");
        Counters_midpi->SetName("Counters_midpi");
        Counters_midp->SetName("Counters_midp");
        Counters_midmu->SetName("Counters_midmu");
        Counters_mideg->SetName("Counters_mideg");
        Counters_midother->SetName("Counters_midother");


        TFile *fout = new TFile(Form("Counters_MC.root"),"RECREATE");
		Counters_el->Write();
		Counters_inel->Write();
		Counters_mcs->Write();
		Counters_midcosmic->Write();
		Counters_midpi->Write();
		Counters_midp->Write();
		Counters_midmu->Write();
		Counters_mideg->Write();
		Counters_midother->Write();
        fout->Close();







}
