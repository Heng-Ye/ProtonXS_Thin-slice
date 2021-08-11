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
#include "./headers/util.h"
#include "./headers/ThinSlice.h"

using namespace std;
using namespace ROOT::Math;

void ProtonThinSlice::Loop() {
	if (fChain == 0) return;

	//various counters ---------------------------------------------------------------------------------------------------------------------------------------------------//
	int n_processmap_error=0; //sansity check: processmap
	int n_true_end=0; //sansity check: true-end label

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
		bool IsIntersection=false; //if any track intersect with our reco track		
		if (primary_truth_Isbeammatched==1) IsBeamMatch=true;
		if (isprimarytrack==1&&isprimaryshower==0) IsPandoraSlice=true; 
		if (!primtrk_hitz->empty()) IsCaloSize=true;
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
		cout<<"primary_truth_byE_origin="<<primary_truth_byE_origin<<""<<endl;
		cout<<"primary_truth_byE_PDG="<<primary_truth_byE_PDG<<""<<endl;
		cout<<"primary_truth_Pdg:"<<primary_truth_Pdg<<endl;
		cout<<"beamtrackPdg:"<<beamtrackPdg<<endl;

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
					//if (intz<0) { //if interaction outside tpc
					//if(strcmp(interactionProcesslist->at(iiii).c_str(),"Transportation")!=0) {
				} //loop over all true interaction hits in this track 
			} //size of interactionProcesslist >=0
*/
		} //hIoni

		//if (IsPureInEL==0&&IsPureEL==0) {
		//	//IsPureMCS=1;
		//}

		////if (strcmp(primary_truth_EndProcess->c_str(),"hIoni")==0) {
			////IsPureEL=true;
		////}
		////if (strcmp(primary_truth_EndProcess->c_str(),"CoulombScat")==0) {
			////IsPureMCS=true;
		////}
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
		cout<<"trueEnd z/y/x:"<<true_endz<<"/"<<true_endy<<"/"<<true_endx<<endl;
		cout<<"trueSt z/y/x:"<<true_stz<<"/"<<true_sty<<"/"<<true_stx<<endl;
		//cout<<"InEL EL MCS:"<<IsPureInEL<<" "<<IsPureEL<<" "<<IsPureMCS<<endl;
		cout<<"InEL EL:"<<IsPureInEL<<" "<<IsPureEL<<" "<<endl;
		cout<<"IsTrueEndOutside:"<<IsTrueEndOutside<<endl;
		if (IsPureInEL==1) cout<<"Summary(TrueEnd, Endoutside, Bm, Orig, EPDG):("<<1<<", "<<IsTrueEndOutside<<", "<<IsBeamMatch<<", "<<primary_truth_byE_origin<<", "<<primary_truth_byE_PDG<<")"<<endl;	
		if (IsPureEL==1) cout<<"Summary(TrueEnd, Endoutside, Bm, Orig, EPDG):("<<2<<", "<<IsTrueEndOutside<<", "<<IsBeamMatch<<", "<<primary_truth_byE_origin<<", "<<primary_truth_byE_PDG<<")"<<endl;	
		//if (IsPureMCS==1) cout<<"Summary(TrueEnd, Endoutside, Bm, Orig, EPDG):("<<3<<", "<<IsTrueEndOutside<<", "<<IsBeamMatch<<", "<<primary_truth_byE_origin<<", "<<primary_truth_byE_PDG<<")"<<endl;	

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
  		}
		//cout<<"ck1"<<endl;
  		for (int iz=key_st+1; iz<(int)beamtrk_z->size(); iz++){
    			if (iz == key_st+1) range_true = 0;
    			range_true += sqrt( pow(beamtrk_x->at(iz)-beamtrk_x->at(iz-1), 2)+
					    pow(beamtrk_y->at(iz)-beamtrk_y->at(iz-1), 2)+	
					    pow(beamtrk_z->at(iz)-beamtrk_z->at(iz-1), 2) );						    	
    			true_trklen_accum[iz] = range_true;
  		}
	
		cout<<"range_true:"<<range_true<<endl;
		cout<<"key_st:"<<key_st<<endl;
		//for (size_t j=0; j<beamtrk_z->size(); ++j) { //MCParticle loop
			//cout<<"beamtrk_z["<<j<<"]:"<<beamtrk_z->at(j)<<" beamtrk_Eng["<<j<<"]:"<<beamtrk_Eng->at(j)<<" true_trklen_accum["<<j<<"]:"<<true_trklen_accum[j]<<endl;
		//} //MCParticle loop
		//Get reco info ----------------------------------------------------------------------------------//

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
		cout<<"kMIDcosmic:"<<kMIDcosmic<<endl;
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

			//reco_startX_sce->Fill(reco_stx);
			//reco_startY_sce->Fill(reco_sty);
			//reco_startZ_sce->Fill(reco_stz);
			
			Fill1DHist(reco_startX_sce, reco_stx);
			Fill1DHist(reco_startY_sce, reco_sty);
			Fill1DHist(reco_startZ_sce, reco_stz);

			double beam_dx=(reco_stx-mean_StartX)/sigma_StartX;
			double beam_dy=(reco_sty-mean_StartY)/sigma_StartY;
			double beam_dz=(reco_stz-mean_StartZ)/sigma_StartZ;
			double beam_dxy=sqrt(pow(beam_dx,2)+pow(beam_dy,2));	
	
			//hdeltaX->Fill(beam_dx);
			//hdeltaY->Fill(beam_dy);
			//hdeltaZ->Fill(beam_dz);
			//hdeltaXY->Fill(beam_dxy);

			Fill1DHist(hdeltaX, beam_dx);
			Fill1DHist(hdeltaY, beam_dy);
			Fill1DHist(hdeltaZ, beam_dz);
			Fill1DHist(hdeltaXY, beam_dxy);

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
		double cosine_beam_spec_primtrk=-99; 
		cosine_beam_spec_primtrk=beamDirx_spec->at(0)*primaryStartDirection[0]+beamDiry_spec->at(0)*primaryStartDirection[1]+beamDirz_spec->at(0)*primaryStartDirection[2]; //cosine between beam_spec and primary trk direction
		if (cosine_beam_spec_primtrk<0) { cosine_beam_spec_primtrk=-1.*cosine_beam_spec_primtrk; }
		//if (cosine_beam_spec_primtrk>cosine_beam_primtrk_min) { IsCosine=true; }
		if (cosine_beam_spec_primtrk>costh_min&&cosine_beam_spec_primtrk<costh_max) { IsCosine=true; }
		//reco_cosineTheta->Fill(cosine_beam_spec_primtrk);
		Fill1DHist(reco_cosineTheta, cosine_beam_spec_primtrk);
		//if (kinel) reco_cosineTheta_inel->Fill(cosine_beam_spec_primtrk);
		//if (kel) reco_cosineTheta_el->Fill(cosine_beam_spec_primtrk);
		//if (kMIDcosmic) reco_cosineTheta_midcosmic->Fill(cosine_beam_spec_primtrk);
		//if (kMIDpi) reco_cosineTheta_midpi->Fill(cosine_beam_spec_primtrk);
		//if (kMIDp) reco_cosineTheta_midp->Fill(cosine_beam_spec_primtrk);
		//if (kMIDmu) reco_cosineTheta_midmu->Fill(cosine_beam_spec_primtrk);
		//if (kMIDeg) reco_cosineTheta_mideg->Fill(cosine_beam_spec_primtrk);
		//if (kMIDother) reco_cosineTheta_midother->Fill(cosine_beam_spec_primtrk);

		if (kinel) Fill1DHist(reco_cosineTheta_inel, cosine_beam_spec_primtrk);
		if (kel) Fill1DHist(reco_cosineTheta_el, cosine_beam_spec_primtrk);
		if (kMIDcosmic) Fill1DHist(reco_cosineTheta_midcosmic, cosine_beam_spec_primtrk);
		if (kMIDpi) Fill1DHist(reco_cosineTheta_midpi, cosine_beam_spec_primtrk);
		if (kMIDp) Fill1DHist(reco_cosineTheta_midp, cosine_beam_spec_primtrk);
		if (kMIDmu) Fill1DHist(reco_cosineTheta_midmu, cosine_beam_spec_primtrk);
		if (kMIDeg) Fill1DHist(reco_cosineTheta_mideg, cosine_beam_spec_primtrk);
		if (kMIDother) Fill1DHist(reco_cosineTheta_midother, cosine_beam_spec_primtrk);

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

			if (IsPureInEL) rangereco_dedxreco_TrueInEL->Fill(range_reco, cali_dedx);
			if (IsPureEL) { 
						rangereco_dedxreco_TrueEL->Fill(range_reco, cali_dedx);
						rr_dedx_truestop->Fill(resrange_reco, cali_dedx);
			}

		  } //loop over reco hits of a given track
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
		cout<<"range_reco:"<<range_reco<<endl;

		//Reco stopping/Inel p cut -----------------------------------------------------------------------------------------------//
		bool IsRecoStop=false;
		bool IsRecoInEL=false;
		double mom_beam_spec=-99; mom_beam_spec=beamMomentum_spec->at(0);
		//double range_reco=-99; if (!primtrk_range->empty()) range_reco=primtrk_range->at(0); //reco primary trklen
		double csda_val_spec=csda_range_vs_mom_sm->Eval(mom_beam_spec);

		if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) IsRecoStop=true;
		if ((range_reco/csda_val_spec)<min_norm_trklen_csda) IsRecoInEL=true;

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

		//some ke calc. -------------------------------------------------------------------------------------//
		//double kereco_calo=0;
		//double kereco_range=0;
		//double kereco_range2=0;
		if (IsCaloSize) { //if calo size not empty
  			//for (int iz=1; iz<(int)zreco_rawindex.size(); iz++) { //calo hit loop
				//double cali_dedx=zreco_dedx[iz].second;
				//double pitch=zreco_dx[iz].second;
				//double len=zreco_lenreco[iz].second;
				//double rr=range_reco-len;

				//kereco_calo+=cali_dedx*pitch;
				//kereco_range+=pitch*dedx_predict(rr);
				//kereco_range2+=pitch*(double)gr_predict_dedx_resrange->Eval(rr);

				//if (IsRecoStop) {
					//rr_dedx_recostop->Fill(rr, cali_dedx);
				//}

				//if (IsPureInEL) rangereco_dedxreco_TrueInEL->Fill(len, cali_dedx);
				//if (IsPureEL) { 
						//rangereco_dedxreco_TrueEL->Fill(len, cali_dedx);
						//rr_dedx_truestop->Fill(rr,cali_dedx);
				//}
				//if (IsPureMCS) rangereco_dedxreco_TrueMCS->Fill(len, cali_dedx);
			//} //calo hit loop

			if (IsRecoStop) { //reco_stop 
				Fill1DHist(KE_ff_recostop, ke_ff);
				Fill1DHist(KE_calo_recostop, kereco_calo);
				Fill1DHist(KE_rrange_recostop, kereco_range);
				Fill1DHist(KE_rrange2_recostop, kereco_range2);
				Fill1DHist(KE_range_recostop, ke_trklen);
				Fill1DHist(KE_simide_recostop, ke_simide);

				Fill1DHist(dKE_range_ff_recostop, ke_trklen-ke_ff);
				Fill1DHist(dKE_calo_ff_recostop, kereco_calo-ke_ff);
				Fill1DHist(dKE_rrange_ff_recostop, kereco_range-ke_ff);
				Fill1DHist(dKE_rrange2_ff_recostop, kereco_range2-ke_ff);

				KE_range_ff_recostop->Fill(kereco_range, ke_ff);
				KE_range_calo_recostop->Fill(kereco_range, kereco_calo);

		  		for (size_t h=0; h<primtrk_dedx->size(); ++h) { //loop over reco hits of a track
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

					rr_dedx_recostop->Fill(resrange_reco, cali_dedx);
				} //loop over reco hits of a track
			} //reco_stop
		} //if calo size not empty

		//countings -------------------------------------------------------//
		if (IsPandoraSlice) n_pan_tot++;
		if (IsPandoraSlice&&IsCaloSize) n_calsz_tot++;
		if (IsPandoraSlice&&IsCaloSize&&IsBQ) n_bq_tot++;
		if (IsPandoraSlice&&IsCaloSize&&IsBQ&&IsRecoInEL) n_recoinel_tot++;

		//[0]pure inel
		if (kinel) { //pure inel
			n_inel++;
			Fill1DHist(trklen_true_inel_NoCut, range_true);
			Fill1DHist(trklen_reco_inel_NoCut, range_reco);
			Fill1DHist(dtrklen_inel_NoCut, range_reco-range_true);
			if (IsPandoraSlice) { //pandora
				n_inel_pan++;
				Fill1DHist(trklen_true_inel_PanS, range_true);
				Fill1DHist(trklen_reco_inel_PanS, range_reco);
				Fill1DHist(dtrklen_inel_PanS,range_reco-range_true);
				if (IsCaloSize) { //calosz
					n_inel_calsz++;
					Fill1DHist(trklen_true_inel_CaloSz, range_true);
					Fill1DHist(trklen_reco_inel_CaloSz, range_reco);
					Fill1DHist(dtrklen_inel_CaloSz, range_reco-range_true);
					if (IsBQ) { //bq
						n_inel_bq++;
						Fill1DHist(trklen_true_inel_BQ, range_true);
						Fill1DHist(trklen_reco_inel_BQ, range_reco);
						Fill1DHist(dtrklen_inel_BQ, range_reco-range_true);
						if (IsRecoInEL) { //reco inel
							n_inel_recoinel++;
							Fill1DHist(trklen_true_inel_RecoInel, range_true);
							Fill1DHist(trklen_reco_inel_RecoInel, range_reco);
							Fill1DHist(dtrklen_inel_RecoInel, range_reco-range_true);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //pure inel		

		//[1]pure el
		if (kel) { //pure el
			n_el++;
			Fill1DHist(trklen_true_el_NoCut, range_true);
			Fill1DHist(trklen_reco_el_NoCut, range_reco);
			Fill1DHist(dtrklen_el_NoCut, range_reco-range_true);
			if (IsPandoraSlice) { //pandora
				n_el_pan++;
				Fill1DHist(trklen_true_el_PanS, range_true);
				Fill1DHist(trklen_reco_el_PanS, range_reco);
				Fill1DHist(dtrklen_el_PanS, range_reco-range_true);
				if (IsCaloSize) { //calosz
					n_el_calsz++;
					Fill1DHist(trklen_true_el_CaloSz, range_true);
					Fill1DHist(trklen_reco_el_CaloSz, range_reco);
					Fill1DHist(dtrklen_el_CaloSz, range_reco-range_true);
					if (IsBQ) { //bq
						n_el_bq++;
						Fill1DHist(trklen_true_el_BQ, range_true);
						Fill1DHist(trklen_reco_el_BQ, range_reco);
						Fill1DHist(dtrklen_el_BQ, range_reco-range_true);
						if (IsRecoInEL) { //reco inel
							n_el_recoinel++;
							Fill1DHist(trklen_true_el_RecoInel, range_true);
							Fill1DHist(trklen_reco_el_RecoInel, range_reco);
							Fill1DHist(dtrklen_el_RecoInel, range_reco-range_true);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //pure el		

		//[3]MID:Cosmic
		if (kMIDcosmic) { //mid:cosmic
			n_midcosmic++;
			Fill1DHist(trklen_true_midcosmic_NoCut, range_true);
			Fill1DHist(trklen_reco_midcosmic_NoCut, range_reco);
			Fill1DHist(dtrklen_midcosmic_NoCut, range_reco-range_true);
			if (IsPandoraSlice) { //pandora
				n_midcosmic_pan++;
				Fill1DHist(trklen_true_midcosmic_PanS, range_true);
				Fill1DHist(trklen_reco_midcosmic_PanS, range_reco);
				Fill1DHist(dtrklen_midcosmic_PanS, range_reco-range_true);
				if (IsCaloSize) { //calosz
					n_midcosmic_calsz++;
					Fill1DHist(trklen_true_midcosmic_CaloSz, range_true);
					Fill1DHist(trklen_reco_midcosmic_CaloSz, range_reco);
					Fill1DHist(dtrklen_midcosmic_CaloSz, range_reco-range_true);
					if (IsBQ) { //bq
						n_midcosmic_bq++;
						Fill1DHist(trklen_true_midcosmic_BQ, range_true);
						Fill1DHist(trklen_reco_midcosmic_BQ, range_reco);
						Fill1DHist(dtrklen_midcosmic_BQ, range_reco-range_true);
						if (IsRecoInEL) { //reco inel
							n_midcosmic_recoinel++;
							Fill1DHist(trklen_true_midcosmic_RecoInel, range_true);
							Fill1DHist(trklen_reco_midcosmic_RecoInel, range_reco);
							Fill1DHist(dtrklen_midcosmic_RecoInel, range_reco-range_true);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //mid:cosmic

		//[4]MID:midpi
		if (kMIDpi) { //mid:pi
			n_midpi++;
			Fill1DHist(trklen_true_midpi_NoCut, range_true);
			Fill1DHist(trklen_reco_midpi_NoCut, range_reco);
			Fill1DHist(dtrklen_midpi_NoCut, range_reco-range_true);
			if (IsPandoraSlice) { //pandora
				n_midpi_pan++;
				Fill1DHist(trklen_true_midpi_PanS, range_true);
				Fill1DHist(trklen_reco_midpi_PanS, range_reco);
				Fill1DHist(dtrklen_midpi_PanS, range_reco-range_true);
				if (IsCaloSize) { //calosz
					n_midpi_calsz++;
					Fill1DHist(trklen_true_midpi_CaloSz, range_true);
					Fill1DHist(trklen_reco_midpi_CaloSz, range_reco);
					Fill1DHist(dtrklen_midpi_CaloSz, range_reco-range_true);
					if (IsBQ) { //bq
						n_midpi_bq++;
						Fill1DHist(trklen_true_midpi_BQ, range_true);
						Fill1DHist(trklen_reco_midpi_BQ, range_reco);
						Fill1DHist(dtrklen_midpi_BQ, range_reco-range_true);
						if (IsRecoInEL) { //reco inel
							n_midpi_recoinel++;
							Fill1DHist(trklen_true_midpi_RecoInel, range_true);
							Fill1DHist(trklen_reco_midpi_RecoInel, range_reco);
							Fill1DHist(dtrklen_midpi_RecoInel, range_reco-range_true);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //mid:pi

		//[5]MID:midp
		if (kMIDp) { //mid:p
			n_midp++;
			Fill1DHist(trklen_true_midp_NoCut, range_true);
			Fill1DHist(trklen_reco_midp_NoCut, range_reco);
			Fill1DHist(dtrklen_midp_NoCut, range_reco-range_true);
			if (IsPandoraSlice) { //pandora
				n_midp_pan++;
				Fill1DHist(trklen_true_midp_PanS, range_true);
				Fill1DHist(trklen_reco_midp_PanS, range_reco);
				Fill1DHist(dtrklen_midp_PanS, range_reco-range_true);
				if (IsCaloSize) { //calosz
					n_midp_calsz++;
					Fill1DHist(trklen_true_midp_CaloSz, range_true);
					Fill1DHist(trklen_reco_midp_CaloSz, range_reco);
					Fill1DHist(dtrklen_midp_CaloSz, range_reco-range_true);
					if (IsBQ) { //bq
						n_midp_bq++;
						Fill1DHist(trklen_true_midp_BQ, range_true);
						Fill1DHist(trklen_reco_midp_BQ, range_reco);
						Fill1DHist(dtrklen_midp_BQ, range_reco-range_true);
						if (IsRecoInEL) { //reco inel
							n_midp_recoinel++;
							Fill1DHist(trklen_true_midp_RecoInel, range_true);
							Fill1DHist(trklen_reco_midp_RecoInel, range_reco);
							Fill1DHist(dtrklen_midp_RecoInel, range_reco-range_true);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //mid:p

		//[6]MID:midmu
		if (kMIDmu) { //mid:mu
			n_midmu++;
			Fill1DHist(trklen_true_midmu_NoCut, range_true);
			Fill1DHist(trklen_reco_midmu_NoCut, range_reco);
			Fill1DHist(dtrklen_midmu_NoCut, range_reco-range_true);
			if (IsPandoraSlice) { //pandora
				n_midmu_pan++;
				Fill1DHist(trklen_true_midmu_PanS, range_true);
				Fill1DHist(trklen_reco_midmu_PanS, range_reco);
				Fill1DHist(dtrklen_midmu_PanS, range_reco-range_true);
				if (IsCaloSize) { //calosz
					n_midmu_calsz++;
					Fill1DHist(trklen_true_midmu_CaloSz, range_true);
					Fill1DHist(trklen_reco_midmu_CaloSz, range_reco);
					Fill1DHist(dtrklen_midmu_CaloSz, range_reco-range_true);
					if (IsBQ) { //bq
						n_midmu_bq++;
						Fill1DHist(trklen_true_midmu_BQ, range_true);
						Fill1DHist(trklen_reco_midmu_BQ, range_reco);
						Fill1DHist(dtrklen_midmu_BQ, range_reco-range_true);
						if (IsRecoInEL) { //reco inel
							n_midmu_recoinel++;
							Fill1DHist(trklen_true_midmu_RecoInel, range_true);
							Fill1DHist(trklen_reco_midmu_RecoInel, range_reco);
							Fill1DHist(dtrklen_midmu_RecoInel, range_reco-range_true);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //mid:mu

		//[7]MID:mideg
		if (kMIDeg) { //mid:eg
			n_mideg++;
			Fill1DHist(trklen_true_mideg_NoCut, range_true);
			Fill1DHist(trklen_reco_mideg_NoCut, range_reco);
			Fill1DHist(dtrklen_mideg_NoCut, range_reco-range_true);
			if (IsPandoraSlice) { //pandora
				n_mideg_pan++;
				Fill1DHist(trklen_true_mideg_PanS, range_true);
				Fill1DHist(trklen_reco_mideg_PanS, range_reco);
				Fill1DHist(dtrklen_mideg_PanS, range_reco-range_true);
				if (IsCaloSize) { //calosz
					n_mideg_calsz++;
					Fill1DHist(trklen_true_mideg_CaloSz, range_true);
					Fill1DHist(trklen_reco_mideg_CaloSz, range_reco);
					Fill1DHist(dtrklen_mideg_CaloSz, range_reco-range_true);
					if (IsBQ) { //bq
						n_mideg_bq++;
						Fill1DHist(trklen_true_mideg_BQ, range_true);
						Fill1DHist(trklen_reco_mideg_BQ, range_reco);
						Fill1DHist(dtrklen_mideg_BQ, range_reco-range_true);
						if (IsRecoInEL) { //reco inel
							n_mideg_recoinel++;
							Fill1DHist(trklen_true_mideg_RecoInel, range_true);
							Fill1DHist(trklen_reco_mideg_RecoInel, range_reco);
							Fill1DHist(dtrklen_mideg_RecoInel, range_reco-range_true);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //mid:eg

		//[8]MID:midother
		if (kMIDother) { //mid:other
			n_midother++;
			Fill1DHist(trklen_true_midother_NoCut, range_true);
			Fill1DHist(trklen_reco_midother_NoCut, range_reco);
			Fill1DHist(dtrklen_midother_NoCut, range_reco-range_true);
			if (IsPandoraSlice) { //pandora
				n_midother_pan++;
				Fill1DHist(trklen_true_midother_PanS, range_true);
				Fill1DHist(trklen_reco_midother_PanS, range_reco);
				Fill1DHist(dtrklen_midother_PanS, range_reco-range_true);
				if (IsCaloSize) { //calosz
					n_midother_calsz++;
					Fill1DHist(trklen_true_midother_CaloSz, range_true);
					Fill1DHist(trklen_reco_midother_CaloSz, range_reco);
					Fill1DHist(dtrklen_midother_CaloSz, range_reco-range_true);
					if (IsBQ) { //bq
						n_midother_bq++;
						Fill1DHist(trklen_true_midother_BQ, range_true);
						Fill1DHist(trklen_reco_midother_BQ, range_reco);
						Fill1DHist(dtrklen_midother_BQ, range_reco-range_true);
						if (IsRecoInEL) { //reco inel
							n_midother_recoinel++;
							Fill1DHist(trklen_true_midother_RecoInel, range_true);
							Fill1DHist(trklen_reco_midother_RecoInel, range_reco);
							Fill1DHist(dtrklen_midother_RecoInel, range_reco-range_true);
						} //reco inel
					} //bq
				} //calosz
			} //pandora
		} //mid:other
		//countings -------------------------------------------------------//

		//ntrklen --------------------------------------------------------------------//
		if (IsPandoraSlice&&IsCaloSize&&IsBQ) {  
			Fill1DHist(ntrklen_BQ, range_reco/csda_val_spec);
			if (kinel) Fill1DHist(ntrklen_inel_BQ, range_reco/csda_val_spec);
			if (kel) Fill1DHist(ntrklen_el_BQ, range_reco/csda_val_spec);
			if (kMIDcosmic) Fill1DHist(ntrklen_midcosmic_BQ, range_reco/csda_val_spec);
			if (kMIDpi) Fill1DHist(ntrklen_midpi_BQ, range_reco/csda_val_spec);
			if (kMIDp) Fill1DHist(ntrklen_midp_BQ, range_reco/csda_val_spec);
			if (kMIDmu) Fill1DHist(ntrklen_midmu_BQ, range_reco/csda_val_spec);
			if (kMIDeg) Fill1DHist(ntrklen_mideg_BQ, range_reco/csda_val_spec);
			if (kMIDother) Fill1DHist(ntrklen_midother_BQ, range_reco/csda_val_spec); 
		}

		//true trklen vs ke ------------------------------------------------------------------//
		if (!beamtrk_Eng->empty()) { //if true container not empty
			for (int hk=0; hk<(int)beamtrk_Eng->size()-1; ++hk) { //loop over true hits
				//double thisZ=beamtrk_z->at(hk);
				double thisLen=true_trklen_accum[hk];
				double this_incE = 1000.*beamtrk_Eng->at(hk); //MeV
				if (kinel) trklen_ke_true_inel->Fill(thisLen, this_incE);
				if (kel) trklen_ke_true_el->Fill(thisLen, this_incE);
			} //loop over true hits
		} //if true container not empty

		//thin-slice method ----------------------------------------------------//
		//true slice ID
		//true_sliceID = int(true_endz/thinslicewidth);
		true_sliceID = int(range_true/thinslicewidth); //HY:Make sure size of true_trk_len vector is !0, otherwise lost truth info
		if (true_sliceID < 0) true_sliceID = -1;
		if (true_endz < 0) true_sliceID = -1; //hy added
		if (true_sliceID >= nthinslices) true_sliceID = nthinslices;
		cout<<"true_endz:"<<true_endz<<"  true_sliceID:"<<true_sliceID<<endl;
		//cout<<"primtrk_range_true->empty():"<<primtrk_range_true->empty()
		   // <<" primtrk_range->empty():"<<primtrk_range->empty()<<endl;

		//evt selection cuts
		bool PassCuts_INT=false; //all bq cut+reco inel cut
		bool PassCuts_INC=false; //all bq cut
		if (IsPandoraSlice&&IsCaloSize) {
			//reco slice ID
			//reco_sliceID = int(primtrk_range->at(0)/thinslicewidth);
			//reco_sliceID = int(reco_endz/thinslicewidth);
			reco_sliceID = int(range_reco/thinslicewidth);
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
				//for (size_t ih=0; ih<zreco_rawindex.size(); ++ih) {
				for (size_t ih=0; ih<primtrk_dedx->size(); ++ih) {
					//double this_calo_z=zreco_widreco[ih].second;
					//int this_sliceID = int(this_calo_z/thinslicewidth);
					//double this_calo_len=zreco_lenreco[ih].second;
					//double this_dE=zreco_de[ih].second;
					int this_sliceID = int(reco_trklen_accum[ih]/thinslicewidth);
					//ke_reco-=this_dE;
					ke_reco-=EDept.at(ih);

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

				//TVector3 pt0(zreco_xreco[0].second,
						//zreco_yreco[0].second,
						//zreco_xreco[0].first); //ST
				//TVector3 pt1(zreco_xreco[-1+zreco_rawindex.size()].second,
						//zreco_yreco[-1+zreco_rawindex.size()].second,
						//zreco_xreco[-1+zreco_rawindex.size()].first);
				TVector3 pt0(primtrk_hitx->at(0), primtrk_hity->at(0), primtrk_hitz->at(0));
				TVector3 pt1(primtrk_hitx->at(-1+primtrk_hitx->size()), primtrk_hity->at(-1+primtrk_hity->size()), primtrk_hitz->at(-1+primtrk_hitz->size()));
				TVector3 dir = pt1 - pt0;
				dir = dir.Unit();
				reco_AngCorr->Fill(dir.Z());
			} //calo & beam_match

			//truth KEs
			if (!beamtrk_Eng->empty()) { //if true container not empty
				std::vector<std::vector<double>> vincE_true(nthinslices);
				for (int hk=0; hk<(int)beamtrk_Eng->size()-1; ++hk) { //loop over true hits
					//double thisZ=beamtrk_z->at(hk);
					//int this_sliceID = int(thisZ/thinslicewidth);
					double thisLen=true_trklen_accum[hk];
					int this_sliceID = int(thisLen/thinslicewidth);
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





	
	} //main entry loop

	//save results -------//
	uf.SaveHistograms();
	CalcXS(uf);
	SaveHistograms();
         
        //counting -- summary -----------------------------------------------------//
	cout<<"\nn_tot:"<<n_tot<<endl;
	cout<<"n_el:"<<n_el<<endl;
	cout<<"n_inel:"<<n_inel<<endl;
	cout<<"n_midcosmic:"<<n_midcosmic<<endl;
	cout<<"n_midpi:"<<n_midpi<<endl;
	cout<<"n_midp:"<<n_midp<<endl;
	cout<<"n_midmu:"<<n_midmu<<endl;
	cout<<"n_mideg:"<<n_mideg<<endl;
	cout<<"n_midother"<<n_midother<<endl;
	cout<<"n_diff:"<<n_tot-(n_el+n_inel+n_midcosmic+n_midpi+n_midp+n_midmu+n_mideg+n_midother)<<endl;

	cout<<"\nn_pan_tot:"<<n_pan_tot<<endl;
	cout<<"n_el_pan:"<<n_el_pan<<endl;
	cout<<"n_inel_pan:"<<n_inel_pan<<endl;
	cout<<"n_midcosmic_pan:"<<n_midcosmic_pan<<endl;
	cout<<"n_midpi_pan:"<<n_midpi_pan<<endl;
	cout<<"n_midp_pan:"<<n_midp_pan<<endl;
	cout<<"n_midmu_pan:"<<n_midmu_pan<<endl;
	cout<<"n_mideg_pan:"<<n_mideg_pan<<endl;
	cout<<"n_midother_pan"<<n_midother_pan<<endl;
	cout<<"n_diff_pan:"<<n_pan_tot-(n_el_pan+n_inel_pan+n_midcosmic_pan+n_midpi_pan+n_midp_pan+n_midmu_pan+n_mideg_pan+n_midother_pan)<<endl;

	cout<<"\nn_calsz_tot:"<<n_calsz_tot<<endl;
	cout<<"n_el_calsz:"<<n_el_calsz<<endl;
	cout<<"n_inel_calsz:"<<n_inel_calsz<<endl;
	cout<<"n_midcosmic_calsz:"<<n_midcosmic_calsz<<endl;
	cout<<"n_midpi_calsz:"<<n_midpi_calsz<<endl;
	cout<<"n_midp_calsz:"<<n_midp_calsz<<endl;
	cout<<"n_midmu_calsz:"<<n_midmu_calsz<<endl;
	cout<<"n_mideg_calsz:"<<n_mideg_calsz<<endl;
	cout<<"n_midother_calsz"<<n_midother_calsz<<endl;
	cout<<"n_diff_calsz:"<<n_calsz_tot-(n_el_calsz+n_inel_calsz+n_midcosmic_calsz+n_midpi_calsz+n_midp_calsz+n_midmu_calsz+n_mideg_calsz+n_midother_calsz)<<endl;

	cout<<"\nn_bq_tot:"<<n_bq_tot<<endl;
	cout<<"n_el_bq:"<<n_el_bq<<endl;
	cout<<"n_inel_bq:"<<n_inel_bq<<endl;
	cout<<"n_midcosmic_bq:"<<n_midcosmic_bq<<endl;
	cout<<"n_midpi_bq:"<<n_midpi_bq<<endl;
	cout<<"n_midp_bq:"<<n_midp_bq<<endl;
	cout<<"n_midmu_bq:"<<n_midmu_bq<<endl;
	cout<<"n_mideg_bq:"<<n_mideg_bq<<endl;
	cout<<"n_midother_bq"<<n_midother_bq<<endl;
	cout<<"n_diff_bq:"<<n_bq_tot-(n_el_bq+n_inel_bq+n_midcosmic_bq+n_midpi_bq+n_midp_bq+n_midmu_bq+n_mideg_bq+n_midother_bq)<<endl;

	cout<<"\nn_recoinel_tot:"<<n_recoinel_tot<<endl;
	cout<<"n_el_recoinel:"<<n_el_recoinel<<endl;
	cout<<"n_inel_recoinel:"<<n_inel_recoinel<<endl;
	cout<<"n_midcosmic_recoinel:"<<n_midcosmic_recoinel<<endl;
	cout<<"n_midpi_recoinel:"<<n_midpi_recoinel<<endl;
	cout<<"n_midp_recoinel:"<<n_midp_recoinel<<endl;
	cout<<"n_midmu_recoinel:"<<n_midmu_recoinel<<endl;
	cout<<"n_mideg_recoinel:"<<n_mideg_recoinel<<endl;
	cout<<"n_midother_recoinel"<<n_midother_recoinel<<endl;
	cout<<"n_diff_recoinel:"<<n_recoinel_tot-(n_el_recoinel+n_inel_recoinel+n_midcosmic_recoinel+n_midpi_recoinel+n_midp_recoinel+n_midmu_recoinel+n_mideg_recoinel+n_midother_recoinel)<<endl;





}
