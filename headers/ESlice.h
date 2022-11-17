#include "TGraphErrors.h"
#include "TVector3.h"
//#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "util.h"
#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "./Unfold.h"

//#include "BetheBloch.h"

//Basic config. -----------------------------------------------------//
std::string fOutputFileName;
TFile *outputFile;
void SetOutputFileName(std::string name){fOutputFileName = name;};
void BookHistograms();
void SaveHistograms();

//XS histograms ----------------------//
//int reco_sliceID;
//int true_sliceID;
TH1D *reco_incE[nthinslices];
TH1D *true_incE[nthinslices];
TH1D *reco_AngCorr;
TH1D *true_AngCorr;

TH1D *h_truesliceid_all;
TH1D *h_true_st_sliceid_all;
//TH1D *h_truesliceid_uf;
TH1D *h_truesliceid_cuts;
TH1D *h_true_st_sliceid_cuts;
TH1D *h_truesliceid_inelastic_all;
//TH1D *h_truesliceid_inelastic_uf;
TH1D *h_truesliceid_inelastic_cuts;

//reco inc
TH1D *h_recosliceid_allevts_cuts;
TH1D *h_recosliceid_allevts_cuts_inel;
TH1D *h_recosliceid_allevts_cuts_el;
TH1D *h_recosliceid_allevts_cuts_midcosmic;
TH1D *h_recosliceid_allevts_cuts_midpi;
TH1D *h_recosliceid_allevts_cuts_midp;
TH1D *h_recosliceid_allevts_cuts_midmu;
TH1D *h_recosliceid_allevts_cuts_mideg;
TH1D *h_recosliceid_allevts_cuts_midother;

//true inc
TH1D *h_truesliceid_allevts_cuts;
TH1D *h_truesliceid_allevts_cuts_inel;
TH1D *h_truesliceid_allevts_cuts_el;
TH1D *h_truesliceid_allevts_cuts_midcosmic;
TH1D *h_truesliceid_allevts_cuts_midpi;
TH1D *h_truesliceid_allevts_cuts_midp;
TH1D *h_truesliceid_allevts_cuts_midmu;
TH1D *h_truesliceid_allevts_cuts_mideg;
TH1D *h_truesliceid_allevts_cuts_midother;


//reco st inc
TH1D *h_reco_st_sliceid_allevts_cuts;
TH1D *h_reco_st_sliceid_allevts_cuts_inel;
TH1D *h_reco_st_sliceid_allevts_cuts_el;
TH1D *h_reco_st_sliceid_allevts_cuts_midcosmic;
TH1D *h_reco_st_sliceid_allevts_cuts_midpi;
TH1D *h_reco_st_sliceid_allevts_cuts_midp;
TH1D *h_reco_st_sliceid_allevts_cuts_midmu;
TH1D *h_reco_st_sliceid_allevts_cuts_mideg;
TH1D *h_reco_st_sliceid_allevts_cuts_midother;

//true st inc
TH1D *h_true_st_sliceid_allevts_cuts;
TH1D *h_true_st_sliceid_allevts_cuts_inel;
TH1D *h_true_st_sliceid_allevts_cuts_el;
TH1D *h_true_st_sliceid_allevts_cuts_midcosmic;
TH1D *h_true_st_sliceid_allevts_cuts_midpi;
TH1D *h_true_st_sliceid_allevts_cuts_midp;
TH1D *h_true_st_sliceid_allevts_cuts_midmu;
TH1D *h_true_st_sliceid_allevts_cuts_mideg;
TH1D *h_true_st_sliceid_allevts_cuts_midother;


//reco int
TH1D *h_recosliceid_cuts;
TH1D *h_reco_st_sliceid_cuts;
TH1D *h_recosliceid_inelastic_cuts;

TH1D *h_recosliceid_recoinelastic_cuts;
TH1D *h_recosliceid_recoinelastic_cuts_inel;
TH1D *h_recosliceid_recoinelastic_cuts_el;
TH1D *h_recosliceid_recoinelastic_cuts_midcosmic;
TH1D *h_recosliceid_recoinelastic_cuts_midpi;
TH1D *h_recosliceid_recoinelastic_cuts_midp;
TH1D *h_recosliceid_recoinelastic_cuts_midmu;
TH1D *h_recosliceid_recoinelastic_cuts_mideg;
TH1D *h_recosliceid_recoinelastic_cuts_midother;

//true int
TH1D *h_truesliceid_recoinelastic_cuts;
TH1D *h_truesliceid_recoinelastic_cuts_inel;
TH1D *h_truesliceid_recoinelastic_cuts_el;
TH1D *h_truesliceid_recoinelastic_cuts_midcosmic;
TH1D *h_truesliceid_recoinelastic_cuts_midpi;
TH1D *h_truesliceid_recoinelastic_cuts_midp;
TH1D *h_truesliceid_recoinelastic_cuts_midmu;
TH1D *h_truesliceid_recoinelastic_cuts_mideg;
TH1D *h_truesliceid_recoinelastic_cuts_midother;


double true_incidents[nthinslices+2];
double true_st_incidents[nthinslices+2];
double true_interactions[nthinslices+2];

//Histograms for basic parameters ---------//
//reco x, y, z [after SCE corr]
TH1D *reco_startX_sce; 
TH1D *reco_startY_sce;
TH1D *reco_startZ_sce;

TH1D *hdeltaX;
TH1D *hdeltaX_inel;
TH1D *hdeltaX_el;
TH1D *hdeltaX_midcosmic;
TH1D *hdeltaX_midpi;
TH1D *hdeltaX_midp;
TH1D *hdeltaX_midmu;
TH1D *hdeltaX_mideg;
TH1D *hdeltaX_midother;

TH1D *hdeltaY;
TH1D *hdeltaY_inel;
TH1D *hdeltaY_el;
TH1D *hdeltaY_midcosmic;
TH1D *hdeltaY_midpi;
TH1D *hdeltaY_midp;
TH1D *hdeltaY_midmu;
TH1D *hdeltaY_mideg;
TH1D *hdeltaY_midother;

TH1D *hdeltaZ;
TH1D *hdeltaZ_inel;
TH1D *hdeltaZ_el;
TH1D *hdeltaZ_midcosmic;
TH1D *hdeltaZ_midpi;
TH1D *hdeltaZ_midp;
TH1D *hdeltaZ_midmu;
TH1D *hdeltaZ_mideg;
TH1D *hdeltaZ_midother;

TH1D *hdeltaXY; //similiar to XY cut
TH1D *hdeltaXY_inel;
TH1D *hdeltaXY_el;
TH1D *hdeltaXY_midcosmic;
TH1D *hdeltaXY_midpi;
TH1D *hdeltaXY_midp;
TH1D *hdeltaXY_midmu;
TH1D *hdeltaXY_mideg;
TH1D *hdeltaXY_midother;







//cosine_theta
TH1D *reco_cosineTheta;
TH1D *reco_cosineTheta_inel;
TH1D *reco_cosineTheta_el;
TH1D *reco_cosineTheta_midcosmic;
TH1D *reco_cosineTheta_midpi;
TH1D *reco_cosineTheta_midp;
TH1D *reco_cosineTheta_midmu;
TH1D *reco_cosineTheta_mideg;
TH1D *reco_cosineTheta_midother;

TH1D *reco_cosineTheta_Pos; //apply position cut
TH1D *reco_cosineTheta_Pos_inel; //apply position cut
TH1D *reco_cosineTheta_Pos_el;
TH1D *reco_cosineTheta_Pos_midcosmic;
TH1D *reco_cosineTheta_Pos_midpi;
TH1D *reco_cosineTheta_Pos_midp;
TH1D *reco_cosineTheta_Pos_midmu;
TH1D *reco_cosineTheta_Pos_mideg;
TH1D *reco_cosineTheta_Pos_midother;


//dE/dx vs rr related histograms ------------//
TH2D *rr_dedx_recostop;
TH2D *rr_dedx_truestop;
TH2D *rangereco_dedxreco_TrueInEL;
TH2D *rangereco_dedxreco_TrueEL;

//KE calc using reco stopping protons ----------------------------------------------------//
TH1D *KE_ff_recostop;
TH1D *KE_calo_recostop;
TH1D *KE_rrange_recostop;
TH1D *KE_rrange2_recostop;
TH1D *KE_range_recostop;
TH1D *KE_simide_recostop;
TH1D *dKE_range_ff_recostop;
TH1D *dKE_calo_ff_recostop;
TH1D *dKE_rrange_ff_recostop;
TH1D *dKE_rrange2_ff_recostop;
TH2D *KE_range_ff_recostop;
TH2D *KE_range_calo_recostop;


//KE Truth info
TH1D *KEtrue_Beam;
TH1D *KEtrue_Beam_inel;
TH1D *KEtrue_ff_inel;
//TH2D *KEtrue_z_inel;
//TH2D *KEtrue_range_inel;
//TH2D *true_z_range_inel;
//TH3D *true_z_range_KEtrue_inel;

//TH2D *KEbb_truetrklen_all;
//TH2D *KEbb_truetrklen_inel;
//TH2D *KEbb_recotrklen_all;
//TH2D *KEbb_recotrklen_inel;

TH2D *KEcalo_truetrklen_all;
TH2D *KEcalo_truetrklen_inel;
TH2D *KEcalo_recotrklen_all;
TH2D *KEcalo_recotrklen_inel;

//True Trklen Patch Dist.
TH1D *true_trklen_patch_all;

//KE reco info
TH2D *KEreco_z_inel;
TH2D *KEreco_range_inel;
TH2D *reco_z_range_inel;

TH2D *dEdx_range_inel;
TH2D *dx_range_inel;
TH2D *dE_range_inel;

//Track length histograms -----------------------------------------------------------------//
//true range
//no cut
TH1D *ke_true_inel_NoCut;
TH1D *ke_true_el_NoCut;
TH1D *ke_true_midcosmic_NoCut;
TH1D *ke_true_midpi_NoCut;
TH1D *ke_true_midp_NoCut;
TH1D *ke_true_midmu_NoCut;
TH1D *ke_true_mideg_NoCut;
TH1D *ke_true_midother_NoCut;

//pandora slice cut
TH1D *ke_true_inel_PanS;
TH1D *ke_true_el_PanS;
TH1D *ke_true_midcosmic_PanS;
TH1D *ke_true_midpi_PanS;
TH1D *ke_true_midp_PanS;
TH1D *ke_true_midmu_PanS;
TH1D *ke_true_mideg_PanS;
TH1D *ke_true_midother_PanS;

//calosz cut
TH1D *ke_true_inel_CaloSz;
TH1D *ke_true_el_CaloSz;
TH1D *ke_true_midcosmic_CaloSz;
TH1D *ke_true_midpi_CaloSz;
TH1D *ke_true_midp_CaloSz;
TH1D *ke_true_midmu_CaloSz;
TH1D *ke_true_mideg_CaloSz;
TH1D *ke_true_midother_CaloSz;

//beam quality cut
TH1D *ke_true_inel_BQ;
TH1D *ke_true_el_BQ;
TH1D *ke_true_midcosmic_BQ;
TH1D *ke_true_midpi_BQ;
TH1D *ke_true_midp_BQ;
TH1D *ke_true_midmu_BQ;
TH1D *ke_true_mideg_BQ;
TH1D *ke_true_midother_BQ;

//RecoInel cut
TH1D *ke_true_inel_RecoInel;
TH1D *ke_true_el_RecoInel;
TH1D *ke_true_midcosmic_RecoInel;
TH1D *ke_true_midpi_RecoInel;
TH1D *ke_true_midp_RecoInel;
TH1D *ke_true_midmu_RecoInel;
TH1D *ke_true_mideg_RecoInel;
TH1D *ke_true_midother_RecoInel;


//RecoEl cut
TH1D *ke_true_inel_RecoEl;
TH1D *ke_true_el_RecoEl;
TH1D *ke_true_midcosmic_RecoEl;
TH1D *ke_true_midpi_RecoEl;
TH1D *ke_true_midp_RecoEl;
TH1D *ke_true_midmu_RecoEl;
TH1D *ke_true_mideg_RecoEl;
TH1D *ke_true_midother_RecoEl;



//reco range
//no cut
TH1D *ke_reco_inel_NoCut;
TH1D *ke_reco_el_NoCut;
TH1D *ke_reco_midcosmic_NoCut;
TH1D *ke_reco_midpi_NoCut;
TH1D *ke_reco_midp_NoCut;
TH1D *ke_reco_midmu_NoCut;
TH1D *ke_reco_mideg_NoCut;
TH1D *ke_reco_midother_NoCut;

//pandora slice cut
TH1D *ke_reco_inel_PanS;
TH1D *ke_reco_el_PanS;
TH1D *ke_reco_midcosmic_PanS;
TH1D *ke_reco_midpi_PanS;
TH1D *ke_reco_midp_PanS;
TH1D *ke_reco_midmu_PanS;
TH1D *ke_reco_mideg_PanS;
TH1D *ke_reco_midother_PanS;

//calosz cut
TH1D *ke_reco_inel_CaloSz;
TH1D *ke_reco_el_CaloSz;
TH1D *ke_reco_midcosmic_CaloSz;
TH1D *ke_reco_midpi_CaloSz;
TH1D *ke_reco_midp_CaloSz;
TH1D *ke_reco_midmu_CaloSz;
TH1D *ke_reco_mideg_CaloSz;
TH1D *ke_reco_midother_CaloSz;

//beam quality cut
//end-point
TH1D *ke_reco_BQ;
TH1D *ke_reco_inel_BQ;
TH1D *ke_reco_el_BQ;
TH1D *dke_reco_el_BQ;
TH1D *kedept_reco_el_BQ;
TH1D *ke_reco_midcosmic_BQ;
TH1D *ke_reco_midpi_BQ;
TH1D *ke_reco_midp_BQ;
TH1D *ke_reco_midmu_BQ;
TH1D *ke_reco_mideg_BQ;
TH1D *ke_reco_midother_BQ;
//start
TH1D *keff_reco_BQ;
TH1D *keff_reco_inel_BQ;
TH1D *keff_reco_el_BQ;
TH1D *keff_reco_midcosmic_BQ;
TH1D *keff_reco_midpi_BQ;
TH1D *keff_reco_midp_BQ;
TH1D *keff_reco_midmu_BQ;
TH1D *keff_reco_mideg_BQ;
TH1D *keff_reco_midother_BQ;



//RecoInel cut
TH1D *ke_reco_RecoInel;
TH1D *ke_reco_inel_RecoInel;
TH1D *ke_reco_el_RecoInel;
TH1D *ke_reco_midcosmic_RecoInel;
TH1D *ke_reco_midpi_RecoInel;
TH1D *ke_reco_midp_RecoInel;
TH1D *ke_reco_midmu_RecoInel;
TH1D *ke_reco_mideg_RecoInel;
TH1D *ke_reco_midother_RecoInel;

//RecoEl cut
TH1D *ke_reco_RecoEl;
TH1D *ke_reco_inel_RecoEl;
TH1D *ke_reco_el_RecoEl;
TH1D *ke_reco_midcosmic_RecoEl;
TH1D *ke_reco_midpi_RecoEl;
TH1D *ke_reco_midp_RecoEl;
TH1D *ke_reco_midmu_RecoEl;
TH1D *ke_reco_mideg_RecoEl;
TH1D *ke_reco_midother_RecoEl;

//MidP-rich cut
TH1D *ke_reco_MidP;
TH1D *ke_reco_inel_MidP;
TH1D *ke_reco_el_MidP;
TH1D *ke_reco_midcosmic_MidP;
TH1D *ke_reco_midpi_MidP;
TH1D *ke_reco_midp_MidP;
TH1D *ke_reco_midmu_MidP;
TH1D *ke_reco_mideg_MidP;
TH1D *ke_reco_midother_MidP;


//(reco-truth)ke
//no cut
TH1D *dke;
TH1D *dke_inel_NoCut;
TH1D *dke_el_NoCut;
TH1D *dke_midcosmic_NoCut;
TH1D *dke_midpi_NoCut;
TH1D *dke_midp_NoCut;
TH1D *dke_midmu_NoCut;
TH1D *dke_mideg_NoCut;
TH1D *dke_midother_NoCut;

//pandora slice cut
TH1D *dke_inel_PanS;
TH1D *dke_el_PanS;
TH1D *dke_midcosmic_PanS;
TH1D *dke_midpi_PanS;
TH1D *dke_midp_PanS;
TH1D *dke_midmu_PanS;
TH1D *dke_mideg_PanS;
TH1D *dke_midother_PanS;

//calosz cut
TH1D *dke_inel_CaloSz;
TH1D *dke_el_CaloSz;
TH1D *dke_midcosmic_CaloSz;
TH1D *dke_midpi_CaloSz;
TH1D *dke_midp_CaloSz;
TH1D *dke_midmu_CaloSz;
TH1D *dke_mideg_CaloSz;
TH1D *dke_midother_CaloSz;

//beam quality cut
TH1D *dke_inel_BQ;
TH1D *dke_el_BQ;
TH1D *dke_midcosmic_BQ;
TH1D *dke_midpi_BQ;
TH1D *dke_midp_BQ;
TH1D *dke_midmu_BQ;
TH1D *dke_mideg_BQ;
TH1D *dke_midother_BQ;

//RecoInel cut
TH1D *dke_inel_RecoInel;
TH1D *dke_el_RecoInel;
TH1D *dke_midcosmic_RecoInel;
TH1D *dke_midpi_RecoInel;
TH1D *dke_midp_RecoInel;
TH1D *dke_midmu_RecoInel;
TH1D *dke_mideg_RecoInel;
TH1D *dke_midother_RecoInel;

//RecoEl cut
TH1D *dke_inel_RecoEl;
TH1D *dke_el_RecoEl;
TH1D *dke_midcosmic_RecoEl;
TH1D *dke_midpi_RecoEl;
TH1D *dke_midp_RecoEl;
TH1D *dke_midmu_RecoEl;
TH1D *dke_mideg_RecoEl;
TH1D *dke_midother_RecoEl;

//ntrklen histograms ----------//
//beam quality cut
TH1D *ntrklen_BQ;
TH1D *ntrklen_inel_BQ;
TH1D *ntrklen_el_BQ;
TH1D *ntrklen_midcosmic_BQ;
TH1D *ntrklen_midpi_BQ;
TH1D *ntrklen_midp_BQ;
TH1D *ntrklen_midmu_BQ;
TH1D *ntrklen_mideg_BQ;
TH1D *ntrklen_midother_BQ;

//trklen vs ke
TH2D *trklen_ke_true_inel;
TH2D *trklen_ke_true_el;

//truth inc/int
TH1D *h_true_incidents;
TH1D *h_true_st_incidents;
TH1D *h_true_interactions;

//Reco E-dept [Inel]
TH2D *reco_dedx_trklen_inel;
TH2D *reco_de_trklen_inel;
TH2D *reco_dx_trklen_inel;

//KEff sansity check
//TH2D *h2d_R_kE1st;

void BookHistograms() { //BookHistograms
	outputFile = TFile::Open(fOutputFileName.c_str(), "recreate"); //open output file

	//XS histograms -------------------------------------------------------------------------------------------------------------------------------------------------------//
	for (int i = 0; i<nthinslices; ++i){
		reco_incE[i] = new TH1D(Form("reco_incE_%d",i),Form("Reco incident energy, %.1f < length < %.1f (cm)",i*thinslicewidth, (i+1)*thinslicewidth), nbinse, 0, 1200.);
		true_incE[i] = new TH1D(Form("true_incE_%d",i),Form("True incident energy, %.1f < length < %.1f (cm)",i*thinslicewidth, (i+1)*thinslicewidth), nbinse, 0, 1200.);
		reco_incE[i]->Sumw2();
		true_incE[i]->Sumw2();
	}

	reco_AngCorr = new TH1D("reco_AngCorr","Reco angle correction", 100, 0, 1.);
	true_AngCorr = new TH1D("true_AngCorr","true angle correction", 100, 0, 1.);
	reco_AngCorr->Sumw2();
	true_AngCorr->Sumw2();

	h_truesliceid_all = new TH1D("h_truesliceid_all","h_truesliceid_all;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_true_st_sliceid_all = new TH1D("h_true_st_sliceid_all","h_true_st_sliceid_all;True Start SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_cuts = new TH1D("h_truesliceid_cuts","h_truesliceid_cuts;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_true_st_sliceid_cuts = new TH1D("h_true_st_sliceid_cuts","h_true_st_sliceid_cuts;True Start SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_inelastic_all = new TH1D("h_truesliceid_inelastic_all","h_truesliceid_inelastic_all;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_inelastic_cuts = new TH1D("h_truesliceid_inelastic_cuts","h_truesliceid_inelastic_cuts;True SliceID", nthinslices + 2, -1, nthinslices + 1);

	h_recosliceid_allevts_cuts = new TH1D("h_recosliceid_allevts_cuts","h_recosliceid_allevts_cuts;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts_inel = new TH1D("h_recosliceid_allevts_cuts_inel","h_recosliceid_allevts_cuts_inel;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts_el = new TH1D("h_recosliceid_allevts_cuts_el","h_recosliceid_allevts_cuts_el;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts_midcosmic = new TH1D("h_recosliceid_allevts_cuts_midcosmic","h_recosliceid_allevts_cuts_midcosmic;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts_midpi = new TH1D("h_recosliceid_allevts_cuts_midpi","h_recosliceid_allevts_cuts_midpi;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts_midp = new TH1D("h_recosliceid_allevts_cuts_midp","h_recosliceid_allevts_cuts_midp;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts_midmu = new TH1D("h_recosliceid_allevts_cuts_midmu","h_recosliceid_allevts_cuts_midmu;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts_mideg = new TH1D("h_recosliceid_allevts_cuts_mideg","h_recosliceid_allevts_cuts_mideg;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts_midother = new TH1D("h_recosliceid_allevts_cuts_midother","h_recosliceid_allevts_cuts_midother;True SliceID", nthinslices + 2, -1, nthinslices + 1);

	h_truesliceid_allevts_cuts = new TH1D("h_truesliceid_allevts_cuts","h_truesliceid_allevts_cuts;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_allevts_cuts_inel = new TH1D("h_truesliceid_allevts_cuts_inel","h_truesliceid_allevts_cuts_inel;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_allevts_cuts_el = new TH1D("h_truesliceid_allevts_cuts_el","h_truesliceid_allevts_cuts_el;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_allevts_cuts_midcosmic = new TH1D("h_truesliceid_allevts_cuts_midcosmic","h_truesliceid_allevts_cuts_midcosmic;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_allevts_cuts_midpi = new TH1D("h_truesliceid_allevts_cuts_midpi","h_truesliceid_allevts_cuts_midpi;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_allevts_cuts_midp = new TH1D("h_truesliceid_allevts_cuts_midp","h_truesliceid_allevts_cuts_midp;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_allevts_cuts_midmu = new TH1D("h_truesliceid_allevts_cuts_midmu","h_truesliceid_allevts_cuts_midmu;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_allevts_cuts_mideg = new TH1D("h_truesliceid_allevts_cuts_mideg","h_truesliceid_allevts_cuts_mideg;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_allevts_cuts_midother = new TH1D("h_truesliceid_allevts_cuts_midother","h_truesliceid_allevts_cuts_midother;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);

	h_reco_st_sliceid_allevts_cuts = new TH1D("h_reco_st_sliceid_allevts_cuts","h_reco_st_sliceid_allevts_cuts;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_reco_st_sliceid_allevts_cuts_inel = new TH1D("h_reco_st_sliceid_allevts_cuts_inel","h_reco_st_sliceid_allevts_cuts_inel;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_reco_st_sliceid_allevts_cuts_el = new TH1D("h_reco_st_sliceid_allevts_cuts_el","h_reco_st_sliceid_allevts_cuts_el;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_reco_st_sliceid_allevts_cuts_midcosmic = new TH1D("h_reco_st_sliceid_allevts_cuts_midcosmic","h_reco_st_sliceid_allevts_cuts_midcosmic;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_reco_st_sliceid_allevts_cuts_midpi = new TH1D("h_reco_st_sliceid_allevts_cuts_midpi","h_reco_st_sliceid_allevts_cuts_midpi;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_reco_st_sliceid_allevts_cuts_midp = new TH1D("h_reco_st_sliceid_allevts_cuts_midp","h_reco_st_sliceid_allevts_cuts_midp;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_reco_st_sliceid_allevts_cuts_midmu = new TH1D("h_reco_st_sliceid_allevts_cuts_midmu","h_reco_st_sliceid_allevts_cuts_midmu;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_reco_st_sliceid_allevts_cuts_mideg = new TH1D("h_reco_st_sliceid_allevts_cuts_mideg","h_reco_st_sliceid_allevts_cuts_mideg;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_reco_st_sliceid_allevts_cuts_midother = new TH1D("h_reco_st_sliceid_allevts_cuts_midother","h_reco_st_sliceid_allevts_cuts_midother;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);

	h_true_st_sliceid_allevts_cuts = new TH1D("h_true_st_sliceid_allevts_cuts","h_true_st_sliceid_allevts_cuts;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_true_st_sliceid_allevts_cuts_inel = new TH1D("h_true_st_sliceid_allevts_cuts_inel","h_true_st_sliceid_allevts_cuts_inel;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_true_st_sliceid_allevts_cuts_el = new TH1D("h_true_st_sliceid_allevts_cuts_el","h_true_st_sliceid_allevts_cuts_el;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_true_st_sliceid_allevts_cuts_midcosmic = new TH1D("h_true_st_sliceid_allevts_cuts_midcosmic","h_true_st_sliceid_allevts_cuts_midcosmic;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_true_st_sliceid_allevts_cuts_midpi = new TH1D("h_true_st_sliceid_allevts_cuts_midpi","h_true_st_sliceid_allevts_cuts_midpi;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_true_st_sliceid_allevts_cuts_midp = new TH1D("h_true_st_sliceid_allevts_cuts_midp","h_true_st_sliceid_allevts_cuts_midp;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_true_st_sliceid_allevts_cuts_midmu = new TH1D("h_true_st_sliceid_allevts_cuts_midmu","h_true_st_sliceid_allevts_cuts_midmu;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_true_st_sliceid_allevts_cuts_mideg = new TH1D("h_true_st_sliceid_allevts_cuts_mideg","h_true_st_sliceid_allevts_cuts_mideg;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_true_st_sliceid_allevts_cuts_midother = new TH1D("h_true_st_sliceid_allevts_cuts_midother","h_true_st_sliceid_allevts_cuts_midother;True SliceID", nthinslices + 2, -1, nthinslices + 1);

	h_recosliceid_cuts = new TH1D("h_recosliceid_cuts","h_recosliceid_cuts;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_reco_st_sliceid_cuts = new TH1D("h_reco_st_sliceid_cuts","h_reco_st_sliceid_cuts;Reco Start SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_inelastic_cuts = new TH1D("h_recosliceid_inelastic_cuts","h_recosliceid_inelastic_cuts;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_recoinelastic_cuts = new TH1D("h_recosliceid_recoinelastic_cuts","h_recosliceid_recoinelastic_cuts;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_recoinelastic_cuts_inel = new TH1D("h_recosliceid_recoinelastic_cuts_inel","h_recosliceid_recoinelastic_cuts_inel;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_recoinelastic_cuts_el = new TH1D("h_recosliceid_recoinelastic_cuts_el","h_recosliceid_recoinelastic_cuts_el;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_recoinelastic_cuts_midcosmic = new TH1D("h_recosliceid_recoinelastic_cuts_midcosmic","h_recosliceid_recoinelastic_cuts_midcosmic;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_recoinelastic_cuts_midpi = new TH1D("h_recosliceid_recoinelastic_cuts_midpi","h_recosliceid_recoinelastic_cuts_midpi;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_recoinelastic_cuts_midp = new TH1D("h_recosliceid_recoinelastic_cuts_midp","h_recosliceid_recoinelastic_cuts_midp;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_recoinelastic_cuts_midmu = new TH1D("h_recosliceid_recoinelastic_cuts_midmu","h_recosliceid_recoinelastic_cuts_midmu;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_recoinelastic_cuts_mideg = new TH1D("h_recosliceid_recoinelastic_cuts_mideg","h_recosliceid_recoinelastic_cuts_mideg;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_recoinelastic_cuts_midother = new TH1D("h_recosliceid_recoinelastic_cuts_midother","h_recosliceid_recoinelastic_cuts_midother;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);

	h_truesliceid_recoinelastic_cuts = new TH1D("h_truesliceid_recoinelastic_cuts","h_truesliceid_recoinelastic_cuts;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_recoinelastic_cuts_inel = new TH1D("h_truesliceid_recoinelastic_cuts_inel","h_truesliceid_recoinelastic_cuts_inel;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_recoinelastic_cuts_el = new TH1D("h_truesliceid_recoinelastic_cuts_el","h_truesliceid_recoinelastic_cuts_el;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_recoinelastic_cuts_midcosmic = new TH1D("h_truesliceid_recoinelastic_cuts_midcosmic","h_truesliceid_recoinelastic_cuts_midcosmic;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_recoinelastic_cuts_midpi = new TH1D("h_truesliceid_recoinelastic_cuts_midpi","h_truesliceid_recoinelastic_cuts_midpi;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_recoinelastic_cuts_midp = new TH1D("h_truesliceid_recoinelastic_cuts_midp","h_truesliceid_recoinelastic_cuts_midp;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_recoinelastic_cuts_midmu = new TH1D("h_truesliceid_recoinelastic_cuts_midmu","h_truesliceid_recoinelastic_cuts_midmu;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_recoinelastic_cuts_mideg = new TH1D("h_truesliceid_recoinelastic_cuts_mideg","h_truesliceid_recoinelastic_cuts_mideg;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_recoinelastic_cuts_midother = new TH1D("h_truesliceid_recoinelastic_cuts_midother","h_truesliceid_recoinelastic_cuts_midother;True SliceID", nthinslices + 2, -1, nthinslices + 1);


	h_truesliceid_all->Sumw2();
	h_true_st_sliceid_all->Sumw2();
	h_truesliceid_cuts->Sumw2();
	h_true_st_sliceid_cuts->Sumw2();
	h_truesliceid_inelastic_all->Sumw2();
	h_truesliceid_inelastic_cuts->Sumw2();

	h_recosliceid_allevts_cuts->Sumw2();
	h_recosliceid_allevts_cuts_inel->Sumw2();
	h_recosliceid_allevts_cuts_el->Sumw2();
	h_recosliceid_allevts_cuts_midcosmic->Sumw2();
	h_recosliceid_allevts_cuts_midpi->Sumw2();
	h_recosliceid_allevts_cuts_midp->Sumw2();
	h_recosliceid_allevts_cuts_midmu->Sumw2();
	h_recosliceid_allevts_cuts_mideg->Sumw2();
	h_recosliceid_allevts_cuts_midother->Sumw2();


	h_truesliceid_allevts_cuts->Sumw2();
	h_truesliceid_allevts_cuts_inel->Sumw2();
	h_truesliceid_allevts_cuts_el->Sumw2();
	h_truesliceid_allevts_cuts_midcosmic->Sumw2();
	h_truesliceid_allevts_cuts_midpi->Sumw2();
	h_truesliceid_allevts_cuts_midp->Sumw2();
	h_truesliceid_allevts_cuts_midmu->Sumw2();
	h_truesliceid_allevts_cuts_mideg->Sumw2();
	h_truesliceid_allevts_cuts_midother->Sumw2();

	h_reco_st_sliceid_allevts_cuts->Sumw2();
	h_reco_st_sliceid_allevts_cuts_inel->Sumw2();
	h_reco_st_sliceid_allevts_cuts_el->Sumw2();
	h_reco_st_sliceid_allevts_cuts_midcosmic->Sumw2();
	h_reco_st_sliceid_allevts_cuts_midpi->Sumw2();
	h_reco_st_sliceid_allevts_cuts_midp->Sumw2();
	h_reco_st_sliceid_allevts_cuts_midmu->Sumw2();
	h_reco_st_sliceid_allevts_cuts_mideg->Sumw2();
	h_reco_st_sliceid_allevts_cuts_midother->Sumw2();

	h_true_st_sliceid_allevts_cuts->Sumw2();
	h_true_st_sliceid_allevts_cuts_inel->Sumw2();
	h_true_st_sliceid_allevts_cuts_el->Sumw2();
	h_true_st_sliceid_allevts_cuts_midcosmic->Sumw2();
	h_true_st_sliceid_allevts_cuts_midpi->Sumw2();
	h_true_st_sliceid_allevts_cuts_midp->Sumw2();
	h_true_st_sliceid_allevts_cuts_midmu->Sumw2();
	h_true_st_sliceid_allevts_cuts_mideg->Sumw2();
	h_true_st_sliceid_allevts_cuts_midother->Sumw2();


	h_recosliceid_cuts->Sumw2();
	h_reco_st_sliceid_cuts->Sumw2();
	h_recosliceid_inelastic_cuts->Sumw2();

	h_recosliceid_recoinelastic_cuts->Sumw2();
	h_recosliceid_recoinelastic_cuts_inel->Sumw2();
	h_recosliceid_recoinelastic_cuts_el->Sumw2();
	h_recosliceid_recoinelastic_cuts_midcosmic->Sumw2();
	h_recosliceid_recoinelastic_cuts_midpi->Sumw2();
	h_recosliceid_recoinelastic_cuts_midp->Sumw2();
	h_recosliceid_recoinelastic_cuts_midmu->Sumw2();
	h_recosliceid_recoinelastic_cuts_mideg->Sumw2();
	h_recosliceid_recoinelastic_cuts_midother->Sumw2();


	h_truesliceid_recoinelastic_cuts->Sumw2();
	h_truesliceid_recoinelastic_cuts_inel->Sumw2();
	h_truesliceid_recoinelastic_cuts_el->Sumw2();
	h_truesliceid_recoinelastic_cuts_midcosmic->Sumw2();
	h_truesliceid_recoinelastic_cuts_midpi->Sumw2();
	h_truesliceid_recoinelastic_cuts_midp->Sumw2();
	h_truesliceid_recoinelastic_cuts_midmu->Sumw2();
	h_truesliceid_recoinelastic_cuts_mideg->Sumw2();
	h_truesliceid_recoinelastic_cuts_midother->Sumw2();


	//for (int i = 0; i<nthinslices; ++i){
	for (int i = 0; i<nthinslices+2; ++i){
		true_interactions[i] = 0;
		true_incidents[i] = 0;
		true_st_incidents[i] = 0;
	}

	//Histograms for basic parameters ----------------------------------------------------------------------------------------------------//
	reco_startX_sce = new TH1D(Form("reco_startX_sce"), Form("reco_startX_sce"), 100, -80, 20);  reco_startX_sce->Sumw2();
	reco_startY_sce = new TH1D(Form("reco_startY_sce"), Form("reco_startY_sce"), 100, 350, 500); reco_startY_sce->Sumw2();
	reco_startZ_sce = new TH1D(Form("reco_startZ_sce"), Form("reco_startZ_sce"), 100, -5, 10);   reco_startZ_sce->Sumw2();

	hdeltaX = new TH1D(Form("hdeltaX"), Form("#Deltax/#sigma_{x}"), 100, -10, 10);  hdeltaX->Sumw2();
	hdeltaX_inel = new TH1D(Form("hdeltaX_inel"), Form("#Deltax/#sigma_{x}"), 100, -10, 10);  hdeltaX_inel->Sumw2();
	hdeltaX_el = new TH1D(Form("hdeltaX_el"), Form("#Deltax/#sigma_{x}"), 100, -10, 10);  hdeltaX_el->Sumw2();
	hdeltaX_midcosmic = new TH1D(Form("hdeltaX_midcosmic"), Form("#Deltax/#sigma_{x}"), 100, -10, 10);  hdeltaX_midcosmic->Sumw2();
	hdeltaX_midpi = new TH1D(Form("hdeltaX_midpi"), Form("#Deltax/#sigma_{x}"), 100, -10, 10);  hdeltaX_midpi->Sumw2();
	hdeltaX_midp = new TH1D(Form("hdeltaX_midp"), Form("#Deltax/#sigma_{x}"), 100, -10, 10);  hdeltaX_midp->Sumw2();
	hdeltaX_midmu = new TH1D(Form("hdeltaX_midmu"), Form("#Deltax/#sigma_{x}"), 100, -10, 10);  hdeltaX_midmu->Sumw2();
	hdeltaX_mideg = new TH1D(Form("hdeltaX_mideg"), Form("#Deltax/#sigma_{x}"), 100, -10, 10);  hdeltaX_mideg->Sumw2();
	hdeltaX_midother = new TH1D(Form("hdeltaX_midother"), Form("#Deltax/#sigma_{x}"), 100, -10, 10);  hdeltaX_midother->Sumw2();

	hdeltaY = new TH1D(Form("hdeltaY"), Form("#Deltay/#sigma_{y}"), 100, -10, 10);  hdeltaY->Sumw2();
	hdeltaY_inel = new TH1D(Form("hdeltaY_inel"), Form("#Deltay/#sigma_{y}"), 100, -10, 10);  hdeltaY_inel->Sumw2();
	hdeltaY_el = new TH1D(Form("hdeltaY_el"), Form("#Deltay/#sigma_{y}"), 100, -10, 10);  hdeltaY_el->Sumw2();
	hdeltaY_midcosmic = new TH1D(Form("hdeltaY_midcosmic"), Form("#Deltay/#sigma_{y}"), 100, -10, 10);  hdeltaY_midcosmic->Sumw2();
	hdeltaY_midpi = new TH1D(Form("hdeltaY_midpi"), Form("#Deltay/#sigma_{y}"), 100, -10, 10);  hdeltaY_midpi->Sumw2();
	hdeltaY_midp = new TH1D(Form("hdeltaY_midp"), Form("#Deltay/#sigma_{y}"), 100, -10, 10);  hdeltaY_midp->Sumw2();
	hdeltaY_midmu = new TH1D(Form("hdeltaY_midmu"), Form("#Deltay/#sigma_{y}"), 100, -10, 10);  hdeltaY_midmu->Sumw2();
	hdeltaY_mideg = new TH1D(Form("hdeltaY_mideg"), Form("#Deltay/#sigma_{y}"), 100, -10, 10);  hdeltaY_mideg->Sumw2();
	hdeltaY_midother = new TH1D(Form("hdeltaY_midother"), Form("#Deltay/#sigma_{y}"), 100, -10, 10);  hdeltaY_midother->Sumw2();

	hdeltaZ = new TH1D(Form("hdeltaZ"), Form("#Deltaz/#sigma_{z}"), 100, -10, 10);  hdeltaZ->Sumw2();
	hdeltaZ_inel = new TH1D(Form("hdeltaZ_inel"), Form("#Deltaz/#sigma_{z}"), 100, -10, 10);  hdeltaZ_inel->Sumw2();
	hdeltaZ_el = new TH1D(Form("hdeltaZ_el"), Form("#Deltaz/#sigma_{z}"), 100, -10, 10);  hdeltaZ_el->Sumw2();
	hdeltaZ_midcosmic = new TH1D(Form("hdeltaZ_midcosmic"), Form("#Deltaz/#sigma_{z}"), 100, -10, 10);  hdeltaZ_midcosmic->Sumw2();
	hdeltaZ_midpi = new TH1D(Form("hdeltaZ_midpi"), Form("#Deltaz/#sigma_{z}"), 100, -10, 10);  hdeltaZ_midpi->Sumw2();
	hdeltaZ_midp = new TH1D(Form("hdeltaZ_midp"), Form("#Deltaz/#sigma_{z}"), 100, -10, 10);  hdeltaZ_midp->Sumw2();
	hdeltaZ_midmu = new TH1D(Form("hdeltaZ_midmu"), Form("#Deltaz/#sigma_{z}"), 100, -10, 10);  hdeltaZ_midmu->Sumw2();
	hdeltaZ_mideg = new TH1D(Form("hdeltaZ_mideg"), Form("#Deltaz/#sigma_{z}"), 100, -10, 10);  hdeltaZ_mideg->Sumw2();
	hdeltaZ_midother = new TH1D(Form("hdeltaZ_midother"), Form("#Deltaz/#sigma_{z}"), 100, -10, 10);  hdeltaZ_midother->Sumw2();

	hdeltaXY = new TH1D(Form("hdeltaXY"), Form("Sqrt((#Deltax/#sigma_{x})^2+(#Deltay/#sigma_{y})^2)"), 100, -10, 10);  hdeltaXY->Sumw2();
	hdeltaXY_inel = new TH1D(Form("hdeltaXY_inel"), Form("Sqrt((#Deltax/#sigma_{x})^2+(#Deltay/#sigma_{y})^2)"), 100, -10, 10);  hdeltaXY_inel->Sumw2();
	hdeltaXY_el = new TH1D(Form("hdeltaXY_el"), Form("Sqrt((#Deltax/#sigma_{x})^2+(#Deltay/#sigma_{y})^2)"), 100, -10, 10);  hdeltaXY_el->Sumw2();
	hdeltaXY_midcosmic = new TH1D(Form("hdeltaXY_midcosmic"), Form("Sqrt((#Deltax/#sigma_{x})^2+(#Deltay/#sigma_{y})^2)"), 100, -10, 10);  hdeltaXY_midcosmic->Sumw2();
	hdeltaXY_midpi = new TH1D(Form("hdeltaXY_midpi"), Form("Sqrt((#Deltax/#sigma_{x})^2+(#Deltay/#sigma_{y})^2)"), 100, -10, 10);  hdeltaXY_midpi->Sumw2();
	hdeltaXY_midp = new TH1D(Form("hdeltaXY_midp"), Form("Sqrt((#Deltax/#sigma_{x})^2+(#Deltay/#sigma_{y})^2)"), 100, -10, 10);  hdeltaXY_midp->Sumw2();
	hdeltaXY_midmu = new TH1D(Form("hdeltaXY_midmu"), Form("Sqrt((#Deltax/#sigma_{x})^2+(#Deltay/#sigma_{y})^2)"), 100, -10, 10);  hdeltaXY_midmu->Sumw2();
	hdeltaXY_mideg = new TH1D(Form("hdeltaXY_mideg"), Form("Sqrt((#Deltax/#sigma_{x})^2+(#Deltay/#sigma_{y})^2)"), 100, -10, 10);  hdeltaXY_mideg->Sumw2();
	hdeltaXY_midother = new TH1D(Form("hdeltaXY_midother"), Form("Sqrt((#Deltax/#sigma_{x})^2+(#Deltay/#sigma_{y})^2)"), 100, -10, 10);  hdeltaXY_midother->Sumw2();


	//Histograms for cosineTheta --------------------------------------------------------------//
	int n_cosine=100;
	//double cosine_min=0.9;
	double cosine_min=0;
	double cosine_max=1.0;
	reco_cosineTheta = new TH1D("reco_cosineTheta","", n_cosine, cosine_min, cosine_max);	reco_cosineTheta->Sumw2();

	reco_cosineTheta_inel=new TH1D("reco_cosineTheta_inel", "", n_cosine, cosine_min, cosine_max); reco_cosineTheta_inel->Sumw2();
	reco_cosineTheta_el=new TH1D("reco_cosineTheta_el", "", n_cosine, cosine_min, cosine_max); reco_cosineTheta_el->Sumw2();
	reco_cosineTheta_midcosmic=new TH1D("reco_cosineTheta_midcosmic", "", n_cosine, cosine_min, cosine_max); reco_cosineTheta_midcosmic->Sumw2();
	reco_cosineTheta_midpi=new TH1D("reco_cosineTheta_midpi", "", n_cosine, cosine_min, cosine_max); reco_cosineTheta_midpi->Sumw2();
	reco_cosineTheta_midp=new TH1D("reco_cosineTheta_midp","", n_cosine, cosine_min, cosine_max); reco_cosineTheta_midp->Sumw2();
	reco_cosineTheta_midmu=new TH1D("reco_cosineTheta_midmu","", n_cosine, cosine_min, cosine_max); reco_cosineTheta_midmu->Sumw2();
	reco_cosineTheta_mideg=new TH1D("reco_cosineTheta_mideg","", n_cosine, cosine_min, cosine_max); reco_cosineTheta_mideg->Sumw2();
	reco_cosineTheta_midother=new TH1D("reco_cosineTheta_midother","", n_cosine, cosine_min, cosine_max); reco_cosineTheta_midother->Sumw2();


	reco_cosineTheta_Pos=new TH1D("reco_cosineTheta_Pos", "", n_cosine, cosine_min, cosine_max); reco_cosineTheta_Pos->Sumw2();
	reco_cosineTheta_Pos_inel=new TH1D("reco_cosineTheta_Pos_inel", "", n_cosine, cosine_min, cosine_max); reco_cosineTheta_Pos_inel->Sumw2();
	reco_cosineTheta_Pos_el=new TH1D("reco_cosineTheta_Pos_el", "", n_cosine, cosine_min, cosine_max); reco_cosineTheta_Pos_el->Sumw2();
	reco_cosineTheta_Pos_midcosmic=new TH1D("reco_cosineTheta_Pos_midcosmic", "", n_cosine, cosine_min, cosine_max); reco_cosineTheta_Pos_midcosmic->Sumw2();
	reco_cosineTheta_Pos_midpi=new TH1D("reco_cosineTheta_Pos_midpi", "", n_cosine, cosine_min, cosine_max); reco_cosineTheta_Pos_midpi->Sumw2();
	reco_cosineTheta_Pos_midp=new TH1D("reco_cosineTheta_Pos_midp","", n_cosine, cosine_min, cosine_max); reco_cosineTheta_Pos_midp->Sumw2();
	reco_cosineTheta_Pos_midmu=new TH1D("reco_cosineTheta_Pos_midmu","", n_cosine, cosine_min, cosine_max); reco_cosineTheta_Pos_midmu->Sumw2();
	reco_cosineTheta_Pos_mideg=new TH1D("reco_cosineTheta_Pos_mideg","", n_cosine, cosine_min, cosine_max); reco_cosineTheta_Pos_mideg->Sumw2();
	reco_cosineTheta_Pos_midother=new TH1D("reco_cosineTheta_Pos_midother","", n_cosine, cosine_min, cosine_max); reco_cosineTheta_Pos_midother->Sumw2();

	//dE/dx vs rr related histograms ---------------------------------------------------------//
	rr_dedx_recostop=new TH2D("rr_dedx_recostop","", 240,0,120, 300,0, 30);
	rr_dedx_truestop=new TH2D("rr_dedx_truestop","", 240,0,120, 300,0, 30);
	rangereco_dedxreco_TrueInEL=new TH2D("rangereco_dedxreco_TrueInEL","",200,0,100,100,0,50);
	rangereco_dedxreco_TrueEL=new TH2D("rangereco_dedxreco_TrueEL","",200,0,100,100,0,50);

	//KE calc using reco stopping protons ----------------------------------------------------//
	int n_ke=160;
	double ke_min=-800;
	double ke_max=800;
	//int n_ke=500;
	//float ke_min=0;
	//float ke_max=1000;
	KE_ff_recostop=new TH1D("KE_ff_recostop","", n_ke, ke_min, ke_max); KE_ff_recostop->Sumw2();
	KE_calo_recostop=new TH1D("KE_calo_recostop","",n_ke, ke_min, ke_max); KE_calo_recostop->Sumw2();
	KE_rrange_recostop=new TH1D("KE_rrange_recostop", "", n_ke, ke_min, ke_max); KE_rrange_recostop->Sumw2();
	KE_rrange2_recostop=new TH1D("KE_rrange2_recostop", "", n_ke, ke_min, ke_max); KE_rrange2_recostop->Sumw2();
	KE_range_recostop=new TH1D("KE_range_recostop", "", n_ke, ke_min, ke_max); KE_range_recostop->Sumw2();
	KE_simide_recostop=new TH1D("KE_simide_recostop","",n_ke, ke_min, ke_max); KE_simide_recostop->Sumw2();
	dKE_range_ff_recostop=new TH1D("dKE_range_ff_recostop","",200,-100,100); dKE_range_ff_recostop->Sumw2();
	dKE_calo_ff_recostop=new TH1D("dKE_calo_ff_recostop","",200,-100,100); dKE_calo_ff_recostop->Sumw2();
	dKE_rrange_ff_recostop=new TH1D("dKE_rrange_ff_recostop","",200,-100,100); dKE_rrange_ff_recostop->Sumw2();
	dKE_rrange2_ff_recostop=new TH1D("dKE_rrange2_ff_recostop","",200,-100,100); dKE_rrange2_ff_recostop->Sumw2();
	KE_range_ff_recostop=new TH2D("KE_range_ff_recostop","",n_ke, ke_min, ke_max, n_ke, ke_min, ke_max); 
	KE_range_ff_recostop->GetXaxis()->SetTitle("KE_{range} [MeV]"); KE_range_ff_recostop->GetYaxis()->SetTitle("KE_{ff} [MeV]"); 
	KE_range_ff_recostop->Sumw2();
	KE_range_calo_recostop=new TH2D("KE_range_calo_recostop","",n_ke, ke_min, ke_max, n_ke, ke_min, ke_max); 
	KE_range_calo_recostop->GetXaxis()->SetTitle("KE_{range} [MeV]"); KE_range_calo_recostop->GetYaxis()->SetTitle("KE_{calo} [MeV]"); 
	KE_range_calo_recostop->Sumw2();


	//KEtruth info
	int n_range=300;
	float range_min=0;
	float range_max=150;
	KEtrue_Beam=new TH1D("KEtrue_Beam","",n_ke, ke_min, ke_max); KEtrue_Beam->Sumw2();
	KEtrue_Beam_inel=new TH1D("KEtrue_Beam_inel","",n_ke, ke_min, ke_max); KEtrue_Beam_inel->Sumw2();
	KEtrue_ff_inel=new TH1D("KEtrue_ff_inel", "", n_ke, ke_min, ke_max); KEtrue_ff_inel->Sumw2();
	//KEtrue_z_inel=new TH2D("KEtrue_z_inel","", n_range, range_min, range_max, n_ke,ke_min,ke_max);
	//KEtrue_range_inel=new TH2D("KEtrue_range_inel","", n_range, range_min, range_max, n_ke,ke_min,ke_max);
	//true_z_range_inel=new TH2D("true_z_range_inel","", n_range, range_min, range_max,  n_range, range_min, range_max);
	//true_z_range_KEtrue_inel=new TH3D("true_z_range_KEtrue_inel","", n_range, range_min, range_max, n_range, range_min, range_max, n_ke,ke_min,ke_max);
	//true_z_range_KEtrue_inel->GetXaxis()->SetTitle("z [cm]");
	//true_z_range_KEtrue_inel->GetYaxis()->SetTitle("range [cm]");
	//true_z_range_KEtrue_inel->GetZaxis()->SetTitle("KE [MeV]");

	//KEreco info
	KEreco_z_inel=new TH2D("KEreco_z_inel","", n_range, range_min, range_max, n_ke,ke_min,ke_max);
	KEreco_range_inel=new TH2D("KEreco_range_inel","", n_range, range_min, range_max, n_ke,ke_min,ke_max);
	reco_z_range_inel=new TH2D("reco_z_range_inel","", n_range, range_min, range_max,  n_range, range_min, range_max);

	dEdx_range_inel=new TH2D("dEdx_range_inel","", n_range, range_min, range_max,1000,0,100);
	dx_range_inel=new TH2D("dx_range_inel","", n_range, range_min, range_max, 10000, 0, 100);
	dE_range_inel=new TH2D("dE_range_inel","", n_range, range_min, range_max, 5000, 0, 500);


	//trklen -----------------------------------------------------------------------------------------------------//
        int n_trklen=34;
        double trklen_min=-4;
        double trklen_max=132;

	//ke -----------------------------------------------------------------------------------------------------//
	int nke=160;
	double kemin=-800;
	double kemax=800;	


	//truth range
	//no cut
	ke_true_inel_NoCut = new TH1D("ke_true_inel_NoCut","",n_ke, ke_min, ke_max); ke_true_inel_NoCut->Sumw2();
	ke_true_el_NoCut = new TH1D("ke_true_el_NoCut","",n_ke, ke_min, ke_max); ke_true_el_NoCut->Sumw2();
	ke_true_midcosmic_NoCut = new TH1D("ke_true_midcosmic_NoCut","",n_ke, ke_min, ke_max); ke_true_midcosmic_NoCut->Sumw2();
	ke_true_midpi_NoCut = new TH1D("ke_true_midpi_NoCut","",n_ke, ke_min, ke_max); ke_true_midpi_NoCut->Sumw2();
	ke_true_midp_NoCut = new TH1D("ke_true_midp_NoCut","",n_ke, ke_min, ke_max); ke_true_midp_NoCut->Sumw2();
	ke_true_midmu_NoCut = new TH1D("ke_true_midmu_NoCut","",n_ke, ke_min, ke_max); ke_true_midmu_NoCut->Sumw2();
	ke_true_mideg_NoCut = new TH1D("ke_true_mideg_NoCut","",n_ke, ke_min, ke_max); ke_true_mideg_NoCut->Sumw2();
	ke_true_midother_NoCut = new TH1D("ke_true_midother_NoCut","",n_ke, ke_min, ke_max); ke_true_midother_NoCut->Sumw2();

	//pandora cut
	ke_true_inel_PanS = new TH1D("ke_true_inel_PanS","",n_ke, ke_min, ke_max); ke_true_inel_PanS->Sumw2();
	ke_true_el_PanS = new TH1D("ke_true_el_PanS","",n_ke, ke_min, ke_max); ke_true_el_PanS->Sumw2();
	ke_true_midcosmic_PanS = new TH1D("ke_true_midcosmic_PanS","",n_ke, ke_min, ke_max); ke_true_midcosmic_PanS->Sumw2();
	ke_true_midpi_PanS = new TH1D("ke_true_midpi_PanS","",n_ke, ke_min, ke_max); ke_true_midpi_PanS->Sumw2();
	ke_true_midp_PanS = new TH1D("ke_true_midp_PanS","",n_ke, ke_min, ke_max); ke_true_midp_PanS->Sumw2();
	ke_true_midmu_PanS = new TH1D("ke_true_midmu_PanS","",n_ke, ke_min, ke_max); ke_true_midmu_PanS->Sumw2();
	ke_true_mideg_PanS = new TH1D("ke_true_mideg_PanS","",n_ke, ke_min, ke_max); ke_true_mideg_PanS->Sumw2();
	ke_true_midother_PanS = new TH1D("ke_true_midother_PanS","",n_ke, ke_min, ke_max); ke_true_midother_PanS->Sumw2();

	//CaloSz
	ke_true_inel_CaloSz = new TH1D("ke_true_inel_CaloSz","",n_ke, ke_min, ke_max); ke_true_inel_CaloSz->Sumw2();
	ke_true_el_CaloSz = new TH1D("ke_true_el_CaloSz","",n_ke, ke_min, ke_max); ke_true_el_CaloSz->Sumw2();
	ke_true_midcosmic_CaloSz = new TH1D("ke_true_midcosmic_CaloSz","",n_ke, ke_min, ke_max); ke_true_midcosmic_CaloSz->Sumw2();
	ke_true_midpi_CaloSz = new TH1D("ke_true_midpi_CaloSz","",n_ke, ke_min, ke_max); ke_true_midpi_CaloSz->Sumw2();
	ke_true_midp_CaloSz = new TH1D("ke_true_midp_CaloSz","",n_ke, ke_min, ke_max); ke_true_midp_CaloSz->Sumw2();
	ke_true_midmu_CaloSz = new TH1D("ke_true_midmu_CaloSz","",n_ke, ke_min, ke_max); ke_true_midmu_CaloSz->Sumw2();
	ke_true_mideg_CaloSz = new TH1D("ke_true_mideg_CaloSz","",n_ke, ke_min, ke_max); ke_true_mideg_CaloSz->Sumw2();
	ke_true_midother_CaloSz = new TH1D("ke_true_midother_CaloSz","",n_ke, ke_min, ke_max); ke_true_midother_CaloSz->Sumw2();

	//beam quality
	ke_true_inel_BQ = new TH1D("ke_true_inel_BQ","",n_ke, ke_min, ke_max); ke_true_inel_BQ->Sumw2();
	ke_true_el_BQ = new TH1D("ke_true_el_BQ","",n_ke, ke_min, ke_max); ke_true_el_BQ->Sumw2();
	ke_true_midcosmic_BQ = new TH1D("ke_true_midcosmic_BQ","",n_ke, ke_min, ke_max); ke_true_midcosmic_BQ->Sumw2();
	ke_true_midpi_BQ = new TH1D("ke_true_midpi_BQ","",n_ke, ke_min, ke_max); ke_true_midpi_BQ->Sumw2();
	ke_true_midp_BQ = new TH1D("ke_true_midp_BQ","",n_ke, ke_min, ke_max); ke_true_midp_BQ->Sumw2();
	ke_true_midmu_BQ = new TH1D("ke_true_midmu_BQ","",n_ke, ke_min, ke_max); ke_true_midmu_BQ->Sumw2();
	ke_true_mideg_BQ = new TH1D("ke_true_mideg_BQ","",n_ke, ke_min, ke_max); ke_true_mideg_BQ->Sumw2();
	ke_true_midother_BQ = new TH1D("ke_true_midother_BQ","",n_ke, ke_min, ke_max); ke_true_midother_BQ->Sumw2();

	//reco inel cut
	ke_true_inel_RecoInel = new TH1D("ke_true_inel_RecoInel","",n_ke, ke_min, ke_max); ke_true_inel_RecoInel->Sumw2();
	ke_true_el_RecoInel = new TH1D("ke_true_el_RecoInel","",n_ke, ke_min, ke_max); ke_true_el_RecoInel->Sumw2();
	ke_true_midcosmic_RecoInel = new TH1D("ke_true_midcosmic_RecoInel","",n_ke, ke_min, ke_max); ke_true_midcosmic_RecoInel->Sumw2();
	ke_true_midpi_RecoInel = new TH1D("ke_true_midpi_RecoInel","",n_ke, ke_min, ke_max); ke_true_midpi_RecoInel->Sumw2();
	ke_true_midp_RecoInel = new TH1D("ke_true_midp_RecoInel","",n_ke, ke_min, ke_max); ke_true_midp_RecoInel->Sumw2();
	ke_true_midmu_RecoInel = new TH1D("ke_true_midmu_RecoInel","",n_ke, ke_min, ke_max); ke_true_midmu_RecoInel->Sumw2();
	ke_true_mideg_RecoInel = new TH1D("ke_true_mideg_RecoInel","",n_ke, ke_min, ke_max); ke_true_mideg_RecoInel->Sumw2();
	ke_true_midother_RecoInel = new TH1D("ke_true_midother_RecoInel","",n_ke, ke_min, ke_max); ke_true_midother_RecoInel->Sumw2();


	//reco el cut
	ke_true_inel_RecoEl = new TH1D("ke_true_inel_RecoEl","",n_ke, ke_min, ke_max); ke_true_inel_RecoEl->Sumw2();
	ke_true_el_RecoEl = new TH1D("ke_true_el_RecoEl","",n_ke, ke_min, ke_max); ke_true_el_RecoEl->Sumw2();
	ke_true_midcosmic_RecoEl = new TH1D("ke_true_midcosmic_RecoEl","",n_ke, ke_min, ke_max); ke_true_midcosmic_RecoEl->Sumw2();
	ke_true_midpi_RecoEl = new TH1D("ke_true_midpi_RecoEl","",n_ke, ke_min, ke_max); ke_true_midpi_RecoEl->Sumw2();
	ke_true_midp_RecoEl = new TH1D("ke_true_midp_RecoEl","",n_ke, ke_min, ke_max); ke_true_midp_RecoEl->Sumw2();
	ke_true_midmu_RecoEl = new TH1D("ke_true_midmu_RecoEl","",n_ke, ke_min, ke_max); ke_true_midmu_RecoEl->Sumw2();
	ke_true_mideg_RecoEl = new TH1D("ke_true_mideg_RecoEl","",n_ke, ke_min, ke_max); ke_true_mideg_RecoEl->Sumw2();
	ke_true_midother_RecoEl = new TH1D("ke_true_midother_RecoEl","",n_ke, ke_min, ke_max); ke_true_midother_RecoEl->Sumw2();

	//reco range
	//no cut
	ke_reco_inel_NoCut = new TH1D("ke_reco_inel_NoCut","",n_ke, ke_min, ke_max); ke_reco_inel_NoCut->Sumw2();
	ke_reco_el_NoCut = new TH1D("ke_reco_el_NoCut","",n_ke, ke_min, ke_max); ke_reco_el_NoCut->Sumw2();
	ke_reco_midcosmic_NoCut = new TH1D("ke_reco_midcosmic_NoCut","",n_ke, ke_min, ke_max); ke_reco_midcosmic_NoCut->Sumw2();
	ke_reco_midpi_NoCut = new TH1D("ke_reco_midpi_NoCut","",n_ke, ke_min, ke_max); ke_reco_midpi_NoCut->Sumw2();
	ke_reco_midp_NoCut = new TH1D("ke_reco_midp_NoCut","",n_ke, ke_min, ke_max); ke_reco_midp_NoCut->Sumw2();
	ke_reco_midmu_NoCut = new TH1D("ke_reco_midmu_NoCut","",n_ke, ke_min, ke_max); ke_reco_midmu_NoCut->Sumw2();
	ke_reco_mideg_NoCut = new TH1D("ke_reco_mideg_NoCut","",n_ke, ke_min, ke_max); ke_reco_mideg_NoCut->Sumw2();
	ke_reco_midother_NoCut = new TH1D("ke_reco_midother_NoCut","",n_ke, ke_min, ke_max); ke_reco_midother_NoCut->Sumw2();

	//pandora cut
	ke_reco_inel_PanS = new TH1D("ke_reco_inel_PanS","",n_ke, ke_min, ke_max); ke_reco_inel_PanS->Sumw2();
	ke_reco_el_PanS = new TH1D("ke_reco_el_PanS","",n_ke, ke_min, ke_max); ke_reco_el_PanS->Sumw2();
	ke_reco_midcosmic_PanS = new TH1D("ke_reco_midcosmic_PanS","",n_ke, ke_min, ke_max); ke_reco_midcosmic_PanS->Sumw2();
	ke_reco_midpi_PanS = new TH1D("ke_reco_midpi_PanS","",n_ke, ke_min, ke_max); ke_reco_midpi_PanS->Sumw2();
	ke_reco_midp_PanS = new TH1D("ke_reco_midp_PanS","",n_ke, ke_min, ke_max); ke_reco_midp_PanS->Sumw2();
	ke_reco_midmu_PanS = new TH1D("ke_reco_midmu_PanS","",n_ke, ke_min, ke_max); ke_reco_midmu_PanS->Sumw2();
	ke_reco_mideg_PanS = new TH1D("ke_reco_mideg_PanS","",n_ke, ke_min, ke_max); ke_reco_mideg_PanS->Sumw2();
	ke_reco_midother_PanS = new TH1D("ke_reco_midother_PanS","",n_ke, ke_min, ke_max); ke_reco_midother_PanS->Sumw2();

	//CaloSz
	ke_reco_inel_CaloSz = new TH1D("ke_reco_inel_CaloSz","",n_ke, ke_min, ke_max); ke_reco_inel_CaloSz->Sumw2();
	ke_reco_el_CaloSz = new TH1D("ke_reco_el_CaloSz","",n_ke, ke_min, ke_max); ke_reco_el_CaloSz->Sumw2();
	ke_reco_midcosmic_CaloSz = new TH1D("ke_reco_midcosmic_CaloSz","",n_ke, ke_min, ke_max); ke_reco_midcosmic_CaloSz->Sumw2();
	ke_reco_midpi_CaloSz = new TH1D("ke_reco_midpi_CaloSz","",n_ke, ke_min, ke_max); ke_reco_midpi_CaloSz->Sumw2();
	ke_reco_midp_CaloSz = new TH1D("ke_reco_midp_CaloSz","",n_ke, ke_min, ke_max); ke_reco_midp_CaloSz->Sumw2();
	ke_reco_midmu_CaloSz = new TH1D("ke_reco_midmu_CaloSz","",n_ke, ke_min, ke_max); ke_reco_midmu_CaloSz->Sumw2();
	ke_reco_mideg_CaloSz = new TH1D("ke_reco_mideg_CaloSz","",n_ke, ke_min, ke_max); ke_reco_mideg_CaloSz->Sumw2();
	ke_reco_midother_CaloSz = new TH1D("ke_reco_midother_CaloSz","",n_ke, ke_min, ke_max); ke_reco_midother_CaloSz->Sumw2();

	//beam quality
	ke_reco_BQ = new TH1D("ke_reco_BQ","",n_ke, ke_min, ke_max); ke_reco_BQ->Sumw2();
	ke_reco_inel_BQ = new TH1D("ke_reco_inel_BQ","",n_ke, ke_min, ke_max); ke_reco_inel_BQ->Sumw2();
	ke_reco_el_BQ = new TH1D("ke_reco_el_BQ","",n_ke, ke_min, ke_max); ke_reco_el_BQ->Sumw2();
	dke_reco_el_BQ = new TH1D("dke_reco_el_BQ","",n_ke, ke_min, ke_max); dke_reco_el_BQ->Sumw2();
	kedept_reco_el_BQ = new TH1D("kedept_reco_el_BQ","",n_ke, ke_min, ke_max); kedept_reco_el_BQ->Sumw2();
	ke_reco_midcosmic_BQ = new TH1D("ke_reco_midcosmic_BQ","",n_ke, ke_min, ke_max); ke_reco_midcosmic_BQ->Sumw2();
	ke_reco_midpi_BQ = new TH1D("ke_reco_midpi_BQ","",n_ke, ke_min, ke_max); ke_reco_midpi_BQ->Sumw2();
	ke_reco_midp_BQ = new TH1D("ke_reco_midp_BQ","",n_ke, ke_min, ke_max); ke_reco_midp_BQ->Sumw2();
	ke_reco_midmu_BQ = new TH1D("ke_reco_midmu_BQ","",n_ke, ke_min, ke_max); ke_reco_midmu_BQ->Sumw2();
	ke_reco_mideg_BQ = new TH1D("ke_reco_mideg_BQ","",n_ke, ke_min, ke_max); ke_reco_mideg_BQ->Sumw2();
	ke_reco_midother_BQ = new TH1D("ke_reco_midother_BQ","",n_ke, ke_min, ke_max); ke_reco_midother_BQ->Sumw2();

	keff_reco_BQ = new TH1D("keff_reco_BQ","",n_ke, ke_min, ke_max); keff_reco_BQ->Sumw2();
	keff_reco_inel_BQ = new TH1D("keff_reco_inel_BQ","",n_ke, ke_min, ke_max); keff_reco_inel_BQ->Sumw2();
	keff_reco_el_BQ = new TH1D("keff_reco_el_BQ","",n_ke, ke_min, ke_max); keff_reco_el_BQ->Sumw2();
	keff_reco_midcosmic_BQ = new TH1D("keff_reco_midcosmic_BQ","",n_ke, ke_min, ke_max); keff_reco_midcosmic_BQ->Sumw2();
	keff_reco_midpi_BQ = new TH1D("keff_reco_midpi_BQ","",n_ke, ke_min, ke_max); keff_reco_midpi_BQ->Sumw2();
	keff_reco_midp_BQ = new TH1D("keff_reco_midp_BQ","",n_ke, ke_min, ke_max); keff_reco_midp_BQ->Sumw2();
	keff_reco_midmu_BQ = new TH1D("keff_reco_midmu_BQ","",n_ke, ke_min, ke_max); keff_reco_midmu_BQ->Sumw2();
	keff_reco_mideg_BQ = new TH1D("keff_reco_mideg_BQ","",n_ke, ke_min, ke_max); keff_reco_mideg_BQ->Sumw2();
	keff_reco_midother_BQ = new TH1D("keff_reco_midother_BQ","",n_ke, ke_min, ke_max); keff_reco_midother_BQ->Sumw2();

	//reco inel cut
	ke_reco_RecoInel = new TH1D("ke_reco_RecoInel","",n_ke, ke_min, ke_max); ke_reco_RecoInel->Sumw2();
	ke_reco_inel_RecoInel = new TH1D("ke_reco_inel_RecoInel","",n_ke, ke_min, ke_max); ke_reco_inel_RecoInel->Sumw2();
	ke_reco_el_RecoInel = new TH1D("ke_reco_el_RecoInel","",n_ke, ke_min, ke_max); ke_reco_el_RecoInel->Sumw2();
	ke_reco_midcosmic_RecoInel = new TH1D("ke_reco_midcosmic_RecoInel","",n_ke, ke_min, ke_max); ke_reco_midcosmic_RecoInel->Sumw2();
	ke_reco_midpi_RecoInel = new TH1D("ke_reco_midpi_RecoInel","",n_ke, ke_min, ke_max); ke_reco_midpi_RecoInel->Sumw2();
	ke_reco_midp_RecoInel = new TH1D("ke_reco_midp_RecoInel","",n_ke, ke_min, ke_max); ke_reco_midp_RecoInel->Sumw2();
	ke_reco_midmu_RecoInel = new TH1D("ke_reco_midmu_RecoInel","",n_ke, ke_min, ke_max); ke_reco_midmu_RecoInel->Sumw2();
	ke_reco_mideg_RecoInel = new TH1D("ke_reco_mideg_RecoInel","",n_ke, ke_min, ke_max); ke_reco_mideg_RecoInel->Sumw2();
	ke_reco_midother_RecoInel = new TH1D("ke_reco_midother_RecoInel","",n_ke, ke_min, ke_max); ke_reco_midother_RecoInel->Sumw2();

	//reco el cut
	ke_reco_RecoEl = new TH1D("ke_reco_RecoEl","",n_ke, ke_min, ke_max); ke_reco_RecoEl->Sumw2();
	ke_reco_inel_RecoEl = new TH1D("ke_reco_inel_RecoEl","",n_ke, ke_min, ke_max); ke_reco_inel_RecoEl->Sumw2();
	ke_reco_el_RecoEl = new TH1D("ke_reco_el_RecoEl","",n_ke, ke_min, ke_max); ke_reco_el_RecoEl->Sumw2();
	ke_reco_midcosmic_RecoEl = new TH1D("ke_reco_midcosmic_RecoEl","",n_ke, ke_min, ke_max); ke_reco_midcosmic_RecoEl->Sumw2();
	ke_reco_midpi_RecoEl = new TH1D("ke_reco_midpi_RecoEl","",n_ke, ke_min, ke_max); ke_reco_midpi_RecoEl->Sumw2();
	ke_reco_midp_RecoEl = new TH1D("ke_reco_midp_RecoEl","",n_ke, ke_min, ke_max); ke_reco_midp_RecoEl->Sumw2();
	ke_reco_midmu_RecoEl = new TH1D("ke_reco_midmu_RecoEl","",n_ke, ke_min, ke_max); ke_reco_midmu_RecoEl->Sumw2();
	ke_reco_mideg_RecoEl = new TH1D("ke_reco_mideg_RecoEl","",n_ke, ke_min, ke_max); ke_reco_mideg_RecoEl->Sumw2();
	ke_reco_midother_RecoEl = new TH1D("ke_reco_midother_RecoEl","",n_ke, ke_min, ke_max); ke_reco_midother_RecoEl->Sumw2();

	//Midp-rich cut
	ke_reco_MidP = new TH1D("ke_reco_MidP","",n_ke, ke_min, ke_max); ke_reco_MidP->Sumw2();
	ke_reco_inel_MidP = new TH1D("ke_reco_inel_MidP","",n_ke, ke_min, ke_max); ke_reco_inel_MidP->Sumw2();
	ke_reco_el_MidP = new TH1D("ke_reco_el_MidP","",n_ke, ke_min, ke_max); ke_reco_el_MidP->Sumw2();
	ke_reco_midcosmic_MidP = new TH1D("ke_reco_midcosmic_MidP","",n_ke, ke_min, ke_max); ke_reco_midcosmic_MidP->Sumw2();
	ke_reco_midpi_MidP = new TH1D("ke_reco_midpi_MidP","",n_ke, ke_min, ke_max); ke_reco_midpi_MidP->Sumw2();
	ke_reco_midp_MidP = new TH1D("ke_reco_midp_MidP","",n_ke, ke_min, ke_max); ke_reco_midp_MidP->Sumw2();
	ke_reco_midmu_MidP = new TH1D("ke_reco_midmu_MidP","",n_ke, ke_min, ke_max); ke_reco_midmu_MidP->Sumw2();
	ke_reco_mideg_MidP = new TH1D("ke_reco_mideg_MidP","",n_ke, ke_min, ke_max); ke_reco_mideg_MidP->Sumw2();
	ke_reco_midother_MidP = new TH1D("ke_reco_midother_MidP","",n_ke, ke_min, ke_max); ke_reco_midother_MidP->Sumw2();


	//dke -----------------------------------------------------------------------------------------------------//
	int n_dke=160;
	float dke_st=-800;
	float dke_end=800; 

	//nocut
	dke_inel_NoCut = new TH1D("dke_inel_NoCut","",n_dke, dke_st, dke_end); dke_inel_NoCut->Sumw2();
	dke_el_NoCut = new TH1D("dke_el_NoCut","",n_dke, dke_st, dke_end); dke_el_NoCut->Sumw2();
	dke_midcosmic_NoCut = new TH1D("dke_midcosmic_NoCut","",n_dke, dke_st, dke_end); dke_midcosmic_NoCut->Sumw2();
	dke_midpi_NoCut = new TH1D("dke_midpi_NoCut","",n_dke, dke_st, dke_end); dke_midpi_NoCut->Sumw2();
	dke_midp_NoCut = new TH1D("dke_midp_NoCut","",n_dke, dke_st, dke_end); dke_midp_NoCut->Sumw2();
	dke_midmu_NoCut = new TH1D("dke_midmu_NoCut","",n_dke, dke_st, dke_end); dke_midmu_NoCut->Sumw2();
	dke_mideg_NoCut = new TH1D("dke_mideg_NoCut","",n_dke, dke_st, dke_end); dke_mideg_NoCut->Sumw2();
	dke_midother_NoCut = new TH1D("dke_midother_NoCut","",n_dke, dke_st, dke_end); dke_midother_NoCut->Sumw2();

	//pandora cut
	dke_inel_PanS = new TH1D("dke_inel_PanS","",n_dke, dke_st, dke_end); dke_inel_PanS->Sumw2();
	dke_el_PanS = new TH1D("dke_el_PanS","",n_dke, dke_st, dke_end); dke_el_PanS->Sumw2();
	dke_midcosmic_PanS = new TH1D("dke_midcosmic_PanS","",n_dke, dke_st, dke_end); dke_midcosmic_PanS->Sumw2();
	dke_midpi_PanS = new TH1D("dke_midpi_PanS","",n_dke, dke_st, dke_end); dke_midpi_PanS->Sumw2();
	dke_midp_PanS = new TH1D("dke_midp_PanS","",n_dke, dke_st, dke_end); dke_midp_PanS->Sumw2();
	dke_midmu_PanS = new TH1D("dke_midmu_PanS","",n_dke, dke_st, dke_end); dke_midmu_PanS->Sumw2();
	dke_mideg_PanS = new TH1D("dke_mideg_PanS","",n_dke, dke_st, dke_end); dke_mideg_PanS->Sumw2();
	dke_midother_PanS = new TH1D("dke_midother_PanS","",n_dke, dke_st, dke_end); dke_midother_PanS->Sumw2();
	
	//calosz
	dke_inel_CaloSz = new TH1D("dke_inel_CaloSz","",n_dke, dke_st, dke_end); dke_inel_CaloSz->Sumw2();
	dke_el_CaloSz = new TH1D("dke_el_CaloSz","",n_dke, dke_st, dke_end); dke_el_CaloSz->Sumw2();
	dke_midcosmic_CaloSz = new TH1D("dke_midcosmic_CaloSz","",n_dke, dke_st, dke_end); dke_midcosmic_CaloSz->Sumw2();
	dke_midpi_CaloSz = new TH1D("dke_midpi_CaloSz","",n_dke, dke_st, dke_end); dke_midpi_CaloSz->Sumw2();
	dke_midp_CaloSz = new TH1D("dke_midp_CaloSz","",n_dke, dke_st, dke_end); dke_midp_CaloSz->Sumw2();
	dke_midmu_CaloSz = new TH1D("dke_midmu_CaloSz","",n_dke, dke_st, dke_end); dke_midmu_CaloSz->Sumw2();
	dke_mideg_CaloSz = new TH1D("dke_mideg_CaloSz","",n_dke, dke_st, dke_end); dke_mideg_CaloSz->Sumw2();
	dke_midother_CaloSz = new TH1D("dke_midother_CaloSz","",n_dke, dke_st, dke_end); dke_midother_CaloSz->Sumw2();

	//beam quality
	dke_inel_BQ = new TH1D("dke_inel_BQ","",n_dke, dke_st, dke_end); dke_inel_BQ->Sumw2();
	dke_el_BQ = new TH1D("dke_el_BQ","",n_dke, dke_st, dke_end); dke_el_BQ->Sumw2();
	dke_midcosmic_BQ = new TH1D("dke_midcosmic_BQ","",n_dke, dke_st, dke_end); dke_midcosmic_BQ->Sumw2();
	dke_midpi_BQ = new TH1D("dke_midpi_BQ","",n_dke, dke_st, dke_end); dke_midpi_BQ->Sumw2();
	dke_midp_BQ = new TH1D("dke_midp_BQ","",n_dke, dke_st, dke_end); dke_midp_BQ->Sumw2();
	dke_midmu_BQ = new TH1D("dke_midmu_BQ","",n_dke, dke_st, dke_end); dke_midmu_BQ->Sumw2();
	dke_mideg_BQ = new TH1D("dke_mideg_BQ","",n_dke, dke_st, dke_end); dke_mideg_BQ->Sumw2();
	dke_midother_BQ = new TH1D("dke_midother_BQ","",n_dke, dke_st, dke_end); dke_midother_BQ->Sumw2();

	//reco inel
	dke_inel_RecoInel = new TH1D("dke_inel_RecoInel","",n_dke, dke_st, dke_end); dke_inel_RecoInel->Sumw2();
	dke_el_RecoInel = new TH1D("dke_el_RecoInel","",n_dke, dke_st, dke_end); dke_el_RecoInel->Sumw2();
	dke_midcosmic_RecoInel = new TH1D("dke_midcosmic_RecoInel","",n_dke, dke_st, dke_end); dke_midcosmic_RecoInel->Sumw2();
	dke_midpi_RecoInel = new TH1D("dke_midpi_RecoInel","",n_dke, dke_st, dke_end); dke_midpi_RecoInel->Sumw2();
	dke_midp_RecoInel = new TH1D("dke_midp_RecoInel","",n_dke, dke_st, dke_end); dke_midp_RecoInel->Sumw2();
	dke_midmu_RecoInel = new TH1D("dke_midmu_RecoInel","",n_dke, dke_st, dke_end); dke_midmu_RecoInel->Sumw2();
	dke_mideg_RecoInel = new TH1D("dke_mideg_RecoInel","",n_dke, dke_st, dke_end); dke_mideg_RecoInel->Sumw2();
	dke_midother_RecoInel = new TH1D("dke_midother_RecoInel","",n_dke, dke_st, dke_end); dke_midother_RecoInel->Sumw2();

	//reco el
	dke_inel_RecoEl = new TH1D("dke_inel_RecoEl","",n_dke, dke_st, dke_end); dke_inel_RecoEl->Sumw2();
	dke_el_RecoEl = new TH1D("dke_el_RecoEl","",n_dke, dke_st, dke_end); dke_el_RecoEl->Sumw2();
	dke_midcosmic_RecoEl = new TH1D("dke_midcosmic_RecoEl","",n_dke, dke_st, dke_end); dke_midcosmic_RecoEl->Sumw2();
	dke_midpi_RecoEl = new TH1D("dke_midpi_RecoEl","",n_dke, dke_st, dke_end); dke_midpi_RecoEl->Sumw2();
	dke_midp_RecoEl = new TH1D("dke_midp_RecoEl","",n_dke, dke_st, dke_end); dke_midp_RecoEl->Sumw2();
	dke_midmu_RecoEl = new TH1D("dke_midmu_RecoEl","",n_dke, dke_st, dke_end); dke_midmu_RecoEl->Sumw2();
	dke_mideg_RecoEl = new TH1D("dke_mideg_RecoEl","",n_dke, dke_st, dke_end); dke_mideg_RecoEl->Sumw2();
	dke_midother_RecoEl = new TH1D("dke_midother_RecoEl","",n_dke, dke_st, dke_end); dke_midother_RecoEl->Sumw2();

	//ntrklen ------------------------------------------------------------------------------------//
	int n_ntrklen=61;
	double st_ntrklen=-0.02;
	double ed_ntrklen=1.2;
	//bq cut
	ntrklen_BQ=new TH1D("ntrklen_BQ", "", n_ntrklen, st_ntrklen, ed_ntrklen);
	ntrklen_inel_BQ=new TH1D("ntrklen_inel_BQ", "", n_ntrklen, st_ntrklen, ed_ntrklen);
	ntrklen_el_BQ=new TH1D("ntrklen_el_BQ", "", n_ntrklen, st_ntrklen, ed_ntrklen);
	ntrklen_midcosmic_BQ=new TH1D("ntrklen_midcosmic_BQ", "", n_ntrklen, st_ntrklen, ed_ntrklen);
	ntrklen_midpi_BQ=new TH1D("ntrklen_midpi_BQ", "", n_ntrklen, st_ntrklen, ed_ntrklen);
	ntrklen_midp_BQ=new TH1D("ntrklen_midp_BQ", "", n_ntrklen, st_ntrklen, ed_ntrklen);
	ntrklen_midmu_BQ=new TH1D("ntrklen_midmu_BQ", "", n_ntrklen, st_ntrklen, ed_ntrklen);
	ntrklen_mideg_BQ=new TH1D("ntrklen_mideg_BQ", "", n_ntrklen, st_ntrklen, ed_ntrklen);
	ntrklen_midother_BQ=new TH1D("ntrklen_midother_BQ", "", n_ntrklen, st_ntrklen, ed_ntrklen);

	//range vs ke
        //int dke=325;
        //float ke_st=-50;
        //float ke_end=600;
	trklen_ke_true_inel=new TH2D("trklen_ke_true_inel","",120,0,120, n_ke, ke_min, ke_max);
	trklen_ke_true_el=new TH2D("trklen_ke_true_el","",120,0,120, n_ke, ke_min, ke_max);

	//truth inc/int
	//h_true_incidents=new TH1D("h_true_incidents","true_incidents", nthinslices, 0, nthinslices-1);
	//h_true_st_incidents=new TH1D("h_true_st_incidents","true_st_incidents", nthinslices, 0, nthinslices-1);
	//h_true_interactions=new TH1D("h_true_interactions","true_interactions", nthinslices, 0, nthinslices-1);
	h_true_incidents=new TH1D("h_true_incidents","true_incidents", nthinslices+2, -1, nthinslices+1);
	h_true_st_incidents=new TH1D("h_true_st_incidents","true_st_incidents", nthinslices+2, -1, nthinslices+1);
	h_true_interactions=new TH1D("h_true_interactions","true_interactions", nthinslices+2, -1, nthinslices+1);
	h_true_incidents->Sumw2();
	h_true_st_incidents->Sumw2();
	h_true_interactions->Sumw2();


	//reco E-dept
	reco_dedx_trklen_inel=new TH2D("reco_dedx_trklen_inel","reco_dedx_trklen_inel",150,0,150, 100,0,100);
	reco_de_trklen_inel=new TH2D("reco_de_trklen_inel","reco_de_trklen_inel",150,0,150, 100,0,100);
	reco_dx_trklen_inel=new TH2D("reco_dx_trklen_inel","reco_dx_trklen_inel",150,0,150,100,0,10);
	reco_dedx_trklen_inel->Sumw2();
	reco_de_trklen_inel->Sumw2();
	reco_dx_trklen_inel->Sumw2();

	//KE(Bethe-Bloch) vs trklen ------------------------------------------------------//
	//KEbb_truetrklen_all=new TH2D("KEbb_truetrklen_all","", 1210,-1,120,800,0,800);
	//KEbb_truetrklen_inel=new TH2D("KEbb_truetrklen_inel","", 1210,-1,120,800,0,800);
	//KEbb_recotrklen_all=new TH2D("KEbb_recotrklen_all","", 1210,-1,120,800,0,800);
	//KEbb_recotrklen_inel=new TH2D("KEbb_recotrklen_inel","", 1210,-1,120,800,0,800);
	//KEbb_truetrklen_all->Sumw2(); //
	//KEbb_truetrklen_inel->Sumw2(); //
	//KEbb_recotrklen_all->Sumw2(); //
	//KEbb_recotrklen_inel->Sumw2(); //


	KEcalo_truetrklen_all=new TH2D("KEcalo_truetrklen_all","", 1210,-1,120,nke,kemin,kemax);
	KEcalo_truetrklen_inel=new TH2D("KEcalo_truetrklen_inel","", 1210,-1,120,nke,kemin,kemax);
	KEcalo_recotrklen_all=new TH2D("KEcalo_recotrklen_all","", 1210,-1,120,nke,kemin,kemax);
	KEcalo_recotrklen_inel=new TH2D("KEcalo_recotrklen_inel","", 1210,-1,120,nke,kemin,kemax);
	KEcalo_truetrklen_all->Sumw2(); //
	KEcalo_truetrklen_inel->Sumw2(); //
	KEcalo_recotrklen_all->Sumw2(); //
	KEcalo_recotrklen_inel->Sumw2(); //


	true_trklen_patch_all=new TH1D("true_trklen_patch_all","",500,0,5);
	true_trklen_patch_all->Sumw2();	

	//R vs KE1st
	//h2d_R_kE1st=new TH2D("h2d_R_kE1st","", 800, 0, 800, 200, 0, 2);
	//h2d_R_kE1st->Sumw2();

} //BookHistograms


void SaveHistograms() { //SaveHistograms
	outputFile->cd();
	outputFile->Write();
	//h_truesliceid_uf->Write("h_truesliceid_uf");
	//h_truesliceid_inelastic_uf->Write("h_truesliceid_inelastic_uf");

} //SaveHistograms

void CalcXS(const Unfold & uf) { //CalcXS

	double slcid[nthinslices] = {0};
	double avg_trueincE[nthinslices] = {0};
	double avg_recoincE[nthinslices] = {0};
	double err_trueincE[nthinslices] = {0};
	double err_recoincE[nthinslices] = {0};
	double reco_trueincE[nthinslices] = {0};
	double err_reco_trueincE[nthinslices] = {0};
	double truexs[nthinslices] = {0};
	double err_truexs[nthinslices] = {0};

	double NA=6.02214076e23;
	double MAr=39.95; //gmol
	double Density = 1.4; // g/cm^3
  	double true_cosangle = 1.;

	for (int i = 0; i<nthinslices; ++i){
		slcid[i] = i;
		avg_trueincE[i] = true_incE[i]->GetMean();
		err_trueincE[i] = true_incE[i]->GetMeanError();
		avg_recoincE[i] = reco_incE[i]->GetMean();
		err_recoincE[i] = reco_incE[i]->GetMeanError();
		reco_trueincE[i] = avg_recoincE[i] - avg_trueincE[i];
		err_reco_trueincE[i] = sqrt(pow(err_trueincE[i],2)+pow(err_recoincE[i],2));
		//std::cout<<i<<" "<<avg_trueincE[i]<<std::endl;
		if (true_incidents[i] && true_interactions[i]){
			//truexs[i] = MAr/(Density*NA*thinslicewidth/true_AngCorr->GetMean())*log(true_incidents[i]/(true_incidents[i]-true_interactions[i]))*1e27;
			//err_truexs[i] = MAr/(Density*NA*thinslicewidth/true_AngCorr->GetMean())*1e27*sqrt(true_interactions[i]+pow(true_interactions[i],2)/true_incidents[i])/true_incidents[i];
      			truexs[i] = MAr/(Density*NA*thinslicewidth/true_cosangle)*log(true_incidents[i]/(true_incidents[i]-true_interactions[i]))*1e27;
      			err_truexs[i] = MAr/(Density*NA*thinslicewidth/true_cosangle)*1e27*sqrt(true_interactions[i]+pow(true_interactions[i],2)/true_incidents[i])/true_incidents[i];
		}
	}

	TGraphErrors *gr_trueincE = new TGraphErrors(nthinslices, &(slcid[0]), &(avg_trueincE[0]), 0, &(err_trueincE[0]));
	TGraphErrors *gr_recoincE = new TGraphErrors(nthinslices, &(slcid[0]), &(avg_recoincE[0]), 0, &(err_recoincE[0]));
	TGraphErrors *gr_reco_trueincE = new TGraphErrors(nthinslices, &(slcid[0]), &(reco_trueincE[0]), 0, &(err_reco_trueincE[0]));

	gr_trueincE->Write("gr_trueincE");
	gr_recoincE->Write("gr_recoincE");
	gr_reco_trueincE->Write("gr_reco_trueincE");

	TGraphErrors *gr_truexs = new TGraphErrors(nthinslices, &(avg_trueincE[0]), &(truexs[0]), 0, &(err_truexs[0]));

	gr_truexs->Write("gr_truexs");

	TH1D *hinc = (TH1D*)h_recosliceid_allevts_cuts->Clone("hinc");
	//TH1D *hint = (TH1D*)h_recosliceid_allevts_cuts->Clone("hint");
	TH1D *hint = (TH1D*)h_recosliceid_recoinelastic_cuts->Clone("hint");
	hinc->Multiply(uf.pur_Inc);
	hint->Multiply(uf.pur_Int);

	//  RooUnfoldBayes   unfold_Inc (&uf.response_SliceID_Inc, uf.pur_num_Inc, 4);
	//  RooUnfoldBayes   unfold_Int (&uf.response_SliceID_Int, uf.pur_num_Int, 4);

	RooUnfoldBayes   unfold_Inc (&uf.response_SliceID_Inc, hinc, 4);
	RooUnfoldBayes   unfold_Int (&uf.response_SliceID_Int, hint, 4);

	//RooUnfoldBayes   unfold_Inc (&uf.response_SliceID_Inc, hinc, 12);
	//RooUnfoldBayes   unfold_Int (&uf.response_SliceID_Int, hint, 12);

	//RooUnfoldSvd     unfold_Inc (&uf.response_SliceID_Inc, hinc, 20);   // OR
	//RooUnfoldSvd     unfold_Int (&uf.response_SliceID_Int, hint, 20);   // OR

	//  RooUnfoldSvd     unfold_Inc (&uf.response_SliceID_Inc, uf.pur_num_Inc, 20);   // OR
	//  RooUnfoldSvd     unfold_Int (&uf.response_SliceID_Int, uf.pur_num_Int, 20);   // OR

	//h_truesliceid_uf = (TH1D*) unfold_Inc.Hreco();
	//h_truesliceid_inelastic_uf = (TH1D*) unfold_Int.Hreco();

} //CalcXS

