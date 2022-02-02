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
//TH1D *h_truesliceid_uf;
TH1D *h_truesliceid_cuts;
TH1D *h_truesliceid_inelastic_all;
//TH1D *h_truesliceid_inelastic_uf;
TH1D *h_truesliceid_inelastic_cuts;
TH1D *h_recosliceid_allevts_cuts;
TH1D *h_recosliceid_allevts_cuts_inel;
TH1D *h_recosliceid_allevts_cuts_el;
TH1D *h_recosliceid_allevts_cuts_midcosmic;
TH1D *h_recosliceid_allevts_cuts_midpi;
TH1D *h_recosliceid_allevts_cuts_midp;
TH1D *h_recosliceid_allevts_cuts_midmu;
TH1D *h_recosliceid_allevts_cuts_mideg;
TH1D *h_recosliceid_allevts_cuts_midother;

TH1D *h_recosliceid_cuts;
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

double true_interactions[nthinslices];
double true_incidents[nthinslices];

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
TH1D *KEtrue_Beam_inel;
TH1D *KEtrue_ff_inel;
TH2D *KEtrue_z_inel;
TH2D *KEtrue_range_inel;
TH2D *true_z_range_inel;
TH3D *true_z_range_KEtrue_inel;


//Track length histograms -----------------------------------------------------------------//
//true range
//no cut
TH1D *trklen_true_inel_NoCut;
TH1D *trklen_true_el_NoCut;
TH1D *trklen_true_midcosmic_NoCut;
TH1D *trklen_true_midpi_NoCut;
TH1D *trklen_true_midp_NoCut;
TH1D *trklen_true_midmu_NoCut;
TH1D *trklen_true_mideg_NoCut;
TH1D *trklen_true_midother_NoCut;

//pandora slice cut
TH1D *trklen_true_inel_PanS;
TH1D *trklen_true_el_PanS;
TH1D *trklen_true_midcosmic_PanS;
TH1D *trklen_true_midpi_PanS;
TH1D *trklen_true_midp_PanS;
TH1D *trklen_true_midmu_PanS;
TH1D *trklen_true_mideg_PanS;
TH1D *trklen_true_midother_PanS;

//calosz cut
TH1D *trklen_true_inel_CaloSz;
TH1D *trklen_true_el_CaloSz;
TH1D *trklen_true_midcosmic_CaloSz;
TH1D *trklen_true_midpi_CaloSz;
TH1D *trklen_true_midp_CaloSz;
TH1D *trklen_true_midmu_CaloSz;
TH1D *trklen_true_mideg_CaloSz;
TH1D *trklen_true_midother_CaloSz;

//beam quality cut
TH1D *trklen_true_inel_BQ;
TH1D *trklen_true_el_BQ;
TH1D *trklen_true_midcosmic_BQ;
TH1D *trklen_true_midpi_BQ;
TH1D *trklen_true_midp_BQ;
TH1D *trklen_true_midmu_BQ;
TH1D *trklen_true_mideg_BQ;
TH1D *trklen_true_midother_BQ;

//RecoInel cut
TH1D *trklen_true_inel_RecoInel;
TH1D *trklen_true_el_RecoInel;
TH1D *trklen_true_midcosmic_RecoInel;
TH1D *trklen_true_midpi_RecoInel;
TH1D *trklen_true_midp_RecoInel;
TH1D *trklen_true_midmu_RecoInel;
TH1D *trklen_true_mideg_RecoInel;
TH1D *trklen_true_midother_RecoInel;

//reco range
//no cut
TH1D *trklen_reco_inel_NoCut;
TH1D *trklen_reco_el_NoCut;
TH1D *trklen_reco_midcosmic_NoCut;
TH1D *trklen_reco_midpi_NoCut;
TH1D *trklen_reco_midp_NoCut;
TH1D *trklen_reco_midmu_NoCut;
TH1D *trklen_reco_mideg_NoCut;
TH1D *trklen_reco_midother_NoCut;

//pandora slice cut
TH1D *trklen_reco_inel_PanS;
TH1D *trklen_reco_el_PanS;
TH1D *trklen_reco_midcosmic_PanS;
TH1D *trklen_reco_midpi_PanS;
TH1D *trklen_reco_midp_PanS;
TH1D *trklen_reco_midmu_PanS;
TH1D *trklen_reco_mideg_PanS;
TH1D *trklen_reco_midother_PanS;

//calosz cut
TH1D *trklen_reco_inel_CaloSz;
TH1D *trklen_reco_el_CaloSz;
TH1D *trklen_reco_midcosmic_CaloSz;
TH1D *trklen_reco_midpi_CaloSz;
TH1D *trklen_reco_midp_CaloSz;
TH1D *trklen_reco_midmu_CaloSz;
TH1D *trklen_reco_mideg_CaloSz;
TH1D *trklen_reco_midother_CaloSz;

//beam quality cut
TH1D *trklen_reco_inel_BQ;
TH1D *trklen_reco_el_BQ;
TH1D *trklen_reco_midcosmic_BQ;
TH1D *trklen_reco_midpi_BQ;
TH1D *trklen_reco_midp_BQ;
TH1D *trklen_reco_midmu_BQ;
TH1D *trklen_reco_mideg_BQ;
TH1D *trklen_reco_midother_BQ;

//RecoInel cut
TH1D *trklen_reco_inel_RecoInel;
TH1D *trklen_reco_el_RecoInel;
TH1D *trklen_reco_midcosmic_RecoInel;
TH1D *trklen_reco_midpi_RecoInel;
TH1D *trklen_reco_midp_RecoInel;
TH1D *trklen_reco_midmu_RecoInel;
TH1D *trklen_reco_mideg_RecoInel;
TH1D *trklen_reco_midother_RecoInel;

//(reco-truth)trklen
//no cut
TH1D *dtrklen;
TH1D *dtrklen_inel_NoCut;
TH1D *dtrklen_el_NoCut;
TH1D *dtrklen_midcosmic_NoCut;
TH1D *dtrklen_midpi_NoCut;
TH1D *dtrklen_midp_NoCut;
TH1D *dtrklen_midmu_NoCut;
TH1D *dtrklen_mideg_NoCut;
TH1D *dtrklen_midother_NoCut;

//pandora slice cut
TH1D *dtrklen_inel_PanS;
TH1D *dtrklen_el_PanS;
TH1D *dtrklen_midcosmic_PanS;
TH1D *dtrklen_midpi_PanS;
TH1D *dtrklen_midp_PanS;
TH1D *dtrklen_midmu_PanS;
TH1D *dtrklen_mideg_PanS;
TH1D *dtrklen_midother_PanS;

//calosz cut
TH1D *dtrklen_inel_CaloSz;
TH1D *dtrklen_el_CaloSz;
TH1D *dtrklen_midcosmic_CaloSz;
TH1D *dtrklen_midpi_CaloSz;
TH1D *dtrklen_midp_CaloSz;
TH1D *dtrklen_midmu_CaloSz;
TH1D *dtrklen_mideg_CaloSz;
TH1D *dtrklen_midother_CaloSz;

//beam quality cut
TH1D *dtrklen_inel_BQ;
TH1D *dtrklen_el_BQ;
TH1D *dtrklen_midcosmic_BQ;
TH1D *dtrklen_midpi_BQ;
TH1D *dtrklen_midp_BQ;
TH1D *dtrklen_midmu_BQ;
TH1D *dtrklen_mideg_BQ;
TH1D *dtrklen_midother_BQ;

//RecoInel cut
TH1D *dtrklen_inel_RecoInel;
TH1D *dtrklen_el_RecoInel;
TH1D *dtrklen_midcosmic_RecoInel;
TH1D *dtrklen_midpi_RecoInel;
TH1D *dtrklen_midp_RecoInel;
TH1D *dtrklen_midmu_RecoInel;
TH1D *dtrklen_mideg_RecoInel;
TH1D *dtrklen_midother_RecoInel;

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
TH1D *h_true_interactions;

//Reco E-dept [Inel]
TH2D *reco_dedx_trklen_inel;
TH2D *reco_de_trklen_inel;
TH2D *reco_dx_trklen_inel;



void BookHistograms() { //BookHistograms
	outputFile = TFile::Open(fOutputFileName.c_str(), "recreate"); //open output file

	//XS histograms -------------------------------------------------------------------------------------------------------------------------------------------------------//
	for (int i = 0; i<nthinslices; ++i){
		reco_incE[i] = new TH1D(Form("reco_incE_%d",i),Form("Reco incident energy, %.1f < z < %.1f (cm)",i*thinslicewidth, (i+1)*thinslicewidth), nbinse, 0, 1200.);
		true_incE[i] = new TH1D(Form("true_incE_%d",i),Form("True incident energy, %.1f < z < %.1f (cm)",i*thinslicewidth, (i+1)*thinslicewidth), nbinse, 0, 1200.);
		reco_incE[i]->Sumw2();
		true_incE[i]->Sumw2();
	}

	reco_AngCorr = new TH1D("reco_AngCorr","Reco angle correction", 100, 0, 1.);
	true_AngCorr = new TH1D("true_AngCorr","true angle correction", 100, 0, 1.);
	reco_AngCorr->Sumw2();
	true_AngCorr->Sumw2();

	h_truesliceid_all = new TH1D("h_truesliceid_all","h_truesliceid_all;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_cuts = new TH1D("h_truesliceid_cuts","h_truesliceid_cuts;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_inelastic_all = new TH1D("h_truesliceid_inelastic_all","h_truesliceid_inelastic_all;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_inelastic_cuts = new TH1D("h_truesliceid_inelastic_cuts","h_truesliceid_inelastic_cuts;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts = new TH1D("h_recosliceid_allevts_cuts","h_recosliceid_allevts_cuts;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts_inel = new TH1D("h_recosliceid_allevts_cuts_inel","h_recosliceid_allevts_cuts_inel;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts_el = new TH1D("h_recosliceid_allevts_cuts_el","h_recosliceid_allevts_cuts_el;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts_midcosmic = new TH1D("h_recosliceid_allevts_cuts_midcosmic","h_recosliceid_allevts_cuts_midcosmic;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts_midpi = new TH1D("h_recosliceid_allevts_cuts_midpi","h_recosliceid_allevts_cuts_midpi;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts_midp = new TH1D("h_recosliceid_allevts_cuts_midp","h_recosliceid_allevts_cuts_midp;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts_midmu = new TH1D("h_recosliceid_allevts_cuts_midmu","h_recosliceid_allevts_cuts_midmu;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts_mideg = new TH1D("h_recosliceid_allevts_cuts_mideg","h_recosliceid_allevts_cuts_mideg;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_allevts_cuts_midother = new TH1D("h_recosliceid_allevts_cuts_midother","h_recosliceid_allevts_cuts_midother;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);


	h_recosliceid_cuts = new TH1D("h_recosliceid_cuts","h_recosliceid_cuts;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
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

	h_truesliceid_all->Sumw2();
	h_truesliceid_cuts->Sumw2();
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


	h_recosliceid_cuts->Sumw2();
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

	for (int i = 0; i<nthinslices; ++i){
		true_interactions[i] = 0;
		true_incidents[i] = 0;
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
	int n_ke=500;
	float ke_min=0;
	float ke_max=1000;
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


	//KEtruth
	KEtrue_Beam_inel=new TH1D("KEtrue_Beam_inel","",n_ke, ke_min, ke_max); KEtrue_Beam_inel->Sumw2();
	KEtrue_ff_inel=new TH1D("KEtrue_ff_inel", "", n_ke, ke_min, ke_max); KEtrue_ff_inel->Sumw2();
	KEtrue_z_inel=new TH2D("KEtrue_z_inel","",300,0,150, 500,0,1000);
	KEtrue_range_inel=new TH2D("KEtrue_range_inel","",300,0,150, 500,0,1000);
	true_z_range_inel=new TH2D("true_z_range_inel","",300,0,150, 300,0,150);
	true_z_range_KEtrue_inel=new TH3D("true_z_range_KEtrue_inel","",300,0,150, 300,0,150, 500,0,1000);
	true_z_range_KEtrue_inel->GetXaxis()->SetTitle("z [cm]");
	true_z_range_KEtrue_inel->GetYaxis()->SetTitle("range [cm]");
	true_z_range_KEtrue_inel->GetZaxis()->SetTitle("KE [MeV]");

	//trklen -----------------------------------------------------------------------------------------------------//
        int n_trklen=34;
        double trklen_min=-4;
        double trklen_max=132;

	//truth range
	//no cut
	trklen_true_inel_NoCut = new TH1D("trklen_true_inel_NoCut","",n_trklen, trklen_min, trklen_max); trklen_true_inel_NoCut->Sumw2();
	trklen_true_el_NoCut = new TH1D("trklen_true_el_NoCut","",n_trklen, trklen_min, trklen_max); trklen_true_el_NoCut->Sumw2();
	trklen_true_midcosmic_NoCut = new TH1D("trklen_true_midcosmic_NoCut","",n_trklen, trklen_min, trklen_max); trklen_true_midcosmic_NoCut->Sumw2();
	trklen_true_midpi_NoCut = new TH1D("trklen_true_midpi_NoCut","",n_trklen, trklen_min, trklen_max); trklen_true_midpi_NoCut->Sumw2();
	trklen_true_midp_NoCut = new TH1D("trklen_true_midp_NoCut","",n_trklen, trklen_min, trklen_max); trklen_true_midp_NoCut->Sumw2();
	trklen_true_midmu_NoCut = new TH1D("trklen_true_midmu_NoCut","",n_trklen, trklen_min, trklen_max); trklen_true_midmu_NoCut->Sumw2();
	trklen_true_mideg_NoCut = new TH1D("trklen_true_mideg_NoCut","",n_trklen, trklen_min, trklen_max); trklen_true_mideg_NoCut->Sumw2();
	trklen_true_midother_NoCut = new TH1D("trklen_true_midother_NoCut","",n_trklen, trklen_min, trklen_max); trklen_true_midother_NoCut->Sumw2();

	//pandora cut
	trklen_true_inel_PanS = new TH1D("trklen_true_inel_PanS","",n_trklen, trklen_min, trklen_max); trklen_true_inel_PanS->Sumw2();
	trklen_true_el_PanS = new TH1D("trklen_true_el_PanS","",n_trklen, trklen_min, trklen_max); trklen_true_el_PanS->Sumw2();
	trklen_true_midcosmic_PanS = new TH1D("trklen_true_midcosmic_PanS","",n_trklen, trklen_min, trklen_max); trklen_true_midcosmic_PanS->Sumw2();
	trklen_true_midpi_PanS = new TH1D("trklen_true_midpi_PanS","",n_trklen, trklen_min, trklen_max); trklen_true_midpi_PanS->Sumw2();
	trklen_true_midp_PanS = new TH1D("trklen_true_midp_PanS","",n_trklen, trklen_min, trklen_max); trklen_true_midp_PanS->Sumw2();
	trklen_true_midmu_PanS = new TH1D("trklen_true_midmu_PanS","",n_trklen, trklen_min, trklen_max); trklen_true_midmu_PanS->Sumw2();
	trklen_true_mideg_PanS = new TH1D("trklen_true_mideg_PanS","",n_trklen, trklen_min, trklen_max); trklen_true_mideg_PanS->Sumw2();
	trklen_true_midother_PanS = new TH1D("trklen_true_midother_PanS","",n_trklen, trklen_min, trklen_max); trklen_true_midother_PanS->Sumw2();

	//CaloSz
	trklen_true_inel_CaloSz = new TH1D("trklen_true_inel_CaloSz","",n_trklen, trklen_min, trklen_max); trklen_true_inel_CaloSz->Sumw2();
	trklen_true_el_CaloSz = new TH1D("trklen_true_el_CaloSz","",n_trklen, trklen_min, trklen_max); trklen_true_el_CaloSz->Sumw2();
	trklen_true_midcosmic_CaloSz = new TH1D("trklen_true_midcosmic_CaloSz","",n_trklen, trklen_min, trklen_max); trklen_true_midcosmic_CaloSz->Sumw2();
	trklen_true_midpi_CaloSz = new TH1D("trklen_true_midpi_CaloSz","",n_trklen, trklen_min, trklen_max); trklen_true_midpi_CaloSz->Sumw2();
	trklen_true_midp_CaloSz = new TH1D("trklen_true_midp_CaloSz","",n_trklen, trklen_min, trklen_max); trklen_true_midp_CaloSz->Sumw2();
	trklen_true_midmu_CaloSz = new TH1D("trklen_true_midmu_CaloSz","",n_trklen, trklen_min, trklen_max); trklen_true_midmu_CaloSz->Sumw2();
	trklen_true_mideg_CaloSz = new TH1D("trklen_true_mideg_CaloSz","",n_trklen, trklen_min, trklen_max); trklen_true_mideg_CaloSz->Sumw2();
	trklen_true_midother_CaloSz = new TH1D("trklen_true_midother_CaloSz","",n_trklen, trklen_min, trklen_max); trklen_true_midother_CaloSz->Sumw2();

	//beam quality
	trklen_true_inel_BQ = new TH1D("trklen_true_inel_BQ","",n_trklen, trklen_min, trklen_max); trklen_true_inel_BQ->Sumw2();
	trklen_true_el_BQ = new TH1D("trklen_true_el_BQ","",n_trklen, trklen_min, trklen_max); trklen_true_el_BQ->Sumw2();
	trklen_true_midcosmic_BQ = new TH1D("trklen_true_midcosmic_BQ","",n_trklen, trklen_min, trklen_max); trklen_true_midcosmic_BQ->Sumw2();
	trklen_true_midpi_BQ = new TH1D("trklen_true_midpi_BQ","",n_trklen, trklen_min, trklen_max); trklen_true_midpi_BQ->Sumw2();
	trklen_true_midp_BQ = new TH1D("trklen_true_midp_BQ","",n_trklen, trklen_min, trklen_max); trklen_true_midp_BQ->Sumw2();
	trklen_true_midmu_BQ = new TH1D("trklen_true_midmu_BQ","",n_trklen, trklen_min, trklen_max); trklen_true_midmu_BQ->Sumw2();
	trklen_true_mideg_BQ = new TH1D("trklen_true_mideg_BQ","",n_trklen, trklen_min, trklen_max); trklen_true_mideg_BQ->Sumw2();
	trklen_true_midother_BQ = new TH1D("trklen_true_midother_BQ","",n_trklen, trklen_min, trklen_max); trklen_true_midother_BQ->Sumw2();

	//reco inel cut
	trklen_true_inel_RecoInel = new TH1D("trklen_true_inel_RecoInel","",n_trklen, trklen_min, trklen_max); trklen_true_inel_RecoInel->Sumw2();
	trklen_true_el_RecoInel = new TH1D("trklen_true_el_RecoInel","",n_trklen, trklen_min, trklen_max); trklen_true_el_RecoInel->Sumw2();
	trklen_true_midcosmic_RecoInel = new TH1D("trklen_true_midcosmic_RecoInel","",n_trklen, trklen_min, trklen_max); trklen_true_midcosmic_RecoInel->Sumw2();
	trklen_true_midpi_RecoInel = new TH1D("trklen_true_midpi_RecoInel","",n_trklen, trklen_min, trklen_max); trklen_true_midpi_RecoInel->Sumw2();
	trklen_true_midp_RecoInel = new TH1D("trklen_true_midp_RecoInel","",n_trklen, trklen_min, trklen_max); trklen_true_midp_RecoInel->Sumw2();
	trklen_true_midmu_RecoInel = new TH1D("trklen_true_midmu_RecoInel","",n_trklen, trklen_min, trklen_max); trklen_true_midmu_RecoInel->Sumw2();
	trklen_true_mideg_RecoInel = new TH1D("trklen_true_mideg_RecoInel","",n_trklen, trklen_min, trklen_max); trklen_true_mideg_RecoInel->Sumw2();
	trklen_true_midother_RecoInel = new TH1D("trklen_true_midother_RecoInel","",n_trklen, trklen_min, trklen_max); trklen_true_midother_RecoInel->Sumw2();

	//reco range
	//no cut
	trklen_reco_inel_NoCut = new TH1D("trklen_reco_inel_NoCut","",n_trklen, trklen_min, trklen_max); trklen_reco_inel_NoCut->Sumw2();
	trklen_reco_el_NoCut = new TH1D("trklen_reco_el_NoCut","",n_trklen, trklen_min, trklen_max); trklen_reco_el_NoCut->Sumw2();
	trklen_reco_midcosmic_NoCut = new TH1D("trklen_reco_midcosmic_NoCut","",n_trklen, trklen_min, trklen_max); trklen_reco_midcosmic_NoCut->Sumw2();
	trklen_reco_midpi_NoCut = new TH1D("trklen_reco_midpi_NoCut","",n_trklen, trklen_min, trklen_max); trklen_reco_midpi_NoCut->Sumw2();
	trklen_reco_midp_NoCut = new TH1D("trklen_reco_midp_NoCut","",n_trklen, trklen_min, trklen_max); trklen_reco_midp_NoCut->Sumw2();
	trklen_reco_midmu_NoCut = new TH1D("trklen_reco_midmu_NoCut","",n_trklen, trklen_min, trklen_max); trklen_reco_midmu_NoCut->Sumw2();
	trklen_reco_mideg_NoCut = new TH1D("trklen_reco_mideg_NoCut","",n_trklen, trklen_min, trklen_max); trklen_reco_mideg_NoCut->Sumw2();
	trklen_reco_midother_NoCut = new TH1D("trklen_reco_midother_NoCut","",n_trklen, trklen_min, trklen_max); trklen_reco_midother_NoCut->Sumw2();

	//pandora cut
	trklen_reco_inel_PanS = new TH1D("trklen_reco_inel_PanS","",n_trklen, trklen_min, trklen_max); trklen_reco_inel_PanS->Sumw2();
	trklen_reco_el_PanS = new TH1D("trklen_reco_el_PanS","",n_trklen, trklen_min, trklen_max); trklen_reco_el_PanS->Sumw2();
	trklen_reco_midcosmic_PanS = new TH1D("trklen_reco_midcosmic_PanS","",n_trklen, trklen_min, trklen_max); trklen_reco_midcosmic_PanS->Sumw2();
	trklen_reco_midpi_PanS = new TH1D("trklen_reco_midpi_PanS","",n_trklen, trklen_min, trklen_max); trklen_reco_midpi_PanS->Sumw2();
	trklen_reco_midp_PanS = new TH1D("trklen_reco_midp_PanS","",n_trklen, trklen_min, trklen_max); trklen_reco_midp_PanS->Sumw2();
	trklen_reco_midmu_PanS = new TH1D("trklen_reco_midmu_PanS","",n_trklen, trklen_min, trklen_max); trklen_reco_midmu_PanS->Sumw2();
	trklen_reco_mideg_PanS = new TH1D("trklen_reco_mideg_PanS","",n_trklen, trklen_min, trklen_max); trklen_reco_mideg_PanS->Sumw2();
	trklen_reco_midother_PanS = new TH1D("trklen_reco_midother_PanS","",n_trklen, trklen_min, trklen_max); trklen_reco_midother_PanS->Sumw2();

	//CaloSz
	trklen_reco_inel_CaloSz = new TH1D("trklen_reco_inel_CaloSz","",n_trklen, trklen_min, trklen_max); trklen_reco_inel_CaloSz->Sumw2();
	trklen_reco_el_CaloSz = new TH1D("trklen_reco_el_CaloSz","",n_trklen, trklen_min, trklen_max); trklen_reco_el_CaloSz->Sumw2();
	trklen_reco_midcosmic_CaloSz = new TH1D("trklen_reco_midcosmic_CaloSz","",n_trklen, trklen_min, trklen_max); trklen_reco_midcosmic_CaloSz->Sumw2();
	trklen_reco_midpi_CaloSz = new TH1D("trklen_reco_midpi_CaloSz","",n_trklen, trklen_min, trklen_max); trklen_reco_midpi_CaloSz->Sumw2();
	trklen_reco_midp_CaloSz = new TH1D("trklen_reco_midp_CaloSz","",n_trklen, trklen_min, trklen_max); trklen_reco_midp_CaloSz->Sumw2();
	trklen_reco_midmu_CaloSz = new TH1D("trklen_reco_midmu_CaloSz","",n_trklen, trklen_min, trklen_max); trklen_reco_midmu_CaloSz->Sumw2();
	trklen_reco_mideg_CaloSz = new TH1D("trklen_reco_mideg_CaloSz","",n_trklen, trklen_min, trklen_max); trklen_reco_mideg_CaloSz->Sumw2();
	trklen_reco_midother_CaloSz = new TH1D("trklen_reco_midother_CaloSz","",n_trklen, trklen_min, trklen_max); trklen_reco_midother_CaloSz->Sumw2();

	//beam quality
	trklen_reco_inel_BQ = new TH1D("trklen_reco_inel_BQ","",n_trklen, trklen_min, trklen_max); trklen_reco_inel_BQ->Sumw2();
	trklen_reco_el_BQ = new TH1D("trklen_reco_el_BQ","",n_trklen, trklen_min, trklen_max); trklen_reco_el_BQ->Sumw2();
	trklen_reco_midcosmic_BQ = new TH1D("trklen_reco_midcosmic_BQ","",n_trklen, trklen_min, trklen_max); trklen_reco_midcosmic_BQ->Sumw2();
	trklen_reco_midpi_BQ = new TH1D("trklen_reco_midpi_BQ","",n_trklen, trklen_min, trklen_max); trklen_reco_midpi_BQ->Sumw2();
	trklen_reco_midp_BQ = new TH1D("trklen_reco_midp_BQ","",n_trklen, trklen_min, trklen_max); trklen_reco_midp_BQ->Sumw2();
	trklen_reco_midmu_BQ = new TH1D("trklen_reco_midmu_BQ","",n_trklen, trklen_min, trklen_max); trklen_reco_midmu_BQ->Sumw2();
	trklen_reco_mideg_BQ = new TH1D("trklen_reco_mideg_BQ","",n_trklen, trklen_min, trklen_max); trklen_reco_mideg_BQ->Sumw2();
	trklen_reco_midother_BQ = new TH1D("trklen_reco_midother_BQ","",n_trklen, trklen_min, trklen_max); trklen_reco_midother_BQ->Sumw2();

	//reco inel cut
	trklen_reco_inel_RecoInel = new TH1D("trklen_reco_inel_RecoInel","",n_trklen, trklen_min, trklen_max); trklen_reco_inel_RecoInel->Sumw2();
	trklen_reco_el_RecoInel = new TH1D("trklen_reco_el_RecoInel","",n_trklen, trklen_min, trklen_max); trklen_reco_el_RecoInel->Sumw2();
	trklen_reco_midcosmic_RecoInel = new TH1D("trklen_reco_midcosmic_RecoInel","",n_trklen, trklen_min, trklen_max); trklen_reco_midcosmic_RecoInel->Sumw2();
	trklen_reco_midpi_RecoInel = new TH1D("trklen_reco_midpi_RecoInel","",n_trklen, trklen_min, trklen_max); trklen_reco_midpi_RecoInel->Sumw2();
	trklen_reco_midp_RecoInel = new TH1D("trklen_reco_midp_RecoInel","",n_trklen, trklen_min, trklen_max); trklen_reco_midp_RecoInel->Sumw2();
	trklen_reco_midmu_RecoInel = new TH1D("trklen_reco_midmu_RecoInel","",n_trklen, trklen_min, trklen_max); trklen_reco_midmu_RecoInel->Sumw2();
	trklen_reco_mideg_RecoInel = new TH1D("trklen_reco_mideg_RecoInel","",n_trklen, trklen_min, trklen_max); trklen_reco_mideg_RecoInel->Sumw2();
	trklen_reco_midother_RecoInel = new TH1D("trklen_reco_midother_RecoInel","",n_trklen, trklen_min, trklen_max); trklen_reco_midother_RecoInel->Sumw2();

	//dtrklen -----------------------------------------------------------------------------------------------------//
	int n_dtrklen=100;
	float dtrklen_st=-100;
	float dtrklen_end=100; 

	//nocut
	dtrklen_inel_NoCut = new TH1D("dtrklen_inel_NoCut","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_inel_NoCut->Sumw2();
	dtrklen_el_NoCut = new TH1D("dtrklen_el_NoCut","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_el_NoCut->Sumw2();
	dtrklen_midcosmic_NoCut = new TH1D("dtrklen_midcosmic_NoCut","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midcosmic_NoCut->Sumw2();
	dtrklen_midpi_NoCut = new TH1D("dtrklen_midpi_NoCut","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midpi_NoCut->Sumw2();
	dtrklen_midp_NoCut = new TH1D("dtrklen_midp_NoCut","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midp_NoCut->Sumw2();
	dtrklen_midmu_NoCut = new TH1D("dtrklen_midmu_NoCut","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midmu_NoCut->Sumw2();
	dtrklen_mideg_NoCut = new TH1D("dtrklen_mideg_NoCut","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_mideg_NoCut->Sumw2();
	dtrklen_midother_NoCut = new TH1D("dtrklen_midother_NoCut","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midother_NoCut->Sumw2();

	//pandora cut
	dtrklen_inel_PanS = new TH1D("dtrklen_inel_PanS","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_inel_PanS->Sumw2();
	dtrklen_el_PanS = new TH1D("dtrklen_el_PanS","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_el_PanS->Sumw2();
	dtrklen_midcosmic_PanS = new TH1D("dtrklen_midcosmic_PanS","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midcosmic_PanS->Sumw2();
	dtrklen_midpi_PanS = new TH1D("dtrklen_midpi_PanS","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midpi_PanS->Sumw2();
	dtrklen_midp_PanS = new TH1D("dtrklen_midp_PanS","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midp_PanS->Sumw2();
	dtrklen_midmu_PanS = new TH1D("dtrklen_midmu_PanS","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midmu_PanS->Sumw2();
	dtrklen_mideg_PanS = new TH1D("dtrklen_mideg_PanS","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_mideg_PanS->Sumw2();
	dtrklen_midother_PanS = new TH1D("dtrklen_midother_PanS","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midother_PanS->Sumw2();
	
	//calosz
	dtrklen_inel_CaloSz = new TH1D("dtrklen_inel_CaloSz","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_inel_CaloSz->Sumw2();
	dtrklen_el_CaloSz = new TH1D("dtrklen_el_CaloSz","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_el_CaloSz->Sumw2();
	dtrklen_midcosmic_CaloSz = new TH1D("dtrklen_midcosmic_CaloSz","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midcosmic_CaloSz->Sumw2();
	dtrklen_midpi_CaloSz = new TH1D("dtrklen_midpi_CaloSz","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midpi_CaloSz->Sumw2();
	dtrklen_midp_CaloSz = new TH1D("dtrklen_midp_CaloSz","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midp_CaloSz->Sumw2();
	dtrklen_midmu_CaloSz = new TH1D("dtrklen_midmu_CaloSz","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midmu_CaloSz->Sumw2();
	dtrklen_mideg_CaloSz = new TH1D("dtrklen_mideg_CaloSz","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_mideg_CaloSz->Sumw2();
	dtrklen_midother_CaloSz = new TH1D("dtrklen_midother_CaloSz","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midother_CaloSz->Sumw2();

	//beam quality
	dtrklen_inel_BQ = new TH1D("dtrklen_inel_BQ","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_inel_BQ->Sumw2();
	dtrklen_el_BQ = new TH1D("dtrklen_el_BQ","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_el_BQ->Sumw2();
	dtrklen_midcosmic_BQ = new TH1D("dtrklen_midcosmic_BQ","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midcosmic_BQ->Sumw2();
	dtrklen_midpi_BQ = new TH1D("dtrklen_midpi_BQ","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midpi_BQ->Sumw2();
	dtrklen_midp_BQ = new TH1D("dtrklen_midp_BQ","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midp_BQ->Sumw2();
	dtrklen_midmu_BQ = new TH1D("dtrklen_midmu_BQ","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midmu_BQ->Sumw2();
	dtrklen_mideg_BQ = new TH1D("dtrklen_mideg_BQ","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_mideg_BQ->Sumw2();
	dtrklen_midother_BQ = new TH1D("dtrklen_midother_BQ","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midother_BQ->Sumw2();

	//reco inel
	dtrklen_inel_RecoInel = new TH1D("dtrklen_inel_RecoInel","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_inel_RecoInel->Sumw2();
	dtrklen_el_RecoInel = new TH1D("dtrklen_el_RecoInel","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_el_RecoInel->Sumw2();
	dtrklen_midcosmic_RecoInel = new TH1D("dtrklen_midcosmic_RecoInel","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midcosmic_RecoInel->Sumw2();
	dtrklen_midpi_RecoInel = new TH1D("dtrklen_midpi_RecoInel","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midpi_RecoInel->Sumw2();
	dtrklen_midp_RecoInel = new TH1D("dtrklen_midp_RecoInel","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midp_RecoInel->Sumw2();
	dtrklen_midmu_RecoInel = new TH1D("dtrklen_midmu_RecoInel","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midmu_RecoInel->Sumw2();
	dtrklen_mideg_RecoInel = new TH1D("dtrklen_mideg_RecoInel","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_mideg_RecoInel->Sumw2();
	dtrklen_midother_RecoInel = new TH1D("dtrklen_midother_RecoInel","",n_dtrklen, dtrklen_st, dtrklen_end); dtrklen_midother_RecoInel->Sumw2();

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
        int dke=325;
        float ke_st=-50;
        float ke_end=600;
	trklen_ke_true_inel=new TH2D("trklen_ke_true_inel","",n_ntrklen, st_ntrklen, ed_ntrklen, dke, ke_st, ke_end);
	trklen_ke_true_el=new TH2D("trklen_ke_true_el","",n_ntrklen, st_ntrklen, ed_ntrklen, dke, ke_st, ke_end);

	//truth inc/int
	h_true_incidents=new TH1D("h_true_incidents","true_incidents", nthinslices, 0, nthinslices-1);
	h_true_interactions=new TH1D("h_true_interactions","true_interactions", nthinslices, 0, nthinslices-1);
	h_true_incidents->Sumw2();
	h_true_interactions->Sumw2();


	//reco E-dept
	reco_dedx_trklen_inel=new TH2D("reco_dedx_trklen_inel","reco_dedx_trklen_inel",150,0,150, 100,0,100);
	reco_de_trklen_inel=new TH2D("reco_de_trklen_inel","reco_de_trklen_inel",150,0,150, 100,0,100);
	reco_dx_trklen_inel=new TH2D("reco_dx_trklen_inel","reco_dx_trklen_inel",150,0,150,100,0,10);
	reco_dedx_trklen_inel->Sumw2();
	reco_de_trklen_inel->Sumw2();
	reco_dx_trklen_inel->Sumw2();

	

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

