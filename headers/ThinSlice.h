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


//class ThinSlice {

	//public:

		int reco_sliceID;
		int true_sliceID;

		TH1D *reco_incE[nthinslices];
		TH1D *true_incE[nthinslices];
		TH1D *reco_AngCorr;
		TH1D *true_AngCorr;

		TH1D *h_truesliceid_all;
		TH1D *h_truesliceid_uf;
		TH1D *h_truesliceid_cuts;
		TH1D *h_truesliceid_inelastic_all;
		TH1D *h_truesliceid_inelastic_uf;
		TH1D *h_truesliceid_inelastic_cuts;
		TH1D *h_recosliceid_allevts_cuts;
		TH1D *h_recosliceid_cuts;
		TH1D *h_recosliceid_inelastic_cuts;


		TH1D *zend_true;
		TH1D *zend_inel_true;
		TH1D *zend_el_true;
		TH1D *zend_mcs_true;
		TH1D *zend_bkg_true;

		TH1D *zend_inel_true_beamQ;
		TH1D *zend_el_true_beamQ;
		TH1D *zend_mcs_true_beamQ;
		TH1D *zend_bkg_true_beamQ;

		TH1D *zend_inel_true_caloSz;
		TH1D *zend_el_true_caloSz;
		TH1D *zend_mcs_true_caloSz;
		TH1D *zend_bkg_true_caloSz;

		TH1D *zend_inel_true_RecoInel;
		TH1D *zend_el_true_RecoInel;
		TH1D *zend_mcs_true_RecoInel;
		TH1D *zend_bkg_true_RecoInel;

		TH1D *zend_inel_true_XY;
		TH1D *zend_el_true_XY;
		TH1D *zend_mcs_true_XY;
		TH1D *zend_bkg_true_XY;

		TH1D *trklen_inel_true_XY;
		TH1D *trklen_el_true_XY;
		TH1D *trklen_mcs_true_XY;
		TH1D *trklen_bkg_true_XY;

		TH1D *dzend_inel;
		TH1D *dzend_el;
		TH1D *dzend_mcs;
		TH1D *dzend_bkg;

		TH1D *dzend_inel_beamQ;
		TH1D *dzend_el_beamQ;
		TH1D *dzend_mcs_beamQ;
		TH1D *dzend_bkg_beamQ;

		TH1D *dzend_inel_caloSz;
		TH1D *dzend_el_caloSz;
		TH1D *dzend_mcs_caloSz;
		TH1D *dzend_bkg_caloSz;

		TH1D *dzend_inel_RecoInel;
		TH1D *dzend_el_RecoInel;
		TH1D *dzend_mcs_RecoInel;
		TH1D *dzend_bkg_RecoInel;

		TH1D *dzend_inel_XY;
		TH1D *dzend_el_XY;
		TH1D *dzend_mcs_XY;
		TH1D *dzend_bkg_XY;

		TH2D *zend_2d_inel;
		TH2D *zend_2d_el;
		TH2D *zend_2d_mcs;
		TH2D *zend_2d_bkg;

		TH2D *zend_2d_inel_XY;
		TH2D *zend_2d_el_XY;
		TH2D *zend_2d_mcs_XY;
		TH2D *zend_2d_bkg_XY;

		TH1D *zend_reco;

		TH1D *trklen_inel_true;
		TH1D *trklen_el_true;
		TH1D *trklen_mcs_true;
		TH1D *trklen_bkg_true;

		TH1D *trklen_inel_reco;
		TH1D *trklen_el_reco;
		TH1D *trklen_mcs_reco;
		TH1D *trklen_bkg_reco;


		TH2D *zend_ke_inel_true_XY;
		TH2D *zend_ke_el_true_XY;
		TH2D *zend_ke_mcs_true_XY;
		TH2D *zend_ke_bkg_true_XY;


	TH2D *ztrue_ketrue_TrueInEL;
	TH2D *ztrue_ketrue_RecoInEL;
	//TH2D *zreco_kereco_RecoInEL;
	//TH2D *zreco_kereco_TrueInEL;

	TH1D *ketrue_TrueInEL;

	TH2D *ztrue_ketrue_TrueMCS;
	TH1D *ketrue_TrueMCS;
	TH1D *dztrue_TrueMCS;

	TH2D *ztrue_ketrue_TrueEL;
	TH1D *ketrue_TrueEL;

	TH2D *ztrue_ketrue;


		//reco presentation with true labels
		TH2D *zreco_kereco_TrueInEL;
		TH2D *zreco_kereco_TrueEL;
		TH2D *zreco_kereco_TrueMCS;

		TH1D *zreco_TrueInEL;
		TH1D *zreco_TrueEL;
		TH1D *zreco_TrueMCS;

		TH1D *kereco_TrueInEL;
		TH1D *kereco_TrueEL;
		TH1D *kereco_TrueMCS;

		//reco presentation with reco labels
		TH2D *zreco_kereco_RecoInEL;
		TH1D *zreco_RecoInEL;
		TH1D *kereco_RecoInEL;


		TH2D *rangereco_dedxreco_TrueInEL;
		TH2D *rangereco_dedxreco_TrueEL;
		TH2D *rangereco_dedxreco_TrueMCS;


		TH1D *KE_ff_recostop;
		TH1D *KE_range_recostop;
		TH1D *KE_calo_recostop;
		TH1D *KE_simide_recostop;
		TH1D *KE_rrange_recostop;
		TH1D *KE_rrange2_recostop;
		//TH1D *KE_truth_recostop;
		
		TH1D *dKE_range_ff_recostop;
		TH1D *dKE_calo_ff_recostop;
		TH1D *dKE_rrange_ff_recostop;
		TH1D *dKE_rrange2_ff_recostop;


                TH1D *KE_ff_recoinel;
                TH1D *KE_range_recoinel;
                TH1D *KE_calo_recoinel;
                TH1D *KE_simide_recoinel;
                TH1D *KE_rrange_recoinel;

		TH2D *rr_dedx_recostop;
		TH2D *rr_ddedx1_recostop;
		TH2D *rr_ddedx2_recostop;

		TH2D *rr_dedx_truestop;


		double true_interactions[nthinslices];
		double true_incidents[nthinslices];

		std::string fOutputFileName;
		TFile *outputFile;
		void SetOutputFileName(std::string name){fOutputFileName = name;};
		void BookHistograms();
		//void FillHistograms(int cut, const HadAna & evt);
		void SaveHistograms();

		//void ProcessEvent(const HadAna & evt, Unfold & uf);
		//void CalcXS(const Unfold & uf);

		//void Run(HadAna & evt, Unfold & uf);

//};


/*
   void FillHistVec1D(TH1D *hist[nParTypes+1], const double &value, const int &partype){
   hist[0]->Fill(value);
   if (partype>=1 && partype < nParTypes+1){
   hist[partype]->Fill(value);
   }
   }

   void FillHistVec2D(TH2D *hist[nParTypes+1], const double &value1, const double &value2, const int &partype){
   hist[0]->Fill(value1, value2);
   if (partype>=1 && partype < nParTypes+1){
   hist[partype]->Fill(value1, value2);
   }
   }
   */

void BookHistograms(){

	outputFile = TFile::Open(fOutputFileName.c_str(), "recreate");

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
	h_recosliceid_cuts = new TH1D("h_recosliceid_cuts","h_recosliceid_cuts;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_inelastic_cuts = new TH1D("h_recosliceid_inelastic_cuts","h_recosliceid_inelastic_cuts;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);

	h_truesliceid_all->Sumw2();
	h_truesliceid_cuts->Sumw2();
	h_truesliceid_inelastic_all->Sumw2();
	h_truesliceid_inelastic_cuts->Sumw2();
	h_recosliceid_allevts_cuts->Sumw2();
	h_recosliceid_cuts->Sumw2();
	h_recosliceid_inelastic_cuts->Sumw2();

	/*
	   for (int i = 0; i < nCuts; ++i){
	   for (int j = 0; j < nParTypes+1; ++j){
	   htrue_beam_endZ[i][j] = new TH1D(Form("htrue_beam_endZ_%d_%d",i,j),Form("true_beam_endZ, %s, %s;true_beam_endZ (cm)", cutName[i], parTypeName[j]), 70, -100, 600);
	   htrue_beam_endZ[i][j]->Sumw2();
	   hreco_beam_endZ[i][j] = new TH1D(Form("hreco_beam_endZ_%d_%d",i,j),Form("reco_beam_endZ, %s, %s;reco_beam_endZ (cm)", cutName[i], parTypeName[j]), 70, -100, 600);
	   hreco_beam_endZ[i][j]->Sumw2();
	   hreco_true_beam_endZ[i][j] = new TH1D(Form("hreco_true_beam_endZ_%d_%d",i,j), Form("reco_true_beam_endZ, %s, %s;reco_beam_endZ - true_beam_endZ (cm)", cutName[i], parTypeName[j]), 100, -100, 100);
	   hreco_true_beam_endZ[i][j]->Sumw2();
	   hreco_vs_true_beam_endZ[i][j]= new TH2D(Form("hreco_vs_true_beam_endZ_%d_%d",i,j), Form("%s, %s;true_beam_endZ (cm);reco_beam_endZ (cm)", cutName[i], parTypeName[j]), 70, -100, 600, 70, -100, 600);
	   hreco_true_vs_true_beam_endZ[i][j]= new TH2D(Form("hreco_true_vs_true_beam_endZ_%d_%d",i,j), Form("%s, %s;true_beam_endZ (cm);reco - true_beam_endZ (cm)", cutName[i], parTypeName[j]), 70, -100, 600, 100, -100, 100);

	   htrue_beam_endZ_SCE[i][j] = new TH1D(Form("htrue_beam_endZ_SCE_%d_%d",i,j),Form("true_beam_endZ_SCE, %s, %s;true_beam_endZ_SCE (cm)", cutName[i], parTypeName[j]), 70, -100, 600);
	   htrue_beam_endZ_SCE[i][j]->Sumw2();
	   hreco_beam_endZ_SCE[i][j] = new TH1D(Form("hreco_beam_endZ_SCE_%d_%d",i,j),Form("reco_beam_endZ_SCE, %s, %s;reco_beam_endZ_SCE (cm)", cutName[i], parTypeName[j]), 70, -100, 600);
	   hreco_beam_endZ_SCE[i][j]->Sumw2();
	   hreco_true_beam_endZ_SCE[i][j] = new TH1D(Form("hreco_true_beam_endZ_SCE_%d_%d",i,j), Form("reco_true_beam_endZ_SCE, %s, %s;reco_beam_endZ_SCE - true_beam_endZ_SCE (cm)", cutName[i], parTypeName[j]), 100, -100, 100);
	   hreco_true_beam_endZ_SCE[i][j]->Sumw2();
	   hreco_vs_true_beam_endZ_SCE[i][j]= new TH2D(Form("hreco_vs_true_beam_endZ_SCE_%d_%d",i,j), Form("%s, %s;true_beam_endZ_SCE (cm);reco_beam_endZ_SCE (cm)", cutName[i], parTypeName[j]), 70, -100, 600, 70, -100, 600);
	   hreco_true_vs_true_beam_endZ_SCE[i][j]= new TH2D(Form("hreco_true_vs_true_beam_endZ_SCE_%d_%d",i,j), Form("%s, %s;true_beam_endZ_SCE (cm);reco - true_beam_endZ_SCE (cm)", cutName[i], parTypeName[j]), 70, -100, 600, 100, -100, 100);


	   htrue_sliceID[i][j] = new TH1D(Form("htrue_sliceID_%d_%d",i,j),Form("true_sliceID, %s, %s;true_sliceID (cm)", cutName[i], parTypeName[j]), 50, -1, 49);
	   htrue_sliceID[i][j]->Sumw2();
	   hreco_sliceID[i][j] = new TH1D(Form("hreco_sliceID_%d_%d",i,j),Form("reco_sliceID, %s, %s;reco_sliceID", cutName[i], parTypeName[j]), 50, -1, 49);
	   hreco_sliceID[i][j]->Sumw2();
	   hreco_true_sliceID[i][j] = new TH1D(Form("hreco_true_sliceID_%d_%d",i,j), Form("reco_true_sliceID, %s, %s;reco_sliceID - true_sliceID", cutName[i], parTypeName[j]), 20, -10, 10);
	   hreco_true_sliceID[i][j]->Sumw2();
	   hreco_vs_true_sliceID[i][j]= new TH2D(Form("hreco_vs_true_sliceID_%d_%d",i,j), Form("%s, %s;true_sliceID;reco_sliceID", cutName[i], parTypeName[j]), 50, -1, 49, 50, -1, 49);
	   hreco_true_vs_true_sliceID[i][j]= new TH2D(Form("hreco_true_vs_true_sliceID_%d_%d",i,j), Form("%s, %s;true_sliceID;reco_sliceID - true_sliceID", cutName[i], parTypeName[j]), 50, -1, 49, 20, -10, 10);

	   hmediandEdx[i][j] = new TH1D(Form("hmediandEdx_%d_%d",i,j), Form("mediandEdx, %s, %s;Median dE/dx (MeV/cm)", cutName[i], parTypeName[j]), 100, 0, 5);
	   hmediandEdx[i][j]->Sumw2();

	   hdaughter_michel_score[i][j] = new TH1D(Form("hdaughter_michel_score_%d_%d",i,j), Form("daughter_michel_score, %s, %s;Michel score", cutName[i], parTypeName[j]), 100, 0, 1);
	   hdaughter_michel_score[i][j]->Sumw2();

	   }
	   }
	   */

	for (int i = 0; i<nthinslices; ++i){
		true_interactions[i] = 0;
		true_incidents[i] = 0;
	}

	//response_SliceID_Pion = new RooUnfoldResponse(nthinslices+2, -1, nthinslices+1, "response_SliceID_Pion");
	//response_SliceID_PionInEl = new RooUnfoldResponse(nthinslices+2, -1, nthinslices+1, "response_SliceID_PionInEl");


	//zend distributions
	int dz=2000;
	float z_st=-50;
	float z_end=150;

	//ke dist
	//int dke=650;
	int dke=325;
	float ke_st=-50;
	float ke_end=600;

	//truth distributions
	zend_true = new TH1D("zend_true","", dz, z_st, z_end);
	zend_inel_true = new TH1D("zend_inel_true","", dz, z_st, z_end);
	zend_el_true = new TH1D("zend_el_true","", dz, z_st, z_end);
	zend_mcs_true = new TH1D("zend_mcs_true","", dz, z_st, z_end);
	zend_bkg_true = new TH1D("zend_bkg_true","", dz, z_st, z_end);

	//(par)_(truthlabel)_(reco or truth par)_(cut)
	zend_inel_true_beamQ = new TH1D("zend_inel_true_beamQ","", dz, z_st, z_end);
	zend_el_true_beamQ = new TH1D("zend_el_true_beamQ","", dz, z_st, z_end);
	zend_mcs_true_beamQ = new TH1D("zend_mcs_true_beamQ","", dz, z_st, z_end);
	zend_bkg_true_beamQ = new TH1D("zend_bkg_true_beamQ","", dz, z_st, z_end);

	zend_inel_true_caloSz = new TH1D("zend_inel_true_caloSz","", dz, z_st, z_end);
	zend_el_true_caloSz = new TH1D("zend_el_true_caloSz","", dz, z_st, z_end);
	zend_mcs_true_caloSz = new TH1D("zend_mcs_true_caloSz","", dz, z_st, z_end);
	zend_bkg_true_caloSz = new TH1D("zend_bkg_true_caloSz","", dz, z_st, z_end);

	zend_inel_true_RecoInel = new TH1D("zend_inel_true_RecoInel","", dz, z_st, z_end);
	zend_el_true_RecoInel = new TH1D("zend_el_true_RecoInel","", dz, z_st, z_end);
	zend_mcs_true_RecoInel = new TH1D("zend_mcs_true_RecoInel","", dz, z_st, z_end);
	zend_bkg_true_RecoInel = new TH1D("zend_bkg_true_RecoInel","", dz, z_st, z_end);

	zend_inel_true_XY = new TH1D("zend_inel_true_XY","", dz, z_st, z_end);
	zend_el_true_XY = new TH1D("zend_el_true_XY","", dz, z_st, z_end);
	zend_mcs_true_XY = new TH1D("zend_mcs_true_XY","", dz, z_st, z_end);
	zend_bkg_true_XY = new TH1D("zend_bkg_true_XY","", dz, z_st, z_end);

	trklen_inel_true_XY = new TH1D("trklen_inel_true_XY","", dz, z_st, z_end);
	trklen_el_true_XY = new TH1D("trklen_el_true_XY","", dz, z_st, z_end);
	trklen_mcs_true_XY = new TH1D("trklen_mcs_true_XY","", dz, z_st, z_end);
	trklen_bkg_true_XY = new TH1D("trklen_bkg_true_XY","", dz, z_st, z_end);

	//dZend (reco-true)
	int ddz=300;
	float dz_st=-150;
	float dz_end=150;

	dzend_inel = new TH1D("dzend_inel","", ddz, dz_st, dz_end);
	dzend_el = new TH1D("dzend_el","", ddz, dz_st, dz_end);
	dzend_mcs = new TH1D("dzend_mcs","", ddz, dz_st, dz_end);
	dzend_bkg = new TH1D("dzend_bkg","", ddz, dz_st, dz_end);

	dzend_inel_beamQ = new TH1D("dzend_inel_beamQ","", ddz, dz_st, dz_end);
	dzend_el_beamQ = new TH1D("dzend_el_beamQ","", ddz, dz_st, dz_end);
	dzend_mcs_beamQ = new TH1D("dzend_mcs_beamQ","", ddz, dz_st, dz_end);
	dzend_bkg_beamQ = new TH1D("dzend_bkg_beamQ","", ddz, dz_st, dz_end);

	dzend_inel_caloSz = new TH1D("dzend_inel_caloSz","", ddz, dz_st, dz_end);
	dzend_el_caloSz = new TH1D("dzend_el_caloSz","", ddz, dz_st, dz_end);
	dzend_mcs_caloSz = new TH1D("dzend_mcs_caloSz","", ddz, dz_st, dz_end);
	dzend_bkg_caloSz = new TH1D("dzend_bkg_caloSz","", ddz, dz_st, dz_end);

	dzend_inel_RecoInel = new TH1D("dzend_inel_RecoInel","", ddz, dz_st, dz_end);
	dzend_el_RecoInel = new TH1D("dzend_el_RecoInel","", ddz, dz_st, dz_end);
	dzend_mcs_RecoInel = new TH1D("dzend_mcs_RecoInel","", ddz, dz_st, dz_end);
	dzend_bkg_RecoInel = new TH1D("dzend_bkg_RecoInel","", ddz, dz_st, dz_end);


	dzend_inel_XY = new TH1D("dzend_inel_XY","", ddz, dz_st, dz_end);
	dzend_el_XY = new TH1D("dzend_el_XY","", ddz, dz_st, dz_end);
	dzend_mcs_XY = new TH1D("dzend_mcs_XY","", ddz, dz_st, dz_end);
	dzend_bkg_XY = new TH1D("dzend_bkg_XY","", ddz, dz_st, dz_end);

	//TH2D *zend_reco_true = new TH2D("zend_reco_true","",dz, z_st, z_end, dz, z_st, z_end);

	zend_2d_inel = new TH2D("zend_2d_inel","",dz, z_st, z_end, dz, z_st, z_end);
	zend_2d_el = new TH2D("zend_2d_el","",dz, z_st, z_end, dz, z_st, z_end);
	zend_2d_mcs = new TH2D("zend_2d_mcs","",dz, z_st, z_end, dz, z_st, z_end);
	zend_2d_bkg = new TH2D("zend_2d_bkg","",dz, z_st, z_end, dz, z_st, z_end);

	zend_2d_inel_XY = new TH2D("zend_2d_inel_XY","",dz, z_st, z_end, dz, z_st, z_end);
	zend_2d_el_XY = new TH2D("zend_2d_el_XY","",dz, z_st, z_end, dz, z_st, z_end);
	zend_2d_mcs_XY = new TH2D("zend_2d_mcs_XY","",dz, z_st, z_end, dz, z_st, z_end);
	zend_2d_bkg_XY = new TH2D("zend_2d_bkg_XY","",dz, z_st, z_end, dz, z_st, z_end);

	zend_reco = new TH1D("zend_reco","", dz, z_st, z_end);
	//dzend = new TH1D("dzend","", 1500, -150, 150);

	//trklen_true = new TH1D("trklen_true","", dz, z_st, z_end);
	trklen_inel_true = new TH1D("trklen_inel_true","", dz, z_st, z_end);
	trklen_el_true = new TH1D("trklen_el_true","", dz, z_st, z_end);
	trklen_mcs_true = new TH1D("trklen_mcs_true","", dz, z_st, z_end);
	trklen_bkg_true = new TH1D("trklen_bkg_true","", dz, z_st, z_end);

	trklen_inel_reco = new TH1D("trklen_inel_reco","", dz, z_st, z_end);
	trklen_el_reco = new TH1D("trklen_el_reco","", dz, z_st, z_end);
	trklen_mcs_reco = new TH1D("trklen_mcs_reco","", dz, z_st, z_end);
	trklen_bkg_reco = new TH1D("trklen_bkg_reco","", dz, z_st, z_end);


	zend_ke_inel_true_XY = new TH2D("zend_ke_inel_true_XY","",dz, z_st, z_end, dke, ke_st, ke_end);
	zend_ke_el_true_XY = new TH2D("zend_ke_el_true_XY","",dz, z_st, z_end, dke, ke_st, ke_end);
	zend_ke_mcs_true_XY = new TH2D("zend_ke_mcs_true_XY","",dz, z_st, z_end, dke, ke_st, ke_end);
	zend_ke_bkg_true_XY = new TH2D("zend_ke_bkg_true_XY","",dz, z_st, z_end, dke, ke_st, ke_end);



	ztrue_ketrue_TrueInEL = new TH2D("ztrue_ketrue_TrueInEL","",dz, z_st, z_end, dke, ke_st, ke_end);
	ztrue_ketrue_RecoInEL = new TH2D("ztrue_ketrue_RecoInEL","",dz, z_st, z_end, dke, ke_st, ke_end);
	zreco_kereco_RecoInEL = new TH2D("zreco_kereco_RecoInEL","",dz, z_st, z_end, dke, ke_st, ke_end);

	zreco_kereco_TrueInEL = new TH2D("zreco_kereco_TrueInEL","",dz, z_st, z_end, dke, ke_st, ke_end);
	zreco_kereco_TrueEL = new TH2D("zreco_kereco_TrueEL","",dz, z_st, z_end, dke, ke_st, ke_end);
	zreco_kereco_TrueMCS = new TH2D("zreco_kereco_TrueMCS","",dz, z_st, z_end, dke, ke_st, ke_end);

	zreco_TrueInEL = new TH1D("zreco_TrueInEL","", dz, z_st, z_end);
	zreco_TrueEL = new TH1D("zreco_TrueEL","", dz, z_st, z_end);
	zreco_TrueMCS = new TH1D("zreco_TrueMCS","", dz, z_st, z_end);

	kereco_TrueInEL = new TH1D("kereco_TrueInEL","", dke, ke_st, ke_end);
	kereco_TrueEL = new TH1D("kereco_TrueEL","", dke, ke_st, ke_end);
	kereco_TrueMCS = new TH1D("kereco_TrueMCS","", dke, ke_st, ke_end);

	zreco_RecoInEL = new TH1D("zreco_RecoInEL","", dz, z_st, z_end);

	kereco_RecoInEL = new TH1D("kereco_RecoInEL","", dke, ke_st, ke_end);


	ketrue_TrueInEL = new TH1D("ketrue_TrueInEL","",dke, ke_st, ke_end);

	ztrue_ketrue_TrueMCS = new TH2D("ztrue_ketrue_TrueMCS","",dz, z_st, z_end, dke, ke_st, ke_end);
	ketrue_TrueMCS = new TH1D("ketrue_TrueMCS","",dke, ke_st, ke_end);
	dztrue_TrueMCS = new TH1D("dztrue_TrueMCS","",200,0,20);


	ztrue_ketrue_TrueEL = new TH2D("ztrue_ketrue_TrueEL","",dz, z_st, z_end,dke, ke_st, ke_end);
	ketrue_TrueEL = new TH1D("ketrue_TrueEL","",dke, ke_st, ke_end);

	ztrue_ketrue = new TH2D("ztrue_ketrue","",dz, z_st, z_end,dke, ke_st, ke_end);

	//zend_true->GetXaxis()->SetTitle("EndZ (truth) [cm]"); zend_true->SetLineColor(2);
	//zend_reco->GetXaxis()->SetTitle("EndZ (reco) [cm]"); zend_reco->SetLineColor(4);
	//dzend->GetXaxis()->SetTitle("#DeltaEndZ (reco-truth) [cm]"); dzend->SetLineColor(4);
	//zend_reco_true->GetXaxis()->SetTitle("EndZ (reco) [cm]");
	//zend_reco_true->GetYaxis()->SetTitle("EndZ (truth) [cm]");

	rangereco_dedxreco_TrueInEL=new TH2D("rangereco_dedxreco_TrueInEL","",200,0,100,100,0,50);
	rangereco_dedxreco_TrueEL=new TH2D("rangereco_dedxreco_TrueEL","",200,0,100,100,0,50);
	rangereco_dedxreco_TrueMCS=new TH2D("rangereco_dedxreco_TrueMCS","",200,0,100,100,0,50);

	int n_ke=140;
	float ke_min=0;
	float ke_max=700;
	KE_ff_recostop=new TH1D("KE_ff_recostop","", n_ke, ke_min, ke_max);
	KE_range_recostop=new TH1D("KE_range_recostop", "", n_ke, ke_min, ke_max);
	KE_calo_recostop=new TH1D("KE_calo_recostop","",n_ke, ke_min, ke_max);
	KE_simide_recostop=new TH1D("KE_simide_recostop","",n_ke, ke_min, ke_max);
	KE_rrange_recostop=new TH1D("KE_rrange_recostop", "", n_ke, ke_min, ke_max);
	KE_rrange2_recostop=new TH1D("KE_rrange2_recostop", "", n_ke, ke_min, ke_max);
	//KE_truth_recostop=new TH1D("KE_truth_recostop","",n_ke, ke_min, ke_max);

	dKE_range_ff_recostop=new TH1D("dKE_range_ff_recostop","",200,-100,100);
	dKE_calo_ff_recostop=new TH1D("dKE_calo_ff_recostop","",200,-100,100);
	dKE_rrange_ff_recostop=new TH1D("dKE_rrange_ff_recostop","",200,-100,100);
	dKE_rrange2_ff_recostop=new TH1D("dKE_rrange2_ff_recostop","",200,-100,100);




        KE_ff_recoinel=new TH1D("KE_ff_recoinel","", n_ke, ke_min, ke_max);
        KE_range_recoinel=new TH1D("KE_range_recoinel", "", n_ke, ke_min, ke_max);
        KE_calo_recoinel=new TH1D("KE_calo_recoinel","",n_ke, ke_min, ke_max);
        KE_simide_recoinel=new TH1D("KE_simide_recoinel","",n_ke, ke_min, ke_max);
        KE_rrange_recoinel=new TH1D("KE_rrange_recoinel", "", n_ke, ke_min, ke_max);

	rr_dedx_recostop=new TH2D("rr_dedx_recostop","", 240,0,120, 300,0, 30);
	rr_dedx_truestop=new TH2D("rr_dedx_truestop","", 240,0,120, 300,0, 30);
	rr_ddedx1_recostop=new TH2D("rr_ddedx1_recostop","", 240,0,120, 600,-30, 30);
	rr_ddedx2_recostop=new TH2D("rr_ddedx2_recostop","", 240,0,120, 600,-30, 30);

}

/*
   void FillHistograms(int cut, const HadAna & evt){

   if (cut>=0 && cut < nCuts){
   FillHistVec1D(htrue_beam_endZ[cut], evt.true_beam_endZ_SCE, evt.partype);
   FillHistVec1D(htrue_beam_endZ_SCE[cut], evt.true_beam_endZ, evt.partype);
   FillHistVec1D(htrue_sliceID[cut], true_sliceID, evt.partype);
   if (!evt.reco_beam_calo_wire->empty()){
   FillHistVec1D(hreco_beam_endZ[cut], evt.reco_beam_endZ, evt.partype);
   FillHistVec1D(hreco_true_beam_endZ[cut], evt.reco_beam_endZ - evt.true_beam_endZ_SCE, evt.partype);
   FillHistVec2D(hreco_vs_true_beam_endZ[cut], evt.true_beam_endZ_SCE, evt.reco_beam_endZ, evt.partype);
   FillHistVec2D(hreco_true_vs_true_beam_endZ[cut], evt.true_beam_endZ_SCE, evt.reco_beam_endZ - evt.true_beam_endZ_SCE, evt.partype);

   FillHistVec1D(hreco_beam_endZ_SCE[cut], evt.reco_beam_calo_endZ, evt.partype);
   FillHistVec1D(hreco_true_beam_endZ_SCE[cut], evt.reco_beam_calo_endZ - evt.true_beam_endZ, evt.partype);
   FillHistVec2D(hreco_vs_true_beam_endZ_SCE[cut], evt.true_beam_endZ, evt.reco_beam_calo_endZ, evt.partype);
   FillHistVec2D(hreco_true_vs_true_beam_endZ_SCE[cut], evt.true_beam_endZ, evt.reco_beam_calo_endZ - evt.true_beam_endZ, evt.partype);

   FillHistVec1D(hreco_sliceID[cut], reco_sliceID, evt.partype);
   FillHistVec1D(hreco_true_sliceID[cut], reco_sliceID - true_sliceID, evt.partype);
   FillHistVec2D(hreco_vs_true_sliceID[cut], true_sliceID, reco_sliceID, evt.partype);
   FillHistVec2D(hreco_true_vs_true_sliceID[cut], true_sliceID, reco_sliceID - true_sliceID, evt.partype);

   FillHistVec1D(hmediandEdx[cut], evt.median_dEdx, evt.partype);
   FillHistVec1D(hdaughter_michel_score[cut], evt.daughter_michel_score, evt.partype);
   }      
   }
   }
   */

void SaveHistograms(){
	outputFile->cd();
	outputFile->Write();
	h_truesliceid_uf->Write("h_truesliceid_uf");
	h_truesliceid_inelastic_uf->Write("h_truesliceid_inelastic_uf");
	//response_SliceID_Pion->Write("response_SliceID_Pion");
	//response_SliceID_PionInEl->Write("response_SliceID_PionInEl");
}

void CalcXS(const Unfold & uf){

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
			truexs[i] = MAr/(Density*NA*thinslicewidth/true_AngCorr->GetMean())*log(true_incidents[i]/(true_incidents[i]-true_interactions[i]))*1e27;
			err_truexs[i] = MAr/(Density*NA*thinslicewidth/true_AngCorr->GetMean())*1e27*sqrt(true_interactions[i]+pow(true_interactions[i],2)/true_incidents[i])/true_incidents[i];
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
	TH1D *hint = (TH1D*)h_recosliceid_allevts_cuts->Clone("hint");
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

	h_truesliceid_uf = (TH1D*) unfold_Inc.Hreco();
	h_truesliceid_inelastic_uf = (TH1D*) unfold_Int.Hreco();

}

