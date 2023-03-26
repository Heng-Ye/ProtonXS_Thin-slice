#include "TGraphErrors.h"
#include "TVector3.h"
//#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "util.h"
#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "./Unfold3D.h"

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
/*
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
*/


/*
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
*/

/*
//reco inc_2D (st, end)
TH2D *h2d_recosliceid_allevts_cuts;
TH2D *h2d_recosliceid_allevts_cuts_inel;
TH2D *h2d_recosliceid_allevts_cuts_el;
TH2D *h2d_recosliceid_allevts_cuts_midcosmic;
TH2D *h2d_recosliceid_allevts_cuts_midpi;
TH2D *h2d_recosliceid_allevts_cuts_midp;
TH2D *h2d_recosliceid_allevts_cuts_midmu;
TH2D *h2d_recosliceid_allevts_cuts_mideg;
TH2D *h2d_recosliceid_allevts_cuts_midother;

//true inc_2D (st, end)
TH2D *h2d_truesliceid_allevts_cuts;
TH2D *h2d_truesliceid_allevts_cuts_inel;
TH2D *h2d_truesliceid_allevts_cuts_el;
TH2D *h2d_truesliceid_allevts_cuts_midcosmic;
TH2D *h2d_truesliceid_allevts_cuts_midpi;
TH2D *h2d_truesliceid_allevts_cuts_midp;
TH2D *h2d_truesliceid_allevts_cuts_midmu;
TH2D *h2d_truesliceid_allevts_cuts_mideg;
TH2D *h2d_truesliceid_allevts_cuts_midother;
*/

//reco inc_3D (st, end, int)
TH3D *h3d_recosliceid_allevts_cuts;
TH3D *h3d_recosliceid_allevts_cuts_inel;
TH3D *h3d_recosliceid_allevts_cuts_el;
TH3D *h3d_recosliceid_allevts_cuts_midcosmic;
TH3D *h3d_recosliceid_allevts_cuts_midpi;
TH3D *h3d_recosliceid_allevts_cuts_midp;
TH3D *h3d_recosliceid_allevts_cuts_midmu;
TH3D *h3d_recosliceid_allevts_cuts_mideg;
TH3D *h3d_recosliceid_allevts_cuts_midother;

//true inc_3D (st, end, int)
//TH3D *h3d_truesliceid_allevts_cuts;
//TH3D *h3d_truesliceid_allevts_cuts_inel;
//TH3D *h3d_truesliceid_allevts_cuts_el;
//TH3D *h3d_truesliceid_allevts_cuts_midp;
//TH3D *h3d_truesliceid_allevts_cuts_midcosmic;
//TH3D *h3d_truesliceid_allevts_cuts_midpi;
//TH3D *h3d_truesliceid_allevts_cuts_midmu;
//TH3D *h3d_truesliceid_allevts_cuts_mideg;
//TH3D *h3d_truesliceid_allevts_cuts_midother;



//double true_incidents[nthinslices+2];
//double true_st_incidents[nthinslices+2];
//double true_interactions[nthinslices+2];



void BookHistograms() { //BookHistograms
	outputFile = TFile::Open(fOutputFileName.c_str(), "recreate"); //open output file

	//XS histograms -------------------------------------------------------------------------------------------------------------------------------------------------------//
/*
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
*/

/*
	h2d_recosliceid_allevts_cuts = new TH2D("h2d_recosliceid_allevts_cuts","h2d_recosliceid_allevts_cuts; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_allevts_cuts_inel = new TH2D("h2d_recosliceid_allevts_cuts_inel","h2d_recosliceid_allevts_cuts_inel; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_allevts_cuts_el = new TH2D("h2d_recosliceid_allevts_cuts_el","h2d_recosliceid_allevts_cuts_el; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_allevts_cuts_midcosmic = new TH2D("h2d_recosliceid_allevts_cuts_midcosmic","h2d_recosliceid_allevts_cuts_midcosmic; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_allevts_cuts_midpi = new TH2D("h2d_recosliceid_allevts_cuts_midpi","h2d_recosliceid_allevts_cuts_midpi; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_allevts_cuts_midp = new TH2D("h2d_recosliceid_allevts_cuts_midp","h2d_recosliceid_allevts_cuts_midp; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_allevts_cuts_midmu = new TH2D("h2d_recosliceid_allevts_cuts_midmu","h2d_recosliceid_allevts_cuts_midmu; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_allevts_cuts_mideg = new TH2D("h2d_recosliceid_allevts_cuts_mideg","h2d_recosliceid_allevts_cuts_mideg; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_allevts_cuts_midother = new TH2D("h2d_recosliceid_allevts_cuts_midother","h2d_recosliceid_allevts_cuts_midother; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
*/

	h3d_recosliceid_allevts_cuts = new TH3D("h3d_recosliceid_allevts_cuts","h3d_recosliceid_allevts_cuts; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h3d_recosliceid_allevts_cuts_inel = new TH3D("h3d_recosliceid_allevts_cuts_inel","h3d_recosliceid_allevts_cuts_inel; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h3d_recosliceid_allevts_cuts_el = new TH3D("h3d_recosliceid_allevts_cuts_el","h3d_recosliceid_allevts_cuts_el; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h3d_recosliceid_allevts_cuts_midcosmic = new TH3D("h3d_recosliceid_allevts_cuts_midcosmic","h3d_recosliceid_allevts_cuts_midcosmic; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h3d_recosliceid_allevts_cuts_midpi = new TH3D("h3d_recosliceid_allevts_cuts_midpi","h3d_recosliceid_allevts_cuts_midpi; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h3d_recosliceid_allevts_cuts_midp = new TH3D("h3d_recosliceid_allevts_cuts_midp","h3d_recosliceid_allevts_cuts_midp; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h3d_recosliceid_allevts_cuts_midmu = new TH3D("h3d_recosliceid_allevts_cuts_midmu","h3d_recosliceid_allevts_cuts_midmu; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h3d_recosliceid_allevts_cuts_mideg = new TH3D("h3d_recosliceid_allevts_cuts_mideg","h3d_recosliceid_allevts_cuts_mideg; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h3d_recosliceid_allevts_cuts_midother = new TH3D("h3d_recosliceid_allevts_cuts_midother","h3d_recosliceid_allevts_cuts_midother; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);


/*
	h2d_truesliceid_allevts_cuts = new TH2D("h2d_truesliceid_allevts_cuts","h2d_truesliceid_allevts_cuts; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_truesliceid_allevts_cuts_inel = new TH2D("h2d_truesliceid_allevts_cuts_inel","h2d_truesliceid_allevts_cuts_inel; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_truesliceid_allevts_cuts_el = new TH2D("h2d_truesliceid_allevts_cuts_el","h2d_truesliceid_allevts_cuts_el; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_truesliceid_allevts_cuts_midcosmic = new TH2D("h2d_truesliceid_allevts_cuts_midcosmic","h2d_truesliceid_allevts_cuts_midcosmic; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_truesliceid_allevts_cuts_midpi = new TH2D("h2d_truesliceid_allevts_cuts_midpi","h2d_truesliceid_allevts_cuts_midpi; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_truesliceid_allevts_cuts_midp = new TH2D("h2d_truesliceid_allevts_cuts_midp","h2d_truesliceid_allevts_cuts_midp; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_truesliceid_allevts_cuts_midmu = new TH2D("h2d_truesliceid_allevts_cuts_midmu","h2d_truesliceid_allevts_cuts_midmu; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_truesliceid_allevts_cuts_mideg = new TH2D("h2d_truesliceid_allevts_cuts_mideg","h2d_truesliceid_allevts_cuts_mideg; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_truesliceid_allevts_cuts_midother = new TH2D("h2d_truesliceid_allevts_cuts_midother","h2d_truesliceid_allevts_cuts_midother; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
*/

	//h3d_truesliceid_allevts_cuts = new TH3D("h3d_truesliceid_allevts_cuts","h3d_truesliceid_allevts_cuts; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	//h3d_truesliceid_allevts_cuts_inel = new TH3D("h3d_truesliceid_allevts_cuts_inel","h3d_truesliceid_allevts_cuts_inel; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	//h3d_truesliceid_allevts_cuts_el = new TH3D("h3d_truesliceid_allevts_cuts_el","h3d_truesliceid_allevts_cuts_el; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	//h3d_truesliceid_allevts_cuts_midp = new TH3D("h3d_truesliceid_allevts_cuts_midp","h3d_truesliceid_allevts_cuts_midp; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	//h3d_truesliceid_allevts_cuts_midcosmic = new TH3D("h3d_truesliceid_allevts_cuts_midcosmic","h3d_truesliceid_allevts_cuts_midcosmic; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	//h3d_truesliceid_allevts_cuts_midpi = new TH3D("h3d_truesliceid_allevts_cuts_midpi","h3d_truesliceid_allevts_cuts_midpi; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	//h3d_truesliceid_allevts_cuts_midmu = new TH3D("h3d_truesliceid_allevts_cuts_midmu","h3d_truesliceid_allevts_cuts_midmu; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	//h3d_truesliceid_allevts_cuts_mideg = new TH3D("h3d_truesliceid_allevts_cuts_mideg","h3d_truesliceid_allevts_cuts_mideg; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	//h3d_truesliceid_allevts_cuts_midother = new TH3D("h3d_truesliceid_allevts_cuts_midother","h3d_truesliceid_allevts_cuts_midother; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);



/*
	h_truesliceid_all->Sumw2();
	h_true_st_sliceid_all->Sumw2();
	h_truesliceid_cuts->Sumw2();
	h_true_st_sliceid_cuts->Sumw2();
	h_truesliceid_inelastic_all->Sumw2();
	h_truesliceid_inelastic_cuts->Sumw2();
*/


/*
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
*/

/*
	h2d_recosliceid_allevts_cuts->Sumw2();
	h2d_recosliceid_allevts_cuts_inel->Sumw2();
	h2d_recosliceid_allevts_cuts_el->Sumw2();
	h2d_recosliceid_allevts_cuts_midcosmic->Sumw2();
	h2d_recosliceid_allevts_cuts_midpi->Sumw2();
	h2d_recosliceid_allevts_cuts_midp->Sumw2();
	h2d_recosliceid_allevts_cuts_midmu->Sumw2();
	h2d_recosliceid_allevts_cuts_mideg->Sumw2();
	h2d_recosliceid_allevts_cuts_midother->Sumw2();
*/

	h3d_recosliceid_allevts_cuts->Sumw2();
	h3d_recosliceid_allevts_cuts_inel->Sumw2();
	h3d_recosliceid_allevts_cuts_el->Sumw2();
	h3d_recosliceid_allevts_cuts_midcosmic->Sumw2();
	h3d_recosliceid_allevts_cuts_midpi->Sumw2();
	h3d_recosliceid_allevts_cuts_midp->Sumw2();
	h3d_recosliceid_allevts_cuts_midmu->Sumw2();
	h3d_recosliceid_allevts_cuts_mideg->Sumw2();
	h3d_recosliceid_allevts_cuts_midother->Sumw2();


/*
	h2d_truesliceid_allevts_cuts->Sumw2();
	h2d_truesliceid_allevts_cuts_inel->Sumw2();
	h2d_truesliceid_allevts_cuts_el->Sumw2();
	h2d_truesliceid_allevts_cuts_midcosmic->Sumw2();
	h2d_truesliceid_allevts_cuts_midpi->Sumw2();
	h2d_truesliceid_allevts_cuts_midp->Sumw2();
	h2d_truesliceid_allevts_cuts_midmu->Sumw2();
	h2d_truesliceid_allevts_cuts_mideg->Sumw2();
	h2d_truesliceid_allevts_cuts_midother->Sumw2();
*/


	//h3d_truesliceid_allevts_cuts->Sumw2();
	//h3d_truesliceid_allevts_cuts_inel->Sumw2();
	//h3d_truesliceid_allevts_cuts_el->Sumw2();
	//h3d_truesliceid_allevts_cuts_midp->Sumw2();
	//h3d_truesliceid_allevts_cuts_midcosmic->Sumw2();
	//h3d_truesliceid_allevts_cuts_midpi->Sumw2();
	//h3d_truesliceid_allevts_cuts_midmu->Sumw2();
	//h3d_truesliceid_allevts_cuts_mideg->Sumw2();
	//h3d_truesliceid_allevts_cuts_midother->Sumw2();

/*
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
*/

/*
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

*/

} //BookHistograms


void SaveHistograms() { //SaveHistograms
	outputFile->cd();
	outputFile->Write();
	//h_truesliceid_uf->Write("h_truesliceid_uf");
	//h_truesliceid_inelastic_uf->Write("h_truesliceid_inelastic_uf");

} //SaveHistograms

