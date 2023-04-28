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
TH1D *h_truesliceid_all;
TH1D *h_true_st_sliceid_all;
TH1D *h_truesliceid_inelastic_all;

TH1D *h_truesliceid_cuts;
TH1D *h_true_st_sliceid_cuts;
TH1D *h_truesliceid_inelastic_cuts;

//TH1D *h_truesliceid_uf;
//TH1D *h_truesliceid_inelastic_uf;

//reco inc
/*
TH1D *h_recosliceid_allevts_cuts;
TH1D *h_recosliceid_allevts_cuts_inel;
TH1D *h_recosliceid_allevts_cuts_el;
TH1D *h_recosliceid_allevts_cuts_midcosmic;
TH1D *h_recosliceid_allevts_cuts_midpi;
TH1D *h_recosliceid_allevts_cuts_midp;
TH1D *h_recosliceid_allevts_cuts_midmu;
TH1D *h_recosliceid_allevts_cuts_mideg;
TH1D *h_recosliceid_allevts_cuts_midother;
*/

//true inc
/*
TH1D *h_truesliceid_allevts_cuts;
TH1D *h_truesliceid_allevts_cuts_inel;
TH1D *h_truesliceid_allevts_cuts_el;
TH1D *h_truesliceid_allevts_cuts_midcosmic;
TH1D *h_truesliceid_allevts_cuts_midpi;
TH1D *h_truesliceid_allevts_cuts_midp;
TH1D *h_truesliceid_allevts_cuts_midmu;
TH1D *h_truesliceid_allevts_cuts_mideg;
TH1D *h_truesliceid_allevts_cuts_midother;
*/

//reco st inc
/*
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
*/

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

//reco int 2D
TH2D *h2d_recosliceid_recoinelastic_cuts;
TH2D *h2d_recosliceid_recoinelastic_cuts_el;
TH2D *h2d_recosliceid_recoinelastic_cuts_midp;
TH2D *h2d_recosliceid_recoinelastic_cuts_inel;
TH2D *h2d_recosliceid_recoinelastic_cuts_midcosmic;
TH2D *h2d_recosliceid_recoinelastic_cuts_midpi;
TH2D *h2d_recosliceid_recoinelastic_cuts_midmu;
TH2D *h2d_recosliceid_recoinelastic_cuts_mideg;
TH2D *h2d_recosliceid_recoinelastic_cuts_midother;

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

//reco inc_2D (st, end)
TH2D *h2d_recosliceid_allevts_cuts;
TH2D *h2d_recosliceid_allevts_cuts_el;
TH2D *h2d_recosliceid_allevts_cuts_midp;
TH2D *h2d_recosliceid_allevts_cuts_inel;
TH2D *h2d_recosliceid_allevts_cuts_midcosmic;
TH2D *h2d_recosliceid_allevts_cuts_midpi;
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

//double true_incidents[nthinslices+2];
//double true_st_incidents[nthinslices+2];
//double true_interactions[nthinslices+2];

//truth inc/int
//TH1D *h_true_incidents;
//TH1D *h_true_st_incidents;
//TH1D *h_true_interactions;

void BookHistograms() { //BookHistograms
	outputFile = TFile::Open(fOutputFileName.c_str(), "recreate"); //open output file

	//std::cout<<"Debug pass [Book1]..."<<std::endl;

	//XS histograms -------------------------------------------------------------------------------------------------------------------------------------------------------//
	h_truesliceid_all = new TH1D("h_truesliceid_all","h_truesliceid_all;True SliceID", p::true_nbins, p::true_bins);
	h_true_st_sliceid_all = new TH1D("h_true_st_sliceid_all","h_true_st_sliceid_all;True Start SliceID", p::true_nbins, p::true_bins);
	h_truesliceid_inelastic_all = new TH1D("h_truesliceid_inelastic_all","h_truesliceid_inelastic_all;True SliceID", p::true_nbins, p::true_bins);


	h_truesliceid_cuts = new TH1D("h_truesliceid_cuts","h_truesliceid_cuts;True SliceID", p::true_nbins, p::true_bins);
	h_true_st_sliceid_cuts = new TH1D("h_true_st_sliceid_cuts","h_true_st_sliceid_cuts;True Start SliceID", p::true_nbins, p::true_bins);
	h_truesliceid_inelastic_cuts = new TH1D("h_truesliceid_inelastic_cuts","h_truesliceid_inelastic_cuts;True SliceID", p::true_nbins, p::true_bins);

/*
	h_truesliceid_all = new TH1D("h_truesliceid_all","h_truesliceid_all;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_true_st_sliceid_all = new TH1D("h_true_st_sliceid_all","h_true_st_sliceid_all;True Start SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_inelastic_all = new TH1D("h_truesliceid_inelastic_all","h_truesliceid_inelastic_all;True SliceID", nthinslices + 2, -1, nthinslices + 1);


	h_truesliceid_cuts = new TH1D("h_truesliceid_cuts","h_truesliceid_cuts;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_true_st_sliceid_cuts = new TH1D("h_true_st_sliceid_cuts","h_true_st_sliceid_cuts;True Start SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_inelastic_cuts = new TH1D("h_truesliceid_inelastic_cuts","h_truesliceid_inelastic_cuts;True SliceID", nthinslices + 2, -1, nthinslices + 1);
*/

/*
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
*/

/*
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
*/

/*
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

	h2d_recosliceid_recoinelastic_cuts = new TH2D("h2d_recosliceid_recoinelastic_cuts","h2d_recosliceid_recoinelastic_cuts; Reco SliceID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_recoinelastic_cuts_el = new TH2D("h2d_recosliceid_recoinelastic_cuts_el","h2d_recosliceid_recoinelastic_cuts_el;Reco SliceID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_recoinelastic_cuts_midp = new TH2D("h2d_recosliceid_recoinelastic_cuts_midp","h2d_recosliceid_recoinelastic_cuts_midp;Reco SliceID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_recoinelastic_cuts_inel = new TH2D("h2d_recosliceid_recoinelastic_cuts_inel","h2d_recosliceid_recoinelastic_cuts_inel;Reco SliceID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_recoinelastic_cuts_midcosmic = new TH2D("h2d_recosliceid_recoinelastic_cuts_midcosmic","h2d_recosliceid_recoinelastic_cuts_midcosmic;Reco SliceID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_recoinelastic_cuts_midpi = new TH2D("h2d_recosliceid_recoinelastic_cuts_midpi","h2d_recosliceid_recoinelastic_cuts_midpi;Reco SliceID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_recoinelastic_cuts_midmu = new TH2D("h2d_recosliceid_recoinelastic_cuts_midmu","h2d_recosliceid_recoinelastic_cuts_midmu;Reco SliceID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_recoinelastic_cuts_mideg = new TH2D("h2d_recosliceid_recoinelastic_cuts_mideg","h2d_recosliceid_recoinelastic_cuts_mideg;Reco SliceID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_recoinelastic_cuts_midother = new TH2D("h2d_recosliceid_recoinelastic_cuts_midother","h2d_recosliceid_recoinelastic_cuts_midother;Reco SliceID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);

	h_truesliceid_recoinelastic_cuts = new TH1D("h_truesliceid_recoinelastic_cuts","h_truesliceid_recoinelastic_cuts;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_recoinelastic_cuts_inel = new TH1D("h_truesliceid_recoinelastic_cuts_inel","h_truesliceid_recoinelastic_cuts_inel;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_recoinelastic_cuts_el = new TH1D("h_truesliceid_recoinelastic_cuts_el","h_truesliceid_recoinelastic_cuts_el;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_recoinelastic_cuts_midcosmic = new TH1D("h_truesliceid_recoinelastic_cuts_midcosmic","h_truesliceid_recoinelastic_cuts_midcosmic;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_recoinelastic_cuts_midpi = new TH1D("h_truesliceid_recoinelastic_cuts_midpi","h_truesliceid_recoinelastic_cuts_midpi;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_recoinelastic_cuts_midp = new TH1D("h_truesliceid_recoinelastic_cuts_midp","h_truesliceid_recoinelastic_cuts_midp;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_recoinelastic_cuts_midmu = new TH1D("h_truesliceid_recoinelastic_cuts_midmu","h_truesliceid_recoinelastic_cuts_midmu;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_recoinelastic_cuts_mideg = new TH1D("h_truesliceid_recoinelastic_cuts_mideg","h_truesliceid_recoinelastic_cuts_mideg;True SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_recoinelastic_cuts_midother = new TH1D("h_truesliceid_recoinelastic_cuts_midother","h_truesliceid_recoinelastic_cuts_midother;True SliceID", nthinslices + 2, -1, nthinslices + 1);

	h2d_recosliceid_allevts_cuts = new TH2D("h2d_recosliceid_allevts_cuts","h2d_recosliceid_allevts_cuts; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_allevts_cuts_el = new TH2D("h2d_recosliceid_allevts_cuts_el","h2d_recosliceid_allevts_cuts_el; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_allevts_cuts_midp = new TH2D("h2d_recosliceid_allevts_cuts_midp","h2d_recosliceid_allevts_cuts_midp; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_allevts_cuts_inel = new TH2D("h2d_recosliceid_allevts_cuts_inel","h2d_recosliceid_allevts_cuts_inel; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_allevts_cuts_midcosmic = new TH2D("h2d_recosliceid_allevts_cuts_midcosmic","h2d_recosliceid_allevts_cuts_midcosmic; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_allevts_cuts_midpi = new TH2D("h2d_recosliceid_allevts_cuts_midpi","h2d_recosliceid_allevts_cuts_midpi; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_allevts_cuts_midmu = new TH2D("h2d_recosliceid_allevts_cuts_midmu","h2d_recosliceid_allevts_cuts_midmu; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_allevts_cuts_mideg = new TH2D("h2d_recosliceid_allevts_cuts_mideg","h2d_recosliceid_allevts_cuts_mideg; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);
	h2d_recosliceid_allevts_cuts_midother = new TH2D("h2d_recosliceid_allevts_cuts_midother","h2d_recosliceid_allevts_cuts_midother; StartID; EndID", nthinslices + 2, -1, nthinslices + 1, nthinslices + 2, -1, nthinslices + 1);

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

	h_recosliceid_cuts = new TH1D("h_recosliceid_cuts","h_recosliceid_cuts;Reco SliceID", p::reco_nbins, p::reco_bins);
	h_reco_st_sliceid_cuts = new TH1D("h_reco_st_sliceid_cuts","h_reco_st_sliceid_cuts;Reco Start SliceID", p::reco_nbins, p::reco_bins);
	h_recosliceid_inelastic_cuts = new TH1D("h_recosliceid_inelastic_cuts","h_recosliceid_inelastic_cuts;Reco SliceID", p::reco_nbins, p::reco_bins);

	h_recosliceid_recoinelastic_cuts = new TH1D("h_recosliceid_recoinelastic_cuts","h_recosliceid_recoinelastic_cuts;Reco SliceID", p::reco_nbins, p::reco_bins);
	h_recosliceid_recoinelastic_cuts_inel = new TH1D("h_recosliceid_recoinelastic_cuts_inel","h_recosliceid_recoinelastic_cuts_inel;Reco SliceID", p::reco_nbins, p::reco_bins);
	h_recosliceid_recoinelastic_cuts_el = new TH1D("h_recosliceid_recoinelastic_cuts_el","h_recosliceid_recoinelastic_cuts_el;Reco SliceID", p::reco_nbins, p::reco_bins);
	h_recosliceid_recoinelastic_cuts_midcosmic = new TH1D("h_recosliceid_recoinelastic_cuts_midcosmic","h_recosliceid_recoinelastic_cuts_midcosmic;Reco SliceID", p::reco_nbins, p::reco_bins);
	h_recosliceid_recoinelastic_cuts_midpi = new TH1D("h_recosliceid_recoinelastic_cuts_midpi","h_recosliceid_recoinelastic_cuts_midpi;Reco SliceID", p::reco_nbins, p::reco_bins);
	h_recosliceid_recoinelastic_cuts_midp = new TH1D("h_recosliceid_recoinelastic_cuts_midp","h_recosliceid_recoinelastic_cuts_midp;Reco SliceID", p::reco_nbins, p::reco_bins);
	h_recosliceid_recoinelastic_cuts_midmu = new TH1D("h_recosliceid_recoinelastic_cuts_midmu","h_recosliceid_recoinelastic_cuts_midmu;Reco SliceID", p::reco_nbins, p::reco_bins);
	h_recosliceid_recoinelastic_cuts_mideg = new TH1D("h_recosliceid_recoinelastic_cuts_mideg","h_recosliceid_recoinelastic_cuts_mideg;Reco SliceID", p::reco_nbins, p::reco_bins);
	h_recosliceid_recoinelastic_cuts_midother = new TH1D("h_recosliceid_recoinelastic_cuts_midother","h_recosliceid_recoinelastic_cuts_midother;Reco SliceID", p::reco_nbins, p::reco_bins);

	h2d_recosliceid_recoinelastic_cuts = new TH2D("h2d_recosliceid_recoinelastic_cuts","h2d_recosliceid_recoinelastic_cuts; Reco SliceID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	h2d_recosliceid_recoinelastic_cuts_el = new TH2D("h2d_recosliceid_recoinelastic_cuts_el","h2d_recosliceid_recoinelastic_cuts_el;Reco SliceID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	h2d_recosliceid_recoinelastic_cuts_midp = new TH2D("h2d_recosliceid_recoinelastic_cuts_midp","h2d_recosliceid_recoinelastic_cuts_midp;Reco SliceID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	h2d_recosliceid_recoinelastic_cuts_inel = new TH2D("h2d_recosliceid_recoinelastic_cuts_inel","h2d_recosliceid_recoinelastic_cuts_inel;Reco SliceID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	h2d_recosliceid_recoinelastic_cuts_midcosmic = new TH2D("h2d_recosliceid_recoinelastic_cuts_midcosmic","h2d_recosliceid_recoinelastic_cuts_midcosmic;Reco SliceID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	h2d_recosliceid_recoinelastic_cuts_midpi = new TH2D("h2d_recosliceid_recoinelastic_cuts_midpi","h2d_recosliceid_recoinelastic_cuts_midpi;Reco SliceID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	h2d_recosliceid_recoinelastic_cuts_midmu = new TH2D("h2d_recosliceid_recoinelastic_cuts_midmu","h2d_recosliceid_recoinelastic_cuts_midmu;Reco SliceID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	h2d_recosliceid_recoinelastic_cuts_mideg = new TH2D("h2d_recosliceid_recoinelastic_cuts_mideg","h2d_recosliceid_recoinelastic_cuts_mideg;Reco SliceID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	h2d_recosliceid_recoinelastic_cuts_midother = new TH2D("h2d_recosliceid_recoinelastic_cuts_midother","h2d_recosliceid_recoinelastic_cuts_midother;Reco SliceID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);

	h_truesliceid_recoinelastic_cuts = new TH1D("h_truesliceid_recoinelastic_cuts","h_truesliceid_recoinelastic_cuts;True SliceID", p::true_nbins, p::true_bins);
	h_truesliceid_recoinelastic_cuts_inel = new TH1D("h_truesliceid_recoinelastic_cuts_inel","h_truesliceid_recoinelastic_cuts_inel;True SliceID", p::true_nbins, p::true_bins);
	h_truesliceid_recoinelastic_cuts_el = new TH1D("h_truesliceid_recoinelastic_cuts_el","h_truesliceid_recoinelastic_cuts_el;True SliceID", p::true_nbins, p::true_bins);
	h_truesliceid_recoinelastic_cuts_midcosmic = new TH1D("h_truesliceid_recoinelastic_cuts_midcosmic","h_truesliceid_recoinelastic_cuts_midcosmic;True SliceID", p::true_nbins, p::true_bins);
	h_truesliceid_recoinelastic_cuts_midpi = new TH1D("h_truesliceid_recoinelastic_cuts_midpi","h_truesliceid_recoinelastic_cuts_midpi;True SliceID", p::true_nbins, p::true_bins);
	h_truesliceid_recoinelastic_cuts_midp = new TH1D("h_truesliceid_recoinelastic_cuts_midp","h_truesliceid_recoinelastic_cuts_midp;True SliceID", p::true_nbins, p::true_bins);
	h_truesliceid_recoinelastic_cuts_midmu = new TH1D("h_truesliceid_recoinelastic_cuts_midmu","h_truesliceid_recoinelastic_cuts_midmu;True SliceID", p::true_nbins, p::true_bins);
	h_truesliceid_recoinelastic_cuts_mideg = new TH1D("h_truesliceid_recoinelastic_cuts_mideg","h_truesliceid_recoinelastic_cuts_mideg;True SliceID", p::true_nbins, p::true_bins);
	h_truesliceid_recoinelastic_cuts_midother = new TH1D("h_truesliceid_recoinelastic_cuts_midother","h_truesliceid_recoinelastic_cuts_midother;True SliceID", p::true_nbins, p::true_bins);

	h2d_recosliceid_allevts_cuts = new TH2D("h2d_recosliceid_allevts_cuts","h2d_recosliceid_allevts_cuts; StartID; EndID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	h2d_recosliceid_allevts_cuts_el = new TH2D("h2d_recosliceid_allevts_cuts_el","h2d_recosliceid_allevts_cuts_el; StartID; EndID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	h2d_recosliceid_allevts_cuts_midp = new TH2D("h2d_recosliceid_allevts_cuts_midp","h2d_recosliceid_allevts_cuts_midp; StartID; EndID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	h2d_recosliceid_allevts_cuts_inel = new TH2D("h2d_recosliceid_allevts_cuts_inel","h2d_recosliceid_allevts_cuts_inel; StartID; EndID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	h2d_recosliceid_allevts_cuts_midcosmic = new TH2D("h2d_recosliceid_allevts_cuts_midcosmic","h2d_recosliceid_allevts_cuts_midcosmic; StartID; EndID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	h2d_recosliceid_allevts_cuts_midpi = new TH2D("h2d_recosliceid_allevts_cuts_midpi","h2d_recosliceid_allevts_cuts_midpi; StartID; EndID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	h2d_recosliceid_allevts_cuts_midmu = new TH2D("h2d_recosliceid_allevts_cuts_midmu","h2d_recosliceid_allevts_cuts_midmu; StartID; EndID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	h2d_recosliceid_allevts_cuts_mideg = new TH2D("h2d_recosliceid_allevts_cuts_mideg","h2d_recosliceid_allevts_cuts_mideg; StartID; EndID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);
	h2d_recosliceid_allevts_cuts_midother = new TH2D("h2d_recosliceid_allevts_cuts_midother","h2d_recosliceid_allevts_cuts_midother; StartID; EndID", p::reco_nbins, p::reco_bins, p::reco_nbins, p::reco_bins);

	h2d_truesliceid_allevts_cuts = new TH2D("h2d_truesliceid_allevts_cuts","h2d_truesliceid_allevts_cuts; StartID; EndID", p::true_nbins, p::true_bins, p::true_nbins, p::true_bins);
	h2d_truesliceid_allevts_cuts_inel = new TH2D("h2d_truesliceid_allevts_cuts_inel","h2d_truesliceid_allevts_cuts_inel; StartID; EndID", p::true_nbins, p::true_bins, p::true_nbins, p::true_bins);
	h2d_truesliceid_allevts_cuts_el = new TH2D("h2d_truesliceid_allevts_cuts_el","h2d_truesliceid_allevts_cuts_el; StartID; EndID", p::true_nbins, p::true_bins, p::true_nbins, p::true_bins);
	h2d_truesliceid_allevts_cuts_midcosmic = new TH2D("h2d_truesliceid_allevts_cuts_midcosmic","h2d_truesliceid_allevts_cuts_midcosmic; StartID; EndID", p::true_nbins, p::true_bins, p::true_nbins, p::true_bins);
	h2d_truesliceid_allevts_cuts_midpi = new TH2D("h2d_truesliceid_allevts_cuts_midpi","h2d_truesliceid_allevts_cuts_midpi; StartID; EndID", p::true_nbins, p::true_bins, p::true_nbins, p::true_bins);
	h2d_truesliceid_allevts_cuts_midp = new TH2D("h2d_truesliceid_allevts_cuts_midp","h2d_truesliceid_allevts_cuts_midp; StartID; EndID", p::true_nbins, p::true_bins, p::true_nbins, p::true_bins);
	h2d_truesliceid_allevts_cuts_midmu = new TH2D("h2d_truesliceid_allevts_cuts_midmu","h2d_truesliceid_allevts_cuts_midmu; StartID; EndID", p::true_nbins, p::true_bins, p::true_nbins, p::true_bins);
	h2d_truesliceid_allevts_cuts_mideg = new TH2D("h2d_truesliceid_allevts_cuts_mideg","h2d_truesliceid_allevts_cuts_mideg; StartID; EndID", p::true_nbins, p::true_bins, p::true_nbins, p::true_bins);
	h2d_truesliceid_allevts_cuts_midother = new TH2D("h2d_truesliceid_allevts_cuts_midother","h2d_truesliceid_allevts_cuts_midother; StartID; EndID", p::true_nbins, p::true_bins, p::true_nbins, p::true_bins);


	h_truesliceid_all->Sumw2();
	h_true_st_sliceid_all->Sumw2();

	h_truesliceid_cuts->Sumw2();
	h_true_st_sliceid_cuts->Sumw2();
	h_truesliceid_inelastic_all->Sumw2();
	h_truesliceid_inelastic_cuts->Sumw2();

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

	h2d_recosliceid_allevts_cuts->Sumw2();
	h2d_recosliceid_allevts_cuts_el->Sumw2();
	h2d_recosliceid_allevts_cuts_midp->Sumw2();
	h2d_recosliceid_allevts_cuts_inel->Sumw2();
	h2d_recosliceid_allevts_cuts_midcosmic->Sumw2();
	h2d_recosliceid_allevts_cuts_midpi->Sumw2();
	h2d_recosliceid_allevts_cuts_midmu->Sumw2();
	h2d_recosliceid_allevts_cuts_mideg->Sumw2();
	h2d_recosliceid_allevts_cuts_midother->Sumw2();

	h2d_truesliceid_allevts_cuts->Sumw2();
	h2d_truesliceid_allevts_cuts_inel->Sumw2();
	h2d_truesliceid_allevts_cuts_el->Sumw2();
	h2d_truesliceid_allevts_cuts_midcosmic->Sumw2();
	h2d_truesliceid_allevts_cuts_midpi->Sumw2();
	h2d_truesliceid_allevts_cuts_midp->Sumw2();
	h2d_truesliceid_allevts_cuts_midmu->Sumw2();
	h2d_truesliceid_allevts_cuts_mideg->Sumw2();
	h2d_truesliceid_allevts_cuts_midother->Sumw2();

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

	h2d_recosliceid_recoinelastic_cuts->Sumw2();
	h2d_recosliceid_recoinelastic_cuts_el->Sumw2();
	h2d_recosliceid_recoinelastic_cuts_midp->Sumw2();
	h2d_recosliceid_recoinelastic_cuts_inel->Sumw2();
	h2d_recosliceid_recoinelastic_cuts_midcosmic->Sumw2();
	h2d_recosliceid_recoinelastic_cuts_midpi->Sumw2();
	h2d_recosliceid_recoinelastic_cuts_midmu->Sumw2();
	h2d_recosliceid_recoinelastic_cuts_mideg->Sumw2();
	h2d_recosliceid_recoinelastic_cuts_midother->Sumw2();

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
	//for (int i = 0; i<nthinslices+2; ++i){
		//true_interactions[i] = 0;
		//true_incidents[i] = 0;
		//true_st_incidents[i] = 0;
	//}

	//KE calc using reco stopping protons ----------------------------------------------------//
	//int n_ke=160;
	//double ke_min=-800;
	//double ke_max=800;

	//dke -----------------------------------------------------------------------------------------------------//
	//int n_dke=160;
	//float dke_st=-800;
	//float dke_end=800; 

	//ntrklen ------------------------------------------------------------------------------------//
	//int n_ntrklen=61;
	//double st_ntrklen=-0.02;
	//double ed_ntrklen=1.2;
	//bq cut

	//range vs ke
        //int dke=325;
        //float ke_st=-50;
        //float ke_end=600;

	//truth inc/int
	//h_true_incidents=new TH1D("h_true_incidents","true_incidents", nthinslices, 0, nthinslices-1);
	//h_true_st_incidents=new TH1D("h_true_st_incidents","true_st_incidents", nthinslices, 0, nthinslices-1);
	//h_true_interactions=new TH1D("h_true_interactions","true_interactions", nthinslices, 0, nthinslices-1);
	//h_true_incidents=new TH1D("h_true_incidents","true_incidents", nthinslices+2, -1, nthinslices+1);
	//h_true_st_incidents=new TH1D("h_true_st_incidents","true_st_incidents", nthinslices+2, -1, nthinslices+1);
	//h_true_interactions=new TH1D("h_true_interactions","true_interactions", nthinslices+2, -1, nthinslices+1);
	//h_true_incidents->Sumw2();
	//h_true_st_incidents->Sumw2();
	//h_true_interactions->Sumw2();

	//std::cout<<"Debug pass [Book2]..."<<std::endl;

} //BookHistograms


void SaveHistograms() { //SaveHistograms
	outputFile->cd();
	outputFile->Write();
	//h_truesliceid_uf->Write("h_truesliceid_uf");
	//h_truesliceid_inelastic_uf->Write("h_truesliceid_inelastic_uf");

} //SaveHistograms

/*
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

	//TH1D *hinc = (TH1D*)h_recosliceid_allevts_cuts->Clone("hinc");
	////TH1D *hint = (TH1D*)h_recosliceid_allevts_cuts->Clone("hint");
	//TH1D *hint = (TH1D*)h_recosliceid_recoinelastic_cuts->Clone("hint");
	//hinc->Multiply(uf.pur_Inc);
	//hint->Multiply(uf.pur_Int);

	//  RooUnfoldBayes   unfold_Inc (&uf.response_SliceID_Inc, uf.pur_num_Inc, 4);
	//  RooUnfoldBayes   unfold_Int (&uf.response_SliceID_Int, uf.pur_num_Int, 4);

	//RooUnfoldBayes   unfold_Inc (&uf.response_SliceID_Inc, hinc, 4);
	//RooUnfoldBayes   unfold_Int (&uf.response_SliceID_Int, hint, 4);

	//RooUnfoldBayes   unfold_Inc (&uf.response_SliceID_Inc, hinc, 12);
	//RooUnfoldBayes   unfold_Int (&uf.response_SliceID_Int, hint, 12);

	//RooUnfoldSvd     unfold_Inc (&uf.response_SliceID_Inc, hinc, 20);   // OR
	//RooUnfoldSvd     unfold_Int (&uf.response_SliceID_Int, hint, 20);   // OR

	//  RooUnfoldSvd     unfold_Inc (&uf.response_SliceID_Inc, uf.pur_num_Inc, 20);   // OR
	//  RooUnfoldSvd     unfold_Int (&uf.response_SliceID_Int, uf.pur_num_Int, 20);   // OR

	//h_truesliceid_uf = (TH1D*) unfold_Inc.Hreco();
	//h_truesliceid_inelastic_uf = (TH1D*) unfold_Int.Hreco();

} //CalcXS
*/

