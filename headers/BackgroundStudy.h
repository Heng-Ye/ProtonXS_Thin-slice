#include "TGraphErrors.h"
#include "TVector3.h"
//#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "util.h"
#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
//#include "./Unfold.h"

//Basic config. -----------------------------------------------------//
std::string fOutputFileName;
TFile *outputFile;
void SetOutputFileName(std::string name){fOutputFileName = name;};
void BookHistograms();
void SaveHistograms();

//reco sliceID histograms
//misID:p-rich sample
TH1D *h_recosliceid_cosLE09;
TH1D *h_recosliceid_cosLE09_inel;
TH1D *h_recosliceid_cosLE09_el;
TH1D *h_recosliceid_cosLE09_midcosmic;
TH1D *h_recosliceid_cosLE09_midpi;
TH1D *h_recosliceid_cosLE09_midp;
TH1D *h_recosliceid_cosLE09_midmu;
TH1D *h_recosliceid_cosLE09_mideg;
TH1D *h_recosliceid_cosLE09_midother;

TH1D *h_recosliceid_cosGT09;
TH1D *h_recosliceid_cosGT09_inel;
TH1D *h_recosliceid_cosGT09_el;
TH1D *h_recosliceid_cosGT09_midcosmic;
TH1D *h_recosliceid_cosGT09_midpi;
TH1D *h_recosliceid_cosGT09_midp;
TH1D *h_recosliceid_cosGT09_midmu;
TH1D *h_recosliceid_cosGT09_mideg;
TH1D *h_recosliceid_cosGT09_midother;

void BookHistograms() { //BookHistograms
	outputFile = TFile::Open(fOutputFileName.c_str(), "recreate"); //open output file

	h_recosliceid_cosLE09 = new TH1D("h_recosliceid_cosLE09","h_recosliceid_cosLE09;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE09_inel = new TH1D("h_recosliceid_cosLE09_inel","h_recosliceid_cosLE09_inel;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE09_el = new TH1D("h_recosliceid_cosLE09_el","h_recosliceid_cosLE09_el;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE09_midcosmic = new TH1D("h_recosliceid_cosLE09_midcosmic","h_recosliceid_cosLE09_midcosmic;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE09_midpi = new TH1D("h_recosliceid_cosLE09_midpi","h_recosliceid_cosLE09_midpi;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE09_midp = new TH1D("h_recosliceid_cosLE09_midp","h_recosliceid_cosLE09_midp;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE09_midmu = new TH1D("h_recosliceid_cosLE09_midmu","h_recosliceid_cosLE09_midmu;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE09_mideg = new TH1D("h_recosliceid_cosLE09_mideg","h_recosliceid_cosLE09_mideg;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE09_midother = new TH1D("h_recosliceid_cosLE09_midother","h_recosliceid_cosLE09_midother;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	
	h_recosliceid_cosLE09->Sumw2();
	h_recosliceid_cosLE09_inel->Sumw2();
	h_recosliceid_cosLE09_el->Sumw2();
	h_recosliceid_cosLE09_midcosmic->Sumw2();
	h_recosliceid_cosLE09_midpi->Sumw2();
	h_recosliceid_cosLE09_midp->Sumw2();
	h_recosliceid_cosLE09_midmu->Sumw2();
	h_recosliceid_cosLE09_mideg->Sumw2();
	h_recosliceid_cosLE09_midother->Sumw2();

	h_recosliceid_cosGT09 = new TH1D("h_recosliceid_cosGT09","h_recosliceid_cosGT09;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosGT09_inel = new TH1D("h_recosliceid_cosGT09_inel","h_recosliceid_cosGT09_inel;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosGT09_el = new TH1D("h_recosliceid_cosGT09_el","h_recosliceid_cosGT09_el;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosGT09_midcosmic = new TH1D("h_recosliceid_cosGT09_midcosmic","h_recosliceid_cosGT09_midcosmic;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosGT09_midpi = new TH1D("h_recosliceid_cosGT09_midpi","h_recosliceid_cosGT09_midpi;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosGT09_midp = new TH1D("h_recosliceid_cosGT09_midp","h_recosliceid_cosGT09_midp;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosGT09_midmu = new TH1D("h_recosliceid_cosGT09_midmu","h_recosliceid_cosGT09_midmu;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosGT09_mideg = new TH1D("h_recosliceid_cosGT09_mideg","h_recosliceid_cosGT09_mideg;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosGT09_midother = new TH1D("h_recosliceid_cosGT09_midother","h_recosliceid_cosGT09_midother;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);

	h_recosliceid_cosGT09->Sumw2();
	h_recosliceid_cosGT09_inel->Sumw2();
	h_recosliceid_cosGT09_el->Sumw2();
	h_recosliceid_cosGT09_midcosmic->Sumw2();
	h_recosliceid_cosGT09_midpi->Sumw2();
	h_recosliceid_cosGT09_midp->Sumw2();
	h_recosliceid_cosGT09_midmu->Sumw2();
	h_recosliceid_cosGT09_mideg->Sumw2();
	h_recosliceid_cosGT09_midother->Sumw2();


} //BookHistograms


void SaveHistograms() { //SaveHistograms
	outputFile->cd();
	outputFile->Write();

} //SaveHistograms


