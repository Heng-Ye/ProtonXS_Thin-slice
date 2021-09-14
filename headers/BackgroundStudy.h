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

TH1D *h_truesliceid_cosLE09;
TH1D *h_truesliceid_cosLE09_inel;
TH1D *h_truesliceid_cosLE09_el;
TH1D *h_truesliceid_cosLE09_midcosmic;
TH1D *h_truesliceid_cosLE09_midpi;
TH1D *h_truesliceid_cosLE09_midp;
TH1D *h_truesliceid_cosLE09_midmu;
TH1D *h_truesliceid_cosLE09_mideg;
TH1D *h_truesliceid_cosLE09_midother;


TH1D *h_recosliceid_cosLE08;
TH1D *h_recosliceid_cosLE08_inel;
TH1D *h_recosliceid_cosLE08_el;
TH1D *h_recosliceid_cosLE08_midcosmic;
TH1D *h_recosliceid_cosLE08_midpi;
TH1D *h_recosliceid_cosLE08_midp;
TH1D *h_recosliceid_cosLE08_midmu;
TH1D *h_recosliceid_cosLE08_mideg;
TH1D *h_recosliceid_cosLE08_midother;


TH1D *h_recosliceid_cosLE07;
TH1D *h_recosliceid_cosLE07_inel;
TH1D *h_recosliceid_cosLE07_el;
TH1D *h_recosliceid_cosLE07_midcosmic;
TH1D *h_recosliceid_cosLE07_midpi;
TH1D *h_recosliceid_cosLE07_midp;
TH1D *h_recosliceid_cosLE07_midmu;
TH1D *h_recosliceid_cosLE07_mideg;
TH1D *h_recosliceid_cosLE07_midother;



TH1D *h_recosliceid_cosGT09;
TH1D *h_recosliceid_cosGT09_inel;
TH1D *h_recosliceid_cosGT09_el;
TH1D *h_recosliceid_cosGT09_midcosmic;
TH1D *h_recosliceid_cosGT09_midpi;
TH1D *h_recosliceid_cosGT09_midp;
TH1D *h_recosliceid_cosGT09_midmu;
TH1D *h_recosliceid_cosGT09_mideg;
TH1D *h_recosliceid_cosGT09_midother;


//cosTheta
TH1D *h_cosTheta_Pos[nthinslices+2];
TH1D *h_cosTheta_Pos_inel[nthinslices+2];
TH1D *h_cosTheta_Pos_el[nthinslices+2];
TH1D *h_cosTheta_Pos_midcosmic[nthinslices+2];
TH1D *h_cosTheta_Pos_midpi[nthinslices+2];
TH1D *h_cosTheta_Pos_midp[nthinslices+2];
TH1D *h_cosTheta_Pos_midmu[nthinslices+2];
TH1D *h_cosTheta_Pos_mideg[nthinslices+2];
TH1D *h_cosTheta_Pos_midother[nthinslices+2];

TH1D *h_cosTheta_Pos_all;
TH1D *h_cosTheta_Pos_all_inel;
TH1D *h_cosTheta_Pos_all_el;
TH1D *h_cosTheta_Pos_all_midcosmic;
TH1D *h_cosTheta_Pos_all_midpi;
TH1D *h_cosTheta_Pos_all_midp;
TH1D *h_cosTheta_Pos_all_midmu;
TH1D *h_cosTheta_Pos_all_mideg;
TH1D *h_cosTheta_Pos_all_midother;

TH1D *h_cosTheta_Pos_all_nosliceidOne;
TH1D *h_cosTheta_Pos_all_nosliceidOne_inel;
TH1D *h_cosTheta_Pos_all_nosliceidOne_el;
TH1D *h_cosTheta_Pos_all_nosliceidOne_midcosmic;
TH1D *h_cosTheta_Pos_all_nosliceidOne_midpi;
TH1D *h_cosTheta_Pos_all_nosliceidOne_midp;
TH1D *h_cosTheta_Pos_all_nosliceidOne_midmu;
TH1D *h_cosTheta_Pos_all_nosliceidOne_mideg;
TH1D *h_cosTheta_Pos_all_nosliceidOne_midother;


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


	h_truesliceid_cosLE09 = new TH1D("h_truesliceid_cosLE09","h_truesliceid_cosLE09;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_cosLE09_inel = new TH1D("h_truesliceid_cosLE09_inel","h_truesliceid_cosLE09_inel;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_cosLE09_el = new TH1D("h_truesliceid_cosLE09_el","h_truesliceid_cosLE09_el;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_cosLE09_midcosmic = new TH1D("h_truesliceid_cosLE09_midcosmic","h_truesliceid_cosLE09_midcosmic;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_cosLE09_midpi = new TH1D("h_truesliceid_cosLE09_midpi","h_truesliceid_cosLE09_midpi;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_cosLE09_midp = new TH1D("h_truesliceid_cosLE09_midp","h_truesliceid_cosLE09_midp;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_cosLE09_midmu = new TH1D("h_truesliceid_cosLE09_midmu","h_truesliceid_cosLE09_midmu;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_cosLE09_mideg = new TH1D("h_truesliceid_cosLE09_mideg","h_truesliceid_cosLE09_mideg;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_truesliceid_cosLE09_midother = new TH1D("h_truesliceid_cosLE09_midother","h_truesliceid_cosLE09_midother;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	
	h_truesliceid_cosLE09->Sumw2();
	h_truesliceid_cosLE09_inel->Sumw2();
	h_truesliceid_cosLE09_el->Sumw2();
	h_truesliceid_cosLE09_midcosmic->Sumw2();
	h_truesliceid_cosLE09_midpi->Sumw2();
	h_truesliceid_cosLE09_midp->Sumw2();
	h_truesliceid_cosLE09_midmu->Sumw2();
	h_truesliceid_cosLE09_mideg->Sumw2();
	h_truesliceid_cosLE09_midother->Sumw2();


	h_recosliceid_cosLE08 = new TH1D("h_recosliceid_cosLE08","h_recosliceid_cosLE08;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE08_inel = new TH1D("h_recosliceid_cosLE08_inel","h_recosliceid_cosLE08_inel;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE08_el = new TH1D("h_recosliceid_cosLE08_el","h_recosliceid_cosLE08_el;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE08_midcosmic = new TH1D("h_recosliceid_cosLE08_midcosmic","h_recosliceid_cosLE08_midcosmic;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE08_midpi = new TH1D("h_recosliceid_cosLE08_midpi","h_recosliceid_cosLE08_midpi;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE08_midp = new TH1D("h_recosliceid_cosLE08_midp","h_recosliceid_cosLE08_midp;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE08_midmu = new TH1D("h_recosliceid_cosLE08_midmu","h_recosliceid_cosLE08_midmu;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE08_mideg = new TH1D("h_recosliceid_cosLE08_mideg","h_recosliceid_cosLE08_mideg;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE08_midother = new TH1D("h_recosliceid_cosLE08_midother","h_recosliceid_cosLE08_midother;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	
	h_recosliceid_cosLE08->Sumw2();
	h_recosliceid_cosLE08_inel->Sumw2();
	h_recosliceid_cosLE08_el->Sumw2();
	h_recosliceid_cosLE08_midcosmic->Sumw2();
	h_recosliceid_cosLE08_midpi->Sumw2();
	h_recosliceid_cosLE08_midp->Sumw2();
	h_recosliceid_cosLE08_midmu->Sumw2();
	h_recosliceid_cosLE08_mideg->Sumw2();
	h_recosliceid_cosLE08_midother->Sumw2();




	h_recosliceid_cosLE07 = new TH1D("h_recosliceid_cosLE07","h_recosliceid_cosLE07;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE07_inel = new TH1D("h_recosliceid_cosLE07_inel","h_recosliceid_cosLE07_inel;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE07_el = new TH1D("h_recosliceid_cosLE07_el","h_recosliceid_cosLE07_el;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE07_midcosmic = new TH1D("h_recosliceid_cosLE07_midcosmic","h_recosliceid_cosLE07_midcosmic;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE07_midpi = new TH1D("h_recosliceid_cosLE07_midpi","h_recosliceid_cosLE07_midpi;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE07_midp = new TH1D("h_recosliceid_cosLE07_midp","h_recosliceid_cosLE07_midp;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE07_midmu = new TH1D("h_recosliceid_cosLE07_midmu","h_recosliceid_cosLE07_midmu;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE07_mideg = new TH1D("h_recosliceid_cosLE07_mideg","h_recosliceid_cosLE07_mideg;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	h_recosliceid_cosLE07_midother = new TH1D("h_recosliceid_cosLE07_midother","h_recosliceid_cosLE07_midother;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
	
	h_recosliceid_cosLE07->Sumw2();
	h_recosliceid_cosLE07_inel->Sumw2();
	h_recosliceid_cosLE07_el->Sumw2();
	h_recosliceid_cosLE07_midcosmic->Sumw2();
	h_recosliceid_cosLE07_midpi->Sumw2();
	h_recosliceid_cosLE07_midp->Sumw2();
	h_recosliceid_cosLE07_midmu->Sumw2();
	h_recosliceid_cosLE07_mideg->Sumw2();
	h_recosliceid_cosLE07_midother->Sumw2();



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

	//cosTheta dists.
        int n_cosine=100;
        double cosine_min=0;
        double cosine_max=1.0;

        int nn_cosine=2; //
	double cosbins[3];
	cosbins[0] = 0;
	cosbins[1] = 0.9;
	cosbins[2] = 1;
	for (int i=0; i<nthinslices+2; ++i) {
		h_cosTheta_Pos[i]=new TH1D(Form("h_cosTheta_Pos_recosliceid_%d", i-1), Form("Pos; Reco SliceID %d",i-1), nn_cosine, cosbins);
		h_cosTheta_Pos_inel[i]=new TH1D(Form("h_cosTheta_Pos_inel_recosliceid_%d", i-1), Form("Pos_inel; Reco SliceID %d",i-1), nn_cosine, cosbins);
		h_cosTheta_Pos_el[i]=new TH1D(Form("h_cosTheta_Pos_el_recosliceid_%d", i-1), Form("Pos_el; Reco SliceID %d",i-1), nn_cosine, cosbins);
		h_cosTheta_Pos_midcosmic[i]=new TH1D(Form("h_cosTheta_Pos_midcosmic_recosliceid_%d", i-1), Form("Pos_midcosmic; Reco SliceID %d",i-1), nn_cosine, cosbins);
		h_cosTheta_Pos_midpi[i]=new TH1D(Form("h_cosTheta_Pos_midpi_recosliceid_%d", i-1), Form("Pos_midpi; Reco SliceID %d",i-1), nn_cosine, cosbins);
		h_cosTheta_Pos_midp[i]=new TH1D(Form("h_cosTheta_Pos_midp_recosliceid_%d", i-1), Form("Pos_midp; Reco SliceID %d",i-1), nn_cosine, cosbins);
		h_cosTheta_Pos_midmu[i]=new TH1D(Form("h_cosTheta_Pos_midmu_recosliceid_%d", i-1), Form("Pos_midmu; Reco SliceID %d",i-1), nn_cosine, cosbins);
		h_cosTheta_Pos_mideg[i]=new TH1D(Form("h_cosTheta_Pos_mideg_recosliceid_%d", i-1), Form("Pos_mideg; Reco SliceID %d",i-1), nn_cosine, cosbins);
		h_cosTheta_Pos_midother[i]=new TH1D(Form("h_cosTheta_Pos_midother_recosliceid_%d", i-1), Form("Pos_midother; Reco SliceID %d",i-1), nn_cosine, cosbins);
		h_cosTheta_Pos[i]->Sumw2();
		h_cosTheta_Pos_inel[i]->Sumw2();
		h_cosTheta_Pos_el[i]->Sumw2();
		h_cosTheta_Pos_midcosmic[i]->Sumw2();
		h_cosTheta_Pos_midpi[i]->Sumw2();
		h_cosTheta_Pos_midp[i]->Sumw2();
		h_cosTheta_Pos_midmu[i]->Sumw2();
		h_cosTheta_Pos_mideg[i]->Sumw2();
		h_cosTheta_Pos_midother[i]->Sumw2();
	}
	h_cosTheta_Pos_all=new TH1D("h_cosTheta_Pos_all","h_cosTheta_Pos_all", n_cosine, cosine_min, cosine_max);
	h_cosTheta_Pos_all_inel=new TH1D("h_cosTheta_Pos_all_inel","h_cosTheta_Pos_all_inel", n_cosine, cosine_min, cosine_max);
	h_cosTheta_Pos_all_el=new TH1D("h_cosTheta_Pos_all_el","h_cosTheta_Pos_all_el", n_cosine, cosine_min, cosine_max);
	h_cosTheta_Pos_all_midcosmic=new TH1D("h_cosTheta_Pos_all_midcosmic","h_cosTheta_Pos_all_midcosmic", n_cosine, cosine_min, cosine_max);
	h_cosTheta_Pos_all_midpi=new TH1D("h_cosTheta_Pos_all_midpi","h_cosTheta_Pos_all_midpi", n_cosine, cosine_min, cosine_max);
	h_cosTheta_Pos_all_midp=new TH1D("h_cosTheta_Pos_all_midp","h_cosTheta_Pos_all_midp", n_cosine, cosine_min, cosine_max);
	h_cosTheta_Pos_all_midmu=new TH1D("h_cosTheta_Pos_all_midmu","h_cosTheta_Pos_all_midmu", n_cosine, cosine_min, cosine_max);
	h_cosTheta_Pos_all_mideg=new TH1D("h_cosTheta_Pos_all_mideg","h_cosTheta_Pos_all_mideg", n_cosine, cosine_min, cosine_max);
	h_cosTheta_Pos_all_midother=new TH1D("h_cosTheta_Pos_all_midother","h_cosTheta_Pos_all_midother", n_cosine, cosine_min, cosine_max);
	h_cosTheta_Pos_all->Sumw2();
	h_cosTheta_Pos_all_inel->Sumw2();
	h_cosTheta_Pos_all_el->Sumw2();
	h_cosTheta_Pos_all_midcosmic->Sumw2();
	h_cosTheta_Pos_all_midpi->Sumw2();
	h_cosTheta_Pos_all_midp->Sumw2();
	h_cosTheta_Pos_all_midmu->Sumw2();
	h_cosTheta_Pos_all_mideg->Sumw2();
	h_cosTheta_Pos_all_midother->Sumw2();


	h_cosTheta_Pos_all_nosliceidOne=new TH1D("h_cosTheta_Pos_all_nosliceidOne","h_cosTheta_Pos_all_nosliceidOne",n_cosine, cosine_min, cosine_max);
	h_cosTheta_Pos_all_nosliceidOne_inel=new TH1D("h_cosTheta_Pos_all_nosliceidOne_inel","h_cosTheta_Pos_all_nosliceidOne_inel", n_cosine, cosine_min, cosine_max);
	h_cosTheta_Pos_all_nosliceidOne_el=new TH1D("h_cosTheta_Pos_all_nosliceidOne_el","h_cosTheta_Pos_all_nosliceidOne_el", n_cosine, cosine_min, cosine_max);
	h_cosTheta_Pos_all_nosliceidOne_midcosmic=new TH1D("h_cosTheta_Pos_all_nosliceidOne_midcosmic","h_cosTheta_Pos_all_nosliceidOne_midcosmic", n_cosine, cosine_min, cosine_max);
	h_cosTheta_Pos_all_nosliceidOne_midpi=new TH1D("h_cosTheta_Pos_all_nosliceidOne_midpi","h_cosTheta_Pos_all_nosliceidOne_midpi", n_cosine, cosine_min, cosine_max);
	h_cosTheta_Pos_all_nosliceidOne_midp=new TH1D("h_cosTheta_Pos_all_nosliceidOne_midp","h_cosTheta_Pos_all_nosliceidOne_midp", n_cosine, cosine_min, cosine_max);
	h_cosTheta_Pos_all_nosliceidOne_midmu=new TH1D("h_cosTheta_Pos_all_nosliceidOne_midmu","h_cosTheta_Pos_all_nosliceidOne_midmu", n_cosine, cosine_min, cosine_max);
	h_cosTheta_Pos_all_nosliceidOne_mideg=new TH1D("h_cosTheta_Pos_all_nosliceidOne_mideg","h_cosTheta_Pos_all_nosliceidOne_mideg", n_cosine, cosine_min, cosine_max);
	h_cosTheta_Pos_all_nosliceidOne_midother=new TH1D("h_cosTheta_Pos_all_nosliceidOne_midother","h_cosTheta_Pos_all_nosliceidOne_midother", n_cosine, cosine_min, cosine_max);


	h_cosTheta_Pos_all_nosliceidOne->Sumw2();
	h_cosTheta_Pos_all_nosliceidOne_inel->Sumw2();
	h_cosTheta_Pos_all_nosliceidOne_el->Sumw2();
	h_cosTheta_Pos_all_nosliceidOne_midcosmic->Sumw2();
	h_cosTheta_Pos_all_nosliceidOne_midpi->Sumw2();
	h_cosTheta_Pos_all_nosliceidOne_midp->Sumw2();
	h_cosTheta_Pos_all_nosliceidOne_midmu->Sumw2();
	h_cosTheta_Pos_all_nosliceidOne_mideg->Sumw2();
	h_cosTheta_Pos_all_nosliceidOne_midother->Sumw2();


} //BookHistograms


void SaveHistograms() { //SaveHistograms
	outputFile->cd();
	outputFile->Write();

} //SaveHistograms


