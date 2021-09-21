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

TH2D *h2d_recotrklen_truetrklen_inel;
TH2D *h2d_recotrklen_truetrklen_el;
TH2D *h2d_recotrklen_truetrklen_misidp;

TH3D *h3d_recotrklen_truetrklen_cosTheta_inel;
TH3D *h3d_recotrklen_truetrklen_cosTheta_el;
TH3D *h3d_recotrklen_truetrklen_cosTheta_misidp;


//cosTheta
TH2D *h2d_ntrklen_chi2_inel;
TH2D *h2d_ntrklen_chi2_el;
TH2D *h2d_ntrklen_chi2_misidp;

TH3D *h3d_ntrklen_chi2_cosTheta_inel;
TH3D *h3d_ntrklen_chi2_cosTheta_el;
TH3D *h3d_ntrklen_chi2_cosTheta_misidp;




void BookHistograms() { //BookHistograms
	outputFile = TFile::Open(fOutputFileName.c_str(), "recreate"); //open output file

	int n_2d=120;
	double trklen_min=0;
	double trklen_max=120;

	int n_cos=100;
	double cos_min=0;
	double cos_max=1;

	h2d_recotrklen_truetrklen_inel=new TH2D("h2d_recotrklen_truetrklen_inel","",n_2d,trklen_min,trklen_max,n_2d,trklen_min,trklen_max);
	h2d_recotrklen_truetrklen_el=new TH2D("h2d_recotrklen_truetrklen_el","",n_2d,trklen_min,trklen_max,n_2d,trklen_min,trklen_max);
	h2d_recotrklen_truetrklen_misidp=new TH2D("h2d_recotrklen_truetrklen_misidp","",n_2d,trklen_min,trklen_max,n_2d,trklen_min,trklen_max);

	h2d_recotrklen_truetrklen_inel->Sumw2();
	h2d_recotrklen_truetrklen_el->Sumw2();
	h2d_recotrklen_truetrklen_misidp->Sumw2();

	h3d_recotrklen_truetrklen_cosTheta_inel=new TH3D("h3d_recotrklen_truetrklen_cosTheta_inel","",n_2d,trklen_min,trklen_max,n_2d,trklen_min,trklen_max,n_cos,cos_min,cos_max);
	h3d_recotrklen_truetrklen_cosTheta_el=new TH3D("h3d_recotrklen_truetrklen_cosTheta_el","",n_2d,trklen_min,trklen_max,n_2d,trklen_min,trklen_max,n_cos,cos_min,cos_max);
	h3d_recotrklen_truetrklen_cosTheta_misidp=new TH3D("h3d_recotrklen_truetrklen_cosTheta_misidp","",n_2d,trklen_min,trklen_max,n_2d,trklen_min,trklen_max,n_cos,cos_min,cos_max);

	h3d_recotrklen_truetrklen_cosTheta_inel->Sumw2();
	h3d_recotrklen_truetrklen_cosTheta_el->Sumw2();
	h3d_recotrklen_truetrklen_cosTheta_misidp->Sumw2();

        int n_chi2=300;
	double chi2_min=0;
	double chi2_max=150;

	int n_ntrklen=120;
	double ntrllen_min=0;
	double ntrklen_max=1.2;
	
	h2d_ntrklen_chi2_inel=new TH2D("h2d_ntrklen_chi2_inel","",n_ntrklen,ntrllen_min,ntrklen_max,n_chi2,chi2_min,ntrklen_max);
	h2d_ntrklen_chi2_el=new TH2D("h2d_ntrklen_chi2_el","",n_ntrklen,ntrllen_min,ntrklen_max,n_chi2,chi2_min,ntrklen_max);
	h2d_ntrklen_chi2_misidp=new TH2D("h2d_ntrklen_chi2_misidp","",n_ntrklen,ntrllen_min,ntrklen_max,n_chi2,chi2_min,ntrklen_max);

	h2d_ntrklen_chi2_inel->Sumw2();
	h2d_ntrklen_chi2_el->Sumw2();
	h2d_ntrklen_chi2_misidp->Sumw2();

	h3d_ntrklen_chi2_cosTheta_inel=new TH2D("h2d_ntrklen_chi2_cosTheta_inel","",n_ntrklen,ntrllen_min,ntrklen_max,n_chi2,chi2_min,ntrklen_max,n_chi2,chi2_min,chi2_max);
	h3d_ntrklen_chi2_cosTheta_el=new TH2D("h2d_ntrklen_chi2_cosTheta_el","",n_ntrklen,ntrllen_min,ntrklen_max,n_chi2,chi2_min,ntrklen_max,n_chi2,chi2_min,chi2_max);
	h3d_ntrklen_chi2_cosTheta_misidp=new TH2D("h2d_ntrklen_chi2_cosTheta_misidp","",n_ntrklen,ntrllen_min,ntrklen_max,n_chi2,chi2_min,ntrklen_max,n_chi2,chi2_min,chi2_max);

	h3d_ntrklen_chi2_cosTheta_inel->Sumw2();
	h3d_ntrklen_chi2_cosTheta_el->Sumw2();
	h3d_ntrklen_chi2_cosTheta_misidp->Sumw2();

} //BookHistograms


void SaveHistograms() { //SaveHistograms
	outputFile->cd();
	outputFile->Write();

} //SaveHistograms


