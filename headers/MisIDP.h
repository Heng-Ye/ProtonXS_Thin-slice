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

TProfile2D *h2d_recotrklen_truetrklen_inel;
TProfile2D *h2d_recotrklen_truetrklen_el;
TProfile2D *h2d_recotrklen_truetrklen_misidp;

TH3D *h3d_recotrklen_truetrklen_cosTheta_inel;
TH3D *h3d_recotrklen_truetrklen_cosTheta_el;
TH3D *h3d_recotrklen_truetrklen_cosTheta_misidp;


//cosTheta
TProfile2D *h2d_ntrklen_chi2_inel;
TProfile2D *h2d_ntrklen_chi2_el;
TProfile2D *h2d_ntrklen_chi2_misidp;

TH3D *h3d_ntrklen_chi2_cosTheta_inel;
TH3D *h3d_ntrklen_chi2_cosTheta_el;
TH3D *h3d_ntrklen_chi2_cosTheta_misidp;


//cosTheta_slicing
const int nn_cos=10;
//track length
TProfile2D *tp2d_range_reco_true_inel[nn_cos];
TProfile2D *tp2d_range_reco_true_el[nn_cos];
TProfile2D *tp2d_range_reco_true_misidp[nn_cos];

//ntrklen vs chi2pid
TProfile2D *tp2d_ntrklen_chi2_inel[nn_cos];
TProfile2D *tp2d_ntrklen_chi2_el[nn_cos];
TProfile2D *tp2d_ntrklen_chi2_misidp[nn_cos];


void BookHistograms() { //BookHistograms
	outputFile = TFile::Open(fOutputFileName.c_str(), "recreate"); //open output file

	int n_2d=120;
	double trklen_min=0;
	double trklen_max=120;

	int n_cos=100;
	double cos_min=0;
	double cos_max=1;

	h2d_recotrklen_truetrklen_inel=new TProfile2D("h2d_recotrklen_truetrklen_inel","",n_2d, trklen_min, trklen_max, n_2d, trklen_min, trklen_max);
	h2d_recotrklen_truetrklen_el=new TProfile2D("h2d_recotrklen_truetrklen_el","",n_2d, trklen_min, trklen_max, n_2d, trklen_min, trklen_max);
	h2d_recotrklen_truetrklen_misidp=new TProfile2D("h2d_recotrklen_truetrklen_misidp","",n_2d, trklen_min, trklen_max, n_2d, trklen_min, trklen_max);

	
	h2d_recotrklen_truetrklen_inel->SetTitle("inel; Reco Track Length[cm]; True Track Length[cm]; cos#Theta");
	h2d_recotrklen_truetrklen_el->SetTitle("El; Reco Track Length[cm]; True Track Length[cm]; cos#Theta");
	h2d_recotrklen_truetrklen_misidp->SetTitle("MisID:p; Reco Track Length[cm]; True Track Length[cm]; cos#Theta");

	h2d_recotrklen_truetrklen_inel->Sumw2();
	h2d_recotrklen_truetrklen_el->Sumw2();
	h2d_recotrklen_truetrklen_misidp->Sumw2();

	h3d_recotrklen_truetrklen_cosTheta_inel=new TH3D("h3d_recotrklen_truetrklen_cosTheta_inel","",n_2d, trklen_min, trklen_max, n_2d, trklen_min, trklen_max, n_cos, cos_min, cos_max);
	h3d_recotrklen_truetrklen_cosTheta_el=new TH3D("h3d_recotrklen_truetrklen_cosTheta_el","",n_2d, trklen_min, trklen_max, n_2d, trklen_min, trklen_max, n_cos, cos_min,cos_max);
	h3d_recotrklen_truetrklen_cosTheta_misidp=new TH3D("h3d_recotrklen_truetrklen_cosTheta_misidp","",n_2d, trklen_min, trklen_max, n_2d, trklen_min, trklen_max, n_cos, cos_min, cos_max);

	h3d_recotrklen_truetrklen_cosTheta_inel->SetTitle("inel; Reco Track Length[cm]; True Track Length[cm]; cos#Theta");
	h3d_recotrklen_truetrklen_cosTheta_el->SetTitle("inel; Reco Track Length[cm]; True Track Length[cm]; cos#Theta");
	h3d_recotrklen_truetrklen_cosTheta_misidp->SetTitle("inel; Reco Track Length[cm]; True Track Length[cm]; cos#Theta");

	h3d_recotrklen_truetrklen_cosTheta_inel->Sumw2();
	h3d_recotrklen_truetrklen_cosTheta_el->Sumw2();
	h3d_recotrklen_truetrklen_cosTheta_misidp->Sumw2();

	float dcos=0.1;
	for (int j=0; j<nn_cos; ++j) {
		float tmp_min=(float)j*dcos;
		float tmp_max=tmp_min+dcos;
		tp2d_range_reco_true_inel[j]=new TProfile2D(Form("tp2d_range_reco_true_inel_%d",j), Form("Cos#Theta:%.1f-%.1f",tmp_min,tmp_max), n_2d, trklen_min, trklen_max, n_2d, trklen_min, trklen_max);
		tp2d_range_reco_true_el[j]=new TProfile2D(Form("tp2d_range_reco_true_el_%d",j), Form("Cos#Theta:%.1f-%.1f",tmp_min,tmp_max), n_2d, trklen_min, trklen_max, n_2d, trklen_min, trklen_max);
		tp2d_range_reco_true_misidp[j]=new TProfile2D(Form("tp2d_range_reco_true_misidp_%d",j), Form("Cos#Theta:%.1f-%.1f",tmp_min,tmp_max), n_2d, trklen_min, trklen_max, n_2d, trklen_min, trklen_max);

		tp2d_range_reco_true_inel[j]->Sumw2();
		tp2d_range_reco_true_el[j]->Sumw2();
		tp2d_range_reco_true_misidp[j]->Sumw2();
	}


        int n_chi2=300;
	double chi2_min=0;
	double chi2_max=150;

	int n_ntrklen=120;
	double ntrllen_min=0;
	double ntrklen_max=1.2;
	
	h2d_ntrklen_chi2_inel=new TProfile2D("h2d_ntrklen_chi2_inel","", n_ntrklen, ntrllen_min, ntrklen_max, n_chi2, chi2_min, chi2_max);
	h2d_ntrklen_chi2_el=new TProfile2D("h2d_ntrklen_chi2_el","", n_ntrklen, ntrllen_min, ntrklen_max, n_chi2, chi2_min, chi2_max);
	h2d_ntrklen_chi2_misidp=new TProfile2D("h2d_ntrklen_chi2_misidp","",n_ntrklen, ntrllen_min, ntrklen_max, n_chi2, chi2_min, chi2_max);

	h2d_ntrklen_chi2_inel->SetTitle("inel; Proton Track Length/CSDA; #chi^{2} PID; cos#Theta");
	h2d_ntrklen_chi2_el->SetTitle("inel; Proton Track Length/CSDA; #chi^{2} PID; cos#Theta");
	h2d_ntrklen_chi2_misidp->SetTitle("inel; Proton Track Length/CSDA; #chi^{2} PID; cos#Theta");

	h2d_ntrklen_chi2_inel->Sumw2();
	h2d_ntrklen_chi2_el->Sumw2();
	h2d_ntrklen_chi2_misidp->Sumw2();

	h3d_ntrklen_chi2_cosTheta_inel=new TH3D("h2d_ntrklen_chi2_cosTheta_inel","", n_ntrklen, ntrllen_min, ntrklen_max, n_chi2, chi2_min, chi2_max, n_cos, cos_min, cos_max);
	h3d_ntrklen_chi2_cosTheta_el=new TH3D("h2d_ntrklen_chi2_cosTheta_el","", n_ntrklen, ntrllen_min, ntrklen_max, n_chi2, chi2_min, chi2_max, n_cos, cos_min, cos_max);
	h3d_ntrklen_chi2_cosTheta_misidp=new TH3D("h2d_ntrklen_chi2_cosTheta_misidp","", n_ntrklen, ntrllen_min, ntrklen_max, n_chi2, chi2_min, chi2_max, n_cos, cos_min, cos_max);

	h3d_ntrklen_chi2_cosTheta_inel->SetTitle("inel; Proton Track Length/CSDA; #chi^{2} PID; cos#Theta");
	h3d_ntrklen_chi2_cosTheta_el->SetTitle("inel; Proton Track Length/CSDA; #chi^{2} PID; cos#Theta");
	h3d_ntrklen_chi2_cosTheta_misidp->SetTitle("inel; Proton Track Length/CSDA; #chi^{2} PID; cos#Theta");

	h3d_ntrklen_chi2_cosTheta_inel->Sumw2();
	h3d_ntrklen_chi2_cosTheta_el->Sumw2();
	h3d_ntrklen_chi2_cosTheta_misidp->Sumw2();

	
	for (int j=0; j<nn_cos; ++j) {
		float tmp_min=(float)j*dcos;
		float tmp_max=tmp_min+dcos;

		tp2d_ntrklen_chi2_inel[j]=new TProfile2D(Form("tp2d_ntrklen_chi2_inel_%d", j), Form("Cos#Theta:%.1f-%.1f",tmp_min,tmp_max), n_ntrklen, ntrllen_min, ntrklen_max, n_chi2, chi2_min, chi2_max);
		tp2d_ntrklen_chi2_el[j]=new TProfile2D(Form("tp2d_ntrklen_chi2_el_%d", j), Form("Cos#Theta:%.1f-%.1f",tmp_min,tmp_max), n_ntrklen, ntrllen_min, ntrklen_max, n_chi2, chi2_min, chi2_max);
		tp2d_ntrklen_chi2_misidp[j]=new TProfile2D(Form("tp2d_ntrklen_chi2_misidp_%d", j), Form("Cos#Theta:%.1f-%.1f",tmp_min,tmp_max), n_ntrklen, ntrllen_min, ntrklen_max, n_chi2, chi2_min, chi2_max);


		tp2d_ntrklen_chi2_inel[j]->Sumw2();
		tp2d_ntrklen_chi2_el[j]->Sumw2();
		tp2d_ntrklen_chi2_misidp[j]->Sumw2();
	}


} //BookHistograms


void SaveHistograms() { //SaveHistograms
	outputFile->cd();
	outputFile->Write();

} //SaveHistograms


