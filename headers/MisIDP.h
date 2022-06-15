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

TH2D *h2d_ntrklen_chi2_misidp_truerangeeq0;
TH2D *h2d_ntrklen_chi2_misidp_truerangegt0;


TH3D *h3d_ntrklen_chi2_cosTheta_inel;
TH3D *h3d_ntrklen_chi2_cosTheta_el;
TH3D *h3d_ntrklen_chi2_cosTheta_misidp;

TH1D *h1d_cosTheta_misidp_truerangeeq0;
TH1D *h1d_cosTheta_misidp_truerangegt0;

//cosTheta_slicing
const int nn_cos=25;
//track length
TProfile2D *tp2d_range_reco_true_inel[nn_cos];
TProfile2D *tp2d_range_reco_true_el[nn_cos];
TProfile2D *tp2d_range_reco_true_misidp[nn_cos];

TProfile2D *tp2d_range_reco_true_chi2pid_inel;
TProfile2D *tp2d_range_reco_true_chi2pid_el;
TProfile2D *tp2d_range_reco_true_chi2pid_misidp;


//ntrklen vs chi2pid
TProfile2D *tp2d_ntrklen_chi2_inel[nn_cos];
TProfile2D *tp2d_ntrklen_chi2_el[nn_cos];
TProfile2D *tp2d_ntrklen_chi2_misidp[nn_cos];


//eff vs reco_range
TProfile2D *tp2d_recorange_eff_inel[nn_cos];
TProfile2D *tp2d_recorange_eff_el[nn_cos];
TProfile2D *tp2d_recorange_eff_misidp[nn_cos];


//dedx vs rr
TH2D *h2d_rr_dedx_el;
TH2D *h2d_rr_dedx_inel;
TH2D *h2d_rr_dedx_misidp;

//eff_vs_range
TH2D *h2d_recolen_eff_inel;
TH2D *h2d_recolen_eff_el;
TH2D *h2d_recolen_eff_misidp;
TH2D *h2d_recolen_eff_all;

//x-y dis
TH2D *h2d_true_xy_upstream_misidp;
TH2D *h2d_reco_xy_upstream_misidp;
TH2D *h2d_true_xy_upstream_inel;
TH2D *h2d_reco_xy_upstream_inel;
TH2D *h2d_true_xy_upstream_el;
TH2D *h2d_reco_xy_upstream_el;

TH2D *h2d_true_xy_misidp;
TH2D *h2d_true_xy_inel;
TH2D *h2d_true_xy_el;

//dxy vs cosine
TH2D *h2d_dxy_cosine_BQ_misidp;
TH2D *h2d_dxy_cosine_BQ_misidp_lenle0;
TH2D *h2d_dxy_cosine_BQ_misidp_lengt0;
TH2D *h2d_dxy_cosine_BQ_inel;
TH2D *h2d_dxy_cosine_BQ_el;

TH1D *h1d_dxy_BQ_misidp;
TH1D *h1d_dxy_BQ_misidp_lenle0;
TH1D *h1d_dxy_BQ_misidp_lengt0;
TH1D *h1d_dxy_BQ_inel;
TH1D *h1d_dxy_BQ_el;


TH2D *h2d_dxy_cosine_Pos_misidp;
TH2D *h2d_dxy_cosine_Pos_misidp_lenle0;
TH2D *h2d_dxy_cosine_Pos_misidp_lengt0;
TH2D *h2d_dxy_cosine_Pos_inel;
TH2D *h2d_dxy_cosine_Pos_el;

TH1D *h1d_dxy_Pos_misidp;
TH1D *h1d_dxy_Pos_misidp_lenle0;
TH1D *h1d_dxy_Pos_misidp_lengt0;
TH1D *h1d_dxy_Pos_inel;
TH1D *h1d_dxy_Pos_el;

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

	float dcos=0.04;
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
	tp2d_range_reco_true_chi2pid_inel=new TProfile2D(Form("tp2d_range_reco_true_chi2pid_inel"), Form(""), n_2d, trklen_min, trklen_max, n_2d, trklen_min, trklen_max);
	tp2d_range_reco_true_chi2pid_el=new TProfile2D(Form("tp2d_range_reco_true_chi2pid_el"), Form(""), n_2d, trklen_min, trklen_max, n_2d, trklen_min, trklen_max);
	tp2d_range_reco_true_chi2pid_misidp=new TProfile2D(Form("tp2d_range_reco_true_chi2pid_misidp"), Form(""), n_2d, trklen_min, trklen_max, n_2d, trklen_min, trklen_max);

	tp2d_range_reco_true_chi2pid_inel->Sumw2();
	tp2d_range_reco_true_chi2pid_el->Sumw2();
	tp2d_range_reco_true_chi2pid_misidp->Sumw2();


        int n_chi2=500;
	double chi2_min=0;
	double chi2_max=250;

	int n_ntrklen=120;
	double ntrllen_min=0;
	double ntrklen_max=1.2;
	
	h2d_ntrklen_chi2_inel=new TProfile2D("h2d_ntrklen_chi2_inel","", n_ntrklen, ntrllen_min, ntrklen_max, n_chi2, chi2_min, chi2_max);
	h2d_ntrklen_chi2_el=new TProfile2D("h2d_ntrklen_chi2_el","", n_ntrklen, ntrllen_min, ntrklen_max, n_chi2, chi2_min, chi2_max);
	h2d_ntrklen_chi2_misidp=new TProfile2D("h2d_ntrklen_chi2_misidp","",n_ntrklen, ntrllen_min, ntrklen_max, n_chi2, chi2_min, chi2_max);

	h2d_ntrklen_chi2_misidp_truerangeeq0=new TH2D("h2d_ntrklen_chi2_misidp_truerangeeq0","",n_ntrklen, ntrllen_min, ntrklen_max, n_chi2, chi2_min, chi2_max);
	h2d_ntrklen_chi2_misidp_truerangegt0=new TH2D("h2d_ntrklen_chi2_misidp_truerangegt0","",n_ntrklen, ntrllen_min, ntrklen_max, n_chi2, chi2_min, chi2_max);

	h2d_ntrklen_chi2_inel->SetTitle("Inel; Proton Track Length/CSDA; #chi^{2} PID; cos#Theta");
	h2d_ntrklen_chi2_el->SetTitle("El.; Proton Track Length/CSDA; #chi^{2} PID; cos#Theta");
	h2d_ntrklen_chi2_misidp->SetTitle("MisID:p; Proton Track Length/CSDA; #chi^{2} PID; cos#Theta");

	h2d_ntrklen_chi2_misidp_truerangeeq0->SetTitle("MisID:p; Proton Track Length/CSDA; #chi^{2} PID");
	h2d_ntrklen_chi2_misidp_truerangegt0->SetTitle("MisID:p; Proton Track Length/CSDA; #chi^{2} PID");

	h2d_ntrklen_chi2_inel->Sumw2();
	h2d_ntrklen_chi2_el->Sumw2();
	h2d_ntrklen_chi2_misidp->Sumw2();

	h2d_ntrklen_chi2_misidp_truerangeeq0->Sumw2();
	h2d_ntrklen_chi2_misidp_truerangegt0->Sumw2();

	h3d_ntrklen_chi2_cosTheta_inel=new TH3D("h2d_ntrklen_chi2_cosTheta_inel","", n_ntrklen, ntrllen_min, ntrklen_max, n_chi2, chi2_min, chi2_max, n_cos, cos_min, cos_max);
	h3d_ntrklen_chi2_cosTheta_el=new TH3D("h2d_ntrklen_chi2_cosTheta_el","", n_ntrklen, ntrllen_min, ntrklen_max, n_chi2, chi2_min, chi2_max, n_cos, cos_min, cos_max);
	h3d_ntrklen_chi2_cosTheta_misidp=new TH3D("h2d_ntrklen_chi2_cosTheta_misidp","", n_ntrklen, ntrllen_min, ntrklen_max, n_chi2, chi2_min, chi2_max, n_cos, cos_min, cos_max);

	h3d_ntrklen_chi2_cosTheta_inel->SetTitle("inel; Proton Track Length/CSDA; #chi^{2} PID; cos#Theta");
	h3d_ntrklen_chi2_cosTheta_el->SetTitle("inel; Proton Track Length/CSDA; #chi^{2} PID; cos#Theta");
	h3d_ntrklen_chi2_cosTheta_misidp->SetTitle("inel; Proton Track Length/CSDA; #chi^{2} PID; cos#Theta");

	h3d_ntrklen_chi2_cosTheta_inel->Sumw2();
	h3d_ntrklen_chi2_cosTheta_el->Sumw2();
	h3d_ntrklen_chi2_cosTheta_misidp->Sumw2();

	h1d_cosTheta_misidp_truerangeeq0=new TH1D("h1d_cosTheta_misidp_truerangeeq0","",n_cos, cos_min, cos_max);
	h1d_cosTheta_misidp_truerangegt0=new TH1D("h1d_cosTheta_misidp_truerangegt0","",n_cos, cos_min, cos_max);
	h1d_cosTheta_misidp_truerangeeq0->Sumw2();
	h1d_cosTheta_misidp_truerangegt0->Sumw2();
	
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

	//dedx vs rr
	h2d_rr_dedx_el=new TH2D("h2d_rr_dedx_el","",240,0,120,90,0,30);
	h2d_rr_dedx_inel=new TH2D("h2d_rr_dedx_inel","",240,0,120,90,0,30);
	h2d_rr_dedx_misidp=new TH2D("h2d_rr_dedx_misidp","",240,0,120,90,0,30);

	h2d_rr_dedx_el->Sumw2();
	h2d_rr_dedx_inel->Sumw2();
	h2d_rr_dedx_misidp->Sumw2();

	//reco_range vs eff
	int n_eff=350;
	double eff_min=-0.5;
	double eff_max=3; 

	h2d_recolen_eff_inel=new TH2D("h2d_recolen_eff_inel","Inel.", n_2d, trklen_min, trklen_max, n_eff, eff_min, eff_max);
	h2d_recolen_eff_el=new TH2D("h2d_recolen_eff_el","El.", n_2d, trklen_min, trklen_max, n_eff, eff_min, eff_max);
	h2d_recolen_eff_misidp=new TH2D("h2d_recolen_eff_misidp","El.", n_2d, trklen_min, trklen_max, n_eff, eff_min, eff_max);
	h2d_recolen_eff_all=new TH2D("h2d_recolen_eff_all","all.", n_2d, trklen_min, trklen_max, n_eff, eff_min, eff_max);
	h2d_recolen_eff_inel->Sumw2();
	h2d_recolen_eff_el->Sumw2();
	h2d_recolen_eff_misidp->Sumw2();
	h2d_recolen_eff_all->Sumw2();

	for (int j=0; j<nn_cos; ++j) {
		float tmp_min=(float)j*dcos;
		float tmp_max=tmp_min+dcos;

		tp2d_recorange_eff_inel[j]=new TProfile2D(Form("tp2d_recorange_eff_inel_%d", j), Form("Cos#Theta:%.1f-%.1f",tmp_min,tmp_max), n_2d, trklen_min, trklen_max, n_eff, eff_min, eff_max);
		tp2d_recorange_eff_el[j]=new TProfile2D(Form("tp2d_recorange_eff_el_%d", j), Form("Cos#Theta:%.1f-%.1f",tmp_min,tmp_max), n_2d, trklen_min, trklen_max, n_eff, eff_min, eff_max);
		tp2d_recorange_eff_misidp[j]=new TProfile2D(Form("tp2d_recorange_eff_misidp_%d", j), Form("Cos#Theta:%.1f-%.1f",tmp_min,tmp_max), n_2d, trklen_min, trklen_max, n_eff, eff_min, eff_max);

		tp2d_recorange_eff_inel[j]->Sumw2();
		tp2d_recorange_eff_el[j]->Sumw2();
		tp2d_recorange_eff_misidp[j]->Sumw2();
	}


	h2d_true_xy_upstream_misidp=new TH2D("h2d_true_xy_upstream_misidp","",70,-60,-10,60,390,450);
	h2d_true_xy_upstream_inel=new TH2D("h2d_true_xy_upstream_inel","",70,-60,-10,60,390,450);
	h2d_true_xy_upstream_el=new TH2D("h2d_true_xy_upstream_el","",70,-60,-10,60,390,450);
	h2d_reco_xy_upstream_misidp=new TH2D("h2d_reco_xy_upstream_misidp","",70,-60,-10,60,390,450);
	h2d_reco_xy_upstream_inel=new TH2D("h2d_reco_xy_upstream_inel","",70,-60,-10,60,390,450);
	h2d_reco_xy_upstream_el=new TH2D("h2d_reco_xy_upstream_el","",70,-60,-10,60,390,450);
	h2d_true_xy_misidp=new TH2D("h2d_true_xy_misidp","",70,-60,-10,60,390,450);
	h2d_true_xy_inel=new TH2D("h2d_true_xy_inel","",70,-60,-10,60,390,450);
	h2d_true_xy_el=new TH2D("h2d_true_xy_el","",70,-60,-10,60,390,450);
	h2d_true_xy_upstream_misidp->Sumw2();
	h2d_true_xy_upstream_inel->Sumw2();
	h2d_true_xy_upstream_el->Sumw2();
	h2d_reco_xy_upstream_misidp->Sumw2();
	h2d_reco_xy_upstream_inel->Sumw2();
	h2d_reco_xy_upstream_el->Sumw2();

	h2d_true_xy_misidp->Sumw2();
	h2d_true_xy_inel->Sumw2();
	h2d_true_xy_el->Sumw2();

	int n_dxy=500;
	double dxy_min=0;
	double dxy_max=50;
	int n_cos1=100;
	double cos1_min=0.;
	double cos1_max=1;
	h2d_dxy_cosine_BQ_misidp=new TH2D("h2d_dxy_cosine_BQ_misidp","",n_dxy,dxy_min,dxy_max,n_cos1,cos1_min,cos1_max);
	h2d_dxy_cosine_BQ_misidp_lenle0=new TH2D("h2d_dxy_cosine_BQ_misidp_lenle0","",n_dxy,dxy_min,dxy_max,n_cos1,cos1_min,cos1_max);
	h2d_dxy_cosine_BQ_misidp_lengt0=new TH2D("h2d_dxy_cosine_BQ_misidp_lengt0","",n_dxy,dxy_min,dxy_max,n_cos1,cos1_min,cos1_max);
	h2d_dxy_cosine_BQ_inel=new TH2D("h2d_dxy_cosine_BQ_inel","",n_dxy,dxy_min,dxy_max,n_cos1,cos1_min,cos1_max);
	h2d_dxy_cosine_BQ_el=new TH2D("h2d_dxy_cosine_BQ_el","",n_dxy,dxy_min,dxy_max,n_cos1,cos1_min,cos1_max);

	h2d_dxy_cosine_BQ_misidp->Sumw2();
	h2d_dxy_cosine_BQ_misidp_lenle0->Sumw2();
	h2d_dxy_cosine_BQ_misidp_lengt0->Sumw2();
	h2d_dxy_cosine_BQ_inel->Sumw2();
	h2d_dxy_cosine_BQ_el->Sumw2();	


	h1d_dxy_BQ_misidp=new TH1D("h1d_dxy_BQ_misidp","",n_dxy,dxy_min,dxy_max);
	h1d_dxy_BQ_misidp_lenle0=new TH1D("h1d_dxy_BQ_misidp_lenle0","",n_dxy,dxy_min,dxy_max);
	h1d_dxy_BQ_misidp_lengt0=new TH1D("h1d_dxy_BQ_misidp_lengt0","",n_dxy,dxy_min,dxy_max);
	h1d_dxy_BQ_inel=new TH1D("h1d_dxy_BQ_inel","",n_dxy,dxy_min,dxy_max);
	h1d_dxy_BQ_el=new TH1D("h1d_dxy_BQ_el","",n_dxy,dxy_min,dxy_max);

	h1d_dxy_BQ_misidp->Sumw2();
	h1d_dxy_BQ_misidp_lenle0->Sumw2();
	h1d_dxy_BQ_misidp_lengt0->Sumw2();
	h1d_dxy_BQ_inel->Sumw2();
	h1d_dxy_BQ_el->Sumw2();	

	//

	h2d_dxy_cosine_Pos_misidp=new TH2D("h2d_dxy_cosine_Pos_misidp","",n_dxy,dxy_min,dxy_max,n_cos1,cos1_min,cos1_max);
	h2d_dxy_cosine_Pos_misidp_lenle0=new TH2D("h2d_dxy_cosine_Pos_misidp_lenle0","",n_dxy,dxy_min,dxy_max,n_cos1,cos1_min,cos1_max);
	h2d_dxy_cosine_Pos_misidp_lengt0=new TH2D("h2d_dxy_cosine_Pos_misidp_lengt0","",n_dxy,dxy_min,dxy_max,n_cos1,cos1_min,cos1_max);
	h2d_dxy_cosine_Pos_inel=new TH2D("h2d_dxy_cosine_Pos_inel","",n_dxy,dxy_min,dxy_max,n_cos1,cos1_min,cos1_max);
	h2d_dxy_cosine_Pos_el=new TH2D("h2d_dxy_cosine_Pos_el","",n_dxy,dxy_min,dxy_max,n_cos1,cos1_min,cos1_max);

	h2d_dxy_cosine_Pos_misidp->Sumw2();
	h2d_dxy_cosine_Pos_misidp_lenle0->Sumw2();
	h2d_dxy_cosine_Pos_misidp_lengt0->Sumw2();
	h2d_dxy_cosine_Pos_inel->Sumw2();
	h2d_dxy_cosine_Pos_el->Sumw2();	


	h1d_dxy_Pos_misidp=new TH1D("h1d_dxy_Pos_misidp","",n_dxy,dxy_min,dxy_max);
	h1d_dxy_Pos_misidp_lenle0=new TH1D("h1d_dxy_Pos_misidp_lenle0","",n_dxy,dxy_min,dxy_max);
	h1d_dxy_Pos_misidp_lengt0=new TH1D("h1d_dxy_Pos_misidp_lengt0","",n_dxy,dxy_min,dxy_max);
	h1d_dxy_Pos_inel=new TH1D("h1d_dxy_Pos_inel","",n_dxy,dxy_min,dxy_max);
	h1d_dxy_Pos_el=new TH1D("h1d_dxy_Pos_el","",n_dxy,dxy_min,dxy_max);

	h1d_dxy_Pos_misidp->Sumw2();
	h1d_dxy_Pos_misidp_lenle0->Sumw2();
	h1d_dxy_Pos_misidp_lengt0->Sumw2();
	h1d_dxy_Pos_inel->Sumw2();
	h1d_dxy_Pos_el->Sumw2();	


} //BookHistograms


void SaveHistograms() { //SaveHistograms
	outputFile->cd();
	outputFile->Write();

} //SaveHistograms


