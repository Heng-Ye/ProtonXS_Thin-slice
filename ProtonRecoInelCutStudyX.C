#define ProtonRecoInelCutStudy_cxx
#include "ProtonRecoInelCutStudy.h"

#include <TH2.h>
#include <TH1.h>
#include "TH2D.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>

#include <TMath.h>
#include <TLine.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH3F.h>
#include <TString.h>
#include <TProfile2D.h>
#include <THStack.h>
#include "TGraph.h" 
#include "TGraphSmooth.h" 
#include "TParameter.h"
#include "TGraphErrors.h"
#include "string"
#include "vector"
#include "TSpline.h"
#include "TH3F.h"
#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>
#include "TVector3.h"

#include <stdio.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "./cali/dedx_function_35ms.h"
#include "./headers/BasicParameters.h"
#include "./headers/BasicFunctions.h"
#include "./headers/util.h"
//#include "./headers/SliceParams.h"
//#include "./headers/ThinSlice.h"

using namespace std;
using namespace ROOT::Math;

void ProtonRecoInelCutStudy::Loop() {
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	Long64_t nbytes = 0, nb = 0;
	bool isTestSample=true;
	//int true_sliceID = -1, reco_sliceID = -1;

	//counters
	int n_tot=0, n_sg=0, n_bkg=0; //no cut
	int n_el=0, n_inel=0, n_midcosmic=0, n_midpi=0, n_midp=0, n_midmu=0, n_mideg=0, n_midother=0;

	int n_pan_tot=0; //pandora cut
	int n_el_pan=0, n_inel_pan=0, n_midcosmic_pan=0, n_midpi_pan=0, n_midp_pan=0, n_midmu_pan=0, n_mideg_pan=0, n_midother_pan=0;

	int n_calsz_tot=0; //calosz cut
	int n_el_calsz=0, n_inel_calsz=0, n_midcosmic_calsz=0, n_midpi_calsz=0, n_midp_calsz=0, n_midmu_calsz=0, n_mideg_calsz=0, n_midother_calsz=0;

	int n_bq_tot=0; //bq cut
	int n_el_bq=0, n_inel_bq=0, n_midcosmic_bq=0, n_midpi_bq=0, n_midp_bq=0, n_midmu_bq=0, n_mideg_bq=0, n_midother_bq=0;

	int n_recoinel_tot=0; //recoinel cut
	int n_el_recoinel=0, n_inel_recoinel=0, n_midcosmic_recoinel=0, n_midpi_recoinel=0, n_midp_recoinel=0, n_midmu_recoinel=0, n_mideg_recoinel=0, n_midother_recoinel=0;

	int n_recoinel2_tot=0; //recoinel2 cut
	int n_el_recoinel2=0, n_inel_recoinel2=0, n_midcosmic_recoinel2=0, n_midpi_recoinel2=0, n_midp_recoinel2=0, n_midmu_recoinel2=0, n_mideg_recoinel2=0, n_midother_recoinel2=0;

	int n_recoinel3_tot=0; //recoinel3 cut
	int n_el_recoinel3=0, n_inel_recoinel3=0, n_midcosmic_recoinel3=0, n_midpi_recoinel3=0, n_midp_recoinel3=0, n_midmu_recoinel3=0, n_mideg_recoinel3=0, n_midother_recoinel3=0;



	//book histograms ------------------------------------------------------------------------------------------------------------------------//
	int n_ntrklen=61;
	double st_ntrklen=-0.02;
	double ed_ntrklen=1.2;

	int n_chi2=750;
	double st_chi2=0;
	double ed_chi2=150;

	//2d ntrklen vs chi2pid
	TH2D *ntrklen_chi2pid_BQ=new TH2D("ntrklen_chi2pid_BQ","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntrklen_chi2pid_BQ_inel=new TH2D("ntrklen_chi2pid_BQ_inel","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntrklen_chi2pid_BQ_el=new TH2D("ntrklen_chi2pid_BQ_el","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntrklen_chi2pid_BQ_midcosmic=new TH2D("ntrklen_chi2pid_BQ_midcosmic","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntrklen_chi2pid_BQ_midpi=new TH2D("ntrklen_chi2pid_BQ_midpi","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntrklen_chi2pid_BQ_midp=new TH2D("ntrklen_chi2pid_BQ_midp","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntrklen_chi2pid_BQ_midmu=new TH2D("ntrklen_chi2pid_BQ_midmu","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntrklen_chi2pid_BQ_mideg=new TH2D("ntrklen_chi2pid_BQ_mideg","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntrklen_chi2pid_BQ_midother=new TH2D("ntrklen_chi2pid_BQ_midother","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntrklen_chi2pid_RecoInel=new TH2D("ntrklen_chi2pid_RecoInel","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntrklen_chi2pid_RecoInel2=new TH2D("ntrklen_chi2pid_RecoInel2","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntrklen_chi2pid_RecoInel3=new TH2D("ntrklen_chi2pid_RecoInel3","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);

	TH2D *ntruetrklen_chi2pid_BQ=new TH2D("ntruetrklen_chi2pid_BQ","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntruetrklen_chi2pid_BQ_inel=new TH2D("ntruetrklen_chi2pid_BQ_inel","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntruetrklen_chi2pid_BQ_el=new TH2D("ntruetrklen_chi2pid_BQ_el","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntruetrklen_chi2pid_BQ_midcosmic=new TH2D("ntruetrklen_chi2pid_BQ_midcosmic","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntruetrklen_chi2pid_BQ_midpi=new TH2D("ntruetrklen_chi2pid_BQ_midpi","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntruetrklen_chi2pid_BQ_midp=new TH2D("ntruetrklen_chi2pid_BQ_midp","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntruetrklen_chi2pid_BQ_midmu=new TH2D("ntruetrklen_chi2pid_BQ_midmu","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntruetrklen_chi2pid_BQ_mideg=new TH2D("ntruetrklen_chi2pid_BQ_mideg","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntruetrklen_chi2pid_BQ_midother=new TH2D("ntruetrklen_chi2pid_BQ_midother","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);

	TH2D *ntrklen_chi2pid_Pos=new TH2D("ntrklen_chi2pid_Pos","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntrklen_chi2pid_Pos_inel=new TH2D("ntrklen_chi2pid_Pos_inel","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntrklen_chi2pid_Pos_el=new TH2D("ntrklen_chi2pid_Pos_el","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);
	TH2D *ntrklen_chi2pid_Pos_midp=new TH2D("ntrklen_chi2pid_Pos_midp","", n_ntrklen, st_ntrklen, ed_ntrklen, n_chi2, st_chi2, ed_chi2);


	//1d pid
	TH1D *h1d_chi2pid_BQ=new TH1D("h1d_chi2pid_BQ","", n_chi2, st_chi2, ed_chi2);
	TH1D *h1d_chi2pid_RecoInel=new TH1D("h1d_chi2pid_RecoInel","", n_chi2, st_chi2, ed_chi2);
	TH1D *h1d_chi2pid_RecoInel2=new TH1D("h1d_chi2pid_RecoInel2","", n_chi2, st_chi2, ed_chi2);
	TH1D *h1d_chi2pid_RecoInel3=new TH1D("h1d_chi2pid_RecoInel3","", n_chi2, st_chi2, ed_chi2);
	TH1D *h1d_chi2pid_RecoStop=new TH1D("h1d_chi2pid_RecoStop","", n_chi2, st_chi2, ed_chi2);
	TH1D *h1d_chi2pid_RecoStop2=new TH1D("h1d_chi2pid_RecoStop2","", n_chi2, st_chi2, ed_chi2);

	TH1D *h1d_chi2pid_BQ_el=new TH1D(Form("h1d_chi2pid_BQ_el"), Form("el"), n_chi2, st_chi2, ed_chi2);
	TH1D *h1d_chi2pid_BQ_inel=new TH1D(Form("h1d_chi2pid_BQ_inel"), Form("inel"), n_chi2, st_chi2, ed_chi2);
	TH1D *h1d_chi2pid_BQ_midcosmic=new TH1D(Form("h1d_chi2pid_BQ_midcosmic"), Form("midcosmic"), n_chi2, st_chi2, ed_chi2);
	TH1D *h1d_chi2pid_BQ_midpi=new TH1D(Form("h1d_chi2pid_BQ_midpi"), Form("midpi"), n_chi2, st_chi2, ed_chi2);
	TH1D *h1d_chi2pid_BQ_midp=new TH1D(Form("h1d_chi2pid_BQ_midp"), Form("midp"), n_chi2, st_chi2, ed_chi2);
	TH1D *h1d_chi2pid_BQ_midmu=new TH1D(Form("h1d_chi2pid_BQ_midmu"), Form("midmu"), n_chi2, st_chi2, ed_chi2);
	TH1D *h1d_chi2pid_BQ_mideg=new TH1D(Form("h1d_chi2pid_BQ_mideg"), Form("mideg"), n_chi2, st_chi2, ed_chi2);
	TH1D *h1d_chi2pid_BQ_midother=new TH1D(Form("h1d_chi2pid_BQ_midother"), Form("midother"), n_chi2, st_chi2, ed_chi2);

	//1d ntrklen
	TH1D *h1d_ntrklen_BQ=new TH1D("h1d_ntrklen_BQ","", n_ntrklen, st_ntrklen, ed_ntrklen);
	TH1D *h1d_ntrklen_RecoInel=new TH1D("h1d_ntrklen_RecoInel","", n_ntrklen, st_ntrklen, ed_ntrklen);
	TH1D *h1d_ntrklen_RecoInel2=new TH1D("h1d_ntrklen_RecoInel2","", n_ntrklen, st_ntrklen, ed_ntrklen);
	TH1D *h1d_ntrklen_RecoInel3=new TH1D("h1d_ntrklen_RecoInel3","", n_ntrklen, st_ntrklen, ed_ntrklen);
	TH1D *h1d_ntrklen_RecoStop=new TH1D("h1d_ntrklen_RecoStop","", n_ntrklen, st_ntrklen, ed_ntrklen);
	TH1D *h1d_ntrklen_RecoStop2=new TH1D("h1d_ntrklen_RecoStop2","", n_ntrklen, st_ntrklen, ed_ntrklen);

	TH1D *h1d_ntrklen_BQ_el=new TH1D(Form("h1d_ntrklen_BQ_el"), Form("el"), n_ntrklen, st_ntrklen, ed_ntrklen);
	TH1D *h1d_ntrklen_BQ_inel=new TH1D(Form("h1d_ntrklen_BQ_inel"), Form("inel"), n_ntrklen, st_ntrklen, ed_ntrklen);
	TH1D *h1d_ntrklen_BQ_midcosmic=new TH1D(Form("h1d_ntrklen_BQ_midcosmic"), Form("midcosmic"), n_ntrklen, st_ntrklen, ed_ntrklen);
	TH1D *h1d_ntrklen_BQ_midpi=new TH1D(Form("h1d_ntrklen_BQ_midpi"), Form("midpi"), n_ntrklen, st_ntrklen, ed_ntrklen);
	TH1D *h1d_ntrklen_BQ_midp=new TH1D(Form("h1d_ntrklen_BQ_midp"), Form("midp"), n_ntrklen, st_ntrklen, ed_ntrklen);
	TH1D *h1d_ntrklen_BQ_midmu=new TH1D(Form("h1d_ntrklen_BQ_midmu"), Form("midmu"), n_ntrklen, st_ntrklen, ed_ntrklen);
	TH1D *h1d_ntrklen_BQ_mideg=new TH1D(Form("h1d_ntrklen_BQ_mideg"), Form("mideg"), n_ntrklen, st_ntrklen, ed_ntrklen);
	TH1D *h1d_ntrklen_BQ_midother=new TH1D(Form("h1d_ntrklen_BQ_midother"), Form("midother"), n_ntrklen, st_ntrklen, ed_ntrklen);

	//1d trklen
	int n_b=30;
	double b_min=0;
	double b_max=150;
	TH1D *h1d_trklen_RecoInel=new TH1D(Form("h1d_trklen_RecoInel"), Form("MC default"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_el=new TH1D(Form("h1d_trklen_RecoInel_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_inel=new TH1D(Form("h1d_trklen_RecoInel_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_midcosmic=new TH1D(Form("h1d_trklen_RecoInel_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_midpi=new TH1D(Form("h1d_trklen_RecoInel_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_midp=new TH1D(Form("h1d_trklen_RecoInel_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_midmu=new TH1D(Form("h1d_trklen_RecoInel_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_mideg=new TH1D(Form("h1d_trklen_RecoInel_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_midother=new TH1D(Form("h1d_trklen_RecoInel_midother"), Form("midother"), n_b, b_min, b_max);

	TH1D *h1d_trklen_RecoInel2=new TH1D(Form("h1d_trklen_RecoInel2"), Form("MC default"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel2_el=new TH1D(Form("h1d_trklen_RecoInel2_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel2_inel=new TH1D(Form("h1d_trklen_RecoInel2_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel2_midcosmic=new TH1D(Form("h1d_trklen_RecoInel2_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel2_midpi=new TH1D(Form("h1d_trklen_RecoInel2_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel2_midp=new TH1D(Form("h1d_trklen_RecoInel2_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel2_midmu=new TH1D(Form("h1d_trklen_RecoInel2_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel2_mideg=new TH1D(Form("h1d_trklen_RecoInel2_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel2_midother=new TH1D(Form("h1d_trklen_RecoInel2_midother"), Form("midother"), n_b, b_min, b_max);

	TH1D *h1d_trklen_RecoInel3=new TH1D(Form("h1d_trklen_RecoInel3"), Form("MC default"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel3_el=new TH1D(Form("h1d_trklen_RecoInel3_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel3_inel=new TH1D(Form("h1d_trklen_RecoInel3_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel3_midcosmic=new TH1D(Form("h1d_trklen_RecoInel3_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel3_midpi=new TH1D(Form("h1d_trklen_RecoInel3_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel3_midp=new TH1D(Form("h1d_trklen_RecoInel3_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel3_midmu=new TH1D(Form("h1d_trklen_RecoInel3_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel3_mideg=new TH1D(Form("h1d_trklen_RecoInel3_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel3_midother=new TH1D(Form("h1d_trklen_RecoInel3_midother"), Form("midother"), n_b, b_min, b_max);

	//other pid parameters
	TH1D *h1d_ntrklen_chi2pid_BQ=new TH1D("h1d_ntrklen_chi2pid_BQ","",1000,0,100);
	TH1D *h1d_ntrklen_chi2pid_BQ_inel=new TH1D("h1d_ntrklen_chi2pid_BQ_inel","",1000,0,100);
	TH1D *h1d_ntrklen_chi2pid_BQ_el=new TH1D("h1d_ntrklen_chi2pid_BQ_el","",1000,0,100);
	TH1D *h1d_ntrklen_chi2pid_RecoInel=new TH1D("h1d_ntrklen_chi2pid_RecoInel","",1000,0,100);

	TH1D *h1d_mediandedx_BQ=new TH1D("h1d_mediandedx_BQ","",100,0,10);
	TH1D *h1d_mediandedx_BQ_el=new TH1D("h1d_mediandedx_BQ_el","",100,0,10);
	TH1D *h1d_mediandedx_BQ_inel=new TH1D("h1d_mediandedx_BQ_inel","",100,0,10);
	TH1D *h1d_mediandedx_BQ_midp=new TH1D("h1d_mediandedx_BQ_midp","",100,0,10);
	TH1D *h1d_mediandedx_BQ_midcosmic=new TH1D("h1d_mediandedx_BQ_midcosmic","",100,0,10);
	TH1D *h1d_mediandedx_BQ_midpi=new TH1D("h1d_mediandedx_BQ_midpi","",100,0,10);
	TH1D *h1d_mediandedx_BQ_midmu=new TH1D("h1d_mediandedx_BQ_midmu","",100,0,10);
	TH1D *h1d_mediandedx_BQ_mideg=new TH1D("h1d_mediandedx_BQ_mideg","",100,0,10);
	TH1D *h1d_mediandedx_BQ_midother=new TH1D("h1d_mediandedx_BQ_midother","",100,0,10);
	TH1D *h1d_mediandedx_RecoInel=new TH1D("h1d_mediandedx_RecoInel","",100,0,10);

	TH1D *h1d_dEdL_BQ=new TH1D("h1d_dEdL_BQ","",100,0,10);
	TH1D *h1d_dEdL_BQ_el=new TH1D("h1d_dEdL_BQ_el","",100,0,10);
	TH1D *h1d_dEdL_BQ_inel=new TH1D("h1d_dEdL_BQ_inel","",100,0,10);
	TH1D *h1d_dEdL_BQ_midp=new TH1D("h1d_dEdL_BQ_midp","",100,0,10);
	TH1D *h1d_dEdL_BQ_midcosmic=new TH1D("h1d_dEdL_BQ_midcosmic","",100,0,10);
	TH1D *h1d_dEdL_BQ_midpi=new TH1D("h1d_dEdL_BQ_midpi","",100,0,10);
	TH1D *h1d_dEdL_BQ_midmu=new TH1D("h1d_dEdL_BQ_midmu","",100,0,10);
	TH1D *h1d_dEdL_BQ_mideg=new TH1D("h1d_dEdL_BQ_mideg","",100,0,10);
	TH1D *h1d_dEdL_BQ_midother=new TH1D("h1d_dEdL_BQ_midother","",100,0,10);
	TH1D *h1d_dEdL_RecoInel=new TH1D("h1d_dEdL_RecoInel","",100,0,10);

	//2D dists after Pos+CaloSize+PandoraSlice cut
	//cosTheta vs ntrklen
        int n_cosine=100;
        double cosine_min=0;
        double cosine_max=1.0;
	TH2D *ntrklen_cosineTheta_Pos_el=new TH2D("ntrklen_cosineTheta_Pos_el","", n_ntrklen, st_ntrklen, ed_ntrklen, n_cosine, cosine_min, cosine_max); 
	TH2D *ntrklen_cosineTheta_Pos_inel=new TH2D("ntrklen_cosineTheta_Pos_inel","", n_ntrklen, st_ntrklen, ed_ntrklen, n_cosine, cosine_min, cosine_max);
	TH2D *ntrklen_cosineTheta_Pos_midp=new TH2D("ntrklen_cosineTheta_Pos_midp","", n_ntrklen, st_ntrklen, ed_ntrklen, n_cosine, cosine_min, cosine_max);
	//cosTheta vs chi^2 PID
        TH2D *chi2pid_cosineTheta_Pos_el=new TH2D("chi2pid_cosineTheta_Pos_el","", n_chi2, st_chi2, ed_chi2, n_cosine, cosine_min, cosine_max);
        TH2D *chi2pid_cosineTheta_Pos_inel=new TH2D("chi2pid_cosineTheta_Pos_inel","", n_chi2, st_chi2, ed_chi2, n_cosine, cosine_min, cosine_max);
        TH2D *chi2pid_cosineTheta_Pos_midp=new TH2D("chi2pid_cosineTheta_Pos_midp","", n_chi2, st_chi2, ed_chi2, n_cosine, cosine_min, cosine_max);

	//cosTheta vs trklen
	TH2D *trklen_cosineTheta_Pos_el=new TH2D(Form("trklen_cosineTheta_Pos_el"), Form(""), n_b, b_min, b_max, n_cosine, cosine_min, cosine_max);
	TH2D *trklen_cosineTheta_Pos_inel=new TH2D(Form("trklen_cosineTheta_Pos_inel"), Form(""), n_b, b_min, b_max, n_cosine, cosine_min, cosine_max);
	TH2D *trklen_cosineTheta_Pos_midp=new TH2D(Form("trklen_cosineTheta_Pos_midp"), Form(""), n_b, b_min, b_max, n_cosine, cosine_min, cosine_max);




	//------------------------------------------------------------------------------------------------------------//

	for (Long64_t jentry=0; jentry<nentries;jentry++) { //main entry loop
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		isTestSample = true;
		if (ientry%2 == 0) isTestSample = false; //Divide MC sample by 2 parts: test+ufold

		//only select protons	
		if (beamtrackPdg!=pdg) continue; //only interested in protons
		n_tot++;


		//Event Selection Cut -- Part 1 ----------------------------------//
		bool IsBeamMatch=false; //if recostructed the right track (recoID=truthID)
		bool IsPandoraSlice=false; //pandora slice cut (can pandora reconstruct this track)
		bool IsCaloSize=false; //if calo size not empty
		bool IsIntersection=false; //if any track intersect with our reco track		
		if (primary_truth_Isbeammatched==1) IsBeamMatch=true;
		if (isprimarytrack==1&&isprimaryshower==0) IsPandoraSlice=true; 
		if (!primtrk_hitz->empty()) IsCaloSize=true;
		if (timeintersection->size()) IsIntersection=true;
		//----------------------------------------------------------------//

		//Truth label of Primarytrack_End ------------------------------------------------------------------------------------------------//
		bool IsPureInEL=false; //inel
		bool IsPureEL=false; //el

		if (strcmp(primary_truth_EndProcess->c_str(),"protonInelastic")==0) { IsPureInEL=true; }
		else { IsPureEL=true; }
		//--------------------------------------------------------------------------------------------------------------------------------//

		//Get true start/end point -----------------------------------------------------------------------//
		double true_endz=primary_truth_EndPosition_MC[2]; 
		double true_endy=primary_truth_EndPosition_MC[1]; 
		double true_endx=primary_truth_EndPosition_MC[0];

		double true_stz=primary_truth_StartPosition_MC[2];
		double true_sty=primary_truth_StartPosition_MC[1];
		double true_stx=primary_truth_StartPosition_MC[0];

		bool IsTrueEndOutside=false;
		if (true_endz<0.) IsTrueEndOutside=true;

		//Get reco info ----------------------------------------------------------------------------------//
		//Evt Classification -----------------------------------------------------------------------------//
		//signal -----------------------------//
		bool kinel=false;
		bool kel=false;
		//bool kmcs=false;
		if (IsBeamMatch) { //beam-match
			if (IsPureInEL) kinel=true;
			if (IsPureEL) kel=true;
			//if (IsPureMCS) kmcs=true;
		} //beam-match

		//background ------------------------------------------------------------------------//
		bool kMIDcosmic=false; //beam or cosmic
		bool kMIDpi=false; //+-pi
		bool kMIDp=false; //p
		bool kMIDmu=false; //mu
		bool kMIDeg=false; //e/gamma
		bool kMIDother=false; //other
		if (!IsBeamMatch) { //!beam-match
			if (primary_truth_byE_origin==2) { 
				kMIDcosmic=true;
			}
			else if (std::abs(primary_truth_byE_PDG)==211) {
				kMIDpi=true;
			}
			else if (primary_truth_byE_PDG==2212) {
				kMIDp=true;
			}
			else if (std::abs(primary_truth_byE_PDG)==13) {
				kMIDmu=true;
			}
			else if (std::abs(primary_truth_byE_PDG)==11 || primary_truth_byE_PDG==22) {
				kMIDeg=true;
			}
			else {
				kMIDother=true;
			}
		} //!beam-match	
		//cout<<"kMIDcosmic:"<<kMIDcosmic<<endl;
		//Evt Classification -----------------------------------------------------------------------------//

		//reco pos info & cut
		double reco_stx=-99, reco_sty=-99, reco_stz=-99;
		double reco_endx=-99, reco_endy=-99, reco_endz=-99;
		//double dx_reco_stx__bposx_ff=-999, dy_reco_sty__bposy_ff=-999, dz_reco_stz__bposz_ff=-999;
		bool IsPos=false;
		if (IsCaloSize) {
			reco_stx=primtrk_hitx->at(0); 
			reco_sty=primtrk_hity->at(0);
			reco_stz=primtrk_hitz->at(0);

			reco_endx=primtrk_hitx->at(primtrk_dedx->size()-1);	
			reco_endy=primtrk_hity->at(primtrk_dedx->size()-1);
			reco_endz=primtrk_hitz->at(primtrk_dedx->size()-1);

			double beam_dx=(reco_stx-mean_StartX)/sigma_StartX;
			double beam_dy=(reco_sty-mean_StartY)/sigma_StartY;
			double beam_dz=(reco_stz-mean_StartZ)/sigma_StartZ;
			double beam_dxy=sqrt(pow(beam_dx,2)+pow(beam_dy,2));	

			if (beam_dx>=dx_min&&beam_dx<=dx_max) { //dx
				if (beam_dy>=dy_min&&beam_dy<=dy_max) { //dy
					if (beam_dz>=dz_min&&beam_dz<=dz_max) { //dz
						if (beam_dxy>=dxy_min&&beam_dxy<=dxy_max) { //dxy
							IsPos=true;
						} //dxy
					} //dz
				} //dy
			} //dx

		}

		//cosine_theta/cut
		bool IsCosine=false;
		double cosine_beam_spec_primtrk=-999; 
		//cosine_beam_spec_primtrk=beamDirx_spec->at(0)*primaryStartDirection[0]+beamDiry_spec->at(0)*primaryStartDirection[1]+beamDirz_spec->at(0)*primaryStartDirection[2]; //cosine between beam_spec and primary trk direction(no SCE corr.)
		TVector3 dir;
		if (IsCaloSize) {	
			//trk direction after SCE corr.
			TVector3 pt0(primtrk_hitx->at(0), primtrk_hity->at(0), primtrk_hitz->at(0));
			TVector3 pt1(primtrk_hitx->at(-1+primtrk_hitx->size()), primtrk_hity->at(-1+primtrk_hity->size()), primtrk_hitz->at(-1+primtrk_hitz->size()));
			//TVector3 dir = pt1 - pt0;
			dir = pt1 - pt0;
			dir = dir.Unit();

			//beam direction
			//TVector3 beamdir(cos(beam_angleX_mc*TMath::Pi()/180), cos(beam_angleY_mc*TMath::Pi()/180), cos(beam_angleZ_mc*TMath::Pi()/180));
			TVector3 beamdir(beamDirx_spec->at(0),beamDiry_spec->at(0),beamDirz_spec->at(0));
			beamdir = beamdir.Unit();
			//beam_costh = dir.Dot(beamdir);
			cosine_beam_spec_primtrk=dir.Dot(beamdir);
		}

		if (cosine_beam_spec_primtrk<0) { cosine_beam_spec_primtrk=-1.*cosine_beam_spec_primtrk; }
		if (cosine_beam_spec_primtrk>cosine_beam_primtrk_min) { IsCosine=true; }

		//xy-cut
		bool IsXY=false;		
		double reco_stx_noSCE=0, reco_sty_noSCE=0, reco_stz_noSCE=0; //start-pos, before sce
		if (primaryEndPosition[2]>primaryStartPosition[2]) { //check if Pandora flip the sign
			reco_stx_noSCE=primaryStartPosition[0];
			reco_sty_noSCE=primaryStartPosition[1];
			reco_stz_noSCE=primaryStartPosition[2];
		} //check if Pandora flip the sign
		else {
			reco_stx_noSCE=primaryEndPosition[0];
			reco_sty_noSCE=primaryEndPosition[1];
			reco_stz_noSCE=primaryEndPosition[2];
		}
		if ((pow(((reco_stx_noSCE-mean_x)/dev_x),2)+pow(((reco_sty_noSCE-mean_y)/dev_y),2))<=1.) IsXY=true;

		//beam quality cut
		bool IsBQ=false;
		if (IsCosine&&IsPos) IsBQ=true;

		int index_reco_endz=0;
		double wid_reco_max=-9999;
		double range_reco=-999;
		vector<double> reco_trklen_accum;
		reco_trklen_accum.reserve(primtrk_hitz->size());
		double kereco_calo=0;
		double kereco_range=0;
		double kereco_range2=0;
		vector<double> EDept;
		if (IsCaloSize) { //if calo size not empty
			for (size_t h=0; h<primtrk_dedx->size(); ++h) { //loop over reco hits of a given track
				double hitx_reco=primtrk_hitx->at(h);
				double hity_reco=primtrk_hity->at(h);
				double hitz_reco=primtrk_hitz->at(h);
				double resrange_reco=primtrk_resrange->at(h);

				double dqdx=primtrk_dqdx->at(h);
				double pitch=primtrk_pitch->at(h);

				int wid_reco=primtrk_wid->at(-1+primtrk_wid->size()-h);
				double pt_reco=primtrk_pt->at(-1+primtrk_wid->size()-h);

				if (wid_reco>wid_reco_max) { 
					wid_reco_max=wid_reco;
					index_reco_endz=(int)-1+primtrk_wid->size()-h;
				}

				double cali_dedx=0.;
				cali_dedx=dedx_function_35ms(dqdx, hitx_reco, hity_reco, hitz_reco);
				EDept.push_back(cali_dedx*pitch);

				if (h==1) range_reco=0;
				if (h>=1) {
					range_reco += sqrt( pow(primtrk_hitx->at(h)-primtrk_hitx->at(h-1), 2)+
							pow(primtrk_hity->at(h)-primtrk_hity->at(h-1), 2)+
							pow(primtrk_hitz->at(h)-primtrk_hitz->at(h-1), 2) );
					reco_trklen_accum[h] = range_reco;
				}

				kereco_calo+=cali_dedx*pitch;
				kereco_range+=pitch*dedx_predict(resrange_reco);
				kereco_range2+=pitch*(double)gr_predict_dedx_resrange->Eval(resrange_reco);

			} //loop over reco hits of a given track
			//range_reco=primtrk_range->at(0);
		} //if calo size not empty


		//Reco stopping/Inel p cut
		bool IsRecoStop=false;
		bool IsRecoInEL=false;
		double mom_beam_spec=-99; mom_beam_spec=beamMomentum_spec->at(0);
		//double range_reco=-99; if (!primtrk_range->empty()) range_reco=primtrk_range->at(0); //reco primary trklen
		double csda_val_spec=csda_range_vs_mom_sm->Eval(mom_beam_spec);
		double norm_trklen=range_reco/csda_val_spec;

		if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) IsRecoStop=true;
		if ((range_reco/csda_val_spec)<min_norm_trklen_csda) IsRecoInEL=true;

		//kinetic energies
		double ke_beam_spec=p2ke(mom_beam_spec); //ke_beam_spec [GeV]
		double ke_beam_spec_MeV=1000.*ke_beam_spec; //ke_beam_spec [MeV]
		double ke_trklen=ke_vs_csda_range_sm->Eval(range_reco); //[unit: GeV]
		double ke_trklen_MeV=1000.*ke_trklen; //[unit: MeV]
		double ke_calo_MeV=0;


		//Get true trklen ---------------------------------------------------------------------------------------//
		double range_true=-999;
		int key_st = 0;
		double tmp_z = 9999;
		//vector<double> true_trklen_accum;
		//cout<<"beamtrk_z->size():"<<beamtrk_z->size()<<endl;
		//true_trklen_accum.reserve(beamtrk_z->size()); // initialize true_trklen_accum
		//cout<<"ck0"<<endl;
		for (int iz=0; iz<(int)beamtrk_z->size(); iz++) {
			if (abs(beamtrk_z->at(iz)) < tmp_z){
				tmp_z = abs(beamtrk_z->at(iz));
				key_st = iz; // find the point where the beam enters the TPC (find the smallest abs(Z))
			}
			//cout<<"ck0/"<<endl;
			//true_trklen_accum[iz] = 0.; // initialize true_trklen_accum
			//cout<<"ck0///"<<endl;
		}
		//cout<<"ck1"<<endl;
		for (int iz=key_st+1; iz<(int)beamtrk_z->size(); iz++){
			if (iz == key_st+1) range_true = 0;
			range_true += sqrt( pow(beamtrk_x->at(iz)-beamtrk_x->at(iz-1), 2)+
					pow(beamtrk_y->at(iz)-beamtrk_y->at(iz-1), 2)+	
					pow(beamtrk_z->at(iz)-beamtrk_z->at(iz-1), 2) );						    	
			//true_trklen_accum[iz] = range_true;
		}

		//cout<<"range_true:"<<range_true<<endl;
		//cout<<"key_st:"<<key_st<<endl;
		//for (size_t j=0; j<beamtrk_z->size(); ++j) { //MCParticle loop
		//cout<<"beamtrk_z["<<j<<"]:"<<beamtrk_z->at(j)<<" beamtrk_Eng["<<j<<"]:"<<beamtrk_Eng->at(j)<<" true_trklen_accum["<<j<<"]:"<<true_trklen_accum[j]<<endl;
		//} //MCParticle loop
		//Get reco info ----------------------------------------------------------------------------------//


		//PID parameter & cut
		double pid=-99; 
		double median_dedx=-99;
		//if (IsBQ&&IsCaloSize&&IsPandoraSlice) { //beam quality cut
		if (IsCaloSize&&IsPandoraSlice) { //beam quality cut
			//calo
			vector<double> trkdedx;
			vector<double> trkres;
			for (size_t h=0; h<primtrk_hitz->size(); ++h) { //loop over reco hits of a given track
				double hitx_reco=primtrk_hitx->at(h);
				double hity_reco=primtrk_hity->at(h);
				double hitz_reco=primtrk_hitz->at(h);
				double rr_reco=primtrk_resrange->at(h);

				double dqdx=primtrk_dqdx->at(h);
				double pitch=primtrk_pitch->at(h);

				int wid_reco=primtrk_wid->at(-1+primtrk_wid->size()-h);
				double pt_reco=primtrk_pt->at(-1+primtrk_wid->size()-h);

				double dedx_reco=0.;
				dedx_reco=dedx_function_35ms(dqdx, hitx_reco, hity_reco, hitz_reco);
				//ke_calo_MeV+=cali_dedx*pitch; //prod4a, pandoracalinoxyzt

				trkdedx.push_back(dedx_reco);
				trkres.push_back(rr_reco);
			} //loop over reco hits of a given track
			pid=chi2pid(trkdedx,trkres); //pid using stopping proton hypothesis

			median_dedx=TMath::Median(trkdedx.size(), &trkdedx.at(0));
		} //beam quality cut

		//define new reco inel cut ----------------------//
		//recoinel2
		bool IsRecoInEL2=false;
		bool IsRecoStop2=false;
		if (IsRecoInEL&&pid>7.5) IsRecoInEL2=true;
		if (IsRecoStop&&pid>10) IsRecoInEL2=true;

		if (IsRecoInEL&&pid<=7.5) IsRecoStop2=true;
		if (IsRecoStop&&pid<=10) IsRecoStop2=true;

		//recoinel3
		bool IsRecoInEL3=false;

		double x1_2d0=0.0328725;
		double y1_2d0=132.511;
		double x2_2d0=0.268832;
		double y2_2d0=55.3046;
		double m_2d0=(y1_2d0-y2_2d0)/(x1_2d0-x2_2d0);
		double b_2d0=y1_2d0-m_2d0*x1_2d0;

		double x1_2d=0.461096;
		double y1_2d=67.9097;
		double x2_2d=0.764785;
		double y2_2d=6.85399;
		double m_2d=(y1_2d-y2_2d)/(x1_2d-x2_2d);
		double b_2d=y1_2d-m_2d*x1_2d;

		if (IsRecoInEL) {
			if (norm_trklen<=0.245) {
				if (pid>=(b_2d0+m_2d0*norm_trklen)) {
					IsRecoInEL3=true;
				}
			}
			if (norm_trklen>0.245&&norm_trklen<=0.472) {
				if (pid>=60) {
					IsRecoInEL3=true;
				}
			}
			if (norm_trklen>0.472) {
				if (pid>=(b_2d+m_2d*norm_trklen)) {
					IsRecoInEL3=true;
				}
			}
		}
		if (IsRecoStop&&pid>10) IsRecoInEL3=true;
		//-----------------------------------------------//

		if (IsPos&&IsCaloSize&&IsPandoraSlice) { //IsPOs
			if (cosine_beam_spec_primtrk>0.5&&cosine_beam_spec_primtrk<0.9) { //misID:P rich(cosine>0.5<0.9)
				ntrklen_chi2pid_Pos->Fill(norm_trklen, pid);
				if (kinel) ntrklen_chi2pid_Pos_inel->Fill(norm_trklen, pid);
				if (kel) ntrklen_chi2pid_Pos_el->Fill(norm_trklen, pid);
				if (kMIDp) ntrklen_chi2pid_Pos_midp->Fill(norm_trklen, pid);
			} //misID:P rich(cosine>0.5<0.9)

			if (kel) { //true el
				ntrklen_cosineTheta_Pos_el->Fill(norm_trklen, cosine_beam_spec_primtrk);
				trklen_cosineTheta_Pos_el->Fill(range_reco, cosine_beam_spec_primtrk);	
				chi2pid_cosineTheta_Pos_el->Fill(pid, cosine_beam_spec_primtrk);
			} //true el

			if (kinel) { //true inel
				ntrklen_cosineTheta_Pos_inel->Fill(norm_trklen, cosine_beam_spec_primtrk);	
				trklen_cosineTheta_Pos_inel->Fill(range_reco, cosine_beam_spec_primtrk);	
				chi2pid_cosineTheta_Pos_inel->Fill(pid, cosine_beam_spec_primtrk);
			} //true inel

			if (kMIDp) { //midP
				ntrklen_cosineTheta_Pos_midp->Fill(norm_trklen, cosine_beam_spec_primtrk);	
				trklen_cosineTheta_Pos_midp->Fill(range_reco, cosine_beam_spec_primtrk);	
				chi2pid_cosineTheta_Pos_midp->Fill(pid, cosine_beam_spec_primtrk);
			} //midP
		} //IsPos



		if (IsBQ&&IsCaloSize&&IsPandoraSlice) { //beam quality cut
			ntrklen_chi2pid_BQ->Fill(norm_trklen, pid);
			Fill1DHist(h1d_chi2pid_BQ, pid);
			Fill1DHist(h1d_ntrklen_BQ, norm_trklen);

			Fill1DHist(h1d_ntrklen_chi2pid_BQ, pid/norm_trklen);
			Fill1DHist(h1d_mediandedx_BQ, median_dedx);
			Fill1DHist(h1d_dEdL_BQ,kereco_calo/range_reco);

			ntruetrklen_chi2pid_BQ->Fill(range_true/csda_val_spec, pid);

			if (kinel) { //true inel
				ntrklen_chi2pid_BQ_inel->Fill(norm_trklen, pid);
				Fill1DHist(h1d_chi2pid_BQ_inel, pid);
				Fill1DHist(h1d_ntrklen_BQ_inel,norm_trklen);

				Fill1DHist(h1d_ntrklen_chi2pid_BQ_inel, pid/norm_trklen);
				Fill1DHist(h1d_mediandedx_BQ_inel, median_dedx);
				Fill1DHist(h1d_dEdL_BQ_inel, kereco_calo/range_reco);

				ntruetrklen_chi2pid_BQ_inel->Fill(range_true/csda_val_spec, pid);	
			} //true inel

			if (kel) { //true el
				ntrklen_chi2pid_BQ_el->Fill(norm_trklen, pid);
				Fill1DHist(h1d_chi2pid_BQ_el, pid);
				Fill1DHist(h1d_ntrklen_BQ_el, norm_trklen);

				Fill1DHist(h1d_ntrklen_chi2pid_BQ_el, pid/norm_trklen);
				Fill1DHist(h1d_mediandedx_BQ_el, median_dedx);
				Fill1DHist(h1d_dEdL_BQ_el, kereco_calo/range_reco);

				ntruetrklen_chi2pid_BQ_el->Fill(range_true/csda_val_spec, pid);	
			} //true el

			if (kMIDcosmic) { //midcosmic
				ntrklen_chi2pid_BQ_midcosmic->Fill(norm_trklen, pid);
				Fill1DHist(h1d_chi2pid_BQ_midcosmic, pid);
				Fill1DHist(h1d_ntrklen_BQ_midcosmic, norm_trklen);

				Fill1DHist(h1d_mediandedx_BQ_midcosmic, median_dedx);
				Fill1DHist(h1d_dEdL_BQ_midcosmic, kereco_calo/range_reco);
				ntruetrklen_chi2pid_BQ_midcosmic->Fill(range_true/csda_val_spec, pid);
			} //midcosmic

			if (kMIDpi) { //midPi
				ntrklen_chi2pid_BQ_midpi->Fill(norm_trklen, pid);
				Fill1DHist(h1d_chi2pid_BQ_midpi, pid);
				Fill1DHist(h1d_ntrklen_BQ_midpi, norm_trklen);

				Fill1DHist(h1d_mediandedx_BQ_midpi, median_dedx);
				Fill1DHist(h1d_dEdL_BQ_midpi, kereco_calo/range_reco);
				ntruetrklen_chi2pid_BQ_midpi->Fill(range_true/csda_val_spec, pid);

			} //midPi

			if (kMIDp) { //midP
				ntrklen_chi2pid_BQ_midp->Fill(norm_trklen, pid);
				Fill1DHist(h1d_chi2pid_BQ_midp, pid);
				Fill1DHist(h1d_ntrklen_BQ_midp, norm_trklen);

				Fill1DHist(h1d_mediandedx_BQ_midp, median_dedx);
				Fill1DHist(h1d_dEdL_BQ_midp, kereco_calo/range_reco);
				ntruetrklen_chi2pid_BQ_midp->Fill(range_true/csda_val_spec, pid);
			} //midP

			if (kMIDmu) { //midmu
				ntrklen_chi2pid_BQ_midmu->Fill(norm_trklen, pid);
				Fill1DHist(h1d_chi2pid_BQ_midmu, pid);
				Fill1DHist(h1d_ntrklen_BQ_midmu, norm_trklen);

				Fill1DHist(h1d_mediandedx_BQ_midmu, median_dedx);
				Fill1DHist(h1d_dEdL_BQ_midmu, kereco_calo/range_reco);
				ntruetrklen_chi2pid_BQ_midmu->Fill(range_true/csda_val_spec, pid);
			} //midmu

			if (kMIDeg) { //mideg
				ntrklen_chi2pid_BQ_mideg->Fill(norm_trklen, pid);
				Fill1DHist(h1d_chi2pid_BQ_mideg, pid);
				Fill1DHist(h1d_ntrklen_BQ_mideg, norm_trklen);

				Fill1DHist(h1d_mediandedx_BQ_mideg, median_dedx);
				Fill1DHist(h1d_dEdL_BQ_mideg, kereco_calo/range_reco);
				ntruetrklen_chi2pid_BQ_mideg->Fill(range_true/csda_val_spec, pid);
			} //mideg
			if (kMIDother) { //midother
				ntrklen_chi2pid_BQ_midother->Fill(norm_trklen, pid);
				Fill1DHist(h1d_chi2pid_BQ_midother, pid);
				Fill1DHist(h1d_ntrklen_BQ_midother, norm_trklen);

				Fill1DHist(h1d_mediandedx_BQ_midother, median_dedx);
				Fill1DHist(h1d_dEdL_BQ_midother, kereco_calo/range_reco);
				ntruetrklen_chi2pid_BQ_midother->Fill(range_true/csda_val_spec, pid);
			} //midother

			if (IsRecoInEL) { //reco inel
				ntrklen_chi2pid_RecoInel->Fill(norm_trklen, pid);
				Fill1DHist(h1d_chi2pid_RecoInel, pid);
				Fill1DHist(h1d_ntrklen_RecoInel, norm_trklen);

				Fill1DHist(h1d_ntrklen_chi2pid_RecoInel, pid/norm_trklen);
				Fill1DHist(h1d_mediandedx_RecoInel, median_dedx);
				Fill1DHist(h1d_dEdL_RecoInel, kereco_calo/range_reco);

				Fill1DHist(h1d_trklen_RecoInel,range_reco);
				if (kel) Fill1DHist(h1d_trklen_RecoInel_el, range_reco);
				if (kinel) Fill1DHist(h1d_trklen_RecoInel_inel, range_reco);
				if (kMIDcosmic) Fill1DHist(h1d_trklen_RecoInel_midcosmic, range_reco);
				if (kMIDpi) Fill1DHist(h1d_trklen_RecoInel_midpi, range_reco);
				if (kMIDp) Fill1DHist(h1d_trklen_RecoInel_midp, range_reco);
				if (kMIDmu) Fill1DHist(h1d_trklen_RecoInel_midmu, range_reco);
				if (kMIDeg) Fill1DHist(h1d_trklen_RecoInel_mideg, range_reco);
				if (kMIDother) Fill1DHist(h1d_trklen_RecoInel_midother, range_reco);

			} //reco inel
			if (IsRecoInEL2) { //reco inel2
				ntrklen_chi2pid_RecoInel2->Fill(norm_trklen, pid);
				Fill1DHist(h1d_chi2pid_RecoInel2, pid);
				Fill1DHist(h1d_ntrklen_RecoInel2, norm_trklen);

				Fill1DHist(h1d_trklen_RecoInel2,range_reco);
				if (kel) Fill1DHist(h1d_trklen_RecoInel2_el, range_reco);
				if (kinel) Fill1DHist(h1d_trklen_RecoInel2_inel, range_reco);
				if (kMIDcosmic) Fill1DHist(h1d_trklen_RecoInel2_midcosmic, range_reco);
				if (kMIDpi) Fill1DHist(h1d_trklen_RecoInel2_midpi, range_reco);
				if (kMIDp) Fill1DHist(h1d_trklen_RecoInel2_midp, range_reco);
				if (kMIDmu) Fill1DHist(h1d_trklen_RecoInel2_midmu, range_reco);
				if (kMIDeg) Fill1DHist(h1d_trklen_RecoInel2_mideg, range_reco);
				if (kMIDother) Fill1DHist(h1d_trklen_RecoInel2_midother, range_reco);

			} //reco inel2

			if (IsRecoInEL3) { //reco inel3
				ntrklen_chi2pid_RecoInel3->Fill(norm_trklen, pid);
				Fill1DHist(h1d_chi2pid_RecoInel3, pid);
				Fill1DHist(h1d_ntrklen_RecoInel3, norm_trklen);

				Fill1DHist(h1d_trklen_RecoInel3,range_reco);
				if (kel) Fill1DHist(h1d_trklen_RecoInel3_el, range_reco);
				if (kinel) Fill1DHist(h1d_trklen_RecoInel3_inel, range_reco);
				if (kMIDcosmic) Fill1DHist(h1d_trklen_RecoInel3_midcosmic, range_reco);
				if (kMIDpi) Fill1DHist(h1d_trklen_RecoInel3_midpi, range_reco);
				if (kMIDp) Fill1DHist(h1d_trklen_RecoInel3_midp, range_reco);
				if (kMIDmu) Fill1DHist(h1d_trklen_RecoInel3_midmu, range_reco);
				if (kMIDeg) Fill1DHist(h1d_trklen_RecoInel3_mideg, range_reco);
				if (kMIDother) Fill1DHist(h1d_trklen_RecoInel3_midother, range_reco);

			} //reco inel3


			if (IsRecoStop) { //reco stop
				Fill1DHist(h1d_chi2pid_RecoStop, pid);
				Fill1DHist(h1d_ntrklen_RecoStop, norm_trklen);
			} //reco stop


			if (IsRecoStop2) { //reco stop 2
				Fill1DHist(h1d_chi2pid_RecoStop2, pid);
				Fill1DHist(h1d_ntrklen_RecoStop2, norm_trklen);

			} //reco stop 2

		} //beam quality cut



		//counting -----------------------------------------------------------//
		if (IsPandoraSlice) n_pan_tot++;
		if (IsPandoraSlice&&IsCaloSize) n_calsz_tot++;
		if (IsPandoraSlice&&IsCaloSize&&IsBQ) n_bq_tot++;
		if (IsPandoraSlice&&IsCaloSize&&IsBQ&&IsRecoInEL) n_recoinel_tot++;
		if (IsPandoraSlice&&IsCaloSize&&IsBQ&&IsRecoInEL2) n_recoinel2_tot++;
		if (IsPandoraSlice&&IsCaloSize&&IsBQ&&IsRecoInEL3) n_recoinel3_tot++;


		//[0]pure inel
		if (kinel) { //pure inel
			n_inel++;
			if (IsPandoraSlice) { //pandora
				n_inel_pan++;
				if (IsCaloSize) { //calosz
					n_inel_calsz++;
					if (IsBQ) { //bq
						n_inel_bq++;
						if (IsRecoInEL) { //reco inel
							n_inel_recoinel++;
						} //reco inel
						if (IsRecoInEL2) { //reco inel2
							n_inel_recoinel2++;
						} //reco inel2
						if (IsRecoInEL3) { //reco inel3
							n_inel_recoinel3++;
						} //reco inel3
					} //bq
				} //calosz
			} //pandora
		} //pure inel		

		//[1]pure el
		if (kel) { //pure el
			n_el++;
			if (IsPandoraSlice) { //pandora
				n_el_pan++;
				if (IsCaloSize) { //calosz
					n_el_calsz++;
					if (IsBQ) { //bq
						n_el_bq++;
						if (IsRecoInEL) { //reco inel
							n_el_recoinel++;
						} //reco inel
						if (IsRecoInEL2) { //reco inel2
							n_el_recoinel2++;
						} //reco inel2
						if (IsRecoInEL3) { //reco inel3
							n_el_recoinel3++;
						} //reco inel3
					} //bq
				} //calosz
			} //pandora
		} //pure el		

		//[3]MID:Cosmic
		if (kMIDcosmic) { //mid:cosmic
			n_midcosmic++;
			if (IsPandoraSlice) { //pandora
				n_midcosmic_pan++;
				if (IsCaloSize) { //calosz
					n_midcosmic_calsz++;
					if (IsBQ) { //bq
						n_midcosmic_bq++;
						if (IsRecoInEL) { //reco inel
							n_midcosmic_recoinel++;
						} //reco inel
						if (IsRecoInEL2) { //reco inel2
							n_midcosmic_recoinel2++;
						} //reco inel2
						if (IsRecoInEL3) { //reco inel3
							n_midcosmic_recoinel3++;
						} //reco inel3
					} //bq
				} //calosz
			} //pandora
		} //mid:cosmic

		//[4]MID:midpi
		if (kMIDpi) { //mid:pi
			n_midpi++;
			if (IsPandoraSlice) { //pandora
				n_midpi_pan++;
				if (IsCaloSize) { //calosz
					n_midpi_calsz++;
					if (IsBQ) { //bq
						n_midpi_bq++;
						if (IsRecoInEL) { //reco inel
							n_midpi_recoinel++;
						} //reco inel
						if (IsRecoInEL2) { //reco inel2
							n_midpi_recoinel2++;
						} //reco inel2
						if (IsRecoInEL3) { //reco inel3
							n_midpi_recoinel3++;
						} //reco inel3
					} //bq
				} //calosz
			} //pandora
		} //mid:pi

		//[5]MID:midp
		if (kMIDp) { //mid:p
			n_midp++;
			if (IsPandoraSlice) { //pandora
				n_midp_pan++;
				if (IsCaloSize) { //calosz
					n_midp_calsz++;
					if (IsBQ) { //bq
						n_midp_bq++;
						if (IsRecoInEL) { //reco inel
							n_midp_recoinel++;
						} //reco inel
						if (IsRecoInEL2) { //reco inel2
							n_midp_recoinel2++;
						} //reco inel2
						if (IsRecoInEL3) { //reco inel3
							n_midp_recoinel3++;
						} //reco inel3
					} //bq
				} //calosz
			} //pandora
		} //mid:p

		//[6]MID:midmu
		if (kMIDmu) { //mid:mu
			n_midmu++;
			if (IsPandoraSlice) { //pandora
				n_midmu_pan++;
				if (IsCaloSize) { //calosz
					n_midmu_calsz++;
					if (IsBQ) { //bq
						n_midmu_bq++;
						if (IsRecoInEL) { //reco inel
							n_midmu_recoinel++;
						} //reco inel
						if (IsRecoInEL2) { //reco inel2
							n_midmu_recoinel2++;
						} //reco inel2
						if (IsRecoInEL3) { //reco inel3
							n_midmu_recoinel3++;
						} //reco inel3
					} //bq
				} //calosz
			} //pandora
		} //mid:mu

		//[7]MID:mideg
		if (kMIDeg) { //mid:eg
			n_mideg++;
			if (IsPandoraSlice) { //pandora
				n_mideg_pan++;
				if (IsCaloSize) { //calosz
					n_mideg_calsz++;
					if (IsBQ) { //bq
						n_mideg_bq++;
						if (IsRecoInEL) { //reco inel
							n_mideg_recoinel++;
						} //reco inel
						if (IsRecoInEL2) { //reco inel2
							n_mideg_recoinel2++;
						} //reco inel2
						if (IsRecoInEL3) { //reco inel3
							n_mideg_recoinel3++;
						} //reco inel3
					} //bq
				} //calosz
			} //pandora
		} //mid:eg

		//[8]MID:midother
		if (kMIDother) { //mid:other
			n_midother++;
			if (IsPandoraSlice) { //pandora
				n_midother_pan++;
				if (IsCaloSize) { //calosz
					n_midother_calsz++;
					if (IsBQ) { //bq
						n_midother_bq++;
						if (IsRecoInEL) { //reco inel
							n_midother_recoinel++;
						} //reco inel
						if (IsRecoInEL2) { //reco inel2
							n_midother_recoinel2++;
						} //reco inel2
						if (IsRecoInEL3) { //reco inel3
							n_midother_recoinel3++;
						} //reco inel3
					} //bq
				} //calosz
			} //pandora
		} //mid:other





	} //main entry loop

	//counting -- summary -----------------------------------------------------//
	cout<<"\nn_tot:"<<n_tot<<endl;
	cout<<"n_el:"<<n_el<<endl;
	cout<<"n_inel:"<<n_inel<<endl;
	cout<<"n_midcosmic:"<<n_midcosmic<<endl;
	cout<<"n_midpi:"<<n_midpi<<endl;
	cout<<"n_midp:"<<n_midp<<endl;
	cout<<"n_midmu:"<<n_midmu<<endl;
	cout<<"n_mideg:"<<n_mideg<<endl;
	cout<<"n_midother"<<n_midother<<endl;
	cout<<"n_diff:"<<n_tot-(n_el+n_inel+n_midcosmic+n_midpi+n_midp+n_midmu+n_mideg+n_midother)<<endl;

	cout<<"\nn_pan_tot:"<<n_pan_tot<<endl;
	cout<<"n_el_pan:"<<n_el_pan<<endl;
	cout<<"n_inel_pan:"<<n_inel_pan<<endl;
	cout<<"n_midcosmic_pan:"<<n_midcosmic_pan<<endl;
	cout<<"n_midpi_pan:"<<n_midpi_pan<<endl;
	cout<<"n_midp_pan:"<<n_midp_pan<<endl;
	cout<<"n_midmu_pan:"<<n_midmu_pan<<endl;
	cout<<"n_mideg_pan:"<<n_mideg_pan<<endl;
	cout<<"n_midother_pan"<<n_midother_pan<<endl;
	cout<<"n_diff_pan:"<<n_pan_tot-(n_el_pan+n_inel_pan+n_midcosmic_pan+n_midpi_pan+n_midp_pan+n_midmu_pan+n_mideg_pan+n_midother_pan)<<endl;

	cout<<"\nn_calsz_tot:"<<n_calsz_tot<<endl;
	cout<<"n_el_calsz:"<<n_el_calsz<<endl;
	cout<<"n_inel_calsz:"<<n_inel_calsz<<endl;
	cout<<"n_midcosmic_calsz:"<<n_midcosmic_calsz<<endl;
	cout<<"n_midpi_calsz:"<<n_midpi_calsz<<endl;
	cout<<"n_midp_calsz:"<<n_midp_calsz<<endl;
	cout<<"n_midmu_calsz:"<<n_midmu_calsz<<endl;
	cout<<"n_mideg_calsz:"<<n_mideg_calsz<<endl;
	cout<<"n_midother_calsz"<<n_midother_calsz<<endl;
	cout<<"n_diff_calsz:"<<n_calsz_tot-(n_el_calsz+n_inel_calsz+n_midcosmic_calsz+n_midpi_calsz+n_midp_calsz+n_midmu_calsz+n_mideg_calsz+n_midother_calsz)<<endl;

	cout<<"\nn_bq_tot:"<<n_bq_tot<<endl;
	cout<<"n_el_bq:"<<n_el_bq<<endl;
	cout<<"n_inel_bq:"<<n_inel_bq<<endl;
	cout<<"n_midcosmic_bq:"<<n_midcosmic_bq<<endl;
	cout<<"n_midpi_bq:"<<n_midpi_bq<<endl;
	cout<<"n_midp_bq:"<<n_midp_bq<<endl;
	cout<<"n_midmu_bq:"<<n_midmu_bq<<endl;
	cout<<"n_mideg_bq:"<<n_mideg_bq<<endl;
	cout<<"n_midother_bq"<<n_midother_bq<<endl;
	cout<<"n_diff_bq:"<<n_bq_tot-(n_el_bq+n_inel_bq+n_midcosmic_bq+n_midpi_bq+n_midp_bq+n_midmu_bq+n_mideg_bq+n_midother_bq)<<endl;

	cout<<"\nn_recoinel_tot:"<<n_recoinel_tot<<endl;
	cout<<"n_el_recoinel:"<<n_el_recoinel<<endl;
	cout<<"n_inel_recoinel:"<<n_inel_recoinel<<endl;
	cout<<"n_midcosmic_recoinel:"<<n_midcosmic_recoinel<<endl;
	cout<<"n_midpi_recoinel:"<<n_midpi_recoinel<<endl;
	cout<<"n_midp_recoinel:"<<n_midp_recoinel<<endl;
	cout<<"n_midmu_recoinel:"<<n_midmu_recoinel<<endl;
	cout<<"n_mideg_recoinel:"<<n_mideg_recoinel<<endl;
	cout<<"n_midother_recoinel"<<n_midother_recoinel<<endl;
	cout<<"n_diff_recoinel:"<<n_recoinel_tot-(n_el_recoinel+n_inel_recoinel+n_midcosmic_recoinel+n_midpi_recoinel+n_midp_recoinel+n_midmu_recoinel+n_mideg_recoinel+n_midother_recoinel)<<endl;

	cout<<"\nn_recoinel2_tot:"<<n_recoinel2_tot<<endl;
	cout<<"n_el_recoinel2:"<<n_el_recoinel2<<endl;
	cout<<"n_inel_recoinel2:"<<n_inel_recoinel2<<endl;
	cout<<"n_midcosmic_recoinel2:"<<n_midcosmic_recoinel2<<endl;
	cout<<"n_midpi_recoinel2:"<<n_midpi_recoinel2<<endl;
	cout<<"n_midp_recoinel2:"<<n_midp_recoinel2<<endl;
	cout<<"n_midmu_recoinel2:"<<n_midmu_recoinel2<<endl;
	cout<<"n_mideg_recoinel2:"<<n_mideg_recoinel2<<endl;
	cout<<"n_midother_recoinel2"<<n_midother_recoinel2<<endl;
	cout<<"n_diff_recoinel2:"<<n_recoinel2_tot-(n_el_recoinel2+n_inel_recoinel2+n_midcosmic_recoinel2+n_midpi_recoinel2+n_midp_recoinel2+n_midmu_recoinel2+n_mideg_recoinel2+n_midother_recoinel2)<<endl;

	cout<<"\nn_recoinel3_tot:"<<n_recoinel3_tot<<endl;
	cout<<"n_el_recoinel3:"<<n_el_recoinel3<<endl;
	cout<<"n_inel_recoinel3:"<<n_inel_recoinel3<<endl;
	cout<<"n_midcosmic_recoinel3:"<<n_midcosmic_recoinel3<<endl;
	cout<<"n_midpi_recoinel3:"<<n_midpi_recoinel3<<endl;
	cout<<"n_midp_recoinel3:"<<n_midp_recoinel3<<endl;
	cout<<"n_midmu_recoinel3:"<<n_midmu_recoinel3<<endl;
	cout<<"n_mideg_recoinel3:"<<n_mideg_recoinel3<<endl;
	cout<<"n_midother_recoinel3"<<n_midother_recoinel3<<endl;
	cout<<"n_diff_recoinel3:"<<n_recoinel3_tot-(n_el_recoinel3+n_inel_recoinel3+n_midcosmic_recoinel3+n_midpi_recoinel3+n_midp_recoinel3+n_midmu_recoinel3+n_mideg_recoinel3+n_midother_recoinel3)<<endl;


	//save results ---------------------------------------------------------//
	TFile *fout = new TFile("mc_proton_studyRecoInelCut.root","RECREATE");
	ntrklen_chi2pid_BQ->Write();
	ntrklen_chi2pid_BQ_inel->Write();
	ntrklen_chi2pid_BQ_el->Write();
	ntrklen_chi2pid_BQ_midcosmic->Write();
	ntrklen_chi2pid_BQ_midpi->Write();
	ntrklen_chi2pid_BQ_midp->Write();
	ntrklen_chi2pid_BQ_midmu->Write();
	ntrklen_chi2pid_BQ_mideg->Write();
	ntrklen_chi2pid_BQ_midother->Write();

	ntrklen_chi2pid_RecoInel->Write();
	ntrklen_chi2pid_RecoInel2->Write();
	ntrklen_chi2pid_RecoInel3->Write();

	ntruetrklen_chi2pid_BQ->Write();
	ntruetrklen_chi2pid_BQ_inel->Write();
	ntruetrklen_chi2pid_BQ_el->Write();
	ntruetrklen_chi2pid_BQ_midcosmic->Write();
	ntruetrklen_chi2pid_BQ_midpi->Write();
	ntruetrklen_chi2pid_BQ_midp->Write();
	ntruetrklen_chi2pid_BQ_midmu->Write();
	ntruetrklen_chi2pid_BQ_mideg->Write();
	ntruetrklen_chi2pid_BQ_midother->Write();

	h1d_chi2pid_BQ->Write();
	h1d_chi2pid_RecoInel->Write();
	h1d_chi2pid_RecoInel2->Write();
	h1d_chi2pid_RecoInel3->Write();
	h1d_chi2pid_RecoStop->Write();
	h1d_chi2pid_RecoStop2->Write();

	h1d_chi2pid_BQ_el->Write();
	h1d_chi2pid_BQ_inel->Write();
	h1d_chi2pid_BQ_midcosmic->Write();
	h1d_chi2pid_BQ_midpi->Write();
	h1d_chi2pid_BQ_midp->Write();
	h1d_chi2pid_BQ_midmu->Write();
	h1d_chi2pid_BQ_mideg->Write();
	h1d_chi2pid_BQ_midother->Write();

	h1d_ntrklen_BQ->Write();
	h1d_ntrklen_RecoInel->Write();
	h1d_ntrklen_RecoInel2->Write();
	h1d_ntrklen_RecoInel3->Write();

	h1d_ntrklen_RecoStop->Write();
	h1d_ntrklen_RecoStop2->Write();

	h1d_ntrklen_BQ_el->Write();
	h1d_ntrklen_BQ_inel->Write();
	h1d_ntrklen_BQ_midcosmic->Write();
	h1d_ntrklen_BQ_midpi->Write();
	h1d_ntrklen_BQ_midp->Write();
	h1d_ntrklen_BQ_midmu->Write();
	h1d_ntrklen_BQ_mideg->Write();
	h1d_ntrklen_BQ_midother->Write();

	h1d_trklen_RecoInel->Write();
	h1d_trklen_RecoInel_el->Write();
	h1d_trklen_RecoInel_inel->Write();
	h1d_trklen_RecoInel_midcosmic->Write();
	h1d_trklen_RecoInel_midpi->Write();
	h1d_trklen_RecoInel_midp->Write();
	h1d_trklen_RecoInel_midmu->Write();
	h1d_trklen_RecoInel_mideg->Write();
	h1d_trklen_RecoInel_midother->Write();


	h1d_trklen_RecoInel2->Write();
	h1d_trklen_RecoInel2_el->Write();
	h1d_trklen_RecoInel2_inel->Write();
	h1d_trklen_RecoInel2_midcosmic->Write();
	h1d_trklen_RecoInel2_midpi->Write();
	h1d_trklen_RecoInel2_midp->Write();
	h1d_trklen_RecoInel2_midmu->Write();
	h1d_trklen_RecoInel2_mideg->Write();
	h1d_trklen_RecoInel2_midother->Write();

	h1d_trklen_RecoInel3->Write();
	h1d_trklen_RecoInel3_el->Write();
	h1d_trklen_RecoInel3_inel->Write();
	h1d_trklen_RecoInel3_midcosmic->Write();
	h1d_trklen_RecoInel3_midpi->Write();
	h1d_trklen_RecoInel3_midp->Write();
	h1d_trklen_RecoInel3_midmu->Write();
	h1d_trklen_RecoInel3_mideg->Write();
	h1d_trklen_RecoInel3_midother->Write();

	h1d_ntrklen_chi2pid_BQ->Write();
	h1d_ntrklen_chi2pid_BQ_inel->Write();
	h1d_ntrklen_chi2pid_BQ_el->Write();
	h1d_ntrklen_chi2pid_RecoInel->Write();


	h1d_mediandedx_BQ->Write();
	h1d_mediandedx_BQ_el->Write();
	h1d_mediandedx_BQ_inel->Write();
	h1d_mediandedx_BQ_midp->Write();
	h1d_mediandedx_BQ_midcosmic->Write();
	h1d_mediandedx_BQ_midpi->Write();
	h1d_mediandedx_BQ_midmu->Write();
	h1d_mediandedx_BQ_mideg->Write();
	h1d_mediandedx_BQ_midother->Write();
	h1d_mediandedx_RecoInel->Write();

	h1d_dEdL_BQ->Write();
	h1d_dEdL_BQ_el->Write();
	h1d_dEdL_BQ_inel->Write();
	h1d_dEdL_BQ_midp->Write();
	h1d_dEdL_BQ_midcosmic->Write();
	h1d_dEdL_BQ_midpi->Write();
	h1d_dEdL_BQ_midmu->Write();
	h1d_dEdL_BQ_mideg->Write();
	h1d_dEdL_BQ_midother->Write();
	h1d_dEdL_RecoInel->Write();

	ntrklen_chi2pid_Pos->Write();
	ntrklen_chi2pid_Pos_inel->Write();
	ntrklen_chi2pid_Pos_el->Write();
	ntrklen_chi2pid_Pos_midp->Write();

	ntrklen_cosineTheta_Pos_el->Write();
	ntrklen_cosineTheta_Pos_inel->Write();
	ntrklen_cosineTheta_Pos_midp->Write();

	trklen_cosineTheta_Pos_el->Write();
	trklen_cosineTheta_Pos_inel->Write();
	trklen_cosineTheta_Pos_midp->Write();

	chi2pid_cosineTheta_Pos_el->Write();
	chi2pid_cosineTheta_Pos_inel->Write();
	chi2pid_cosineTheta_Pos_midp->Write();

	fout->Close();



}
