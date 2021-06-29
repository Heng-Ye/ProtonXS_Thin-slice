
double p2ke(double p) {
	double ke=m_proton*(-1.+sqrt(1.+(p*p)/(m_proton*m_proton)));
	return ke;
}

double ke2p(double ke) { //input ke unit: GeV
	double p=m_proton*sqrt(-1+pow((1+ke/m_proton),2));
	return p;
}

bool myComparison(const pair<double,int> &a,const pair<double,int> &b) {
	return a.first<b.first;
}

//Read file of dedx versus kinetic energy -----------------------//
TString conv_path="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/conversion/";
TFile *fke_dedx=new TFile(Form("%sproton_dedx_ke_MeV.root", conv_path.Data()));
TGraph *dedx_vs_ke_sm=(TGraph *)fke_dedx->Get("dedx_vs_ke_sm");

//Read file of csda range versus momentum --------------------------------------//
TFile *fmom_csda=new TFile(Form("%sproton_mom_csda_converter.root", conv_path.Data()));
TGraph *csda_range_vs_mom_sm=(TGraph *)fmom_csda->Get("csda_range_vs_mom_sm");
TGraph *mom_vs_csda_range_sm=(TGraph *)fmom_csda->Get("mom_vs_csda_range_sm");

//Read file of csda range versus kinetic energy---------------------------------//
TFile *fke_csda=new TFile(Form("%sproton_ke_csda_converter_reduction.root", conv_path.Data()));
TGraph *csda_range_vs_ke_sm=(TGraph *)fke_csda->Get("csda_range_vs_ke_sm");
TGraph *ke_vs_csda_range_sm=(TGraph *)fke_csda->Get("ke_vs_csda_range_sm_rd");

Double_t fitg(Double_t* x,Double_t *par) {
	double m=par[0];
	double s=par[1];
	double n=1;

	double g=n*TMath::Exp(-(x[0]-m)*(x[0]-m)/(2*s*s));
	//Double_t g=n/(s*sqrt(2*3.14159))*TMath::Exp(-(x[0]-m)*(x[0]-m)/(2*s*s));
	return g;
}

Double_t govg(Double_t* x,Double_t *par) {
	//g1
	double m1=par[0];
	double s1=par[1];
	//double n1=1.;
	double g1=-(x[0]-m1)*(x[0]-m1)/(2*s1*s1);

	//g2
	double m2=par[2];
	double s2=par[3];
	//double n2=1.;
	double g2=-(x[0]-m2)*(x[0]-m2)/(2*s2*s2);

	//g2/g1
	double g_ov_g=0; 
	g_ov_g=TMath::Exp(g2-g1);
	if (m1==m2&&s1==s2) g_ov_g=1;

	return g_ov_g;
}

double cutAPA3_Z = 226.;
bool endAPA3(double reco_beam_endZ){
	return(reco_beam_endZ < cutAPA3_Z);
}

Double_t dedx_predict(double rr) {
	double a=17.;
	double b=-0.42;

	return a*pow(rr,b);
}

//read the predicted dE/dx vs residual range
TFile *fdedx_rr=new TFile(Form("%sdedx_rr.root",conv_path.Data()));
TGraph *gr_predict_dedx_resrange=(TGraph *)fdedx_rr->Get("gr_predict_dedx_resrange");
TGraph *gr_wq_dedx_resrange=(TGraph *)fdedx_rr->Get("gr_wq_dedx_resrange");


