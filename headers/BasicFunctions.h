#include "TMath.h"
#include "betheBloch.h"

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

//Read file of dedx versus kinetic energy -----------------------------------------------------//
TString conv_path="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/conversion/";
/*
   TFile *fke_dedx=new TFile(Form("%sproton_dedx_ke_MeV.root", conv_path.Data()));
   TGraph *dedx_vs_ke_sm=(TGraph *)fke_dedx->Get("dedx_vs_ke_sm");

//Load dE/dx vs KE from NIST data base
//TFile *fKE_dEdx=new TFile(Form("%sproton_ydedx_xke.root", conv_path.Data()));
//TGraph *dEdx_vs_KE_sm=(TGraph *)fKE_dEdx->Get("dEdx_vs_KE_sm"); //x:KE in MeV; y:dE/dx in MeV/cm
*/


//Read file of csda range versus momentum --------------------------------------//
TFile *fmom_csda=new TFile(Form("%sproton_mom_csda_converter.root", conv_path.Data()));
TGraph *csda_range_vs_mom_sm=(TGraph *)fmom_csda->Get("csda_range_vs_mom_sm");
TGraph *mom_vs_csda_range_sm=(TGraph *)fmom_csda->Get("mom_vs_csda_range_sm");

//Read file of csda range versus kinetic energy---------------------------------//
TFile *fke_csda=new TFile(Form("%sproton_ke_csda_converter_reduction.root", conv_path.Data()));
TGraph *csda_range_vs_ke_sm=(TGraph *)fke_csda->Get("csda_range_vs_ke_sm");
TGraph *ke_vs_csda_range_sm=(TGraph *)fke_csda->Get("ke_vs_csda_range_sm_rd");

//Function to convert trklen to Edept -----------------------------------------------------//
/*
   void hist_NIST(double E_init, TH1D* h_bethe){
   for(int i=1; i <= h_bethe->GetNbinsX(); i++){
   h_bethe->SetBinContent( i, dEdx_vs_KE_sm->Eval(E_init));
   h_bethe->SetBinError(i, 0.001 );
   E_init = E_init - dEdx_vs_KE_sm->Eval(E_init);
   if(E_init <= 0) return;
   };
   };
   */

/*
   void hist_bethe_mean_distance(double E_init, double mass_particle, TH1D* h_bethe ) { 
   for(int i=1; i <= h_bethe->GetNbinsX(); i++){
   h_bethe->SetBinContent( i, betheBloch(E_init, mass_particle));
   h_bethe->SetBinError(i, 0.001 );
   E_init = E_init - betheBloch(E_init, mass_particle);
   if(E_init <= 0) return;
   };
   };

   class LEN2E {
   public:
   void setmap(double); //input:E_ini [in MeV]
   double E(double); //input: length [in cm]
   LEN2E();
   ~LEN2E();

   private:
   double E_init; //E_ini

//map size [convert trklen to Edept]
int n_len=300;
double len_min=0;  
double len_max=300;

TH1* cumulative;
};

void LEN2E::setmap(double E0) {
E_init=E0;

//create the map to convert trklen to Edept
TH1D* dEdx;
dEdx = new TH1D("dEdx", "", n_len, len_min, len_max);
//hist_NIST(E_init, dEdx); //loading in dE/dx map based on E_int
hist_bethe_mean_distance(E_init, mass_particle, dEdx); //loading in dE/dx map based on E_int
cumulative = dEdx->GetCumulative();

delete dEdx;
}

double LEN2E::E(double len) {
double Edept_len=-9999;

int bin_cen=0;
//int bin_b=0;
//int bin_a=0;

//get Edept
//int n_len=cumulative->GetNbinsX();
bin_cen=cumulative->GetXaxis()->FindBin(len);
//bin_b=bin_cen-1;
//bin_a=bin_cen+1;
if (bin_cen>n_len) bin_cen=n_len;
if (bin_cen<0) bin_cen=0;

//bin_cen=cumulative->GetXaxis()->FindBin(len);
//bin_b=bin_cen-1;
//bin_a=bin_cen+1;
//if (bin_a>n_len) bin_a=n_len;
//if (bin_b<0) bin_b=0;

//double dept_cen=cumulative->GetBinContent(bin_cen);
double E_len=cumulative->GetBinContent(bin_cen);

//remove interpolation calc
//double dept_b=cumulative->GetBinContent(bin_b);
//double dept_a=cumulative->GetBinContent(bin_a);

//double m_dept=(dept_a-dept_b)/(cumulative->GetBinCenter(bin_a)-cumulative->GetBinCenter(bin_b));
//double b_dept=dept_a-m_dept*cumulative->GetBinCenter(bin_a);
//double E_len=b_dept+m_dept*len; 


if (E_init>0) { //E_init>0 
	if (len==0) Edept_len=0;
	if (len<0) Edept_len=-99999; //not-possible
	if (len>0) Edept_len=E_len;

	if (Edept_len>E_init) {
		std::cout<<"\nWARNING!! E_len>E_init! :: E_len="<<Edept_len<<" MeV;  E_init="<<E_init<<" MeV\n"<<endl;
		//Edept_len=E_init;
	}
} //E_init>0
else { //E_init<=0 [up-stream INT::KE_ff:=0]
	Edept_len=-99999;
} //E_init<=0 [up-stream INT::KE_ff:=0]

return Edept_len;
}

LEN2E::LEN2E(void) {}
LEN2E::~LEN2E(void) {}
*/
//Function to convert trklen to Edept -----------------------------------------------//




Double_t fitg(Double_t* x,Double_t *par) {
	double m=par[0];
	double s=par[1];
	double n=par[2];

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

Double_t agovg(Double_t* x,Double_t *par) {
	//g1
	double m1=par[0];
	double s1=par[1];
	double a1=par[2];
	double g1=-(x[0]-m1)*(x[0]-m1)/(2*s1*s1);

	//g2
	double m2=par[3];
	double s2=par[4];
	double a2=par[5];	
	double g2=-(x[0]-m2)*(x[0]-m2)/(2*s2*s2);

	//g2/g1
	double g_ov_g=0; 
	g_ov_g=(a1/a2)*TMath::Exp(g2-g1);
	if (m1==m2&&s1==s2&&a1==a2) g_ov_g=1;

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

//PID using stopping proton hypothesis
TFile f_pid("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/PDSPProd2/with_SCE_cali/dEdxrestemplates.root");
TProfile* dedx_range_pro = (TProfile*)f_pid.Get("dedx_range_pro");
double chi2pid(std::vector<double> &trkdedx, std::vector<double> &trkres) {
	int npt = 0;
	double chi2pro = 0;
	for (size_t i = 0; i<trkdedx.size(); ++i){ //hits
		//ignore the first and the last point
		if (i==0 || i==trkdedx.size()-1) continue;
		if (trkdedx[i]>1000) continue; //protect against large pulse height

		int bin = dedx_range_pro->FindBin(trkres[i]);
		//    std::cout<<"bin proton "<<bin<<std::endl; 
		if (bin>=1&&bin<=dedx_range_pro->GetNbinsX()){
			double bincpro = dedx_range_pro->GetBinContent(bin);
			if (bincpro<1e-6){//for 0 bin content, using neighboring bins
				bincpro = (dedx_range_pro->GetBinContent(bin-1)+dedx_range_pro->GetBinContent(bin+1))/2;
			}
			double binepro = dedx_range_pro->GetBinError(bin);
			if (binepro<1e-6){
				binepro = (dedx_range_pro->GetBinError(bin-1)+dedx_range_pro->GetBinError(bin+1))/2;
			}
			double errdedx = 0.04231+0.0001783*trkdedx[i]*trkdedx[i]; //resolution on dE/dx
			errdedx *= trkdedx[i];
			chi2pro += pow((trkdedx[i]-bincpro)/std::sqrt(pow(binepro,2)+pow(errdedx,2)),2);
			++npt;
		}
	}
	if (npt>0) return (chi2pro/npt);
	else return 9999;
}


//Gaussian fit
TF1* VFit(TH1D* h, Int_t col) {
	//pre-fit parameters
	float pre_mean=h->GetBinCenter(h->GetMaximumBin());
	float pre_max=h->GetBinContent(h->GetMaximumBin());
	float pre_rms=h->GetRMS();
	cout<<"mean: "<<pre_mean<<endl;
	cout<<"rms: "<<pre_rms<<endl;
	cout<<"max: "<<pre_max<<endl;
	cout<<""<<endl;

	//1st fitting ---------------------------------------------------------//
	TF1 *gg=new TF1("gg", fitg, pre_mean-3*pre_rms, pre_mean+3*pre_rms, 3);
	gg->SetParameter(0,pre_mean);
	gg->SetParameter(1,pre_rms);
	gg->SetParameter(2,pre_max);
	//if (pre_rms>1.0e+06) { gg->SetParLimits(1,0,100); }

	//gg->SetLineColor(col);
	//gg->SetLineStyle(2);
	h->Fit("gg","remn");

	//2nd fitting -----------------------------------------------------------------------------------------------------------//
	TF1 *g=new TF1("g", fitg, gg->GetParameter(0)-3.*gg->GetParameter(1), gg->GetParameter(0)+3.*gg->GetParameter(1), 3);
	//TF1 *g=new TF1("g",fitg,0.3,0.5,3);

	//TF1 *g=new TF1("g",fitg,gg->GetParameter(0)-1,gg->GetParameter(0)+.5,3);
	g->SetParameter(0, gg->GetParameter(0));
	g->SetParameter(1, gg->GetParameter(1));
	g->SetParameter(2, gg->GetParameter(2));

	//g->SetParLimits(0,gg->GetParameter(0)-1*gg->GetParameter(1), gg->GetParameter(0)+1*gg->GetParameter(1));
	//double sss=gg->GetParameter(1); if (sss<0) sss=-sss;
	//g->SetParLimits(1,0,5.*sss);
	//g->SetParLimits(2,0,100.*sqrt(pre_max));

	g->SetLineColor(col);
	g->SetLineStyle(2);
	g->SetLineWidth(2);

	h->Fit("g","remn");
	return g;
}

TF1* VNFit(TH1D* h, float pre_mean, float n_sigma) {
	//pre-fit parameters
	//float pre_mean=h->GetBinCenter(h->GetMaximumBin());
	float pre_max=h->GetMaximum();
	float pre_rms=h->GetRMS();
	cout<<"mean: "<<pre_mean<<endl;
	cout<<"rms: "<<pre_rms<<endl;
	cout<<"max: "<<pre_max<<endl;
	cout<<""<<endl;

	//1st fitting
	TF1 *gg=new TF1("gg",fitg,pre_mean-n_sigma*pre_rms,pre_mean+n_sigma*pre_rms,3);
	gg->SetParameter(0,pre_mean);
	gg->SetParameter(1,pre_rms);
	gg->SetParameter(2,pre_max);
	//if (pre_rms>1.0e+06) { gg->SetParLimits(1,0,100); }

	//gg->SetLineColor(col);
	gg->SetLineStyle(2);
	h->Fit("gg","remn");

	//2nd fitting
	TF1 *g=new TF1("g",fitg,gg->GetParameter(0)-n_sigma*gg->GetParameter(1),gg->GetParameter(0)+n_sigma*gg->GetParameter(1),3);
	//TF1 *g=new TF1("g",fitg,gg->GetParameter(0)-1,gg->GetParameter(0)+.5,3);
	g->SetParameter(0,gg->GetParameter(0));
	g->SetParameter(1,gg->GetParameter(1));
	g->SetParameter(2,gg->GetParameter(2));

	//g->SetParLimits(0,gg->GetParameter(0)-3*gg->GetParameter(1), gg->GetParameter(0)+3*gg->GetParameter(1));
	//double sss=gg->GetParameter(1); if (sss<0) sss=-sss;
	//g->SetParLimits(1,0,5.*sss);
	//g->SetParLimits(2,0,10.*sqrt(pre_max));

	//g->SetLineColor(col);
	g->SetLineStyle(2);
	g->SetLineWidth(2);

	h->Fit("g","remn");
	return g;
}

//beam momentum reweighting --------------------------------------------------------------------------------------------------------------------------------------------//
TString fpath_bmrw=Form("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/rw/");
TString file_bmrw=Form("bmrw_bestfit.root");
TFile *f_bmrw=new TFile(Form("%s%s",fpath_bmrw.Data(),file_bmrw.Data()));
TF1 *bmrw_func=(TF1 *)f_bmrw->Get("bmrw_minchi2");

//efficiency calc.-------------------------------------------------------------------//
double err_p(double nom, double denom) {
	double ef=nom/denom;

	double err_p=0;
	double err_m=0;
	if (ef>0.&&ef<1.) {
		err_p=(1.0/(sqrt(denom)))*sqrt((ef*(1.0-ef)));
		err_m=err_p;
	}
	else {
		double inverse_denom=1./denom;
		double err=1-pow(0.72,inverse_denom);

		if (ef==0.) {
			err_p=err;
			err_m=0.;
		}
		if (ef==1.) {
			err_p=0.;
			err_m=err;
		}
	}
	return err_p;
}

double err_m(double nom, double denom) {
	double ef=nom/denom;

	double err_p=0;
	double err_m=0;
	if (ef>0.&&ef<1.) {
		err_p=(1.0/(sqrt(denom)))*sqrt((ef*(1.0-ef)));
		err_m=err_p;
	}
	else {
		double inverse_denom=1./denom;
		double err=1-pow(0.72,inverse_denom);

		if (ef==0.) {
			err_p=err;
			err_m=0.;
		}
		if (ef==1.) {
			err_p=0.;
			err_m=err;
		}
	}
	return err_m;
}


double Convert_Proton_KE_Spectrometer_to_KE_ff(double KE_RecoBeam, TString key, int syst){
	double out = KE_RecoBeam;
	double delta_E = 0.;
	double p0 = 0., p1 = 0., p2 = 0.;

	if(key == "ElasTrue"){
		if(syst == 0){
			p0 = 39.78;
			p1 = -0.2396;
			p2 = 0.0004498;
		}
		else if(syst == -1){
			p0 = 17.92;
			p1 = -0.1334;
			p2 = 0.0003075;
		}
		else if(syst == 1){
			p0 = 60.52;
			p1 = -0.3404;
			p2 = 0.0005857;
		}
		else{
			return out;
		}
	}
	else if(key == "AllTrue"){
		if(syst == 0){
			p0 = 51.33;
			p1 = -0.2954;
			p2 = 0.0005164;
		}
		else if(syst == -1){
			p0 = 28.91;
			p1 = -0.1827;
			p2 = 0.0003602;
		}
		else if(syst == 1){
			p0 = 72.64;
			p1 = -0.4028;
			p2 = 0.0006661;
		}
		else{
			return out;
		}
	}
	else if(key == "ElasFitted"){
		if(syst == 0){
			p0 = 27.68;
			p1 = -0.1727;
			p2 = 0.0003779;
		}
		else if(syst == -1){
			p0 = -6.135;
			p1 = -0.01069;
			p2 = 0.0001675;
		}
		else if(syst == 1){
			p0 = 60.19;
			p1 = -0.3284;
			p2 = 0.0005808;
		}
		else{
			return out;
		}
	}
	else if(key == "AllFitted"){
		if(syst == 0){
			p0 = 37.57;
			p1 = -0.2144;
			p2 = 0.0004282;
		}
		else if(syst == -1){
			p0 = 12.66;
			p1 = -0.09451;
			p2 = 0.0002699;
		}
		else if(syst == 1){
			p0 = 61.69;
			p1 = -0.3305;
			p2 = 0.0005820;
		}
		else{
			return out;
		}
	}
	else if(key == "Data"){
		if(syst == 0){
			p0 = 47.35;
			p1 = -0.1748;
			p2 = 0.0003067;
		}
		else if(syst == -1){
			p0 = 28.27;
			p1 = -0.07865;
			p2 = 0.0001744;
		}
		else if(syst == 1){
			p0 = 65.99;
			p1 = -0.2688;
			p2 = 0.0004365;
		}
		else{
			return out;
		}
	}
	else {
		return out;
	}

	delta_E = p0 + p1 * out + p2 * out * out;

	return out - delta_E;

}

//chi^2 functions ------------------------------------------------------------------------------------------------------//
//neyman_chi2
double neyman_chi2_data_mc(vector<double> D /*data*/, vector<double> er_D, vector<double> E, vector<double> er_E) {
	double chi2=0.;
	for (int i=0; i<(int)D.size(); i++) {
		double nom=pow(D.at(i)-E.at(i),2);
		double denom=pow(er_D.at(i),2);	

		if (denom) chi2+=nom/denom;	
	} 
	return chi2;
}

//pearson_chi2
double pearson_chi2_data_mc(vector<double> D /*data*/, vector<double> er_D, vector<double> E, vector<double> er_E) {
	double chi2=0.;
	for (int i=0; i<(int)D.size(); i++) {
		double nom=pow(D.at(i)-E.at(i),2);
		double denom=pow(er_E.at(i),2);	

		if (denom) chi2+=nom/denom;	
	} 
	return chi2;
}

double ml_data_mc(vector<double> D /*data*/, vector<double> er_D, vector<double> E /*MC*/, vector<double> er_E) {
	double ml=0.;
	for (int i=0; i<(int)D.size(); i++) { //sum over all bins
		double mui=E.at(i);
		double ni=D.at(i);

		double first=mui-ni;
		double second=0.;	

		if (ni>0&&mui>0) { 	
			second=ni*TMath::Log(ni/mui);
		}
		else { 	
			second=0;
		}
		ml+=first+second;
	} //sum over all bins

	return 2.*ml;
}

//Covariance matrix =========================================================================================//
//Matrix element calculation using Combined Neyman-Pearson(CNP)
double cov_cnp(double m /*data(measurement)*/, double p /*prediction from mc*/) {
	double cnp=0;
	if(m!=0&&p!=0) { 
		cnp=3.*m*p/(p+2.*m);
	}
	if (m==0&&p!=0) { //measurement bin=0
		cnp=p/2.; //Poisson approx.
	}
	if (p==0&&m!=0) { //prediction from mc=0
		cnp=m; //Neyman approx.
	}
	//if m=p=0, no info, set to zero

	return cnp;
}


