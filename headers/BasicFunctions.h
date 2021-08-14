#include "TMath.h"

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

       h->Fit("g","rem+");
       return g;
}


//beam momentum reweighting --------------------------------------------------------------------------------------------------------------------------------------------//
TString fpath_bmrw=Form("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/rw/");
TString file_bmrw=Form("bmrw_bestfit.root");
TFile *f_bmrw=new TFile(Form("%s%s",fpath_bmrw.Data(),file_bmrw.Data()));
TF1 *bmrw_func=(TF1 *)f_bmrw->Get("bmrw_minchi2");


