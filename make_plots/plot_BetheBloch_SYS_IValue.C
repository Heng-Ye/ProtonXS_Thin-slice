//#include "./cali/dedx_function_35ms.h"
#include "../headers/BasicParameters.h"
#include "../headers/BasicFunctions.h"
//#include "../headers/ESliceParams.h"
//#include "../headers/util.h"
//#include "../headers/ESlice.h"
#include "../headers/BetheBloch.h"
#include "../headers/BetheBloch_SYS_I.h"

using namespace std;
using namespace ROOT::Math;

void plot_BetheBloch_SYS_IValue() {
	//Output plot ------------------------------
	TString fout=Form("./plot_BetheBloch_IV/");

        //plot style --------------------------------------------------------------------//
        gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
        gROOT->SetStyle("protoDUNEStyle");
        gROOT->ForceStyle();
        gStyle->SetTitleX(0.5);
        gStyle->SetTitleAlign(23);
        gStyle->SetOptStat(0);

	//Basic configure ------------------//
	BetheBloch BB;
	BB.SetPdgCode(2212);
	BetheBloch_SYS_I BB_SYS_I; 
	BB_SYS_I.SetPdgCode(2212);

	double mu_denom_data=411.05145837595467; //new with event-by-event R corr (const E-loss using stopping protons)
	double sg_denom_data=47.48714821962207; //new with event-by-event R corr (const E-loss using stopping protons)

	double mu_nom_mc=416.224743039812; //for mc(KEbeam-const) with R=1 (R=Ratio of KEff(Fit)/(KEbeam-dE))
	double sg_nom_mc=42.786018962508784; //


	//[0]Toy model to generate KEff based on expected shape 
	//[1]Given KEff
	double keff_reco=0;
        keff_reco=mu_denom_data;
	//[2]Expected Length (stopping protons) based on stopping protons
	double range_reco=0;
	double range_sys_reco=0;
	range_reco=BB.RangeFromKESpline(keff_reco);
	range_sys_reco=BB_SYS_I.RangeFromKESpline(keff_reco);
	cout<<"keff_reco="<<keff_reco<<" range_reco="<<range_reco<<endl;	
	cout<<"keff_reco="<<keff_reco<<" range_reco="<<range_sys_reco<<endl;	

	//[3]KE vs Length -------------------------------------------------------------------//
	int n_seg=100;
	double dl_BB=range_reco/(double)n_seg;

	vector<double> LEN;
	vector<double> KE_bb;
	vector<double> KE_bb_sys_I;
	vector<double> KE_RATIO;
	double len=0;
	int cnt=0;
	for (int i=0; i<n_seg-1; ++i) {
		if (len>=range_reco) { 
			len=range_reco;
			cnt++;
			if (cnt>1) break;
		}

		LEN.push_back(len);
		
		double kebb=BB.KEAtLength(keff_reco,len); 
		double kebb_sys_i=BB_SYS_I.KEAtLength(keff_reco,len); 
		if (kebb<0||kebb_sys_i<0) break;

		KE_bb.push_back(kebb);
		KE_bb_sys_I.push_back(kebb_sys_i);
		double tmp_ratio=(kebb_sys_i-kebb)/kebb;
		KE_RATIO.push_back(100.*tmp_ratio);
		cout<<"len="<<len<<" kebb="<<kebb<<" kebb_sys_i="<<kebb_sys_i<<" tmp_ratio="<<tmp_ratio<<endl;

		len+=(double)i*dl_BB;

	}

	std::cout<<"KE_end="<<KE_bb.at(KE_bb.size()-1)<<std::endl;
	std::cout<<"KE_end2="<<KE_bb_sys_I.at(KE_bb_sys_I.size()-1)<<std::endl;
	std::cout<<"KE_RATIO="<<KE_RATIO.at(KE_RATIO.size()-1)<<std::endl;


	TGraph *gr_len_kebb = new TGraph(LEN.size(), &LEN.at(0), &KE_bb.at(0));
	TGraph *gr_len_kebb_sys_I = new TGraph(LEN.size(), &LEN.at(0), &KE_bb_sys_I.at(0));
	TGraph *gr_len_ratio = new TGraph(LEN.size(), &LEN.at(0), &KE_RATIO.at(0));
	
	
	//[4]Plot ------------------------------------------------------------------------//
	TCanvas *c_=new TCanvas("c_","");
	c_->Divide(1,2);
	TH2D* f2d=new TH2D("f2d","", 100, 0, 100, 500, -10, 500);
	f2d->SetTitle("Bethe-Bloch KE; Track Length[cm]; Proton Kinetic Energy[MeV]"); 
	c_->cd(1);
	f2d->Draw();
	gr_len_kebb->Draw("pl same");
	gr_len_kebb_sys_I->SetMarkerColor(2);
	gr_len_kebb_sys_I->SetLineColor(2);
	gr_len_kebb_sys_I->Draw("pl same");

	TLegend *leg = new TLegend(0.5,0.6,0.85,0.9);
	leg->SetFillStyle(0);
	leg->AddEntry(gr_len_kebb, "I=189 eV [theory]","lp");
	leg->AddEntry(gr_len_kebb_sys_I, "I=205 eV [sys.]","lp");
	leg->Draw();

	c_->cd(2);
	gr_len_ratio->SetTitle("; Track Length[cm]; KE Ratio [a.u.]");
	TH2D* f2d2=new TH2D("f2d2","", 100, 0, 100, 40, -1, 10);
	f2d2->SetTitle("; Track Length[cm]; 100*(1-#DeltaKE/KE) [%]"); 
	f2d2->Draw();
	gr_len_ratio->Draw("pl same");
	c_->Print(Form("%ske_ke.eps",fout.Data()));







}
