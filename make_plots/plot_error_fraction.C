#include <iostream>
#include <vector>
#include <utility>  // for std::pair
#include <string>  // for std::string
#include "THStack.h"

void plot_error_fraction() {

	//Outputfig -------------------------------------------------------------------------//
	TString out_fig=Form("./plots_sys_percentage/");

	//load xs histograms ----------------------------------------------------------------------------//
	//statistical uncertainty only
	//TString str_cen="./xs_files/xs_Eslice_dE20MeV_40slcs_bmrw_data_cen.root"; //stat. only file
	TString str_cen="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_data_cen_newslcid.root"; //stat. only file
	TFile *f_cen = TFile::Open(str_cen.Data()); //xs_central
	TGraphErrors* xs_cen=(TGraphErrors* )f_cen->Get("gr_recoxs"); //xs_central

	//sys1(xs shift due to keff shift)
	//TString str_bmrw="./xs_files/xs_Eslice_dE20MeV_40slcs_bmrw_data_cen_sysbmrw.root"; //sys only
	TString str_bmrw="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_SYS_BMRW_newslcid.root"; //sys only
	TFile *f_bmrw = TFile::Open(str_bmrw.Data()); //bmrw only
	TGraphAsymmErrors* xs_sys_bmrw=(TGraphAsymmErrors* )f_bmrw->Get("reco_xs_sysonly"); //sys_only

	//sys2(xs shiftt due to KE theory)
	//TString str_bb="./xs_files/xs_Eslice_dE20MeV_40slcs_bmrw_data_cen_sysBB.root"; //sys only
	TString str_bb="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_SYS_BB_newslcid.root"; //sys only
	TFile *f_bb = TFile::Open(str_bb.Data()); //bb only
	TGraphAsymmErrors* xs_sys_bb=(TGraphAsymmErrors* )f_bb->Get("reco_xs_sysonly"); //sys_only

	//sys3(El-bks:simple scaling from bkg-rich channel)
	//TString str_el="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_SYS_ElScale_newslcid.root"; //sys only
	//TFile *f_el = TFile::Open(str_el.Data()); //el-scale only
	//TGraphAsymmErrors* xs_sys_el=(TGraphAsymmErrors* )f_el->Get("reco_xs_sysonly"); //sys_only
	
	//sys3(xs shift due to Bayesian unfolding process) 


	//Fill error composition histogram --------------------------------------------------//
	//xs central values(stat.)
	Double_t *ke=xs_cen->GetX();
	Double_t *err_ke=xs_cen->GetEX();
	Double_t *cen_xs=xs_cen->GetY();
	Double_t *err_cen_xs=xs_cen->GetEY();

	//sys1(xs shift due to keff shift)
	Double_t *err_h_ke=xs_sys_bmrw->GetEXhigh();
	Double_t *err_l_ke=xs_sys_bmrw->GetEXlow();
	Double_t *err_h_sys_bmrw=xs_sys_bmrw->GetEYhigh();
	Double_t *err_l_sys_bmrw=xs_sys_bmrw->GetEYlow();

	//sys2(xs shiftt due to KE theory)
	Double_t *err_h_sys_bb=xs_sys_bb->GetEYhigh();
	Double_t *err_l_sys_bb=xs_sys_bb->GetEYlow();

	//sys3(xs shiftt due to el-scaling)
	//Double_t *err_h_sys_el=xs_sys_el->GetEYhigh();
	//Double_t *err_l_sys_el=xs_sys_el->GetEYlow();

	int n=xs_cen->GetN();
	int n_bin=err_ke[0]; //+-10 MeV, bin_size=20 MeV
	double ke_min=ke[n-1]-(double)n_bin;
	double ke_max=ke[0]+(double)n_bin;

	TH1D* erry_CEN=new TH1D("erry_CEN","", n, ke_min, ke_max); erry_CEN->Sumw2();
	TH1D* erry_SYS_BMRW=new TH1D("erry_SYS_BMRW","", n, ke_min, ke_max); erry_SYS_BMRW->Sumw2();
	TH1D* erry_SYS_BB=new TH1D("erry_SYS_BB","", n, ke_min, ke_max); erry_SYS_BB->Sumw2();
	//TH1D* erry_SYS_EL=new TH1D("erry_SYS_EL","", n, ke_min, ke_max); erry_SYS_EL->Sumw2();

	//construct individual error histogram -------------------------------------------//
	for (int i=n-1; i>=0; --i) {
		int j=n-i;
		cout<<"i="<<i<<" j="<<j<<endl;

		//std::cout<<"["<<i<<"] "<<"ke:"<<ke[i]<<" +- "<<err_ke[i]<<std::endl;
		//systematic u
		double erry_h_cen=err_cen_xs[i];
		double erry_l_cen=erry_h_cen;
		double erry_cen=pow(erry_h_cen,2)+pow(erry_l_cen,2);

		//sys: bmrw
		double erry_h_sys_bmrw=err_h_sys_bmrw[i];
		double erry_l_sys_bmrw=err_l_sys_bmrw[i];
		double erry_bmrw=pow(erry_h_sys_bmrw,2)+pow(erry_l_sys_bmrw,2);

		//sys: bb
		double erry_h_sys_bb=err_h_sys_bb[i];
		double erry_l_sys_bb=err_l_sys_bb[i];
		double erry_bb=pow(erry_h_sys_bb,2)+pow(erry_l_sys_bb,2);

		//sys: el
		//double erry_h_sys_el=err_h_sys_el[i];
		//double erry_l_sys_el=err_l_sys_el[i];
		//double erry_el=pow(erry_h_sys_el,2)+pow(erry_l_sys_el,2);

		//total
		double erry_all=erry_cen+erry_bmrw+erry_bb;
		//double erry_all=erry_cen+erry_bmrw+erry_bb+erry_el;
		double frac_cen=100.*erry_cen/erry_all;
		double frac_erry_bmrw=100.*erry_bmrw/erry_all;
		double frac_erry_bb=100.*erry_bb/erry_all;
		double frac_erry_el=100.*erry_el/erry_all;

		if (erry_all>0) {
			erry_CEN->SetBinContent(j, frac_cen);
			erry_SYS_BMRW->SetBinContent(j, frac_erry_bmrw);
			erry_SYS_BB->SetBinContent(j, frac_erry_bb);
			//erry_SYS_EL->SetBinContent(j, frac_erry_el);

			std::cout<<frac_cen<<" | "<<frac_erry_bmrw<<" | "<<frac_erry_bb<<std::endl;
		}	
	}

        //plot style --------------------------------------------------------------------//
        gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
        gROOT->SetStyle("protoDUNEStyle");
        gROOT->ForceStyle();
        gStyle->SetTitleX(0.5);
        gStyle->SetTitleAlign(23);
        gStyle->SetOptStat(0);

	TCanvas *c_=new TCanvas(Form("c_"),"",900, 600);
	c_->Divide(1,1);
	erry_CEN->SetFillColor(kAzure+1);
	erry_SYS_BMRW->SetFillColor(3);
	erry_SYS_BB->SetFillColor(6);
	erry_SYS_EL->SetFillColor(kOrange-3);


	THStack* hs=new THStack("hs","");
	hs->Add(erry_CEN);
	hs->Add(erry_SYS_BMRW);
	hs->Add(erry_SYS_BB);
	hs->Add(erry_SYS_EL);

	TH2D *f2d=new TH2D("f2d","",440,0,340,130,0,130);
	f2d->SetTitle("; Proton Kinetic Energy [MeV]; Systematic Percentage [%]");
	f2d->Draw();
	hs->Draw("same");

	//TLegend
	TLegend *leg = new TLegend(0.12,0.75,0.86,0.9);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->AddEntry(erry_CEN, "Stat.","f");
	leg->AddEntry(erry_SYS_BMRW, "Sys. KE_{ff} ","f");
	leg->AddEntry(erry_SYS_BB, "Sys. KE (Bethe-Bloch)","f");
        //leg->AddEntry((TObject*)0, "", "");
        //leg->AddEntry((TObject*)0, "", "");
	//leg->AddEntry(erry_SYS_EL, "Sys. Elastic-Scatt. Bkg.","f");
	leg->SetNColumns(3);
	leg->Draw();

	//Preliminary
        TLatex **txt_preliminary=new TLatex*[1];
        txt_preliminary[0]=new TLatex(150, 40, Form("Preliminary"));
        txt_preliminary[0]->SetTextColor(0);
        txt_preliminary[0]->SetTextSize(0.09);
        txt_preliminary[0]->Draw();

	c_->Print(Form("%ssys_percentage.eps",out_fig.Data()));


	/*
	//[1]fill error in one histogram
	//[2]sum all histogram in quardidure and take square root
	//[3]normalize each histogram (each/sum)
	//plot results

	std::vector<std::pair<std::string, int>> vec;
	// Add some elements to the vector
	vec.push_back(std::make_pair("apple", 3));
	vec.push_back(std::make_pair("banana", 5));
	vec.push_back(std::make_pair("orange", 4));

	// Print the elements in the vector
	for (const auto& p : vec) {
	std::cout << p.first << ": " << p.second << std::endl;
	}

	//return 0;
	*/

}
