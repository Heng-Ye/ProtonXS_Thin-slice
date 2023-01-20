#include <iostream>
#include <vector>
#include <utility>  // for std::pair
#include <string>  // for std::string
#include "THStack.h"

void plot_dataXS_withSys() {

        //plot style --------------------------------------------------------------------//
        gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
        gROOT->SetStyle("protoDUNEStyle");
        gROOT->ForceStyle();
        gStyle->SetTitleX(0.5);
        gStyle->SetTitleAlign(23);
        gStyle->SetOptStat(0);

        gStyle->SetPadRightMargin(0.0);
 	gStyle->SetPadLeftMargin(0.16);
    	gStyle->SetPadRightMargin(0.04);//0.02

	//Outputfig -------------------------------------------------------------------------//
	TString out_fig=Form("./plots_sys_xs/");

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
	TString str_el="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_SYS_ElScale_newslcid.root"; //sys only
	TFile *f_el = TFile::Open(str_el.Data()); //el-scale only
	TGraphAsymmErrors* xs_sys_el=(TGraphAsymmErrors* )f_el->Get("reco_xs_sysonly"); //sys_only

	//sys4(MisIDP-bks:simple scaling from bkg-rich channel)
	TString str_misidp="./xs_files_newslcid/xs_Eslice_dE20MeV_40slcs_SYS_MisIDPScale_newslcid.root"; //sys only
	TFile *f_misidp = TFile::Open(str_misidp.Data()); //misidp-scale only
	TGraphAsymmErrors* xs_sys_misidp=(TGraphAsymmErrors* )f_misidp->Get("reco_xs_sysonly"); //sys_only

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
	Double_t *err_h_sys_el=xs_sys_el->GetEYhigh();
	Double_t *err_l_sys_el=xs_sys_el->GetEYlow();

	//sys4(xs shiftt due to misidp-scaling)
	Double_t *err_h_sys_misidp=xs_sys_misidp->GetEYhigh();
	Double_t *err_l_sys_misidp=xs_sys_misidp->GetEYlow();

	vector<double> KE_CEN;
	vector<double> err_KE_CEN;

	vector<double> KE_SYS;
	vector<double> err_h_KE_SYS;
	vector<double> err_l_KE_SYS;

	vector<double> XS_CEN;
	vector<double> err_XS_CEN;

	//vector<double> XS_SYS;
	vector<double> err_h_XS_SYS;
	vector<double> err_l_XS_SYS;


	int n=xs_cen->GetN();
	int n_bin=err_ke[0]; //+-10 MeV, bin_size=20 MeV
	double ke_min=ke[n-1]-(double)n_bin;
	double ke_max=ke[0]+(double)n_bin;

	//construct individual error histogram -------------------------------------------//
	for (int i=n-1; i>=0; --i) {
		int j=n-i;
		//cout<<"i="<<i<<" j="<<j<<endl;

		//ke
		double errx_cen=err_ke[i];
		double errx_h_ke=err_h_ke[i];
		double errx_l_ke=err_l_ke[i];

		//stat.
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
		double erry_h_sys_el=err_h_sys_el[i];
		double erry_l_sys_el=err_l_sys_el[i];
		double erry_el=pow(erry_h_sys_el,2)+pow(erry_l_sys_el,2);

		//sys: misid:p
		double erry_h_sys_misidp=err_h_sys_misidp[i];
		double erry_l_sys_misidp=err_l_sys_misidp[i];
		double erry_misidp=pow(erry_h_sys_misidp,2)+pow(erry_l_sys_misidp,2);

		//total
		//double erry_all=erry_cen+erry_bmrw+erry_bb;
		//double erry_all=erry_cen+erry_bmrw+erry_bb+erry_el;
		double erry_all=erry_cen+erry_bmrw+erry_bb+erry_el+erry_misidp;

		double frac_cen=100.*erry_cen/erry_all;
		double frac_erry_bmrw=100.*erry_bmrw/erry_all;
		double frac_erry_bb=100.*erry_bb/erry_all;
		double frac_erry_el=100.*erry_el/erry_all;
		double frac_erry_misidp=100.*erry_misidp/erry_all;


		//double errx_h_ke_tot=sqrt(pow(errx_cen,2)+pow(errx_h_ke,2));
		//double errx_l_ke_tot=sqrt(pow(errx_cen,2)+pow(errx_l_ke,2));
		//double erry_h_xs_tot=sqrt(pow(erry_h_cen,2)+pow(erry_h_sys_bmrw,2)+pow(erry_h_sys_bb,2));
		//double erry_l_xs_tot=sqrt(pow(erry_l_cen,2)+pow(erry_l_sys_bmrw,2)+pow(erry_l_sys_bb,2));

		//only sys errors
		double errx_h_ke_tot=sqrt(0+pow(errx_h_ke,2));
		double errx_l_ke_tot=sqrt(0+pow(errx_l_ke,2));
		//double erry_h_xs_tot=sqrt(0+pow(erry_h_sys_bmrw,2)+pow(erry_h_sys_bb,2));
		//double erry_l_xs_tot=sqrt(0+pow(erry_l_sys_bmrw,2)+pow(erry_l_sys_bb,2));
		//double erry_h_xs_tot=sqrt(0+pow(erry_h_sys_bmrw,2)+pow(erry_h_sys_bb,2)+pow(erry_h_sys_el,2));
		//double erry_l_xs_tot=sqrt(0+pow(erry_l_sys_bmrw,2)+pow(erry_l_sys_bb,2)+pow(erry_l_sys_el,2));
		double erry_h_xs_tot=sqrt(0+pow(erry_h_sys_bmrw,2)+pow(erry_h_sys_bb,2)+pow(erry_h_sys_el,2)+pow(erry_h_sys_misidp,2));
		double erry_l_xs_tot=sqrt(0+pow(erry_l_sys_bmrw,2)+pow(erry_l_sys_bb,2)+pow(erry_l_sys_el,2)+pow(erry_l_sys_misidp,2));


		if (erry_all>0) {
			//std::cout<<ke[i]-err_ke[i]<<"-"<<ke[i]+err_ke[i]<<" & "<<cen_xs[i]<<" & \\pm"<<erry_h_cen<<" & $^{+"<<erry_h_sys_bmrw<<"}_{-"<<erry_l_sys_bmrw<<"}$ & $^{+"<<erry_h_sys_bb<<"}_{-"<<erry_l_sys_bb<<"}$ & $^{+"<<erry_h_sys_el<<"}_{-"<<erry_l_sys_el<<"}$ & $^{+"<<erry_h_sys_misidp<<"}_{-"<<erry_l_sys_misidp<<"}$ \\\\"<<endl;
			
			printf("%.0f-%.0f & %.1f & $\\pm$%.1f(%.1f\\,\\\%) & $^{+%.1f(%.1f\\,\\\%)}_{-%.1f(%.1f\\,\\\%)}$ & $^{+%.1f(%.1f\\,\\\%)}_{-%.1f(%.1f\\,\\\%)}$ & $^{+%.1f(%.1f\\,\\\%)}_{-%.1f(%.1f\\,\\\%)}$ & $^{+%.1f(%.1f\\,\\\%)}_{-%.1f(%.1f\\,\\\%)}$ \\\\ \n",ke[i]-err_ke[i],ke[i]+err_ke[i],cen_xs[i],erry_h_cen,100.*erry_h_cen/cen_xs[i],erry_h_sys_bmrw,100.*erry_h_sys_bmrw/cen_xs[i], erry_l_sys_bmrw, 100.*erry_l_sys_bmrw/cen_xs[i],erry_h_sys_bb,100.*erry_h_sys_bb/cen_xs[i],erry_l_sys_bb,100.*erry_l_sys_bb/cen_xs[i],erry_h_sys_el,100.*erry_h_sys_el/cen_xs[i],erry_l_sys_el,100.*erry_l_sys_el/cen_xs[i],erry_h_sys_misidp,100.*erry_h_sys_misidp/cen_xs[i],erry_l_sys_misidp,100.*erry_l_sys_misidp/cen_xs[i]);

			//erry_CEN->SetBinContent(j, frac_cen);
			//erry_SYS_BMRW->SetBinContent(j, frac_erry_bmrw);
			//erry_SYS_BB->SetBinContent(j, frac_erry_bb);
			//std::cout<<frac_cen<<" | "<<frac_erry_bmrw<<" | "<<frac_erry_bb<<std::endl;

			KE_CEN.push_back(ke[i]);
			err_KE_CEN.push_back(err_ke[i]);

			err_h_KE_SYS.push_back(errx_h_ke_tot);
			err_l_KE_SYS.push_back(errx_l_ke_tot);

			XS_CEN.push_back(cen_xs[i]);
			err_XS_CEN.push_back(err_cen_xs[i]);

			err_h_XS_SYS.push_back(erry_h_xs_tot);
			err_l_XS_SYS.push_back(erry_l_xs_tot);

			//cout<<"KE_CEN:"<<ke[i]<<" "
		}	
	}
	TGraphErrors* xs_reco_stat=new TGraphErrors(KE_CEN.size(), &KE_CEN.at(0), &XS_CEN.at(0), &err_KE_CEN.at(0), &err_XS_CEN.at(0));
	TGraphAsymmErrors* xs_reco_sys=new TGraphAsymmErrors(KE_CEN.size(), &KE_CEN.at(0), &XS_CEN.at(0), &err_l_KE_SYS.at(0), &err_h_KE_SYS.at(0), &err_l_XS_SYS.at(0), &err_h_XS_SYS.at(0));

	//Load other xs models
        //Geant4 -------------------------------------------------------------------------------------//
        TFile f_xs("/dune/data2/users/hyliao/GeantReweight/xs_cascade/proton_cross_section.root");
        TGraph *total_inel_KE = (TGraph*)f_xs.Get("inel_KE");

	//Neutrino Gen -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
        TFile f_other("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/make_plots/model_prediction/model_pred.root");
        TGraph *xs_GENIE_ha2018 = (TGraph*)f_other.Get("xs_GENIE_ha2018");
        TGraph *xs_NEUT_2019 = (TGraph*)f_other.Get("xs_NEUT_2019");
        TGraph *xs_NuWRO_2019 = (TGraph*)f_other.Get("xs_NuWRO_2019");
        TGraph *xs_GENIE_hN2018 = (TGraph*)f_other.Get("xs_GENIE_hN2018");
        TGraph *xs_GENIE_INCL_pp = (TGraph*)f_other.Get("xs_GENIE_INCL_pp");
	


        TCanvas *c_xs = new TCanvas("c_xs", "c_xs", 1400, 900);
	c_xs->Divide(1,1);
	c_xs->cd(1);

        float ymax=1600;
        //float xmax=400;
        float xmax=200;
        float xmin=0;

        TH2D *f2d_xs=new TH2D("f2d_xs","",450,xmin,xmax,ymax,0,ymax);
        f2d_xs->GetXaxis()->SetTitle("Proton Kinetic Energy [MeV]");
        f2d_xs->GetYaxis()->SetTitleOffset(1.3);
        f2d_xs->GetYaxis()->SetTitle("P-Ar inelastic cross section [mb]");
        f2d_xs->Draw("");

        //gr_truexs->SetLineWidth(2);
        //gr_truexs->SetMarkerColor(3);
        //gr_truexs->SetLineColor(3);
        //gr_truexs->SetMarkerStyle(21);
        //gr_truexs->SetMarkerSize(1.5);


	//reco_xs_all->SetMarkerColor(2);
	//reco_xs_all->SetLineColor(2);
	//reco_xs_all->SetMarkerStyle(20);
	//reco_xs_all->Draw("p same");
	//
	xs_reco_stat->SetMarkerStyle(20);
	xs_reco_sys->SetMarkerStyle(20);
	xs_reco_stat->SetLineColor(1);
	xs_reco_sys->SetLineColor(1);
	xs_reco_sys->SetMarkerColor(1);

	//xs_GENIE_ha2018->SetMarkerColor(1); xs_GENIE_ha2018->SetLineColor(1); 
	xs_GENIE_INCL_pp->SetMarkerColor(6); xs_GENIE_INCL_pp->SetLineColor(6);
	//xs_GENIE_INCL_pp->SetMarkerColor(6); xs_GENIE_INCL_pp->SetLineColor(6);
	//xs_GENIE_hN2018->SetMarkerColor(8); xs_GENIE_hN2018->SetLineColor(8);
	//xs_NEUT_2019->SetMarkerColor(7); xs_NEUT_2019->SetLineColor(7);
	xs_NEUT_2019->SetMarkerColor(3); xs_NEUT_2019->SetLineColor(3);
	xs_GENIE_hN2018->SetMarkerColor(kOrange-3); xs_GENIE_hN2018->SetLineColor(kOrange-3);
	//xs_NEUT_2019->SetMarkerColor(kOrange-3); xs_NEUT_2019->SetLineColor(kOrange-3);
	//xs_GENIE_hN2018->SetMarkerColor(3); xs_GENIE_hN2018->SetLineColor(3);

	//xs_NEUT_2019->SetMarkerColor(kOrange-3); xs_NEUT_2019->SetLineColor(kOrange-3);
	//xs_GENIE_hN2018->SetMarkerColor(7); xs_GENIE_hN2018->SetLineColor(7);

	xs_NuWRO_2019->SetMarkerColor(4); xs_NuWRO_2019->SetLineColor(4);	

	xs_GENIE_ha2018->SetLineStyle(8);
	xs_NEUT_2019->SetLineStyle(10);
	xs_GENIE_INCL_pp->SetLineStyle(2);
	xs_GENIE_hN2018->SetLineStyle(2);
	xs_NuWRO_2019->SetLineStyle(6);

        total_inel_KE->SetLineColor(2);

        total_inel_KE->Draw("c same");
	xs_GENIE_ha2018->Draw("l same");
	xs_NuWRO_2019->Draw("l same");
	xs_NEUT_2019->Draw("l same");
	xs_GENIE_INCL_pp->Draw("l same");
	xs_GENIE_hN2018->Draw("l same");

	xs_reco_sys->SetLineWidth(2);
	xs_reco_stat->SetLineWidth(2);
	xs_reco_stat->SetMarkerSize(1.3);
	xs_reco_sys->Draw("p same");
	xs_reco_stat->Draw("pz same"); //with 'z' option means get rid of small horizontal bars

        //TLegend *leg_xs = new TLegend(0.35,0.6,0.95,0.85);
        TLegend *leg_xs = new TLegend(0.2,0.7,0.94,0.86);

        leg_xs->SetFillStyle(0);
	leg_xs->SetNColumns(3);
        //leg_xs->AddEntry(xs_reco_stat, "Data(stat.)", "pe");
        //leg_xs->AddEntry(xs_reco_sys, "Data(stat.+sys.)", "pe");
        leg_xs->AddEntry(xs_reco_stat, "Data(stat.+sys.)", "pe");
        //leg_xs->AddEntry(gr_truexs, "MC Truth", "pe");
        leg_xs->AddEntry(total_inel_KE, "Geant4", "l");

        leg_xs->AddEntry(xs_GENIE_INCL_pp, "GENIE INCL++", "l");
        leg_xs->AddEntry(xs_GENIE_ha2018, "GENIE hA2018", "l");
        leg_xs->AddEntry(xs_GENIE_hN2018, "GENIE hN2018", "l");
        leg_xs->AddEntry(xs_NEUT_2019, "NEUT 2019", "l");
        leg_xs->AddEntry(xs_NuWRO_2019, "NuWRO 2019", "l");

        leg_xs->Draw();




   	//pDUNE Logo ----------------------------------------------------------------------
        float logo_y=ymax+15;
        TLatex **txt_pdune1=new TLatex*[1];
        txt_pdune1[0]=new TLatex(xmin, logo_y, Form("#bf{DUNE:ProtoDUNE-SP}"));
        txt_pdune1[0]->SetTextColor(1);
        txt_pdune1[0]->Draw();
        //
        //Beam Logo ---------------------------------------------------------------------
        TLatex **txt_p1=new TLatex*[1];
        //txt_p1[0]=new TLatex(xmax-6.3, logo_y, Form("Protons (1 GeV/c)")); //x:0-20
        //txt_p1[0]=new TLatex(xmax-120, logo_y, Form("Protons (1 GeV/c)")); //x:0-40
        //txt_p1[0]=new TLatex(142, logo_y, Form("Protons (1 GeV/c)")); //x:0-40
        txt_p1[0]=new TLatex(xmax-58, logo_y, Form("Protons (1 GeV/c)")); //x:0-40

        txt_p1[0]->SetTextColor(1);
        txt_p1[0]->SetTextSize(0.05);
        txt_p1[0]->Draw();

	//Preliminary
        TLatex **txt_preliminary=new TLatex*[1];
        txt_preliminary[0]=new TLatex(70, 260, Form("Preliminary"));
        txt_preliminary[0]->SetTextColor(14);
        txt_preliminary[0]->SetTextSize(0.09);
        txt_preliminary[0]->Draw();


        //c_xs->Print(Form("%sxs_reco_stat_sys.eps",out_fig.Data()));
        c_xs->Print(Form("%sxs_reco_stat_sys_zoom.eps",out_fig.Data()));



	//larger range
        TH2D *f2d_xs2=new TH2D("f2d_xs2","",340,xmin,340,ymax,0,ymax);
        f2d_xs2->GetXaxis()->SetTitle("Proton Kinetic Energy [MeV]");
        f2d_xs2->GetYaxis()->SetTitleOffset(1.3);
        f2d_xs2->GetYaxis()->SetTitle("P-Ar inelastic cross section [mb]");
        f2d_xs2->Draw("");
        total_inel_KE->Draw("c same");
	xs_GENIE_ha2018->Draw("l same");
	xs_NuWRO_2019->Draw("l same");
	xs_NEUT_2019->Draw("l same");
	xs_GENIE_INCL_pp->Draw("l same");
	xs_GENIE_hN2018->Draw("l same");
	xs_reco_sys->Draw("p same");
	xs_reco_stat->Draw("pz same"); //with 'z' option means get rid of small horizontal bars
        leg_xs->Draw();
        txt_pdune1[0]->Draw();

        //Beam Logo ---------------------------------------------------------------------
        TLatex **txt_p2=new TLatex*[1];
        txt_p2[0]=new TLatex(340-100, logo_y, Form("Protons (1 GeV/c)")); //x:0-40
        txt_p2[0]->SetTextColor(1);
        txt_p2[0]->SetTextSize(0.05);
        txt_p2[0]->Draw();

	//Preliminary
        TLatex **txt_preliminary2=new TLatex*[1];
        txt_preliminary2[0]=new TLatex(100, 260, Form("Preliminary"));
        txt_preliminary2[0]->SetTextColor(14);
        txt_preliminary2[0]->SetTextSize(0.09);
        txt_preliminary2[0]->Draw();

        c_xs->Print(Form("%sxs_reco_stat_sys.eps",out_fig.Data()));





	//TCanvas *c_=new TCanvas(Form("c_"),"",900, 600);
	//c_->Divide(1,1);




}
