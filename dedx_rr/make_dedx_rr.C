#include "TGraph.h"
#include "TGraphSmooth.h"
#include "TParameter.h"

#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;

void make_dedx_rr() {
		//Read file of csda range versus momentum
		vector<double> predict_resrange;
		vector<double> predict_dedx;
		vector<double> predict_data;
		vector<double> ex;
		double buffer_predict;
		//ifstream f_predict_in("../../proton_mp_dedx_vs_range_0.50_2GeV_KE.txt");
		std::ifstream f_predict_in("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/mcc12_validation/new_alg/proton_mp_dedx_vs_range_1GeV.txt");
		while (f_predict_in>>buffer_predict) { predict_data.push_back(buffer_predict); }
		f_predict_in.close();

		size_t n_predict=2; //2 columns
		for (size_t kk=0; kk<predict_data.size(); kk+=n_predict) {
		predict_resrange.push_back(predict_data[kk]);
		predict_dedx.push_back(predict_data[kk+1]);
		ex.push_back(0);	
		} 
		//TGraphErrors *gr_predict_dedx_resrange = new TGraphErrors(predict_dedx.size(),&predict_resrange[0],&predict_dedx[0],&ex[0],&ex[0]);
		TGraph *gr_predict_dedx_resrange = new TGraph(predict_dedx.size(),&predict_resrange[0],&predict_dedx[0]);
		gr_predict_dedx_resrange->SetName("gr_predict_dedx_resrange");
		gr_predict_dedx_resrange->SetMarkerStyle(20);
		gr_predict_dedx_resrange->SetMarkerSize(0.7);
		gr_predict_dedx_resrange->SetMarkerColor(2);
		gr_predict_dedx_resrange->SetLineColor(2);

		//Read WQ's csda range versus dedx (max. amp of dedx in G4)
		vector<double> wq_resrange;
		vector<double> wq_dedx;
		vector<double> wq_data;
		vector<double> ex_wq;
		double buffer_wq;
		//ifstream f_predict_in("proton_mp_dedx_vs_range_0.50_2GeV_KE.txt");
		std::ifstream f_wq_in("/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/nosce/prod3_test/wenqiang_g4_dedx_rr/proton_mp_dedx_vs_range_geant4_extrc.txt");
		while (f_wq_in>>buffer_wq) { wq_data.push_back(buffer_wq); }
		f_wq_in.close();

		//size_t n_predict=2; //2 columns
		for (size_t kk=0; kk<wq_data.size(); kk+=n_predict) {
			wq_resrange.push_back(0.25+wq_data[kk]);
			wq_dedx.push_back(wq_data[kk+1]);
			ex_wq.push_back(0);
		}
		//TGraphErrors *gr_predict_dedx_resrange = new TGraphErrors(predict_dedx.size(),&predict_resrange[0],&predict_dedx[0],&ex[0],&ex[0]);
		TGraph *gr_wq_dedx_resrange = new TGraph(wq_dedx.size(),&wq_resrange[0],&wq_dedx[0]);
		gr_wq_dedx_resrange->SetName("gr_wq_dedx_resrange");
		gr_wq_dedx_resrange->SetMarkerStyle(20);
		gr_wq_dedx_resrange->SetMarkerSize(0.7);
		gr_wq_dedx_resrange->SetMarkerColor(6);
		gr_wq_dedx_resrange->SetLineColor(6);



		TFile *fout = new TFile("dedx_rr.root","RECREATE");	
		  gr_predict_dedx_resrange->Write();
		  gr_wq_dedx_resrange->Write();
		fout->Close();

}
