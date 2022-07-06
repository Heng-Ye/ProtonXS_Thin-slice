#include "../headers/BasicParameters.h"
#include "/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/headers/BasicAnaFunc.h"

void plotBMRW(TString fdata, TString fmc, TString outpath) {
	//TString rep="trklen";
	//TString x_axis_label="Proton Track Length [cm]";

	TString rep="pcalo";
	TString x_axis_label="Proton Momentum [MeV/c]";

	//plot style --------------------------------------------------------------------//
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	gStyle->SetOptStat(0);

	//read data ---------------------------------------------------------------------//
	TFile *f_data = TFile::Open(fdata.Data());
	TH1D *pbeam_data=(TH1D*)f_data->Get("h1d_pbeam");
	TF1 *pbeam_fit_data=(TF1*)f_data->Get("pbeam_fit");
	TH1D* h1d=(TH1D*)f_data->Get(Form("h1d_%s_stop",rep.Data()));

	TH1D *pbeam_stop_data=(TH1D*)f_data->Get("h1d_pbeam_stop");
	TF1 *pbeam_stop_fit_data=(TF1*)f_data->Get("pbeam_stop_fit");
	int n_data=pbeam_data->Integral(); 
	int n_stop_data=pbeam_stop_data->Integral(); 
	cout<<"n_data:"<<n_data<<endl;
	cout<<"n_stop_data:"<<n_stop_data<<endl;

	pbeam_data->SetLineColor(1); 	      pbeam_data->SetMarkerColor(1);
	pbeam_fit_data->SetLineColor(1);      pbeam_fit_data->SetMarkerColor(1);
	pbeam_stop_data->SetLineColor(1);     pbeam_stop_data->SetMarkerColor(1);
	pbeam_stop_fit_data->SetLineColor(1); pbeam_stop_fit_data->SetMarkerColor(1);

	//read MC -----------------------------------------------------------------------//
	TFile *f_mc = TFile::Open(fmc.Data());
	TH1D *pbeam_mc=(TH1D*)f_mc->Get("h1d_pbeam");
	TF1 *pbeam_fit_mc=(TF1*)f_mc->Get("pbeam_fit");

	TH1D *pbeam_stop_mc=(TH1D*)f_mc->Get("h1d_pbeam_stop");
	TF1 *pbeam_stop_fit_mc=(TF1*)f_mc->Get("pbeam_stop_fit");

	TH1D *pff_stop_mc=(TH1D*)f_mc->Get("h1d_pff_stop");
	TF1 *pff_stop_fit_mc=(TF1*)f_mc->Get("pff_stop_fit");

	int n_mc=pbeam_mc->Integral(); 
	int n_stop_mc=pbeam_stop_mc->Integral(); 
	int n_ff_mc=pff_stop_mc->Integral(); 

	pbeam_mc->Scale((double)n_data/(double)n_mc);
	pff_stop_mc->Scale((double)n_data/(double)n_mc);
	cout<<"pbeam_mc:"<<pbeam_mc->Integral()<<endl;

	pbeam_mc->SetLineColor(2); 	     pbeam_mc->SetMarkerColor(2);
	pbeam_fit_mc->SetLineColor(2);       pbeam_fit_mc->SetMarkerColor(2);
	pbeam_stop_mc->SetLineColor(2);      pbeam_stop_mc->SetMarkerColor(2);
	pbeam_stop_fit_mc->SetLineColor(2);  pbeam_stop_fit_mc->SetMarkerColor(2);
	pff_stop_mc->SetLineColor(2);        pff_stop_mc->SetMarkerColor(2);
	pff_stop_fit_mc->SetLineColor(2);    pff_stop_fit_mc->SetMarkerColor(2);

	//beam momentum reweighting -------------------------------------------------------//
	//get mu & sigma info
	//mu
	TParameter<Int_t> *bm_nmu=(TParameter<Int_t>*)f_mc->Get("bm_nmu");
	Int_t nmu=bm_nmu->GetVal();
	TParameter<Double_t>* bm_dmu=(TParameter<Double_t>*)f_mc->Get("bm_dmu");
	Double_t dmu=bm_dmu->GetVal();
	TParameter<Double_t>* bm_mu_st=(TParameter<Double_t>*)f_mc->Get("bm_mu_st");
	Double_t mu_st=bm_mu_st->GetVal();
	cout<<"mu:"<<nmu<<" "<<dmu<<" "<<mu_st<<endl;

	//sigma
	TParameter<Int_t> *bm_nsigma=(TParameter<Int_t>*)f_mc->Get("bm_nsigma");
	Int_t nsigma=bm_nsigma->GetVal();
	TParameter<Double_t>* bm_dsigma=(TParameter<Double_t>*)f_mc->Get("bm_dsigma");
	Double_t dsigma=bm_dsigma->GetVal();
	TParameter<Double_t>* bm_sigma_st=(TParameter<Double_t>*)f_mc->Get("bm_sigma_st");
	Double_t sigma_st=bm_sigma_st->GetVal();
	cout<<"sigma:"<<nsigma<<" "<<dsigma<<" "<<sigma_st<<endl;





/*
	//Proton Momentum --------------------------------------------------------------//
	TCanvas *c0=new TCanvas("c0","");
	c0->Divide(1,1);
	c0->cd(1);
	TH2D* frame2d=new TH2D("frame2d","", 600, 600, 1400, 600, 0, 600); //zend_2d
	frame2d->SetTitle(";Proton Momentum [MeV/c];");
	frame2d->GetXaxis()->CenterTitle();
	frame2d->Draw();
	pbeam_data->Draw("ep same");
	pbeam_mc->Draw("hist same");

	TLegend *leg0 = new TLegend(0.14,0.65,.6,0.85);
	leg0->SetFillStyle(0);
	leg0->AddEntry(pbeam_data, "Data", "ep");
	leg0->AddEntry(pbeam_mc, "MC", "l");
	leg0->Draw();

	c0->Print(Form("%s/pbeam_data_mc.eps",outpath.Data()));
*/

	//read rw histograms -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
	Int_t n_1d=nmu*nsigma;
	const int n_array=(const int)n_1d;
	TH1D *h1d_rw[n_array];
	vector<double> frac_mom;
	vector<double> mean_mom;
	vector<double> sigma_mom;
	vector<int> selection;

	int cnt_array=0;
	//double m1=pbeam_fit_mc->GetParameter(0); //MC prod p4a [spec]
	//double s1=pbeam_fit_mc->GetParameter(1); //MC prod p4a [spec]
        //double m1=1007.1482; //MC prod4a [spec]
        //double s1=60.703307; //MC prod4a [spec]
        double m1=997.969; //MC prod4a [truth]
        double s1=54.4602; //MC prod4a [truth]

	cout<<"m1:"<<m1<<"|  s1:"<<s1<<endl;
	//cout<<"amp:"<<pbeam_fit_data->GetParameter(2)<<endl;

	//data mu, sigma
	double m1_data=pbeam_fit_data->GetParameter(0); //Data prod4 reco2
	double s1_data=pbeam_fit_data->GetParameter(1); //Data prod4 reco2
	cout<<"m1_data:"<<m1_data<<"|  s1_data:"<<s1_data<<endl;

	//th2f to save best-fit n-sigm region
	TH2F *mu_sigma_chi2_bestfit=new TH2F("mu_sigma_chi2_bestfit","",nmu, mu_st-(float)(nmu-1)*dmu, mu_st-(float)0*dmu,  nsigma, sigma_st-(float)(nsigma-1)*dsigma, sigma_st-(float)(0)*dsigma);
	mu_sigma_chi2_bestfit->GetXaxis()->SetTitle("#mu/#mu_{0}");	
	mu_sigma_chi2_bestfit->GetYaxis()->SetTitle("s/s_{0}");	

	int index_original=0;
	int index_data=0;
	for (int imu=0; imu<nmu; ++imu){ //mu loop
		double frac_mu=mu_st-(double)imu*dmu;
		double mu=m1*frac_mu;
		//cout<<"frac_mu["<<imu<<"]:"<<frac_mu<<endl;

		for (int isigma=0; isigma<nsigma; ++isigma){ //sigma loop
			double frac_sigma=sigma_st-(double)isigma*dsigma;
			double sigma=s1*frac_sigma;

			if (std::abs(mu-m1)<0.0001&&std::abs(sigma-s1)<0.0001) { //no rw
				index_original=cnt_array;
				mu=m1;
				sigma=s1;
				cout<<"index_original:"<<index_original<<endl;
			} //no rw

			//if (std::abs(mu-m1_data)<0.4&&std::abs(sigma-s1_data)<0.4) { //data gaussian
			if (std::abs(mu-m1_data)<0.35&&std::abs(sigma-s1_data)<0.35) { //data gaussian
				index_data=cnt_array;
				cout<<"index_data:"<<index_data<<endl;
			} //data gaussian

			//mean_mom.push_back(mu);
			//sigma_mom.push_back(sigma);
			mean_mom.push_back(frac_mu);
			sigma_mom.push_back(frac_sigma);

			//rw histograms
			h1d_rw[cnt_array]=(TH1D*)f_mc->Get(Form("h1d_%s_rw_%d",rep.Data(),cnt_array));
			cnt_array++;
		} //sigma loop
	} //mu loop
	//read rw histograms -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

	//plot figure -----------------------------------------------------------------------//
	TCanvas *c_1=new TCanvas("c_1",""); c_1->SetName("c_1");
	c_1->Divide(1,1);
	c_1->cd(1);
	//TH2D* frame2dx=new TH2D("frame2dx"," ", 150, 0, 150, 800, 0, 800); //for z2, trklen
	TH2D* frame2dx=new TH2D("frame2dx"," ", 800, 600, 1400, 800, 0, 100+h1d->GetBinContent(h1d->GetMaximumBin())); //for momentum
	frame2dx->GetXaxis()->SetTitle(x_axis_label.Data());
	frame2dx->GetYaxis()->SetTitle(Form("Counts"));
	frame2dx->GetXaxis()->CenterTitle();
	frame2dx->GetYaxis()->CenterTitle();

	//data
	frame2dx->Draw(); 
	h1d->Draw("ep same");

	//chi2 distributions -------------------------------------------------//
	//preparation for data and mc
	vector<double> D; //data
	vector<double> er_D; //error of data
	//h1d->Sumw2();
	for (int k=1; k<=h1d->GetNbinsX(); k++){
		D.push_back(h1d->GetBinContent(k));
		er_D.push_back(sqrt(h1d->GetBinContent(k)));
	}

	//vector<double> chi2_ndf;
	vector<double> chi2;
	vector<double> Dchi2;
	vector<double> ndf;
	double min_chi2=9999999999999;
	double max_chi2=-9999999999999;
	size_t min_index=-1;
	size_t max_index=-1;
	for (size_t i=0; i<(size_t)n_1d; i++) { //loop over all array histograms
		//h1d_rw[i]->Sumw2();
		//normalization constant
		double norm=(double)n_stop_data/(double)(h1d_rw[i]->Integral());

		//chi2/ndf calculation
		vector<double> MC; //MC
		vector<double> er_MC; //error of MC
		for (int k=1; k<=h1d_rw[i]->GetNbinsX(); k++) { //loop over all bin contents of each histogram
			MC.push_back(norm*h1d_rw[i]->GetBinContent(k));
			//er_MC.push_back(norm*h1d_rw[i]->GetBinError(k));
			er_MC.push_back(norm*sqrt(h1d_rw[i]->GetBinContent(k)));
		} //loop over all bin contents of each histogram

		//normalize histogram
		h1d_rw[i]->Scale(norm);
		//h1d_rw[i]->SetLineColor(38);
		h1d_rw[i]->SetLineColor(3);

		//for (int k=1; k<=h1d_rw[i]->GetNbinsX(); k++) { //loop over all bin contents of each histogram
		//MC.push_back(h1d_rw[i]->GetBinContent(k));
		//er_MC.push_back(h1d_rw[i]->GetBinError(k));
		//} //loop over all bin contents of each histogram

		//double chi2_=chi2ndf_data_mc(D, er_D, 1, MC, er_MC);
		//double chi2_=neyman_chi2_data_mc(D, er_D, MC, er_MC);
		//double chi2_=pearson_chi2_data_mc(D, er_D, MC, er_MC);
		double chi2_=ml_data_mc(D, er_D, MC, er_MC);
		//double chi2ndf=chi2/(double)D.size();
		chi2.push_back(chi2_);
		//std::cout<<"index:"<<i<<" chi2:"<<chi2_<<std::endl;
		if (min_chi2>chi2_) {
			min_chi2=chi2_;
			min_index=i;
		}
		if (max_chi2<chi2_) {
			max_chi2=chi2_;
			max_index=i;
		}
	} //loop over all array histograms

	//Delta chi2
	for (size_t t=0; t<chi2.size(); ++t) {
		Dchi2.push_back(chi2.at(t)-chi2.at(min_index));	
	}

	//2 1D projections to find best-fit n-sigma region --------------------------------------------------------------------------//
	//minimum chi2 as func. of x
	vector<double> x_mom;
	vector<double> x_chi2;
	vector<int> index_x_chi2;

	int cnt_=0;
	for (int imu=0; imu<nmu; ++imu){ //mu loop
		double frac_mu=mu_st-(double)imu*dmu;

		double min_chi2_=99999;
		int tmp_index=0;		
		for (int isigma=0; isigma<nsigma; ++isigma){ //sigma loop
			double frac_sigma=sigma_st-(double)isigma*dsigma;

			if (mean_mom.at(cnt_)==frac_mu) {
				if (chi2.at(cnt_)<min_chi2_) {
					min_chi2_=chi2.at(cnt_);
					tmp_index=cnt_;
				}
			}
			cnt_++;
		} //sigma loop
		x_chi2.push_back(min_chi2_);
		index_x_chi2.push_back(tmp_index);
		x_mom.push_back(frac_mu);
	} //mu loop

	//setup best-fit region (1-s, 2-s, 3-s)
	int index_min_x_chi2=std::distance(x_chi2.begin(), std::min_element(x_chi2.begin(), x_chi2.end())); //index of min. chi2 in x
	double min_x_chi2=*min_element(x_chi2.begin(), x_chi2.end()); //global min. chi2 in x
	double min_x_mom=x_mom.at(index_min_x_chi2); //x(global min. chi2)
	double max_x_chi2=*max_element(x_chi2.begin(), x_chi2.end()); //max. chi2 in x
	cout<<"Min chi^2 value and index of (x-axis):"<<min_x_chi2<<" ["<<index_min_x_chi2<<"]="<<min_x_mom<<endl;

	const int n_sl=3; //total number of sigma level
	double dchi2[n_sl]={1, 4, 9}; //associtated delta chi^2 values in 1-sigma, 2-sigma, 3-sigma

	double set_x_chi2[n_sl]; //best-fit chi2 
	double get_xl_chi2[n_sl];
	double get_xr_chi2[n_sl];
	double get_xl[n_sl]; 
	double get_xr[n_sl]; 

	for (int ni=0; ni<n_sl; ni++) { //chi2-n-sigma loop
		set_x_chi2[ni]=min_x_chi2+dchi2[ni];

		get_xl_chi2[ni]=min_x_chi2; //set min chi2 as default
		get_xr_chi2[ni]=min_x_chi2; //set min chi2 as default

		get_xl[ni]=min_x_mom; //set min chi2 as default
		get_xr[ni]=min_x_mom; //set min chi2 as default
	} //chi2-n-sigma loop

	//started from minimum, find associtated points in the given range 
	//get points on the right
	for (int ni=0; ni<n_sl; ni++) { //n-sigma loop
		for (int ir=index_min_x_chi2; ir>=0; ir--) { //sweep all the points on the right
			double xi_chi2=x_chi2.at(ir);

			if (set_x_chi2[ni]<xi_chi2) {
				get_xr_chi2[ni]=xi_chi2;
				get_xr[ni]=x_mom.at(ir);
				break;
			}
		} //sweep all the point on the left
	} //n-sigma loop
	cout<<"x 1-sigma set at "<<set_x_chi2[0]<<"-->  get right (chi2, mom):("<<get_xr_chi2[0]<<", "<<get_xr[0]<<")"<<endl;	
	cout<<"x 2-sigma set at "<<set_x_chi2[1]<<"-->  get right (chi2, mom):("<<get_xr_chi2[1]<<", "<<get_xr[1]<<")"<<endl;	
	cout<<"x 3-sigma set at "<<set_x_chi2[2]<<"-->  get right (chi2, mom):("<<get_xr_chi2[2]<<", "<<get_xr[2]<<")"<<endl;	

	//get points on the left
	for (int ni=0; ni<n_sl; ni++) { //n-sigma loop
		for (int il=index_min_x_chi2; il<x_chi2.size(); il++) { //sweep all the points on the left
			double xi_chi2=x_chi2.at(il);

			if (set_x_chi2[ni]<xi_chi2) {
				get_xl_chi2[ni]=xi_chi2;
				get_xl[ni]=x_mom.at(il);
				break;
			}
		} //sweep all the point on the left
	} //n-sigma loop
	cout<<"x 1-sigma set at "<<set_x_chi2[0]<<"-->  get left (chi2, mom):("<<get_xl_chi2[0]<<", "<<get_xl[0]<<")"<<endl;	
	cout<<"x 2-sigma set at "<<set_x_chi2[1]<<"-->  get left (chi2, mom):("<<get_xl_chi2[1]<<", "<<get_xl[1]<<")"<<endl;	
	cout<<"x 3-sigma set at "<<set_x_chi2[2]<<"-->  get left (chi2, mom):("<<get_xl_chi2[2]<<", "<<get_xl[2]<<")"<<endl;	

	cout<<"========================================================================"<<endl;
	cout<<"Best-fit n-sigma in x:"<<endl;
	cout<<"x_mom (1-sigma):"<<min_x_mom<<" -"<<(min_x_mom-get_xl[0])<<" +"<<(get_xr[0]-min_x_mom)<<endl;
	cout<<"x_mom (2-sigma):"<<min_x_mom<<" -"<<(min_x_mom-get_xl[1])<<" +"<<(get_xr[1]-min_x_mom)<<endl;
	cout<<"x_mom (3-sigma):"<<min_x_mom<<" -"<<(min_x_mom-get_xl[2])<<" +"<<(get_xr[2]-min_x_mom)<<endl;
	cout<<"========================================================================"<<endl;

	//plot figure
	float fx2d_ymin=0;
	float fx2d_ymax=max_x_chi2+20;
	int ny_fx2d=100;

	float fx2d_xmin=0.975;
	float fx2d_xmax=1.01;
	int nx_fx2d=100;

	TH2D* f2d_x=new TH2D("f2d_x","",nx_fx2d, fx2d_xmin, fx2d_xmax, ny_fx2d, fx2d_ymin,fx2d_ymax);
	f2d_x->GetXaxis()->SetTitle("#mu/#mu_{0}");
	f2d_x->GetXaxis()->CenterTitle();
	f2d_x->GetYaxis()->SetTitle("Minimum #chi^{2}");
	f2d_x->GetYaxis()->CenterTitle();	

	TCanvas *c_chi2_x=new TCanvas("c_chi2_x",""); c_chi2_x->SetName("c_chi2_x");
	c_chi2_x->Divide(1,1);
	c_chi2_x->cd(1);
	f2d_x->Draw();

	TGraph *gr_chi2_x = new TGraph(x_mom.size(), &x_mom.at(0), &x_chi2.at(0));
	gr_chi2_x->Draw("psame");

	TLine **x_hor=new TLine*[n_sl]; //horizontal
	TLine **x_lv=new TLine*[n_sl]; //left
	TLine **x_rv=new TLine*[n_sl]; //right
	TLine* x_cen=new TLine(min_x_mom, fx2d_ymin, min_x_mom, min_x_chi2); //central
	x_cen->SetLineStyle(1); 
	x_cen->SetLineColor(2);
	x_cen->SetLineWidth(2);

	c_chi2_x->cd(1);
	x_cen->Draw("same");

	for (int ll=0; ll<n_sl; ll++) {
		x_hor[ll]=new TLine(get_xl[ll], set_x_chi2[ll], get_xr[ll], set_x_chi2[ll]);
		x_lv[ll]=new TLine(get_xl[ll], fx2d_ymin, get_xl[ll], get_xl_chi2[ll]);
		x_rv[ll]=new TLine(get_xr[ll], fx2d_ymin, get_xr[ll], get_xr_chi2[ll]);

		x_hor[ll]->SetLineStyle(2);
		x_lv[ll]->SetLineStyle(2);
		x_rv[ll]->SetLineStyle(2);

		x_hor[ll]->SetLineWidth(2);
		x_lv[ll]->SetLineWidth(2);
		x_rv[ll]->SetLineWidth(2);

		if (ll==0) {
			x_hor[ll]->SetLineColor(2);
			x_lv[ll]->SetLineColor(2);
			x_rv[ll]->SetLineColor(2);
		}

		if (ll==1) {
			x_hor[ll]->SetLineColor(4);
			x_lv[ll]->SetLineColor(4);
			x_rv[ll]->SetLineColor(4);
		}

		if (ll==2) {
			x_hor[ll]->SetLineColor(3);
			x_lv[ll]->SetLineColor(3);
			x_rv[ll]->SetLineColor(3);
		}

		c_chi2_x->cd(1);
		x_hor[ll]->Draw("same");
		x_lv[ll]->Draw("same");
		x_rv[ll]->Draw("same");
	}

	TLegend *legx=new TLegend(0.16,0.45,0.85,0.85);
	legx->SetFillColor(0);
	legx->SetFillStyle(0);
	TLegendEntry* llx[10];
	llx[0]=legx->AddEntry((TObject*)0, Form("Central value (#chi^{2}_{min}): %.3f", min_x_mom),"");
	for (int j=0; j<n_sl; ++j) {
		llx[j]=legx->AddEntry((TObject*)0, Form("%d-#sigma: ^{+%.3f}_{-%.3f} (#chi^{2}_{min}+%.0f)", j+1, get_xr[j]-min_x_mom, min_x_mom-get_xl[j], dchi2[j]),"");
		if (j==0) llx[j]->SetTextColor(2);
		if (j==1) llx[j]->SetTextColor(4);
		if (j==2) llx[j]->SetTextColor(3);
	}
	c_chi2_x->cd(1);
	legx->Draw();

	cout<<"\n\n\n"<<endl;


	//minimum chi2 as func. of y
	vector<double> y_mom;
	vector<double> y_chi2;
	vector<int> index_y_chi2;
	for (int ksigma=0; ksigma<nsigma; ++ksigma){ //sigma loop
		double frac_ksigma=sigma_st-(double)ksigma*dsigma;
		double min_chi2_=99999;
		int tmp_index=0;
		//sweap through all the points in the array --------------//
		int cnt_y=0;
		for (int imu=0; imu<nmu; ++imu){ //mu loop
			double frac_mu=mu_st-(double)imu*dmu;
			for (int isigma=0; isigma<nsigma; ++isigma){ //sigma loop
				double frac_sigma=sigma_st-(double)isigma*dsigma;

				if (frac_ksigma==frac_sigma) {
					if (chi2.at(cnt_y)<min_chi2_) {
						min_chi2_=chi2.at(cnt_y);
						tmp_index=cnt_y;
					}
				}
				cnt_y++;
			} //sigma loop
		} //mu loop
		//sweap through all the points in the array --------------//
		y_chi2.push_back(min_chi2_);
		index_y_chi2.push_back(tmp_index);
		y_mom.push_back(frac_ksigma);
	} //sigma loop

	//setup best-fit region (1-s, 2-s, 3-s)
	int index_min_y_chi2=std::distance(y_chi2.begin(), std::min_element(y_chi2.begin(), y_chi2.end())); //index of min. chi2 in y
	double min_y_chi2=*min_element(y_chi2.begin(), y_chi2.end()); //global min. chi2 in y
	double min_y_mom=y_mom.at(index_min_y_chi2); //y(global min. chi2)
	double max_y_chi2=*max_element(y_chi2.begin(), y_chi2.end()); //max. chi2 in y
	cout<<"Min chi^2 value and index of (y-axis):"<<min_y_chi2<<" ["<<index_min_y_chi2<<"]="<<min_y_mom<<endl;

	double set_y_chi2[n_sl]; //best-fit chi2 
	double get_yl_chi2[n_sl];
	double get_yr_chi2[n_sl];
	double get_yl[n_sl]; 
	double get_yr[n_sl]; 

	for (int ni=0; ni<n_sl; ni++) { //chi2-n-sigma loop
		set_y_chi2[ni]=min_y_chi2+dchi2[ni];

		get_yl_chi2[ni]=min_y_chi2; //set min chi2 as default
		get_yr_chi2[ni]=min_y_chi2; //set min chi2 as default

		get_yl[ni]=min_y_mom; //set min chi2 as default
		get_yr[ni]=min_y_mom; //set min chi2 as default
	} //chi2-n-sigma loop

	//started from minimum, find associtated points in the given range 
	//get points on the right
	for (int ni=0; ni<n_sl; ni++) { //n-sigma loop
		for (int ir=index_min_y_chi2; ir>=0; ir--) { //sweep all the points on the right
			double yi_chi2=y_chi2.at(ir);
			if (ni==0) cout<<"ir:"<<ir<<"yi_chi2:"<<yi_chi2<<" ; set_y_chi2:"<<set_y_chi2[ni]<<endl;

			if (set_y_chi2[ni]<yi_chi2) {
				get_yr_chi2[ni]=yi_chi2;
				get_yr[ni]=y_mom.at(ir);
				break;
			}
		} //sweep all the point on the left

		if (get_yr[ni]==min_y_mom) { 
			get_yr[ni]=y_mom.at(0);
		}
	} //n-sigma loop

	cout<<"y 1-sigma set at "<<set_y_chi2[0]<<"-->  get right (chi2, mom):("<<get_yr_chi2[0]<<", "<<get_yr[0]<<")"<<endl;	
	cout<<"y 2-sigma set at "<<set_y_chi2[1]<<"-->  get right (chi2, mom):("<<get_yr_chi2[1]<<", "<<get_yr[1]<<")"<<endl;	
	cout<<"y 3-sigma set at "<<set_y_chi2[2]<<"-->  get right (chi2, mom):("<<get_yr_chi2[2]<<", "<<get_yr[2]<<")"<<endl;	

	//get points on the left
	for (int ni=0; ni<n_sl; ni++) { //n-sigma loop
		for (int il=index_min_y_chi2; il<y_chi2.size(); il++) { //sweep all the points on the left
			double yi_chi2=y_chi2.at(il);

			if (set_y_chi2[ni]<yi_chi2) {
				get_yl_chi2[ni]=yi_chi2;
				get_yl[ni]=y_mom.at(il);
				break;
			}
		} //sweep all the point on the left
	} //n-sigma loop
	cout<<"y 1-sigma set at "<<set_y_chi2[0]<<"-->  get left (chi2, mom):("<<get_yl_chi2[0]<<", "<<get_yl[0]<<")"<<endl;	
	cout<<"y 2-sigma set at "<<set_y_chi2[1]<<"-->  get left (chi2, mom):("<<get_yl_chi2[1]<<", "<<get_yl[1]<<")"<<endl;	
	cout<<"y 3-sigma set at "<<set_y_chi2[2]<<"-->  get left (chi2, mom):("<<get_yl_chi2[2]<<", "<<get_yl[2]<<")"<<endl;	

	cout<<"========================================================================"<<endl;
	cout<<"Best-fit n-sigma in y:"<<endl;
	cout<<"y_mom (1-sigma):"<<min_y_mom<<" -"<<(min_y_mom-get_yl[0])<<" +"<<(get_yr[0]-min_y_mom)<<endl;
	cout<<"y_mom (2-sigma):"<<min_y_mom<<" -"<<(min_y_mom-get_yl[1])<<" +"<<(get_yr[1]-min_y_mom)<<endl;
	cout<<"y_mom (3-sigma):"<<min_y_mom<<" -"<<(min_y_mom-get_yl[2])<<" +"<<(get_yr[2]-min_y_mom)<<endl;
	cout<<"========================================================================"<<endl;

	//plot figure
	float fy2d_ymin=0;
	float fy2d_ymax=max_y_chi2+20;
	int ny_fy2d=50;

	//float fy2d_xmin=0.9;
	//float fy2d_xmax=1.4;
	float fy2d_xmin=0.9;
	float fy2d_xmax=1.6;
	int nx_fy2d=nsigma;

	TH2D* f2d_y=new TH2D("f2d_y","",nx_fy2d, fy2d_xmin, fy2d_xmax, ny_fy2d, fy2d_ymin,fy2d_ymax);
	f2d_y->GetXaxis()->SetTitle("#sigma/#sigma_{0}");
	f2d_y->GetXaxis()->CenterTitle();
	f2d_y->GetYaxis()->SetTitle("Minimum #chi^{2}");
	f2d_y->GetYaxis()->CenterTitle();

	TCanvas *c_chi2_y=new TCanvas("c_chi2_y",""); c_chi2_y->SetName("c_chi2_y");
	c_chi2_y->Divide(1,1);
	c_chi2_y->cd(1);
	f2d_y->Draw();

	TGraph *gr_chi2_y = new TGraph(y_mom.size(), &y_mom.at(0), &y_chi2.at(0));
	gr_chi2_y->Draw("p same");

	TLine **y_hor=new TLine*[n_sl]; //horizontal
	TLine **y_lv=new TLine*[n_sl]; //left
	TLine **y_rv=new TLine*[n_sl]; //right
	TLine* y_cen=new TLine(min_y_mom, fy2d_ymin, min_y_mom, min_y_chi2); //central
	y_cen->SetLineStyle(1); 
	y_cen->SetLineColor(2);
	y_cen->SetLineWidth(2);

	c_chi2_y->cd(1);
	y_cen->Draw("same");

	for (int ll=0; ll<n_sl; ll++) {
		y_hor[ll]=new TLine(get_yl[ll], set_y_chi2[ll], get_yr[ll], set_y_chi2[ll]);
		y_lv[ll]=new TLine(get_yl[ll], fy2d_ymin, get_yl[ll], get_yl_chi2[ll]);
		y_rv[ll]=new TLine(get_yr[ll], fy2d_ymin, get_yr[ll], get_yr_chi2[ll]);

		y_hor[ll]->SetLineStyle(2);
		y_lv[ll]->SetLineStyle(2);
		y_rv[ll]->SetLineStyle(2);

		y_hor[ll]->SetLineWidth(2);
		y_lv[ll]->SetLineWidth(2);
		y_rv[ll]->SetLineWidth(2);

		if (ll==0) {
			y_hor[ll]->SetLineColor(2);
			y_lv[ll]->SetLineColor(2);
			y_rv[ll]->SetLineColor(2);
		}

		if (ll==1) {
			y_hor[ll]->SetLineColor(4);
			y_lv[ll]->SetLineColor(4);
			y_rv[ll]->SetLineColor(4);
		}

		if (ll==2) {
			y_hor[ll]->SetLineColor(3);
			y_lv[ll]->SetLineColor(3);
			y_rv[ll]->SetLineColor(3);
		}

		c_chi2_y->cd(1);
		y_hor[ll]->Draw("same");
		y_lv[ll]->Draw("same");
		y_rv[ll]->Draw("same");
	}

	TLegend *legy=new TLegend(0.12,0.4,0.85,0.85);
	legy->SetFillColor(0);
	legy->SetFillStyle(0);
	TLegendEntry* lly[10];
	lly[0]=legy->AddEntry((TObject*)0, Form("Central value (#chi^{2}_{min}): %.3f", min_y_mom),"");
	for (int j=0; j<n_sl; ++j) {
		lly[j]=legy->AddEntry((TObject*)0, Form("%d-#sigma: ^{+%.3f}_{-%.3f} (#chi^{2}_{min}+%.0f)", j+1, get_yr[j]-min_y_mom, min_y_mom-get_yl[j], dchi2[j]),"");
		if (j==0) lly[j]->SetTextColor(2);
		if (j==1) lly[j]->SetTextColor(4);
		if (j==2) lly[j]->SetTextColor(3);
	}
	c_chi2_y->cd(1);
	legy->Draw();
	//2 1D projections to find best-fit n-sigma region --------------------------------------------------------------------------//

	//label best-fit area ---------------------------------------------------------//
	vector<int> selection_bestfit;
	int n_sigma=1;
	//int n_sigma=1;
	vector<double> axis_elast_bestfit;
	vector<double> axis_react_bestfit;
	vector<double> chi2_bestfit;
	int index_boundary=0;
	for (size_t i=0; i<(size_t)n_1d; i++) { //loop over all array histograms
		double chi2__=chi2.at(i);

		//best-fit n-sigma region
		if (chi2__<=(min_chi2+(double)n_sigma)) { //best-fit n-sigma
			selection_bestfit.push_back(1);

			//h1d_rw[i]->Draw("hist same");

			//axis_elast_bestfit.push_back(axis_elast.at(i));
			//axis_react_bestfit.push_back(axis_react.at(i));
			chi2_bestfit.push_back(chi2.at(i));
			index_boundary=i;

			mu_sigma_chi2_bestfit->Fill(mean_mom.at(i), sigma_mom.at(i), chi2.at(i));
		}
		else {
			selection_bestfit.push_back(0);
			//mu_sigma_chi2_bestfit->Fill(mean_mom.at(i), sigma_mom.at(i), chi2.at(i));
		}
	}
	c_1->cd(1);

	h1d_rw[index_original]->SetLineColor(4);
	h1d_rw[index_original]->SetLineWidth(4);

	h1d_rw[min_index]->SetLineColor(2);
	h1d_rw[min_index]->SetLineWidth(4);

	h1d_rw[index_data]->SetLineColor(3);
	//h1d_rw[index_data]->Draw("hist same"); //index data
	h1d_rw[index_original]->Draw("hist same"); //no rw
	h1d_rw[min_index]->Draw("hist same"); //min. chi2
	//h1d_rw[5782]->SetLineColor(kViolet+1); //min. chi2 from z2-z0	
	//h1d_rw[5782]->Draw("hist same"); //min. chi2 from z2-z0
	//h1d_rw[max_index]->Draw("hist same"); //max. chi2
	//h1d->Draw("ep same");

	TLegend *leg=new TLegend(0.128677,0.669353,0.839161,0.881581); //z2, trklen
	//TLegend *leg=new TLegend(0.13,0.58,0.84,0.88); //z2, trklen
	//TLegend *leg=new TLegend(0.18,0.5,0.9,0.79); //b2, liny
	//TLegend *leg=new TLegend(0.25,0.6,0.9,0.85); //b2, logy
	leg->SetFillColor(0);
	leg->SetFillStyle(0);

	TLegendEntry* ll[10];
	int NDF=0;
	for (int uu=1; uu<(int)h1d->GetNbinsX(); uu++) {
		int tmp_n=h1d->GetBinContent(uu);
		if (tmp_n>0) NDF++; 
	}
	ll[0]=leg->AddEntry(h1d,"Data","ep"); 
	ll[1]=leg->AddEntry(h1d_rw[index_original],Form("MC [Default]: #chi^{2}/ndf=%.1f/%d",chi2.at(index_original),NDF),"l"); 
	//ll[1]=leg->AddEntry(h1d_rw[index_original],Form("MC [No RW]: #chi^{2}/ndf=%.1f/%d",chi2.at(index_original),h1d->GetNbinsX()),"l"); 
	//ll[2]=leg->AddEntry(h1d_rw[index_data],Form("MC [RW using data (#mu,#sigma)]: #chi^{2}/ndf=%.1f/%d",chi2.at(index_data),h1d->GetNbinsX()-2),"l"); 
	//ll[3]=leg->AddEntry(h1d_rw[min_index],Form("MC [RW with min. #chi^{2}]: #chi^{2}/ndf=%.1f/%d",chi2.at(min_index),h1d->GetNbinsX()-2),"l"); 
	ll[3]=leg->AddEntry(h1d_rw[min_index],Form("MC [Weighted]: #chi^{2}/ndf=%.1f/%d",chi2.at(min_index),NDF-2),"l"); 
	//ll[4]=leg->AddEntry(h1d_rw[5782],Form("MC [RW with min. #chi^{2} from (z_{2}-z_{0})/cos#theta]: #chi^{2}/ndf=%.1f/%d",chi2.at(5782),h1d->GetNbinsX()-2),"l"); 
	c_1->cd(1);
	leg->Draw();
	//h1d_tot->SetFillColor(3);
	//h1d_tot->SetLineColor(3);
	//ll[4]=leg->AddEntry(h1d_rw[index_boundary-1],Form("MC(RW #chi^{2} between min. and max.)"),"f");
	//ll[4]=leg->AddEntry(h1d_tot,Form("MC(RW #chi^{2} between min. and max.)"),"f");
	//ll[4]=leg->AddEntry(h1d_tot,Form("MC(RW #chi^{2} within min. #chi^{2}+%d)",n_sigma),"f");

	cout<<"min_index:"<<min_index<<endl;

	//pDUNE Logo
	TLatex **txt_pdune=new TLatex*[1];
	txt_pdune[0]=new TLatex(0.002, 806, Form("#bf{DUNE:ProtoDUNE-SP}"));
	txt_pdune[0]->SetTextColor(1);
	txt_pdune[0]->Draw(); 

	//Beam Logo
	TLatex **txt_p=new TLatex*[1];
	txt_p[0]=new TLatex(77.1929,806, Form("Stopping Protons (1 GeV/c)"));
	txt_p[0]->SetTextColor(1);
	//txt_pdune[0]->SetTextSize(0.07);
	txt_p[0]->Draw();

	//2d (mu, sigma) space, define contour ----------------------------------------------//
	vector<double> x_mu_contour;	
	vector<double> y_sigma_contour;	
	vector<double> chi2_contour;	
	int cnt_i=0;
	//nmu,start_rw-dmu*(double)int,start_rw,nsigma,sigma_st-(double)nsigma*dsigma,sigma_st  
	//TProfile *prof = new TProfile("prof", "Prof",nsigma,start_rw-dmu*(double)dmu,start_rw,nsigma,sigma_st-(double)nsigma*dsigma,sigma_st); //Profile
	for (int imu=0; imu<nmu; ++imu){ //mu loop
		//double frac_mu=mu_st-(double)imu*dmu;
		//double mu=m1*frac_mu;
		//cout<<"frac_mu["<<imu<<"]:"<<frac_mu<<endl;

		double bc_max_sigma=-999;
		double bc_min_sigma=999;
		int index_bc_max=0;
		int index_bc_min=0;
		int cnt_slice=0;
		for (int isigma=0; isigma<nsigma; ++isigma){ //sigma loop
			//double frac_sigma=sigma_st-(double)isigma*dsigma;
			//double sigma=s1*frac_sigma;
			double sigma_=sigma_mom.at(cnt_i);

			if (selection_bestfit.at(cnt_i)==1) { //best-fit area
				cnt_slice++;	
				if (sigma_>bc_max_sigma) { 
					bc_max_sigma=sigma_;
					index_bc_max=cnt_i;
				}
				if (sigma_<bc_min_sigma) { 
					bc_min_sigma=sigma_;
					index_bc_min=cnt_i;
				}
			} //best-fit area
			cnt_i++;
		} //sigma loop

		if (cnt_slice>0) { //min. chi2 within
			if (bc_max_sigma!=-999&&bc_min_sigma!=999) {
				if (bc_max_sigma==bc_min_sigma) {
					//x_mu_contour.push_back(mean_mom.at(cnt_i));
					//y_sigma_contour.push_back(sigma_mom.at(cnt_i));
				}
				if (bc_max_sigma!=bc_min_sigma) {
					//x_mu_contour.push_back(mean_mom.at(index_bc_max));
					//y_sigma_contour.push_back(sigma_mom.at(index_bc_max));

					//x_mu_contour.push_back(mean_mom.at(index_bc_min));
					//y_sigma_contour.push_back(sigma_mom.at(index_bc_min));
					//prof->Fill(mean_mom.at(index_bc_max),sigma_mom.at(index_bc_max));
					//prof->Fill(mean_mom.at(index_bc_min),sigma_mom.at(index_bc_min));
				}
			}
		} //min chi2 within
	} //mu loop
	//2d (mu, sigma) space, define contour ----------------------------------------------//

	//2d chi^2 ------------------------------------------------------------------------------------------//
	TCanvas *c_2=new TCanvas("c_2",""); c_2->SetName("c_2");
	//gStyle->SetPalette(53);
	c_2->Divide(1,1);
	c_2->cd(1);
	c_2->cd(1)->SetLogz();
	gStyle->SetPalette(104);
	TGraph2D *gr_chi2 = new TGraph2D(mean_mom.size(), &mean_mom.at(0), &sigma_mom.at(0), &Dchi2.at(0));
	//gr_chi2->SetTitle("#chi^{2}; #mu [MeV/c]; #sigma [MeV/c]");
	//gr_chi2->SetTitle("#chi^{2}; #mu/#mu_{0}; #sigma/#sigma_{0}");
	//gr_chi2->SetTitle("#Delta#chi^{2}=-2#Deltaln#lambda; #mu/#mu_{0}; #sigma/#sigma_{0}");
	gr_chi2->SetTitle("#Delta#chi^{2}=-2#Deltaln#lambda; #mu/#mu_{0}; s/s_{0}");

	gr_chi2->GetHistogram()->GetXaxis()->CenterTitle();
	gr_chi2->GetHistogram()->GetYaxis()->CenterTitle();
	gr_chi2->GetHistogram()->GetZaxis()->CenterTitle();

	//gr_chi2->GetHistogram()->GetXaxis()->SetTitleOffset(1.5);
	//gr_chi2->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
	//gr_chi2->GetHistogram()->GetZaxis()->SetTitleOffset(1.1);

	//gr_chi2->Draw("surf1");
	gr_chi2->Draw("colz");

	//minimum chi^2
	//TMarker *minchi2 = new TMarker(mean_mom.at(min_index), sigma_mom.at(min_index), 43);
	TMarker *minchi2 = new TMarker(mean_mom.at(min_index), sigma_mom.at(min_index), 47);
	minchi2->SetMarkerColor(2);
	minchi2->Draw();

	//original
	TMarker *origchi2 = new TMarker(mean_mom.at(index_original), sigma_mom.at(index_original), 41);
	origchi2->SetMarkerColor(4);
	//origchi2->Draw();

	//aussume data mu, sigma
	TMarker *datachi2 = new TMarker(mean_mom.at(index_data), sigma_mom.at(index_data), 71);
	datachi2->SetMarkerColor(3);
	//datachi2->Draw();

	//selected area
	//TGraph *gr_margin = new TGraph(x_mu_contour.size(), &x_mu_contour.at(0), &y_sigma_contour.at(0));
	//gr_margin->Draw("al same");
	//x_mu_contour
	//prof->Draw("HIST L SAME");	

	std::cout<<"\n\nmin chi2 index:"<<min_index<<"  (mu,sigma):("<<mean_mom.at(min_index)<<","<<sigma_mom.at(min_index)<<")"<<std::endl;
	std::cout<<"original chi2 index:"<<index_original<<"  (mu,sigma):("<<mean_mom.at(index_original)<<","<<sigma_mom.at(index_original)<<")"<<std::endl;
	std::cout<<"data chi2 index:"<<index_data<<"  (mu,sigma):("<<mean_mom.at(index_data)<<","<<sigma_mom.at(index_data)<<")"<<std::endl;

	gPad->Update();

	//plot figure -----------------------------------------------------------------------//
	TCanvas *c_ns=new TCanvas("c_ns",""); c_ns->SetName("c_ns");
	c_ns->Divide(1,1);
	c_ns->cd(1);

	//float contour_level=min_chi2+(float)n_sigma;
	//Double_t contour_level[1]={min_chi2+(float)n_sigma};
	//mu_sigma_chi2_bestfit->SetContourLevel(1, contour_level[0]); //total number of contour level, value
	//mu_sigma_chi2_bestfit->Draw("colz");

	const int n_level=2;
	double contours[n_level];
	//contours[0] = min_chi2;
	contours[0] = min_chi2;
	contours[1] = min_chi2+(float)n_sigma;
	//contours[2] = 900;
	mu_sigma_chi2_bestfit->SetContour(2, contours);

	//Int_t colors[n_level] = {kRed};
	//gStyle->SetPalette(1,colors);
	//gStyle->SetNumberContours(1);
	mu_sigma_chi2_bestfit->SetLineColor(2);
	mu_sigma_chi2_bestfit->Draw("cont3");

	//c_2->cd(1);	
	//mu_sigma_chi2_bestfit->Draw("cont3 same");

	//best-fit lines 
	float xmin_mu=mu_st-(float)(nmu-1)*dmu;
	float xmax_mu=mu_st;

	TLine **hor_dn=new TLine*[1]; //dn
	TLine **hor_up=new TLine*[1]; //dn
	hor_dn[0]=new TLine(xmin_mu,get_yl[0],xmax_mu,get_yl[0]);
	hor_up[0]=new TLine(xmin_mu,get_yr[0],xmax_mu,get_yr[0]);
	hor_dn[0]->SetLineColor(2);
	hor_up[0]->SetLineColor(2);
	hor_dn[0]->SetLineStyle(2);
	hor_up[0]->SetLineStyle(2);
	hor_dn[0]->SetLineWidth(2);
	hor_up[0]->SetLineWidth(2);

	//hor_dn[0]->Draw(" same");
	//hor_up[0]->Draw(" same");


	float ymin_sigma=sigma_st-(float)(nsigma-1)*dsigma;
	float ymax_sigma=sigma_st;

	TLine **ver_l=new TLine*[1]; //dn
	TLine **ver_r=new TLine*[1]; //dn
	ver_l[0]=new TLine(get_xl[0],ymin_sigma,get_xl[0],ymax_sigma);
	ver_r[0]=new TLine(get_xr[0],ymin_sigma,get_xr[0],ymax_sigma);
	ver_l[0]->SetLineColor(2);
	ver_r[0]->SetLineColor(2);
	ver_l[0]->SetLineStyle(2);
	ver_r[0]->SetLineStyle(2);
	ver_l[0]->SetLineWidth(2);
	ver_r[0]->SetLineWidth(2);

	//ver_l[0]->Draw(" same");
	//ver_r[0]->Draw(" same");


	//Pick up some data points for systematic study ----------------------------------------//
	double mean_mom_minchi2=mean_mom.at(min_index);
	double sigma_mom_minchi2=sigma_mom.at(min_index);
	double chi2_min=chi2.at(min_index);
	vector<int> key_sel0;
	vector<double> mean_sel0;
	vector<double> sigma_sel0;
	for (size_t ii=0; ii<mean_mom.size(); ii++) {
		if (chi2.at(ii)>chi2_min+1.) continue; //only care about 1-sigma region

		double mean_mom_i=mean_mom.at(ii);
		double sigma_mom_i=sigma_mom.at(ii);
		//if (mean_mom_i==mean_mom_minchi2) {
		if (sigma_mom_i==sigma_mom_minchi2) {
			//cout<<"index:"<<ii<<" mean_mom:"<<mean_mom_i<<endl;
			key_sel0.push_back(ii);
			mean_sel0.push_back(mean_mom.at(ii));
			sigma_sel0.push_back(sigma_mom.at(ii));
		}
	}

	//skip some points
	int cnt_jump=0;
	int n_skip=4;
	vector<int> key_sel;
	vector<double> mean_sel;
	vector<double> sigma_sel;
	for (size_t ii=0; ii<mean_sel0.size(); ii++) {
		cnt_jump++;
		if (cnt_jump>=n_skip) cnt_jump=0;
		//if (cnt_jump==n_skip-1) {
		key_sel.push_back(key_sel0.at(ii));
		mean_sel.push_back(mean_sel0.at(ii));
		sigma_sel.push_back(sigma_sel0.at(ii));

		//cout<<"key_list.push_back("<<key_sel0.at(ii)<<");"<<endl;
		//}
	}
	cout<<"ich habe "<<mean_sel.size()<<" Punkte..."<<endl;
	TGraph *gr_sel=new TGraph(mean_sel.size(), &mean_sel.at(0), &sigma_sel.at(0));
	gr_sel->SetMarkerStyle(8);
	gr_sel->SetMarkerSize(.5);
	gr_sel->SetMarkerColor(0);

	gr_sel->Draw("p same");


        //c0->Print(Form("%s/pbeam_data_mc.eps",outpath.Data()));
	c_1->Print(Form("%s/bmrw_trklen.eps",outpath.Data()));

	c_chi2_x->Print(Form("%s/bmrw_chi2x.eps",outpath.Data()));
	c_chi2_y->Print(Form("%s/bmrw_chi2y.eps",outpath.Data()));
	c_2->Print(Form("%s/mu_sigma_chi2.eps",outpath.Data()));
	c_ns->Print(Form("%s/mu_sigma_chi2_withcontour.eps",outpath.Data()));


	//save the config. & best-fit parameter --------------------------------------------------------------------//
	//min chi^2 index
        TParameter<Int_t>* key_minchi2=new TParameter<Int_t>("key_minchi2",0.); 
	key_minchi2->SetVal(min_index);

	//weighting func of minchi^2
	TF1 *gn_default;
        TF1 *gn_minchi2;
        TF1 *gng_minchi2;
	double mu_minchi2=0;
	double sigma_minchi2=0;
        double xmin=0.; //pmin [MeV/c]
        double xmax=2000.; //pmax [MeV/c]
	int key=0;
        for (int imu=0; imu<nmu; ++imu){ //mu loop
                double frac_mu=mu_st-(double)imu*dmu;
                double mu=m1*frac_mu;
                for (int isigma=0; isigma<nsigma; ++isigma){ //sigma loop
                        double frac_sigma=sigma_st-(double)isigma*dsigma;
                        double sigma=s1*frac_sigma;

			if (key==min_index) { //key of minchi2 	
				mu_minchi2=mu;
				sigma_minchi2=sigma;
                        	//Gaussian with changed mean and sigma
                        	gn_minchi2=new TF1(Form("gn_%d",key),fitg,xmin,xmax,3);
                        	gn_minchi2->SetParameter(0,mu);
                        	gn_minchi2->SetParameter(1,sigma);
                        	gn_minchi2->SetParameter(2,1);

                        	//weighting func. (beam mom)
                        	gng_minchi2=new TF1(Form("gng_%d",key),govg,xmin,xmax,4);
                                gng_minchi2->SetParameter(0,m1);
                        	gng_minchi2->SetParameter(1,s1);
                        	gng_minchi2->SetParameter(2,mu);
                        	gng_minchi2->SetParameter(3,sigma);
			} //key of minchi2

                        key++;
                 } //sigma loop
        } //mu loop
	gn_default=new TF1(Form("gn_default"),fitg,xmin,xmax,3);
	gn_default->SetParameter(0,m1);
        gn_default->SetParameter(1,s1);
        gn_default->SetParameter(2,1);

	gn_default->SetName("gaus_default");
	gn_minchi2->SetName("gaus_minchi2");
	gng_minchi2->SetName("bmrw_minchi2");

	TF1 *gn_data=new TF1(Form("gn_data"),fitg,xmin,xmax,3);
	gn_data->SetParameter(0,m1_data);
        gn_data->SetParameter(1,s1_data);
        gn_data->SetParameter(2,1);


	
	//plot figure -----------------------------------------------------------------------//
	TCanvas *c_g=new TCanvas("c_g",""); 
	c_g->Divide(1,1);
	c_g->cd(1);
	c_g->cd(1)->SetGridx();
	c_g->cd(1)->SetGridy();
	gn_default->SetLineColor(1);
	gn_minchi2->SetLineColor(2);
	gn_data->SetLineColor(4);

	TH2D* f2d=new TH2D("f2d","", 800, 600, 1400, 150, 0, 1.5);
	f2d->SetTitle(";Proton Momentum [MeV/c];");
	f2d->GetXaxis()->CenterTitle();
	f2d->Draw();
	gn_default->Draw("c same");
	gn_minchi2->Draw("c same");
	gn_data->Draw("c same");

	TLegend *leg0g = new TLegend(0.14,0.75,.75,0.88);
	leg0g->SetFillStyle(0);
	leg0g->AddEntry(gn_default, Form("MC Default - (#mu,#sigma):(%.2f,%.2f) MeV/c",m1,s1), "l");
	leg0g->AddEntry(gn_minchi2, Form("MC Weighted - (#mu,#sigma):(%.2f,%.2f) MeV/c",mu_minchi2,sigma_minchi2), "l");
	leg0g->AddEntry(gn_data, Form("Data - (#mu,#sigma):(%.2f,%.2f) MeV/c",m1_data,s1_data), "l");

	leg0g->Draw();

	//another pad for weighting func.
	TPad *overlay = new TPad("overlay","",0,0,1,1);
   	overlay->SetFillStyle(4000); //will be transparent
   	overlay->SetFrameFillStyle(0);
	overlay->Draw();
	overlay->cd();

        //Double_t pxmin = overlay->GetUxmin();
   	//Double_t pxmax = overlay->GetUxmax();
        Double_t pxmin = xmin;
   	Double_t pxmax = xmax;
	Double_t pymin = 0;
   	Double_t pymax = 40;
   	TH1F *hframe = overlay->DrawFrame(pxmin,pymin,pxmax,pymax);
	//hframe->GetXaxis()->SetLabelOffset(99);
	hframe->GetYaxis()->SetLabelOffset(99);
	hframe->GetYaxis()->SetLabelSize(0);
	hframe->GetXaxis()->SetLabelSize(0);
	hframe->GetXaxis()->SetTickLength(0.0); 
	hframe->GetYaxis()->SetTickLength(0.0); 
	
	gng_minchi2->SetLineColor(8);
	gng_minchi2->SetLineStyle(2);
   	gng_minchi2->Draw("same");

   	//Draw an axis on the right side
   	TGaxis *axis = new TGaxis(pxmax,pymin,pxmax, pymax,pymin,pymax,510,"+L");
   	axis->SetLineColor(8);
   	axis->SetLabelColor(8);
	axis->SetTitleColor(8); 
	axis->SetTitle("Weighting function");
   	axis->Draw();
	c_g->Print(Form("%s/gaus_beforeafter_bmrw.eps",outpath.Data()));


	//TFile *f_out = TFile::Open("bmrw_bestfit.root","recreate");
	//TFile *f_out = TFile::Open("bmrw_bestfit_calo.root","recreate");
	TFile *f_out = TFile::Open("bmrw_bestfit_calo_rmtrack.root","recreate");
		key_minchi2->Write();

		bm_nmu->Write();
		bm_dmu->Write();
		bm_mu_st->Write();

		bm_nsigma->Write();
		bm_dsigma->Write();
		bm_sigma_st->Write();

		gn_default->Write();
		gn_minchi2->Write();
		gng_minchi2->Write();
	f_out->Close();








}
