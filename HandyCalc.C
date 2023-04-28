#include "./headers/BasicParameters.h"
#include "./headers/BasicFunctions.h"
//#include "./headers/SliceParams.h"
//#include "./headers/ThinSlice.h"
#include "./headers/util.h"
#include "./headers/BetheBloch.h"
//#include "./headers/ESliceParams.h"

void HandyCalc(){

/*
double mom_beam=2; //unit:GeV/c
double ke_beam_MeV=1000.*p2ke(mom_beam); //ke_beam [MeV]
//double ke_beam_spec_MeV=1000.*ke_beam_spec; //ke_beam_spec [MeV]
double csda=csda_range_vs_mom_sm->Eval(mom_beam);

std::cout<<"mom_beam:"<<mom_beam<<" [GeV/C]"<<std::endl;
std::cout<<"ke_beam_MeV:"<<ke_beam_MeV<<" [MeV]"<<std::endl;
std::cout<<"csda:"<<csda<<" [cm]"<<std::endl;
*/

/*
//systematic uncertainties ----------------------------------------------------------------//
double ke=400; //MeV
double ke_0=1000.*p2ke(Pbeam_sys(ke2p(ke/1000.), 0.3/100., 0)); //central value
double ke_up=1000.*p2ke(Pbeam_sys(ke2p(ke/1000.), 0.3/100., 1)); //up value
double ke_dn=1000.*p2ke(Pbeam_sys(ke2p(ke/1000.), 0.3/100., -1)); //dn value
cout<<"ke:"<<ke<<endl;
cout<<"ke_0:"<<ke_0<<endl;
cout<<"ke_up:"<<ke_up<<endl;
cout<<"ke_dn:"<<ke_dn<<endl;
cout<<"dke_up:"<<ke_up-ke_0<<endl;
cout<<"dke_dn:"<<ke_0-ke_dn<<endl;
*/

/*
int ke=Emax;
for (int i=0;i<nthinslices;++i) {
	cout<<"["<<i<<"]: "<<ke<<endl;
	ke-=thinslicewidth;

}
*/

/*
const int nthinslices=40; //total # of slices
const double thinslicewidth=20; //MeV
const int name_thinslicewidth=20; //MeV
const double Emin=0; //emin
const double Emax=800; //emax

double KE_ff_reco=404.1;
int reco_st_sliceID=int(ceil((Emax-KE_ff_reco)/thinslicewidth));

double KEend_reco=50.3;
int reco_sliceID = int(floor((Emax-KEend_reco)/thinslicewidth));
std::cout<<"reco_st_sliceID="<<reco_st_sliceID<<std::endl;
std::cout<<"(Emax-KE_ff_reco)/thinslicewidth="<<(Emax-KE_ff_reco)/thinslicewidth<<std::endl;
std::cout<<"ceil((Emax-KE_ff_reco)/thinslicewidth)="<<ceil((Emax-KE_ff_reco)/thinslicewidth)<<std::endl;
std::cout<<"reco_sliceID="<<reco_sliceID<<std::endl;

int another_reco_st_sliceID;
for (int i=0; i<nthinslices; ++i) {
   double ke_cen=Emax-((double)i+0.5)*thinslicewidth;
   double ke_min=ke_cen-10;
   double ke_max=ke_cen+10;
   //cout<<"ID:"<<i<<"ke_min-ke_max:"<<ke_min<<" - "<<ke_max<<endl;

   if (KE_ff_reco > ke_min) { 
	another_reco_st_sliceID=i;	
	break;
   }	
	
}

int another_reco_sliceID;
for (int i=0; i<nthinslices; ++i) {
   double ke_cen=Emax-((double)i+0.5)*thinslicewidth;
   double ke_min=ke_cen-10;
   double ke_max=ke_cen+10;
   //cout<<"ID:"<<i<<"ke_min-ke_max:"<<ke_min<<" - "<<ke_max<<endl;

   if (KEend_reco > ke_min) { 
	another_reco_sliceID=i;	
	break;
   }	
	
}

std::cout<<"\n another_reco_st_sliceID="<<another_reco_st_sliceID<<std::endl;
std::cout<<"\n another_reco_sliceID="<<another_reco_sliceID<<std::endl;
*/

/*
		const int n_eloss=14; //14 energy slicing in total
		int de_eloss=50; //50 MeV Slice
		int eloss_min=0;
		vector<int> Eloss;
		Eloss.push_back(eloss_min);

		for (int ii=0; ii<n_eloss; ++ii) {
			cout<<eloss_min<<" "<<eloss_min+de_eloss<<endl;

			eloss_min+=de_eloss;
			Eloss.push_back(eloss_min);
		}

		int n=Eloss.size();
		cout<<"n="<<n<<endl;

		for (int ii=0; ii<n; ++ii) {
			cout<<Eloss.at(ii)<<endl;
		}
		for (int ii=0; ii<n-1; ++ii) {
			cout<<Eloss.at(ii)<<" "<<Eloss.at(ii+1)<<endl;
		}
		cout<<"\n"<<endl;


		const int n_kebeam_slice=14; //14 energy slicing in total
		int d_kebeam=50; //50 MeV Slice
		int kebeam_min=0;
		vector<int> KEbeam_slice;
		KEbeam_slice.push_back(kebeam_min);

		//per histogram
		int n_edept=60;
		double edept_min=-20;
		double edept_max=100;
		
		for (int ii=0; ii<n_kebeam_slice; ++ii) {
			kebeam_min+=d_kebeam;
			KEbeam_slice.push_back(kebeam_min);
		}

		for (int ii=0; ii<(int)KEbeam_slice.size()-1; ++ii) {
			cout<<KEbeam_slice.at(ii)<<" "<<KEbeam_slice.at(ii+1)<<endl;
		}

*/

	//Simulation on SliceID ---------------------------------------------------------//
	//smaller bin at higher KE ---------------------------------
        const int nthinslices_small=91; //total # of slices
        const double thinslicewidth_small=2; //MeV
        const double name_thinslicewidth_small=22; //MeV
        const double Emin_small=420; //emin
        const double Emax_small=600; //emax

	int bin_st=0;
	vector< pair <double,double> > bin_ke;
	for (int i=0; i<nthinslices_small; ++i) {
		int tmp_KE=Emax_small-i*thinslicewidth_small;

		bin_ke.push_back(std::make_pair(bin_st, tmp_KE));
		cout<<bin_st<<", "<<tmp_KE<<endl;
		bin_st++;
	}

	//larger bin ------------------------------------------------
        const int nthinslices=21; //total # of slices
        const double thinslicewidth=20; //MeV
        const double name_thinslicewidth=20; //MeV
        const double Emin=0; //emin
        const double Emax=Emin_small; //emax

	for (int i=0; i<nthinslices; ++i) {
		int tmp_KE=Emax-thinslicewidth-i*thinslicewidth;

		bin_ke.push_back(std::make_pair(bin_st, tmp_KE));
		cout<<bin_st<<", "<<tmp_KE<<endl;
		bin_st++;
	}

	TString reco_KE;
	TString true_KE;
	TString reco_bins;
	TString true_bins;
	
	reco_KE+="const double reco_KE[reco_nbins] = {";
	true_KE+="const double true_KE[true_nbins] = {";

	reco_bins+="const int reco_bins[reco_nbins+1] = {-1, ";
	true_bins+="const int true_bins[true_nbins+1] = {-1, ";

	for (int i=0; i<bin_ke.size(); ++i) {
		reco_KE+=bin_ke[i].second;
		true_KE+=bin_ke[i].second;

		reco_bins+=bin_ke[i].first;
		true_bins+=bin_ke[i].first;

		if (i<(int)bin_ke.size()-1) { 
			reco_KE+=", ";
			true_KE+=", ";

			reco_bins+=", ";
			true_bins+=", ";
		}
		if (i==(int)bin_ke.size()-1) { 
			reco_KE+="";
			true_KE+="";

			reco_bins+="";
			true_bins+="";
		}
	}
	reco_KE+="};";
	true_KE+="};";

	reco_bins+="};";
	true_bins+="};";

	TString reco_nbins;
	TString true_nbins;
	cout<<"const int reco_nbins="<<bin_ke.size()<<";"<<endl;
	cout<<"const int true_nbins="<<bin_ke.size()<<";"<<endl;
	
	cout<<reco_KE<<endl;
	cout<<true_KE<<endl;
	//cout<<""<<endl;
	cout<<reco_bins<<endl;
	cout<<true_bins<<endl;
	

	//for (int i=0; i<nthinslices; ++i) {
	





/*
	//smaller bin
	const int nthinslices_fine=90; //total # of slices
        const double thinslicewidth_fine=2; //MeV
        const double name_thinslicewidth_fine=2; //MeV
        const double Emin_fine=Emax; //emin
        const double Emax_fine=600; //emax

	int bin_st=0;
	//vector< pair <double,double> > bin_ke;
	for (int i=0; i<nthinslices+1; ++i) {
		bin_ke.push_back(std::make_pair(bin_st, int(Emax-i*thinslicewidth)));
		bin_st++;
	}






	vector<int> reco_KE;
	vector<int> reco_bins;
	TString str_E;
	TString str_BIN;	
	str_E+=Form("{");
	double mm1=1007.1482; //MC prod4a [spec]
	double ss1=60.703307; //MC prod4a [spec]
	double mu_min=mm1-3.*ss1;
	double mu_max=mm1+3.*ss1;
	str_BIN+=Form("{");
	for (int i=0; i<bin_ke.size(); ++i) {
		std::cout << bin_ke[i].first << ": " << bin_ke[i].second << std::endl;
		str_E+=Form("%.0f", bin_ke[i].second);
		str_BIN+=Form("%.0f", bin_ke[i].first);
		if (i!=bin_ke.size()-1) { 
			str_E+=Form(",");
			str_BIN+=Form(",");
		}
	}
	str_E+=Form("};");
	str_BIN+=Form("};");
	cout<<bin_ke.size()<<endl;

	std::cout<<str_BIN<<std::endl;
	std::cout<<str_E<<std::endl;


	//double KE_ff_reco=386.45;
	//double KEend_reco=384.084;
	double KE_ff_reco=445.998;
	double KEend_reco=426.7;
	int reco_st_sliceID=int(ceil((Emax-KE_ff_reco)/thinslicewidth));
	int reco_sliceID = int(floor((Emax-KEend_reco)/thinslicewidth));

	std::cout<<"reco_st_sliceID="<<reco_st_sliceID<<"=ceil("<<((Emax-KE_ff_reco)/thinslicewidth)<<")"<<std::endl;
	std::cout<<"reco_sliceID="<<reco_sliceID<<"=floor("<<((Emax-KEend_reco)/thinslicewidth)<<")"<<std::endl;

*/


/*
	const int N=41;
	const double bins[N+1] = {-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40};
	const double KE[N] = {800,780,760,740,720,700,680,660,640,620,600,580,560,540,520,500,480,460,440,420,400,380,360,340,320,300,280,260,240,220,200,180,160,140,120,100,80,60,40,20,0};

	//double KE_ff_reco=433.6; //MeV

	double KE_ff_scan=500.;
	double dKE_scan=0.1;

	vector<int> dSLICE;
	vector<int> KEff;
	vector<int> ex;
	for (int k=0; k<2000; ++k) {
		double KE_ff_reco=KE_ff_scan;
		KEff.push_back(KE_ff_reco);	

		//Proton's approach -----------------------------------------------//
		//int reco_st_sliceID=int(ceil((Emax-KE_ff_reco)/thinslicewidth)+0.5);
		//int reco_st_sliceID=int(ceil((Emax-KE_ff_reco)/thinslicewidth)+1);
		int reco_st_sliceID=int((Emax-KE_ff_reco)/thinslicewidth+0.5);

		//std::cout<<"[HY] reco_st_sliceID="<<reco_st_sliceID<<std::endl;
	
		//Pion's approach ------------------------------------------------------//
		int reco_ini_sliceID=-99;
        	for (reco_ini_sliceID=0; reco_ini_sliceID<N-2; ++reco_ini_sliceID) {
          		if (KE_ff_reco > KE[reco_ini_sliceID]) break;
        	}
		//std::cout<<"[Pion] reco_st_sliceID="<<reco_ini_sliceID<<std::endl;
		if ((reco_ini_sliceID-reco_st_sliceID)!=0) {
			std::cout<<"reco_ini_sliceID(pion)-reco_st_sliceID(proton)="<<(reco_ini_sliceID-reco_st_sliceID)<<"="<<reco_ini_sliceID<<"-"<<reco_st_sliceID<<" KE="<<KE_ff_reco<<std::endl;
		}

		dSLICE.push_back((double)(reco_ini_sliceID-reco_st_sliceID));
		ex.push_back(0);

		KE_ff_scan-=dKE_scan;
	}

	TCanvas *c_=new TCanvas(Form("c"),"",900, 600);
	TH2D *f2d=new TH2D("f2d","",300,250,550, 10,-5,5);
	f2d->SetTitle(" ;KE[MeV]; SliceID(Proton)-Slice(Pion) [a.u.]");
	f2d->Draw();
	TGraph *gr = new TGraph(KEff.size(), &KEff.at(0), &dSLICE.at(0));
	gr->SetMarkerSize(1);
	gr->SetMarkerStyle(20);
	gr->Draw("pl same");
	
	c_->Print("dslc_0.5.png");
*/


	double mm1=1007.1482; //MC prod4a [spec]
	double ss1=60.703307; //MC prod4a [spec]
	double mu_min=mm1-3.*ss1;
	double mu_max=mm1+3.*ss1;

	double ke_beam_min=1000*p2ke(mu_min/1000); //ke_beam_spec [GeV]
	double ke_beam_max=1000*p2ke(mu_max/1000); //ke_beam_spec [GeV]

	std::cout<<"ke_beam_min="<<ke_beam_min<<" - "<<ke_beam_max<<std::endl;




}
