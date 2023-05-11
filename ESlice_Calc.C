#include "./headers/BasicParameters.h"
#include "./headers/BasicFunctions.h"
//#include "./headers/SliceParams.h"
//#include "./headers/ThinSlice.h"
#include "./headers/util.h"
#include "./headers/BetheBloch.h"
//#include "./headers/ESliceParams.h"

void ESlice_Calc() {

	//** SliceID Calc **************************************************************************************//
	//[1]Smaller bin at higher KE ----------------------------
        //const int nthinslices_small=91; //total # of slices
        const double thinslicewidth_small=20; //MeV
        const double Emin_small=420; //emin
        const double Emax_small=600; //emax
        const int nthinslices_small=1+(int)((Emax_small-Emin_small)/thinslicewidth_small); //total # of slices

	cout<<"nthinslices_small:"<<nthinslices_small<<endl;

	int bin_st=0;
	vector< pair <double,double> > bin_ke;
	for (int i=0; i<nthinslices_small; ++i) {
		int tmp_KE=Emax_small-i*thinslicewidth_small;

		bin_ke.push_back(std::make_pair(bin_st, tmp_KE));
		cout<<bin_st<<", "<<tmp_KE<<endl;
		bin_st++;
	}

	//[2]Larger bin at lower KE -------------------------------
        //const int nthinslices=21; //total # of slices
        const double thinslicewidth=20; //MeV
        const double Emin=0; //emin
        const double Emax=Emin_small; //emax
        const int nthinslices=(Emax-Emin)/thinslicewidth;
	cout<<"nthinslices:"<<nthinslices<<"\n\n"<<endl;

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
	cout<<"\n\n\n**Copy and paste the following strings ...***"<<endl;
	cout<<"const int reco_nbins="<<bin_ke.size()<<";"<<endl;
	cout<<"const int true_nbins="<<bin_ke.size()<<";"<<endl;
	
	cout<<reco_KE<<endl;
	cout<<true_KE<<endl;
	//cout<<""<<endl;
	cout<<reco_bins<<endl;
	cout<<true_bins<<endl;



}
