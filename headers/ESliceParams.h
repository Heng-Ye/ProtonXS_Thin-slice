
	//const int nthinslices=25; //total # of slices
	//const double thinslicewidth=20; //MeV
	//const int name_thinslicewidth=20; //MeV
	//const double Emin=0; //emin
	//const double Emax=500; //emax

	//const int nthinslices=30; //total # of slices
	//const double thinslicewidth=20; //MeV
	//const int name_thinslicewidth=20; //MeV
	//const double Emin=0; //emin
	//const double Emax=600; //emax

	//const int nthinslices=600; //total # of slices
	//const double thinslicewidth=1; //MeV
	//const int name_thinslicewidth=1; //MeV
	//const double Emin=0; //emin
	//const double Emax=600; //emax

namespace p{

	const Int_t reco_nbins=31;
	const Int_t true_nbins=31;
	const Double_t reco_KE[reco_nbins] = {600, 580, 560, 540, 520, 500, 480, 460, 440, 420, 400, 380, 360, 340, 320, 300, 280, 260, 240, 220, 200, 180, 160, 140, 120, 100, 80, 60, 40, 20, 0};
	const Double_t true_KE[true_nbins] = {600, 580, 560, 540, 520, 500, 480, 460, 440, 420, 400, 380, 360, 340, 320, 300, 280, 260, 240, 220, 200, 180, 160, 140, 120, 100, 80, 60, 40, 20, 0};
	const Double_t reco_bins[reco_nbins+1] = {-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
	const Double_t true_bins[true_nbins+1] = {-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
}

        //const int nthinslices=120; //total # of slices
        //const double thinslicewidth=5; //MeV
        //const int name_thinslicewidth=5; //MeV
        //const double Emin=0; //emin
        //const double Emax=600; //emax

	//const int nthinslices=60; //total # of slices
	//const double thinslicewidth=10; //MeV
	//const int name_thinslicewidth=10; //MeV
	//const double Emin=0; //emin
	//const double Emax=600; //emax

	//const int nthinslices=600; //total # of slices
	//const double thinslicewidth=1; //MeV
	//const int name_thinslicewidth=1; //MeV
	//const double Emin=0; //emin
	//const double Emax=600; //emax

        //const int nthinslices=40; //total # of slices
        //const double thinslicewidth=20; //MeV
        //const int name_thinslicewidth=20; //MeV
        //const double Emin=0; //emin
        //const double Emax=800; //emax

	//Best config
        //const int nthinslices=30; //total # of slices
        //const double thinslicewidth=20; //MeV
        //const int name_thinslicewidth=20; //MeV
        //const double Emin=0; //emin
        //const double Emax=600; //emax

        //const int nthinslices=24; //total # of slices
        //const double thinslicewidth=25; //MeV
        //const int name_thinslicewidth=25; //MeV
        //const double Emin=0; //emin
        //const double Emax=600; //emax

        //const int nthinslices=12; //total # of slices
        //const double thinslicewidth=50; //MeV
        //const int name_thinslicewidth=50; //MeV
        //const double Emin=0; //emin
        //const double Emax=600; //emax


        //const int nthinslices=10; //total # of slices
        //const double thinslicewidth=80; //MeV
        //const int name_thinslicewidth=80; //MeV
        //const double Emin=0; //emin
        //const double Emax=800; //emax

        //const int nthinslices=20; //total # of slices
        //const double thinslicewidth=30; //MeV
        //const int name_thinslicewidth=30; //MeV
        //const double Emin=0; //emin
        //const double Emax=600; //emax

        //const int nthinslices=40; //total # of slices
        //const double thinslicewidth=15; //MeV
        //const int name_thinslicewidth=15; //MeV
        //const double Emin=0; //emin
        //const double Emax=600; //emax


	Int_t nbinse=120; //nbins for KE 
	//Int_t nbinse=12; //nbins for KE 
	Int_t nbinsthickness = 1200; //nbins for trk pitch

