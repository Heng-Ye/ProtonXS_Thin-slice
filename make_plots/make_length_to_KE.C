//#include "../headers/betheBloch.h"
#include "../headers/BasicParameters.h"
#include "../headers/BasicFunctions.h"
#include "../headers/BetheBloch.h"
#include "./len2ke.h"

double mass_proton = 938.272046; //proton [unit: MeV/c^2]
double P_in_proton = 1000.; //MeV
//double KE_in_proton=sqrt(pow(P_in_proton,2) + pow(mass_proton,2) ) - mass_proton;
//double KE_in_proton=412.65;
double KE_in_proton=399.365;

/*
void hist_bethe_mean_distance(double E_init, double mass_particle, TH1D* h_bethe ){

	for(int i=1; i <= h_bethe->GetNbinsX(); i++){
		h_bethe->SetBinContent( i, betheBloch(E_init, mass_particle));
		h_bethe->SetBinError(i, 0.001 );
		E_init = E_init - betheBloch(E_init, mass_particle);
		if(E_init <= 0) return;

	};
};


void hist_NIST(double E_init, TH1D* h_bethe){
	for(int i=1; i <= h_bethe->GetNbinsX(); i++){
		h_bethe->SetBinContent( i, dEdx_vs_KE_sm->Eval(E_init));
		h_bethe->SetBinError(i, 0.001 );
		E_init = E_init - dEdx_vs_KE_sm->Eval(E_init);
		if(E_init <= 0) return;

	};
};
*/


double Len2KE(double len, double ke_ini) {

	float len_min=0;
	float len_max=300;
	int n_len=300;

	//create the map to convert trklen to Edept
	TH1D* dEdx = new TH1D("dEdx", "", n_len, len_min, len_max);
	//hist_NIST(ke_ini, dEdx); //loading in dE/dx map based on KE_int
	hist_bethe_mean_distance(ke_ini, mass_proton, dEdx); //loading in dE/dx map based on KE_int
	TH1* cumulative = dEdx->GetCumulative();

	int bin_cen=0;
	int bin_b=0;
	int bin_a=0;
	//int n_len=cumulative->GetNbinsX();
	
	std::cout<<"\n\nlen="<<len<<std::endl;
	std::cout<<"n_len="<<n_len<<std::endl;

	//get Edept
	bin_cen=cumulative->GetXaxis()->FindBin(len);
	bin_b=bin_cen-1;
	bin_a=bin_cen+1;
	if (bin_a>n_len) bin_a=n_len;
	if (bin_b<0) bin_b=0;

	bin_cen=cumulative->GetXaxis()->FindBin(len);
	bin_b=bin_cen-1;
	bin_a=bin_cen+1;
	if (bin_a>n_len) bin_a=n_len;
	if (bin_b<0) bin_b=0;

        double dept_cen=cumulative->GetBinContent(bin_cen);
	double dept_b=cumulative->GetBinContent(bin_b);
	double dept_a=cumulative->GetBinContent(bin_a);
	
	double m_dept=(dept_a-dept_b)/(cumulative->GetBinCenter(bin_a)-cumulative->GetBinCenter(bin_b));
	double b_dept=dept_a-m_dept*cumulative->GetBinCenter(bin_a);
	double ke_len=b_dept+m_dept*len; 

	if (ke_len>ke_ini) {
	  std::cout<<"NOT Possible conversion!"<<std::endl;
	  ke_len=ke_ini;
	}

	cout<<"\n\nTest!:: dept_cen="<<dept_cen<<" MeV"<<""<<endl;
	cout<<"            ke_len="<<ke_len<<" MeV"<<""<<endl;
	cout<<"            dept_b="<<dept_b<<" MeV"<<""<<endl;
	cout<<"            dept_a="<<dept_a<<" MeV"<<"\n\n"<<endl;

	return ke_len;

}


void make_length_to_KE() {

	//plot style --------------------------------------------------------------------//
	gROOT->LoadMacro(" ~/protoDUNEStyle.C"); //load pDUNE style
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	gStyle->SetOptStat(0);

	float len_min=0;
	float len_max=100;
	int n_len=100;

	TH1D* dEdx3 = new TH1D("dEdx3", "", n_len, len_min, len_max);
	hist_bethe_mean_distance( KE_in_proton, mass_proton, dEdx3);

	//TH1D* dEdx2 = new TH1D("dEdx2", "", n_len, len_min, len_max);
	//hist_NIST(KE_in_proton, dEdx2); //loading in dE/dx map based on KE_int

	TH1* cumulative = dEdx3->GetCumulative();
	cumulative->SetLineColor(3);
	cumulative->SetLineStyle(2);

	//TH1* cumulative2 = dEdx2->GetCumulative();
	//cumulative2->Smooth();

	//Edept calc -------------------------//
	int bin_cen=0;
	int bin_b=0;
	int bin_a=0;

	double len=50.1; //cm
	double trueKE_len=Len2KE(len, KE_in_proton);
	cout<<"trueKE_len1:"<<trueKE_len<<std::endl;

	LEN2E Len2KE_;
        Len2KE_.setmap(KE_in_proton);
        double trueKE_len2=Len2KE_.E(len);
	cout<<"trueKE_len2:"<<trueKE_len2<<std::endl;
	cout<<"len3=20.; trueKE_len3:"<<Len2KE_.E(20.)<<std::endl;
	cout<<"len4=10.; trueKE_len4:"<<Len2KE_.E(10.)<<std::endl;
	cout<<"len5=100.; trueKE_len5:"<<Len2KE_.E(100.)<<std::endl;
	cout<<"len5=-10.; trueKE_len6:"<<Len2KE_.E(-10.)<<std::endl;
    


	//test TJ's class///
	BetheBloch BB;
	BB.SetPdgCode(2212);
	//BB.meandEdx(KE_in_proton); //input: KE_ff
	//BB.KEAtLength(KE_in_proton, 20.);
	
	double ke_len3=KE_in_proton-BB.KEAtLength(KE_in_proton, 20.);
 	cout<<"\n len3=20; KE="<<ke_len3<<endl;
	cout<<"meandEdx:"<<BB.meandEdx(KE_in_proton)<<endl;
	cout<<"RangeFromKE:"<<BB.RangeFromKE(KE_in_proton)<<endl;
	cout<<"RangeFromKESpline:"<<BB.RangeFromKESpline(KE_in_proton)<<endl;
	//cout<<"KEFromRangeSpline:"<<KEFromRangeSpline(



/*
	bin_cen=cumulative2->GetXaxis()->FindBin(len);
	bin_b=bin_cen-1;
	bin_a=bin_cen+1;
	if (bin_a>n_len) bin_a=n_len;
	if (bin_b<0) bin_b=0;
        double dept_cen=cumulative2->GetBinContent(bin_cen);
	double dept_b=cumulative2->GetBinContent(bin_b);
	double dept_a=cumulative2->GetBinContent(bin_a);
	
	double m_dept=(dept_a-dept_b)/(cumulative2->GetBinCenter(bin_a)-cumulative2->GetBinCenter(bin_b));
	double b_dept=dept_a-m_dept*cumulative2->GetBinCenter(bin_a);
	double ke_len=b_dept+m_dept*len; 
*/


	//cumulative2->GetXaxis()->SetTitle("True Track Length [cm]"); cumulative->GetYaxis()->SetTitle("True Deposited Energy [MeV]");   
	cumulative->SetTitle(Form("Proton Kinetic Energy:%.2f MeV ; True Track Length [cm]; True Deposited Energy [MeV]", KE_in_proton));   
	cumulative->GetXaxis()->CenterTitle();
	cumulative->SetLineColor(4);
	cumulative->SetLineWidth(3);

	TCanvas *c = new TCanvas("c", "", 600,400);
	cumulative->Draw("");
	//cumulative->Draw("same");

	c->Update();

	//scale to pad coordinates
	Float_t rightmax = 1.1*dEdx3->GetMaximum();
	Float_t scale = gPad->GetUymax()/rightmax;
	dEdx3->SetLineColor(kRed);
	dEdx3->SetMarkerColor(kRed);
	dEdx3->SetMarkerSize(0.4);
	dEdx3->Scale(scale);
	//dEdx->GetYaxis()->SetRangeUser(1.5,4);
	dEdx3->Draw("C SAME");


	TGaxis *axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(), 0, rightmax, 510, "+L");
	axis->SetLineColor(kRed);
	axis->SetLabelColor(kRed);
	axis->SetTitleColor(kRed);
	axis->SetLabelSize(0.03);
	axis->SetTitle("dEdx [MeV / cm]");
	axis->Draw();

	//TLegend *leg = new TLegend(0.76, 0.45, 0.86, 0.55);
	TLegend *leg = new TLegend(0.155518,0.638298,0.471572,0.875);
	leg->SetFillColor(0);
	leg->SetTextSize(0.036);
	//leg->SetHeader("1GeV Momentum Proton", "C");
	//leg->AddEntry(dEdx, "dE/dx (NIST Data Base)", "L");
	//leg->AddEntry(cumulative2, "Deposited energy (NIST)", "L");
	//leg->AddEntry(cumulative, "Deposited energy (Bethe-Bloch formula)", "L");


	TLegendEntry* l1[3];
	//l1[0]=leg->AddEntry(dEdx3, "dE/dx (NIST Data Base)", "L"); l1[0]->SetTextColor(1);
	l1[0]=leg->AddEntry(dEdx3, "dE/dx (Bethe-Bloch)", "L"); l1[0]->SetTextColor(1);
	//l1[1]=leg->AddEntry(cumulative2, "Deposited energy (NIST)", "L"); l1[1]->SetTextColor(1);
	l1[2]=leg->AddEntry(cumulative, "Deposited energy (Bethe-Bloch formula)", "L"); l1[2]->SetTextColor(1);


	leg->Draw();


	c->Print(Form("Edept_len_KE%.2fMeV.eps",KE_in_proton));


}
