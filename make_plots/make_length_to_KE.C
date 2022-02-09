#include "../headers/betheBloch.h"
#include "../headers/BasicParameters.h"
#include "../headers/BasicFunctions.h"

double mass_proton = 938.272046; //proton [unit: MeV/c^2]
double P_in_proton = 1000.; //MeV
//double KE_in_proton=sqrt(pow(P_in_proton,2) + pow(mass_proton,2) ) - mass_proton;
double KE_in_proton=700;

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


double Len2KE(double len, double ke_ini) {

	float len_min=0;
	float len_max=300;
	int n_len=300;

	//create the map to convert trklen to Edept
	TH1D* dEdx = new TH1D("dEdx", "", n_len, len_min, len_max);
	hist_NIST(ke_ini, dEdx); //loading in dE/dx map based on KE_int
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
	float len_max=300;
	int n_len=300;

	TH1D* dEdx = new TH1D("dEdx", "", n_len, len_min, len_max);
	hist_bethe_mean_distance( KE_in_proton, mass_proton, dEdx);

	TH1D* dEdx2 = new TH1D("dEdx2", "", n_len, len_min, len_max);
	hist_NIST(KE_in_proton, dEdx2); //loading in dE/dx map based on KE_int

	TH1* cumulative = dEdx->GetCumulative();
	cumulative->SetLineColor(3);
	cumulative->SetLineStyle(2);

	TH1* cumulative2 = dEdx2->GetCumulative();
	//cumulative2->Smooth();

	//Edept calc -------------------------//
	int bin_cen=0;
	int bin_b=0;
	int bin_a=0;

	double len=150.1; //cm
	double trueKE_len=Len2KE(len, KE_in_proton);
	cout<<"trueKE_len:"<<trueKE_len<<std::endl;

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
	cumulative2->SetTitle(Form("Proton Kinetic Energy:%.2f MeV ; True Track Length [cm]; True Deposited Energy [MeV]", KE_in_proton));   
	cumulative2->GetXaxis()->CenterTitle();
	cumulative2->SetLineColor(4);
	cumulative2->SetLineWidth(3);

	TCanvas *c = new TCanvas("c", "", 600,400);
	cumulative2->Draw("");
	cumulative->Draw("same");

	c->Update();

	//scale to pad coordinates
	Float_t rightmax = 1.1*dEdx->GetMaximum();
	Float_t scale = gPad->GetUymax()/rightmax;
	dEdx2->SetLineColor(kRed);
	dEdx2->SetMarkerColor(kRed);
	dEdx2->SetMarkerSize(0.4);
	dEdx2->Scale(scale);
	//dEdx->GetYaxis()->SetRangeUser(1.5,4);
	dEdx2->Draw("C SAME");


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
	l1[0]=leg->AddEntry(dEdx2, "dE/dx (NIST Data Base)", "L"); l1[0]->SetTextColor(1);
	l1[1]=leg->AddEntry(cumulative2, "Deposited energy (NIST)", "L"); l1[1]->SetTextColor(1);
	l1[2]=leg->AddEntry(cumulative, "Deposited energy (Bethe-Bloch formula)", "L"); l1[2]->SetTextColor(1);


	leg->Draw();


	c->Print(Form("Edept_len_KE%.2fMeV.eps",KE_in_proton));


}
