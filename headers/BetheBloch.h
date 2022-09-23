#include "TF1.h"
#include "TSpline.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <map>

#include "Math/VavilovAccurate.h"

ROOT::Math::VavilovAccurate vav;

using namespace std;

class TSpline3;

class BetheBloch {

 public:

  BetheBloch();
  BetheBloch(int pdg);

  void SetPdgCode(int pdg);

  int GetPdgCode(){ return pdgcode;};

  double meandEdx(double KE);

  double MPVdEdx(double KE, double pitch);

  double IntegratedEdx(double KE0, double KE1, int n = 10000);

  double RangeFromKE(double KE);

  double RangeFromKESpline(double KE);

  double KEFromRangeSpline(double range);

  double Best_Chi2=-1;

  double KEAtLength(double KE0, double tracklength);

  void CreateSplineAtKE(int iKE);

  double Fit_dEdx_Residual_Length(const vector<double> dEdx, const vector<double> ResRange, int PID, bool save_graph); //Sungbin's dE/dx vs RR Fitter  

  double Landau_xi(double KE, double pitch);

  double Get_Wmax(double KE);
 
  double dEdx_PDF(double KE, double pitch, double dEdx); //Sungbin's dEdx_PDF

  double Fit_Proton_Residual_Length_Likelihood(const vector<double> dEdx, const vector<double> ResRange, int PID, bool save_graph); //Sungbin's dE/dx vs RR Fitter using likelihood func.
 
 private:

  int pdgcode;
  double mass;
  int charge;

  TSpline3 *sp_KE_range;
  TSpline3 *sp_range_KE;

  map<int, TSpline3*> spmap;
  
  double densityEffect(double beta, double gamma);

  double betaGamma(double KE);

  void CreateSplines(int np = 1000, double minke = .01, double maxke = 2e5);

};


BetheBloch::BetheBloch()
  : pdgcode(0)
  , mass(0)
  , charge(0)
  , sp_KE_range(0)
  , sp_range_KE(0){
}

BetheBloch::BetheBloch(int pdg)
  : pdgcode(0)
  , mass(0)
  , charge(0)
  , sp_KE_range(0)
  , sp_range_KE(0){
  SetPdgCode(pdg);
}

void BetheBloch::SetPdgCode(int pdg){

  pdgcode = pdg;

  if (abs(pdgcode) == 13){//muon
    mass = 105.6583755;
    charge = 1;
  }
  else if (abs(pdgcode) == 211){//pion
    mass = 139.57039;
    charge = 1;
  }
  else if (abs(pdgcode) == 321){//kaon
    mass = 493.677;
    charge = 1;
  }
  else if (pdgcode == 2212){//proton
    mass = 938.27208816;
    charge = 1;
  }
  else{
    cout<<"Unknown pdg code "<<pdgcode<<endl;
    exit(1);
  }

  CreateSplines();

}

double BetheBloch::densityEffect(double beta, double gamma){

  double lar_C = 5.215, lar_x0 = 0.201, lar_x1 = 3, lar_a = 0.196, lar_k = 3;
  double x = log10(beta * gamma);
  
  if( x >= lar_x1 ){
    return 2*log(10)*x - lar_C;
  }
  else if ( lar_x0 <= x && x < lar_x1){
    return 2*log(10)*x - lar_C + lar_a * pow(( lar_x1 - x ) , lar_k );
  }
  else{
    return 0; //if x < lar_x0
  }
}

double BetheBloch::betaGamma(double KE){

  double gamma, beta;
  gamma = (KE + mass) / mass;
  beta = sqrt( 1 - 1/pow(gamma,2));
   
  return beta*gamma;

}

double BetheBloch::meandEdx(double KE){

  //KE is kinetic energy in MeV
  
  double K = 0.307;
  double rho = 1.396;
  double Z = 18;
  double A = 39.948;
  double I = pow(10,-6)*10.5*18; //MeV
  double me = 0.511; //MeV me*c^2
  
  double gamma = (KE + mass) / mass;
  double beta = sqrt( 1 - 1/pow(gamma,2));
  
  double wmax = 2*me*pow(beta,2)*pow(gamma,2)/(1+2*gamma*me/mass + pow(me,2)/pow(mass,2));

  //cout<<wmax<<" "<<rho*K*Z*pow(charge,2)<<" "<<A*pow(beta,2)<<" "<<0.5*log(2*me*pow(gamma,2)*pow(beta,2)*wmax/pow(I,2))<<" "<<pow(beta,2)<<" "<<densityEffect( beta, gamma )/2<<endl;
  double dEdX = (rho*K*Z*pow(charge,2))/(A*pow(beta,2))*(0.5*log(2*me*pow(gamma,2)*pow(beta,2)*wmax/pow(I,2)) - pow(beta,2) - densityEffect( beta, gamma )/2 );

  return dEdX;
}

double BetheBloch::MPVdEdx(double KE, double pitch){

  //KE is kinetic energy in MeV
  //pitch is in cm

  double K = 0.307;
  double rho = 1.396;
  double Z = 18;
  double A = 39.948;
  double I = pow(10,-6)*10.5*18; //MeV
  double  me = 0.511; //MeV me*c^2

  double gamma = (KE + mass) / mass;
  double beta = sqrt( 1 - 1/pow(gamma,2));

  double xi = ( K/2 )*( Z/A )* ( pitch * rho / pow(beta,2));
  
  double eloss_mpv = xi*(log( 2*me*pow(gamma,2)*pow(beta,2) / I ) + log( xi / I ) + 0.2 - pow(beta,2) - densityEffect( beta, gamma ) )/pitch;

  return eloss_mpv;
}

double BetheBloch::IntegratedEdx(double KE0, double KE1, int n){

  if (KE0>KE1) swap(KE0, KE1);

  double step = (KE1-KE0)/n;

  double area = 0;
  
  for (int i = 0; i<n; ++i){
    double dEdx = meandEdx(KE0 + (i+0.5)*step);
    if (dEdx)
      area += 1/dEdx*step;
  }
  return area;
}

double BetheBloch::RangeFromKE(double KE){

  return IntegratedEdx(0, KE);
}

void BetheBloch::CreateSplines(int np, double minke, double maxke){

  if (sp_KE_range) delete sp_KE_range;
  if (sp_range_KE) delete sp_range_KE;

  for (const auto & x : spmap){
    if (x.second) delete x.second;
  }
  spmap.clear();
  
  double *KE = new double[np];
  double *Range = new double[np];
//  vector<double> KE(np);
//  vector<double> Range(np);

  for (int i = 0; i<np; ++i){
    double ke = pow(10, log10(minke)+i*log10(maxke/minke)/np);
//    KE.push_back(ke);
//    Range.push_back(RangeFromKE(ke));
    KE[i] = ke;
    Range[i] = RangeFromKE(ke);
    //if (pdgcode ==2212) cout<<KE[i]<<" "<<Range[i]<<" "<<meandEdx(KE[i])<<endl;
  }

  sp_KE_range = new TSpline3("sp_KE_range", KE, Range, np, "b2e2", 0, 0);
  sp_range_KE = new TSpline3("sp_range_KE", Range, KE, np, "b2e2", 0, 0);
  //cout<<sp_KE_range->Eval(10)<<endl;
  delete[] KE;
  delete[] Range;
  cout<<"Done creating splines for particle with pdgcode "<<pdgcode<<endl;
}

double BetheBloch::RangeFromKESpline(double KE){
  if (!sp_KE_range){
    cout<<"Spline does not exist."<<endl;
    exit(1);
  }
  return sp_KE_range->Eval(KE);
}

double BetheBloch::KEFromRangeSpline(double range){
  if (!sp_range_KE){
    cout<<"Spline does not exit."<<endl;
    exit(1);
  }
  return sp_range_KE->Eval(range);
}

double BetheBloch::KEAtLength(double KE0, double tracklength){
  double tmp_ke_at_len=0;
  int iKE = int(KE0);

  if (spmap.find(iKE)==spmap.end()){
    CreateSplineAtKE(iKE);
  }

  double deltaE = spmap[iKE]->Eval(tracklength);

  tmp_ke_at_len=KE0 - deltaE;
  if (deltaE < 0) { 
	//cout<<"Negative delta E: "<<deltaE<<endl; //meaning in the region of KE<500 keV
	tmp_ke_at_len=-700;
  }	
  //if (tmp_ke_at_len < 0) { 
        //tmp_ke_at_len=-10;
	//cout<<"Negative KE: "<<KE0 - deltaE<<" --> set KE to zero"<<endl;
  //}
  
  return tmp_ke_at_len;

}

void BetheBloch::CreateSplineAtKE(int iKE){

  double KE0 = iKE;

  // Sample every 10 MeV
  int np = int(KE0/10);
  double *deltaE;
  double *trklength;
  if (np>1){
    deltaE = new double[np];
    trklength = new double[np];
    for (int i = 0; i<np; ++i){
      double KE = KE0 - i*10;
      deltaE[i] = KE0 - KE;
      trklength[i] = IntegratedEdx(KE, KE0);
      //cout<<iKE<<" "<<KE<<" "<<i<<"/"<<np<<" "<<trklength[i]<<endl;
    }
  }
  else{
    //cout<<"KE too low: "<<iKE<<endl;
    np = 2;
    deltaE = new double[np];
    trklength = new double[np];
    deltaE[0] = 0;
    trklength[0] = 0;
    deltaE[1] = KE0;
    trklength[1] = RangeFromKE(KE0);
  }

  spmap[iKE] = new TSpline3(Form("KE %d",iKE), trklength, deltaE, np, "b2e2", 0, 0);
  delete[] trklength;
  delete[] deltaE;

}


double BetheBloch::Fit_dEdx_Residual_Length(const vector<double> dEdx, const vector<double> ResRange, int PID, bool save_graph) {

  bool this_is_beam = true;
  int N_max = 200; // == Maximum number of hits used for the Bethe-Bloch fitting

  // == PID input : mass hypothesis, valid only for muons, charged pions, and protons
  int abs_PID = abs(PID);
  if(!(abs(PID) == 13 || PID == 2212 || abs(PID) == 211)){
    //cout << "[HadAna::Fit_dEdx_Residual_Length] Not a valid PID!" << endl;
    return -9999.;
  }

  double best_additional_res_length = -0.1;
  double best_chi2 = 99999.;
  double min_additional_res_length = 0.; // == [cm]
  double max_additional_res_length = 200.; // == [cm]
  //double res_length_step = 0.5; // == [cm] //HY::Spicky structure appears if set to 0.5 cm step size
  double res_length_step = 0.105; // == [cm //HY:Set step size to 0.105 cm
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;
  int this_N_calo = dEdx.size();
  if(this_N_calo <= 15){
    //cout << "[HadAna::Fit_dEdx_Residual_Length] Too small number of hits!" << endl;
    return -9999.; // == Too small number of hits
  }
  int this_N_hits = TMath::Min(this_N_calo, N_max); // == Use how many hits
  int i_bestfit = -1;
  vector<double> chi2_vector;
  vector<double> additional_res_legnth_vector;
  for(int i = 0; i < res_length_trial; i++){
    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_chi2 = 0.;
    for(int j = 5; j < this_N_hits - 5; j++){ // == Do not use first and last 5 hits
      int this_index = this_N_calo - 1 - j;
      if(this_is_beam) this_index = j;
      double this_res_length = ResRange.at(this_index) - ResRange.at(this_N_calo - this_N_hits) + this_additional_res_length;
      if(this_is_beam) this_res_length = ResRange.at(this_index) - ResRange.at(this_N_hits - 1) + this_additional_res_length;
      
      //cout << Form("[HadAna::Fit_dEdx_Residual_Length] ResRange.at(%d) = %f, dEdx.at(%d) = %f", j, ResRange.at(j), j, dEdx.at(j)) << endl;

      //double this_KE = ResLength_to_KE_BB(this_res_length, this_mass);
      //double this_KE = map_BB[abs_PID]->KEFromRangeSpline(this_res_length);
      double this_KE = KEFromRangeSpline(this_res_length);
      //double dEdx_theory = dEdx_Bethe_Bloch(this_KE, this_mass);
      //double dEdx_theory = map_BB[abs_PID]->meandEdx(this_KE);
      double dEdx_theory = meandEdx(this_KE);
      double dEdx_measured = dEdx.at(this_index);
      if(dEdx_measured < 0.5 || dEdx_measured > 20.0) continue; // == Truncate, it should be modified to consider protons

      // == Gaussian approx.
      //double dEdx_theory_err = dEdx_theory * 0.02;
      this_chi2 += pow(dEdx_measured - dEdx_theory, 2);
    }
    this_chi2 = this_chi2 / (this_N_hits + 0.); // == chi2 / n.d.f
    if(this_chi2 < best_chi2){
      best_chi2 = this_chi2;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }
    Best_Chi2=best_chi2;

    if(save_graph){
      // == Save vectors for graphes
      chi2_vector.push_back(this_chi2);
      additional_res_legnth_vector.push_back(this_additional_res_length);
    }
  }

  if(save_graph){
/*
    // == Vectors for graphes
    vector<double> range_original;
    vector<double> range_bestfit;
    vector<double> range_reco;
    vector<double> dEdx_ordered;
    for(int i = 5; i < this_N_hits - 5; i++){
      int this_index = this_N_calo - 1 - i;
      if(this_is_beam) this_index = i;
      double this_range_original = ResRange.at(this_index) - ResRange.at(this_N_calo - this_N_hits);
      if(this_is_beam) this_range_original = ResRange.at(this_index) - ResRange.at(this_N_hits - 1);
      range_original.push_back(this_range_original);

      double this_range_bestfit = ResRange.at(this_index) - ResRange.at(this_N_calo - this_N_hits) + best_additional_res_length;
      if(this_is_beam) this_range_bestfit = ResRange.at(this_index) - ResRange.at(this_N_hits - 1) + best_additional_res_length;
      range_bestfit.push_back(this_range_bestfit);
      range_reco.push_back(ResRange.at(this_index));
      dEdx_ordered.push_back(dEdx.at(this_index));
    }
    TGraph *dEdx_gr = new TGraph(this_N_hits - 10, &range_original[0], &dEdx_ordered[0]);
    dEdx_gr -> SetName(Form("dEdx_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits - 1));
    dEdx_gr -> Write();
    delete dEdx_gr;

    TGraph *dEdx_bestfit_gr = new TGraph(this_N_hits - 10,&range_bestfit[0], &dEdx_ordered[0]);
    dEdx_bestfit_gr -> SetName(Form("dEdx_bestfit_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    dEdx_bestfit_gr -> Write();
    delete dEdx_bestfit_gr;

    TGraph *dEdx_reco_gr = new TGraph(this_N_hits - 10,&range_reco[0], &dEdx_ordered[0]);
    dEdx_reco_gr -> SetName(Form("dEdx_reco_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    dEdx_reco_gr -> Write();
    delete dEdx_reco_gr;

    TGraph *chi2_gr = new TGraph(additional_res_legnth_vector.size(), &additional_res_legnth_vector[0], &chi2_vector[0]);
    chi2_gr -> SetName(Form("Chi2_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    chi2_gr -> Write();
    chi2_vector.clear();
    additional_res_legnth_vector.clear();
    delete chi2_gr;
*/
  }

  double original_res_length = ResRange.at(this_N_calo - 1) - ResRange.at(this_N_calo - this_N_hits); // == [cm]
  if(this_is_beam) original_res_length = ResRange.at(0) - ResRange.at(this_N_hits - 1); // == [cm]
  double best_total_res_length = best_additional_res_length + original_res_length;
  //double best_KE = map_BB[abs_PID]->KEFromRangeSpline(best_total_res_length);
  double best_KE = KEFromRangeSpline(best_total_res_length);
  //double best_mom = map_BB[abs_PID]->KEtoMomentum(best_KE);
  //double best_mom = KEtoMomentum(best_KE);

  // == Define fitting failed cases
  if(i_bestfit == res_length_trial - 1){
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : no mimumum" << endl;
    return -9999.;
  }
  else if(best_chi2 > 99990.){
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : best_chi2 > 99990." << endl;
    return -9999.;
  }
  else if(best_chi2 < 1.0e-11){
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : best_chi2 < 1.0e-11" << endl;
    return -9999.;
  }

  //return best_mom;
  return best_total_res_length;


}


double BetheBloch::Landau_xi(double KE, double pitch){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double xi = rho * pitch * 0.5 * K * (Z / A) * pow(1. / beta, 2);
  return xi;
}

double BetheBloch::Get_Wmax(double KE){
  double me = 0.511; //MeV me*c^2
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * me * pow(beta * gamma, 2)) / (1.0 + 2.0 * me * (gamma / mass) + pow((me / mass),2));

  return Wmax;
}



double dEdx_PDF_fuction(double *x, double *par){
  // == par[5] = {kappa, beta^2, xi, <dE/dx>BB, width}
  double a = par[2] / par[4];
  double b = (0.422784 + par[1] + log(par[0])) * par[2] / par[4] + par[3];
  double y = (x[0] - b) / a;

  double this_vav = 0.;

  if(par[0] < 0.01){ // == Landau
    this_vav = TMath::Landau(y);
    this_vav =  this_vav / a;
  }
  else if(par[0] > 10.){ // == Gaussian
    double mu = vav.Mean(par[0], par[1]);
    double sigma = sqrt(vav.Variance(par[0], par[1]));
    this_vav =  TMath::Gaus(y, mu, sigma);
  }
  else{ // == Vavilov
    this_vav =  vav.Pdf(y, par[0], par[1]);
    this_vav =  this_vav / a;
  }

  // == Vavilov PDF only - out of range for very low kappa values
  //this_vav =  vav.Pdf(y, par[0], par[1]);
  //this_vav =  this_vav / a;

  return this_vav;
}


double BetheBloch::dEdx_PDF(double KE, double pitch, double dEdx){

  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double this_xi = Landau_xi(KE, pitch);
  double this_Wmax = Get_Wmax(KE);
  double this_kappa = this_xi / this_Wmax;
  double this_dEdx_BB = meandEdx(KE);
  double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, pitch};
  
  TF1 *PDF = new TF1("", dEdx_PDF_fuction, -100., 1000., 5);
  PDF -> SetParameters(par[0], par[1], par[2], par[3], par[4]);

  double out = PDF -> Eval(dEdx);
  delete PDF;
  return out;
}


double BetheBloch::Fit_Proton_Residual_Length_Likelihood(const vector<double> dEdx, const vector<double> ResRange, int PID, bool save_graph){
  // == only for protons
  int abs_PID = abs(PID);
  //if(!abs(PID) == 2212){
  if(PID != 2212){
    return -9999.;
  }
  double best_additional_res_length = -0.1;
  double best_m2lnL = 99999.;
  double min_additional_res_length = 0.; // == [cm]
  double max_additional_res_length = 450.; // == [cm]
  //double res_length_step = 1.0; // == [cm]
  double res_length_step = 0.105; // == [cm //HY:Set step size to 0.105 cm
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;
  int N_skip = 3;
  int this_N_calo = dEdx.size();
  if(this_N_calo <= 15){
    return -9999.; // == Too small number of hits 
  }
  int this_N_hits = this_N_calo; // == Use how many hits
  int i_bestfit = -1;
  double dEdx_truncate_upper = 5.;
  double dEdx_truncate_bellow = 0.5;
  vector<double> m2lnL_vector;
  vector<double> additional_res_legnth_vector;
  for(int i = 0; i < res_length_trial; i++){

    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_m2lnL = 0.;
    for(int j = N_skip; j < this_N_hits - N_skip; j++){ // == Do not use first and last N_skip hits 
      double this_res_length = ResRange.at(j) + this_additional_res_length;
      //double this_KE = map_BB[abs_PID]->KEFromRangeSpline(this_res_length);
      double this_KE = KEFromRangeSpline(this_res_length);
      //double dEdx_theory = map_BB[abs_PID]->meandEdx(this_KE);
      double dEdx_theory = meandEdx(this_KE);
      double dEdx_measured = dEdx.at(j);
      if(dEdx_measured < dEdx_truncate_bellow || dEdx_measured > dEdx_truncate_upper) continue; // == Truncate

      // == Likelihood
      double this_pitch = fabs(ResRange.at(j - 1) - ResRange.at(j + 1)) / 2.0;
      //cout << "[HadAna::Fit_Pion_Residual_Length_Likelihood] " << j << ", ResRange.at(j) : " << ResRange.at(j) << ", this_KE : " << this_KE << ", this_pitch : " << this_pitch << ", dEdx_measured : " << dEdx_measured << endl;
      //double this_likelihood = map_BB[abs_PID] -> dEdx_PDF(this_KE, this_pitch, dEdx_measured);
      double this_likelihood = dEdx_PDF(this_KE, this_pitch, dEdx_measured);
      //cout << "[HadAna::Fit_Pion_Residual_Length_Likelihood] dEdx_measured : " << dEdx_measured << ", dEdx_PDF : " << this_likelihood << ", log(this_likelihood) : " << log(this_likelihood) << endl;
      this_m2lnL += (-2.0) * log(this_likelihood);
    }
    //this_chi2 = this_chi2 / (this_N_hits + 0.); // == chi2 / n.d.f                                                                                                                                       
    if(this_m2lnL < best_m2lnL){
      best_m2lnL = this_m2lnL;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }
    if(save_graph){
      // == Save vectors for graphes
      m2lnL_vector.push_back(this_m2lnL);
      additional_res_legnth_vector.push_back(this_additional_res_length);
    }
  }

  if(save_graph){
/*
    vector<double> range_original;
    vector<double> range_bestfit;
    vector<double> range_reco;
    vector<double> dEdx_ordered;
    for(int i = N_skip; i < this_N_hits - N_skip; i++){
      double this_range_original = ResRange.at(i);
      range_original.push_back(this_range_original);

      double this_range_bestfit = ResRange.at(i) + best_additional_res_length;
      range_bestfit.push_back(this_range_bestfit);
      range_reco.push_back(ResRange.at(i));
      dEdx_ordered.push_back(dEdx.at(i));
    }
    TGraph *dEdx_gr = new TGraph(this_N_hits - 10, &range_original[0], &dEdx_ordered[0]);
    dEdx_gr -> SetName(Form("dEdx_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits - 1));
    dEdx_gr -> Write();
    delete dEdx_gr;

    TGraph *dEdx_bestfit_gr = new TGraph(this_N_hits - 10,&range_bestfit[0], &dEdx_ordered[0]);
    dEdx_bestfit_gr -> SetName(Form("dEdx_bestfit_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    dEdx_bestfit_gr -> Write();
    delete dEdx_bestfit_gr;

    TGraph *dEdx_reco_gr = new TGraph(this_N_hits - 10,&range_reco[0], &dEdx_ordered[0]);
    dEdx_reco_gr -> SetName(Form("dEdx_reco_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    dEdx_reco_gr -> Write();
    delete dEdx_reco_gr;

    TGraph *chi2_gr = new TGraph(additional_res_legnth_vector.size(), &additional_res_legnth_vector[0], &m2lnL_vector[0]);
    chi2_gr -> SetName(Form("Chi2_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    chi2_gr -> Write();
    m2lnL_vector.clear();
    additional_res_legnth_vector.clear();
    delete chi2_gr;
*/
  }

  double original_res_length = ResRange.at(this_N_calo - 1); // == [cm]
  double best_total_res_length = best_additional_res_length + original_res_length;
  //double best_KE = map_BB[abs_PID]->KEFromRangeSpline(best_total_res_length);
  double best_KE = KEFromRangeSpline(best_total_res_length);
  //double best_mom = map_BB[abs_PID]->KEtoMomentum(best_KE);
  //double best_mom = KEtoMomentum(best_KE);

  // == Define fitting failed cases
  if(i_bestfit == res_length_trial - 1){
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : no mimumum" << endl; 
    m2lnL_vector.clear();
    additional_res_legnth_vector.clear();
    return -9999.;
  }
  else if(best_m2lnL > 99990.){
    m2lnL_vector.clear();
    additional_res_legnth_vector.clear();
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : best_chi2 > 99990." << endl;
    return -9999.;
  }
  else if(best_m2lnL < 1.0e-11){
    m2lnL_vector.clear();
    additional_res_legnth_vector.clear();
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : best_chi2 < 1.0e-11" << endl;
    return -9999.;
  }
  m2lnL_vector.clear();
  additional_res_legnth_vector.clear();
  return best_total_res_length;
}


