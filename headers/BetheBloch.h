#include "TF1.h"
#include "TSpline.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <map>

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

  double KEAtLength(double KE0, double tracklength);

  void CreateSplineAtKE(int iKE);

  
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
	cout<<"Negative delta E: "<<deltaE<<endl; //meaning in the region of KE<500 keV
	tmp_ke_at_len=0;
  }	
  if (tmp_ke_at_len < 0) { 
        tmp_ke_at_len=0;
	cout<<"Negative KE: "<<KE0 - deltaE<<" --> set KE to zero"<<endl;
  }
  
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
    cout<<"KE too low: "<<iKE<<endl;
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
