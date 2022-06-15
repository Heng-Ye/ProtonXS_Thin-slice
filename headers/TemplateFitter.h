#include "TMinuit.h"
#include "TH1D.h"
#include <iostream>

class TH1D;

class TemplateFitter {
  
 public:
  
  TemplateFitter();
  
  void SetHistograms(TH1D *hist0, TH1D *hist1, TH1D *hist2);
  static void fcn(int &npar, double *gin, double &f, double *par, int iflag);
  void SetFitRange(int imin, int imax);

  void Fit();
  
  double GetPar();
  double GetParError();

  static TH1D *h0, *h1, *h2;
  static int i0, i1;

  TMinuit *gMinuit;

 private:

  bool fitsuccess;
  
};

TH1D* TemplateFitter::h0;
TH1D* TemplateFitter::h1;
TH1D* TemplateFitter::h2;

int TemplateFitter::i0;
int TemplateFitter::i1;


TemplateFitter::TemplateFitter(){
  gMinuit = new TMinuit(1);
  gMinuit->SetFCN(fcn);
  double arglist[10];
  int ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR",arglist,1,ierflg);
}

void TemplateFitter::SetHistograms(TH1D *hist0, TH1D *hist1, TH1D *hist2){
  h0 = hist0;
  h1 = hist1;
  h2 = hist2;

  i0 = 1;
  i1 = h0->GetNbinsX();

}

void TemplateFitter::SetFitRange(int imin, int imax){
  if (imin>=0) i0 = imin;
  if (imax>=0) i1 = imax;
}

void TemplateFitter::fcn(int &npar, double *gin, double &f, double *par, int iflag){
  
  double chisq = 0;

  for (int i = i0; i<=i1; ++i){
    double x = h0->GetBinContent(i);
    double y = h1->GetBinContent(i);
    double z = h2->GetBinContent(i);
    
    double ex = h0->GetBinError(i);
    double ey = h1->GetBinError(i);
    double ez = h2->GetBinError(i);

    if (!ex) ex = 1;
    chisq += pow(x-y-z*par[0],2)/(pow(ex,2)+pow(ey,2)+pow(ez*par[0],2));
  }

  f = chisq;

}

void TemplateFitter::Fit(){

  double arglist[10];
  int ierflg = 0;

  gMinuit->mncler();
  double vstart = 1;
  double step = 0.01;
  gMinuit->mnparm(0,"corr_fact",vstart,step,0,10,ierflg);

  fitsuccess = false;

  if (h0->Integral(i0,i1) && h2->Integral(i0,i1)){
    arglist[0] = 500;
    arglist[1] = 1;
    gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
    
    double par, epar;
    gMinuit->GetParameter(0,par,epar);

    TString test =  gMinuit->fCstatu.Data(); 
    if (test.EqualTo("CONVERGED ")){
      std::cout<<"Best fit = "<<par<<" error = "<<epar<<std::endl;
      fitsuccess = true;
    }
  }
  else{
    std::cout<<"No fit was done because data and/or template are empty."<<std::endl;
  }
}

double TemplateFitter::GetPar(){
  if (fitsuccess){
    double par, epar;
    gMinuit->GetParameter(0,par,epar);
    return par;
  }
  else{
    return 1;
  }
}

double TemplateFitter::GetParError(){
  if (fitsuccess){
    double par, epar;
    gMinuit->GetParameter(0,par,epar);
    return epar;
  }
  else{
    return 2;
  }
}


