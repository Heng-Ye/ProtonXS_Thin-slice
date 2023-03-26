//#ifndef UNFOLD_H
//#define UNFOLD_H

R__LOAD_LIBRARY(libRooUnfold.so) //load share lib
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"

class Unfold {
  
 public:

  //Unfold(int nb, double xlo, double xhi, TH2D* h2d_reco, TH2D* h2d_true);
  Unfold(int nb, double xlo, double xhi, TH2D* h2d_reco, TH2D* h2d_true, TH2D* h2d_reco2, TH2D* h2d_true2);

  RooUnfoldResponse response_SliceID_Int;  //Interaction
  RooUnfoldResponse response_SliceID_Inc;  //Incident
  RooUnfoldResponse response_st_SliceID_Inc;  //Incident(start sliceID)
  RooUnfoldResponse response_SliceID_2D;   //2D response matrix for INC: INC_ST(x-axis), INC_Int(y-axis)
  RooUnfoldResponse response_SliceID_Int_2D;  //2D response matrix for INT: INC_Int(x-axis), INT(y-axis)

  //RooUnfoldResponse response_SliceID_shape_Int;  //Interaction [by shape]
  //RooUnfoldResponse response_SliceID_shape_Inc;  //Incident [by shape]

  TH1D *eff_num_Int; //Interaction efficiency numerator
  TH1D *eff_den_Int; //Interaction efficiency denominator

  TH1D *eff_num_Inc; //Incident efficiency numerator
  TH1D *eff_den_Inc; //Incident efficiency denominator

  TH1D *eff_num_st_Inc; //Incident efficiency numerator
  TH1D *eff_den_st_Inc; //Incident efficiency denominator


  TH1D *pur_num_Int; //Interaction purity numerator
  TH1D *pur_num_Inc; //Incident purity numerator
  TH1D *pur_num_st_Inc; //Incident purity numerator

  TH1D *pur_den;     //Interaction/incident purity denominator
  TH1D *pur_den_Inc;     //Incident purity denominator
  TH1D *pur_den_st_Inc;     //Incident purity denominator
  TH1D *pur_den_Int;     //Interaction purity denominator

  TH1D *eff_Int;
  TH1D *eff_Inc;
  TH1D *eff_st_Inc;
  TH1D *pur_Int;
  TH1D *pur_Inc;
  TH1D *pur_st_Inc;

  TH1D *res_Inc_reco;	
  TH1D *res_st_Inc_reco;	
  TH1D *res_Inc_truth;	
  TH1D *res_st_Inc_truth;	
  TH1D *res_Int_reco;	
  TH1D *res_Int_truth;	

  void SaveHistograms();

};


//Unfold::Unfold(int nb, double xlo, double xhi)
//Unfold::Unfold(int nb, double xlo, double xhi, TH2D* h2d_reco, TH2D* h2d_true)
Unfold::Unfold(int nb, double xlo, double xhi, TH2D* h2d_reco, TH2D* h2d_true, TH2D* h2d_reco2, TH2D* h2d_true2)
  : response_SliceID_Int(nb, xlo, xhi)
  , response_SliceID_Inc(nb, xlo, xhi)
  , response_st_SliceID_Inc(nb, xlo, xhi)
  , response_SliceID_2D(h2d_reco, h2d_true)
  , response_SliceID_Int_2D(h2d_reco2, h2d_true2)
{

  response_SliceID_Int.UseOverflow(false);
  response_SliceID_Inc.UseOverflow(false);
  response_st_SliceID_Inc.UseOverflow(false);
  response_SliceID_2D.UseOverflow(false);
  response_SliceID_Int_2D.UseOverflow(false);

  eff_num_Int = new TH1D("eff_num_Int", "eff_num_Int", nb, xlo, xhi);
  eff_den_Int = new TH1D("eff_den_Int", "eff_den_Int", nb, xlo, xhi);
  eff_num_Inc = new TH1D("eff_num_Inc", "eff_num_Inc", nb, xlo, xhi);
  eff_den_Inc = new TH1D("eff_den_Inc", "eff_den_Inc", nb, xlo, xhi);

  eff_num_st_Inc = new TH1D("eff_num_st_Inc", "eff_num_st_Inc", nb, xlo, xhi);
  eff_den_st_Inc = new TH1D("eff_den_st_Inc", "eff_den_st_Inc", nb, xlo, xhi);
  pur_num_Int = new TH1D("pur_num_Int", "pur_num_Int", nb, xlo, xhi);
  pur_num_Inc = new TH1D("pur_num_Inc", "pur_num_Inc", nb, xlo, xhi);
  pur_num_st_Inc = new TH1D("pur_num_st_Inc", "pur_num_st_Inc", nb, xlo, xhi);
  pur_den     = new TH1D("pur_den",     "pur_den",     nb, xlo, xhi);
  pur_den_Inc = new TH1D("pur_den_Inc", "pur_den_Inc", nb, xlo, xhi);
  pur_den_st_Inc = new TH1D("pur_den_st_Inc", "pur_den_st_Inc", nb, xlo, xhi);
  pur_den_Int = new TH1D("pur_den_Int", "pur_den_Int", nb, xlo, xhi);

  res_Inc_reco = new TH1D("res_Inc_reco", "res_Inc_reco", nb, xlo, xhi);
  res_st_Inc_reco = new TH1D("res_st_Inc_reco", "res_st_Inc_reco", nb, xlo, xhi);
  res_Inc_truth = new TH1D("res_Inc_truth", "res_Inc_truth", nb, xlo, xhi);	
  res_st_Inc_truth = new TH1D("res_st_Inc_truth", "res_st_Inc_truth", nb, xlo, xhi);	
  res_Int_reco = new TH1D("res_Int_reco", "res_Int_reco", nb, xlo, xhi);
  res_Int_truth = new TH1D("res_Int_truth", "res_Int_truth", nb, xlo, xhi);	

  eff_num_Int->Sumw2();
  eff_den_Int->Sumw2();
  eff_num_Inc->Sumw2();
  eff_num_st_Inc->Sumw2();
  eff_den_Inc->Sumw2();
  eff_den_st_Inc->Sumw2();
  pur_num_Int->Sumw2();
  pur_num_Inc->Sumw2();
  pur_num_st_Inc->Sumw2();
  pur_den->Sumw2();
  pur_den_Inc->Sumw2();
  pur_den_st_Inc->Sumw2();
  pur_den_Int->Sumw2();

  res_Inc_reco->Sumw2();
  res_st_Inc_reco->Sumw2();
  res_Inc_truth->Sumw2();
  res_st_Inc_truth->Sumw2();
  res_Int_reco->Sumw2();
  res_Int_truth->Sumw2();	


}  

void Unfold::SaveHistograms(){

  response_SliceID_Int.Write("response_SliceID_Int");
  response_SliceID_Inc.Write("response_SliceID_Inc");
  response_st_SliceID_Inc.Write("response_st_SliceID_Inc");

  response_SliceID_2D.Write("response_SliceID_2D");
  response_SliceID_Int_2D.Write("response_SliceID_Int_2D");

  eff_num_Int->Write("eff_num_Int");
  eff_den_Int->Write("eff_den_Int");
  eff_num_Inc->Write("eff_num_Inc");
  eff_num_st_Inc->Write("eff_num_st_Inc");
  eff_den_Inc->Write("eff_den_Inc");
  eff_den_st_Inc->Write("eff_den_st_Inc");
  pur_num_Int->Write("pur_num_Int");
  pur_num_Inc->Write("pur_num_Inc");
  pur_num_st_Inc->Write("pur_num_st_Inc");
  pur_den->Write("pur_den");
  pur_den_Inc->Write("pur_den_Inc");
  pur_den_st_Inc->Write("pur_den_st_Inc");
  pur_den_Int->Write("pur_den_Int");

  eff_Int = (TH1D*)eff_num_Int->Clone("eff_Int");
  eff_Int->Divide(eff_den_Int);
  eff_Int->Write("eff_Int");

  eff_Inc = (TH1D*)eff_num_Inc->Clone("eff_Inc");
  eff_Inc->Divide(eff_den_Inc);
  eff_Inc->Write("eff_Inc");

  eff_st_Inc = (TH1D*)eff_num_st_Inc->Clone("eff_st_Inc");
  eff_st_Inc->Divide(eff_den_st_Inc);
  eff_st_Inc->Write("eff_st_Inc");

  pur_Int = (TH1D*)pur_num_Int->Clone("pur_Int");
  pur_Int->Divide(pur_den_Int);
  pur_Int->Write("pur_Int");

  pur_Inc = (TH1D*)pur_num_Inc->Clone("pur_Inc");
  //pur_Inc->Divide(pur_den);
  pur_Inc->Divide(pur_den_Inc);
  pur_Inc->Write("pur_Inc");

  pur_st_Inc = (TH1D*)pur_num_st_Inc->Clone("pur_st_Inc");
  pur_st_Inc->Divide(pur_den_st_Inc);
  pur_st_Inc->Write("pur_st_Inc");


  //TH2D *hint = (TH2D*)response_SliceID_Int.Hresponse();
  //hint->SetTitle("Proton Inelastic Scatterings;Reco Slice ID;True Slice ID");
  //hint->Write("response_SliceID_Int");

  //TH2D *hinc = (TH2D*)response_SliceID_Inc.Hresponse();
  //hinc->SetTitle("All Protons; Reco Slice ID; True Slice ID");
  //hinc->Write("response_SliceID_Inc");


  res_Inc_reco->Write("res_Inc_reco");
  res_st_Inc_reco->Write("res_st_Inc_reco");
  res_Inc_truth->Write("res_Inc_truth");
  res_st_Inc_truth->Write("res_st_Inc_truth");
  res_Int_reco->Write("res_Int_reco");
  res_Int_truth->Write("res_Int_truth");	

}

//#endif
