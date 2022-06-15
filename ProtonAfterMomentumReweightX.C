#define ProtonAfterMomentumReweight_cxx
#include "ProtonAfterMomentumReweight.h"

#include <TH2.h>
#include <TH1.h>
#include "TH2D.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>

#include <TMath.h>
#include <TLine.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH3F.h>
#include <TString.h>
#include <TProfile2D.h>
#include <THStack.h>
#include "TGraph.h" 
#include "TGraphSmooth.h" 
#include "TParameter.h"
#include "TGraphErrors.h"
#include "string"
#include "vector"
#include "TSpline.h"
#include "TH3F.h"
#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>
#include "TVector3.h"

#include <stdio.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "./cali/dedx_function_35ms.h"
#include "./headers/BasicParameters.h"
#include "./headers/BasicFunctions.h"
#include "./headers/util.h"
//#include "./headers/SliceParams.h"
//#include "./headers/ThinSlice.h"
//#include "./headers/sce_map.h"

using namespace std;
using namespace ROOT::Math;

   TSpline3 *spline_dx_fwd_neg[31][37];
   TSpline3 *spline_dy_fwd_neg[19][37];
   TSpline3 *spline_dz_fwd_neg[19][31];

   TSpline3 *spline_dx_bkwd_neg[31][37];
   TSpline3 *spline_dy_bkwd_neg[19][37];
   TSpline3 *spline_dz_bkwd_neg[19][31];

   TSpline3 *spline_dEx_neg[31][37];
   TSpline3 *spline_dEy_neg[19][37];
   TSpline3 *spline_dEz_neg[19][31];

   TSpline3 *spline_dx_fwd_pos[31][37];
   TSpline3 *spline_dy_fwd_pos[19][37];
   TSpline3 *spline_dz_fwd_pos[19][31];

   TSpline3 *spline_dx_bkwd_pos[31][37];
   TSpline3 *spline_dy_bkwd_pos[19][37];
   TSpline3 *spline_dz_bkwd_pos[19][31];

   TSpline3 *spline_dEx_pos[31][37];
   TSpline3 *spline_dEy_pos[19][37];
   TSpline3 *spline_dEz_pos[19][31];

TSpline3* MakeSpline(TH3F* spline_hist, int dim1, int dim2_bin, int dim3_bin, int maptype, int driftvol) {
  TSpline3 *spline = 0;
  
  if(dim1 == 1)
  {
    double a[19];
    double b[19];
    for(int x = 1; x <= 19; x++)
    {
      a[x-1] = spline_hist->GetXaxis()->GetBinCenter(x);
      b[x-1] = spline_hist->GetBinContent(x,dim2_bin,dim3_bin);
    }

    if(maptype == 1)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,19,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 2)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,19,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 3)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,19,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
  }
  else if(dim1 == 2)
  {
    double a[31];
    double b[31];
    for(int y = 1; y <= 31; y++)
    {
      a[y-1] = spline_hist->GetYaxis()->GetBinCenter(y);
      b[y-1] = spline_hist->GetBinContent(dim2_bin,y,dim3_bin);
    }

    if(maptype == 1)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,31,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 2)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,31,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 3)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,31,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
  }
  else if(dim1 == 3)
  {
    double a[37];
    double b[37];
    for(int z = 1; z <= 37; z++)
    {
      a[z-1] = spline_hist->GetZaxis()->GetBinCenter(z);
      b[z-1] = spline_hist->GetBinContent(dim2_bin,dim3_bin,z);
    }

    if(maptype == 1)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,37,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 2)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,37,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
    else if(maptype == 3)
    {
      spline = new TSpline3(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol),a,b,37,"b2e2",0,0);
      spline->SetName(Form("spline_%d_%d_%d_%d_%d",dim1,dim2_bin,dim3_bin,maptype,driftvol));
    }
  }

  return spline;
}


double InterpolateSplines(TH3F* interp_hist, double xVal, double yVal, double zVal, int dim, int maptype, int driftvol) {
  int bin_x = interp_hist->GetXaxis()->FindBin(xVal);
  int bin_y = interp_hist->GetYaxis()->FindBin(yVal);
  int bin_z = interp_hist->GetZaxis()->FindBin(zVal);

  int bincenter_x = interp_hist->GetXaxis()->GetBinCenter(bin_x);
  int bincenter_y = interp_hist->GetYaxis()->GetBinCenter(bin_y);
  int bincenter_z = interp_hist->GetZaxis()->GetBinCenter(bin_z);

  int max_x = interp_hist->GetNbinsX();
  int max_y = interp_hist->GetNbinsY();
  int max_z = interp_hist->GetNbinsZ();
  
  int low_x;
  int high_x;
  if(bin_x <= 1)
  {
    low_x = 1;
    high_x = 2;
  }
  else if(bin_x >= max_x)
  {
    low_x = max_x-1;
    high_x = max_x;
  }
  else if(xVal > bincenter_x)
  {
    low_x = bin_x;
    high_x = bin_x+1;
  }
  else
  {
    low_x = bin_x-1;
    high_x = bin_x;
  }

  int low_y;
  int high_y;
  if(bin_y <= 1)
  {
    low_y = 1;
    high_y = 2;
  }
  else if(bin_y >= max_y)
  {
    low_y = max_y-1;
    high_y = max_y;
  }
  else if(yVal > bincenter_y)
  {
    low_y = bin_y;
    high_y = bin_y+1;
  }
  else
  {
    low_y = bin_y-1;
    high_y = bin_y;
  }

  int low_z;
  int high_z;
  if(bin_z <= 1)
  {
    low_z = 1;
    high_z = 2;
  }
  else if(bin_z >= max_z)
  {
    low_z = max_z-1;
    high_z = max_z;
  }
  else if(zVal > bincenter_z)
  {
    low_z = bin_z;
    high_z = bin_z+1;
  }
  else
  {
    low_z = bin_z-1;
    high_z = bin_z;
  }

  double interp_val = 0.0;
  
  if(dim == 1)
  {
    double a_1 = interp_hist->GetYaxis()->GetBinCenter(low_y);
    double a_2 = interp_hist->GetYaxis()->GetBinCenter(high_y);

    double b_1 = interp_hist->GetZaxis()->GetBinCenter(low_z);
    double b_2 = interp_hist->GetZaxis()->GetBinCenter(high_z);

    double f_11 = 0.0;
    double f_12 = 0.0;
    double f_21 = 0.0;
    double f_22 = 0.0;
    if(driftvol == 1)
    {
      if(maptype == 1)
      {
        f_11 = spline_dx_fwd_neg[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_fwd_neg[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_fwd_neg[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_fwd_neg[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dx_bkwd_neg[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_bkwd_neg[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_bkwd_neg[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_bkwd_neg[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEx_neg[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dEx_neg[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dEx_neg[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dEx_neg[high_y-1][high_z-1]->Eval(xVal);
      }
    }
    else if(driftvol == 2)
    {
      if(maptype == 1)
      {
        f_11 = spline_dx_fwd_pos[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_fwd_pos[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_fwd_pos[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_fwd_pos[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dx_bkwd_pos[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dx_bkwd_pos[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dx_bkwd_pos[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dx_bkwd_pos[high_y-1][high_z-1]->Eval(xVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEx_pos[low_y-1][low_z-1]->Eval(xVal);
        f_12 = spline_dEx_pos[low_y-1][high_z-1]->Eval(xVal);
        f_21 = spline_dEx_pos[high_y-1][low_z-1]->Eval(xVal);
        f_22 = spline_dEx_pos[high_y-1][high_z-1]->Eval(xVal);
      }
    }

    interp_val = (f_11*(a_2-yVal)*(b_2-zVal) + f_21*(yVal-a_1)*(b_2-zVal) + f_12*(a_2-yVal)*(zVal-b_1) + f_22*(yVal-a_1)*(zVal-b_1))/((a_2-a_1)*(b_2-b_1));
  }
  else if(dim == 2)
  {
    double a_1 = interp_hist->GetXaxis()->GetBinCenter(low_x);
    double a_2 = interp_hist->GetXaxis()->GetBinCenter(high_x);

    double b_1 = interp_hist->GetZaxis()->GetBinCenter(low_z);
    double b_2 = interp_hist->GetZaxis()->GetBinCenter(high_z);

    double f_11 = 0.0;
    double f_12 = 0.0;
    double f_21 = 0.0;
    double f_22 = 0.0;
    if(driftvol == 1)
    {
      if(maptype == 1)
      {
        f_11 = spline_dy_fwd_neg[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_fwd_neg[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_fwd_neg[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_fwd_neg[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dy_bkwd_neg[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_bkwd_neg[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_bkwd_neg[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_bkwd_neg[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEy_neg[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dEy_neg[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dEy_neg[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dEy_neg[high_x-1][high_z-1]->Eval(yVal);
      }
    }
    else if(driftvol == 2)
    {
      if(maptype == 1)
      {
        f_11 = spline_dy_fwd_pos[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_fwd_pos[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_fwd_pos[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_fwd_pos[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dy_bkwd_pos[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dy_bkwd_pos[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dy_bkwd_pos[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dy_bkwd_pos[high_x-1][high_z-1]->Eval(yVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEy_pos[low_x-1][low_z-1]->Eval(yVal);
        f_12 = spline_dEy_pos[low_x-1][high_z-1]->Eval(yVal);
        f_21 = spline_dEy_pos[high_x-1][low_z-1]->Eval(yVal);
        f_22 = spline_dEy_pos[high_x-1][high_z-1]->Eval(yVal);
      }
    }

    interp_val = (f_11*(a_2-xVal)*(b_2-zVal) + f_21*(xVal-a_1)*(b_2-zVal) + f_12*(a_2-xVal)*(zVal-b_1) + f_22*(xVal-a_1)*(zVal-b_1))/((a_2-a_1)*(b_2-b_1));
  }
  else if(dim == 3)
  {
    double a_1 = interp_hist->GetXaxis()->GetBinCenter(low_x);
    double a_2 = interp_hist->GetXaxis()->GetBinCenter(high_x);

    double b_1 = interp_hist->GetYaxis()->GetBinCenter(low_y);
    double b_2 = interp_hist->GetYaxis()->GetBinCenter(high_y);

    double f_11 = 0.0;
    double f_12 = 0.0;
    double f_21 = 0.0;
    double f_22 = 0.0;
    if(driftvol == 1)
    {
      if(maptype == 1)
      {
        f_11 = spline_dz_fwd_neg[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_fwd_neg[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_fwd_neg[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_fwd_neg[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dz_bkwd_neg[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_bkwd_neg[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_bkwd_neg[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_bkwd_neg[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEz_neg[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dEz_neg[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dEz_neg[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dEz_neg[high_x-1][high_y-1]->Eval(zVal);
      }
    }
    else if(driftvol == 2)
    {
      if(maptype == 1)
      {
        f_11 = spline_dz_fwd_pos[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_fwd_pos[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_fwd_pos[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_fwd_pos[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 2)
      {
        f_11 = spline_dz_bkwd_pos[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dz_bkwd_pos[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dz_bkwd_pos[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dz_bkwd_pos[high_x-1][high_y-1]->Eval(zVal);
      }
      else if(maptype == 3)
      {
        f_11 = spline_dEz_pos[low_x-1][low_y-1]->Eval(zVal);
        f_12 = spline_dEz_pos[low_x-1][high_y-1]->Eval(zVal);
        f_21 = spline_dEz_pos[high_x-1][low_y-1]->Eval(zVal);
        f_22 = spline_dEz_pos[high_x-1][high_y-1]->Eval(zVal);
      }
    }

    interp_val = (f_11*(a_2-xVal)*(b_2-yVal) + f_21*(xVal-a_1)*(b_2-yVal) + f_12*(a_2-xVal)*(yVal-b_1) + f_22*(xVal-a_1)*(yVal-b_1))/((a_2-a_1)*(b_2-b_1));
  }

  return interp_val;
}


void LoadHist(){








  TFile *infile = TFile::Open("/cvmfs/dune.opensciencegrid.org/products/dune/dune_pardata/v01_66_00/SpaceChargeProtoDUNE/SCE_DataDriven_180kV_v4.root");


  //Load in files
  TH3F* hDx_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_X_Pos");
  TH3F* hDy_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Y_Pos");
  TH3F* hDz_sim_pos_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Z_Pos");
  TH3F* hEx_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_X_Pos");
  TH3F* hEy_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_Y_Pos");
  TH3F* hEz_sim_pos_orig = (TH3F*)infile->Get("Reco_ElecField_Z_Pos");
  
  TH3F* hDx_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_X_Neg");
  TH3F* hDy_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Y_Neg");
  TH3F* hDz_sim_neg_orig = (TH3F*)infile->Get("RecoFwd_Displacement_Z_Neg");
  TH3F* hEx_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_X_Neg");
  TH3F* hEy_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_Y_Neg");
  TH3F* hEz_sim_neg_orig = (TH3F*)infile->Get("Reco_ElecField_Z_Neg");
  
  TH3F* hDx_sim_pos = (TH3F*)hDx_sim_pos_orig->Clone("hDx_pos");
  TH3F* hDy_sim_pos = (TH3F*)hDy_sim_pos_orig->Clone("hDy_pos");
  TH3F* hDz_sim_pos = (TH3F*)hDz_sim_pos_orig->Clone("hDz_pos");
  TH3F* hEx_sim_pos = (TH3F*)hEx_sim_pos_orig->Clone("hEx_pos");
  TH3F* hEy_sim_pos = (TH3F*)hEy_sim_pos_orig->Clone("hEy_pos");
  TH3F* hEz_sim_pos = (TH3F*)hEz_sim_pos_orig->Clone("hEz_pos");
  
  TH3F* hDx_sim_neg = (TH3F*)hDx_sim_neg_orig->Clone("hDx_neg");
  TH3F* hDy_sim_neg = (TH3F*)hDy_sim_neg_orig->Clone("hDy_neg");
  TH3F* hDz_sim_neg = (TH3F*)hDz_sim_neg_orig->Clone("hDz_neg");
  TH3F* hEx_sim_neg = (TH3F*)hEx_sim_neg_orig->Clone("hEx_neg");
  TH3F* hEy_sim_neg = (TH3F*)hEy_sim_neg_orig->Clone("hEy_neg");
  TH3F* hEz_sim_neg = (TH3F*)hEz_sim_neg_orig->Clone("hEz_neg");
  
  hDx_sim_pos->SetDirectory(0);
  hDy_sim_pos->SetDirectory(0);
  hDz_sim_pos->SetDirectory(0);
  hEx_sim_pos->SetDirectory(0);
  hEy_sim_pos->SetDirectory(0);
  hEz_sim_pos->SetDirectory(0);
  
  hDx_sim_neg->SetDirectory(0);
  hDy_sim_neg->SetDirectory(0);
  hDz_sim_neg->SetDirectory(0);
  hEx_sim_neg->SetDirectory(0);
  hEy_sim_neg->SetDirectory(0);
  hEz_sim_neg->SetDirectory(0);

  for(int y = 1; y <= 31; y++){
    for(int z = 1; z <= 37; z++){
      spline_dx_fwd_neg[y-1][z-1] = MakeSpline(hDx_sim_neg,1,y,z,1,1);
      spline_dx_fwd_pos[y-1][z-1] = MakeSpline(hDx_sim_pos,1,y,z,1,2);
      spline_dEx_neg[y-1][z-1] = MakeSpline(hEx_sim_neg,1,y,z,3,1);
      spline_dEx_pos[y-1][z-1] = MakeSpline(hEx_sim_pos,1,y,z,3,2);
    }
  }
  for(int x = 1; x <= 19; x++){
    for(int z = 1; z <= 37; z++){
      spline_dy_fwd_neg[x-1][z-1] = MakeSpline(hDy_sim_neg,2,x,z,1,1);
      spline_dy_fwd_pos[x-1][z-1] = MakeSpline(hDy_sim_pos,2,x,z,1,2);
      spline_dEy_neg[x-1][z-1] = MakeSpline(hEy_sim_neg,2,x,z,3,1);
      spline_dEy_pos[x-1][z-1] = MakeSpline(hEy_sim_pos,2,x,z,3,2);
    }
  }
  for(int x = 1; x <= 19; x++){
    for(int y = 1; y <= 31; y++){
      spline_dz_fwd_neg[x-1][y-1] = MakeSpline(hDz_sim_neg,3,x,y,1,1);
      spline_dz_fwd_pos[x-1][y-1] = MakeSpline(hDz_sim_pos,3,x,y,1,2);
      spline_dEz_neg[x-1][y-1] = MakeSpline(hEz_sim_neg,3,x,y,3,1);
      spline_dEz_pos[x-1][y-1] = MakeSpline(hEz_sim_pos,3,x,y,3,2);
    }
  }
}



void ProtonAfterMomentumReweight::Loop() {
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntries();
	std::cout<<"nentries: "<<nentries<<std::endl;

	Long64_t nbytes = 0, nb = 0;
	bool isTestSample=true;
	//int true_sliceID = -1, reco_sliceID = -1;

   	LoadHist();
 	TFile f2("/cvmfs/dune.opensciencegrid.org/products/dune/dune_pardata/v01_66_00/SpaceChargeProtoDUNE/SCE_DataDriven_180kV_v4.root");
   	TH3F *RecoFwd_Displacement_Z_Neg = (TH3F*)f2.Get("RecoFwd_Displacement_Z_Neg");
   	TH3F *RecoFwd_Displacement_Z_Pos = (TH3F*)f2.Get("RecoFwd_Displacement_Z_Pos");	


	//book histograms --------------------------------------------------------------------------------------------//
	//int n_b=150;
	int n_b=30;
	double b_min=0;
	double b_max=150;

	TH1D *h1d_trklen_CaloSz=new TH1D(Form("h1d_trklen_CaloSz"), Form("MC CaloSz"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_el=new TH1D(Form("h1d_trklen_CaloSz_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_inel=new TH1D(Form("h1d_trklen_CaloSz_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_midcosmic=new TH1D(Form("h1d_trklen_CaloSz_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_midpi=new TH1D(Form("h1d_trklen_CaloSz_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_midp=new TH1D(Form("h1d_trklen_CaloSz_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_midmu=new TH1D(Form("h1d_trklen_CaloSz_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_mideg=new TH1D(Form("h1d_trklen_CaloSz_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_midother=new TH1D(Form("h1d_trklen_CaloSz_midother"), Form("midother"), n_b, b_min, b_max);

	TH1D *h1d_trklen_CaloSz_bmrw=new TH1D(Form("h1d_trklen_CaloSz_bmrw"), Form("MC CaloSz_bmrw"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_bmrw_el=new TH1D(Form("h1d_trklen_CaloSz_bmrw_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_bmrw_inel=new TH1D(Form("h1d_trklen_CaloSz_bmrw_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_bmrw_midcosmic=new TH1D(Form("h1d_trklen_CaloSz_bmrw_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_bmrw_midpi=new TH1D(Form("h1d_trklen_CaloSz_bmrw_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_bmrw_midp=new TH1D(Form("h1d_trklen_CaloSz_bmrw_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_bmrw_midmu=new TH1D(Form("h1d_trklen_CaloSz_bmrw_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_bmrw_mideg=new TH1D(Form("h1d_trklen_CaloSz_bmrw_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_bmrw_midother=new TH1D(Form("h1d_trklen_CaloSz_bmrw_midother"), Form("midother"), n_b, b_min, b_max);

	TH1D *h1d_trklen_Pos=new TH1D(Form("h1d_trklen_Pos"), Form("MC default"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_el=new TH1D(Form("h1d_trklen_Pos_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_inel=new TH1D(Form("h1d_trklen_Pos_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_midcosmic=new TH1D(Form("h1d_trklen_Pos_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_midpi=new TH1D(Form("h1d_trklen_Pos_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_midp=new TH1D(Form("h1d_trklen_Pos_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_midmu=new TH1D(Form("h1d_trklen_Pos_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_mideg=new TH1D(Form("h1d_trklen_Pos_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_midother=new TH1D(Form("h1d_trklen_Pos_midother"), Form("midother"), n_b, b_min, b_max);

	TH1D *h1d_trklen_Pos_bmrw=new TH1D(Form("h1d_trklen_Pos_bmrw"), Form("MC default"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_bmrw_el=new TH1D(Form("h1d_trklen_Pos_bmrw_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_bmrw_inel=new TH1D(Form("h1d_trklen_Pos_bmrw_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_bmrw_midcosmic=new TH1D(Form("h1d_trklen_Pos_bmrw_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_bmrw_midpi=new TH1D(Form("h1d_trklen_Pos_bmrw_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_bmrw_midp=new TH1D(Form("h1d_trklen_Pos_bmrw_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_bmrw_midmu=new TH1D(Form("h1d_trklen_Pos_bmrw_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_bmrw_mideg=new TH1D(Form("h1d_trklen_Pos_bmrw_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_bmrw_midother=new TH1D(Form("h1d_trklen_Pos_bmrw_midother"), Form("midother"), n_b, b_min, b_max);


	TH1D *h1d_trklen_BQ=new TH1D(Form("h1d_trklen_BQ"), Form("MC default"), n_b, b_min, b_max);
	TH1D *h1d_trklen_BQ_el=new TH1D(Form("h1d_trklen_BQ_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_BQ_inel=new TH1D(Form("h1d_trklen_BQ_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_BQ_midcosmic=new TH1D(Form("h1d_trklen_BQ_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_BQ_midpi=new TH1D(Form("h1d_trklen_BQ_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_BQ_midp=new TH1D(Form("h1d_trklen_BQ_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_BQ_midmu=new TH1D(Form("h1d_trklen_BQ_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_BQ_mideg=new TH1D(Form("h1d_trklen_BQ_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_BQ_midother=new TH1D(Form("h1d_trklen_BQ_midother"), Form("midother"), n_b, b_min, b_max);

	TH1D *h1d_trklen_BQ_bmrw=new TH1D(Form("h1d_trklen_BQ_bmrw"), Form("MC weighted"), n_b, b_min, b_max);
	TH1D *h1d_trklen_BQ_bmrw_el=new TH1D(Form("h1d_trklen_BQ_bmrw_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_BQ_bmrw_inel=new TH1D(Form("h1d_trklen_BQ_bmrw_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_BQ_bmrw_midcosmic=new TH1D(Form("h1d_trklen_BQ_bmrw_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_BQ_bmrw_midpi=new TH1D(Form("h1d_trklen_BQ_bmrw_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_BQ_bmrw_midp=new TH1D(Form("h1d_trklen_BQ_bmrw_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_BQ_bmrw_midmu=new TH1D(Form("h1d_trklen_BQ_bmrw_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_BQ_bmrw_mideg=new TH1D(Form("h1d_trklen_BQ_bmrw_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_BQ_bmrw_midother=new TH1D(Form("h1d_trklen_BQ_bmrw_midother"), Form("midother"), n_b, b_min, b_max);


	TH1D *h1d_trklen_RecoInel=new TH1D(Form("h1d_trklen_RecoInel"), Form("MC default"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_el=new TH1D(Form("h1d_trklen_RecoInel_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_inel=new TH1D(Form("h1d_trklen_RecoInel_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_midcosmic=new TH1D(Form("h1d_trklen_RecoInel_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_midpi=new TH1D(Form("h1d_trklen_RecoInel_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_midp=new TH1D(Form("h1d_trklen_RecoInel_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_midmu=new TH1D(Form("h1d_trklen_RecoInel_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_mideg=new TH1D(Form("h1d_trklen_RecoInel_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_midother=new TH1D(Form("h1d_trklen_RecoInel_midother"), Form("midother"), n_b, b_min, b_max);

	TH1D *h1d_trklen_RecoInel_bmrw=new TH1D(Form("h1d_trklen_RecoInel_bmrw"), Form("MC weighted"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_bmrw_el=new TH1D(Form("h1d_trklen_RecoInel_bmrw_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_bmrw_inel=new TH1D(Form("h1d_trklen_RecoInel_bmrw_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_bmrw_midcosmic=new TH1D(Form("h1d_trklen_RecoInel_bmrw_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_bmrw_midpi=new TH1D(Form("h1d_trklen_RecoInel_bmrw_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_bmrw_midp=new TH1D(Form("h1d_trklen_RecoInel_bmrw_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_bmrw_midmu=new TH1D(Form("h1d_trklen_RecoInel_bmrw_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_bmrw_mideg=new TH1D(Form("h1d_trklen_RecoInel_bmrw_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_RecoInel_bmrw_midother=new TH1D(Form("h1d_trklen_RecoInel_bmrw_midother"), Form("midother"), n_b, b_min, b_max);

	//MisID:P-rich
	TH1D *h1d_trklen_MidpRich=new TH1D(Form("h1d_trklen_MidpRich"), Form("MC default"), n_b, b_min, b_max);
	TH1D *h1d_trklen_MidpRich_el=new TH1D(Form("h1d_trklen_MidpRich_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_MidpRich_inel=new TH1D(Form("h1d_trklen_MidpRich_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_MidpRich_midcosmic=new TH1D(Form("h1d_trklen_MidpRich_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_MidpRich_midpi=new TH1D(Form("h1d_trklen_MidpRich_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_MidpRich_midp=new TH1D(Form("h1d_trklen_MidpRich_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_MidpRich_midmu=new TH1D(Form("h1d_trklen_MidpRich_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_MidpRich_mideg=new TH1D(Form("h1d_trklen_MidpRich_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_MidpRich_midother=new TH1D(Form("h1d_trklen_MidpRich_midother"), Form("midother"), n_b, b_min, b_max);

	TH1D *h1d_trklen_MidpRich_bmrw=new TH1D(Form("h1d_trklen_MidpRich_bmrw"), Form("MC weighted"), n_b, b_min, b_max);
	TH1D *h1d_trklen_MidpRich_bmrw_el=new TH1D(Form("h1d_trklen_MidpRich_bmrw_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_MidpRich_bmrw_inel=new TH1D(Form("h1d_trklen_MidpRich_bmrw_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_MidpRich_bmrw_midcosmic=new TH1D(Form("h1d_trklen_MidpRich_bmrw_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_MidpRich_bmrw_midpi=new TH1D(Form("h1d_trklen_MidpRich_bmrw_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_MidpRich_bmrw_midp=new TH1D(Form("h1d_trklen_MidpRich_bmrw_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_MidpRich_bmrw_midmu=new TH1D(Form("h1d_trklen_MidpRich_bmrw_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_MidpRich_bmrw_mideg=new TH1D(Form("h1d_trklen_MidpRich_bmrw_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_MidpRich_bmrw_midother=new TH1D(Form("h1d_trklen_MidpRich_bmrw_midother"), Form("midother"), n_b, b_min, b_max);

	//El-rich
	//bq
	TH1D *h1d_trklen_ElRich=new TH1D(Form("h1d_trklen_ElRich"), Form("MC default"), n_b, b_min, b_max);
	TH1D *h1d_trklen_ElRich_el=new TH1D(Form("h1d_trklen_ElRich_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_ElRich_inel=new TH1D(Form("h1d_trklen_ElRich_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_ElRich_midcosmic=new TH1D(Form("h1d_trklen_ElRich_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_ElRich_midpi=new TH1D(Form("h1d_trklen_ElRich_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_ElRich_midp=new TH1D(Form("h1d_trklen_ElRich_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_ElRich_midmu=new TH1D(Form("h1d_trklen_ElRich_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_ElRich_mideg=new TH1D(Form("h1d_trklen_ElRich_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_ElRich_midother=new TH1D(Form("h1d_trklen_ElRich_midother"), Form("midother"), n_b, b_min, b_max);

	TH1D *h1d_trklen_ElRich_bmrw=new TH1D(Form("h1d_trklen_ElRich_bmrw"), Form("MC weighted"), n_b, b_min, b_max);
	TH1D *h1d_trklen_ElRich_bmrw_el=new TH1D(Form("h1d_trklen_ElRich_bmrw_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_ElRich_bmrw_inel=new TH1D(Form("h1d_trklen_ElRich_bmrw_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_ElRich_bmrw_midcosmic=new TH1D(Form("h1d_trklen_ElRich_bmrw_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_ElRich_bmrw_midpi=new TH1D(Form("h1d_trklen_ElRich_bmrw_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_ElRich_bmrw_midp=new TH1D(Form("h1d_trklen_ElRich_bmrw_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_ElRich_bmrw_midmu=new TH1D(Form("h1d_trklen_ElRich_bmrw_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_ElRich_bmrw_mideg=new TH1D(Form("h1d_trklen_ElRich_bmrw_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_ElRich_bmrw_midother=new TH1D(Form("h1d_trklen_ElRich_bmrw_midother"), Form("midother"), n_b, b_min, b_max);

	//calosz
	TH1D *h1d_trklen_CaloSz_ElRich=new TH1D(Form("h1d_trklen_CaloSz_ElRich"), Form("MC default"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_ElRich_el=new TH1D(Form("h1d_trklen_CaloSz_ElRich_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_ElRich_inel=new TH1D(Form("h1d_trklen_CaloSz_ElRich_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_ElRich_midcosmic=new TH1D(Form("h1d_trklen_CaloSz_ElRich_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_ElRich_midpi=new TH1D(Form("h1d_trklen_CaloSz_ElRich_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_ElRich_midp=new TH1D(Form("h1d_trklen_CaloSz_ElRich_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_ElRich_midmu=new TH1D(Form("h1d_trklen_CaloSz_ElRich_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_ElRich_mideg=new TH1D(Form("h1d_trklen_CaloSz_ElRich_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_ElRich_midother=new TH1D(Form("h1d_trklen_CaloSz_ElRich_midother"), Form("midother"), n_b, b_min, b_max);

	TH1D *h1d_trklen_CaloSz_ElRich_bmrw=new TH1D(Form("h1d_trklen_CaloSz_ElRich_bmrw"), Form("MC weighted"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_ElRich_bmrw_el=new TH1D(Form("h1d_trklen_CaloSz_ElRich_bmrw_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_ElRich_bmrw_inel=new TH1D(Form("h1d_trklen_CaloSz_ElRich_bmrw_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_ElRich_bmrw_midcosmic=new TH1D(Form("h1d_trklen_CaloSz_ElRich_bmrw_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_ElRich_bmrw_midpi=new TH1D(Form("h1d_trklen_CaloSz_ElRich_bmrw_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_ElRich_bmrw_midp=new TH1D(Form("h1d_trklen_CaloSz_ElRich_bmrw_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_ElRich_bmrw_midmu=new TH1D(Form("h1d_trklen_CaloSz_ElRich_bmrw_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_ElRich_bmrw_mideg=new TH1D(Form("h1d_trklen_CaloSz_ElRich_bmrw_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_CaloSz_ElRich_bmrw_midother=new TH1D(Form("h1d_trklen_CaloSz_ElRich_bmrw_midother"), Form("midother"), n_b, b_min, b_max);

	//pos
	TH1D *h1d_trklen_Pos_ElRich=new TH1D(Form("h1d_trklen_Pos_ElRich"), Form("MC default"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_ElRich_el=new TH1D(Form("h1d_trklen_Pos_ElRich_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_ElRich_inel=new TH1D(Form("h1d_trklen_Pos_ElRich_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_ElRich_midcosmic=new TH1D(Form("h1d_trklen_Pos_ElRich_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_ElRich_midpi=new TH1D(Form("h1d_trklen_Pos_ElRich_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_ElRich_midp=new TH1D(Form("h1d_trklen_Pos_ElRich_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_ElRich_midmu=new TH1D(Form("h1d_trklen_Pos_ElRich_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_ElRich_mideg=new TH1D(Form("h1d_trklen_Pos_ElRich_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_ElRich_midother=new TH1D(Form("h1d_trklen_Pos_ElRich_midother"), Form("midother"), n_b, b_min, b_max);

	TH1D *h1d_trklen_Pos_ElRich_bmrw=new TH1D(Form("h1d_trklen_Pos_ElRich_bmrw"), Form("MC weighted"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_ElRich_bmrw_el=new TH1D(Form("h1d_trklen_Pos_ElRich_bmrw_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_ElRich_bmrw_inel=new TH1D(Form("h1d_trklen_Pos_ElRich_bmrw_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_ElRich_bmrw_midcosmic=new TH1D(Form("h1d_trklen_Pos_ElRich_bmrw_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_ElRich_bmrw_midpi=new TH1D(Form("h1d_trklen_Pos_ElRich_bmrw_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_ElRich_bmrw_midp=new TH1D(Form("h1d_trklen_Pos_ElRich_bmrw_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_ElRich_bmrw_midmu=new TH1D(Form("h1d_trklen_Pos_ElRich_bmrw_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_ElRich_bmrw_mideg=new TH1D(Form("h1d_trklen_Pos_ElRich_bmrw_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_trklen_Pos_ElRich_bmrw_midother=new TH1D(Form("h1d_trklen_Pos_ElRich_bmrw_midother"), Form("midother"), n_b, b_min, b_max);



	//endz
/*
	TH1D *h1d_endz_CaloSz=new TH1D(Form("h1d_endz_CaloSz"), Form("MC CaloSz"), n_b, b_min, b_max);
	TH1D *h1d_endz_CaloSz_el=new TH1D(Form("h1d_endz_CaloSz_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_endz_CaloSz_inel=new TH1D(Form("h1d_endz_CaloSz_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_endz_CaloSz_midcosmic=new TH1D(Form("h1d_endz_CaloSz_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_endz_CaloSz_midpi=new TH1D(Form("h1d_endz_CaloSz_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_endz_CaloSz_midp=new TH1D(Form("h1d_endz_CaloSz_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_endz_CaloSz_midmu=new TH1D(Form("h1d_endz_CaloSz_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_endz_CaloSz_mideg=new TH1D(Form("h1d_endz_CaloSz_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_endz_CaloSz_midother=new TH1D(Form("h1d_endz_CaloSz_midother"), Form("midother"), n_b, b_min, b_max);

	TH1D *h1d_endz_Pos=new TH1D(Form("h1d_endz_Pos"), Form("MC Pos"), n_b, b_min, b_max);
	TH1D *h1d_endz_Pos_el=new TH1D(Form("h1d_endz_Pos_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_endz_Pos_inel=new TH1D(Form("h1d_endz_Pos_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_endz_Pos_midcosmic=new TH1D(Form("h1d_endz_Pos_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_endz_Pos_midpi=new TH1D(Form("h1d_endz_Pos_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_endz_Pos_midp=new TH1D(Form("h1d_endz_Pos_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_endz_Pos_midmu=new TH1D(Form("h1d_endz_Pos_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_endz_Pos_mideg=new TH1D(Form("h1d_endz_Pos_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_endz_Pos_midother=new TH1D(Form("h1d_endz_Pos_midother"), Form("midother"), n_b, b_min, b_max);

	TH1D *h1d_endz_BQ=new TH1D(Form("h1d_endz_BQ"), Form("MC BQ"), n_b, b_min, b_max);
	TH1D *h1d_endz_BQ_el=new TH1D(Form("h1d_endz_BQ_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_endz_BQ_inel=new TH1D(Form("h1d_endz_BQ_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_endz_BQ_midcosmic=new TH1D(Form("h1d_endz_BQ_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_endz_BQ_midpi=new TH1D(Form("h1d_endz_BQ_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_endz_BQ_midp=new TH1D(Form("h1d_endz_BQ_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_endz_BQ_midmu=new TH1D(Form("h1d_endz_BQ_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_endz_BQ_mideg=new TH1D(Form("h1d_endz_BQ_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_endz_BQ_midother=new TH1D(Form("h1d_endz_BQ_midother"), Form("midother"), n_b, b_min, b_max);


	//trklen vs endz
	TH2D *h2d_trklen_endz_CaloSz=new TH2D(Form("h2d_trklen_endz_CaloSz"), Form("MC CaloSz"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_CaloSz_el=new TH2D(Form("h2d_trklen_endz_CaloSz_el"), Form("el"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_CaloSz_inel=new TH2D(Form("h2d_trklen_endz_CaloSz_inel"), Form("inel"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_CaloSz_midcosmic=new TH2D(Form("h2d_trklen_endz_CaloSz_midcosmic"), Form("midcosmic"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_CaloSz_midpi=new TH2D(Form("h2d_trklen_endz_CaloSz_midpi"), Form("midpi"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_CaloSz_midp=new TH2D(Form("h2d_trklen_endz_CaloSz_midp"), Form("midp"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_CaloSz_midmu=new TH2D(Form("h2d_trklen_endz_CaloSz_midmu"), Form("midmu"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_CaloSz_mideg=new TH2D(Form("h2d_trklen_endz_CaloSz_mideg"), Form("mideg"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_CaloSz_midother=new TH2D(Form("h2d_trklen_endz_CaloSz_midother"), Form("midother"), n_b, b_min, b_max, n_b, b_min, b_max);

	TH2D *h2d_trklen_endz_Pos=new TH2D(Form("h2d_trklen_endz_Pos"), Form("MC Pos"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_Pos_el=new TH2D(Form("h2d_trklen_endz_Pos_el"), Form("el"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_Pos_inel=new TH2D(Form("h2d_trklen_endz_Pos_inel"), Form("inel"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_Pos_midcosmic=new TH2D(Form("h2d_trklen_endz_Pos_midcosmic"), Form("midcosmic"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_Pos_midpi=new TH2D(Form("h2d_trklen_endz_Pos_midpi"), Form("midpi"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_Pos_midp=new TH2D(Form("h2d_trklen_endz_Pos_midp"), Form("midp"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_Pos_midmu=new TH2D(Form("h2d_trklen_endz_Pos_midmu"), Form("midmu"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_Pos_mideg=new TH2D(Form("h2d_trklen_endz_Pos_mideg"), Form("mideg"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_Pos_midother=new TH2D(Form("h2d_trklen_endz_Pos_midother"), Form("midother"), n_b, b_min, b_max, n_b, b_min, b_max);

	TH2D *h2d_trklen_endz_BQ=new TH2D(Form("h2d_trklen_endz_BQ"), Form("MC BQ"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_BQ_el=new TH2D(Form("h2d_trklen_endz_BQ_el"), Form("el"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_BQ_inel=new TH2D(Form("h2d_trklen_endz_BQ_inel"), Form("inel"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_BQ_midcosmic=new TH2D(Form("h2d_trklen_endz_BQ_midcosmic"), Form("midcosmic"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_BQ_midpi=new TH2D(Form("h2d_trklen_endz_BQ_midpi"), Form("midpi"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_BQ_midp=new TH2D(Form("h2d_trklen_endz_BQ_midp"), Form("midp"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_BQ_midmu=new TH2D(Form("h2d_trklen_endz_BQ_midmu"), Form("midmu"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_BQ_mideg=new TH2D(Form("h2d_trklen_endz_BQ_mideg"), Form("mideg"), n_b, b_min, b_max, n_b, b_min, b_max);
	TH2D *h2d_trklen_endz_BQ_midother=new TH2D(Form("h2d_trklen_endz_BQ_midother"), Form("midother"), n_b, b_min, b_max, n_b, b_min, b_max);
*/

	//true len
	TH1D *h1d_truetrklen_CaloSz=new TH1D(Form("h1d_truetrklen_CaloSz"), Form("MC CaloSz"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_CaloSz_el=new TH1D(Form("h1d_truetrklen_CaloSz_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_CaloSz_inel=new TH1D(Form("h1d_truetrklen_CaloSz_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_CaloSz_midcosmic=new TH1D(Form("h1d_truetrklen_CaloSz_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_CaloSz_midpi=new TH1D(Form("h1d_truetrklen_CaloSz_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_CaloSz_midp=new TH1D(Form("h1d_truetrklen_CaloSz_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_CaloSz_midmu=new TH1D(Form("h1d_truetrklen_CaloSz_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_CaloSz_mideg=new TH1D(Form("h1d_truetrklen_CaloSz_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_CaloSz_midother=new TH1D(Form("h1d_truetrklen_CaloSz_midother"), Form("midother"), n_b, b_min, b_max);

	TH1D *h1d_truetrklen_Pos=new TH1D(Form("h1d_truetrklen_Pos"), Form("MC Pos"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_Pos_el=new TH1D(Form("h1d_truetrklen_Pos_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_Pos_inel=new TH1D(Form("h1d_truetrklen_Pos_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_Pos_midcosmic=new TH1D(Form("h1d_truetrklen_Pos_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_Pos_midpi=new TH1D(Form("h1d_truetrklen_Pos_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_Pos_midp=new TH1D(Form("h1d_truetrklen_Pos_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_Pos_midmu=new TH1D(Form("h1d_truetrklen_Pos_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_Pos_mideg=new TH1D(Form("h1d_truetrklen_Pos_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_Pos_midother=new TH1D(Form("h1d_truetrklen_Pos_midother"), Form("midother"), n_b, b_min, b_max);

	TH1D *h1d_truetrklen_BQ=new TH1D(Form("h1d_truetrklen_BQ"), Form("MC BQ"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_BQ_el=new TH1D(Form("h1d_truetrklen_BQ_el"), Form("el"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_BQ_inel=new TH1D(Form("h1d_truetrklen_BQ_inel"), Form("inel"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_BQ_midcosmic=new TH1D(Form("h1d_truetrklen_BQ_midcosmic"), Form("midcosmic"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_BQ_midpi=new TH1D(Form("h1d_truetrklen_BQ_midpi"), Form("midpi"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_BQ_midp=new TH1D(Form("h1d_truetrklen_BQ_midp"), Form("midp"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_BQ_midmu=new TH1D(Form("h1d_truetrklen_BQ_midmu"), Form("midmu"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_BQ_mideg=new TH1D(Form("h1d_truetrklen_BQ_mideg"), Form("mideg"), n_b, b_min, b_max);
	TH1D *h1d_truetrklen_BQ_midother=new TH1D(Form("h1d_truetrklen_BQ_midother"), Form("midother"), n_b, b_min, b_max);


	//true endz
	//TH1D *h1d_trueendz_CaloSz=new TH1D(Form("h1d_trueendz_CaloSz"), Form("MC CaloSz"), n_b, b_min, b_max);
	//TH1D *h1d_trueendz_Pos=new TH1D(Form("h1d_trueendz_Pos"), Form("MC Pos"), n_b, b_min, b_max);
	//TH1D *h1d_trueendz_BQ=new TH1D(Form("h1d_trueendz_BQ"), Form("MC BQ"), n_b, b_min, b_max);

	//TH1D *h1d_trueendzSCE_CaloSz=new TH1D(Form("h1d_trueendzSCE_CaloSz"), Form("MC CaloSz"), n_b, b_min, b_max);
	//TH1D *h1d_trueendzSCE_Pos=new TH1D(Form("h1d_trueendzSCE_Pos"), Form("MC Pos"), n_b, b_min, b_max);
	//TH1D *h1d_trueendzSCE_BQ=new TH1D(Form("h1d_trueendzSCE_BQ"), Form("MC BQ"), n_b, b_min, b_max);

	//------------------------------------------------------------------------------------------------------------//

        //MC Beam Mom Gaussian 
        double m1=1007.1482; //MC prod4a [spec]
        double s1=60.703307; //MC prod4a [spec]

        //momentum cut range    
        double mu_min=m1-3.*s1;
        double mu_max=m1+3.*s1;



	for (Long64_t jentry=0; jentry<nentries;jentry++) { //main entry loop
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		isTestSample = true;
		if (ientry%2 == 0) isTestSample = false; //Divide MC sample by 2 parts: test+ufold

		//only select protons	
		if (beamtrackPdg!=pdg) continue; //only interested in protons

		//Event Selection Cut -- Part 1 ----------------------------------//
		bool IsBeamMatch=false; //if recostructed the right track (recoID=truthID)
		bool IsPandoraSlice=false; //pandora slice cut (can pandora reconstruct this track)
		bool IsCaloSize=false; //if calo size not empty
		bool IsIntersection=false; //if any track intersect with our reco track		
		if (primary_truth_Isbeammatched==1) IsBeamMatch=true;
		if (isprimarytrack==1&&isprimaryshower==0) IsPandoraSlice=true; 
		if (!primtrk_hitz->empty()) IsCaloSize=true;
		if (timeintersection->size()) IsIntersection=true;
		//----------------------------------------------------------------//

		//Truth label of Primarytrack_End ------------------------------------------------------------------------------------------------//
		bool IsPureInEL=false; //inel
		bool IsPureEL=false; //el

		if (strcmp(primary_truth_EndProcess->c_str(),"protonInelastic")==0) { IsPureInEL=true; }
		else { IsPureEL=true; }
		//--------------------------------------------------------------------------------------------------------------------------------//

		//Get true start/end point -----------------------------------------------------------------------//
		double true_endz=primary_truth_EndPosition_MC[2]; 
		double true_endy=primary_truth_EndPosition_MC[1]; 
		double true_endx=primary_truth_EndPosition_MC[0];

		double true_stz=primary_truth_StartPosition_MC[2];
		double true_sty=primary_truth_StartPosition_MC[1];
		double true_stx=primary_truth_StartPosition_MC[0];

		double offset_z = 0;
/*
		double true_endz_sce=true_endz;


		TVector3 true_Endpoint(true_endx, true_endy, true_endz);
     		if (true_Endpoint.X()>0){
       			//true_endz += RecoFwd_Displacement_Z_Pos->GetBinContent(RecoFwd_Displacement_Z_Pos->FindBin(true_beam_endX, true_beam_endY, true_beam_endZ));
       			offset_z = InterpolateSplines(RecoFwd_Displacement_Z_Pos, true_Endpoint.X(), true_Endpoint.Y(), true_Endpoint.Z(), 3, 1, 2);
     		}
     		else{
       			//true_endz += RecoFwd_Displacement_Z_Neg->GetBinContent(RecoFwd_Displacement_Z_Neg->FindBin(true_beam_endX, true_beam_endY, true_beam_endZ));
       			offset_z = InterpolateSplines(RecoFwd_Displacement_Z_Neg, true_Endpoint.X(), true_Endpoint.Y(), true_Endpoint.Z(), 3, 1, 1);
     		}
     		//std::cout<<true_beam_endX<<" "<<true_beam_endY<<" "<<true_beam_endZ<<" "<<offset_z<<std::endl;
     		true_endz_sce += offset_z;
*/

		bool IsTrueEndOutside=false;
		if (true_endz<0.) IsTrueEndOutside=true;

		//Get reco info ----------------------------------------------------------------------------------//
		//Evt Classification -----------------------------------------------------------------------------//
		//signal -----------------------------//
		bool kinel=false;
		bool kel=false;
		//bool kmcs=false;
		if (IsBeamMatch) { //beam-match
			if (IsPureInEL) kinel=true;
			if (IsPureEL) kel=true;
			//if (IsPureMCS) kmcs=true;
		} //beam-match

		//background ------------------------------------------------------------------------//
		bool kMIDcosmic=false; //beam or cosmic
		bool kMIDpi=false; //+-pi
		bool kMIDp=false; //p
		bool kMIDmu=false; //mu
		bool kMIDeg=false; //e/gamma
		bool kMIDother=false; //other
		if (!IsBeamMatch) { //!beam-match
			if (primary_truth_byE_origin==2) { 
				kMIDcosmic=true;
			}
			else if (std::abs(primary_truth_byE_PDG)==211) {
				kMIDpi=true;
			}
			else if (primary_truth_byE_PDG==2212) {
				kMIDp=true;
			}
			else if (std::abs(primary_truth_byE_PDG)==13) {
				kMIDmu=true;
			}
			else if (std::abs(primary_truth_byE_PDG)==11 || primary_truth_byE_PDG==22) {
				kMIDeg=true;
			}
			else {
				kMIDother=true;
			}
		} //!beam-match	
		//cout<<"kMIDcosmic:"<<kMIDcosmic<<endl;
		//Evt Classification -----------------------------------------------------------------------------//

		//reco pos info & cut
		double reco_stx=-99, reco_sty=-99, reco_stz=-99;
		double reco_endx=-99, reco_endy=-99, reco_endz=-99;
		//double dx_reco_stx__bposx_ff=-999, dy_reco_sty__bposy_ff=-999, dz_reco_stz__bposz_ff=-999;
		bool IsPos=false;
		if (IsCaloSize) {
			reco_stx=primtrk_hitx->at(0); 
			reco_sty=primtrk_hity->at(0);
			reco_stz=primtrk_hitz->at(0);

			reco_endx=primtrk_hitx->at(primtrk_dedx->size()-1);	
			reco_endy=primtrk_hity->at(primtrk_dedx->size()-1);
			reco_endz=primtrk_hitz->at(primtrk_dedx->size()-1);

			if (reco_stz>reco_endz) {
				reco_endx=primtrk_hitx->at(0); 
				reco_endy=primtrk_hity->at(0);
				reco_endz=primtrk_hitz->at(0);

				reco_stx=primtrk_hitx->at(primtrk_dedx->size()-1);	
				reco_sty=primtrk_hity->at(primtrk_dedx->size()-1);
				reco_stz=primtrk_hitz->at(primtrk_dedx->size()-1);
			}

			double beam_dx=(reco_stx-mean_StartX)/sigma_StartX;
			double beam_dy=(reco_sty-mean_StartY)/sigma_StartY;
			double beam_dz=(reco_stz-mean_StartZ)/sigma_StartZ;
			double beam_dxy=sqrt(pow(beam_dx,2)+pow(beam_dy,2));	

			if (beam_dx>=dx_min&&beam_dx<=dx_max) { //dx
				if (beam_dy>=dy_min&&beam_dy<=dy_max) { //dy
					if (beam_dz>=dz_min&&beam_dz<=dz_max) { //dz
						if (beam_dxy>=dxy_min&&beam_dxy<=dxy_max) { //dxy
							IsPos=true;
						} //dxy
					} //dz
				} //dy
			} //dx

		}

		//cosine_theta/cut
		bool IsCosine=false;
		double cosine_beam_spec_primtrk=-999; 
		//cosine_beam_spec_primtrk=beamDirx_spec->at(0)*primaryStartDirection[0]+beamDiry_spec->at(0)*primaryStartDirection[1]+beamDirz_spec->at(0)*primaryStartDirection[2]; //cosine between beam_spec and primary trk direction(no SCE corr.)
		TVector3 dir;
		if (IsCaloSize) {	
			//trk direction after SCE corr.
			TVector3 pt0(primtrk_hitx->at(0), primtrk_hity->at(0), primtrk_hitz->at(0));
			TVector3 pt1(primtrk_hitx->at(-1+primtrk_hitx->size()), primtrk_hity->at(-1+primtrk_hity->size()), primtrk_hitz->at(-1+primtrk_hitz->size()));
			//TVector3 dir = pt1 - pt0;
			dir = pt1 - pt0;
			dir = dir.Unit();

			//beam direction
      			//TVector3 beamdir(cos(beam_angleX_mc*TMath::Pi()/180), cos(beam_angleY_mc*TMath::Pi()/180), cos(beam_angleZ_mc*TMath::Pi()/180));
      			TVector3 beamdir(beamDirx_spec->at(0),beamDiry_spec->at(0),beamDirz_spec->at(0));
      			beamdir = beamdir.Unit();
      			//beam_costh = dir.Dot(beamdir);
      			cosine_beam_spec_primtrk=dir.Dot(beamdir);
		}

		if (cosine_beam_spec_primtrk<0) { cosine_beam_spec_primtrk=-1.*cosine_beam_spec_primtrk; }
		if (cosine_beam_spec_primtrk>cosine_beam_primtrk_min) { IsCosine=true; }

		//xy-cut
		bool IsXY=false;		
		double reco_stx_noSCE=0, reco_sty_noSCE=0, reco_stz_noSCE=0; //start-pos, before sce
		if (primaryEndPosition[2]>primaryStartPosition[2]) { //check if Pandora flip the sign
			reco_stx_noSCE=primaryStartPosition[0];
			reco_sty_noSCE=primaryStartPosition[1];
			reco_stz_noSCE=primaryStartPosition[2];
		} //check if Pandora flip the sign
		else {
			reco_stx_noSCE=primaryEndPosition[0];
			reco_sty_noSCE=primaryEndPosition[1];
			reco_stz_noSCE=primaryEndPosition[2];
		}
		if ((pow(((reco_stx_noSCE-mean_x)/dev_x),2)+pow(((reco_sty_noSCE-mean_y)/dev_y),2))<=1.) IsXY=true;

		//beam quality cut
		bool IsBQ=false;
		if (IsCosine&&IsPos) IsBQ=true;

		int index_reco_endz=0;
		double wid_reco_max=-9999;
		double range_reco=-999;
		vector<double> reco_trklen_accum;
  		reco_trklen_accum.reserve(primtrk_hitz->size());
		double kereco_calo=0;
		double kereco_range=0;
		double kereco_range2=0;
		vector<double> EDept;
		double pid=-99;
		if (IsCaloSize) { //if calo size not empty
		  vector<double> trkdedx;
		  vector<double> trkres;
		  for (size_t h=0; h<primtrk_dedx->size(); ++h) { //loop over reco hits of a given track
			double hitx_reco=primtrk_hitx->at(h);
			double hity_reco=primtrk_hity->at(h);
			double hitz_reco=primtrk_hitz->at(h);
			double resrange_reco=primtrk_resrange->at(h);

			double dqdx=primtrk_dqdx->at(h);
			double pitch=primtrk_pitch->at(h);

			int wid_reco=primtrk_wid->at(-1+primtrk_wid->size()-h);
			double pt_reco=primtrk_pt->at(-1+primtrk_wid->size()-h);

			if (wid_reco>wid_reco_max) { 
				wid_reco_max=wid_reco;
				index_reco_endz=(int)-1+primtrk_wid->size()-h;
			}

			double cali_dedx=0.;
			cali_dedx=dedx_function_35ms(dqdx, hitx_reco, hity_reco, hitz_reco);
			EDept.push_back(cali_dedx*pitch);

			if (h==1) range_reco=0;
			if (h>=1) {
    					range_reco += sqrt( pow(primtrk_hitx->at(h)-primtrk_hitx->at(h-1), 2)+
					    		    pow(primtrk_hity->at(h)-primtrk_hity->at(h-1), 2)+
					    		    pow(primtrk_hitz->at(h)-primtrk_hitz->at(h-1), 2) );
					reco_trklen_accum[h] = range_reco;
			}

			kereco_calo+=cali_dedx*pitch;
			kereco_range+=pitch*dedx_predict(resrange_reco);
			kereco_range2+=pitch*(double)gr_predict_dedx_resrange->Eval(resrange_reco);

			trkdedx.push_back(cali_dedx);
			trkres.push_back(resrange_reco);

		  } //loop over reco hits of a given track
		  //range_reco=primtrk_range->at(0);

		  pid=chi2pid(trkdedx,trkres); //pid using stopping proton hypothesis
		} //if calo size not empty


		//Reco stopping/Inel p cut
		bool IsRecoStop=false;
		bool IsRecoInEL=false;
		double mom_beam_spec=-99; mom_beam_spec=beamMomentum_spec->at(0);
		//double range_reco=-99; if (!primtrk_range->empty()) range_reco=primtrk_range->at(0); //reco primary trklen
		double csda_val_spec=csda_range_vs_mom_sm->Eval(mom_beam_spec);

		//if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) IsRecoStop=true;
		//if ((range_reco/csda_val_spec)<min_norm_trklen_csda) IsRecoInEL=true;

		if ((range_reco/csda_val_spec)<min_norm_trklen_csda) { //inel region
			if (pid>pid_1) IsRecoInEL=true; 
			if (pid<=pid_1) IsRecoStop=true; 
		} //inel region
		if ((range_reco/csda_val_spec)>=min_norm_trklen_csda&&(range_reco/csda_val_spec)<max_norm_trklen_csda) { //stopping p region
			if (pid>pid_2) IsRecoInEL=true; 
			if (pid<=pid_2) IsRecoStop=true;
		} //stopping p region


		//kinetic energies
		double ke_beam_spec=p2ke(mom_beam_spec); //ke_beam_spec [GeV]
		double ke_beam_spec_MeV=1000.*ke_beam_spec; //ke_beam_spec [MeV]
		double ke_trklen=ke_vs_csda_range_sm->Eval(range_reco); //[unit: GeV]
		double ke_trklen_MeV=1000.*ke_trklen; //[unit: MeV]
		double ke_calo_MeV=0;

		//Get true trklen ---------------------------------------------------------------------------------------//
		double range_true=-999;
		int key_st = 0;
		double tmp_z = 9999;
		vector<double> true_trklen_accum;
		//cout<<"beamtrk_z->size():"<<beamtrk_z->size()<<endl;
		true_trklen_accum.reserve(beamtrk_z->size()); // initialize true_trklen_accum
		//cout<<"ck0"<<endl;
		for (int iz=0; iz<(int)beamtrk_z->size(); iz++) {
			if (abs(beamtrk_z->at(iz)) < tmp_z){
				tmp_z = abs(beamtrk_z->at(iz));
				key_st = iz; // find the point where the beam enters the TPC (find the smallest abs(Z))
			}
			//cout<<"ck0/"<<endl;
			true_trklen_accum[iz] = 0.; // initialize true_trklen_accum
			//cout<<"ck0///"<<endl;
		}
		//cout<<"ck1"<<endl;
		for (int iz=key_st+1; iz<(int)beamtrk_z->size(); iz++){
			if (iz == key_st+1) range_true = 0;
			range_true += sqrt( pow(beamtrk_x->at(iz)-beamtrk_x->at(iz-1), 2)+
					pow(beamtrk_y->at(iz)-beamtrk_y->at(iz-1), 2)+	
					pow(beamtrk_z->at(iz)-beamtrk_z->at(iz-1), 2) );						    	
			true_trklen_accum[iz] = range_true;
		}




		if (IsPandoraSlice&&IsCaloSize) { //calosz cut
			//before bmrw
			Fill1DHist(h1d_trklen_CaloSz, range_reco);
			Fill1DHist(h1d_truetrklen_CaloSz, range_true);
			//Fill1DHist(h1d_trueendz_CaloSz, true_endz);
			//Fill1DHist(h1d_trueendzSCE_CaloSz, true_endz_sce);

			if (kinel) Fill1DHist(h1d_trklen_CaloSz_inel, range_reco); 
			if (kel) Fill1DHist(h1d_trklen_CaloSz_el, range_reco); 
			if (kMIDcosmic) Fill1DHist(h1d_trklen_CaloSz_midcosmic, range_reco); 
			if (kMIDpi) Fill1DHist(h1d_trklen_CaloSz_midpi, range_reco);
			if (kMIDp) Fill1DHist(h1d_trklen_CaloSz_midp, range_reco);
			if (kMIDmu) Fill1DHist(h1d_trklen_CaloSz_midmu, range_reco);
			if (kMIDeg) Fill1DHist(h1d_trklen_CaloSz_mideg, range_reco);
			if (kMIDother) Fill1DHist(h1d_trklen_CaloSz_midother, range_reco);

			if (kinel) Fill1DHist(h1d_truetrklen_CaloSz_inel, range_true); 
			if (kel) Fill1DHist(h1d_truetrklen_CaloSz_el, range_true); 
			if (kMIDcosmic) Fill1DHist(h1d_truetrklen_CaloSz_midcosmic, range_true); 
			if (kMIDpi) Fill1DHist(h1d_truetrklen_CaloSz_midpi, range_true);
			if (kMIDp) Fill1DHist(h1d_truetrklen_CaloSz_midp, range_true);
			if (kMIDmu) Fill1DHist(h1d_truetrklen_CaloSz_midmu, range_true);
			if (kMIDeg) Fill1DHist(h1d_truetrklen_CaloSz_mideg, range_true);
			if (kMIDother) Fill1DHist(h1d_truetrklen_CaloSz_midother, range_true);

/*
			//before bmrw
			//Fill1DHist(h1d_endz_CaloSz, reco_endz);

			//if (kinel) Fill1DHist(h1d_endz_CaloSz_inel, reco_endz); 
			//if (kel) Fill1DHist(h1d_endz_CaloSz_el, reco_endz); 
			//if (kMIDcosmic) Fill1DHist(h1d_endz_CaloSz_midcosmic, reco_endz); 
			//if (kMIDpi) Fill1DHist(h1d_endz_CaloSz_midpi, reco_endz);
			//if (kMIDp) Fill1DHist(h1d_endz_CaloSz_midp, reco_endz);
			//if (kMIDmu) Fill1DHist(h1d_endz_CaloSz_midmu, reco_endz);
			//if (kMIDeg) Fill1DHist(h1d_endz_CaloSz_mideg, reco_endz);
			//if (kMIDother) Fill1DHist(h1d_endz_CaloSz_midother, reco_endz);


			//before bmrw
			h2d_trklen_endz_CaloSz->Fill(range_reco, reco_endz);

			if (kinel) h2d_trklen_endz_CaloSz_inel->Fill(range_reco, reco_endz);
			if (kel) h2d_trklen_endz_CaloSz_el->Fill(range_reco, reco_endz); 
			if (kMIDcosmic) h2d_trklen_endz_CaloSz_midcosmic->Fill(range_reco, reco_endz);
			if (kMIDpi) h2d_trklen_endz_CaloSz_midpi->Fill(range_reco, reco_endz);
			if (kMIDp) h2d_trklen_endz_CaloSz_midp->Fill(range_reco, reco_endz);
			if (kMIDmu) h2d_trklen_endz_CaloSz_midmu->Fill(range_reco, reco_endz);
			if (kMIDeg) h2d_trklen_endz_CaloSz_mideg->Fill(range_reco, reco_endz);
			if (kMIDother) h2d_trklen_endz_CaloSz_midother->Fill(range_reco, reco_endz);
*/

			double mom_rw_minchi2=1.;
			if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) { //beam-mom cut (within 3-sigma)
				mom_rw_minchi2=bmrw_func->Eval(mom_beam_spec*1000.); //bmrw
				Fill1DWHist(h1d_trklen_CaloSz_bmrw, range_reco, mom_rw_minchi2);

				if (kinel) Fill1DWHist(h1d_trklen_CaloSz_bmrw_inel, range_reco, mom_rw_minchi2); 
				if (kel) Fill1DWHist(h1d_trklen_CaloSz_bmrw_el, range_reco, mom_rw_minchi2); 
				if (kMIDcosmic) Fill1DWHist(h1d_trklen_CaloSz_bmrw_midcosmic, range_reco, mom_rw_minchi2); 
				if (kMIDpi) Fill1DWHist(h1d_trklen_CaloSz_bmrw_midpi, range_reco, mom_rw_minchi2);
				if (kMIDp) Fill1DWHist(h1d_trklen_CaloSz_bmrw_midp, range_reco, mom_rw_minchi2);
				if (kMIDmu) Fill1DWHist(h1d_trklen_CaloSz_bmrw_midmu, range_reco, mom_rw_minchi2);
				if (kMIDeg) Fill1DWHist(h1d_trklen_CaloSz_bmrw_mideg, range_reco, mom_rw_minchi2);
				if (kMIDother) Fill1DWHist(h1d_trklen_CaloSz_bmrw_midother, range_reco, mom_rw_minchi2);

			} //beam-mom cut (within 3-sigma)
			else { //tail of beam
				Fill1DWHist(h1d_trklen_CaloSz_bmrw, range_reco, 1); //set weight to one

				if (kinel) Fill1DWHist(h1d_trklen_CaloSz_bmrw_inel, range_reco, 1); 
				if (kel) Fill1DWHist(h1d_trklen_CaloSz_bmrw_el, range_reco, 1); 
				if (kMIDcosmic) Fill1DWHist(h1d_trklen_CaloSz_bmrw_midcosmic, range_reco, 1); 
				if (kMIDpi) Fill1DWHist(h1d_trklen_CaloSz_bmrw_midpi, range_reco, 1);
				if (kMIDp) Fill1DWHist(h1d_trklen_CaloSz_bmrw_midp, range_reco, 1);
				if (kMIDmu) Fill1DWHist(h1d_trklen_CaloSz_bmrw_midmu, range_reco, 1);
				if (kMIDeg) Fill1DWHist(h1d_trklen_CaloSz_bmrw_mideg, range_reco, 1);
				if (kMIDother) Fill1DWHist(h1d_trklen_CaloSz_bmrw_midother, range_reco, 1);

			} //tail of beam			
		} //calosz cut

		if (IsPos&&IsPandoraSlice&&IsCaloSize) { //pos cut
			//before bmrw
			Fill1DHist(h1d_trklen_Pos, range_reco);
			Fill1DHist(h1d_truetrklen_Pos, range_true);
			//Fill1DHist(h1d_trueendz_Pos, true_endz);
			//Fill1DHist(h1d_trueendzSCE_Pos, true_endz_sce);

			if (kinel) Fill1DHist(h1d_trklen_Pos_inel, range_reco); 
			if (kel) Fill1DHist(h1d_trklen_Pos_el, range_reco); 
			if (kMIDcosmic) Fill1DHist(h1d_trklen_Pos_midcosmic, range_reco); 
			if (kMIDpi) Fill1DHist(h1d_trklen_Pos_midpi, range_reco);
			if (kMIDp) Fill1DHist(h1d_trklen_Pos_midp, range_reco);
			if (kMIDmu) Fill1DHist(h1d_trklen_Pos_midmu, range_reco);
			if (kMIDeg) Fill1DHist(h1d_trklen_Pos_mideg, range_reco);
			if (kMIDother) Fill1DHist(h1d_trklen_Pos_midother, range_reco);

			if (kinel) Fill1DHist(h1d_truetrklen_Pos_inel, range_true); 
			if (kel) Fill1DHist(h1d_truetrklen_Pos_el, range_true); 
			if (kMIDcosmic) Fill1DHist(h1d_truetrklen_Pos_midcosmic, range_true); 
			if (kMIDpi) Fill1DHist(h1d_truetrklen_Pos_midpi, range_true);
			if (kMIDp) Fill1DHist(h1d_truetrklen_Pos_midp, range_true);
			if (kMIDmu) Fill1DHist(h1d_truetrklen_Pos_midmu, range_true);
			if (kMIDeg) Fill1DHist(h1d_truetrklen_Pos_mideg, range_true);
			if (kMIDother) Fill1DHist(h1d_truetrklen_Pos_midother, range_true);

/*
			//before bmrw
			Fill1DHist(h1d_endz_Pos, reco_endz);

			if (kinel) Fill1DHist(h1d_endz_Pos_inel, reco_endz); 
			if (kel) Fill1DHist(h1d_endz_Pos_el, reco_endz); 
			if (kMIDcosmic) Fill1DHist(h1d_endz_Pos_midcosmic, reco_endz); 
			if (kMIDpi) Fill1DHist(h1d_endz_Pos_midpi, reco_endz);
			if (kMIDp) Fill1DHist(h1d_endz_Pos_midp, reco_endz);
			if (kMIDmu) Fill1DHist(h1d_endz_Pos_midmu, reco_endz);
			if (kMIDeg) Fill1DHist(h1d_endz_Pos_mideg, reco_endz);
			if (kMIDother) Fill1DHist(h1d_endz_Pos_midother, reco_endz);


			//before bmrw
			h2d_trklen_endz_Pos->Fill(range_reco, reco_endz);

			if (kinel) h2d_trklen_endz_Pos_inel->Fill(range_reco, reco_endz);
			if (kel) h2d_trklen_endz_Pos_el->Fill(range_reco, reco_endz); 
			if (kMIDcosmic) h2d_trklen_endz_Pos_midcosmic->Fill(range_reco, reco_endz);
			if (kMIDpi) h2d_trklen_endz_Pos_midpi->Fill(range_reco, reco_endz);
			if (kMIDp) h2d_trklen_endz_Pos_midp->Fill(range_reco, reco_endz);
			if (kMIDmu) h2d_trklen_endz_Pos_midmu->Fill(range_reco, reco_endz);
			if (kMIDeg) h2d_trklen_endz_Pos_mideg->Fill(range_reco, reco_endz);
			if (kMIDother) h2d_trklen_endz_Pos_midother->Fill(range_reco, reco_endz);
*/

			double mom_rw_minchi2=1.;
			if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) { //beam-mom cut (within 3-sigma)
				mom_rw_minchi2=bmrw_func->Eval(mom_beam_spec*1000.); //bmrw
				Fill1DWHist(h1d_trklen_Pos_bmrw, range_reco, mom_rw_minchi2);

				if (kinel) Fill1DWHist(h1d_trklen_Pos_bmrw_inel, range_reco, mom_rw_minchi2); 
				if (kel) Fill1DWHist(h1d_trklen_Pos_bmrw_el, range_reco, mom_rw_minchi2); 
				if (kMIDcosmic) Fill1DWHist(h1d_trklen_Pos_bmrw_midcosmic, range_reco, mom_rw_minchi2); 
				if (kMIDpi) Fill1DWHist(h1d_trklen_Pos_bmrw_midpi, range_reco, mom_rw_minchi2);
				if (kMIDp) Fill1DWHist(h1d_trklen_Pos_bmrw_midp, range_reco, mom_rw_minchi2);
				if (kMIDmu) Fill1DWHist(h1d_trklen_Pos_bmrw_midmu, range_reco, mom_rw_minchi2);
				if (kMIDeg) Fill1DWHist(h1d_trklen_Pos_bmrw_mideg, range_reco, mom_rw_minchi2);
				if (kMIDother) Fill1DWHist(h1d_trklen_Pos_bmrw_midother, range_reco, mom_rw_minchi2);

			} //beam-mom cut (within 3-sigma)
			else { //tail of beam
				Fill1DWHist(h1d_trklen_Pos_bmrw, range_reco, 1); //set weight to one

				if (kinel) Fill1DWHist(h1d_trklen_Pos_bmrw_inel, range_reco, 1); 
				if (kel) Fill1DWHist(h1d_trklen_Pos_bmrw_el, range_reco, 1); 
				if (kMIDcosmic) Fill1DWHist(h1d_trklen_Pos_bmrw_midcosmic, range_reco, 1); 
				if (kMIDpi) Fill1DWHist(h1d_trklen_Pos_bmrw_midpi, range_reco, 1);
				if (kMIDp) Fill1DWHist(h1d_trklen_Pos_bmrw_midp, range_reco, 1);
				if (kMIDmu) Fill1DWHist(h1d_trklen_Pos_bmrw_midmu, range_reco, 1);
				if (kMIDeg) Fill1DWHist(h1d_trklen_Pos_bmrw_mideg, range_reco, 1);
				if (kMIDother) Fill1DWHist(h1d_trklen_Pos_bmrw_midother, range_reco, 1);

			} //tail of beam			
		} //pos cut




		if (IsPandoraSlice&&IsBQ&&IsCaloSize) { //basic cuts
			//before bmrw
			Fill1DHist(h1d_trklen_BQ, range_reco);
			Fill1DHist(h1d_truetrklen_BQ, range_true);
			//Fill1DHist(h1d_trueendz_BQ, true_endz);
			//Fill1DHist(h1d_trueendzSCE_BQ, true_endz_sce);

			if (kinel) Fill1DHist(h1d_trklen_BQ_inel, range_reco); 
			if (kel) Fill1DHist(h1d_trklen_BQ_el, range_reco); 
			if (kMIDcosmic) Fill1DHist(h1d_trklen_BQ_midcosmic, range_reco); 
			if (kMIDpi) Fill1DHist(h1d_trklen_BQ_midpi, range_reco);
			if (kMIDp) Fill1DHist(h1d_trklen_BQ_midp, range_reco);
			if (kMIDmu) Fill1DHist(h1d_trklen_BQ_midmu, range_reco);
			if (kMIDeg) Fill1DHist(h1d_trklen_BQ_mideg, range_reco);
			if (kMIDother) Fill1DHist(h1d_trklen_BQ_midother, range_reco);

			if (kinel) Fill1DHist(h1d_truetrklen_BQ_inel, range_true); 
			if (kel) Fill1DHist(h1d_truetrklen_BQ_el, range_true); 
			if (kMIDcosmic) Fill1DHist(h1d_truetrklen_BQ_midcosmic, range_true); 
			if (kMIDpi) Fill1DHist(h1d_truetrklen_BQ_midpi, range_true);
			if (kMIDp) Fill1DHist(h1d_truetrklen_BQ_midp, range_true);
			if (kMIDmu) Fill1DHist(h1d_truetrklen_BQ_midmu, range_true);
			if (kMIDeg) Fill1DHist(h1d_truetrklen_BQ_mideg, range_true);
			if (kMIDother) Fill1DHist(h1d_truetrklen_BQ_midother, range_true);


/*
			//before bmrw
			Fill1DHist(h1d_endz_BQ, reco_endz);

			if (kinel) Fill1DHist(h1d_endz_BQ_inel, reco_endz); 
			if (kel) Fill1DHist(h1d_endz_BQ_el, reco_endz); 
			if (kMIDcosmic) Fill1DHist(h1d_endz_BQ_midcosmic, reco_endz); 
			if (kMIDpi) Fill1DHist(h1d_endz_BQ_midpi, reco_endz);
			if (kMIDp) Fill1DHist(h1d_endz_BQ_midp, reco_endz);
			if (kMIDmu) Fill1DHist(h1d_endz_BQ_midmu, reco_endz);
			if (kMIDeg) Fill1DHist(h1d_endz_BQ_mideg, reco_endz);
			if (kMIDother) Fill1DHist(h1d_endz_BQ_midother, reco_endz);

			//before bmrw
			h2d_trklen_endz_BQ->Fill(range_reco, reco_endz);

			if (kinel) h2d_trklen_endz_BQ_inel->Fill(range_reco, reco_endz);
			if (kel) h2d_trklen_endz_BQ_el->Fill(range_reco, reco_endz); 
			if (kMIDcosmic) h2d_trklen_endz_BQ_midcosmic->Fill(range_reco, reco_endz);
			if (kMIDpi) h2d_trklen_endz_BQ_midpi->Fill(range_reco, reco_endz);
			if (kMIDp) h2d_trklen_endz_BQ_midp->Fill(range_reco, reco_endz);
			if (kMIDmu) h2d_trklen_endz_BQ_midmu->Fill(range_reco, reco_endz);
			if (kMIDeg) h2d_trklen_endz_BQ_mideg->Fill(range_reco, reco_endz);
			if (kMIDother) h2d_trklen_endz_BQ_midother->Fill(range_reco, reco_endz);
*/

			double mom_rw_minchi2=1.;
			if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) { //beam-mom cut (within 3-sigma)
				mom_rw_minchi2=bmrw_func->Eval(mom_beam_spec*1000.); //bmrw
				Fill1DWHist(h1d_trklen_BQ_bmrw, range_reco, mom_rw_minchi2);

				if (kinel) Fill1DWHist(h1d_trklen_BQ_bmrw_inel, range_reco, mom_rw_minchi2); 
				if (kel) Fill1DWHist(h1d_trklen_BQ_bmrw_el, range_reco, mom_rw_minchi2); 
				if (kMIDcosmic) Fill1DWHist(h1d_trklen_BQ_bmrw_midcosmic, range_reco, mom_rw_minchi2); 
				if (kMIDpi) Fill1DWHist(h1d_trklen_BQ_bmrw_midpi, range_reco, mom_rw_minchi2);
				if (kMIDp) Fill1DWHist(h1d_trklen_BQ_bmrw_midp, range_reco, mom_rw_minchi2);
				if (kMIDmu) Fill1DWHist(h1d_trklen_BQ_bmrw_midmu, range_reco, mom_rw_minchi2);
				if (kMIDeg) Fill1DWHist(h1d_trklen_BQ_bmrw_mideg, range_reco, mom_rw_minchi2);
				if (kMIDother) Fill1DWHist(h1d_trklen_BQ_bmrw_midother, range_reco, mom_rw_minchi2);

			} //beam-mom cut (within 3-sigma)
			else { //tail of beam
				Fill1DWHist(h1d_trklen_BQ_bmrw, range_reco, 1); //set weight to one

				if (kinel) Fill1DWHist(h1d_trklen_BQ_bmrw_inel, range_reco, 1); 
				if (kel) Fill1DWHist(h1d_trklen_BQ_bmrw_el, range_reco, 1); 
				if (kMIDcosmic) Fill1DWHist(h1d_trklen_BQ_bmrw_midcosmic, range_reco, 1); 
				if (kMIDpi) Fill1DWHist(h1d_trklen_BQ_bmrw_midpi, range_reco, 1);
				if (kMIDp) Fill1DWHist(h1d_trklen_BQ_bmrw_midp, range_reco, 1);
				if (kMIDmu) Fill1DWHist(h1d_trklen_BQ_bmrw_midmu, range_reco, 1);
				if (kMIDeg) Fill1DWHist(h1d_trklen_BQ_bmrw_mideg, range_reco, 1);
				if (kMIDother) Fill1DWHist(h1d_trklen_BQ_bmrw_midother, range_reco, 1);

			} //tail of beam

			
			
		} //basic cuts


		if (IsPandoraSlice&&IsBQ&&IsCaloSize&&IsRecoInEL) { //basic cuts+RecoInel
			//before bmrw
			Fill1DHist(h1d_trklen_RecoInel, range_reco);

			if (kinel) Fill1DHist(h1d_trklen_RecoInel_inel, range_reco); 
			if (kel) Fill1DHist(h1d_trklen_RecoInel_el, range_reco); 
			if (kMIDcosmic) Fill1DHist(h1d_trklen_RecoInel_midcosmic, range_reco); 
			if (kMIDpi) Fill1DHist(h1d_trklen_RecoInel_midpi, range_reco);
			if (kMIDp) Fill1DHist(h1d_trklen_RecoInel_midp, range_reco);
			if (kMIDmu) Fill1DHist(h1d_trklen_RecoInel_midmu, range_reco);
			if (kMIDeg) Fill1DHist(h1d_trklen_RecoInel_mideg, range_reco);
			if (kMIDother) Fill1DHist(h1d_trklen_RecoInel_midother, range_reco);

			double mom_rw_minchi2=1.;
			if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) { //beam-mom cut (within 3-sigma)
				mom_rw_minchi2=bmrw_func->Eval(mom_beam_spec*1000.); //bmrw
				Fill1DWHist(h1d_trklen_RecoInel_bmrw, range_reco, mom_rw_minchi2);

				if (kinel) Fill1DWHist(h1d_trklen_RecoInel_bmrw_inel, range_reco, mom_rw_minchi2); 
				if (kel) Fill1DWHist(h1d_trklen_RecoInel_bmrw_el, range_reco, mom_rw_minchi2); 
				if (kMIDcosmic) Fill1DWHist(h1d_trklen_RecoInel_bmrw_midcosmic, range_reco, mom_rw_minchi2); 
				if (kMIDpi) Fill1DWHist(h1d_trklen_RecoInel_bmrw_midpi, range_reco, mom_rw_minchi2);
				if (kMIDp) Fill1DWHist(h1d_trklen_RecoInel_bmrw_midp, range_reco, mom_rw_minchi2);
				if (kMIDmu) Fill1DWHist(h1d_trklen_RecoInel_bmrw_midmu, range_reco, mom_rw_minchi2);
				if (kMIDeg) Fill1DWHist(h1d_trklen_RecoInel_bmrw_mideg, range_reco, mom_rw_minchi2);
				if (kMIDother) Fill1DWHist(h1d_trklen_RecoInel_bmrw_midother, range_reco, mom_rw_minchi2);

			} //beam-mom cut (within 3-sigma)
			else { //tail of beam
				Fill1DWHist(h1d_trklen_RecoInel_bmrw, range_reco, 1); //set weight to one

				if (kinel) Fill1DWHist(h1d_trklen_RecoInel_bmrw_inel, range_reco, 1); 
				if (kel) Fill1DWHist(h1d_trklen_RecoInel_bmrw_el, range_reco, 1); 
				if (kMIDcosmic) Fill1DWHist(h1d_trklen_RecoInel_bmrw_midcosmic, range_reco, 1); 
				if (kMIDpi) Fill1DWHist(h1d_trklen_RecoInel_bmrw_midpi, range_reco, 1);
				if (kMIDp) Fill1DWHist(h1d_trklen_RecoInel_bmrw_midp, range_reco, 1);
				if (kMIDmu) Fill1DWHist(h1d_trklen_RecoInel_bmrw_midmu, range_reco, 1);
				if (kMIDeg) Fill1DWHist(h1d_trklen_RecoInel_bmrw_mideg, range_reco, 1);
				if (kMIDother) Fill1DWHist(h1d_trklen_RecoInel_bmrw_midother, range_reco, 1);

			} //tail of beam			
		} //basic cuts+RecoInel

		//misid:p-rich
		if (cosine_beam_spec_primtrk<=0.9&&IsPos&&IsRecoInEL&&IsPandoraSlice&&IsCaloSize) { //cos<=0.9+pos+RecoInel
			//before bmrw
			Fill1DHist(h1d_trklen_MidpRich, range_reco);

			if (kinel) Fill1DHist(h1d_trklen_MidpRich_inel, range_reco); 
			if (kel) Fill1DHist(h1d_trklen_MidpRich_el, range_reco); 
			if (kMIDcosmic) Fill1DHist(h1d_trklen_MidpRich_midcosmic, range_reco); 
			if (kMIDpi) Fill1DHist(h1d_trklen_MidpRich_midpi, range_reco);
			if (kMIDp) Fill1DHist(h1d_trklen_MidpRich_midp, range_reco);
			if (kMIDmu) Fill1DHist(h1d_trklen_MidpRich_midmu, range_reco);
			if (kMIDeg) Fill1DHist(h1d_trklen_MidpRich_mideg, range_reco);
			if (kMIDother) Fill1DHist(h1d_trklen_MidpRich_midother, range_reco);

			double mom_rw_minchi2=1.;
			if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) { //beam-mom cut (within 3-sigma)
				mom_rw_minchi2=bmrw_func->Eval(mom_beam_spec*1000.); //bmrw
				Fill1DWHist(h1d_trklen_MidpRich_bmrw, range_reco, mom_rw_minchi2);

				if (kinel) Fill1DWHist(h1d_trklen_MidpRich_bmrw_inel, range_reco, mom_rw_minchi2); 
				if (kel) Fill1DWHist(h1d_trklen_MidpRich_bmrw_el, range_reco, mom_rw_minchi2); 
				if (kMIDcosmic) Fill1DWHist(h1d_trklen_MidpRich_bmrw_midcosmic, range_reco, mom_rw_minchi2); 
				if (kMIDpi) Fill1DWHist(h1d_trklen_MidpRich_bmrw_midpi, range_reco, mom_rw_minchi2);
				if (kMIDp) Fill1DWHist(h1d_trklen_MidpRich_bmrw_midp, range_reco, mom_rw_minchi2);
				if (kMIDmu) Fill1DWHist(h1d_trklen_MidpRich_bmrw_midmu, range_reco, mom_rw_minchi2);
				if (kMIDeg) Fill1DWHist(h1d_trklen_MidpRich_bmrw_mideg, range_reco, mom_rw_minchi2);
				if (kMIDother) Fill1DWHist(h1d_trklen_MidpRich_bmrw_midother, range_reco, mom_rw_minchi2);

			} //beam-mom cut (within 3-sigma)
			else { //tail of beam
				Fill1DWHist(h1d_trklen_MidpRich_bmrw, range_reco, 1); //set weight to one

				if (kinel) Fill1DWHist(h1d_trklen_MidpRich_bmrw_inel, range_reco, 1); 
				if (kel) Fill1DWHist(h1d_trklen_MidpRich_bmrw_el, range_reco, 1); 
				if (kMIDcosmic) Fill1DWHist(h1d_trklen_MidpRich_bmrw_midcosmic, range_reco, 1); 
				if (kMIDpi) Fill1DWHist(h1d_trklen_MidpRich_bmrw_midpi, range_reco, 1);
				if (kMIDp) Fill1DWHist(h1d_trklen_MidpRich_bmrw_midp, range_reco, 1);
				if (kMIDmu) Fill1DWHist(h1d_trklen_MidpRich_bmrw_midmu, range_reco, 1);
				if (kMIDeg) Fill1DWHist(h1d_trklen_MidpRich_bmrw_mideg, range_reco, 1);
				if (kMIDother) Fill1DWHist(h1d_trklen_MidpRich_bmrw_midother, range_reco, 1);

			} //tail of beam			
		} //cos<=0.9+pos+RecoInel


		//El-rich
		//bq
		if (IsRecoStop&&IsBQ&&IsPandoraSlice&&IsCaloSize) { //BQ+El-rich
			//before bmrw
			Fill1DHist(h1d_trklen_ElRich, range_reco);

			if (kinel) Fill1DHist(h1d_trklen_ElRich_inel, range_reco); 
			if (kel) Fill1DHist(h1d_trklen_ElRich_el, range_reco); 
			if (kMIDcosmic) Fill1DHist(h1d_trklen_ElRich_midcosmic, range_reco); 
			if (kMIDpi) Fill1DHist(h1d_trklen_ElRich_midpi, range_reco);
			if (kMIDp) Fill1DHist(h1d_trklen_ElRich_midp, range_reco);
			if (kMIDmu) Fill1DHist(h1d_trklen_ElRich_midmu, range_reco);
			if (kMIDeg) Fill1DHist(h1d_trklen_ElRich_mideg, range_reco);
			if (kMIDother) Fill1DHist(h1d_trklen_ElRich_midother, range_reco);

			double mom_rw_minchi2=1.;
			if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) { //beam-mom cut (within 3-sigma)
				mom_rw_minchi2=bmrw_func->Eval(mom_beam_spec*1000.); //bmrw
				Fill1DWHist(h1d_trklen_ElRich_bmrw, range_reco, mom_rw_minchi2);

				if (kinel) Fill1DWHist(h1d_trklen_ElRich_bmrw_inel, range_reco, mom_rw_minchi2); 
				if (kel) Fill1DWHist(h1d_trklen_ElRich_bmrw_el, range_reco, mom_rw_minchi2); 
				if (kMIDcosmic) Fill1DWHist(h1d_trklen_ElRich_bmrw_midcosmic, range_reco, mom_rw_minchi2); 
				if (kMIDpi) Fill1DWHist(h1d_trklen_ElRich_bmrw_midpi, range_reco, mom_rw_minchi2);
				if (kMIDp) Fill1DWHist(h1d_trklen_ElRich_bmrw_midp, range_reco, mom_rw_minchi2);
				if (kMIDmu) Fill1DWHist(h1d_trklen_ElRich_bmrw_midmu, range_reco, mom_rw_minchi2);
				if (kMIDeg) Fill1DWHist(h1d_trklen_ElRich_bmrw_mideg, range_reco, mom_rw_minchi2);
				if (kMIDother) Fill1DWHist(h1d_trklen_ElRich_bmrw_midother, range_reco, mom_rw_minchi2);

			} //beam-mom cut (within 3-sigma)
			else { //tail of beam
				Fill1DWHist(h1d_trklen_ElRich_bmrw, range_reco, 1); //set weight to one

				if (kinel) Fill1DWHist(h1d_trklen_ElRich_bmrw_inel, range_reco, 1); 
				if (kel) Fill1DWHist(h1d_trklen_ElRich_bmrw_el, range_reco, 1); 
				if (kMIDcosmic) Fill1DWHist(h1d_trklen_ElRich_bmrw_midcosmic, range_reco, 1); 
				if (kMIDpi) Fill1DWHist(h1d_trklen_ElRich_bmrw_midpi, range_reco, 1);
				if (kMIDp) Fill1DWHist(h1d_trklen_ElRich_bmrw_midp, range_reco, 1);
				if (kMIDmu) Fill1DWHist(h1d_trklen_ElRich_bmrw_midmu, range_reco, 1);
				if (kMIDeg) Fill1DWHist(h1d_trklen_ElRich_bmrw_mideg, range_reco, 1);
				if (kMIDother) Fill1DWHist(h1d_trklen_ElRich_bmrw_midother, range_reco, 1);

			} //tail of beam			
		} //BQ+El-rich

		//calosz
		if (IsRecoStop&&IsPandoraSlice&&IsCaloSize) { //calosz+El-rich
			//before bmrw
			Fill1DHist(h1d_trklen_CaloSz_ElRich, range_reco);

			if (kinel) Fill1DHist(h1d_trklen_CaloSz_ElRich_inel, range_reco); 
			if (kel) Fill1DHist(h1d_trklen_CaloSz_ElRich_el, range_reco); 
			if (kMIDcosmic) Fill1DHist(h1d_trklen_CaloSz_ElRich_midcosmic, range_reco); 
			if (kMIDpi) Fill1DHist(h1d_trklen_CaloSz_ElRich_midpi, range_reco);
			if (kMIDp) Fill1DHist(h1d_trklen_CaloSz_ElRich_midp, range_reco);
			if (kMIDmu) Fill1DHist(h1d_trklen_CaloSz_ElRich_midmu, range_reco);
			if (kMIDeg) Fill1DHist(h1d_trklen_CaloSz_ElRich_mideg, range_reco);
			if (kMIDother) Fill1DHist(h1d_trklen_CaloSz_ElRich_midother, range_reco);

			double mom_rw_minchi2=1.;
			if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) { //beam-mom cut (within 3-sigma)
				mom_rw_minchi2=bmrw_func->Eval(mom_beam_spec*1000.); //bmrw
				Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw, range_reco, mom_rw_minchi2);

				if (kinel) Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw_inel, range_reco, mom_rw_minchi2); 
				if (kel) Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw_el, range_reco, mom_rw_minchi2); 
				if (kMIDcosmic) Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw_midcosmic, range_reco, mom_rw_minchi2); 
				if (kMIDpi) Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw_midpi, range_reco, mom_rw_minchi2);
				if (kMIDp) Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw_midp, range_reco, mom_rw_minchi2);
				if (kMIDmu) Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw_midmu, range_reco, mom_rw_minchi2);
				if (kMIDeg) Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw_mideg, range_reco, mom_rw_minchi2);
				if (kMIDother) Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw_midother, range_reco, mom_rw_minchi2);

			} //beam-mom cut (within 3-sigma)
			else { //tail of beam
				Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw, range_reco, 1); //set weight to one

				if (kinel) Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw_inel, range_reco, 1); 
				if (kel) Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw_el, range_reco, 1); 
				if (kMIDcosmic) Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw_midcosmic, range_reco, 1); 
				if (kMIDpi) Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw_midpi, range_reco, 1);
				if (kMIDp) Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw_midp, range_reco, 1);
				if (kMIDmu) Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw_midmu, range_reco, 1);
				if (kMIDeg) Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw_mideg, range_reco, 1);
				if (kMIDother) Fill1DWHist(h1d_trklen_CaloSz_ElRich_bmrw_midother, range_reco, 1);

			} //tail of beam			
		} //calosz+El-rich

		//pos
		if (IsRecoStop&&IsPos&&IsPandoraSlice&&IsCaloSize) { //pos+El-rich
			//before bmrw
			Fill1DHist(h1d_trklen_Pos_ElRich, range_reco);

			if (kinel) Fill1DHist(h1d_trklen_Pos_ElRich_inel, range_reco); 
			if (kel) Fill1DHist(h1d_trklen_Pos_ElRich_el, range_reco); 
			if (kMIDcosmic) Fill1DHist(h1d_trklen_Pos_ElRich_midcosmic, range_reco); 
			if (kMIDpi) Fill1DHist(h1d_trklen_Pos_ElRich_midpi, range_reco);
			if (kMIDp) Fill1DHist(h1d_trklen_Pos_ElRich_midp, range_reco);
			if (kMIDmu) Fill1DHist(h1d_trklen_Pos_ElRich_midmu, range_reco);
			if (kMIDeg) Fill1DHist(h1d_trklen_Pos_ElRich_mideg, range_reco);
			if (kMIDother) Fill1DHist(h1d_trklen_Pos_ElRich_midother, range_reco);

			double mom_rw_minchi2=1.;
			if ((mom_beam_spec*1000.)>=mu_min&&(mom_beam_spec*1000.)<=mu_max) { //beam-mom cut (within 3-sigma)
				mom_rw_minchi2=bmrw_func->Eval(mom_beam_spec*1000.); //bmrw
				Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw, range_reco, mom_rw_minchi2);

				if (kinel) Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw_inel, range_reco, mom_rw_minchi2); 
				if (kel) Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw_el, range_reco, mom_rw_minchi2); 
				if (kMIDcosmic) Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw_midcosmic, range_reco, mom_rw_minchi2); 
				if (kMIDpi) Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw_midpi, range_reco, mom_rw_minchi2);
				if (kMIDp) Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw_midp, range_reco, mom_rw_minchi2);
				if (kMIDmu) Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw_midmu, range_reco, mom_rw_minchi2);
				if (kMIDeg) Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw_mideg, range_reco, mom_rw_minchi2);
				if (kMIDother) Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw_midother, range_reco, mom_rw_minchi2);

			} //beam-mom cut (within 3-sigma)
			else { //tail of beam
				Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw, range_reco, 1); //set weight to one

				if (kinel) Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw_inel, range_reco, 1); 
				if (kel) Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw_el, range_reco, 1); 
				if (kMIDcosmic) Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw_midcosmic, range_reco, 1); 
				if (kMIDpi) Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw_midpi, range_reco, 1);
				if (kMIDp) Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw_midp, range_reco, 1);
				if (kMIDmu) Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw_midmu, range_reco, 1);
				if (kMIDeg) Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw_mideg, range_reco, 1);
				if (kMIDother) Fill1DWHist(h1d_trklen_Pos_ElRich_bmrw_midother, range_reco, 1);

			} //tail of beam			
		} //pos+El-rich


	} //main entry loop

	//save results ---------------------------------------------------------//
   	//TFile *fout = new TFile("mc_proton_after_bmrw.root","RECREATE");
   	//TFile *fout = new TFile("mc_proton_tightxy_after_bmrw.root","RECREATE");
   	//TFile *fout = new TFile("mc_proton_after_bmrw_HD.root","RECREATE");
   	TFile *fout = new TFile("mc_proton_after_bmrw.root","RECREATE");

		h1d_trklen_CaloSz->Write();
		h1d_trklen_CaloSz_inel->Write();
		h1d_trklen_CaloSz_el->Write();
		h1d_trklen_CaloSz_midcosmic->Write();
		h1d_trklen_CaloSz_midpi->Write();
		h1d_trklen_CaloSz_midp->Write();
		h1d_trklen_CaloSz_midmu->Write();
		h1d_trklen_CaloSz_mideg->Write();
		h1d_trklen_CaloSz_midother->Write();

		h1d_trklen_CaloSz_bmrw->Write();
		h1d_trklen_CaloSz_bmrw_inel->Write();
		h1d_trklen_CaloSz_bmrw_el->Write();
		h1d_trklen_CaloSz_bmrw_midcosmic->Write();
		h1d_trklen_CaloSz_bmrw_midpi->Write();
		h1d_trklen_CaloSz_bmrw_midp->Write();
		h1d_trklen_CaloSz_bmrw_midmu->Write();
		h1d_trklen_CaloSz_bmrw_mideg->Write();
		h1d_trklen_CaloSz_bmrw_midother->Write();



		h1d_trklen_Pos->Write();
		h1d_trklen_Pos_inel->Write();
		h1d_trklen_Pos_el->Write();
		h1d_trklen_Pos_midcosmic->Write();
		h1d_trklen_Pos_midpi->Write();
		h1d_trklen_Pos_midp->Write();
		h1d_trklen_Pos_midmu->Write();
		h1d_trklen_Pos_mideg->Write();
		h1d_trklen_Pos_midother->Write();

		h1d_trklen_Pos_bmrw->Write();
		h1d_trklen_Pos_bmrw_inel->Write();
		h1d_trklen_Pos_bmrw_el->Write();
		h1d_trklen_Pos_bmrw_midcosmic->Write();
		h1d_trklen_Pos_bmrw_midpi->Write();
		h1d_trklen_Pos_bmrw_midp->Write();
		h1d_trklen_Pos_bmrw_midmu->Write();
		h1d_trklen_Pos_bmrw_mideg->Write();
		h1d_trklen_Pos_bmrw_midother->Write();



		h1d_trklen_BQ->Write();
		h1d_trklen_BQ_inel->Write();
		h1d_trklen_BQ_el->Write();
		h1d_trklen_BQ_midcosmic->Write();
		h1d_trklen_BQ_midpi->Write();
		h1d_trklen_BQ_midp->Write();
		h1d_trklen_BQ_midmu->Write();
		h1d_trklen_BQ_mideg->Write();
		h1d_trklen_BQ_midother->Write();

		h1d_trklen_BQ_bmrw->Write();
		h1d_trklen_BQ_bmrw_inel->Write();
		h1d_trklen_BQ_bmrw_el->Write();
		h1d_trklen_BQ_bmrw_midcosmic->Write();
		h1d_trklen_BQ_bmrw_midpi->Write();
		h1d_trklen_BQ_bmrw_midp->Write();
		h1d_trklen_BQ_bmrw_midmu->Write();
		h1d_trklen_BQ_bmrw_mideg->Write();
		h1d_trklen_BQ_bmrw_midother->Write();


		h1d_trklen_RecoInel->Write();
		h1d_trklen_RecoInel_inel->Write();
		h1d_trklen_RecoInel_el->Write();
		h1d_trklen_RecoInel_midcosmic->Write();
		h1d_trklen_RecoInel_midpi->Write();
		h1d_trklen_RecoInel_midp->Write();
		h1d_trklen_RecoInel_midmu->Write();
		h1d_trklen_RecoInel_mideg->Write();
		h1d_trklen_RecoInel_midother->Write();

		h1d_trklen_RecoInel_bmrw->Write();
		h1d_trklen_RecoInel_bmrw_inel->Write();
		h1d_trklen_RecoInel_bmrw_el->Write();
		h1d_trklen_RecoInel_bmrw_midcosmic->Write();
		h1d_trklen_RecoInel_bmrw_midpi->Write();
		h1d_trklen_RecoInel_bmrw_midp->Write();
		h1d_trklen_RecoInel_bmrw_midmu->Write();
		h1d_trklen_RecoInel_bmrw_mideg->Write();
		h1d_trklen_RecoInel_bmrw_midother->Write();

		h1d_trklen_MidpRich->Write();
		h1d_trklen_MidpRich_inel->Write();
		h1d_trklen_MidpRich_el->Write();
		h1d_trklen_MidpRich_midcosmic->Write();
		h1d_trklen_MidpRich_midpi->Write();
		h1d_trklen_MidpRich_midp->Write();
		h1d_trklen_MidpRich_midmu->Write();
		h1d_trklen_MidpRich_mideg->Write();
		h1d_trklen_MidpRich_midother->Write();

		h1d_trklen_MidpRich_bmrw->Write();
		h1d_trklen_MidpRich_bmrw_inel->Write();
		h1d_trklen_MidpRich_bmrw_el->Write();
		h1d_trklen_MidpRich_bmrw_midcosmic->Write();
		h1d_trklen_MidpRich_bmrw_midpi->Write();
		h1d_trklen_MidpRich_bmrw_midp->Write();
		h1d_trklen_MidpRich_bmrw_midmu->Write();
		h1d_trklen_MidpRich_bmrw_mideg->Write();
		h1d_trklen_MidpRich_bmrw_midother->Write();



		h1d_trklen_ElRich->Write();
		h1d_trklen_ElRich_inel->Write();
		h1d_trklen_ElRich_el->Write();
		h1d_trklen_ElRich_midcosmic->Write();
		h1d_trklen_ElRich_midpi->Write();
		h1d_trklen_ElRich_midp->Write();
		h1d_trklen_ElRich_midmu->Write();
		h1d_trklen_ElRich_mideg->Write();
		h1d_trklen_ElRich_midother->Write();

		h1d_trklen_ElRich_bmrw->Write();
		h1d_trklen_ElRich_bmrw_inel->Write();
		h1d_trklen_ElRich_bmrw_el->Write();
		h1d_trklen_ElRich_bmrw_midcosmic->Write();
		h1d_trklen_ElRich_bmrw_midpi->Write();
		h1d_trklen_ElRich_bmrw_midp->Write();
		h1d_trklen_ElRich_bmrw_midmu->Write();
		h1d_trklen_ElRich_bmrw_mideg->Write();
		h1d_trklen_ElRich_bmrw_midother->Write();



		h1d_trklen_CaloSz_ElRich->Write();
		h1d_trklen_CaloSz_ElRich_inel->Write();
		h1d_trklen_CaloSz_ElRich_el->Write();
		h1d_trklen_CaloSz_ElRich_midcosmic->Write();
		h1d_trklen_CaloSz_ElRich_midpi->Write();
		h1d_trklen_CaloSz_ElRich_midp->Write();
		h1d_trklen_CaloSz_ElRich_midmu->Write();
		h1d_trklen_CaloSz_ElRich_mideg->Write();
		h1d_trklen_CaloSz_ElRich_midother->Write();

		h1d_trklen_CaloSz_ElRich_bmrw->Write();
		h1d_trklen_CaloSz_ElRich_bmrw_inel->Write();
		h1d_trklen_CaloSz_ElRich_bmrw_el->Write();
		h1d_trklen_CaloSz_ElRich_bmrw_midcosmic->Write();
		h1d_trklen_CaloSz_ElRich_bmrw_midpi->Write();
		h1d_trklen_CaloSz_ElRich_bmrw_midp->Write();
		h1d_trklen_CaloSz_ElRich_bmrw_midmu->Write();
		h1d_trklen_CaloSz_ElRich_bmrw_mideg->Write();
		h1d_trklen_CaloSz_ElRich_bmrw_midother->Write();




		h1d_trklen_Pos_ElRich->Write();
		h1d_trklen_Pos_ElRich_inel->Write();
		h1d_trklen_Pos_ElRich_el->Write();
		h1d_trklen_Pos_ElRich_midcosmic->Write();
		h1d_trklen_Pos_ElRich_midpi->Write();
		h1d_trklen_Pos_ElRich_midp->Write();
		h1d_trklen_Pos_ElRich_midmu->Write();
		h1d_trklen_Pos_ElRich_mideg->Write();
		h1d_trklen_Pos_ElRich_midother->Write();

		h1d_trklen_Pos_ElRich_bmrw->Write();
		h1d_trklen_Pos_ElRich_bmrw_inel->Write();
		h1d_trklen_Pos_ElRich_bmrw_el->Write();
		h1d_trklen_Pos_ElRich_bmrw_midcosmic->Write();
		h1d_trklen_Pos_ElRich_bmrw_midpi->Write();
		h1d_trklen_Pos_ElRich_bmrw_midp->Write();
		h1d_trklen_Pos_ElRich_bmrw_midmu->Write();
		h1d_trklen_Pos_ElRich_bmrw_mideg->Write();
		h1d_trklen_Pos_ElRich_bmrw_midother->Write();

/*
		h1d_endz_CaloSz->Write();
		h1d_endz_CaloSz_inel->Write();
		h1d_endz_CaloSz_el->Write();
		h1d_endz_CaloSz_midcosmic->Write();
		h1d_endz_CaloSz_midpi->Write();
		h1d_endz_CaloSz_midp->Write();
		h1d_endz_CaloSz_midmu->Write();
		h1d_endz_CaloSz_mideg->Write();
		h1d_endz_CaloSz_midother->Write();

		h1d_endz_Pos->Write();
		h1d_endz_Pos_inel->Write();
		h1d_endz_Pos_el->Write();
		h1d_endz_Pos_midcosmic->Write();
		h1d_endz_Pos_midpi->Write();
		h1d_endz_Pos_midp->Write();
		h1d_endz_Pos_midmu->Write();
		h1d_endz_Pos_mideg->Write();
		h1d_endz_Pos_midother->Write();

		h1d_endz_BQ->Write();
		h1d_endz_BQ_inel->Write();
		h1d_endz_BQ_el->Write();
		h1d_endz_BQ_midcosmic->Write();
		h1d_endz_BQ_midpi->Write();
		h1d_endz_BQ_midp->Write();
		h1d_endz_BQ_midmu->Write();
		h1d_endz_BQ_mideg->Write();
		h1d_endz_BQ_midother->Write();

		h2d_trklen_endz_CaloSz->Write();
		h2d_trklen_endz_CaloSz_inel->Write();
		h2d_trklen_endz_CaloSz_el->Write();
		h2d_trklen_endz_CaloSz_midcosmic->Write();
		h2d_trklen_endz_CaloSz_midpi->Write();
		h2d_trklen_endz_CaloSz_midp->Write();
		h2d_trklen_endz_CaloSz_midmu->Write();
		h2d_trklen_endz_CaloSz_mideg->Write();
		h2d_trklen_endz_CaloSz_midother->Write();

		h2d_trklen_endz_Pos->Write();
		h2d_trklen_endz_Pos_inel->Write();
		h2d_trklen_endz_Pos_el->Write();
		h2d_trklen_endz_Pos_midcosmic->Write();
		h2d_trklen_endz_Pos_midpi->Write();
		h2d_trklen_endz_Pos_midp->Write();
		h2d_trklen_endz_Pos_midmu->Write();
		h2d_trklen_endz_Pos_mideg->Write();
		h2d_trklen_endz_Pos_midother->Write();

		h2d_trklen_endz_BQ->Write();
		h2d_trklen_endz_BQ_inel->Write();
		h2d_trklen_endz_BQ_el->Write();
		h2d_trklen_endz_BQ_midcosmic->Write();
		h2d_trklen_endz_BQ_midpi->Write();
		h2d_trklen_endz_BQ_midp->Write();
		h2d_trklen_endz_BQ_midmu->Write();
		h2d_trklen_endz_BQ_mideg->Write();
		h2d_trklen_endz_BQ_midother->Write();
*/

		h1d_truetrklen_CaloSz->Write();
		h1d_truetrklen_CaloSz_inel->Write();
		h1d_truetrklen_CaloSz_el->Write();
		h1d_truetrklen_CaloSz_midcosmic->Write();
		h1d_truetrklen_CaloSz_midpi->Write();
		h1d_truetrklen_CaloSz_midp->Write();
		h1d_truetrklen_CaloSz_midmu->Write();
		h1d_truetrklen_CaloSz_mideg->Write();
		h1d_truetrklen_CaloSz_midother->Write();

		h1d_truetrklen_Pos->Write();
		h1d_truetrklen_Pos_inel->Write();
		h1d_truetrklen_Pos_el->Write();
		h1d_truetrklen_Pos_midcosmic->Write();
		h1d_truetrklen_Pos_midpi->Write();
		h1d_truetrklen_Pos_midp->Write();
		h1d_truetrklen_Pos_midmu->Write();
		h1d_truetrklen_Pos_mideg->Write();
		h1d_truetrklen_Pos_midother->Write();

		h1d_truetrklen_BQ->Write();
		h1d_truetrklen_BQ_inel->Write();
		h1d_truetrklen_BQ_el->Write();
		h1d_truetrklen_BQ_midcosmic->Write();
		h1d_truetrklen_BQ_midpi->Write();
		h1d_truetrklen_BQ_midp->Write();
		h1d_truetrklen_BQ_midmu->Write();
		h1d_truetrklen_BQ_mideg->Write();
		h1d_truetrklen_BQ_midother->Write();


/*
		h1d_trueendz_CaloSz->Write();
		h1d_trueendz_Pos->Write();
		h1d_trueendz_BQ->Write();

		h1d_trueendzSCE_CaloSz->Write();
		h1d_trueendzSCE_Pos->Write();
		h1d_trueendzSCE_BQ->Write();
*/


	fout->Close();



}
