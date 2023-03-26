//#ifndef UNFOLD_H
//#define UNFOLD_H

R__LOAD_LIBRARY(libRooUnfold.so) //load share lib
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"

class Unfold {
  
 public:

  Unfold(TH3D* h3d_reco, TH3D* h3d_true);
  RooUnfoldResponse response_SliceID_3D;   //3D response matrix for all: INC_ST(x-axis), INC_Int(y-axis), INT(z-axis)

  void SaveHistograms();

};


Unfold::Unfold(TH3D* h3d_reco, TH3D* h3d_true)
  : response_SliceID_3D(h3d_reco, h3d_true)
{
  response_SliceID_3D.UseOverflow(false);
}  

void Unfold::SaveHistograms(){
  response_SliceID_3D.Write("response_SliceID_3D");
}

//#endif
