#include "./headers/BasicParameters.h"
#include "./headers/BasicFunctions.h"
//#include "./headers/SliceParams.h"
//#include "./headers/ThinSlice.h"
#include "./headers/util.h"
#include "./headers/BetheBloch.h"

void HandyCalc(){

double mom_beam=2; //unit:GeV/c
double ke_beam_MeV=1000.*p2ke(mom_beam); //ke_beam [MeV]
//double ke_beam_spec_MeV=1000.*ke_beam_spec; //ke_beam_spec [MeV]
double csda=csda_range_vs_mom_sm->Eval(mom_beam);

std::cout<<"mom_beam:"<<mom_beam<<" [GeV/C]"<<std::endl;
std::cout<<"ke_beam_MeV:"<<ke_beam_MeV<<" [MeV]"<<std::endl;
std::cout<<"csda:"<<csda<<" [cm]"<<std::endl;

}
