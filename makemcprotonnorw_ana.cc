#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>

#include "TString.h"
#include "TFile.h"
#include "TApplication.h"
#include "TTree.h"
#include "TChain.h"
#include "./headers/color_text.h"

using namespace std;

int main(int argc, char *argv[], char *envp[]) {
//int main() {
 //assert(argc == 2);
 //TString str_in(argv[1]);
 TString fgoodlist(argv[1]);
 TString fclsname(argv[2]);
 //TString fgoodlist="goodfile_list.txt";

 //TString str_out(argv[2]);
 //TApplication theApp("App", &argc, argv); //necessary for tree

 //setup output file names ---------------------//
 //std::cout<<"Input root file path : "<<str_in<<std::endl;
 //std::cout<<"Output root file : "<<str_out<<std::endl;
 
 ifstream f_string(fgoodlist.Data());
 string buffer_string;
 vector<string> filename;
 while (f_string>>buffer_string) { filename.push_back(buffer_string); }
  f_string.close();

 std::cout<<"Looping over the good beam data list : "<<red<<fgoodlist.Data()<<std::endl;

 //TChain *chain=new TChain("protonanalysis/PandoraBeam");
 TChain *chain=new TChain("protonmcnorw/PandoraBeam");
 for (size_t ii=0; ii<filename.size();ii++) { //loop over all the selected ana files
   //TString inputfilename(filename[ii]);
   std::cout<<green<<"--> reading beam data:  "<<filename[ii]<<std::endl;

   chain->Add(filename[ii].c_str());
   
 } //loop over all the selected ana files

 std::cout<<reset<<"\n"<<std::endl;


 chain->MakeClass(Form("%s",fclsname.Data()));
 //chain->MakeClass("ProtonMCSelector");
 //chain->MakeClass("PionMCSelector");




}                   
