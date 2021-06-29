
void RunAna(TString class_name){
 //load class
 TString tmp_name=Form("ANAMC");
 gROOT->ProcessLine(Form(".L %s.C+", class_name.Data()));
 gROOT->ProcessLine(Form("%s %s",class_name.Data(), tmp_name.Data()));
 gROOT->ProcessLine(Form("%s.Loop()", tmp_name.Data())); 

 //gROOT->ProcessLine("ANA.Show(0)"); 

}
