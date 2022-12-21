void model_root_productor() {

  string line;  int n=2; //2 columns
  //GENIE_ha2018 ---------------------------------------------------------------------------------------------//
  vector<double> x_GENIE_ha2018;
  vector<double> y_GENIE_ha2018;
  vector<double> data_GENIE_ha2018;

  //ke reac sig
  double buffer_GENIE_ha2018;
  ifstream f_in_GENIE_ha2018("GENIE_ha2018.dat");
  while (f_in_GENIE_ha2018>>buffer_GENIE_ha2018) { data_GENIE_ha2018.push_back(buffer_GENIE_ha2018); }
  f_in_GENIE_ha2018.close();  for (int k=0; k<data_GENIE_ha2018.size(); k+=3) {
    x_GENIE_ha2018.push_back(1000.*data_GENIE_ha2018[k]); //KE in MeV
    y_GENIE_ha2018.push_back(data_GENIE_ha2018[k+1]); //[mb]
  }

  TGraph *xs_GENIE_ha2018=new TGraph(x_GENIE_ha2018.size(), &x_GENIE_ha2018.at(0), &y_GENIE_ha2018.at(0));
  xs_GENIE_ha2018->SetName("xs_GENIE_ha2018");

  //NEUT_2019 ---------------------------------------------------------------------------------------------//
  vector<double> x_NEUT_2019;
  vector<double> y_NEUT_2019;
  vector<double> data_NEUT_2019;

  double buffer_NEUT_2019;
  ifstream f_in_NEUT_2019("NEUT_2019.dat");
  while (f_in_NEUT_2019>>buffer_NEUT_2019) { data_NEUT_2019.push_back(buffer_NEUT_2019); }
  f_in_NEUT_2019.close();  for (int k=0; k<data_NEUT_2019.size(); k+=n) {
    x_NEUT_2019.push_back(1000.*data_NEUT_2019[k]); //KE in MeV
    y_NEUT_2019.push_back(data_NEUT_2019[k+1]); //[mb]
  }

  TGraph *xs_NEUT_2019=new TGraph(x_NEUT_2019.size(), &x_NEUT_2019.at(0), &y_NEUT_2019.at(0));
  xs_NEUT_2019->SetName("xs_NEUT_2019");

  //NuWRO_2019 ---------------------------------------------------------------------------------------------//
  vector<double> x_NuWRO_2019;
  vector<double> y_NuWRO_2019;
  vector<double> data_NuWRO_2019;

  double buffer_NuWRO_2019;
  ifstream f_in_NuWRO_2019("NuWRO_2019.dat");
  while (f_in_NuWRO_2019>>buffer_NuWRO_2019) { data_NuWRO_2019.push_back(buffer_NuWRO_2019); }
  f_in_NuWRO_2019.close();  for (int k=0; k<data_NuWRO_2019.size(); k+=n) {
    x_NuWRO_2019.push_back(data_NuWRO_2019[k]); //KE in MeV
    y_NuWRO_2019.push_back(data_NuWRO_2019[k+1]); //[mb]
  }

  TGraph *xs_NuWRO_2019=new TGraph(x_NuWRO_2019.size(), &x_NuWRO_2019.at(0), &y_NuWRO_2019.at(0));
  xs_NuWRO_2019->SetName("xs_NuWRO_2019");

  //GENIE_hN2018 ---------------------------------------------------------------------------------------------//
  vector<double> x_GENIE_hN2018;
  vector<double> y_GENIE_hN2018;
  vector<double> data_GENIE_hN2018;

  //ke reac sig
  double buffer_GENIE_hN2018;
  ifstream f_in_GENIE_hN2018("GENIE_hN2018.dat");
  while (f_in_GENIE_hN2018>>buffer_GENIE_hN2018) { data_GENIE_hN2018.push_back(buffer_GENIE_hN2018); }
  f_in_GENIE_hN2018.close();  for (int k=0; k<data_GENIE_hN2018.size(); k+=3) {
    x_GENIE_hN2018.push_back(1000.*data_GENIE_hN2018[k]); //KE in MeV
    y_GENIE_hN2018.push_back(data_GENIE_hN2018[k+1]); //[mb]
  }

  TGraph *xs_GENIE_hN2018=new TGraph(x_GENIE_hN2018.size(), &x_GENIE_hN2018.at(0), &y_GENIE_hN2018.at(0));
  xs_GENIE_hN2018->SetName("xs_GENIE_hN2018");

  //GENIE_INCL_pp ---------------------------------------------------------------------------------------------//
  vector<double> x_GENIE_INCL_pp;
  vector<double> y_GENIE_INCL_pp;
  vector<double> data_GENIE_INCL_pp;

  //ke reac sig
  double buffer_GENIE_INCL_pp;
  ifstream f_in_GENIE_INCL_pp("GENIE_INCL_pp.dat");
  while (f_in_GENIE_INCL_pp>>buffer_GENIE_INCL_pp) { data_GENIE_INCL_pp.push_back(buffer_GENIE_INCL_pp); }
  f_in_GENIE_INCL_pp.close();  for (int k=0; k<data_GENIE_INCL_pp.size(); k+=3) {
    x_GENIE_INCL_pp.push_back(1000.*data_GENIE_INCL_pp[k]); //KE in MeV
    y_GENIE_INCL_pp.push_back(data_GENIE_INCL_pp[k+1]); //[mb]
  }

  TGraph *xs_GENIE_INCL_pp=new TGraph(x_GENIE_INCL_pp.size(), &x_GENIE_INCL_pp.at(0), &y_GENIE_INCL_pp.at(0));
  xs_GENIE_INCL_pp->SetName("xs_GENIE_INCL_pp");




  TFile *fout=new TFile("model_pred.root","create");
    xs_GENIE_ha2018->Write();
    xs_NEUT_2019->Write();
    xs_NuWRO_2019->Write();
    xs_GENIE_hN2018->Write();
    xs_GENIE_INCL_pp->Write();

  fout->Close();




}
