//Basic parameters ------------------------------------------//
//PDG
//const int pdg=211; //pi+
const int pdg=2212; //proton
const double m_proton=0.938272046; //proton_mass, unit:GeV/c2

//TPC boundary
double minX =  -360.0;
double maxX = 360.0;
double minY =0.0;
double maxY = 600.0;
double minZ =  0.0; // G10 at z=1.2
double maxZ = 695.0;

//some constants
double NA=6.02214076e23;
double MAr=39.95; //gmol
double Density = 1.39; // g/cm^3

//Misc. Cut values
//XY-Cut
//double mean_x=-29.73; //prod4 mc
//double mean_y=422.4;
//double dev_x=1.5*4.046;
//double dev_y=1.5*3.679;

//double mean_x_data=-26.58; //prod4 data
//double mean_y_data=423.5; //prod4 data
//double dev_x_data=1.5*3.753; //prod4 data
//double dev_y_data=1.5*4.354; //prod4 data

double mean_x=-29.25; //prod4a mc
double mean_y=422.1; //prod4a mc
double dev_x=1.5*4.4; //prod4a mc
double dev_y=1.5*3.883; //prod4a mc

double mean_x_data=-26.59; //prod4a data
double mean_y_data=423.5; //prod4a data
double dev_x_data=1.5*3.744; //prod4a data
double dev_y_data=1.5*4.364; //prod4a data


//dedx cut
double dedx_min=30.;

//beam quality cut
double min1_dx=0.; //new p3
double min2_dx=2.0; //new p3
double min1_dy=-2.0; //new p3
double min2_dy=2.0; //new p3
double min1_dz=-.5; //new p3
double min2_dz=1.5; //new p3

double min1_z=5.10816e-02-3.*2.13366e-01; //p4 
double min2_z=5.10816e-02+3.*2.13366e-01; //p4
double min1_y=4.21863e+02-3.*4.11359e+00; //p4
double min2_y=4.21863e+02+3.*4.11359e+00; //p4
double min1_x=-3.05895e+01-3.*4.69242e+00; //p4  
double min2_x=-3.05895e+01+3.*4.69242e+00; //p4

double cosine_beam_primtrk_min=(9.92194e-01)-4.*(3.96921e-03); //new p4 [spec]


//stopping proton cut
double mean_norm_trklen_csda=9.32064e-01; //prod4 spec
double sigma_norm_trklen_csda=1.32800e-02; //prod4 spec
double min_norm_trklen_csda=mean_norm_trklen_csda-2.*sigma_norm_trklen_csda;
double max_norm_trklen_csda=mean_norm_trklen_csda+3.*sigma_norm_trklen_csda;



