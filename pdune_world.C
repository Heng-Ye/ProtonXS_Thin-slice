
void pdune_world() {
  gSystem->Load("libGeom");
  gSystem->Load("libGdml");
  
  //TGeoManager* geom = TGeoManager::Import("/dune/app/users/hyliao/WORK/larsoft_mydev/v09_32_00_e20_prof_test/srcs/dunetpc/dune/Geometry/gdml/protodune_v8_refactored_nowires.gdml");
  //TGeoManager* geom = TGeoManager::Import("/dune/app/users/hyliao/WORK/larsoft_mydev/v09_32_00_e20_prof_test/srcs/dunetpc/dune/Geometry/gdml/protodune_v7_refactored_nowires.gdml");
  TGeoManager* geom = TGeoManager::Import("/dune/app/users/hyliao/WORK/larsoft_mydev/v09_32_00_e20_prof_version_for_reco2_and_prod4a/srcs/dunetpc/dune/Geometry/gdml/protodune_v7_refactored_nowires.gdml");
  //TGeoManager* geom = TGeoManager::Import("/dune/app/users/hyliao/WORK/larsoft_mydev/v09_32_00_e20_prof_test/srcs/dunetpc/dune/Geometry/gdml/protodune_v8_refactored.gdml");
  //TBrowser*       b = new TBrowser();
  //geom->GetTopVolume()->Draw();
  //geom->SetVisLevel(4);
  //geom->GetTopVolume()->Draw("ogl");
  //
  //

   // Creating a view
   //TView3D *view = (TView3D*) TView::CreateView(1);
   //view->SetRange(5,5,5,25,25,25);

  //geom->GetTopVolume()->Draw("");


 TGeoIterator next(gGeoManager->GetTopVolume());
 TGeoNode *node = 0;

 gGeoManager->GetVolume("volSteelShell")->SetLineColor(19);
 gGeoManager->GetVolume("volSteelShell")->SetVisibility(1);
 //gGeoManager->GetVolume("volSteelShell")->SetTransparency(90);
 gGeoManager->GetVolume("volSteelShell")->SetTransparency(100);

 gGeoManager->GetVolume("volGaseousArgon")->SetLineColor(kYellow-7);
 gGeoManager->GetVolume("volGaseousArgon")->SetVisibility(1);
 gGeoManager->GetVolume("volGaseousArgon")->SetTransparency(85);


 while ( (node=(TGeoNode*)next()) ){
   const char* nm = node->GetName();

   if( (strncmp(nm, "volCathode", 10) == 0) ){
     //node->GetVolume()->SetLineColor(kOrange+3); 
     //node->GetVolume()->SetVisibility(1);
     //node->GetVolume()->SetTransparency(30);
     node->GetVolume()->SetLineColor(0); 
     node->GetVolume()->SetTransparency(100);
     node->GetVolume()->SetVisibility(0);

   }
   if( (strncmp(nm, "volTPCActive", 12) == 0) ){
     node->GetVolume()->SetLineColor(kGreen-7); 
     //node->GetVolume()->SetVisibility(1);
     //node->GetVolume()->SetTransparency(80);
     node->GetVolume()->SetTransparency(100);
     node->GetVolume()->SetVisibility(0);

   }
   if( (strncmp(nm, "volAPAFrame", 11) == 0) ){
     node->GetVolume()->SetLineColor(kGray); 
     //node->GetVolume()->SetVisibility(1);
     //node->GetVolume()->SetTransparency(20);

     node->GetVolume()->SetTransparency(100);
     node->GetVolume()->SetVisibility(0);
   }
   if( (strncmp(nm, "volG10Board", 11) == 0) ){
     node->GetVolume()->SetLineColor(kMagenta-10); 
     //node->GetVolume()->SetVisibility(1);
     //node->GetVolume()->SetTransparency(40);
     node->GetVolume()->SetVisibility(0);
     node->GetVolume()->SetTransparency(100);
   }
   if( (strncmp(nm, "volTPCPlane", 11) == 0) ){
     node->GetVolume()->SetLineColor(kBlue-9); 
     //node->GetVolume()->SetVisibility(1);
     //node->GetVolume()->SetTransparency(50);
     node->GetVolume()->SetTransparency(100);
     node->GetVolume()->SetVisibility(0);
   }
   if( (strncmp(nm, "volTPCInner", 11) == 0) ){
     node->GetVolume()->SetLineColor(kWhite); 
     node->GetVolume()->SetVisibility(1);
     node->GetVolume()->SetTransparency(100);
     node->GetVolume()->SetVisibility(0);
   }
   if( (strncmp(nm, "volTPCOuter", 11) == 0) ){
     node->GetVolume()->SetLineColor(kWhite); 
     node->GetVolume()->SetVisibility(1);
     node->GetVolume()->SetTransparency(100);
     node->GetVolume()->SetVisibility(0);
   }
   if( (strncmp(nm, "volOpDetSensitive", 17) == 0) ){
     node->GetVolume()->SetLineColor(kRed-4); 
     //node->GetVolume()->SetVisibility(1);
     //node->GetVolume()->SetTransparency(10);
     node->GetVolume()->SetTransparency(100);
     node->GetVolume()->SetVisibility(0);
   }
   if( (strncmp(nm, "volWorld", 8) == 0) ){
     //node->GetVolume()->SetLineColor(kOrange-7); 
     //node->GetVolume()->SetVisibility(1);
     node->GetVolume()->SetLineColor(0); 
     node->GetVolume()->SetTransparency(100);
     node->GetVolume()->SetVisibility(0);
   }
   if( (strncmp(nm, "volFoamPadding", 14) == 0) ){
     //node->GetVolume()->SetLineColor(kCyan-10); 
     //node->GetVolume()->SetVisibility(1);
     //node->GetVolume()->SetTransparency(95);
     node->GetVolume()->SetLineColor(0); 
     node->GetVolume()->SetTransparency(100);
     node->GetVolume()->SetVisibility(0);
	
   }
   if( (strncmp(nm, "volSteelSupport", 15) == 0) ){
     //node->GetVolume()->SetLineColor(kGray); 
     //node->GetVolume()->SetVisibility(1);
     //node->GetVolume()->SetTransparency(95);
     node->GetVolume()->SetLineColor(0);      
     node->GetVolume()->SetTransparency(100);
     node->GetVolume()->SetVisibility(0);
   }

 }

  gGeoManager->GetTopNode();
  //gGeoManager->CheckOverlaps(1e-5);
  //gGeoManager->PrintOverlaps();
  gGeoManager->SetMaxVisNodes(70000);

  //gGeoManager->GetTopVolume()->Draw("ogl");
  //gGeoManager->FindVolumeFast("volWorld")->Draw("ogl");
  gGeoManager->FindVolumeFast("volWorld")->Draw("ogl");
  //gGeoManager->FindVolumeFast("volTPCPlaneUInner")->Draw("ogl");
  //if ( ! volName.IsNull() ) gGeoManager->FindVolumeFast(volName)->Draw("ogl");
//  gGeoManager->FindVolumeFast("volCryostat")->Draw("X3D");



  //TFile *tf = new TFile("protodune_v7.root", "RECREATE"); 
    //gGeoManager->Write();
  //tf->Close();

}
