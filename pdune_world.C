
void pdune_world() {
  //TGeoManager* geom = TGeoManager::Import("/dune/app/users/hyliao/WORK/larsoft_mydev/v09_32_00_e20_prof_test/srcs/dunetpc/dune/Geometry/gdml/protodune_v8_refactored_nowires.gdml");
  TGeoManager* geom = TGeoManager::Import("/dune/app/users/hyliao/WORK/larsoft_mydev/v09_32_00_e20_prof_test/srcs/dunetpc/dune/Geometry/gdml/protodune_v8_refactored.gdml");
  TBrowser*       b = new TBrowser();
  geom->GetTopVolume()->Draw();
  //geom->SetVisLevel(4);
  //geom->GetTopVolume()->Draw("ogle");


}
