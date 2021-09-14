#!/bin/bash

#file_in="../prod4a_thinslice_dx4cm_24slcs.root"
file_in="../prod4a_thinslice_dx4cm_25slcs.root"
#file_in="../prod4a_thinslice_dx4cm_25slcs_1stHitKEff.root"
#file_in="../prod4a_thinslice_dx5cm_20slcs_1stHitKEff.root"
file_data_in="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_run5387_prod4a_thinslice_dx4cm_25slcs.root"
fout_path="./plots"
#fout_path="./plots_KEff_1stHit"
#file_in="../prod4a_thinslice_dx4cm_25slcs.root"
#fout_path="./plots_25slices"

#for bkg study
file_in2="../mc_proton_studyRecoInelCut.root"
file_data_in2="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_proton_bmrw.root"

#for bkg fit study
file_in3="../prod4a_bkgstudy_dx4cm_25slcs.root"
file_data_in3="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_run5387_prod4reco2_bkgstudy_dx4cm_25slcs.root"

exe0_str="root -b -q 'plotslices.C(\""$fout_path\"")'"
exe1_str="root -b -q 'makeXS.C(\""$file_in\"", \""$fout_path\"")'"
exe2_str="root -b -q 'plotincE.C(\""$file_in\"", \""$fout_path\"")'"
exe3_str="root -b -q 'ploteffpur.C(\""$file_in\"", \""$fout_path\"")'"
exe4_str="root -b -q 'plotAngCorr.C(\""$file_in\"", \""$fout_path\"")'"
exe5_str="root -b -q 'plot_trklen_ke_true.C(\""$file_in\"", \""$fout_path\"")'"
exe6_str="root -b -q 'plot_KEs.C(\""$file_in\"", \""$fout_path\"")'"
exe7_str="root -b -q 'plot_trklentrue.C(\""$file_in\"", \""$fout_path\"", \""NoCut\"")'"
exe8_str="root -b -q 'plot_trklentrue.C(\""$file_in\"", \""$fout_path\"", \""PanS\"")'"
exe9_str="root -b -q 'plot_trklentrue.C(\""$file_in\"", \""$fout_path\"", \""CaloSz\"")'"
exe10_str="root -b -q 'plot_trklentrue.C(\""$file_in\"", \""$fout_path\"", \""BQ\"")'"
exe11_str="root -b -q 'plot_trklentrue.C(\""$file_in\"", \""$fout_path\"", \""RecoInel\"")'"

exe12_str="root -b -q 'plot_trklenreco.C(\""$file_in\"", \""$file_data_in\"", \""$fout_path\"", \""NoCut\"")'"
exe13_str="root -b -q 'plot_trklenreco.C(\""$file_in\"", \""$file_data_in\"", \""$fout_path\"", \""PanS\"")'"
exe14_str="root -b -q 'plot_trklenreco.C(\""$file_in\"", \""$file_data_in\"", ,\""$fout_path\"", \""CaloSz\"")'"
exe15_str="root -b -q 'plot_trklenreco.C(\""$file_in\"", \""$file_data_in\"", \""$fout_path\"", \""BQ\"")'"
exe16_str="root -b -q 'plot_trklenreco.C(\""$file_in\"", \""$file_data_in\"", \""$fout_path\"", \""RecoInel\"")'"

exe17_str="root -b -q 'plot_dtrklen.C(\""$file_in\"", \""$fout_path\"", \""NoCut\"")'"
exe18_str="root -b -q 'plot_dtrklen.C(\""$file_in\"", \""$fout_path\"", \""PanS\"")'"
exe19_str="root -b -q 'plot_dtrklen.C(\""$file_in\"", \""$fout_path\"", \""CaloSz\"")'"
exe20_str="root -b -q 'plot_dtrklen.C(\""$file_in\"", \""$fout_path\"", \""BQ\"")'"
exe21_str="root -b -q 'plot_dtrklen.C(\""$file_in\"", \""$fout_path\"", \""RecoInel\"")'"

exe22_str="root -b -q 'plot_ntrklenreco.C(\""$file_in\"", \""$file_data_in\"", \""$fout_path\"", \""BQ\"")'"
exe23_str="root -b -q 'plot_reco_cosineTheta.C(\""$file_in\"", \""$file_data_in\"", \""reco_cosineTheta\"", \""$fout_path\"")'"
exe23_str1="root -b -q 'plot_reco_cosineTheta.C(\""$file_in\"", \""$file_data_in\"", \""reco_cosineTheta_Pos\"", \""$fout_path\"")'"

exe24_str="root -b -q 'plot_hdeltaXYZ.C(\""$file_in\"", \""$file_data_in\"", \""$fout_path\"", \""X\"")'"
exe25_str="root -b -q 'plot_hdeltaXYZ.C(\""$file_in\"", \""$file_data_in\"", \""$fout_path\"", \""Y\"")'"
exe26_str="root -b -q 'plot_hdeltaXYZ.C(\""$file_in\"", \""$file_data_in\"", \""$fout_path\"", \""Z\"")'"
exe27_str="root -b -q 'plot_hdeltaXYZ.C(\""$file_in\"", \""$file_data_in\"", \""$fout_path\"", \""XY\"")'"

exe28_str="root -b -q 'plot_startXYZ.C(\""$file_in\"", \""$file_data_in\"", \""$fout_path\"")'"

exe29_str1="root -b -q 'study_bkg_using_cosineTheta.C(\""$file_in\"", \""$file_data_in\"", \""reco_cosineTheta_Pos\"", \""$fout_path\"")'"
exe29_str2="root -b -q 'study_bkg_using_chi2pid.C(\""$file_in2\"", \""$file_data_in2\"", \""chi2pid_BQ\"", \""$fout_path\"")'"

#exe30_str1="root -b -q 'plot_recosliceid.C(\""$file_in\"", \""$fout_path\"", \""recoid_all\"", \""h_recosliceid_allevts_cuts\"", \""All\ Protons\"")'"
#exe30_str2="root -b -q 'plot_recosliceid.C(\""$file_in\"", \""$fout_path\"", \""recoid_inelastic\"", \""h_recosliceid_recoinelastic_cuts\"", \""Proton\ Inelastic\ Scatterings\"")'"

exe30_str1="root -b -q 'plot_recosliceid.C(\""$file_in\"", \""$fout_path\"", \""recoid_all\"", \""h_recosliceid_allevts_cuts\"", \""\"")'"
exe30_str2="root -b -q 'plot_recosliceid.C(\""$file_in\"", \""$fout_path\"", \""recoid_inelastic\"", \""h_recosliceid_recoinelastic_cuts\"", \""\"")'"

#exe31_str="root -b -q 'plot_bkgfit_recosliceid.C(\""$file_in3\"", \""$file_data_in3\"", \""$fout_path\"", \""recosliceid_misidp_rich\"", \""h_recosliceid_cosLE09\"", \""cos\#Theta\<\=0.9\"")'"
#exe31_str="root -b -q 'plot_bkgfit_recosliceid.C(\""$file_in3\"", \""$file_data_in3\"", \""$fout_path\"", \""recosliceid_misidp_rich\"", \""h_recosliceid_cosLE08\"", \""cos\#Theta\<\=0.8\"")'"
exe31_str1="root -b -q 'plot_bkgfit_recosliceid.C(\""$file_in3\"", \""$file_data_in3\"", \""$fout_path\"", \""recosliceid_misidp_rich_LE07\"", \""h_recosliceid_cosLE07\"", \""cos\#Theta\<\=0.7\"")'"
exe31_str2="root -b -q 'plot_bkgfit_recosliceid.C(\""$file_in3\"", \""$file_data_in3\"", \""$fout_path\"", \""recosliceid_misidp_rich_LE08\"", \""h_recosliceid_cosLE08\"", \""cos\#Theta\<\=0.8\"")'"
exe31_str3="root -b -q 'plot_bkgfit_recosliceid.C(\""$file_in3\"", \""$file_data_in3\"", \""$fout_path\"", \""recosliceid_misidp_rich_LE09\"", \""h_recosliceid_cosLE09\"", \""cos\#Theta\<\=0.9\"")'"

#echo $exe29_str1" ......"
#eval $exe29_str1

#echo $exe29_str2" ......"
#eval $exe29_str2

#eval $exe30_str1
#eval $exe30_str2

eval $exe31_str3
