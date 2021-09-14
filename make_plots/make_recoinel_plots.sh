#!/bin/bash

#file_mc_in="../mc_proton_studyRecoInelCut_old.root"
file_mc_in="../mc_proton_studyRecoInelCut.root"
#file_data_in="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_run5387_prod4a_thinslice_dx4cm_25slcs.root"
file_data_in="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/data_proton_bmrw.root"
fout_path="./plots_recoinel"

exe0_str="root -b -q 'plot_recoinel_cut_study.C(\""$file_mc_in\"", \""$file_data_in\"", \""$fout_path\"", \""BQ\"")'"
exe1_str="root -b -q 'plot_recoinel_cut_study.C(\""$file_mc_in\"", \""$file_data_in\"", \""$fout_path\"", \""trklen\"",\""RecoInel\"")'"
exe2_str="root -b -q 'plot_recoinel_cut_study.C(\""$file_mc_in\"", \""$file_data_in\"", \""$fout_path\"", \""trklen\"",\""RecoInel2\"")'"
exe3_str="root -b -q 'plot_recoinel_cut_study.C(\""$file_mc_in\"", \""$file_data_in\"", \""$fout_path\"", \""trklen\"",\""RecoInel2\"")'"
exe4_str="root -b -q 'plot_recoinel_cut_study.C(\""$file_mc_in\"", \""$file_data_in\"", \""$fout_path\"", \""chi2pid\"",\""BQ\"")'"
exe5_str="root -b -q 'plot_recoinel_cut_study.C(\""$file_mc_in\"", \""$file_data_in\"", \""$fout_path\"", \""mediandedx\"",\""BQ\"")'"
exe6_str="root -b -q 'plot_recoinel_cut_study.C(\""$file_mc_in\"", \""$file_data_in\"", \""$fout_path\"", \""dEdL\"",\""BQ\"")'"


#echo $exe0_str" ......"
#eval $exe0_str

echo $exe1_str" ......"
#eval $exe1_str

echo $exe2_str" ......"
#eval $exe2_str

echo $exe3_str" ......"
#eval $exe3_str


#eval $exe4_str

#eval $exe5_str

eval $exe6_str

