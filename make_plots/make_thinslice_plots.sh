#!/bin/bash

#file_in="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4_MC_1GeV_reco1_sce_datadriven/thinslice_width4cm_25slices_Etruth.root"
#file_in="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4_MC_1GeV_reco1_sce_datadriven/thinslice_width4cm_25slices.root"
#file_in="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4_MC_1GeV_reco1_sce_datadriven/thinslice_thickness4cm_25slices.root"
file_in="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/prod4a_thinslice_dx4cm_25slcs.root"
fout_path="./plots"

exe1_str="root -b -q 'makeXS.C(\""$file_in\"", \""$fout_path\"")'"
exe2_str="root -b -q 'plotincE.C(\""$file_in\"", \""$fout_path\"")'"
exe3_str="root -b -q 'ploteffpur.C(\""$file_in\"", \""$fout_path\"")'"
exe4_str="root -b -q 'plotAngCorr.C(\""$file_in\"", \""$fout_path\"")'"
exe5_str="root -b -q 'plot_z_ke_true.C(\""$file_in\"", \""$fout_path\"")'"
exe6_str="root -b -q 'plot_KEs.C(\""$file_in\"", \""$fout_path\"")'"


echo $exe1_str" ......"
eval $exe1_str

echo $exe2_str" ......"
eval $exe2_str

echo $exe3_str" ......"
eval $exe3_str

echo $exe4_str" ......"
eval $exe4_str

echo $exe5_str" ......"
eval $exe5_str

echo $exe6_str" ......"
eval $exe6_str

