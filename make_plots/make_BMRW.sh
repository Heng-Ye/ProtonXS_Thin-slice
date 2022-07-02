#!/bin/bash

file_bmrw_data="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/rw/data_proton_beamxy_beammom_bmrw.root"
file_bmrw_mc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/rw/mc_proton_beamxy_beammom_bmrw.root"
fout_path="./plots_beamxy_bmrw"

#exe1_str="root -b -q 'plotBMRW.C(\""$file_bmrw_data\"", \""$file_bmrw_mc\"", \""$fout_path\"")'"
#exe1_str="root -b -q 'plotDataMCKE.C(\""$file_bmrw_data\"", \""$file_bmrw_mc\"", \""$fout_path\"")'"


#draw mean and sigma for data and MC 
exe1_str="root -l 'plotDataMCKE.C(\""$file_bmrw_data\"", \""$file_bmrw_mc\"", \""$fout_path\"")'"

#calculate bmrw
exe2_str="root -b -q 'plotBMRW.C(\""$file_bmrw_data\"", \""$file_bmrw_mc\"", \""$fout_path\"")'"




#eval $exe1_str
eval $exe2_str
