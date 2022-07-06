#!/bin/bash
#data sample
file_bmrw_data="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/realdata/p1gev/code_timedep_trkpos/protodune-sp_runset_5387_reco2/rw/data_proton_beamxy_beammom_bmrw.root"

#mc sample [before bmrw, range]
#file_mc="/dune/app/users/hyliao/WORK/analysis/protodune/proton/analysis/mcdata/sce/MC_PDSPProd4a_MC_1GeV_reco1_sce_datadriven_v1/xs_thinslice/rw/mc_proton_beamxy_beammom_bmrw.root"

#mc sample [before bmrw, calo]
#file_mc="../rw/mc_proton_beamxy_beammom_bmrw_calo.root"
file_mc="../rw/mc_proton_beamxy_beammom_bmrw_old_usebeamspec_as_default.root"

#mc sample [after bmrw] ---------------------------------------------------//
#use truth momentum as baseline for bmrw
#file_mc_bmrw="../mc_proton_beamxy_beammom_afterbmrw.root"
#fout_path="./plots_beamxy_bmrw"

#use spec as baseline for bmrw
#file_mc_bmrw="../mc_proton_beamxy_beammom_afterbmrw_usespecasinput.root"
#fout_path="./plots_beamxy_bmrw_old"

#use truth as baseline for bmrw
#file_mc_bmrw="../mc_proton_beamxy_beammom_calo_afterbmrw.root"
#fout_path="./plots_beamxy_bmrw_calo"

file_mc_bmrw="../rw/mc_proton_beamxy_beammom_afterbmrw_usespecasinput.root"
fout_path="./plots_beamxy_bmrw_old"


#exe1_str="root -b -q 'plotBMRW.C(\""$file_bmrw_data\"", \""$file_bmrw_mc\"", \""$fout_path\"")'"
#exe1_str="root -b -q 'plotDataMCKE.C(\""$file_bmrw_data\"", \""$file_bmrw_mc\"", \""$fout_path\"")'"


#draw mean and sigma for data and MC(before bmrw) 
exe1_str="root -b -q 'plotDataMCKE.C(\""$file_bmrw_data\"", \""$file_mc\"", \""$fout_path\"")'"

#calculate bmrw
#note: Edit plotBMRW.C before running the code
exe2_str="root -b -q 'plotBMRW.C(\""$file_bmrw_data\"", \""$file_mc\"", \""$fout_path\"")'"


#calculate bmrw(after bmrw)
exe3_str="root -l 'plotDataMCKE_BMRW.C(\""$file_bmrw_data\"", \""$file_mc\"", \""$file_mc_bmrw\"", \""$fout_path\"")'"



#eval $exe1_str
#eval $exe2_str
eval $exe3_str

