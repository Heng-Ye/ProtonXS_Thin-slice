#!/bin/bash

#file_in="../prod4a_thinslice_dx4cm_24slcs.root"
file_in="../prod4a_thinslice_dx4cm_25slcs.root"
fout_path="./plots"
#file_in="../prod4a_thinslice_dx4cm_25slcs.root"
#fout_path="./plots_25slices"

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

exe12_str="root -b -q 'plot_trklenreco.C(\""$file_in\"", \""$fout_path\"", \""NoCut\"")'"
exe13_str="root -b -q 'plot_trklenreco.C(\""$file_in\"", \""$fout_path\"", \""PanS\"")'"
exe14_str="root -b -q 'plot_trklenreco.C(\""$file_in\"", \""$fout_path\"", \""CaloSz\"")'"
exe15_str="root -b -q 'plot_trklenreco.C(\""$file_in\"", \""$fout_path\"", \""BQ\"")'"
exe16_str="root -b -q 'plot_trklenreco.C(\""$file_in\"", \""$fout_path\"", \""RecoInel\"")'"

exe17_str="root -b -q 'plot_dtrklen.C(\""$file_in\"", \""$fout_path\"", \""NoCut\"")'"
exe18_str="root -b -q 'plot_dtrklen.C(\""$file_in\"", \""$fout_path\"", \""PanS\"")'"
exe19_str="root -b -q 'plot_dtrklen.C(\""$file_in\"", \""$fout_path\"", \""CaloSz\"")'"
exe20_str="root -b -q 'plot_dtrklen.C(\""$file_in\"", \""$fout_path\"", \""BQ\"")'"
exe21_str="root -b -q 'plot_dtrklen.C(\""$file_in\"", \""$fout_path\"", \""RecoInel\"")'"

exe22_str="root -b -q 'plot_ntrklenreco.C(\""$file_in\"", \""$fout_path\"", \""BQ\"")'"
exe23_str="root -b -q 'plot_reco_cosineTheta.C(\""$file_in\"", \""$fout_path\"")'"

echo $exe0_str" ......"
eval $exe0_str

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

echo $exe7_str" ......"
eval $exe7_str

echo $exe8_str" ......"
eval $exe8_str

echo $exe9_str" ......"
eval $exe9_str

echo $exe10_str" ......"
eval $exe10_str

echo $exe11_str" ......"
eval $exe11_str

echo $exe12_str" ......"
eval $exe12_str

echo $exe13_str" ......"
eval $exe13_str

echo $exe14_str" ......"
eval $exe14_str

echo $exe15_str" ......"
eval $exe15_str

echo $exe16_str" ......"
eval $exe16_str

echo $exe17_str" ......"
eval $exe17_str

echo $exe18_str" ......"
eval $exe18_str

echo $exe19_str" ......"
eval $exe19_str

echo $exe20_str" ......"
eval $exe20_str

echo $exe21_str" ......"
eval $exe21_str

echo $exe22_str" ......"
eval $exe22_str

echo $exe23_str" ......"
eval $exe23_str


