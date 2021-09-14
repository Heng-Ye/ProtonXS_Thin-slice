#!/bin/bash

#[1]Generate the file list for analysis
#class_name="ProtonCounters"
#class_name="ProtonThinSlice"
#class_name="ProtonMomentumReweight"
#class_name="ProtonAfterMomentumReweight"
class_name="ProtonBackgroundFit"
#class_name="ProtonRecoInelCutStudy"
#class_name="ProtonEvtClassification"
#class_name="ProtonEvtCounter"
class_namex=$class_name"X"
echo $class_namex


#[2]Generate the ana module [to get the data structure of selected trees]
g++ makemcprotonnorw_ana.cc `root-config --libs --cflags` -o makeproton_ana

#[3]Run the ana module (input can be changable if needed but still need compile to loop over the selected files)
#./makeproton_ana files_prod4a.txt $class_name
./makeproton_ana files_prod4a_new.txt $class_name

#[4]Fix bugs in the generated makeclass module & create analysis file based on the template
sed '/Init(tree)\;/i if (tree-\>InheritsFrom(\"TChain\")) ((TChain\*)tree)-\>LoadTree(0);' $class_name".h" > $class_name"_0.h"
mv $class_name"_0.h" $class_name".h"
cp -prv $class_namex".C" $class_name".C"

#[5]Run analysis code
root_exe_str="root -b -q 'RunAna.C(\""$class_name\"")'"
echo $root_exe_str" ......"
eval $root_exe_str
