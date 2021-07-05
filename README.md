# Proton Cross-section Analysis using Thin Slice Alg. <br/>

Run procedure:<br/>
*Edit the main analysis code: #ProtonThinSliceX.C<br/>
*Edit the input file path in  #files_prod4a.txt<br/>
*Execute the thin-slice analysis<br/>
./compile.sh

*Reference:
Use RooUnfold for unfolding: https://gitlab.cern.ch/RooUnfold/RooUnfold<br/>
-RooUnfold Installation:<br/>
git clone https://gitlab.cern.ch/RooUnfold/RooUnfold.git<br/>
mkdir RooUnfold/build<br/>
cd RooUnfold/build<br/>
cmake ..<br/>
make -j4<br/>
