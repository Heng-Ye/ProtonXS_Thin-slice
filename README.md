# Proton Cross-section Analysis using Thin-slice Algorithm <br/>

Run procedure:<br/>
- Edit the main analysis code: **ProtonThinSliceX.C**
- Edit the input file path in  **files_prod4a.txt**
- Execute the thin-slice analysis\
**./compile.sh**

*Reference about unfolding:
- Use RooUnfold for unfolding: https://gitlab.cern.ch/RooUnfold/RooUnfold<br/>
- RooUnfold Installation:<br/>
git clone https://gitlab.cern.ch/RooUnfold/RooUnfold.git<br/>
mkdir RooUnfold/build<br/>
cd RooUnfold/build<br/>
cmake ..<br/>
make -j4<br/>
