=============================================
Procedure of the Beam momentum systematics:
Python code author: Jake Calcutt
---------------------------------------------

[1] Create tree file for analysis
py3_reco_momentum_tuples.py '/pnfs/dune/persistent/dunepro/beam_data/simulation/mcc10/H4_v34b_1GeV_-27.7_10M_*root' 1 reco_momentum.root '100, -1., 1.'

[2] Estimation of systematics: 
 -Calculate the beam momentum shift/smearing due to the fiber position shift [1.45+-0.18 mm]
 Run as:
 [2.1]tree->Draw("(minus_shift - reco_p)/reco_p>>h(1000)","PDG == 2212")
 [2.2]tree->Draw("(plus_shift - reco_p)/reco_p>>h(1000)","PDG == 2212")
 -> Get the mean values of these 2 histograms
 For protons, similiar result like pions, ~0.7%

[2.3] reco/true difference:
tree->Draw("(reco_p*1.e3 - NP04front_p)/NP04front_p>>h(1000)","PDG == 2212")

[3]Overall systematics:
(a) 1% on B-field -> 1% on p (Info from beamline experts)
    since p=(299.7924/theta)*int_{0}^{L_{mag}} (B*dl)

(b) 0.7% calculated from [2]
Overall systematic uncewrtainty=sqrt(1%^2+0.7%^2)
