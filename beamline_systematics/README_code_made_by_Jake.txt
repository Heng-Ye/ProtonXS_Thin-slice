Procedure of the Beam momentum systematics:
Python code author: Jake Calcutt

[1]Create tree file for analysis
py3_reco_momentum_tuples.py '/pnfs/dune/persistent/dunepro/beam_data/simulation/mcc10/H4_v34b_1GeV_-27.7_10M_*root' 1 reco_momentum.root '100, -1., 1.'

[2]Estimation of systematics: 
-Calculation of beam momentum shift due to the fibers was shift in the third fiber plane 
Run as:
[2.1]tree->Draw("(minus_shift - reco_p)/reco_p>>h(1000)","PDG == 2212")
[2.2]tree->Draw("(plus_shift - reco_p)/reco_p>>h(1000)","PDG == 2212")
-> Get the mean values of these 2 histograms

[2.3]reco/true difference:
tree->Draw("(reco_p*1.e3 - NP04front_p)/NP04front_p>>h(1000)","PDG == 2212")


