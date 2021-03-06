import ROOT as RT
import sys
from math import floor, acos, cos, sqrt, tan, ceil
from glob import glob as ls
from array import array

RT.gROOT.SetBatch(1)

mag_P1 = 5.82044830e-3;
mag_P3 = -4.68880000e-6;
mag_P4 = 324.573967;

#f = TFile(sys.argv[1])
#t = f.Get("NTuples/GoodParticle")
files = ls(sys.argv[1]) 
t = RT.TChain("NTuples/GoodParticle")

if( len(sys.argv) > 5 ): max_files = int(sys.argv[5])
else: max_files = 1000


nFiles = 0
for f in files: 

  if(nFiles > max_files): break

  print("Adding", f)
  t.Add(f)
  nFiles = nFiles + 1

def shift_x(x):
  return ( floor(x) + .5 )

def to_f(x):
  if x > 0.: return int(ceil(x))
  else: return int(floor(x))

def momentum_costheta(x1,x2,x3):
  L1 = 1.98
  L2 = 1.69472
  L3 = 2.11666
  fBeamBend = .12003

  a = (x2*L3 - x3*L2)*cos(fBeamBend)/(L3-L2);

 
  numTerm = (a - x1)*( (L3 - L2)*tan(fBeamBend) + (x3 - x2)*cos(fBeamBend) ) + L1*( L3 - L2 );

  denomTerm1 = sqrt( L1*L1 + (a - x1)*(a - x1) );
  denomTerm2 = sqrt( ( (L3 - L2)*tan(fBeamBend) + (x3 - x2)*cos(fBeamBend) )**2 + ( (L3 - L2) )**2 );
  denom = denomTerm1 * denomTerm2;

  cosTheta = numTerm/denom;  
  return cosTheta;


fout = RT.TFile(sys.argv[3], "RECREATE")
h = RT.TH1D("h", "", 150,0.5,1.5)
hTrue = RT.TH1D("hTrue", "", 150,0.5,1.5)
true_vs_reco = RT.TH2D("true_vs_reco", "", 150,0.5,1.5,150,0.5,1.5)

r = RT.TH1D("r", "", 200,-1.,1.)
r_vs_true = RT.TH2D("r_vs_true", "", 150,0.5,1.5, 200,-1.,1.)
r_vs_p = RT.TH2D("r_vs_p", "", 150,0.5,1.5, 200,-1.,1.)

hF1 = RT.TH1D("hF1", "", 192, 0, 192)
hF2 = RT.TH1D("hF2", "", 192, 0, 192)
hF3 = RT.TH1D("hF3", "", 192, 0, 192)

f1 = array("i", [0])
f2 = array("i", [0])
f3 = array("i", [0])

x1 = array("d", [0])
x2 = array("d", [0])
x3 = array("d", [0])

true_p = array("d", [0.])
preSpec_p = array("d", [0.])
TOF1_p = array("d", [0.])
postSpec_p = array("d", [0.])
NP04front_p = array("d", [0.])
NP04FieldCage_p = array("d", [0.])
reco_p = array("d", [0.])
reco_tof = array("d", [0.])
flip_p = array("d", [0.])
plus_p_1 = array("d", [0.])
minus_p_1 = array("d", [0.])
plus_p_2 = array("d", [0.])
minus_p_2 = array("d", [0.])
plus_p_3 = array("d", [0.])
minus_p_3 = array("d", [0.])

plus_shift = array('d', [0.])
minus_shift = array('d', [0.])

PDG = array("i", [0])
fcPDG = array("i", [0])

TOF1_TrackID = array("i", [0])
TRIG2_TrackID = array("i", [0])
TOF1_ParentID = array("i", [0])
TRIG2_ParentID = array("i", [0])
TOF1_EventID = array("i", [0])
TRIG2_EventID = array("i", [0])

outtree = RT.TTree("tree","")
outtree.Branch("f1", f1, "f1/I")
outtree.Branch("f2", f2, "f2/I")
outtree.Branch("f3", f3, "f3/I")
outtree.Branch("x1", x1, "x1/D")
outtree.Branch("x2", x2, "x2/D")
outtree.Branch("x3", x3, "x3/D")

outtree.Branch("true_p", true_p, "true_p/D")
outtree.Branch("preSpec_p", preSpec_p, "preSpec_p/D")
outtree.Branch("TOF1_p", TOF1_p, "TOF1_p/D")
outtree.Branch("postSpec_p", postSpec_p, "postSpec_p/D")
outtree.Branch("NP04front_p", NP04front_p, "NP04front_p/D")
outtree.Branch("NP04FieldCage_p", NP04FieldCage_p, "NP04FieldCage_p/D")

outtree.Branch("TOF1_TrackID", TOF1_TrackID, "TOF1_TrackID/I")
outtree.Branch("TRIG2_TrackID", TRIG2_TrackID, "TRIG2_TrackID/I")
outtree.Branch("TOF1_ParentID", TOF1_ParentID, "TOF1_ParentID/I")
outtree.Branch("TRIG2_ParentID", TRIG2_ParentID, "TRIG2_ParentID/I")
outtree.Branch("TOF1_EventID", TOF1_EventID, "TOF1_EventID/I")
outtree.Branch("TRIG2_EventID", TRIG2_EventID, "TRIG2_EventID/I")

outtree.Branch("reco_p", reco_p, "reco_p/D")
outtree.Branch("reco_tof", reco_tof, "reco_tof/D")
outtree.Branch("flip_p", flip_p, "flip_p/D")
outtree.Branch("plus_p_1", plus_p_1, "plus_p_1/D")
outtree.Branch("minus_p_1", minus_p_1, "minus_p_1/D")
outtree.Branch("plus_p_2", plus_p_2, "plus_p_2/D")
outtree.Branch("minus_p_2", minus_p_2, "minus_p_2/D")
outtree.Branch("plus_p_3", plus_p_3, "plus_p_3/D")
outtree.Branch("minus_p_3", minus_p_3, "minus_p_3/D")
outtree.Branch("plus_shift", plus_shift, "plus_shift/D")
outtree.Branch("minus_shift", minus_shift, "minus_shift/D")

outtree.Branch("PDG", PDG, "PDG/I")
outtree.Branch("fcPDG", fcPDG, "fcPDG/I")

def momentum(x1,x2,x3):
  return 299792458*LB/(1.E9 * acos(momentum_costheta(x1,x2,x3)))


count = 0
print("Tree has", t.GetEntries(), "entries")
for e in t:

  PDG[0] = int(e.NP04front_PDGid)
  fcPDG[0] = int(e.NP04FieldCage_PDGid)

  TOF1_TrackID[0] = int(e.TOF1_TrackID)
  TRIG2_TrackID[0] = int(e.TRIG2_TrackID)
  TOF1_ParentID[0] = int(e.TOF1_ParentID)
  TRIG2_ParentID[0] = int(e.TRIG2_ParentID)
  TOF1_EventID[0] = int(e.TOF1_EventID)
  TRIG2_EventID[0] = int(e.TRIG2_EventID)
  
  if not (count % 1000): print(count)
  count = count + 1
  px = e.AfterTarget_Px
  py = e.AfterTarget_Py
  pz = e.AfterTarget_Pz

  true_p[0] = sqrt( px**2 + py**2 + pz**2 )/1.e3
  #if true_p[0] < 1.e-5: continue

  x1[0] = -1.e-3*shift_x(e.BPROF1_x)
  x2[0] = -1.e-3*shift_x(e.BPROF2_x)
  x3[0] = -1.e-3*shift_x(e.BPROF3_x)

  f1[0] = int( 96 - e.BPROF1_x)
  f2[0] = int( 96 - e.BPROF2_x)
  f3[0] = int( 96 - e.BPROF3_x)

  preSpec_p[0] = sqrt( e.BPROF1_Px**2 + e.BPROF1_Py**2 + e.BPROF1_Pz**2 )
  TOF1_p[0] = sqrt( e.TOF1_Px**2 + e.TOF1_Py**2 + e.TOF1_Pz**2 )
  postSpec_p[0] = sqrt( e.BPROF3_Px**2 + e.BPROF3_Py**2 + e.BPROF3_Pz**2 )
  NP04front_p[0] = sqrt( e.NP04front_Px**2 + e.NP04front_Py**2 + e.NP04front_Pz**2 )
  NP04FieldCage_p[0] = sqrt( e.NP04FieldCage_Px**2 + e.NP04FieldCage_Py**2 + e.NP04FieldCage_Pz**2 )

  if NP04front_p[0] < 1.e-5: continue

  hTrue.Fill(true_p[0])
  hF1.Fill( f1[0] )
  hF2.Fill( f2[0] )
  hF3.Fill( f3[0] )

  reco_tof[0] = (e.TRIG2_t - e.TOF1_t)

  #print momentum_costheta(x1,x2,x3), acos(momentum_costheta(x1,x2,x3))

  #1GeV: 68.8 
  #2GeV: 137.5
  #3GeV: 206.2

  if(   int(sys.argv[2]) == 1 ): current = 68.8
  elif( int(sys.argv[2]) == 2 ): current = 137.5 
  elif( int(sys.argv[2]) == 3 ): current = 206.2 
  elif( int(sys.argv[2]) == 6 ): current = 419.7

  LB = mag_P1 * current 
  deltaI = current  - mag_P4
  if deltaI > 0: LB = LB + mag_P3 * deltaI**2
  #print 299792458*LB/(1.E9 * acos(momentum_costheta(x1,x2,x3)));
  #reco_p[0] = 299792458*LB/(1.E9 * acos(momentum_costheta(x1,x2,x3)))
  reco_p[0] = momentum(x1[0], x2[0], x3[0] )
  flip_p[0] = momentum(-1.*x1[0], -1.*x2[0], -1.*x3[0] )
  
  h.Fill(reco_p[0])
  #r.Fill(1. - reco_p[0] / true_p[0] )
  #r_vs_true.Fill(true_p[0], 1. - reco_p[0] / true_p[0] )
  #r_vs_p.Fill(reco_p[0], 1. - reco_p[0] / true_p[0] )
  true_vs_reco.Fill( reco_p[0], true_p[0] )


  plus_p_1[0] =  momentum(x1[0]-.5e-3, x2[0] ,x3[0])
  minus_p_1[0] = momentum(x1[0]+.5e-3, x2[0] ,x3[0])

  plus_p_2[0] =  momentum(x1[0], x2[0]-.5e-3 ,x3[0])
  minus_p_2[0] = momentum(x1[0], x2[0]+.5e-3 ,x3[0])

  plus_p_3[0] =  momentum(x1[0], x2[0], x3[0]-.5e-3)
  minus_p_3[0] = momentum(x1[0], x2[0], x3[0]+.5e-3)

  plus_shift[0] =  momentum(x1[0], x2[0], x3[0]-.18e-3)
  minus_shift[0] = momentum(x1[0], x2[0], x3[0]+.18e-3)

  outtree.Fill()

fout.cd()
h.Write()
hTrue.Write()
r.Write()
r_vs_true.Write()
r_vs_p.Write()

hF1.Write()
hF2.Write()
hF3.Write()

true_vs_reco.Write()

outtree.Write()

##Making the resolution histograms
binning = sys.argv[4]
outtree.Draw("(reco_p*1.e3 - NP04front_p)/NP04front_p>>hPiRes(" + binning + ")", "PDG == 211 || PDG == -13")
outtree.Draw("(reco_p*1.e3 - NP04front_p)/NP04front_p>>hProtRes(" + binning + ")", "PDG == 2212")
outtree.Draw("(reco_p*1.e3 - NP04front_p)/NP04front_p>>hElecRes(" + binning + ")", "PDG == -11")

hPiRes = RT.gDirectory.Get("hPiRes")
hProtRes = RT.gDirectory.Get("hProtRes")
hElecRes = RT.gDirectory.Get("hElecRes")

hPiRes.Write()
hProtRes.Write()
hElecRes.Write()

fout.Close()
