# Output files 
In this example, we used map EMD-2566 and part of its aligned structure(chain 9, PDB-ID: 3J6B). In the released, we only output DAQ(AA) score since it performed best among all 3 scores.

1. 2566_3J6B_9.mrc: Input map file in mrc format.
2. 2566_3J6B_9_new.mrc: Resized map file with grid_size=1.
3. 2566_3J6B_9.trimmap: an input file for the network that saves the normalized box input.
4. prediction.txt: a prediction output file by the network, which includes all our predicted probabilities for all input boxes of the map.
   For each line, colum are corresponding to 0-2:xyz, 3-22:AA type Probabilities, 23-28: Atom Probabilities, 29-31: Secondary structure Probabilities.
6. daq_raw_score.pdb: an output pdb with raw score, where the score is directly computed by the network without window average. For each ATOM line, the Amino-Acid Type scores (AAscores) are saved in b-factor column. The "#ATOM=" lines indicate the C-alpha ATOM scores and AA scores. The "#ProbAA" lines indicate the probabilities of 20 Amino-acid types. The "#DAQAA" lines indicate the DAQ(AA) scores of 20 Amino-acid types.
```
ATOM      0  CA  VAL 9  29     375.211 218.122 165.067  1.00 -2.37
#ATOM= 0.440 AA= -2.370
#ProbAA COORDS:375.211,218.122,165.067,ALA:0.004,VAL:0.004,PHE:0.033,PRO:0.022,MET:0.016,ILE:0.007,LEU:0.019,ASP:0.028,GLU:0.052,LYS:0.240,ARG:0.120,SER:0.012,THR:0.012,TYR:0.052,HIS:0.132,CYS:0.001,ASN:0.090,TRP:0.009,GLN:0.143,GLY:0.004
#DAQAA COORDS:375.211,218.122,165.067,ALA:-2.363,VAL:-2.370,PHE:-0.486,PRO:-0.462,MET:-0.632,ILE:-2.024,LEU:-1.647,ASP:-0.672,GLU:-0.258,LYS:1.152,ARG:0.460,SER:-1.430,THR:-1.368,TYR:-0.119,HIS:1.410,CYS:-2.066,ASN:0.395,TRP:-0.755,GLN:0.877,GLY:-1.933
ATOM      1  CA  ILE 9  30     374.127 214.846 166.700  1.00 -1.00
#ATOM= 0.110 AA= -1.000
#ProbAA COORDS:374.127,214.846,166.700,ALA:0.007,VAL:0.013,PHE:0.032,PRO:0.064,MET:0.010,ILE:0.020,LEU:0.040,ASP:0.043,GLU:0.061,LYS:0.149,ARG:0.041,SER:0.023,THR:0.033,TYR:0.037,HIS:0.085,CYS:0.001,ASN:0.194,TRP:0.020,GLN:0.122,GLY:0.005
#DAQAA COORDS:374.127,214.846,166.700,ALA:-1.949,VAL:-1.261,PHE:-0.520,PRO:0.622,MET:-1.054,ILE:-1.000,LEU:-0.890,ASP:-0.236,GLU:-0.096,LYS:0.675,ARG:-0.615,SER:-0.786,THR:-0.341,TYR:-0.456,HIS:0.963,CYS:-1.745,ASN:1.163,TRP:-0.007,GLN:0.722,GLY:-1.707
ATOM      2  CA  TYR 9  31     377.109 214.570 169.056  1.00  1.75
#ATOM= 1.463 AA= 1.745
#ProbAA COORDS:377.109,214.570,169.056,ALA:0.005,VAL:0.002,PHE:0.235,PRO:0.008,MET:0.012,ILE:0.005,LEU:0.014,ASP:0.006,GLU:0.003,LYS:0.016,ARG:0.128,SER:0.011,THR:0.005,TYR:0.338,HIS:0.093,CYS:0.004,ASN:0.014,TRP:0.075,GLN:0.019,GLY:0.007
#DAQAA COORDS:377.109,214.570,169.056,ALA:-2.227,VAL:-2.933,PHE:1.483,PRO:-1.399,MET:-0.888,ILE:-2.502,LEU:-1.947,ASP:-2.272,GLU:-3.056,LYS:-1.550,ARG:0.528,SER:-1.512,THR:-2.298,TYR:1.745,HIS:1.057,CYS:-0.691,ASN:-1.477,TRP:1.309,GLN:-1.120,GLY:-1.385
```
6. daq_score_w9.pdb: a final overall output pdb with window averaged score. Here the scores are saved in b-factor column.
7. daq_score_w9_9.pdb： a chain-based pdb for window averaged daq score, here it's for the chain 9.
8. visualization.png: visualization result of the daq score.  

# Example Analysis
This example is discussed in the paper “Residue-Wise Local Quality Estimation for Protein Models from Cryo-EM Maps” by Genki Terashi,  Xiao Wang, Sai Raghavendra Maddhuri Venkata Subramaniya, John J. G. Tesmer, and Daisuke Kihara (in submission, 2022) in Fig. 4c.

This protein is Chain 9 of the Ribonuclease III domain from PDB entries 3J6B. When 3J6B is compared with another PDB entry 5MRF, the sequence identity of these pairs is 99.5%. However, residue 227 to 241 are shifted between the two structures. 

In the 3.2 Å resolution EM map (EMD-2566, 3J6B), density for the terminal helix region (His227-Lys241) is not of high enough quality on its own to confirm choices on residue or atom identity. Positions in this model truncated as alanine underscore the lack of interpretability in this region. There are, however, a few telltale signs that there are misalignments, such as hydrophobic side-chains exposed to solvent (3J6B Ile232, which aligns with 5MRF Asn234) and polar residues packed into the core of the fold (3J6B His227 vs. 5MRF Leu229, 3J6B Ser228 vs. 5MRF Val230, and 3J6B Ser219 vs. 5MRF Leu220). At these positions, side chains in 5MRF fit the environment well. For these reasons, 5MRF-9 seems to be a better fit for both EM maps, consistent with the DAQ(AA) score, and the inconsistent regions in 3J6B (227 to 241), which correspond to the red regions in the associated image, may have errors.

To have a more clear understanding of the region, we attached the comparison in visualization.png. The left panel shows the DAQ(AA) score for this structure (EMD-2566, 3J6B chain 9). The misalignment region is shown in red by DAQ, which indicates lower scores for this region. In the right panel, we compared the assignment of 3J6B-9 (shown in green) and 5MRF-9 (shown in cyan) in the helix region where DAQ outputs low scores. We can see the difference of the residue assignment in 3J6B-9 from 5MRF-9.
![](https://github.com/kiharalab/DAQ/blob/main/result/visualization.png)
