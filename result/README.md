# Output files 
In this example, we used map EMD-2566 and part of its aligned structure(chain 9, PDB-ID: 3J6B). In the released, we only output DAQ(AA) score since it performed best among all 3 scores.

1. 2566_3J6B_9.mrc: Input map file in mrc format.
2. 2566_3J6B_9_new.mrc: Resized map file with grid_size=1.
3. 2566_3J6B_9.trimmap: an input file for the network that saves the normalized box input.
4. prediction.txt: a prediction output file by the network, which includes all our predicted probabilities for all input boxes of the map.
5. daq_raw_score.pdb: an output pdb with raw score, where the score is directly computed by the network without window average. For each ATOM line, the Amino-Acid Type scores (AAscores) are saved in b-factor column. The "#ATOM=" lines indicate the C-alpha ATOM scores and AA scores.
```
ATOM      0  CA  VAL 9  29     375.211 218.122 165.067  1.00 -2.37
#ATOM= 0.440 AA= -2.370
ATOM      1  CA  ILE 9  30     374.127 214.846 166.700  1.00 -1.00
#ATOM= 0.110 AA= -1.000
ATOM      2  CA  TYR 9  31     377.109 214.570 169.056  1.00  1.75
#ATOM= 1.463 AA= 1.745
ATOM      3  CA  LEU 9  32     377.281 211.233 170.917  1.00  0.37
#ATOM= 1.057 AA= 0.366
ATOM      4  CA  HIS 9  33     380.233 211.908 173.279  1.00  0.69
#ATOM= 0.502 AA= 0.688
ATOM      5  CA  LYS 9  34     382.772 209.030 173.465  1.00  1.21
#ATOM= 1.562 AA= 1.207
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
