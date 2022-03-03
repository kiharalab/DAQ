# Output files 
In this example, we used map EMD-2566 and part of its aligned structure chain 9 of structure with PDB-ID 3J6B.

1. 2566_3J6B_9.mrc: Input map file in mrc format.
2. 2566_3J6B_9_new.mrc: Resized map file with grid_size=1
3. 2566_3J6B_9.trimmap: input file for the network that saves the normalized box input
4. prediction.txt: prediction output by the network, which includes all our predicted probabilities for all input boxes of the map.
5. dqa_raw_score.pdb: an output pdb with raw score, where the score is directly computed by the network without window average. Here the scores are saved in b-factor column.
6. dqa_score_w9.pdb: a final overall output pdb with window averaged score. Here the scores are saved in b-factor column.
7. dqa_score_w9_9.pdb： a chain-based pdb for window averaged daq score, here it's for the chain 9.
8. visualization.png: visualization result of the daq score.  
