# DAQ

<a href="https://github.com/marktext/marktext/releases/latest">
   <img src="https://img.shields.io/badge/DAQ-v1.0.0-green">
   <img src="https://img.shields.io/badge/platform-Linux%20%7C%20Mac%20-green">
   <img src="https://img.shields.io/badge/Language-python3-green">
   <img src="https://img.shields.io/badge/Language-C-green">
   <img src="https://img.shields.io/badge/dependencies-tested-green">
   <img src="https://img.shields.io/badge/licence-GNU-green">
</a>      <br>
DAQ is a computational tool using deep learning that can estimate the residue-wise local quality fro protein models from cryo-Electron Microscopy (EM) maps.  

Copyright (C) 2021 Genki Terashi* , Xiao Wang*, Sai Raghavendra Maddhuri Venkata Subramaniya, John J. G. Tesmer, and Daisuke Kihara, and Purdue University. 

License: GPL v3 for academic use. (For commercial use, please contact us for different licensing.)

Contact: Daisuke Kihara (dkihara@purdue.edu)

## Citation:
Genki Terashi* , Xiao Wang*, Sai Raghavendra Maddhuri Venkata Subramaniya, John J. G. Tesmer & Daisuke Kihara. Residue-Wise Local Quality Estimation for Protein Models from Cryo-EM Maps. (2021). [PDF](Submitted)
```
@article{genki2021DAQ,   
  title={Residue-Wise Local Quality Estimation for Protein Models from Cryo-EM Maps},   
  author={Genki Terashi, Xiao Wang, Sai Raghavendra Maddhuri Venkata Subramaniya, John J. G. Tesmer, and Daisuke Kihara},    
  journal={(Submitted)},    
  year={2021}    
}   
```

## Colab Website (Online platform): https://bit.ly/daq-score
**All the functions in this github are available here. Related instructions are included in the Colab website.**

## Introduction
An increasing number of protein structures are determined by cryogenic electron microscopy (cryo-EM). Although the resolution of determined cryo-EM density maps is improving in general, there are still many cases where amino acids of a protein are assigned with different levels of confidence, including those assigned with relatively high ambiguity. Here, we developed a method that identifies potential misassignment of residues in the map, including residue shifts along an otherwise correct main-chain trace. The score, named DAQ, computes the likelihood that the local density corresponds to different amino acids, atoms, and secondary structures from the map density distribution and assesses how well amino acids in the reconstructed model structure agree with the likelihood. DAQ is complementary to existing model validation scores for cryo-EM that examine local density gradient in the map or stereochemical geometry of the structure model. When DAQ was applied to different versions of model structure entries in PDB that were derived from the same density maps, a clear improvement of DAQ-score was observed in the newer versions of the models. The DAQ-score also found potential misassignment errors in a substantial number of over 4400 deposited protein structure models built into cryo-EM maps.


## Overall Protocol
![protocol](https://user-images.githubusercontent.com/50850224/142276557-c79df306-5cf9-40f9-a0b8-f7ef08176a7a.jpeg)

## Overview of DAQ Score
![image](https://user-images.githubusercontent.com/50850224/167702710-5350bcd7-2acb-4054-b424-e751bdf6aee0.png)
where aa(i) is the amino acid type of residue i, P_aa(i)(i) is the computed probability for amino acid type aa(i) for the nearest grid point to the Cα atom of residue i. As shown in the equation, the probability is normalized by the average probability of amino acid type aa(i) across over all atom positions in the protein model.<br>
*If the assignment is correct, DAQ will be positive, and negative if the assignment may be incorrect.<br>
*If a position in the map does not have distinct density pattern for the assigned amino acid (or secondary structure, Calpha atom), DAQ will be close to 0.



## Pre-required software
Python 3 : https://www.python.org/downloads/    
Pymol(for visualization): https://pymol.org/2/   

## Installation  
### 1. [`Install git`](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) 
### 2. Clone the repository in your computer 
```
git clone https://github.com/kiharalab/DAQ && cd DAQ
```

### 3. Build dependencies.   
You have two options to install dependency on your computer:
#### 3.1 Install with pip and python(Ver 3.6.9).
##### 3.1.1[`install pip`](https://pip.pypa.io/en/stable/installing/).
##### 3.1.2  Install dependency in command line.
```
pip3 install -r requirements.txt --user
```
If you encounter any errors, you can install each library one by one:
```
!pip install mrcfile==1.2.0 
!pip install numpy>=1.19.4
!pip install numba>=0.52.0
!pip install torch>=1.6.0
!pip install scipy>=1.6.0
```

#### 3.2 Install with anaconda
##### 3.2.1 [`install conda`](https://bit.ly/daq-score). 
##### 3.2.2 Install dependency in command line
```
conda create -n daq python=3.8.5
conda activate daq
pip install -r requirements.txt 
```
Each time when you want to run my code, simply activate the environment by
```
conda activate daq
conda deactivate(If you want to exit) 
```

## Usage
```
python3 main.py -h:
  -h, --help            show this help message and exit
  -F F                  Map file path
  -M M                  QA deep learning model path, default:"best_model/qa_model/Multimodel.pth"
  -P P                  PDB file path
  --mode MODE           Running Mode
  --stride STRIDE       Stride size for scanning maps (default:1)
  --voxel_size          input voxel size (default:11)
  --gpu GPU             specify the gpu to use
  --batch_size          batch size for inference (default:256)
  --cardinality         ResNeXt cardinality
  --window WINDOW       half window size to smooth the score for output (default:9)
```

## 1. Run DAQ 
Since DAQ(AA) yields the best score 
```
python main.py --mode=0 -F [Map_path]  -P [Structure_path] --window [half_window_size] --stride [stride_size] 
```
**Please Run the script under "DAQ" directory**, otherwise it may raise errors because of the complilation failure. 

Here [Map_path] is the cryo-EM map file path in your computerl [Structure_path] is the protein structure in pdb format; [half_window_size] is half of the window size that used for smoothing the residue-wise score based on a sliding window scanning the entire sequence, here half_window_size=(window_size-1)/2; [stride_size]  is the stride step to scan the maps.<br>
Output will be saved in "Predict_Result_WithPDB/[Input_Map_Name]". 
### Running Example
```
python main.py --mode=0 -F example/2566_3J6B_9.mrc -P example/3J6B_9.pdb --window 9 --stride 2
```
Results of this example is saved in [2566_Result](https://github.com/kiharalab/DAQ/tree/main/result)
### Preparing the input map
If the cryo-EM map grid spacing is not 1, it typically takes longer time to resample the map to have grid spacing 1 by our script. Hence, you can also use [ChimeraX](https://www.rbvi.ucsf.edu/chimerax/) to accelerate the speed by providing the script a resampled map:
```
1 open your map via chimeraX. 
2 In the bottom command line to type command: vol resample #1 spacing 1.0
3 In the bottom command line to type command: save newmap.mrc model #2
4 Then you can use the resampled map to upload
```

## 2. Visualization Result
In Pymol, open "dqa_score_w9.pdb" file, please type the following command line:
```
spectrum b, red_white_blue,  all, -1,1
```
![](https://github.com/kiharalab/DAQ/blob/main/result/visualization.png)
Here blue region means the quality is acceptable while red region means the quality is not so good. Here we put DAQ(AA) score in the b-factor columns. Detailed explanation is [here](https://github.com/kiharalab/DAQ/blob/main/result/README.md).

## Output file
The detailed instructions are [here](https://github.com/kiharalab/DAQ/blob/main/result/README.md).

## Q&A
For possible errors and solutions, please check [QA.md](https://github.com/kiharalab/DAQ/blob/main/QA.md)
