# DAQ

<a href="https://github.com/marktext/marktext/releases/latest">
   <img src="https://img.shields.io/badge/DAQ-v1.0.0-green">
   <img src="https://img.shields.io/badge/platform-Linux%20%7C%20Mac%20-green">
   <img src="https://img.shields.io/badge/Language-python3-green">
   <img src="https://img.shields.io/badge/Language-C-green">
   <img src="https://img.shields.io/badge/dependencies-tested-green">
   <img src="https://img.shields.io/badge/licence-GNU-green">
</a>      
DAQ is a computational tool using deep learning that can estimate the residue-wise local quality in cryo-Electron Microscopy (EM) maps.  

Copyright (C) 2021 Genki Terashi* , Xiao Wang*, Sai Raghavendra Maddhuri Venkata Subramaniya, John J. G. Tesmer, and Daisuke Kihara, and Purdue University. 

License: GPL v3 for academic use. (For commercial use, please contact us for different licensing.)

Contact: Daisuke Kihara (dkihara@purdue.edu)

## Citation:
Genki Terashi* , Xiao Wang*, Sai Raghavendra Maddhuri Venkata Subramaniya, John J. G. Tesmer & Daisuke Kihara. Residue-Wise Local Quality Estimation for Protein Models from Cryo-EM Maps.  bioArxiv (2021). [PDF]()
```
@article{genki2021DAQ,   
  title={Residue-Wise Local Quality Estimation for Protein Models from Cryo-EM Maps},   
  author={Genki Terashi, Xiao Wang, Sai Raghavendra Maddhuri Venkata Subramaniya, John J. G. Tesmer, and Daisuke Kihara},    
  journal={BioArxiv},    
  year={2021}    
}   
```

## Colab Website (Online platform): https://bit.ly/daq-score

## Introduction
An increasing number of protein structures are determined by cryogenic electron microscopy (cryo-EM). Although the resolution of determined cryo-EM density maps is improving in general, there are still many cases where amino acids of a protein are assigned with different levels of confidence, including those assigned with relatively high ambiguity, to the map. Here, we developed a method that identifies potential miss-assignment residues in the map, including residue shifts along a main-chain trace. The method named DAQ computes the likelihood that each local density corresponds to different amino acids, atoms, and secondary structures from the map density distribution and assesses how well amnio acids in the reconstructed model structure agree with the likelihood. DAQ is complementary to existing model validation scores for cryo-EM that examine local density gradient in the map or stereochemical geometry of the structure model. When DAQ was applied to different versions of model structure entries in PDB that were derived from the same density maps, a clear improvement of DAQ-score was observed in the newer version of the models. We also applied DAQ to a larger dataset of protein structure models from cryo-EM and suggested potential misassignment of residues.

## Overall Protocol


## Pre-required software
Python 3 : https://www.python.org/downloads/    
Pymol(for visualization): https://pymol.org/2/   

## Installation  
### 1. [`Install git`](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) 
### 2. Clone the repository in your computer 
```
git clone git@github.com:kiharalab/DeepMainMast.git && cd DeepMainMast
```

### 3. Build dependencies.   
You have two options to install dependency on your computer:
#### 3.1 Install with pip and python(Ver 3.6.9).
##### 3.1.1[`install pip`](https://pip.pypa.io/en/stable/installing/).
##### 3.1.2  Install dependency in command line.
```
pip3 install -r requirement.txt --user
```
If you encounter any errors, you can install each library one by one:
```
pip3 install mrcfile==1.2.0
pip3 install numpy==1.19.4
pip3 install numba==0.52.0
pip3 install torch==1.6.0
pip3 install scipy==1.6.0
```

#### 3.2 Install with anaconda
##### 3.2.1 [`install conda`](https://bit.ly/daq-score). 
##### 3.2.2 Install dependency in command line
```
conda create -n deepmainamst python=3.8.5
conda activate deepmainmast
pip install -r requirement.txt 
```
Each time when you want to run my code, simply activate the environment by
```
conda activate deepmainmast
conda deactivate(If you want to exit) 
```

## Usage
```
python3 main.py -h:
  -h, --help            show this help message and exit
  -F F                  Map file path
  -M1 M1                Phase1 model path
  -M2 M2                Phase2 model path
  --mode MODE           Running Mode
  --contour CONTOUR     Map contour level
  --stride STRIDE       Stride size for scanning maps (default:1)
  --voxel_size1 VOXEL_SIZE1
                        Phase1 input voxel size
  --voxel_size2 VOXEL_SIZE2
                        Phase2 input voxel size
  --gpu GPU             specify the gpu we will use
  --batch_size BATCH_SIZE
                        batch size for inference
  --cardinality CARDINALITY
                        ResNeXt cardinality
  --workers WORKERS     number of workers for dataloader
```

### 1. Build Protein Structures with EM maps
```
python3 main.py --mode=0 -F=[Map_path] --gpu=0 --contour=[contour_level]
```
Here [Map_path] is the cryo-EM map file path in your computer. [contour_level] is optional, if not please just use default 0.
```
