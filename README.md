# Modesto_LIF_3D

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/ScreenNeuroPharm/3DNeuronalModel/blob/master/LICENSE)

> The repository contains the data and the functions needed to reproduce the analysis reported in the article "Modeling the three-dimensional connectivity of in vitro cortical ensembles coupled to Micro-Electrode Arrays"


## Details
The code comprehends:
- a folder with the Python code necessary for the computational model in each of the examined configuration. The name of the folder corresponds to the name of the simulated data in the Data folder.
- a folder with the matlab code code necessary for the analysis. All uploaded scripts for the analyses work with a .mat format. To reproduce our analysis is necessary to convert the ```.txt``` format file in ```.mat``` format file using the function ```TxT2Mat.m``` in the Conversion folder. ```TxT2Mat.m``` function allows obtaining for each electrode the peak train .mat file. Peak_train file is a sparse vector that reports the spike occur.

All recordings are sampled at 10 KHz. 


### Code folder architecture:

#### Computational Model folder:
- Gaussian folder: contains the .py functions to obtain the 2D and 3D Gaussian configuration
- Gaussian_DifferentSources folder: contains the .py functions to obtain the 3D Gaussian configuration with different projecting sources
- Gaussian_modularity: contains the .py functions to obtain two 3D modules interconnected at the bottom layer
- SF-all folder:contains the .py functions to obtain the 2D and 3D Scale-Free configuration
- SF-L0 folder:contains the .py functions to obtain the configuration where the scale-free connectivity is implemented only in readout layer
- EternalFunctions: .py folder that contains supplementary functions
- requirements.txt

#### Data analysis folder:
- Conversion folder:
    * Txt2Mat: function to convert ```.txt``` format file in ```.mat``` format file
    * SplitLayers_3D_all: function to divide the peak train of the different layers in the different folders
    * SplitLayers_3D_all_modular: function to divide the peak train of the different modules and of the different layers in the different folders

- BurstAnalysis folder:
    * burstFeaturesAnalysis: functions to obtain the burst feautures
    * StringMethod: functions to detect the burst

- NetworkBurstDetection folder: 
    * IBEi: functions to extract the Inter Burst Event interval (threshold used to detect the Network Burst)
    * NetBurstDetection: functions for the Network Burst Detection
    * STH: function to compute the Time Spike Histogram
	
- Similarityfolder:
    * Similarity_test: for the extraction of the network burst initial events and the computation of their similarity
    * VP_compute_normalized_dist: function to compute the normalized victor purpura distance

- Utilities folder: supplementary functions