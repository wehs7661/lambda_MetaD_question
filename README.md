# Questions about alchemical metadynamics
In this repostiroy, I've described two problems we encountered before wrapping up the paper. Below is a summary for the questions.

## Problem 1
### System of interest
A molecule (solvated in water) composed of 4 vdW sites with 0 net charges. In this simple system (solvated in water), the slowest degree of freedom is the dihedral 1-2-3-4. The force constant of the dihedral was greatly increased such that sampling only in the alchemical space (e.g. expanded ensemble or 1D alchemical metadynamics) can not estimate the free energy difference between the coupled and uncoupled states accurately. (This is shown by the difference in the estimated free energy differences obtained from expanded ensemble simulations starting from different torsional state of the molecule.)
### Description of the problems
Our hypothesis is that in a 2D alchemical metadynamics simulation where the dihedral is also biased, the free energy difference can be recovered correctly due to sufficient sampling in the configurational space. Therefore, we performed a 2D alchemical metadynamics starting from each of the torsional states of the molecule. However, it was shown that the two estimated values were inconsistent to each other, as shown below:
- 2D alchemical metadynamics starting from state A (dihedral around 180 degrees): -2.518 +/- 0.043 kT
- 2D alchemical metadynamics starting from state B (dihedral around 0 degrees): -4.994 +/- 0.032 kT

### Relevant folder
Note that some large simulation outputs can be downloaded from [this link](https://drive.google.com/drive/folders/19mCLDtWa1L9jtyh13_DHYnLnhJqpXaN8?usp=sharing).
```
Problem_1
├── analysis_results
├── Problem_1.ipynb
├── state_A
└── state_B
```

## Problem 2
### System of interest
CB8-G3 host-guest binding complex (solvated in water) from SAMPL6 SAMPLing challenge. The slowest degree of freedom of this binding complex is assumed to be the water molecules entering and exiting the binding cavity. The corresponding CV is the number is therefore the number of water molecules in the binding cavity. 

### Description of the problems
- The problems have been summarized in the notebook `Problem_2.ipynb`


### Relevant folder
Note that some large simulation outputs can be downloaded from [this link](https://drive.google.com/drive/folders/19mCLDtWa1L9jtyh13_DHYnLnhJqpXaN8?usp=sharing).
```
Problem_2
├── CV_hist.png
├── plumed.dat
├── Problem_2.ipynb
└── time_series.png
```
