# Questions about alchemical metadynamics

In this repository, we present our methods for calculating the free energy surface and the corresponding uncertainty from alchemical metadynamics. For each method, we have one or two questions about the uncertainty assessment, as elaborated in the Jupyter notebooks (`Method_1.ipynb` and `Method_2.ipynb`) in folders `Method_1` and `Method_2`.

## 1. Directory structure
- In the folder `Method_1`, we present our first method, which is the one taught in [PLUMED masterclass 21-2]((https://www.plumed.org/doc-v2.7/user-doc/html/masterclass-21-2.html)) (Exercise 9). 
- In the folder `Method_2`, we used the second method to calculate the uncertainty of the free energy surface, which was Equation 10 in the paper [Using metadynamics to explore complex free-energy landscapes](https://www.nature.com/articles/s42254-020-0153-0).
- In the folder `input_files`, we stored the simulation outputs (including `COLVAR` and `HILLS_LAMBDA`) from a 5ns 1D alchemical metadynamics, which are the inputs for our data analysis in `Method_1` and `Method_2`. In the alchemical metadynamics simulation, we defined 6 alchemical states to disappear an argon atom from water and the goal was to calculate the free energy difference and its uncertainty between the coupled and the uncoupled states using `Method_1` and `Method_2`. 

## 2. Results of free energy calculations
- As a reference, the benchmark of the free energy difference was obtained from a 5 ns expanded ensemble simulation, which was about -3.137 +/- 0.135 kT. 
- As a result of the assessment of the influence of the block size on the uncertainty, 20 ps might be a reasonable block size. 
- Using 20 ps as the block size, the free energy difference estimated by Method 1 was -3.20068 +/- 0.05761 kT.
- Using 20 ps as the block size, the free energy difference estimated by Method 2 was -3.20068 +/- 0.06523 kT.

## 3. Our questions in summary
### Method 1
- As demonstrated in `Method_1.ipynb`, Method 1 underestimated the uncertainty. We wonder if you could provide some insights into the reason for that. 
- As mentioned in the notebook, the average free energy difference calculated from 20 repetitions of alchemical metadynamics was -3.120 +/- 0.415 kT. However, we expect that 1D alchemical metadynamics should have roughly the same performance as expanded ensemble simulations of the same length, which led to a much smaller uncertainty (-3.137 +/- 0.135 kT). We wonder if our way to calculate the free energy difference is correct or if our implementation of alchemical metadynamics is problematic. 

### Method 2
- We wonder if our way of implementing Equation 10 in the paper is correct. If so, could you help us identify the reason that it underestimated the true uncertainty? Here we provide the code of Method 2, `calculate_free_energy.py` in case that you are interested in taking a look at it. 
- As mentioned `Method_2.ipynb`, our implementation of Method 2 was extremely computationally expensive compared to Method 1. We wonder if there is a way to accelerate the calculation of the cumulative bias B(s, t_fill + i * L_b), which was the bottleneck of the method. 
