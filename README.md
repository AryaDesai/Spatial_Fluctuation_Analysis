# Overview
This repo contains code for a biophysics project. The goal is to infer the number of bound proteins from images of cells. 

Each imaged cell has a certain 'mean' intensity, which as the name suggests is the intensity in each pixel divided by the number of pixels the cell occupies. 
We can identify the location of bound particles by measuring the variance in intensity. Areas of the cell where proteins are bound to operators will have a higher intensity than the mean, which will show up as a significantly high variance. 
Collecting data for different concentrations of proteins in cells gives us an array of different mean-variance pairs, the goal is to use this data to infer how many bound proteins there are, along with other properties of interest such as degree of cooperativity in binding and the dissociation constant.

# Types of files 
This repo has .py files which contain the logic to generate synthetic data based on a specific binding scheme, and jupyter notebooks which analyze either the synthetic data or real data from experiment.
It also contains .m files containing the matlab equivalent of the code.
