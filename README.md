# Large-scale constrained Gaussian processes for shape-restriced function estimation

This repository contains R codes and functions for implementing the large-scale constrained Gaussian processes for shape restricted function estimate based on Maatouk et al. (2023a), Maatouk et al. (2023b), Ray et al. (2020), and Maatouk and Bay (2017).

# Description of the associated R files:
1. 'all_base_functions.R' file contains all the required base functions to generate the large-scale constrained GP models, including functions like
Fast.LS (draws a very large Gaussian vector prior based on Maatouk et al. [2023b]), LS.KLE (draws a very large Gaussian vector prior based on Maatouk et al. [2023a]),
samp.WC (draws a sample based on sampling scheme developed in Wood and Chan [1994]), ESS (draws a sample using elliptical slice sampler by Murray et al. [2010]), nu.MH2 (draws posterior samples on the hyperparameters using Metropolis–Hastings and inverse Cholesky factor). One can find detailed description of each functions in the file.
2. 'all_models.R' file contains all models that implemented the large-scale constrained Gaussian processes for shape-restricted function estimation with and without hyperparameters updates.
3. 'Illustrative-examples.R' file contains the code of all numerical examples shown in the paper.

   For more details on the codes or the functions, refer to the associated R files.

# Note:
Part of this work was conducted with the support of the consortium in Applied Mathematics CIROQUO, gathering partners in technological research (BRGM, CEA, IFPEN, IRSN, Safran, Storengy) and academia (CNRS, Ecole Centrale de Lyon, Mines Saint-Etienne, University of Grenoble, University of Nice, University of Toulouse) around advanced methods for Computer Experiments. [link]( https://doi.org/10.5281/zenodo.65812).

# Authors:
Hassan Maatouk (CYU Cergy-Paris Université) and Didier Rullière (École des Mines de St-Étienne).

# Maintainer: 
Hassan Maatouk, hassan.maatouk@cyu.fr

# References
Maatouk, H. and Bay, X. (2017). "Gaussian process emulators for computer experiments with inequality constraints". Mathematical Geosciences, 49(5): 557-582. [doi](https://link.springer.com/article/10.1007/s11004-017-9673-2)

Maatouk, H. and Rullière, D. and Bay, X. (2023a). "Sampling large hyperplane‐truncated multivariate normal distributions", To appear in Computational Statistics. [doi](https://link.springer.com/article/10.1007/s00180-023-01416-7)

Maatouk, H. and Rullière, D. and Bay, X. (2023b). "Large-scale constrained Gaussian processes for shape-restricted function estimation". [preprint](https://hal.science/hal-04348962/file/LS-CGP.pdf)

Ray, P. and D., Pati and A. Bhattacharya (2020). "Efficient Bayesian shape-restricted function estimation with constrained Gaussian process priors". Statistics Computing 30, 839–853. [doi](https://link.springer.com/article/10.1007/s11222-020-09922-0) 

Wood, A. and Chan, G. (1994). "Simulation of Stationary Gaussian Processes in [0,1]^d". Journal of Computational and Graphical Statistics. [doi](https://www.jstor.org/stable/1390903)
