# Scaling_up_uncertainties_in_tree_allometry

This project introduces four different approaches for estimating the uncertainty in allometric models, 
including analytical approach, slope-intercept sampling, bootstrapping, and bayesian approach: 
(1) the analytical approach is based on the equations in Snedecor and Cochran (1989) and Yanai et al. (2010); 
(2) the slope-intercept sampling approach is based on the covariance structure between the slope and intercept; 
(3) the bootstrapping approach is based on resampling the entire number of allometric trees with replacement, 
and then refitting a regression equation to the sample; and 
(4) the bayesian approach is adapted from the paper Rejou-Mechain et al. (2017) and the R package "Biomass". See the detailed description of these four approaches in Lin et al. (2023). 

In the paper “Lin et al. (2023)”, we apply the above four approaches to four different sites with 
different tree characteristics: Plantation (Hawaii), Young tropical (Yucatán), Subtropical (Paranaense), and Semi-arid (Chaco). 
See the detailed site description in Lin et al. (2023). 
