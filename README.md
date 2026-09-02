# RoyleYoung
Royle Young models for area or transect searches where you detect individuals (not their sign). Density and within home range resource selection covariates.  Data simulators and files to fit models in nimble. Known and latent ID versions. 
Latent ID versions consider 1) only spatial information, 2) additional Poisson pairwise match scores, 3) zero-inflated negative binomial pairwise match scored. 
These models use N-prior data augmentation: https://github.com/benaug/SCR-N-Prior-Data-Augmentation

For a similar model for sign transect or area searches with a Poisson observation model, see here;
https://github.com/benaug/SCR_Sign_Search

The within home range space use model is the same as the between primary period relocation model used in the Jolly-Seber N-Prior 
Data Augmentation repository. It is an RSF weighted BVN movement model. The factored representation used to fit the model 
here is the same as used in Jolly-Seber and is described in the movement model section of the readme for that repository:
https://github.com/benaug/Jolly-Seber-N-Prior-DA


Royle-Young paper:
https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/07-0601.1

Royle, J. Andrew, and Kevin V. Young. "A hierarchical model for spatial capture–recapture data." Ecology 89.8 (2008): 2281-2289.