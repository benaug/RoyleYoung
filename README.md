# Royle-Young

Royle-Young models for area or transect searches where individuals, rather than their sign, are detected. These models allow density and within-home-range resource-selection covariates. Data simulators and files for fitting the models in NIMBLE are provided for both known and latent individual identity.

The latent-ID versions consider three sources of identity information:

1. spatial information only;
2. additional Poisson pairwise match scores; and
3. zero-inflated negative binomial pairwise match scores.

The match scores here are counts, which might be the number of match keypoints/features from HotSpotter (SIFT) or a similar algorithm. Other match score distributions could be used for other data types.

These models use N-prior data augmentation:

https://github.com/benaug/SCR-N-Prior-Data-Augmentation

For a related model for sign transect or area searches with a Poisson observation model, see:

https://github.com/benaug/SCR_Sign_Search

The within home range space use model is closely related to the between primary period relocation model used in the Jolly-Seber N-Prior Data Augmentation repository:

https://github.com/benaug/Jolly-Seber-N-Prior-DA

The original Royle-Young model is described in:

Royle, J. A., and K. V. Young. 2008. A hierarchical model for spatial capture-recapture data. *Ecology* 89:2281-2289.

https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/07-0601.1

## Within home range space use

The observation model can be viewed as an extension of the spatial capture-recapture model of Royle and Young (2008).
In the original Royle-Young formulation, an individual's location during a sampling occasion is modeled as a realization 
from a bivariate normal distribution centered on its activity center. Detection is possible only when that realized location
falls within the observation window (surveyed area), and an individual within the observation window is detected with a fixed
probability $p$.
Therefore, unlike conventional SCR models in which detection probability is specified directly as a function of distance between
an activity center and a detector, the Royle-Young model explicitly represents the individual's realized continuous location
during sampling and separates the probability of being available for detection within the searched area from the probability of detection
given availability. An interesting consequence of this model is that detection probability is exactly proportional to space use.

The implementation here extends this idea from a single observation window to a regular grid of observation windows covering
the searched area. On each sampling occasion, an individual occupies one realized continuous location, denoted
$\mathbf{u} _{ik}$, and therefore can occur in only one grid cell. The underlying availability distribution for
$\mathbf{u} _{ik}$ is a bivariate normal distribution centered on the individual's activity center $\mathbf{s}_i$.
The model additionally allows spatial density covariates and within home range resource selection through a 
cell-specific resource selection function (RSF). Conditional on the selected cell, the realized location remains
continuous and follows the bivariate normal distribution truncated to the boundaries of that cell.

It is useful to distinguish three components of the within home range space use model. The availability
distribution describes where an individual could be located during a survey based only on its activity center and
the spatial scale of movement, and is given here
by the bivariate normal kernel centered on $\mathbf{s}_i$. The resource selection function describes the relative
preference for locations with different habitat characteristics. Then, the use distribution is proportional to the product
of the availability distribution and the RSF and is normalized over the spatial domain to sum to one.
Thus, $\sigma$ controls the spatial pattern of availability around the activity center, whereas
$\boldsymbol{\beta}^{\mathrm{RSF}}$ controls how that availability is redistributed among locations according to
habitat.

<p align="center">
  <img src="TripPlot.png" width="900"><br>
  <em>Example availability distribution, resource-selection function, and resulting within-home-range use distribution.</em>
</p>

### Space-use formulation

For individual $i$, let $\mathbf{s}_i$ denote its activity center and let cell $c$ have resource
selection weight

$$
r_c = \exp(\mathbf{x}_c^\top\boldsymbol{\beta}^{\mathrm{RSF}}).
$$

The availability of cell $c$ is obtained by integrating the bivariate normal kernel centered on $\mathbf{s}_i$ over
the area of the cell. Denote this cell probability by

$$
a_{ic} = P(\mathbf{u} _{ik} \in c \mid \mathbf{s}_i,\sigma),
$$

where $\sigma$ controls the spatial scale of availability around the activity center. Resource selection then reweights this cell availability
distribution, giving the probability that the realized location occurs in cell $c$,

$$
P(\mathbf{u} _{ik}\in c\mid\mathbf{s}_i,\sigma,
\boldsymbol{\beta}^{\mathrm{RSF}}) =
\frac{r_c a_{ic}}
{\sum_{h=1}^{C}r_h a_{ih}}.
$$

Conditional on $\mathbf{u} _{ik}$ occurring in cell $c$, its continuous location within that cell follows the bivariate
normal availability distribution truncated to the cell boundaries,

$$
f(\mathbf{u} _{ik}\mid\mathbf{u} _{ik}\in c,\mathbf{s}_i,\sigma) =
\frac{
f_{\mathrm{BVN}}(\mathbf{u} _{ik}\mid\mathbf{s}_i,\sigma)
}{
a_{ic}
},
\qquad \mathbf{u} _{ik}\in c.
$$

Combining the RSF-weighted probability of using cell $c$ with the conditional within cell distribution gives
the continuous space use density

$$
f(\mathbf{u} _{ik}\mid\mathbf{s}_i,\sigma,
\boldsymbol{\beta}^{\mathrm{RSF}}) =
\frac{
r_c f_{\mathrm{BVN}}(\mathbf{u} _{ik}\mid\mathbf{s}_i,\sigma)
}{
\sum_{h=1}^{C}r_h a_{ih}
},
\qquad \mathbf{u} _{ik}\in c.
$$

Thus, although habitat selection is represented at the grid-cell level, the resulting space-use distribution remains
continuous within cells; the grid determines only where the RSF weight changes.

A direct implementation of the Royle-Young formulation would retain a latent realized location $\mathbf{u} _{ik}$
for every individual and sampling occasion. These latent locations can lead to poor MCMC mixing because each
$\mathbf{u} _{ik}$ must be updated while remaining consistent with the individual's activity center. Similarly,
the activity centers are updated conditional on the current latent locations. In this situation, both the activity centers and 
latent locations can only be moved very small distances on each update, requiring a potentially very large number of MCMC 
iterations to achieve an acceptable posterior effective sample size.

The approach used here therefore integrates unobserved $\mathbf{u} _{ik}$ values out of the likelihood. For a detection,
$\mathbf{u} _{ik}$ is observed as the individual's continuous location on that occasion and contributes directly to the
continuous space use likelihood. For a nondetection, $\mathbf{u} _{ik}$ is latent and the likelihood integrates over all
possible continuous locations at which the individual could have occurred without being detected. Because the RSF weight
and detection process are defined at the cell level, this continuous space integral can be evaluated using the integrated
bivariate normal probability for each cell. The models here therefore do not retain or update latent $\mathbf{u} _{ik}$
values for nondetections, which removes the strong dependence between these latent locations and the activity centers and
substantially improves MCMC mixing and effective sample size per unit time.


### Factored representation

A direct implementation of the marginal model could calculate and store $a_{ic}$ for every individual and every two-dimensional
habitat cell. This becomes expensive because the same bivariate normal cell probabilities must be repeatedly calculated during
MCMC. For the isotropic bivariate normal kernel used here, the $x$ and $y$ dimensions are independent. Thus, the probability mass
of a rectangular grid cell factors exactly into one-dimensional probabilities. If cell $c$ corresponds to x-cell $j(c)$
and y-cell $k(c)$,

$$
a_{ic} = a^x_{i,j(c)}a^y_{i,k(c)},
$$

where

$$
a^x_{ij} = P(x\in\text{x-cell }j\mid s_{i,x},\sigma)
$$

and

$$
a^y_{ik} = P(y\in\text{y-cell }k\mid s_{i,y},\sigma).
$$

Each one-dimensional probability is calculated from differences in normal CDF values at the corresponding cell boundaries.
Instead of storing a vector of length equal to the total number of two-dimensional cells for each individual, the fitted model
therefore stores only the two one-dimensional vectors `avail.x` and `avail.y`.

Using these factored one-dimensional availability probabilities, the RSF-weighted normalizing constant for individual $i$ is

$$
D_i = \sum_{c=1}^{C}r_c a^x_{i,j(c)} a^y_{i,k(c)}.
$$

The probability that the continuous realized location occurs in any particular cell can be reconstructed when needed as

$$
P(\mathbf{u}_{ik}\in c\mid\mathbf{s}_i) =
\frac{r_c a^x_{i,j(c)} a^y_{i,k(c)}}{D_i}.
$$

For a detected individual, the exact continuous location $\mathbf{u}_{ik}$ is observed. Within its cell, the continuous
bivariate normal density factors into independent x- and y-normal densities. When the cell-use probability and the
conditional within cell density are combined, the cell availability term cancels with the truncation normalization.
The continuous likelihood contribution for a detected location in cell $c$ can therefore be written as

$$
f(\mathbf{u}_{ik}\mid\mathbf{s}_i) =
\frac{r_c \phi(u_{ik,x}\mid s_{i,x},\sigma)\phi(u_{ik,y}\mid s_{i,y},\sigma)}{D_i},
$$

where $\phi(\cdot)$ denotes a univariate normal density.

For a nondetection, the continuous location $\mathbf{u}_{ik}$ is integrated out. Because the detection probability is
constant within an observation cell, this integral depends on the RSF-weighted probability mass of the relevant cells rather
than on the unknown exact within cell location. Thus, both detected and nondetected observations can be evaluated from the
same factored availability representation without retaining latent continuous locations for nondetections.

The factorization itself is exact for a rectangular grid and independent x- and y-components of the bivariate normal model.
A different approach is required for availability distributions that do not factor into independent x- and y-components.

#### Trimming and grid resolution

For computational efficiency, the one-dimensional availability distributions are trimmed before the normal CDF differences
are evaluated. For a user-specified trimming value `avail.z`, cell probabilities are calculated only for cells intersecting
the interval $s_{i,d} \pm \texttt{avail.z}\,\sigma$ for dimension $d$, with probabilities for cells outside this interval set
to zero. Larger values of `avail.z` retain more of the Gaussian tails and give a closer approximation to the untrimmed
availability distribution, whereas smaller values reduce computation by evaluating fewer cells. Thus, `avail.z` controls
a user-adjustable tradeoff between numerical approximation and computational efficiency.

This approach is similar to the local evaluation strategy of Turek et al. (2021), who restrict SCR detection calculations
to detectors within a fixed maximum distance, $d_{\max}$, of an individual's activity center and precompute the relevant local
detectors for each habitat cell. Here, rather than using a fixed absolute distance and a precomputed spatial lookup table,
the limits are defined by a fixed number of standard deviations from the activity center. The absolute trimming distance
therefore changes with the current value of $\sigma$ and can also vary among individuals when individual-specific space use
scales are modeled. Because the habitat grid is regular, the range of cells falling within these limits can be determined
directly from the activity center, $\sigma$, and `avail.z`, without storing a lookup table.

The grid resolution, `res`, introduces a second computational tradeoff. The Gaussian availability probability for each grid cell is
calculated exactly from normal CDF differences, conditional on the chosen grid, so reducing `res` is not required to approximate
the Gaussian integral itself. Rather, `res` determines the spatial resolution at which RSF weights, search effort, density
covariates, and other cell-level quantities can vary. Finer grids provide a closer approximation to continuously varying
spatial processes but increase the number of cells involved in the space-use calculations. In applications considered here,
a cell width of approximately $\sigma$ has provided a good compromise between spatial resolution and computational cost.
This should be treated as a practical starting point rather than a general requirement. When $\sigma$ varies among individuals,
grid resolution should be considered relative to the smallest value of $\sigma$ that is substantively important.

The fitted models therefore maintain three quantities for each active individual:

- `avail.x`: one-dimensional x availability probabilities;
- `avail.y`: one-dimensional y availability probabilities; and
- `use.denom`: the RSF-weighted normalizing constant $D_i$.

When an activity center or $\sigma$ changes, the appropriate one-dimensional availability distributions and normalizing
constant are updated. When the RSF coefficient changes, the availability distributions do not change, so only the RSF-weighted
normalizing constant needs to be recomputed. This separation substantially reduces repeated normal CDF calculations during MCMC.

## Gating activity centers and space use by population membership

Unlike typical data augmentation approaches, the activity center and within home range space use components are gated by the
data augmentation indicator $z_i$. When $z_i=0$, the individual's activity center, one-dimensional availability vectors, and
use normalizing constant are set to zero and contribute no likelihood. When an augmented individual is proposed to be turned on,
its activity center and associated space use quantities are proposed jointly with $z_i$. This avoids maintaining and updating
spatial latent variables for individuals that are not currently part of the population, reducing unnecessary computation and
improving MCMC efficiency. A secondary benefit of this approach is that reversible jump MCMC can be used more reliably for density covariates because the
locations of individuals not currently in the population are not considered when proposing to turn spatial predictors on or off.


## Joint abundance and density intercept update

These models use N-prior data augmentation rather than conventional binomial data augmentation, with abundance updates
that add or subtract one individual at a time. In sparse data settings, particularly with large $N$, strong posterior dependence
between realized abundance and the density intercept can substantially reduce mixing. I therefore include an optional joint
$N/z/D_0$ update designed to move along this posterior ridge.

Expected abundance is determined by the spatial density model. For example, with a single density covariate,

$$
D_c = D_0\exp(\beta^{D}x_c),
$$

and expected abundance over the state space is

$$
\lambda = D_0\sum_{c=1}^{C}A_c I_c\exp(\beta^{D}x_c),
$$

where $A_c$ is cell area and $I_c$ indicates whether the cell is part of the state space. Realized abundance then follows

$$
N\sim\mathrm{Poisson}(\lambda).
$$

As a result, $N$ and the density intercept $D_0$ can have strongly correlated posteriors. When the conventional birth or death
update changes $N$, it leaves $D_0$, and therefore expected abundance, unchanged. This can constrain movement of $N$ because
proposed changes in realized abundance are not accompanied by corresponding changes in the density model. The problem can be
particularly severe when the posterior correlation between $N$ and $D_0$ is strong.

The optional joint update addresses this by changing $N$ and $D_0$ together. When a birth is proposed, $D_0$ is also increased,
and when a death is proposed, $D_0$ is decreased. This allows $N$ to move more freely while remaining consistent with the
expected abundance implied by the density model. For a proposed birth,

$$
N' = N+1
$$

and

$$
D_0' = D_0\frac{N+1}{N}.
$$

For a proposed death,

$$
N' = N-1
$$

and

$$
D_0' = D_0\frac{N-1}{N}.
$$

Because $\lambda$ is proportional to $D_0$, this produces

$$
\lambda' = \lambda\frac{N'}{N},
$$

and therefore preserves

$$
\frac{\lambda'}{N'} = \frac{\lambda}{N}.
$$

The proposal thus moves along the principal relationship between expected and realized abundance rather than changing $N$
while holding expected abundance fixed.

Because the proposal for $D_0$ is a deterministic multiplicative transformation, the Metropolis-Hastings ratio includes the
corresponding Jacobian,

$$
\left|\frac{\partial D_0'}{\partial D_0}\right| = \frac{N'}{N}.
$$

Thus the log acceptance ratio includes

$$
\log\left(\frac{N'}{N}\right)
$$

in addition to the usual target density and birth/death proposal terms.

The equivalent update can be written on the log density scale. If

$$
D_{\beta0}=\log(D_0),
$$

then

$$
D_{\beta0}' = D_{\beta0} + \log\left(\frac{N'}{N}\right).
$$

On this scale the deterministic transformation is a translation and its Jacobian is one. The two updates therefore describe
the same ridge move expressed on different parameter scales. This deterministic joint proposal is similar to centered proposal
constructions for improving between-state moves in reversible jump MCMC (Brooks et al. 2003). In the application that motivated
this approach, the effective sample size for realized abundance approximately doubled.

The joint update does not need to replace all conventional $N/z$ updates. In models where the relationship between $N$ and
the density intercept is not extremely strong, for example when density covariates introduce substantial posterior dependence
between the density intercept and slope, a mixture of ordinary $N/z$ proposals and joint $N/z/D_0$ proposals can provide 
better overall mixing.

Brooks, S. P., P. Giudici, and G. O. Roberts. 2003. Efficient construction of reversible jump Markov chain Monte Carlo proposal
distributions. *Journal of the Royal Statistical Society: Series B* 65:3-39.