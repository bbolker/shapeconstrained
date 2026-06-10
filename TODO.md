## to do

##  cosmetic

* tweak RF fits pic: leave out RTMB_holling, rename RTMB2_mpd .. adjust _ to / 
* turn off SD prediction to speed up redeye_odo fits
* improve reflist CSL, e.g. "Pya and Wood" vs "Pya et al"
* pix, natural history background for reed frogs/red-eyed tree frogs (pred, prey)
* improve captions for plotly figs
* troubleshoot semimech fit??
* McMaster logo?
* improve README

## technical notes/future to-do

* instability of Laplace approx for small, noisy data sets. `inner.control` helps a little bit ...
* possibly related to rank-deficiency of penalty matrix (== precision matrix)?  AFAICT the Wood 2004 appendix A trick to split off the rank-deficient components is incompatible with the positivity constraints (i.e. we would like to do Laplace approximation only over the 'random' components but ...)
* GCV vs ML/REML selection
* are there other Wood/Pya computational tricks we can use in an RTMB/TMB context? 
* regularization?
* `scam`: better knot selection (e.g. based on *unique* x values)? fix binomial/N>1 bug?
* experiment with constrained opt. vs exponentiation
* bases/constraints corresponding to unimodal smooths?
* fart around with optimization alternatives, e.g. `DEoptim` ???
* `tmbstan` (with priors/regularization)?

## non-cosmetic

* ecdf/AIC for RTMB? (test with `scam_binom_test.R`, `scam_mpd_test.R` examples)

???

## reed frogs

* RTMB_mpd: penalization too strong? Is this (RE)ML vs GCV?
 * compare `m_scam_mpd[c("trA", "aic", "sp", "edf")]`
 * works OK for binom_test (better data!)
 
* get lambda, ecdf from RTMB_mpd, scam_mpd
   * `sp` is multiplied by 
* compare AIC values?
* tmbstan
* constrained optimization?
* bad knots? too many, wrong place? (hard to adjust ...)

## waterbugs

* fit!
* pix
* AIC tables for various parametric fits
* scam, RTMB (semimech)?? one smooth for mpd, 


larvae of red-eyed treefrogs Agalychnis callidryas and two species of aquatic invertebrate predators, adult predatory water bugs (Belostoma sp. Belostomatidae) and dragonfly nymphs (Pantala flavescens Libellulidae)


## other

* hypothesis tests (goodness-of-fit/AIC, monotonicity, concavity, etc.)
* dynamical sensitivity???

## technical bits

* plotting? `gratia`?
* figure out weights for scam: are weights precision weights or analytic weights (sensu Lumley https://notstatschat.rbind.io/2020/08/04/weights-in-statistics/) ?

  are weights invariant under scaling (precision weights) or not (frequency weights)?

## waterbugs

* what did we do originally?
* functions for attack(size): Ricker, power-Ricker, logistic, hyperbolic, exponential,
* functions for handling(size): exponential, linear, proportional, independent

prob = 1/((1/a) + h*init_dens)

* gam, scam, glmmTMB, RTMB, JAGS ?
*  x simulated data, reedfrog, McCoy data
*  x Gaussian, binomial
*  x GAM, scGAM, parametric

regularization??

## talk outline

* semipar models
* Levins
* shape-constrained models
* `gam`, `scam`
* bases and constraints

* list of options?
* test 

scam smooth codes:
*  m = monotonic
*  p = p-spline
*  i/d = increasing/decreasing
*  cv = concavity

unimodal splines? (uniReg package)
