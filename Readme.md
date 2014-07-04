# equi: Equipercentile Equating R Library

*equi* is an R library for equipercentile equating in Equivalent Groups (EG), Single Group (SG) and Nonequivalent groups with Anchor Test with Chained Equating (NEAT-CE) designs.

## CHANGES in v. 0.6


**Changes and bug fixes:**

* 'kern' now it works with continous variables
* changes in 'smoothtab': while using postsmoothing=TRUE it produces smoothed CDF with N points, while N is given by 'grid' parameter (default is 1000)
* in equi* functions df parameter changes it's default value to "BIC" and uses automatic df estimation
* some minor bug fixes
* 'conttab' if asked prints probabilities instead of counts.

**New:**

* 'kernsmooth' function for kernel smoothing of distributions. It can be used for demonstrational pourposes as lower level 'kern' is used internally for postsmoothing and there is no need for using 'kernsmooth' manually.
* 'interp' function for linear interpolation with estimation of NA's by spline regression (see: nslm)
* 'nslm' function for natural cubic spline regression

**To do:**

* generic methods: print, summary, predict for 'equi', 'smoothtab' and 'nslm' classes