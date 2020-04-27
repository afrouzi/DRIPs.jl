## A Software Package for LQG Dynamic Rational Inattention Problems (DRIPs)

[![Build Status](https://travis-ci.com/afrouzi/DRIP.jl.svg?branch=master)](https://travis-ci.com/afrouzi/DRIP.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/afrouzi/DRIP.jl?svg=true)](https://ci.appveyor.com/project/afrouzi/DRIP-jl)
[![Codecov](https://codecov.io/gh/afrouzi/DRIP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/afrouzi/DRIP.jl)
[![Coveralls](https://coveralls.io/repos/github/afrouzi/DRIP.jl/badge.svg?branch=master)](https://coveralls.io/github/afrouzi/DRIP.jl?branch=master)

This package provides a fast and robust method for solving LQG Dynamic Ratinoal Inattnetion Problems based on methods developed in [Afrouzi and Yang (2019): Dynamic Rational Inattention and the Phillips Curve](http://www.afrouzi.com/dynamic_inattention.pdf). It contains:
  1.  `solve_drip`: a function for solving the **steady state distribution** of actions and shocks.
  2.  `solve_trip`: a function for solving the **transition dynamics** to the steady state distribution from an aribitrary initial distribution.
  3.  `dripirfs`: a function for generating the impulse response functions of actions and beliefs both in the **steady state** and **transition path** of the optimal information structure.
### Resources
* [Afrouzi and Yang (2019): Dynamic Rational Inattention and the Phillips Curve](http://www.afrouzi.com/dynamic_inattention.pdf)
* [GitHub Repository for Matlab code, Jupyter notebooks and earlier versions](https://github.com/choongryulyang/dynamic_multivariate_RI)

### Examples and Replications
 Examples and replications include:
* Pricing with Rational Inattention without Feedback
* Pricing with Rational Inattention with Endogenous Feedback
* Replication of [Mackowiack and Wiederholt (2009)](https://www.aeaweb.org/articles?id=10.1257/aer.99.3.769)
* Replication of [Sims (2011)](http://sims.princeton.edu/yftp/RIMP/handbookChapterRI2.pdf) 

All examples are also available as Jupyter notebooks at [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/choongryulyang/dynamic_multivariate_RI/master) (no software is needed on the local machine).
