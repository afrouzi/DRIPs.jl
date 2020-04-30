## A Software Package for LQG Dynamic Rational Inattention Problems (DRIPs)

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://afrouzi.github.io/DRIPs.jl/dev)
[![Build Status](https://travis-ci.com/afrouzi/DRIPs.jl.svg?branch=master)](https://travis-ci.com/afrouzi/DRIPs.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/afrouzi/DRIPs.jl?svg=true)](https://ci.appveyor.com/project/afrouzi/DRIPs-jl)
[![Codecov](https://codecov.io/gh/afrouzi/DRIPs.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/afrouzi/DRIPs.jl)
[![Coveralls](https://coveralls.io/repos/github/afrouzi/DRIPs.jl/badge.svg?branch=master)](https://coveralls.io/github/afrouzi/DRIPs.jl?branch=master)

This package provides a fast and robust method for solving LQG Dynamic Ratinoal Inattnetion Problems based on methods developed in [Afrouzi and Yang (2019): Dynamic Rational Inattention and the Phillips Curve](http://www.afrouzi.com/dynamic_inattention.pdf). For detailed documentation of the software, see the [documentation page](https://afrouzi.github.io/DRIPs.jl/dev).

### Installation
The package is registered in the [General](https://github.com/JuliaRegistries/General) registry and can be installed at the REPL with `] add DRIPs`.

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
