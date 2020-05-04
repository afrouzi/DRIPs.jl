## A Software Package for LQG Dynamic Rational Inattention Problems (DRIPs)

<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://afrouzi.github.io/DRIPs.jl/stable) -->
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://afrouzi.github.io/DRIPs.jl/dev)
[![Build Status](https://travis-ci.com/afrouzi/DRIPs.jl.svg?branch=master)](https://travis-ci.com/afrouzi/DRIPs.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/afrouzi/DRIPs.jl?svg=true)](https://ci.appveyor.com/project/afrouzi/DRIPs-jl)
[![Codecov](https://codecov.io/gh/afrouzi/DRIPs.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/afrouzi/DRIPs.jl)
<!-- [![Coveralls](https://coveralls.io/repos/github/afrouzi/DRIPs.jl/badge.svg?branch=master)](https://coveralls.io/github/afrouzi/DRIPs.jl?branch=master)
 -->
This package provides a fast and robust method for solving LQG Dynamic Ratinoal Inattention Problems based on methods developed in [Afrouzi and Yang (2019): Dynamic Rational Inattention and the Phillips Curve](http://www.afrouzi.com/dynamic_inattention.pdf). For detailed documentation of the software, see the [documentation page](https://afrouzi.github.io/DRIPs.jl/dev).

### Installation
The package is registered in the [General](https://github.com/JuliaRegistries/General) registry and can be installed at the REPL with `] add DRIPs`.

### Resources
* [Afrouzi and Yang (2019): Dynamic Rational Inattention and the Phillips Curve](http://www.afrouzi.com/dynamic_inattention.pdf)
* [GitHub Repository for Matlab code, Jupyter notebooks and earlier versions](https://github.com/choongryulyang/dynamic_multivariate_RI)

### Examples and Replications
 See the examples folder for:
* Pricing with Rational Inattention without Feedback [[HTML]](https://afrouzi.github.io/DRIPs.jl/dev/examples/ex1_pricing_nofeedback/ex1_pricing_pe_nofeedback/) [[Jupyter notebook]](https://github.com/afrouzi/DRIPs.jl/blob/master/examples/ex1_pricing_pe_nofeedback.ipynb) [[.jl file]](https://github.com/afrouzi/DRIPs.jl/blob/master/examples/ex1_pricing_pe_nofeedback.jl)
* Pricing with Rational Inattention with Endogenous Feedback [[HTML]](https://afrouzi.github.io/DRIPs.jl/dev/examples/ex2_pricing_wfeedback/ex2_pricing_pe_with_feedback/) [[Jupyter notebook]](https://github.com/afrouzi/DRIPs.jl/blob/master/examples/ex2_pricing_pe_with_feedback.ipynb) [[.jl file]](https://github.com/afrouzi/DRIPs.jl/blob/master/examples/ex2_pricing_pe_with_feedback.jl)
* Replication of [Mackowiack and Wiederholt (2009)](https://www.aeaweb.org/articles?id=10.1257/aer.99.3.769) [[HTML]](https://afrouzi.github.io/DRIPs.jl/dev/examples/ex3_mw2009/ex3_Mackowiak_Wiederholt_2009/) [[Jupyter notebook]](https://github.com/afrouzi/DRIPs.jl/blob/master/examples/ex3_Mackowiak_Wiederholt_2009.ipynb) [[.jl file]](https://github.com/afrouzi/DRIPs.jl/blob/master/examples/ex3_Mackowiak_Wiederholt_2009.jl)
* Replication of [Sims (2011)](http://sims.princeton.edu/yftp/RIMP/handbookChapterRI2.pdf) [[HTML]](https://afrouzi.github.io/DRIPs.jl/dev/examples/ex4_sims2011/ex4_Sims_2011/) [[Jupyter notebook]](https://github.com/afrouzi/DRIPs.jl/blob/master/examples/ex4_Sims_2011.ipynb) [[.jl file]](https://github.com/afrouzi/DRIPs.jl/blob/master/examples/ex4_Sims_2011.jl)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/choongryulyang/dynamic_multivariate_RI/master) to run and modify the examples online (no software is needed on the local machine).
