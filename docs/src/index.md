# Dynamic Rational Inattention Problems (DRIPs)

`DRIPs.jl` is a Julia software package that provides a fast and robust method for solving LQG Dynamic Rational Inattention models using the methods developed by [Afrouzi and Yang (2020)](http://www.afrouzi.com/dynamic_inattention.pdf).
## Installation
To add the package, execute the following in Julia REPL:

```@julia
using Pkg; Pkg.add("DRIPs");
```
To import and use the package, execute:

```julia
using DRIPs;
```

## Overview
```@contents
Pages = ["overview.md"]
Depth = 2
```
## Syntax
```@contents
Pages = ["syntax.md",]
Depth = 3
```
## Examples and Replications
```@contents
Pages = [
"examples/ex1_pricing_nofeedback/ex1_pricing_pe_nofeedback.md",
"examples/ex2_pricing_wfeedback/ex2_pricing_pe_with_feedback.md",
"examples/ex3_mw2009/ex3_Mackowiak_Wiederholt_2009.md",
"examples/ex4_sims2010/ex4_Sims_2010.md",
"examples/ex5_mmw2018/ex5_Mackowiak_Matejka_Wiederholt_2018.md",
"examples/ex6_ay2020/ex6_Afrouzi_Yang_2020.md"]
Depth = 1
```
