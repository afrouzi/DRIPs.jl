# Dynamic Rational Inattention Problems

This package provides a fast and robust method for solving LQG Dynamic Rational Inattention models using the methods developed by [Afrouzi and Yang (2019)](http://www.afrouzi.com/dynamic_inattention.pdf).

```@contents
```
## Installation
To add the package, execute:

```@julia
using Pkg; Pkg.add("DRIPs");
```
To import and use the package, execute:

```julia
using DRIP;
```

## Steady State of Dynamic Rational Inattention Problems (DRIPs)
```@autodocs
Modules = [DRIPs]
Pages   = ["drip_methods.jl"]
```

## Transition dynamics of Rational Inattention Problems (TRIPs)
```@autodocs
Modules = [DRIPs]
Pages   = ["trip_methods.jl"]
```

## Impulse Response Functions
```@autodocs
Modules = [DRIPs]
Pages   = ["dripirfs_methods.jl"]
```

## Aux. Functions
```@autodocs
Modules = [DRIPs]
Pages   = ["aux_funcs.jl"]
```
