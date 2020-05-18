"""
# A Package for Solving Dynamic Rational Inattention Problems (DRIPs)
Copyright © 2020 Hassan Afrouzi, Choongryul Yang and Miguel Acosta

When using this package for your work, cite
[Afrouzi and Yang (2019)](http://www.afrouzi.com/dynamic_inattention.pdf)
for reference.
"""
module DRIPs

using LinearAlgebra, Random

include("drip_methods.jl")
include("trip_methods.jl")
include("irfs_methods.jl")
include("simulations_methods.jl")
include("aux_funcs.jl")

export Drip, Trip, irfs, simulate

end # module
