"""
# A Package for Solving Dynamic Rational Inattention Problems (DRIPs)
Copyright © 2020 Hassan Afrouzi, Choongryul Yang and Miguel Acosta

When using this package for your work, cite
[Afrouzi and Yang (2019)](http://www.afrouzi.com/dynamic_inattention.pdf)
for reference.
"""
module DRIPs

using LinearAlgebra

include("structs.jl")
include("functions.jl")

export Drip, Trip, Dripirfs, Signal,
	   solve_drip,
	   solve_trip,
	   dripirfs

end # module
