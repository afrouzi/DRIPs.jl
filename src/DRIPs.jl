## A Software Package for Solving Dynamic Rational Inattention Problems (D.R.I.P.)
## Copyright © 2020 Hassan Afrouzi, Choongryul Yang and Miguel Acosta
## When using this package for your work, cite Afrouzi and Yang (2019) for reference.
## Julia code written by Miguel Acosta / Modified by Hassan Afrouzi

module DRIPs

using LinearAlgebra

include("structs.jl")
include("functions.jl")

export Drip, Trip, Dripirfs, Signal,
	   solve_drip,
	   solve_trip,
	   dripirfs 
		
end # module
