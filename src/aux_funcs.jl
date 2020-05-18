"""
        capacity(P::Drip;      # Drip structure
                 unit = "bit"  # optional: unit of capacity (bit or nat).
                 )
Returns the amount of information processes per unit of time in the steady state
of the DRIP `P`.
"""
function capacity(P::Drip;      # Drip structure
                  unit = "bit"  # optional: unit of capacity (bit or nat).
                  )
    if unit == "bit"
        κ = 0.5*log(det(P.ss.Σ_1)/det(P.ss.Σ_p))/log(2); #returns capacity in bits
    elseif unit == "nat"
        κ = 0.5*log(det(P.ss.Σ_1)/det(P.ss.Σ_p));        #returns capacity in nats
    else
        println("Invalid input for unit! Capacity is reported in bits.")
        κ = 0.5*log(det(P.ss.Σ_1)/det(P.ss.Σ_p))/log(2);
    end
    return(κ)
end

"""
    infinitesum(func; tol = 1e-6,maxit = 1000,start=0)
Returns the infinite sum `Σₓfunc(x)` starting from `x = start` up to tolderance
`tol` or max iteration `maxit`.
"""
function infinitesum(func; tol = 1e-6,maxit = 1000,start=0)
    diff  = 1.0
    infsum = func(start)
    it    = start + 1
    while (diff > tol) & (it < maxit)
        func_it = func(it)
        infsum += func_it
        diff = maximum(func_it)
        it += 1
    end
    return(infsum)
end
