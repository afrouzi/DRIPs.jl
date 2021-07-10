## A General Structure for D.R.I.P.
"""
# Summary
    A type srtucture for storing the steady state solution of a Drip.

# Fields
    err    : Convergence error for the solution
    Σ_p    : Steady-state posterior covariance matrix under the solution
    Σ_1    : Steady-state prior covariance matrix under the solution
    D      : Marginal values of information
    K      : Kalman gain matrix
    Y      : Weight vector for evolution of actions
    Σ_z    : Covariance matrix of the rational inattention error
    Ω      : Dynamic benefit matrix
"""
struct SteadyState
    err::Float64
    Ω::Array{Float64,2}
    Σ_1::Array{Float64,2}
    Σ_p::Array{Float64,2}
    D::Array{Float64,1}
    K::Array{Float64,2}
    Y::Array{Float64,2}
    Σ_z::Array{Float64,2}
    SteadyState(err, Ω, Σ_1) = new(err, Ω, Σ_1)
    SteadyState(err, Ω, Σ_1, Σ_p, D, K, Y, Σ_z) = new(err, Ω, Σ_1, Σ_p, D, K, Y, Σ_z)
end

"""
# Summary
  A Type Structure for LQG Dynamic Rational Inattention Problems (DRIPs)

# Fields
## Primitives of the DRIP
    ω      : Cost of information
    β      : Discount factor
    A      : Transition matrix: x=Ax+Qu
    Q      : Std. Dev. matrix: x=Ax+Qu
    H      : Mapping of shocks to actions: v=-0.5(a'-x'H)(a-H'x)

## Solution of the DRIP in the Steady State
    ss     : Steady State Solution as a SteadyState type
    See also: [`SteadyState`](@ref)
"""
struct Drip
# primitives
    ω::Float64
    β::Float64
    A::Array{Float64,2}
    Q::Array{Float64,2}
    H::Array{Float64,2}
# steady state
    ss::SteadyState
end

"""
    Drip(ω,β,A,Q,H; kwargs...) -> Drip

Solves for the steady state of a Dynamic Rational Inattention Problem (DRIP)
    defined by the arguments and stores the solution in a Drip type.
    See [Afrouzi and  Yang (2019)](http://afrouzi.com/dynamic_inattention.pdf)
    for details.
# Arguments
The function takes the primitives of the Drip as arguments:

    * ω      : Cost of information
    * β      : Discount factor
    * A      : Transition matrix: x=Ax+Qu
    * Q      : Std. Dev. matrix: x=Ax+Qu
    * H      : Mapping of shocks to actions: v=-0.5(a'-x'H)(a-H'x)
## Optional Arguments (kwargs...)
Default values are set unless specified otherwise by user.

    * fcap  [= false]    [if `true` then solves the problem with fixed capacity = ω bits]
    * Ω0    [= H*H']     [initial guess for steady state information matrix]
    * Σ0    [= A*A'+Q*Q'][initial guess for steady state prior]
    * w     [= 1]        [updating weight on the new guess in iteration]
    * tol   [= 1e-4]     [tolerance level for convergence]
    * maxit [= 10000]    [maximum number of iterations]
# Outputs
The function returns a `Drip` structure with the primitives and the solution objects:

    * Y      : Weight vector for evolution of actions
    * Σ_z    : Covariance matrix of the rational inattention error
    * K      : Kalman gain matrix
    * Σ_1    : Steady-state prior covariance matrix under the solution
    * Σ_p    : Steady-state posterior covariance matrix under the solution
    * Ω      : Dynamic benefit matrix
# Examples
```julia-repl
julia> P = Drip(ω,β,A,Q,H)
```
"""
function Drip(ω,β,A,Q,H;         # primitives of the D.R.I.P.
              fcap=false,     # optional: if true then solves the problem with fixed capacity κ = ω.
              Ω0=H * H',      # optional: initial guess for steady state information matrix
              Σ0=A * A' + Q * Q', # optional: initial guess for steady state prior
              w=1,         # optional: updating weight in iteration
              tol=1e-4,      # optional: tolerance level for convergence
              maxit=10000,     # optional: maximum number of iterations
              fast=false      # optional: does not compute anything other than the covariance matrices 
              )

    A, Q, H = collect(A)[:,:], collect(Q)[:, :], collect(H)[:, :];

    ## initialize
    (n, m) = size(H)
    # n: dimension of state, m: number of actions
    eye   = Matrix{Float64}(I, n, n);
    err   = 1
    iter  = 0
    SqRΣ  = sqrt(Σ0);
    Σ1    = Matrix{Float64}(I, n, n)
    Ω1    = Matrix{Float64}(I, n, n)
    Σp    = Matrix{Float64}(I, n, n)
    κ     = ω
    Ω_c   = H * H';
    Σq    = Q * Q';
    # iterate
    while (err > tol) & (iter < maxit)
        D, U    = eigen!(Symmetric(SqRΣ * Ω0 * SqRΣ));
        D       = diagm(D);
        U       = real.(U);

        if fcap == true
            ω = (2^(2 * κ) / det(max.(ω * eye, D)))^(-1 / n);
        end
        Σp      = Symmetric(ω * SqRΣ * U / (max.(D, ω * eye)) * U' * SqRΣ);

        SqRΣ    = real.(sqrt(Σ0));
        invSqRΣ = real.(inv(SqRΣ));

        Σ1      = Symmetric(A * Σp * A' + Σq);
        Ω1      = Ω_c .+ β * A' * invSqRΣ * U * (min.(D, ω * eye)) * U' * invSqRΣ * A
        err     = norm(Σ1 - Σ0, 1) / norm(Σ0, 1) + norm(Ω1 - Ω0, 1) / norm(Ω0, 1);

        Σ0      = w * Σ1 + (1 - w) * Σ0;
        Ω0      = w * Ω1 + (1 - w) * Ω0;
        Ω0      = (abs.(Ω0) .> 1e-10) .* Ω0;

        iter   += 1
    end

    if iter == maxit
        print("RI Code hit maxit\n")
    end
    Σ_1  = collect(Σ1)[:,:];
    Ω    = collect(Ω1)[:,:];
    if fast == true
        return(Drip(ω, β, A, Q, H, SteadyState(err, Ω, Σ_1)))
    else
        Σ_p  = collect(Σp)[:,:];
        inv_Σ1 = inv(Σ1) ;
        Y    = collect((eye - Σ_p * inv_Σ1)' * H)[:,:];
        Σ_z  = collect(H' * (Σ_p - Σ_p * inv_Σ1 * Σ_p) * H)[:,:];
        D    = collect(diag(D));
        K    = collect(Σ1 * Y * inv(Y' * Σ1 * Y .+ Σ_z))[:,:];
        return(Drip(ω, β, A, Q, H, SteadyState(err, Ω, Σ_1, Σ_p, D, K, Y, Σ_z)))
    end
end
