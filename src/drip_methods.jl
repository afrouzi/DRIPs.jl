## A General Structure for D.R.I.P.
"""
# Summary
  A [Mutable] Structure for LQG Dynamic Rational Inattention Problems (DRIPs)

# Fields
## Primitives of the DRIP
    ω      : Cost of information
    β      : Discount factor
    A      : Transition matrix: x=Ax+Qu
    Q      : Std. Dev. matrix: x=Ax+Qu
    H      : Mapping of shocks to actions: v=-0.5(a'-x'H)(a-H'x)

## Solution of the DRIP in the Steady State
    K      : Kalman gain matrix
    Y      : Weight vector for evolution of actions
    Σ_z    : Covariance matrix of the rational inattention error
    Σ_p    : Steady-state posterior covariance matrix under the solution
    Σ_1    : Steady-state prior covariance matrix under the solution
    Ω      : Dynamic benefit matrix
    err    : Convergence error for the solution
"""
mutable struct Drip
    ω; β; A; Q; H;                    # primitives
    K; Y; Σ_z; Σ_p; Σ_1; Ω;           # solution
    err;                              # convergence err
    Drip()          = new();
    Drip(ω,β,A,Q,H) = new(ω,β,A,Q,H);
end

"""
    solve_drip(ω,β,A,Q,H;
               fcap  = false,
               Ω0    = H*H',
               Σ0    = A*A'+Q*Q',
               w     = 1,
               tol   = 1e-4,
               maxit = 10000) -> Drip

Solves for the steady state of a Dynamic Rational Inattention Problem (DRIP)
    defined by the arguments. See [Afrouzi and  Yang (2019)](http://afrouzi.com/dynamic_inattention.pdf)
    for details.
# Arguments
The function takes the primitives of the Drip as arguments:

    * ω      : Cost of information
    * β      : Discount factor
    * A      : Transition matrix: x=Ax+Qu
    * Q      : Std. Dev. matrix: x=Ax+Qu
    * H      : Mapping of shocks to actions: v=-0.5(a'-x'H)(a-H'x)
## Optional Arguments
Default values are set unless specified otherwise by user.

    * fcap  = false    [if `true` then solves the problem with fixed capacity = ω bits]
    * Ω0    = H*H'     [initial guess for steady state information matrix]
    * Σ0    = A*A'+Q*Q'[initial guess for steady state prior]
    * w     = 1        [updating weight on the new guess in iteration]
    * tol   = 1e-4     [tolerance level for convergence]
    * maxit = 10000    [maximum number of iterations]
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
julia> P = solve_drip(ω,β,A,Q,H)
```
"""
function solve_drip(ω,β,A,Q,H;              # primitives of the D.R.I.P.
                    fcap::Bool = false,     # optional: if true then solves the problem with fixed capacity κ = ω.
                    Ω0         = H*H',      # optional: initial guess for steady state information matrix
                    Σ0         = A*A'+Q*Q', # optional: initial guess for steady state prior
                    w          = 1,         # optional: updating weight in iteration
                    tol        = 1e-4,      # optional: tolerance level for convergence
                    maxit      = 10000)     # optional: maximum number of iterations

    ## initialize
    (n,m) = length(size(H)) == 2 ? size(H) : (size(H,1),1);
    # n: dimension of state, m: number of actions
    eye   = Matrix{Float64}(I,n,n);
    err   = 1
    iter  = 0
    SqRΣ  = sqrt(Σ0)
    Ω_c   = H*H';
    Σ1    = Matrix{Float64}(I,n,n)
    Σ_p   = Matrix{Float64}(I,n,n)
    Λ     = Matrix{Float64}(I,n,n)
    κ     = ω
    # iterate
    while (err > tol) & (iter < maxit)
        D, U    = eigen!(Symmetric(SqRΣ*Ω0*SqRΣ));
        D       = diagm(D);
        U       = getreal(U);

        if fcap == true
            ω = (2^(2*κ)/det(max.(ω*eye,D)))^(-1/n);
        end
        Σ_p     = Symmetric(ω*SqRΣ*U/(max.(D,ω*eye))*U'*SqRΣ);

        Σ1      = Symmetric(A*Σ_p*A' + Q*Q');
        err     = norm(Σ1 - Σ0,2)/norm(Σ0,2);

        Σ0      = w*Σ1 + (1-w)*Σ0

        SqRΣ    = getreal(sqrt(Σ0));
        invSqRΣ = getreal(inv(SqRΣ));

        Ω0      = w*(Ω_c .+ β*A'*invSqRΣ*U*(min.(D,ω*eye))*U'*invSqRΣ*A)
                +(1-w)*Ω0;
        Ω0      = (abs.(Ω0).>1e-10).*Ω0;

        iter   += 1
    end

    if iter == maxit
        print("RI Code hit maxit\n")
    end
    P      = Drip(ω,β,A,Q,H);
    Σ1     = (abs.(Σ1).>1e-10).*Σ1 + diagm(abs.(diag(Σ1)).<= 1e-10)*1e-8
    inv_Σ1 = pinv(Σ1) ;
    P.Σ_1  = Σ1;
    P.Σ_p  = Σ_p;
    P.Y    = (eye - Σ_p*inv_Σ1)'*H;
    P.Σ_z  = H'*(Σ_p - Σ_p*inv_Σ1*Σ_p)*H;
    P.K    = Σ1*P.Y*pinv(P.Y'*Σ1*P.Y + P.Σ_z);
    P.Ω    = Ω0;
    P.err  = err;

    return(P)
end

"""
    solve_drip(P::Drip;...) -> Drip
Same as above but infers `ω,β,A,Q` and `H` from `P` and returns a `Drip` structure
with the primitives and the solution.
# Examples
```julia-repl
julia> P = Drip(ω,β,A,Q,H)
julia> P = solve_drip(P)
```
"""
function solve_drip(P   ::Drip;             # D.R.I.P. to be solved
                    fcap::Bool = false,     # optional: if true then solves the problem with fixed capacity κ = ω.
                    Ω0         = P.H*P.H',  # optional: initial guess for steady state information matrix
                    Σ0         = P.A*P.A'+P.Q*P.Q', # optional: initial guess for steady state prior
                    w          = 1,         # optional: updating weight in iteration
                    tol        = 1e-4,      # optional: tolerance level for convergence
                    maxit      = 10000)     # optional: maximum number of iterations
    P = solve_drip(P.ω,P.β,P.A,P.Q,P.H;
                    fcap  = fcap,
                    Ω0    = Ω0,
                    Σ0    = Σ0,
                    w     = w,
                    tol   = tol,
                    maxit = maxit);
    return(P)
end
