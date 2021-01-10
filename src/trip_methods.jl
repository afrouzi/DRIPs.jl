"""
# Summary
  A Structure for the Transition dynamics of Rational Inattention Problems (TRIPs)

# Fields
    T       : length of TRIP
    Σ_1s    : sequence of prior covariance matrices
    Σ_ps    : sequence of posterior covariance matrices
    Ωs      : sequence of information benefit matrices
    Ds      : eigenvalues of Σ_t^(0.5)Ω_tΣ_t^(0.5) over time (marginal values of information)
    err     : convergence error for shooting algorithm
    con_err : convergence error for the steady state prior 
"""
struct Trip
# Drip
    p       :: Drip
# length of T.R.I.P.
    T       :: Int32
# priors, posteriors and benefit matrices
    Σ_1s    :: Array{Float64,3}
    Σ_ps    :: Array{Float64,3}
    Ωs      :: Array{Float64,3}
# eigenvalues of Σ_t^(0.5)Ω_tΣ_t^(0.5)
    Ds      :: Array{Float64,2}
# convergence err
    err     :: Float64
    con_err :: Float64
end

"""
# Summary
  A Signal Structure for Information Treatments in DRIPs. The structure encodes
  the signal `S = L'*x+z, z₀~N(0,Σ_z)`.
# Fields

* L    : loading of the signal on `x`
* Σ_z  : variance covariance matrix of the noise
"""
struct Signal
    L   :: Array{Float64,2}
    Σ_z :: Array{Float64,2}
    function Signal(L, Σ_z)
        L, Σ_z = collect(L)[:,:], collect(Σ_z)[:,:]
        return new(L,Σ_z)
    end
end

"""
         Trip(p::Drip,              # steady state of D.R.I.P.
              Σ0::Array{Float64,2}; # initial prior matrix
              T     = 100,          # optional: guess for time until convergence to steady state
              tol   = 1e-4,         # optional: tolerance for convergence
              maxit = 1000          # optional: max iterations
              ) -> Trip
Solves for the transition dynamics of the optimal information structure starting
    from the initial prior distribution with covariance matrix `Σ0`.
    See [Afrouzi and  Yang (2019)](http://afrouzi.com/dynamic_inattention.pdf)
    for details.
# Outputs
Returns a `Trip` structure with the steady state and transition path of
    the optimal information structure.

# Examples
```julia-repl
julia> p = solve_drip(ω,β,A,Q,H)
julia> Σ0 = 0.1*p.Σ_1;
julia> Pt = solve_trip(p,Σ0);
```
"""
function Trip(p::Drip,Σ0;
              T     = 100,  # optional: time until convergence to steady state
              tol   = 1e-4, # optional: tolerance for convergence
              maxit = 1000, # optional: max iterations
              )
    Σ0 = collect(Σ0)[:,:]

    ## Initialize
    Ωs  = repeat(p.ss.Ω, inner = [1,1,T]);
    Ωsp = Ωs;

    Σ_1s  = repeat(p.ss.Σ_1, inner = [1,1,T]); Σ_1s[:,:,1] = Σ0;
    Σ_1sp = Σ_1s;
    (n,m) = size(p.H)
    # n: dimension of state, m: number of actions

    Σ_ps= repeat(p.ss.Σ_p, inner = [1,1,T]); # initialize posteriors
    Ds  = zeros(n,T);                     # initialize eigenvalues

    iter    = 0;
    err     = 1;
    con_err = 1;
    eye  = Matrix(I,n,n);
    while (err > tol) & (iter <= maxit)
        ## Given Ωs, find Sigmas using the law of motion for priors
        for i in 1:1:T-1
            SqSigma        = real.(sqrt(Σ_1s[:,:,i]));
            D, U           = eigen!(Symmetric(SqSigma*Ωs[:,:,i]*SqSigma));
            Ds[:,i]        = real.(D);
            U              = real.(U);
            D              = diagm(Ds[:,i]);
            Σ_ps[:,:,i]    = Symmetric(p.ω*SqSigma*U/(max.(D,p.ω*eye))*U'*SqSigma);
            Σ_1sp[:,:,i+1] = Symmetric(p.Q*p.Q'+p.A*Σ_ps[:,:,i]*p.A');
        end
        # Given Σ_1s, Find Omegas using the Euler equation
        for i = T-1:-1:1
            SqSigma      = sqrt(Σ_1sp[:,:,i+1]);
            invSqSigma   = real.(pinv(SqSigma));
            D, U         = eigen(SqSigma*Ωs[:,:,i+1]*SqSigma);
            D            = diagm(real.(D));
            U            = real.(U);
            Ωsp[:,:,i]   = Symmetric(p.H*p.H' .+
                           p.β*p.A'*invSqSigma*U*
                           min.(D,p.ω*eye)*U'*invSqSigma*p.A);
        end
        err    = norm(Σ_1sp-Σ_1s)/norm(Σ_1s);
        Σ_1s   = real.(Σ_1sp);
        Ωs     = real.(Ωsp);
        iter  += 1;
    end
    # Throw error if T was too short for convergence to steady state
    con_err = norm(p.Q*p.Q'.+p.A*Σ_ps[:,:,end-1]*p.A'.-p.ss.Σ_1)/norm(p.ss.Σ_1);
    # if (con_err  > tol)
    #     error("T was too short for convergence of Trip. Try larger T.");
    # end
    # Finally, store the eigenvalues of steady state
    SqSigma      = real.(sqrt(Σ_1s[:,:,end]));
    D, U         = eigen!(Symmetric(SqSigma*Ωs[:,:,end]*SqSigma));
    Ds[:,end]    = real.(D);
    # return the Trip structure
    return(Trip(p,T,Σ_1s,Σ_ps,Ωs,Ds,err,con_err))
end

## solve_Trip: this function solves for the transition dynamics given an initial signal for information treatment
"""
         Trip(p::Drip,             # steady state of D.R.I.P.
              S::Signal;            # information treatment in the steady state
              T     = 100,          # optional: guess for time until convergence to steady state
              tol   = 1e-4,         # optional: tolerance for convergence
              maxit = 1000          # optional: max iterations
              ) -> Trip
Solves for the transition dynamics of the optimal information structure starting
    from a one time treatment with a signal `S` in the steady state.
    See [Afrouzi and  Yang (2019)](http://afrouzi.com/dynamic_inattention.pdf)
    for details.
# Outputs
Returns a `Trip` structure with the steady state and transition path of
    the optimal information structure.

# Examples
```julia-repl
julia> p  = Drip(ω,β,A,Q,H)
julia> S  = Signal(L,Σ_z);
julia> pt = Trip(p,S);
```
"""
function Trip(p::Drip,         # D.R.I.P.
              S::Signal;       # initial Signal for information treatment
              T     = 100,     # optional: time until convergence to steady state
              tol   = 1e-4,    # optional: tolerance for convergence
              maxit = 1000     # optional: max iterations
              )
    K0   = p.ss.Σ_1*S.L/(S.L'*p.ss.Σ_1*S.L.+S.Σ_z); # Get prior covariance matrix from new signal
    Σ0   = p.ss.Σ_1.-K0*S.L'*p.ss.Σ_1;              # Get Kalman gain from new signal
    trip = Trip(p,Σ0;T=T,tol=tol,maxit=maxit);
    return(trip)
end
