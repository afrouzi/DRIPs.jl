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

## A General Structure for Transition dynamics of Rational Inattention Problems (T.R.I.P.)
"""
# Summary
  A Structure for the Transition dynamics of Rational Inattention Problems (TRIPs)

# Fields
    P    : a DRIP structure with its primitives and steady state solution
    T    : length of TRIP
    Σ_1s : sequence of prior covariance matrices
    Σ_ps : sequence of posterior covariance matrices
    Ωs   : sequence of information benefit matrices
    Ds   : eigenvalues of Σ_t^(0.5)Ω_tΣ_t^(0.5) over time (marginal values of information)
    err  : convergence err
"""
struct Trip
    P::Drip;                # problem and solution in steady state
    T::Int;                 # length of T.R.I.P.
    Σ_1s; Σ_ps; Ωs;         # priors, posteriors and benefit matrices
    Ds;                     # eigenvalues of Σ_t^(0.5)Ω_tΣ_t^(0.5)
    err;                    # convergence err
end

## A Structure for the impulse responses of D.R.I.P. for the steady state posterior
"""
# Summary
  A Structure for the impulse response functions of DRIPs

# Fields
    T     : length of IRFs
    x     : IRFs of the fundamental shocks
    x_hat : IRFs of beliefs
    a     : IRFs of actions

In particular, if `n` is the
dimension of `x`, `m` is the dimension of `a` and `k` is the number of
structural shocks, then

* `x` has dimension `n*k*T` where `x(i,j,:)` is the impulse response function of
    the `i`'th dimension of `x` to the `j`'th structural shock.
* `x_hat` has dimension `n*k*T` where `x_hat(i,j,:)` is the impulse response
    function of the agent's average belief about the `i`'th dimension of `x` to
    the `j`'th structural shock.
* `a` has dimension `m*k*T` where `a(i,j,:)` is the impulse response function of
    the `i`'th action to the `j`'th structural shock.
"""
struct Dripirfs
    T     :: Int            # length of IRFs
    x     :: Array{Float64} # IRFs of the fundamental shocks
    x_hat :: Array{Float64} # IRFs of beliefs
    a     :: Array{Float64} # IRFs of actions
end

## A structure for additional signals for treatment
# S₀ = L'x⃗₀+z₀, z₀~N(0,Σ_z)
"""
# Summary
  A Signal Structure for Information Treatments in DRIPs. The structure encodes
  the signal `S = L'*x+z, z₀~N(0,Σ_z)`.
# Fields
    L    : loading on `x`
    Σ_z  : variance covariance matrix of the noise


"""
struct Signal
    L; Σ_z;
end
