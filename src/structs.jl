## A General Structure for D.R.I.P.
mutable struct Drip 
    ω; β; A; Q; H;                    # primitives
    K; Y; Σ_z; Σ_p; Σ_1; Ω;           # solution
    err;                              # convergence err
    Drip()          = new();
    Drip(ω,β,A,Q,H) = new(ω,β,A,Q,H);
end

## A General Structure for Transition dynamics of Rational Inattention Problems (T.R.I.P.)
struct Trip 
    P::Drip;                # problem and solution in steady state
    T::Int;                 # length of T.R.I.P.
    Σ_1s; Σ_ps; Ωs;         # priors, posteriors and benefit matrices 
    Ds;                     # eigenvalues of Σ_t^(0.5)Ω_tΣ_t^(0.5)
    err;                    # convergence err
end 

## A Structure for the impulse responses of D.R.I.P. for the steady state posterior
struct Dripirfs
    T     :: Int            # length of IRFs
    x     :: Array{Float64} # IRFs of the fundamental shocks
    x_hat :: Array{Float64} # IRFs of beliefs
    a     :: Array{Float64} # IRFs of actions
end

## A structure for additional signals for treatment 
## S₀ = L'x⃗₀+z₀, z₀~N(0,Σ_z) 
struct Signal  
    L; Σ_z;
end