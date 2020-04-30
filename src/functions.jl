## solve_Drip: this function solves for the steady state of of the D.R.I.P.
## Inputs
# ω      : Cost of information
# β      : Discount factor
# A      : Transition matrix: x=Ax+Qu
# Q      : Std. Dev. matrix: x=Ax+Qu
# H      : Mapping of shocks to actions: v=-0.5(a'-x'H)(a-H'x)

## Outputs: a Drip structure with
# Σ_1    : Steady-state Prior uncertainty
# Σ_p    : Steady-state posterior uncertainty
# Λ      : Shadow matrix on the no-forgetting constraint -- Λ*(Σ-Σ_p) = 0
# Ω      : Dynamic benefit matrix
# Y      : Weight vector for evolution of actions
# Σ_z    : Covariance matrix of the rational inattention error
# K      : Kalman gain matrix

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

        Ω0      = w*(Ω_c + β*A'*invSqRΣ*U*(min.(D,ω*eye))*U'*invSqRΣ*A)+(1-w)*Ω0;
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

## solve_Trip: this function solves for the transition dynamics given an initial prior covariance matrix

"""
         solve_trip(Ss::Drip,             # steady state of D.R.I.P.
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
julia> Ss = solve_drip(ω,β,A,Q,H)
julia> Σ0 = 0.1*Ss.Σ_1;
julia> Pt = solve_trip(Ss,Σ0);
```
"""
function solve_trip(Ss::Drip,             # D.R.I.P. steady state
                    Σ0::Array{Float64};   # initial prior matrix
                    T     = 100,          # optional: time until convergence to steady state
                    tol   = 1e-4,         # optional: tolerance for convergence
                    maxit = 1000          # optional: max iterations
                    )
    ## Initialize
    Ωs  = repeat(Ss.Ω, inner = [1,1,T]);
    Ωsp = Ωs;

    Σ_1s  = repeat(Ss.Σ_1, inner = [1,1,T]); Σ_1s[:,:,1] = Σ0;
    Σ_1sp = Σ_1s;
    (n,m) = length(size(Ss.H)) == 2 ? size(Ss.H) : (size(Ss.H,1),1);
    # n: dimension of state, m: number of actions

    Σ_ps= repeat(Ss.Σ_p, inner = [1,1,T]); # initialize posteriors
    Ds  = zeros(n,T);                          # initialize eigenvalues

    iter = 0;
    err  = 1;
    eye  = Matrix(I,n,n);
    while (err > tol) & (iter <= maxit)
        ## Given Ωs, find Sigmas using the law of motion for priors
        for i in 1:1:T-1
            SqSigma        = real.(sqrt(Σ_1s[:,:,i]));
            D, U           = eigen(SqSigma*Ωs[:,:,i]*SqSigma);
            Ds[:,i]        = real.(D);
            U              = real.(U);
            D              = diagm(Ds[:,i]);
            Σ_ps[:,:,i]    = Ss.ω*SqSigma*U/(max.(D,Ss.ω*eye))*U'*SqSigma;
            Σ_1sp[:,:,i+1] = Ss.Q*Ss.Q'+Ss.A*Σ_ps[:,:,i]*Ss.A';
        end
        # Given Σ_1s, Find Omegas using the Euler equation
        for i = T-1:-1:1
            SqSigma      = sqrt(Σ_1sp[:,:,i+1]);
            invSqSigma   = real.(pinv(SqSigma));
            D, U         = eigen(SqSigma*Ωs[:,:,i+1]*SqSigma);
            D            = diagm(getreal(D));
            U            = real.(U);
            Ωsp[:,:,i]   = Ss.H*Ss.H'+Ss.β*Ss.A'*invSqSigma*U*min.(D,Ss.ω*eye)*U'*invSqSigma*Ss.A;
        end
        err    = norm(Σ_1sp-Σ_1s)/norm(Σ_1s);
        Σ_1s   = real.(Σ_1sp);
        Ωs     = real.(Ωsp);
        iter  += 1;
    end
    # Throw error if T was too short for convergence to steady state
    con_err = norm(Ss.Q*Ss.Q'+Ss.A*Σ_ps[:,:,end-1]*Ss.A'-Ss.Σ_1)/norm(Ss.Σ_1);
    if (con_err  > tol)
        error("T was too short for convergence of Trip. Try larger T.");
    end
    # Finally, store the eigenvalues of steady state
    SqSigma      = real.(sqrt(Σ_1s[:,:,end]));
    D, U         = eigen(SqSigma*Ωs[:,:,end]*SqSigma);
    Ds[:,end]    = real.(D);
    # return the Trip structure
    return(Trip(Ss,T,Σ_1s,Σ_ps,Ωs,Ds,err))
end

## solve_Trip: this function solves for the transition dynamics given an initial signal for information treatment
"""
         solve_trip(Ss::Drip,             # steady state of D.R.I.P.
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
julia> Ss = solve_drip(ω,β,A,Q,H)
julia> S  = Signal(L,Σ_z);
julia> Pt = solve_trip(Ss,S);
```
"""
function solve_trip(Ss::Drip,        # D.R.I.P.
                    S::Signal;       # initial Signal for information treatment
                    T     = 100,     # optional: time until convergence to steady state
                    tol   = 1e-4,    # optional: tolerance for convergence
                    maxit = 1000     # optional: max iterations
                    )
    K0   = Ss.Σ_1*S.L/(S.L'*Ss.Σ_1*S.L+S.Σ_z); # Get prior covariance matrix from new signal
    Σ0   = Ss.Σ_1-K0*S.L'*Ss.Σ_1;              # Get Kalman gain from new signal
    trip = solve_trip(Ss,Σ0;T=T,tol=tol,maxit=maxit);
    return(trip)
end

## get steady state irfs
"""
         dripirfs(P::Drip;  # Steady state of the DRIP
                  T = 40    # Optional: length of impulse response functions
                  ) -> Dripirfs
Returns a `Dripirfs` structure with the impulse response functions of the fundamental (`x`), beliefs (`x_hat`)
    beliefs (`x_hat`) and actions (`a`) to all the structural shocks
    **under the steady state information strucutre**. In particular, if `n` is the
    dimension of `x`, `m` is the dimension of `a` and `k` is the number of
    structural shocks, then

* `x` has dimension `n*k*T` where `x(i,j,:)` is the impulse response function
    of the `i`'th dimension of `x` to the `j`'th structural shock.
* `x_hat` has dimension `n*k*T` where `x_hat(i,j,:)` is the impulse response function
    of the agent's average belief about the `i`'th dimension of `x` to the `j`'th
    structural shock.
* `a` has dimension `m*k*T` where `a(i,j,:)` is the impulse response function
    of the `i`'th action to the `j`'th structural shock.
"""
function dripirfs(P::Drip;
                  T = 40)
    (n,m) = length(size(P.H)) == 2 ? size(P.H) : (size(P.H,1),1)
    (_,k) = length(size(P.Q)) == 2 ? size(P.Q) : (size(P.Q,1),1)
    x     = zeros(n,k,T);
    x_hat = zeros(n,k,T);
    a     = zeros(m,k,T);
    for kk in 1:k
        e_k = zeros(k,1); e_k[kk] = 1;
        for ii in 1:T
            if ii==1
                x[:,kk,ii]     = P.Q*e_k;
                x_hat[:,kk,ii] = (P.K*P.Y')*(x[:,kk,ii]);
            else
                x[:,kk,ii]     =P.A*x[:,kk,ii-1];
                x_hat[:,kk,ii] =P.A*x_hat[:,kk,ii-1]+(P.K*P.Y')*(x[:,kk,ii]-P.A*x_hat[:,kk,ii-1]);
            end
            a[:,kk,ii]  .= P.H'*x_hat[:,kk,ii];
        end
    end
    return(Dripirfs(T,x,x_hat,a))
end

## get irfs with optimal signals acquired under Pt::Trip
"""
         dripirfs(P::Trip;  # Transition dynamics of the DRIP
                  T = 40    # Optional: length of impulse response functions
                  ) -> Dripirfs
Returns a `Dripirfs` structure with the impulse response functions of the fundamental (`x`), beliefs (`x_hat`)
    beliefs (`x_hat`) and actions (`a`) to all the structural shocks
    **under the information structure implied by** `P`. In particular, if `n` is the
    dimension of `x`, `m` is the dimension of `a` and `k` is the number of
    structural shocks, then

* `x` has dimension `n*k*T` where `x(i,j,:)` is the impulse response function
    of the `i`'th dimension of `x` to the `j`'th structural shock.
* `x_hat` has dimension `n*k*T` where `x_hat(i,j,:)` is the impulse response function
    of the agent's average belief about the `i`'th dimension of `x` to the `j`'th
    structural shock.
* `a` has dimension `m*k*T` where `a(i,j,:)` is the impulse response function
    of the `i`'th action to the `j`'th structural shock.
"""
function dripirfs(Pt::Trip; T = 40)
    Ss    = Pt.P; # steady state Drip
    (n,m) = length(size(Ss.H)) == 2 ? size(Ss.H) : (size(Ss.H,1),1) # dimensions of state and actions
    (_,k) = length(size(Ss.Q)) == 2 ? size(Ss.Q) : (size(Ss.Q,1),1) # number of structural shocks
    L     = size(Pt.Ds,1) # length of Trip

    eye   = Matrix(I,n,n);
    x     = zeros(n,k,T); # initialize IRFs of state
    x_hat = zeros(n,k,T); # initialize IRFs of beliefs
    a     = zeros(m,k,T); # initialize IRFs of actions
    for kk in 1:k
        e_k = zeros(k,1); e_k[kk] = 1;
        for ii in 1:T
            if ii==1
                x[:,kk,ii]     = Ss.Q*e_k;
                x_hat[:,kk,ii] = (eye-Pt.Σ_ps[:,:,ii]*pinv(Pt.Σ_1s[:,:,ii]))*(x[:,kk,ii]);
            elseif ii <= L
                x[:,kk,ii]     =Ss.A*x[:,kk,ii-1];
                x_hat[:,kk,ii] =Ss.A*x_hat[:,kk,ii-1]+(eye-Pt.Σ_ps[:,:,ii]*pinv(Pt.Σ_1s[:,:,ii]))*(x[:,kk,ii]-Ss.A*x_hat[:,kk,ii-1]);
            else
                x[:,kk,ii]     =Ss.A*x[:,kk,ii-1];
                x_hat[:,kk,ii] =Ss.A*x_hat[:,kk,ii-1]+(eye-Pt.Σ_ps[:,:,end]*pinv(Pt.Σ_1s[:,:,end]))*(x[:,kk,ii]-Ss.A*x_hat[:,kk,ii-1]);
            end
            a[:,kk,ii]  .= Ss.H'*x_hat[:,kk,ii];
        end
    end
    return(Dripirfs(T,x,x_hat,a))
end

## get irfs with information treatment S at time zero.
"""
         dripirfs(Ss::Drip,          # Steady state of the DRIP (when treatment happens)
                  S::Signal;         # Signal for treatment at time 0
                  T = 40,            # optional: length of irfs
                  reoptimize = true, # optional: if true gives the irfs with reoptimized signals, if false with steady state signals
                  trip       = false # optional: if false solves for the optimal trip, if = P::trip then takes P as the transition dynamics after treatment
                  ) -> Dripirfs
Returns a `Dripirfs` structure with the impulse response functions of the fundamental (`x`), beliefs (`x_hat`)
    beliefs (`x_hat`) and actions (`a`) to all the structural shocks
    **under the information structure implied by a one time information treatment
    with** `S` ** in the steady state of the DRIP** `P`. In particular, if `n` is the
    dimension of `x`, `m` is the dimension of `a` and `k` is the number of
    structural shocks, then

* `x` has dimension `n*k*T` where `x(i,j,:)` is the impulse response function
    of the `i`'th dimension of `x` to the `j`'th structural shock.
* `x_hat` has dimension `n*k*T` where `x_hat(i,j,:)` is the impulse response function
    of the agent's average belief about the `i`'th dimension of `x` to the `j`'th
    structural shock.
* `a` has dimension `m*k*T` where `a(i,j,:)` is the impulse response function
    of the `i`'th action to the `j`'th structural shock.
"""
function dripirfs(Ss::Drip,          # Steady state of the DRIP (when treatment happens)
                  S::Signal;         # Signal for treatment
                  T = 40,            # optional: length of irfs
                  reoptimize = true, # if true gives the irfs with reoptimized signals, if false with steady state signals
                  trip       = false # if false solves the trip, if = P::trip then takes P as the trip
                  )

    (n,m) = length(size(Ss.H)) == 2 ? size(Ss.H) : (size(Ss.H,1),1) # dimensions of state and actions
    (_,k) = length(size(Ss.Q)) == 2 ? size(Ss.Q) : (size(Ss.Q,1),1) # number of structural shocks
    eye   = Matrix(I,n,n);

    if reoptimize == true
        if trip == false
            trip = solve_trip(Ss,S;T = T);
        end
        L     = size(trip.Ds,1) # length of Trips
    end

    x     = zeros(n,k,T); # initialize IRFs of state
    x_hat = zeros(n,k,T); # initialize IRFs of beliefs
    a     = zeros(m,k,T); # initialize IRFs of actions
    for kk in 1:k
        K0   = Ss.Σ_1*S.L/(S.L'*Ss.Σ_1*S.L+S.Σ_z); # Get prior covariance matrix from new signal
        Σ0   = Ss.Σ_1-K0*S.L'*Ss.Σ_1;
        e_k = zeros(k,1); e_k[kk] = 1;
        for ii in 1:T
            if ii==1
                x[:,kk,ii]     = Ss.Q*e_k;
                x_hat_p        = K0*S.L'*(x[:,kk,ii]);
                if reoptimize == true
                    x_hat[:,kk,ii] = x_hat_p + (eye-trip.Σ_ps[:,:,ii]*pinv(trip.Σ_1s[:,:,ii]))*(x[:,kk,ii]-x_hat_p);
                else
                    K0 = Σ0*Ss.Y/(Ss.Y'*Σ0*Ss.Y+Ss.Σ_z);
                    Σ0 = Ss.Q*Ss.Q' + Ss.A*(Σ0-K0*Ss.Y'*Σ0)*Ss.A';
                    x_hat[:,kk,ii] = x_hat_p + (K0*Ss.Y')*(x[:,kk,ii]-x_hat_p);
                end
            else
                x[:,kk,ii]     =Ss.A*x[:,kk,ii-1];
                if reoptimize  == true
                    if ii <= L
                        x_hat[:,kk,ii] =Ss.A*x_hat[:,kk,ii-1]+(eye-trip.Σ_ps[:,:,ii]*pinv(trip.Σ_1s[:,:,ii]))*(x[:,kk,ii]-Ss.A*x_hat[:,kk,ii-1]);
                    else
                        x_hat[:,kk,ii] =Ss.A*x_hat[:,kk,ii-1]+(eye-trip.Σ_ps[:,:,end]*pinv(trip.Σ_1s[:,:,end]))*(x[:,kk,ii]-Ss.A*x_hat[:,kk,ii-1]);
                    end
                else
                    K0 = Σ0*Ss.Y/(Ss.Y'*Σ0*Ss.Y+Ss.Σ_z);
                    Σ0 = Ss.Q*Ss.Q' + Ss.A*(Σ0-K0*Ss.Y'*Σ0)*Ss.A';
                    x_hat[:,kk,ii] =Ss.A*x_hat[:,kk,ii-1]+(K0*Ss.Y')*(x[:,kk,ii]-Ss.A*x_hat[:,kk,ii-1]);
                end
            end
            a[:,kk,ii]  .= Ss.H'*x_hat[:,kk,ii];
        end
    end
    return(Dripirfs(T,x,x_hat,a))
end


########## Aux. Functions ############

"""
    getreal(M)
Returns the real part of `M` (Same as `real.(M)`).
"""
function getreal(M)
    #if maximum(abs.(imag.(M))) > 1e-6
    #    print("Warning: Matrix decomposition returned complex numbers larger than 1e-6")
    #end
    return(real.(M))
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

## calculate the amount of information acquired in bits/nats for a Drip
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
        κ = 0.5*log(det(P.Σ_1)/det(P.Σ_p))/log(2); #returns capacity in bits
    elseif unit == "nat"
        κ = 0.5*log(det(P.Σ_1)/det(P.Σ_p));        #returns capacity in nats
    else
        println("Invalid input for unit! Capacity is reported in bits.")
        κ = 0.5*log(det(P.Σ_1)/det(P.Σ_p))/log(2);
    end
    return(κ)
end
