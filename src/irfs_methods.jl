"""
# Summary
  A Structure for the irfs/simulations of DRIPs

# Fields
    T     : length of IRFs/simulation
    x     : IRFs/simulation of the fundamental shocks
    x_hat : IRFs/simulation of beliefs
    a     : IRFs/simulation of actions

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
struct Path
    T     :: Int            # length of IRFs
    x     :: Array{Float64} # IRFs of the fundamental shocks
    x_hat :: Array{Float64} # IRFs of beliefs
    a     :: Array{Float64} # IRFs of actions
end


"""
              irfs(p :: Drip;        # Steady state of the DRIP
                   T :: Int64 = 40   # Optional: length of impulse response functions
                   ) -> Path
Returns a [`Path`](@ref) structure with the impulse response functions of the fundamental (`x`), beliefs (`x_hat`)
    and actions (`a`) to all the structural shocks
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
function irfs(p :: Drip;        # Steady state of the DRIP
              T :: Int64 = 40   # Optional: length of impulse response functions
              )
    (n,m) = size(p.H)
    (_,k) = size(p.Q)
    x     = zeros(n,k,T);
    x_hat = zeros(n,k,T);
    a     = zeros(m,k,T);
    for kk in 1:k
        e_k = zeros(k,1); e_k[kk] = 1;
        for ii in 1:T
            if ii==1
                x[:,kk,ii]     = p.Q*e_k;
                x_hat[:,kk,ii] = (p.ss.K*p.ss.Y')*(x[:,kk,ii]);
            else
                x[:,kk,ii]     =p.A*x[:,kk,ii-1];
                x_hat[:,kk,ii] =p.A*x_hat[:,kk,ii-1]+(p.ss.K*p.ss.Y')*(x[:,kk,ii]-p.A*x_hat[:,kk,ii-1]);
            end
            a[:,kk,ii]  .= p.H'*x_hat[:,kk,ii];
        end
    end
    return(Path(T,x,x_hat,a))
end

## get irfs with optimal signals acquired under pt::Trip
"""
          irfs(pt :: Trip;       # Transition dynamics of the DRIP
               T  :: Int64 = 40  # Optional: length of impulse response functions
               ) -> Path
Returns a [`Path`](@ref) structure with the impulse response functions of the fundamental (`x`), beliefs (`x_hat`)
    and actions (`a`) to all the structural shocks
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
function irfs(pt :: Trip;       # Transition dynamics of the DRIP
              T  :: Int64 = 40  # Optional: length of impulse response functions
              )
    p     = pt.p; # steady state Drip
    (n,m) = size(p.H)
    (_,k) = size(p.Q)
    L     = size(pt.Ds,1) # length of Trip

    eye   = Matrix(I,n,n);
    x     = zeros(n,k,T); # initialize IRFs of state
    x_hat = zeros(n,k,T); # initialize IRFs of beliefs
    a     = zeros(m,k,T); # initialize IRFs of actions
    for kk in 1:k
        e_k = zeros(k,1); e_k[kk] = 1;
        for ii in 1:T
            if ii==1
                x[:,kk,ii]     = p.Q*e_k;
                x_hat[:,kk,ii] = (eye-pt.Σ_ps[:,:,ii]*pinv(pt.Σ_1s[:,:,ii]))*(x[:,kk,ii]);
            elseif ii <= L
                x[:,kk,ii]     =p.A*x[:,kk,ii-1];
                x_hat[:,kk,ii] =p.A*x_hat[:,kk,ii-1]+(eye-pt.Σ_ps[:,:,ii]*pinv(pt.Σ_1s[:,:,ii]))*(x[:,kk,ii]-p.A*x_hat[:,kk,ii-1]);
            else
                x[:,kk,ii]     =p.A*x[:,kk,ii-1];
                x_hat[:,kk,ii] =p.A*x_hat[:,kk,ii-1]+(eye-pt.Σ_ps[:,:,end]*pinv(pt.Σ_1s[:,:,end]))*(x[:,kk,ii]-p.A*x_hat[:,kk,ii-1]);
            end
            a[:,kk,ii]  .= p.H'*x_hat[:,kk,ii];
        end
    end
    return(Path(T,x,x_hat,a))
end

## get irfs with information treatment S at time zero.
"""
         irfs(p          :: Drip,                          # Steady state of the DRIP (when treatment happens)
              S          :: Signal;                        # Signal for treatment
              T          :: Int64 = 40,                    # optional: length of irfs
              reoptimize :: Bool = true,                   # if true gives the irfs with reoptimized signals, if false with steady state signals
              trip       :: Union{Trip, Nothing} = nothing # if false solves the trip, if = P::trip then takes P as the trip
              )
Returns a [`Path`](@ref) structure with the impulse response functions of the fundamental (`x`), beliefs (`x_hat`)
    and actions (`a`) to all the structural shocks
    **under the information structure implied by a one time information treatment
    with** `S` **in the steady state of the DRIP** `P`. In particular, if `n` is the
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
function irfs(p          :: Drip,                          # Steady state of the DRIP (when treatment happens)
              S          :: Signal;                        # Signal for treatment
              T          :: Int64 = 40,                    # optional: length of irfs
              reoptimize :: Bool = true,                   # if true gives the irfs with reoptimized signals, if false with steady state signals
              trip       :: Union{Trip, Nothing} = nothing # if false solves the trip, if = P::trip then takes P as the trip
              )

    (n,m) = size(p.H)
    (_,k) = size(p.Q)
    eye   = Matrix(I,n,n);

    if reoptimize == true
        if trip == nothing
            trip = Trip(p,S;T = T);
        end
        L     = size(trip.Ds,1) # length of Trips
    end

    x     = zeros(n,k,T); # initialize IRFs of state
    x_hat = zeros(n,k,T); # initialize IRFs of beliefs
    a     = zeros(m,k,T); # initialize IRFs of actions
    for kk in 1:k
        K0   = p.ss.Σ_1*S.L/(S.L'*p.ss.Σ_1*S.L+S.Σ_z); # Get prior covariance matrix from new signal
        Σ0   = p.ss.Σ_1-K0*S.L'*p.ss.Σ_1;
        e_k = zeros(k,1); e_k[kk] = 1;
        for ii in 1:T
            if ii==1
                x[:,kk,ii]     = p.Q*e_k;
                x_hat_p        = K0*S.L'*(x[:,kk,ii]);
                if reoptimize == true
                    x_hat[:,kk,ii] = x_hat_p + (eye-trip.Σ_ps[:,:,ii]*pinv(trip.Σ_1s[:,:,ii]))*(x[:,kk,ii]-x_hat_p);
                else
                    K0 = Σ0*p.ss.Y/(p.ss.Y'*Σ0*p.ss.Y+p.ss.Σ_z);
                    Σ0 = p.Q*p.Q' + p.A*(Σ0-K0*p.ss.Y'*Σ0)*p.A';
                    x_hat[:,kk,ii] = x_hat_p + (K0*p.ss.Y')*(x[:,kk,ii]-x_hat_p);
                end
            else
                x[:,kk,ii]     =p.A*x[:,kk,ii-1];
                if reoptimize  == true
                    if ii <= L
                        x_hat[:,kk,ii] =p.A*x_hat[:,kk,ii-1]+(eye-trip.Σ_ps[:,:,ii]*pinv(trip.Σ_1s[:,:,ii]))*(x[:,kk,ii]-p.A*x_hat[:,kk,ii-1]);
                    else
                        x_hat[:,kk,ii] =p.A*x_hat[:,kk,ii-1]+(eye-trip.Σ_ps[:,:,end]*pinv(trip.Σ_1s[:,:,end]))*(x[:,kk,ii]-p.A*x_hat[:,kk,ii-1]);
                    end
                else
                    K0 = Σ0*p.ss.Y/(p.ss.Y'*Σ0*p.ss.Y+p.ss.Σ_z);
                    Σ0 = p.Q*p.Q' + p.A*(Σ0-K0*p.ss.Y'*Σ0)*p.A';
                    x_hat[:,kk,ii] =p.A*x_hat[:,kk,ii-1]+(K0*p.ss.Y')*(x[:,kk,ii]-p.A*x_hat[:,kk,ii-1]);
                end
            end
            a[:,kk,ii]  .= p.H'*x_hat[:,kk,ii];
        end
    end
    return(Path(T,x,x_hat,a))
end
