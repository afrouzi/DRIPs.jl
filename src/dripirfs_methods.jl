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
