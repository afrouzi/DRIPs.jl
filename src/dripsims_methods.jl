"""
         dripsims(P::Drip;      # Steady state of the DRIP
                  T    = 500,   # Optional: length of simulation
                  N    = 1,     # Optional: number of simulated agents
                  burn = 100,   # Optional: length of initial burn (in addition to T)
                  seed = true   # Optional: seed number for fundamental shocks
                  ) -> Dripirfs
Returns a `Dripirfs` structure with a simulated path of the fundamental (`x`), beliefs (`x_hat`)
    and actions (`a`)
    **under the steady state information structure**. In particular, if `n` is the
    dimension of `x` and `m` is the dimension of `a`, then

* `x` has dimension `n*T` where `x(:,t)` is the simulated value of `x` at time `t`.
* `x_hat` has dimension `n*N*T` where `x_hat(:,i,t)` is the simulated value of `x_hat` of agent `i` at time `t`.
* `a` has dimension `m*N*T` where `a(:,i,t)` is the the simulated value of `x_hat` of agent `i` at time `t`.
"""
function dripsims(P::Drip;
                  T    = 500,
                  N    = 1,
                  burn = 100,
                  seed = true)
    (n,m)    = length(size(P.H)) == 2 ? size(P.H) : (size(P.H,1),1)
    (_,k)    = length(size(P.Q)) == 2 ? size(P.Q) : (size(P.Q,1),1)
    if seed == true
        x_Shocks = randn(k,burn+T);
    else
        Random.seed!(seed);
        x_Shocks = randn(k,burn+T);
    end
    Random.seed!();
    a_Shocks = randn(m,N,burn+T);
    x        = zeros(n,burn+T);
    x_hat    = zeros(n,N,burn+T);
    a        = zeros(m,N,burn+T);
    SqΣz     = sqrt(P.Σ_z);
    for t in 2:T+burn
        x[:,t] = P.A*x[:,t-1]+P.Q[:,:]*x_Shocks[:,t];
        for i in 1:N
            x_hat[:,i,t] = P.A*x_hat[:,i,t-1]+(P.K*P.Y')*(x[:,t]-P.A*x_hat[:,i,t-1])+P.K[:,:]*SqΣz*a_Shocks[:,i,t];
            a[:,i,t]    .= P.H'*x_hat[:,i,t];
        end
    end
    x     = x[:,burn+1:T+burn];
    x_hat = x_hat[:,:,burn+1:T+burn];
    a     = a[:,:,burn+1:T+burn];
    return(Dripirfs(T,x,x_hat,a))
end
