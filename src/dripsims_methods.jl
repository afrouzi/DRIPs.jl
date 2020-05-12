"""
         dripsims(P::Drip;      # Steady state of the DRIP
                  T    = 500,   # Optional: length of simulation
                  burn = 100,   # Optional: length of initial burn
                  seed = true   # Optional: seed number for fundamental shocks
                  ) -> Dripirfs
Returns a `Dripirfs` structure with a simulated path of the fundamental (`x`), beliefs (`x_hat`)
    and actions (`a`)
    **under the steady state information structure**. In particular, if `n` is the
    dimension of `x` and `m` is the dimension of `a`, then

* `x` has dimension `n*T` where `x(:,t)` is the simulated value of `x` at time `t`.
* `x_hat` has dimension `n*T` where `x_hat(:,t)` is the simulated value of `x_hat` at time `t`.
* `a` has dimension `m*T` where `a(:,t)` is the the simulated value of `x_hat` at time `t`.
"""
function dripsims(P::Drip;
                  T    = 500,
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
    a_Shocks = sqrt(P.Σ_z)*randn(m,burn+T);
    x        = zeros(n,burn+T);
    x_hat    = zeros(n,burn+T);
    a        = zeros(m,burn+T);
    for ii in 1:T+burn
        if ii==1
            x[:,ii]     = P.Q*x_Shocks[:,ii];
            x_hat[:,ii] = (P.K*P.Y')*(x[:,ii])+P.K[:,:]*a_Shocks[:,ii];
        else
            x[:,ii]     = P.A*x[:,ii-1]+P.Q*x_Shocks[:,ii];
            x_hat[:,ii] = P.A*x_hat[:,ii-1]+(P.K*P.Y')*(x[:,ii]-P.A*x_hat[:,ii-1])+P.K[:,:]*a_Shocks[:,ii];
        end
        a[:,ii]  .= P.H'*x_hat[:,ii];
    end
    x     = x[:,burn+1:T+burn];
    x_hat = x_hat[:,burn+1:T+burn];
    a     = a[:,burn+1:T+burn];
    return(Dripirfs(T,x,x_hat,a))
end
