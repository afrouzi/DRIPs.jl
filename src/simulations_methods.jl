"""
         simulate(p    :: Drip;
                  T    :: Int64 = 500,                    # Optional: length of simulation
                  burn :: Int64 = 100,                    # Optional: length of initial burn (in addition to T)
                  N    :: Union{Int4,Nothing} = nothing,  # Optional: number of simulated agents (returns the average beliefs of a large set of agents by default)
                  seed :: Union{Int64,Nothing} = nothing  # Optional: seed number for fundamental shocks
                  ) -> Path
Returns a `Path` structure with a simulated path of the fundamental (`x`), beliefs (`x_hat`)
    and actions (`a`)
    **under the steady state information structure**. In particular, if `n` is the
    dimension of `x` and `m` is the dimension of `a`, then

* `x` has dimension `n*T` where `x(:,t)` is the simulated value of `x` at time `t`.
* `x_hat` has dimension `n*T*N` where `x_hat(:,i,t)` is the simulated value of `x_hat` of agent `i` at time `t`.
* `a` has dimension `m*T*N` where `a(:,i,t)` is the the simulated value of `x_hat` of agent `i` at time `t`.
"""
function simulate(p    :: Drip;
                  T    :: Int64 = 500,                    # Optional: length of simulation
                  burn :: Int64 = 100,                    # Optional: length of initial burn (in addition to T)
                  N    :: Union{Int64,Nothing} = nothing, # Optional: number of simulated agents (returns the average beliefs of a large set of agents by default)
                  seed :: Union{Int64,Nothing} = nothing  # Optional: seed number for fundamental shocks
                  )
    # create flag for if N is nothing and change N to 1
    average_belief = (N == nothing)
    if average_belief; N = 1; end
    # get dimensions
    (n,m)    = size(p.H)
    (_,k)    = size(p.Q)
    # set seed for fundamental shocks
    seed != nothing ? Random.seed!(seed) : Random.seed!()
    x_Shocks = randn(k,burn+T);

    Random.seed!()
    a_Shocks = average_belief ? zeros(m,burn+T,N) : randn(m,burn+T,N);
    x        = zeros(n,burn+T);
    x_hat    = zeros(n,burn+T,N);
    a        = zeros(m,burn+T,N);

    SqΣz = average_belief ? 0 : sqrt(p.ss.Σ_z);

    for t in 2:T+burn
        x[:,t] = p.A*x[:,t-1]+p.Q*x_Shocks[:,t];
        for i in 1:N
            x_hat[:,t,i] = p.A*x_hat[:,t-1,i]+(p.ss.K*p.ss.Y')*(x[:,t]-p.A*x_hat[:,t-1,i])+p.ss.K*SqΣz*a_Shocks[:,t,i];
            a[:,t,i]    .= p.H'*x_hat[:,t,i];
        end
    end
    x     = x[:,burn+1:T+burn];
    x_hat = average_belief ? x_hat[:,burn+1:T+burn] : x_hat[:,burn+1:T+burn,:];
    a = average_belief ? a[:,burn+1:T+burn] : a[:,burn+1:T+burn,:];
    return(Path(T,x,x_hat,a))
end
