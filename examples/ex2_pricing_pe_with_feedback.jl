using DRIPs, Plots, LaTeXStrings, BenchmarkTools, Suppressor; pyplot();

ρ   = 0.6;        #persistence of money growth
σ_u = 0.1;        #std. deviation of shocks to money growth
α   = 0.8;        #degree of strategic complementarity
L   = 40;         #length of truncation
Hq  = ρ.^(0:L-1); #state-space rep. of Δq 

## specifying the primitives of the drip
ω   = 0.2;       
β   = 0.99;  
A   = [1 zeros(1,L-2) 0; Matrix(I,L-1,L-1) zeros(L-1,1)];
Q   = [σ_u; zeros(L-1,1)];

function ge_drip(ω,β,A,Q,          #primitives of drip except for H because H is endogenous
                 α,                #strategic complementarity
                 Hq,               #state space rep. of Δq
                 L;                #length of truncation
                 H0       = Hq,    #optional: initial guess for H (Hq is the true solution when α=0)
                 maxit    = 200,   #optional: max number of iterations for GE code 
                 tol      = 1e-4)  #optional: tolerance for iterations
    # set primitives
    err   = 1;  
    iter  = 0;
    M     = [zeros(1,L-1) 0; Matrix(I,L-1,L-1) zeros(L-1,1)];
    # iterate on GE
    while (err > tol) & (iter < maxit)
            if iter == 0 
                global ge  = solve_drip(ω,β,A,Q,H0, w = 0.9); 
            else 
                global ge  = solve_drip(ω,β,A,Q,H0;Ω0 = ge.Ω ,Σ0 = ge.Σ_1,maxit=15);
            end

            XFUN(jj) = ((I-ge.K*ge.Y')*ge.A)^jj * (ge.K*ge.Y') * (M')^jj
            X = infinitesum(XFUN; maxit=L, start = 0);  #E[x⃗]=X×x⃗

            XpFUN(jj) = α^jj * X^(jj)
            Xp = infinitesum(XpFUN; maxit=L, start = 0);

            H1 = (1-α)*Xp'*Hq;
            err= 0.5*norm(H1-H0,2)/norm(H0)+0.5*err;
            H0 = H1;
        
            iter += 1;
            if iter == maxit
                print("GE loop hit maxit\n")
            end
            println("Iteration $iter. Difference: $err")
    end
    return(ge)
end;

@suppress @benchmark ge_drip(ω,β,A,Q,α,Hq,L) setup = (ω=rand()) # solves and times the fixedp point for different values of ω 

@time @suppress  ge = ge_drip(ω,β,A,Q,α,Hq,L); # remove suppress to see convergence log

T      = L;
geirfs = dripirfs(ge,T);

M  = [zeros(1,T-1) 0; Matrix(I,T-1,T-1) zeros(T-1,1)];
dq = diagm(Hq)*geirfs.x[1,1,:];
Pi = (I-M)*geirfs.a[1,1,:];
y  = inv(I-M)*(dq-Pi);

p1 = plot(1:T,[dq,Pi],
     label = [L"Agg. Demand Growth ($\Delta q$)" L"Inflation ($\pi$)"]);

p2 = plot(1:T,y,
     label  = L"Output ($y$)");

plot(p1,p2,
    layout     = (1,2),
    xlim       = (1,20),
    lw         = 3,
    legend     = :topright,
    legendfont = font(12),
    tickfont   = font(12),
    size       = (900,370),
    framestyle = :box)
