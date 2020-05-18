using DRIPs
using Test
using LinearAlgebra

@testset "MW2009.jl" begin
	## Parameters
	ρ  = 0.95;
	σq = 0.01;
	σz = 11.8*σq;
	κ  = 3;
	ξ  = 1;
	α  = 1 - 0.15;

	## Primitives of drip
	L  = 21; # length of trunction
	A  = [zeros(1,L);[Matrix(I,L-1,L-1);zeros(1,L-1)]']; # MW truncate the state space with linear irfs of length 20
	Qq = zeros(L,1); Qq[1]=σq;
	Qz = zeros(L,1); Qz[1]=σz;
	H  = zeros(L,1); H[1:21] = Array(1:-1/20:0);

	## Function for solving the fixed point of the aggregate shock given ω
	# Note that for a certain range of ω the fixed point is not unique:
	# there is a fixed point with zero capacity and one with positive capacity.
	# In this range, we would like to find the fixed point in which the firms are
	# acqurining a positive amount of information.

	function agg_drip(ω,A,Qq,        #primitives of drip except for H because H is endogenous
	                  α,             #strategic complementarity
	                  H;             #state space rep. of q
	                  β     = 1,     #optional: discount factor, MW's parameterization implies β = 1
	                  H0    = H,     #optional: initial guess for HΔ (H is the true solution when α=0)
	                  maxit = 10000, #optional: max number of iterations for GE code
	                  tol   = 1e-4,  #optional: tolerance for iterations
	                  w     = 1)     #optional: update weight for RI
	    # set primitives
	    errmin= 1;
	    err   = 1;
	    iter  = 0;
	    L     = length(H);
	    while (err > tol) & (iter < maxit)
	            if iter == 0
	                global agg  = Drip(ω,β,A,Qq,H0;w = w);
	            else
	                global agg  = Drip(ω,β,A,Qq,H0;Ω0 = agg.Ω , Σ0 = agg.Σ_1,w = w);
	            end

	            XFUN(jj) = ((I-agg.K*agg.Y')*agg.A)^jj * (agg.K*agg.Y') * (agg.A')^jj
	            X = DRIPs.infinitesum(XFUN; maxit=200, start = 0);  #E[x⃗]=X×x⃗

	            XpFUN(jj) = α^jj * X^(jj)
	            Xp = DRIPs.infinitesum(XpFUN; maxit=200, start = 0);

	            H1 = (1-α)*Xp'*H;

	            err= 0.5*norm(H1-H0,2)/norm(H0)+0.5*err;
	            if DRIPs.capacity(agg) < 1e-2 # perturb the initial guess if solution is the zero capacity one
	                H0 = H0+rand(L).*(H-H0);
	            else # store the solution if it has positive capacity
	                H0 = H1;
	                if err < errmin
	                    global aggmin = agg;
	                    errmin = err;
	                end
	            end
	            iter += 1;
	    end
	    return(aggmin,min(err,errmin))
	end;
	println("Solving Mackowiack and Wiederholt (2009) ...")
	## Solve for κ = 3
	agg, err = agg_drip(2.4σq^2,A,Qq,α,H; H0 = rand(L), maxit = 500, w = 0.95)
	@test agg.err < 1e-4
	agg, err = agg_drip(2.82σq^2,A,Qq,α,H; H0 = rand(L), maxit = 500, w = 0.95)
	@test agg.err < 1e-4
	agg, err = agg_drip(0.01*σq^2,A,Qq,α,H; H0 = rand(L), maxit = 500, w = 0.95)
	@test agg.err < 1e-4
	agg, err = agg_drip(10*σq^2,A,Qq,α,H; H0 = rand(L), maxit = 500, w = 0.95)
	@test agg.err < 1e-4
	idi  = Drip(2.4σq^2,1,A,Qz,H,w = 0.9);
	@test idi.err < 1e-4

end
