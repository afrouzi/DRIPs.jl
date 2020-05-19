using DRIPs
using Test

@testset "Sims2011.jl" begin
    # Write your own tests here.
    println("Running tests ...")
	    β = 0.9;
		ω = 1.0;
		A = [0.95 0.0; 0.0 0.4];
		Q = [√0.0975 0.0; 0.0 √0.86];
		H = [1.0; 1.0];
	println("Solving Sims (2011) ...")
	    p   = Drip(ω,β,A,Q,H)
	    @test p.ss.err < 1e-4
    println("Solving Sims (2011) with fixed capacity ...")
	    cap = DRIPs.capacity(p)
	    capn= DRIPs.capacity(p; unit = "nat")
	    pf  = Drip(cap,β,A,Q,H;fcap=true)
	    @test pf.ss.err < 1e-4
	    @test pf.ss.Σ_p ≈ p.ss.Σ_p atol = 1e-3
    println("Solving transition dynamics for Sims (2011) ...")
	    S   = DRIPs.Signal([1 0; 0 1],[0 0; 0 0])
	    pt  = Trip(p,S;T=30)
	    @test pt.err < 1e-4;
    println("Checking IRFs ...")
		pirfs     = irfs(p,T = 15)
		psims1    = simulate(p,T = 500, burn = 100, seed = 0)
		psims2    = simulate(p,T = 500, burn = 100, N=100, seed = 0)
		@test psims1.x == psims2.x
		psims     = simulate(p,T = 500, burn = 100, N=1, seed = 0)
		ptss      = Trip(p,p.ss.Σ_1;T = 15)
		ptssirfs  = irfs(ptss,T = 15)
		Sp        = DRIPs.Signal([0 0; 0 0], [1 0; 0 1]);
		ptssirfsp = irfs(p,Sp, T = 30)
		ptssirfsp = irfs(p,Sp, T = 10)
		ptssirfsp = irfs(p,Sp, T = 15)
		ptssirfspn= irfs(p,Sp; reoptimize = false, T = 15)
		@test pirfs.x_hat ≈ ptssirfs.x_hat atol = 1e-3
		@test pirfs.x_hat ≈ ptssirfsp.x_hat atol = 1e-3
		@test pirfs.x_hat ≈ ptssirfspn.x_hat atol = 1e-3
	println("Final tests on internal functions ...")
		x(j) = 0.5^j
		sum  = DRIPs.infinitesum(x)
		@test sum ≈ 2 atol = 1e-4
end
