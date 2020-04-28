using DRIPs
using Test

@testset "DRIP.jl" begin
    # Write your own tests here.
    println("Running tests ...")
	    β = 0.9;
		ω = 1.0;
		A = [0.95 0.0; 0.0 0.4];
		Q = [√0.0975 0.0; 0.0 √0.86];
		H = [1.0; 1.0];
	println("Solving Sims (2011) ...")
		p   = Drip();
	    p   = Drip(ω,β,A,Q,H);
	    p   = solve_drip(p);
	    pp  = solve_drip(ω,β,A,Q,H)
	    @test p.err < 1e-4
	    @test p.Σ_p ≈ pp.Σ_p atol = 1e-4
    println("Solving Sims (2011) with fixed capacity ...")
	    cap = DRIPs.capacity(p);
	    capn= DRIPs.capacity(p; unit = "nat");
	    pf  = p;
	    pf  = Drip(cap,β,A,Q,H);
	    pf  = solve_drip(pf;fcap=true);
	    @test pf.err < 1e-4
	    @test pf.Σ_p ≈ p.Σ_p atol = 1e-2
    println("Solving transition dynamics for Sims (2011) ...")
	    S   = Signal([1 0; 0 1],[0 0;0 0]);
	    pt  = solve_trip(p,S;T=30);
	    @test pt.err < 1e-4;
    println("Checking IRFs ...")
		pirfs     = dripirfs(p,15);
		ptss      = solve_trip(p,p.Σ_1;T = 15);
		ptssirfs  = dripirfs(ptss,15);
		Sp        = Signal([0 0; 0 0], [1 0; 0 1]);
		ptssirfsp = dripirfs(p,15,Sp);
		ptssirfspn= dripirfs(p,15,Sp; reoptimize = false);
		@test pirfs.x_hat ≈ ptssirfs.x_hat atol = 1e-4
		@test pirfs.x_hat ≈ ptssirfsp.x_hat atol = 1e-4
		@test pirfs.x_hat ≈ ptssirfspn.x_hat atol = 1e-4
	println("Final tests on internal functions ...")
		x(j) = 0.5^j
		sum  = DRIPs.infinitesum(x)
		@test sum ≈ 2 atol = 1e-4
end
