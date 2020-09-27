using SafeTestsets

@safetestset "Sims 2011" begin include("sims2010.jl") end
@safetestset "Mackowiack and Wiederholt 2009" begin include("MW2009.jl") end
@safetestset "Endogenous Feedback" begin include("EndFeedback.jl") end
