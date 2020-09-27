using SafeTestsets

@safetestset "Sims 2010" begin include("Sims2010.jl") end
@safetestset "Mackowiack and Wiederholt 2009" begin include("MW2009.jl") end
@safetestset "Endogenous Feedback" begin include("EndFeedback.jl") end
