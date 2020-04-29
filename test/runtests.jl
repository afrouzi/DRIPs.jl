using SafeTestsets

@safetestset "Sims 2011" begin include("Sims2011.jl") end
@safetestset "Mackowiack and Wiederholt 2009" begin include("MW2009.jl") end
@safetestset "Endogenous Feedback" begin include("EndFeedback.jl") end 
