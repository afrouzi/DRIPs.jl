# compile markdowns
Literate.markdown("examples/ex1_pricing_pe_nofeedback.jl", "docs/src/examples/ex1_pricing_nofeedback"; documenter=true, execute = true, codefence = "```julia" => "```")
Literate.markdown("examples/ex2_pricing_pe_with_feedback.jl", "docs/src/examples/ex2_pricing_wfeedback"; documenter=true, execute = true, codefence = "```julia" => "```")
Literate.markdown("examples/ex3_Mackowiak_Wiederholt_2009.jl", "docs/src/examples/ex3_mw2009"; documenter=true, execute = true, codefence = "```julia" => "```")
Literate.markdown("examples/ex4_Sims_2011.jl", "docs/src/examples/ex4_sims2011"; documenter=true, execute = true, codefence = "```julia" => "```")

# compile jupyter notebooks
Literate.notebook("examples/ex1_pricing_pe_nofeedback.jl", "examples/"; documenter = true)
Literate.notebook("examples/ex2_pricing_pe_with_feedback.jl", "examples/"; documenter = true)
Literate.notebook("examples/ex4_Sims_2011.jl", "examples/"; documenter = true)
Literate.notebook("examples/ex3_Mackowiak_Wiederholt_2009.jl","examples/"; documenter = true)