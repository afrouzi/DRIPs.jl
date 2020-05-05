# compile markdowns
Literate.markdown("examples/src/ex1_pricing_pe_nofeedback.jl", "docs/src/examples/ex1_pricing_nofeedback"; documenter=true, execute = true, codefence = "```julia" => "```")
Literate.markdown("examples/src/ex2_pricing_pe_with_feedback.jl", "docs/src/examples/ex2_pricing_wfeedback"; documenter=true, execute = true, codefence = "```julia" => "```")
Literate.markdown("examples/src/ex3_Mackowiak_Wiederholt_2009.jl", "docs/src/examples/ex3_mw2009"; documenter=true, execute = true, codefence = "```julia" => "```")
Literate.markdown("examples/src/ex4_Sims_2011.jl", "docs/src/examples/ex4_sims2011"; documenter=true, execute = true, codefence = "```julia" => "```")

function update_notebooks(content)
    content = replace(content, "# [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/afrouzi/DRIPs.jl/binder?filepath=examples) to run and modify the following code (no software is needed on the local machine).

" => "")
    return content
end

# compile jupyter notebooks
Literate.notebook("examples/src/ex1_pricing_pe_nofeedback.jl", "examples/notebooks/"; documenter = true, preprocess = update_notebooks)
Literate.notebook("examples/src/ex2_pricing_pe_with_feedback.jl", "examples/notebooks/"; documenter = true, preprocess = update_notebooks)
Literate.notebook("examples/src/ex4_Sims_2011.jl", "examples/notebooks/"; documenter = true, preprocess = update_notebooks)
Literate.notebook("examples/src/ex3_Mackowiak_Wiederholt_2009.jl","examples/notebooks/"; documenter = true, preprocess = update_notebooks)
