using Documenter, DRIPs

# uncomment the following line to re-compile examples
# using Literate; include("compile_examples.jl");

makedocs(
    modules  = [DRIPs],
    format = Documenter.HTML(prettyurls = "--local" ∉ ARGS),
    sitename = "DRIPs.jl",
    pages    = [
         "Home" => "index.md",
         "Overview" => "overview.md",
         "Syntax" => "syntax.md",
         "Examples and Replications" => Any[
            "examples/ex1_pricing_nofeedback/ex1_pricing_pe_nofeedback.md",
            "examples/ex2_pricing_wfeedback/ex2_pricing_pe_with_feedback.md",
            "examples/ex3_mw2009/ex3_Mackowiak_Wiederholt_2009.md",
            "examples/ex4_sims2010/ex4_Sims_2011.md",
            "examples/ex5_ay2020/ex5_Afrouzi_Yang_2019.md",
            ]
         ]
    )


const devurl = "dev"

deploydocs(
    target = "build",
    repo   = "github.com/afrouzi/DRIPs.jl.git",
    branch = "gh-pages",
    deps   = nothing,
    make   = nothing,
    devbranch = "master",
    devurl = devurl,
    versions = ["v#.#", devurl => devurl],
    push_preview = false
)
