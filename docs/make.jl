using Documenter, DRIPs

makedocs(
    modules  = [DRIPs],
    format = Documenter.HTML(prettyurls = "--local" ∉ ARGS),
    sitename = "DRIPs.jl",
    pages    = [
         "Home" => "index.md",
         "Overview" => "overview.md",
         "Syntax" => "syntax.md",
         "Examples and Replications" => Any["examples/sims2011/sims2011.md"]
         ]
    )
#makedocs()
deploydocs(
    target = "build",
    repo   = "github.com/afrouzi/DRIPs.jl.git",
    branch = "gh-pages",
    deps   = nothing,
    make   = nothing,
    devbranch = "master",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#"],
    push_preview = false
)
