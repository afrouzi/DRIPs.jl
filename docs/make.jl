using Documenter, DRIPs

# makedocs(
#     modules  = [DRIPs],
#     format = Documenter.HTML(prettyurls = "--local" ∉ ARGS),
#     sitename = "DRIPs.jl",
#     pages    = Any[
#         "Home" => "index.md",
#        ]
#     )
makedocs()

deploydocs(
    repo   = "github.com/afrouzi/DRIPs.jl.git",
)
