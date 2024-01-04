cd(@__DIR__)
ENV["JULIA_DEBUG"] = "Documenter"
using ParticleTracking
using Documenter
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
import Literate
using GLMakie # to avoid precompilation within scripts

# compile tutorial to markdown
Literate.markdown(
    joinpath(@__DIR__, "src", "tutorial.jl"),
    joinpath(@__DIR__, "src");
)

pages = [
    "Home" => "index.md",
    "Tutorial" => "tutorial.md",
    #"API" => "api.md"
]

makedocs(
    sitename = "ParticleTracking.jl",
    authors = "Riccardo Foffi",
    modules = [ParticleTracking],
    pages = pages,
    expandfirst = ["index.md"],
    format = Documenter.HTML(
        prettyurls = CI,
    ),
    warnonly = true,
)

if CI
    deploydocs(;
        repo = "github.com/mastrof/ParticleTracking.jl",
        target = "build",
        push_preview = true
    )
end
