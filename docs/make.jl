import Pkg
Pkg.add("Documenter")
using Documenter
using CELFI
push!(LOAD_PATH,"../src/")
makedocs(
    modules=[CELFI],
    sitename="CELFI.jl",
    pages = [
        "Home" => "index.md",
    ]
)
deploydocs(
    repo = "github.com/m3g/CELFI.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#" ],
)
