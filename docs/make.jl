import Pkg
Pkg.add("Documenter")
using Documenter
using FundamentosDMC
push!(LOAD_PATH,"../src/")
makedocs(
    modules=[FundamentosDMC],
    sitename="FundamentosDMC.jl",
    pages = [
        "Início" => "index.md",
        "Sistema simulado" => "sistema.md",
        "Simulación microcanónica" => "simple.md",
        "Termostato isocinético" => "isokinetic.md",
        "Termostato de Berendsen" => "berendsen.md",
        "Termostato de Langevin" => "langevin.md",
        "Monte-Carlo" => "montecarlo.md",
        "Comparando MD y MC" => "md_vs_mc.md",
        "Ejecución en 3D" => "running3D.md",
        "Help entries" => "help.md",
    ]
)
deploydocs(
    repo = "github.com/m3g/FundamentosDMC.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#" ],
)
