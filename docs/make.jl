using Documenter
using FundamentosDMC
makedocs(
    modules=[FundamentosDMC],
    sitename="FundamentosDMC.jl",
    top_menu = [
        "Português" => [
            "Inicio" => "pt/index.md",
            "Sistema simulado" => "pt/sistema.md",
            "Simulação microcanônica" => "pt/simple.md",
            "Termostato isocinético" => "pt/isokinetic.md",
            "Termostato de Berendsen" => "pt/berendsen.md",
            "Termostato de Langevin" => "pt/langevin.md",
            "Monte-Carlo" => "pt/montecarlo.md",
            "Comparando MD y MC" => "pt/md_vs_mc.md",
            "Execução em 3D" => "pt/running3D.md",
            "Distribução de velocidades" => "pt/velocities.md",
        ],
        "Español" => [
            "Início" => "es/index.md",
            "Sistema simulado" => "es/sistema.md",
            "Simulación microcanónica" => "es/simple.md",
            "Termostato isocinético" => "es/isokinetic.md",
            "Termostato de Berendsen" => "es/berendsen.md",
            "Termostato de Langevin" => "es/langevin.md",
            "Monte-Carlo" => "es/montecarlo.md",
            "Comparando MD y MC" => "es/md_vs_mc.md",
            "Ejecución en 3D" => "es/running3D.md",
            "Distribución de velocidades" => "es/velocities.md",
        ],
        "English" => [
            "Start" => "en/index.md",
            "Simulated system" => "en/sistema.md",
            "Microcanonical simulation" => "en/simple.md",
            "Isokinetic thermostat" => "en/isokinetic.md",
            "Berendsen thermostat" => "en/berendsen.md",
            "Langevin thermostat" => "en/langevin.md",
            "Monte Carlo" => "en/montecarlo.md",
            "Comparing MD and MC" => "en/md_vs_mc.md",
            "Running in 3D" => "en/running3D.md",
            "Velocity distribution" => "en/velocities.md",
        ],
        "Help entries" => ["help.md",],
    ]
)
deploydocs(
    repo = "github.com/m3g/FundamentosDMC.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#" ],
)
