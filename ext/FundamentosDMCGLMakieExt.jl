module FundamentosDMCGLMakieExt

using GLMakie
using FundamentosDMC
using FundamentosDMC: Point2D, System, Options
import FundamentosDMC: potential, forces!, kinetic, remove_drift!, image, dim

# Simulation kinds exposed in the interface, and their display labels.
const KIND_OPTIONS = [
    "Microcanonical (NVE)" => :md,
    "Isokinetic bath (NVT)" => :md_isokinetic,
    "Berendsen bath (NVT)" => :md_berendsen,
    "Langevin bath (NVT)" => :md_langevin,
    "Monte Carlo" => :mc,
]
const KIND_LABELS = Dict(v => k for (k, v) in KIND_OPTIONS)

const VELOCITY_OPTIONS = ["Normal" => :normal, "Flat" => :flat, "Zero" => :zero]
const VELOCITY_LABELS = Dict(v => k for (k, v) in VELOCITY_OPTIONS)

#
# Mutable state of the interactive simulation: the current parameters (which
# mirror `System` and `Options`), the current physical state (positions,
# velocities, forces), and the history of the logged quantities used for the
# energy/temperature plots. A fresh `SimState` is built every time a
# parameter is changed, which is how the interface implements "restart on
# option change".
#
mutable struct SimState
    # Parameters (mirror System/Options)
    kind::Symbol
    n::Int
    Lx::Float64
    Ly::Float64
    dt::Float64
    nsteps::Int
    eps::Float64
    sig::Float64
    initial_velocities::Symbol
    kT::Float64
    ibath::Int
    iequil::Int
    tau::Int
    lambda::Float64
    alpha::Float64
    minimize_first::Bool

    # Derived System/Options
    sys::System{Point2D}
    opt::Options

    # Physical state
    x::Vector{Point2D}
    v::Vector{Point2D}
    f::Vector{Point2D}
    flast::Vector{Point2D}
    xtrial::Vector{Point2D}
    ucurrent::Float64
    naccepted::Int

    # Bookkeeping and logs
    step::Int
    time::Float64
    steps_history::Vector{Int}
    potential_history::Vector{Float64}
    kinetic_history::Vector{Float64}
    total_history::Vector{Float64}
    temperature_history::Vector{Float64}

    # Run control
    running::Bool
    stop::Bool
end

function SimState(;
    kind::Symbol=:md,
    n::Int=100,
    Lx::Float64=100.0,
    Ly::Float64=100.0,
    dt::Float64=0.05,
    nsteps::Int=2000,
    eps::Float64=1.0,
    sig::Float64=2.0,
    initial_velocities::Symbol=:normal,
    kT::Float64=0.6,
    ibath::Int=10,
    iequil::Int=10,
    tau::Int=10,
    lambda::Float64=0.1,
    alpha::Float64=0.1,
    minimize_first::Bool=true,
)
    n = max(n, 2)
    nsteps = max(nsteps, 1)
    Lx = max(Lx, 1.0)
    Ly = max(Ly, 1.0)

    sys = System(n=n, sides=[Lx, Ly])
    opt = Options(; dt, nsteps, eps, sig, initial_velocities, kT, ibath, iequil, tau, lambda, alpha)

    if minimize_first
        minimize!(sys, opt)
    end

    x = copy(sys.x0)
    v = kind == :mc ? zeros(Point2D, n) : init_velocities(sys, opt)
    f = zeros(Point2D, n)
    flast = zeros(Point2D, n)
    xtrial = copy(x)

    u0 = potential(x, sys, opt)
    k0 = kinetic(v)
    if kind != :mc
        forces!(f, x, sys, opt)
        flast .= f
    end

    return SimState(
        kind, n, Lx, Ly, dt, nsteps, eps, sig, initial_velocities, kT, ibath, iequil, tau, lambda, alpha, minimize_first,
        sys, opt,
        x, v, f, flast, xtrial, u0, 0,
        0, 0.0,
        [0], [u0], [k0], [u0 + k0], [k0 / n],
        false, false,
    )
end

#
# Performs a single Velocity-Verlet integration step, reproducing exactly
# the algorithms in md-simple.jl, md-isokinetic.jl, md-berendsen.jl and
# md-langevin.jl, selected by `state.kind`. Returns `false` if the
# simulation exploded (as in the original `md` functions).
#
function md_step!(state::SimState)
    state.step += 1
    istep = state.step
    x, v, f, flast = state.x, state.v, state.f, state.flast
    sys, opt = state.sys, state.opt
    dt = opt.dt
    T = Point2D

    @. x = x + v * dt + 0.5 * f * dt^2
    @. flast = f
    forces!(f, x, sys, opt)

    if state.kind == :md_langevin
        @. f = f - opt.lambda * v
        @. v = v + 0.5 * (f + flast) * dt + sqrt(2 * opt.lambda * opt.kT * opt.dt) * randn(T)
    else
        @. v = v + 0.5 * (f + flast) * dt
    end
    remove_drift!(v)

    ustep = potential(x, sys, opt)
    kstep = kinetic(v)
    kavg = kstep / sys.n

    if ustep > 1e10
        return false
    end

    if state.kind == :md_isokinetic && istep <= opt.iequil && mod(istep, opt.ibath) == 0
        @. v = v * sqrt((dim(T) * opt.kT / 2) / kavg)
    elseif state.kind == :md_berendsen && istep <= opt.iequil
        lam = sqrt(1 + (opt.dt / opt.tau) * ((dim(T) * opt.kT / 2) / kavg - 1))
        @. v = v * lam
    end

    state.time += dt
    push!(state.steps_history, istep)
    push!(state.potential_history, ustep)
    push!(state.kinetic_history, kstep)
    push!(state.total_history, ustep + kstep)
    push!(state.temperature_history, kavg)
    return true
end

# Performs a single Monte Carlo trial move, reproducing mc.jl.
function mc_step!(state::SimState)
    state.step += 1
    x, xtrial = state.x, state.xtrial
    sys, opt = state.sys, state.opt
    T = Point2D

    @. xtrial = x + opt.alpha * randn(T)
    utrial = potential(xtrial, sys, opt)

    if (utrial < state.ucurrent) || (exp(-(utrial - state.ucurrent) / opt.kT) > rand())
        state.ucurrent = utrial
        @. x = xtrial
        state.naccepted += 1
    end

    push!(state.steps_history, state.step)
    push!(state.potential_history, state.ucurrent)
    push!(state.kinetic_history, 0.0)
    push!(state.total_history, state.ucurrent)
    push!(state.temperature_history, 0.0)
    return true
end

# Runs the simulation to completion (or until stopped), notifying `obs`
# after every step so the interface redraws. Meant to be `@async`ed.
#
# Every write to `obs` is guarded by an identity check against `state`: if a
# parameter was changed while this loop was running, `restart!` will have
# already replaced `obs[]` with a brand new `SimState`, and this (now stale)
# loop must not clobber it back.
function run!(obs::Observable{SimState})
    state = obs[]
    state.running = true
    state.stop = false
    obs[] = state

    step! = state.kind == :mc ? mc_step! : md_step!
    exploded = false
    while state.step < state.nsteps && !state.stop
        exploded = !step!(state)
        obs[] === state && (obs[] = state)
        exploded && break
        sleep(1 / 60)
    end
    exploded && @warn "Simulation exploded (potential energy too large). Stopping."

    state.running = false
    obs[] === state && (obs[] = state)
    return nothing
end

# Rebuilds the simulation state with a single field changed, which
# restarts the simulation with the new parameters (or a fresh random
# configuration, if `field` is set to its current value).
function restart!(obs::Observable{SimState}, field::Symbol, value)
    old = obs[]
    old.stop = true
    kwargs = Dict{Symbol,Any}(
        :kind => old.kind, :n => old.n, :Lx => old.Lx, :Ly => old.Ly,
        :dt => old.dt, :nsteps => old.nsteps, :eps => old.eps, :sig => old.sig,
        :initial_velocities => old.initial_velocities, :kT => old.kT,
        :ibath => old.ibath, :iequil => old.iequil, :tau => old.tau,
        :lambda => old.lambda, :alpha => old.alpha, :minimize_first => old.minimize_first,
    )
    kwargs[field] = value
    obs[] = SimState(; kwargs...)
    return nothing
end

function particles_title(s::SimState)
    label = KIND_LABELS[s.kind]
    if s.kind == :mc
        ar = s.step == 0 ? 0.0 : 100 * s.naccepted / s.step
        return "$label  |  step $(s.step)/$(s.nsteps)  |  acceptance = $(round(ar, digits=1))%"
    else
        return "$label  |  step $(s.step)/$(s.nsteps)  |  t = $(round(s.time, digits=2))"
    end
end

function FundamentosDMC.simulate_gui(; n::Int=100, sides=(100.0, 100.0), kind::Symbol=:md)
    state = SimState(; kind, n, Lx=Float64(sides[1]), Ly=Float64(sides[2]))
    obs = Observable(state)

    GLMakie.activate!(title="FundamentosDMC - Interactive simulation")
    fig = Figure(size=(1500, 820))

    controls = fig[1, 1] = GridLayout(tellwidth=false, valign=:top)

    row = 0
    next_row!() = (row += 1; row)

    #
    # Simulation type and initial velocities menus
    #
    r = next_row!()
    Label(controls[r, 1], "Simulation type", halign=:right)
    kind_menu = Menu(controls[r, 2], options=first.(KIND_OPTIONS), default=KIND_LABELS[state.kind])
    on(kind_menu.selection) do s
        restart!(obs, :kind, Dict(KIND_OPTIONS)[s])
    end

    r = next_row!()
    Label(controls[r, 1], "Initial velocities", halign=:right)
    iv_menu = Menu(controls[r, 2], options=first.(VELOCITY_OPTIONS), default=VELOCITY_LABELS[state.initial_velocities])
    on(iv_menu.selection) do s
        restart!(obs, :initial_velocities, Dict(VELOCITY_OPTIONS)[s])
    end

    #
    # Numeric parameters
    #
    function add_numeric_row!(label, field::Symbol, valtype::Type)
        r = next_row!()
        Label(controls[r, 1], label, halign=:right)
        tb = Textbox(
            controls[r, 2];
            placeholder=@lift(string(getfield($obs, field))),
            validator=valtype,
            reset_on_defocus=true,
        )
        on(tb.stored_string) do s
            restart!(obs, field, parse(valtype, s))
        end
        return tb
    end

    add_numeric_row!("Number of particles (n)", :n, Int)
    add_numeric_row!("Box side Lx", :Lx, Float64)
    add_numeric_row!("Box side Ly", :Ly, Float64)
    add_numeric_row!("Time step (dt)", :dt, Float64)
    add_numeric_row!("Number of steps (nsteps)", :nsteps, Int)
    add_numeric_row!("LJ eps", :eps, Float64)
    add_numeric_row!("LJ sig", :sig, Float64)
    add_numeric_row!("Target temperature (kT)", :kT, Float64)
    add_numeric_row!("Isokinetic bath frequency (ibath)", :ibath, Int)
    add_numeric_row!("Equilibration steps (iequil)", :iequil, Int)
    add_numeric_row!("Berendsen relaxation time (tau)", :tau, Int)
    add_numeric_row!("Langevin friction (lambda)", :lambda, Float64)
    add_numeric_row!("MC trial displacement (alpha)", :alpha, Float64)

    r = next_row!()
    Label(controls[r, 1], "Minimize before run", halign=:right)
    minimize_cb = Checkbox(controls[r, 2], checked=state.minimize_first)
    on(minimize_cb.checked) do checked
        restart!(obs, :minimize_first, checked)
    end

    #
    # Run / Stop / Reset buttons
    #
    r = next_row!()
    buttons = controls[r, 1:3] = [Button(fig, label="Run"), Button(fig, label="Stop"), Button(fig, label="Reset")]
    on(buttons[1].clicks) do _
        obs[].running || @async run!(obs)
    end
    on(buttons[2].clicks) do _
        obs[].stop = true
    end
    on(buttons[3].clicks) do _
        restart!(obs, :n, obs[].n)
    end

    #
    # Visualization: particle motion, energies, and temperature
    #
    viz = fig[1, 2] = GridLayout()

    ax_particles = Axis(
        viz[1:2, 1],
        aspect=DataAspect(),
        title=@lift(particles_title($obs)),
    )
    boxpts = @lift(let sides = $(obs).sys.sides
        hx, hy = sides[1] / 2, sides[2] / 2
        Point2f[(-hx, -hy), (hx, -hy), (hx, hy), (-hx, hy), (-hx, -hy)]
    end)
    lines!(ax_particles, boxpts, color=:black)
    positions = @lift([Point2f(image(p, $(obs).sys.sides)) for p in $(obs).x])
    scatter!(
        ax_particles, positions;
        markersize=@lift(Float32($(obs).sig)),
        markerspace=:data,
        color=(:dodgerblue, 0.85),
        strokewidth=1,
        strokecolor=:black,
    )

    last_box = Ref((NaN, NaN))
    on(obs) do s
        b = (s.sys.sides[1], s.sys.sides[2])
        if b != last_box[]
            last_box[] = b
            hx, hy = b[1] / 2, b[2] / 2
            pad = 0.05 * max(b[1], b[2])
            limits!(ax_particles, -hx - pad, hx + pad, -hy - pad, hy + pad)
        end
    end

    ax_energy = Axis(viz[1, 2], xlabel="step", ylabel="Energy", title="Potential, kinetic and total energy")
    lines!(ax_energy, @lift(Point2f.($(obs).steps_history, $(obs).potential_history)), color=:royalblue, label="Potential")
    lines!(ax_energy, @lift(Point2f.($(obs).steps_history, $(obs).kinetic_history)),
        color=:seagreen, label="Kinetic", visible=@lift($(obs).kind != :mc))
    lines!(ax_energy, @lift(Point2f.($(obs).steps_history, $(obs).total_history)),
        color=:firebrick, label="Total", visible=@lift($(obs).kind != :mc))
    axislegend(ax_energy, position=:rb)

    ax_temp = Axis(viz[2, 2], xlabel="step", ylabel="Temperature", title="Temperature (average kinetic energy/particle)")
    lines!(ax_temp, @lift(Point2f.($(obs).steps_history, $(obs).temperature_history)),
        color=:purple, visible=@lift($(obs).kind != :mc))
    hlines!(ax_temp, @lift(Float32($(obs).kT)), color=:gray, linestyle=:dash)

    on(obs) do _
        autolimits!(ax_energy)
        autolimits!(ax_temp)
    end

    colsize!(viz, 1, Relative(0.45))
    colsize!(viz, 2, Relative(0.55))
    colsize!(fig.layout, 1, Fixed(300))
    colsize!(fig.layout, 2, Auto())

    return fig
end

end # module
