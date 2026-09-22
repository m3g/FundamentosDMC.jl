module FundamentosDMCGLMakieExt

using GLMakie
using Printf: @sprintf
using FundamentosDMC
using FundamentosDMC: Point2D, System, Options
import FundamentosDMC: potential, forces!, kinetic, remove_drift!, image, dim

# Simulation kinds exposed in the interface, and their display labels.
const KIND_OPTIONS = [
    "NVE" => :md,
    "NVT - Isokinetic" => :md_isokinetic,
    "NVT - Berendsen" => :md_berendsen,
    "NVT - Langevin" => :md_langevin,
    "Monte Carlo" => :mc,
]
const KIND_LABELS = Dict(v => k for (k, v) in KIND_OPTIONS)

const VELOCITY_OPTIONS = ["Normal" => :normal, "Flat" => :flat, "Zero" => :zero]
const VELOCITY_LABELS = Dict(v => k for (k, v) in VELOCITY_OPTIONS)

# Compact control-panel styling, so labels/boxes don't overflow the figure.
const CTRL_FONTSIZE = 11
const CTRL_HEIGHT = 19
const CTRL_TEXTBOX_WIDTH = 140

#
# Mutable state of the interactive simulation: the current parameters (which
# mirror `System` and `Options`), the current physical state (positions,
# velocities, forces), and the history of the logged quantities used for the
# energy/temperature plots. A fresh `SimState` is built every time a
# parameter is changed, which is how the interface implements "restart on
# option change". Minimization is *not* performed automatically here: it is
# a separate, explicit action triggered from the interface (see `minimize_now!`).
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
    run_start_step::Int
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
    iequil::Int=1000,
    tau::Int=10,
    lambda::Float64=0.1,
    alpha::Float64=0.05,
    x0::Union{Nothing,Vector{Point2D}}=nothing,
    v0::Union{Nothing,Vector{Point2D}}=nothing,
)
    n = max(n, 2)
    nsteps = max(nsteps, 1)
    Lx = max(Lx, 1.0)
    Ly = max(Ly, 1.0)

    # Reuse the given coordinates/velocities when they are compatible with
    # the (possibly new) number of particles; otherwise fall back to a
    # fresh random configuration. This is how a parameter change continues
    # the current trajectory instead of jumping to a new configuration.
    reuse_x = x0 !== nothing && length(x0) == n
    sys = reuse_x ? System(n=n, x0=copy(x0), sides=[Lx, Ly]) : System(n=n, sides=[Lx, Ly])
    opt = Options(; dt, nsteps, eps, sig, initial_velocities, kT, ibath, iequil, tau, lambda, alpha)

    x = copy(sys.x0)
    v = if kind == :mc
        zeros(Point2D, n)
    elseif v0 !== nothing && length(v0) == n
        copy(v0)
    else
        init_velocities(sys, opt)
    end
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
        kind, n, Lx, Ly, dt, nsteps, eps, sig, initial_velocities, kT, ibath, iequil, tau, lambda, alpha,
        sys, opt,
        x, v, f, flast, xtrial, u0, 0,
        0, 0, 0.0,
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
# `opt.iequil`/`opt.ibath` gate the isokinetic/Berendsen rescaling relative
# to `state.run_start_step` (the step count at which the current "Run"
# click started, set in `run!`), not relative to absolute step 0 — so
# clicking "Run" again always re-applies equilibration to the *next*
# `iequil` steps of that run, rather than only ever to the first `iequil`
# steps of the whole trajectory.
function md_step!(state::SimState)
    state.step += 1
    istep = state.step
    rel_step = istep - state.run_start_step
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

    if state.kind == :md_isokinetic && rel_step <= opt.iequil && mod(rel_step, opt.ibath) == 0
        @. v = v * sqrt((dim(T) * opt.kT / 2) / kavg)
    elseif state.kind == :md_berendsen && rel_step <= opt.iequil
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

# Appends `msg` (possibly several lines) to the on-screen status log,
# keeping only the last 4 lines. Kept as a plain Observable{Vector{String}},
# separate from `SimState`, so it persists across restarts instead of being
# wiped out every time the simulation parameters change.
function log!(log_obs::Observable{Vector{String}}, msg::AbstractString)
    lines = [String(strip(l)) for l in split(msg, '\n') if !isempty(strip(l))]
    isempty(lines) && return nothing
    newlog = vcat(log_obs[], lines)
    log_obs[] = newlog[max(1, end - 3):end]
    return nothing
end

# Runs the simulation to completion (or until stopped), notifying `obs`
# after every step so the interface redraws. Meant to be `@async`ed.
#
# `nsteps` (and, through `run_start_step`, `iequil`/`ibath` in `md_step!`)
# count from the step the run *starts at*, not from absolute step 0:
# clicking "Run" always advances the simulation by `nsteps` more steps from
# wherever it currently is, and re-applies equilibration to the next
# `iequil` steps of this run.
#
# Every write to `obs` is guarded by an identity check against `state`: it is
# only ever swapped out from under a running loop by `reset!` (a genuinely
# fresh restart); `update_params!` mutates this very state object in place,
# so a running loop simply keeps going, picking up new parameters (and even
# a new `kind`, hence re-selecting `step!` on every iteration below) on the
# very next step.
function run!(obs::Observable{SimState}, log_obs::Observable{Vector{String}})
    state = obs[]
    state.run_start_step = state.step
    state.running = true
    state.stop = false
    obs[] = state
    log!(log_obs, "Running $(KIND_LABELS[state.kind]) from step $(state.step) to $(state.run_start_step + state.nsteps)...")

    exploded = false
    while state.step < state.run_start_step + state.nsteps && !state.stop
        step! = state.kind == :mc ? mc_step! : md_step!
        exploded = !step!(state)
        obs[] === state && (obs[] = state)
        exploded && break
        sleep(1 / 60)
    end
    if exploded
        log!(log_obs, "Simulation exploded (potential energy too large). Stopping.")
    else
        log!(log_obs, "Stopped at step $(state.step)/$(state.run_start_step + state.nsteps).")
    end

    state.running = false
    obs[] === state && (obs[] = state)
    return nothing
end

# Applies a single changed parameter (or a new `kind`) directly to the
# *current* SimState, in place: the step counter, elapsed time, and the
# whole plot history are left untouched, so tweaking a parameter (say,
# `kT`) or switching between NVE/NVT/MC continues the running trajectory
# and its plots instead of restarting them. A running simulation loop is
# not interrupted either (see `run!`), since it reads `state.sys`,
# `state.opt` and `state.kind` fresh on every step.
#
# Coordinates/arrays are only regenerated when `n` changes (the array
# length wouldn't match anymore). Velocities are only regenerated when
# `initial_velocities` itself is the field being changed (a new
# distribution was explicitly requested), or when switching away from
# Monte Carlo (which does not carry meaningful velocities). An actual
# fresh random configuration is only ever produced at startup or via
# `reset!` (the "Reset" button).
function update_params!(obs::Observable{SimState}, field::Symbol, value)
    s = obs[]
    old_kind = s.kind
    setfield!(s, field, value)

    s.opt = Options(;
        dt=s.dt, nsteps=s.nsteps, eps=s.eps, sig=s.sig,
        initial_velocities=s.initial_velocities, kT=s.kT,
        ibath=s.ibath, iequil=s.iequil, tau=s.tau, lambda=s.lambda, alpha=s.alpha,
    )

    if field == :n
        n = max(s.n, 2)
        s.n = n
        s.sys = System(n=n, sides=[s.Lx, s.Ly])
        s.x = copy(s.sys.x0)
        s.xtrial = copy(s.x)
        s.f = zeros(Point2D, n)
        s.flast = zeros(Point2D, n)
        s.v = s.kind == :mc ? zeros(Point2D, n) : init_velocities(s.sys, s.opt)
        s.naccepted = 0
    elseif field in (:Lx, :Ly)
        s.sys = System(n=s.n, x0=copy(s.x), sides=[s.Lx, s.Ly])
    end

    if field == :initial_velocities || (field == :kind && old_kind == :mc && s.kind != :mc)
        s.v = init_velocities(s.sys, s.opt)
    elseif field == :kind && s.kind == :mc
        s.v = zeros(Point2D, s.n)
    end

    if s.kind == :mc
        s.ucurrent = potential(s.x, s.sys, s.opt)
    else
        forces!(s.f, s.x, s.sys, s.opt)
        s.flast .= s.f
    end

    obs[] = s
    return nothing
end

# Generates a genuinely fresh random configuration, keeping the current
# parameters. This is the explicit user action (the "Reset" button) that
# `update_params!` deliberately does not perform on its own.
function reset!(obs::Observable{SimState})
    old = obs[]
    old.stop = true
    obs[] = SimState(;
        kind=old.kind, n=old.n, Lx=old.Lx, Ly=old.Ly,
        dt=old.dt, nsteps=old.nsteps, eps=old.eps, sig=old.sig,
        initial_velocities=old.initial_velocities, kT=old.kT,
        ibath=old.ibath, iequil=old.iequil, tau=old.tau,
        lambda=old.lambda, alpha=old.alpha,
    )
    return nothing
end

# Minimizes the *current* configuration in place (explicit user action, not
# run automatically). Pauses any running loop first, then recomputes the
# forces (or the MC energy) and resets the step counter and the logged
# history, so the simulation (and its plots) restart cleanly from the
# relaxed configuration, exactly as `minimize!` followed by `md`/`mc` does
# in the tutorial. The energies before/after (the same figures `minimize!`
# itself prints to the terminal) are mirrored into the on-screen status log.
function minimize_now!(obs::Observable{SimState}, log_obs::Observable{Vector{String}})
    s = obs[]
    s.stop = true
    ubefore = potential(s.x, s.sys, s.opt)
    tmp_sys = System(n=s.sys.n, x0=s.x, sides=s.sys.sides)
    minimize!(tmp_sys, s.opt)
    uafter = potential(s.x, s.sys, s.opt)
    log!(log_obs, "Energy before minimization: $ubefore")
    log!(log_obs, "Energy after minimization: $uafter")
    if s.kind == :mc
        s.ucurrent = uafter
        u0, k0 = uafter, 0.0
    else
        forces!(s.f, s.x, s.sys, s.opt)
        s.flast .= s.f
        u0, k0 = uafter, kinetic(s.v)
    end
    s.step = 0
    s.run_start_step = 0
    s.time = 0.0
    empty!(s.steps_history)
    empty!(s.potential_history)
    empty!(s.kinetic_history)
    empty!(s.total_history)
    empty!(s.temperature_history)
    push!(s.steps_history, 0)
    push!(s.potential_history, u0)
    push!(s.kinetic_history, k0)
    push!(s.total_history, u0 + k0)
    push!(s.temperature_history, k0 / s.sys.n)
    obs[] === s && (obs[] = s)
    return nothing
end

function particles_title(s::SimState)
    label = KIND_LABELS[s.kind]
    # Target step of the current (or most recent) run, not just `nsteps` on
    # its own, since `nsteps` counts from `run_start_step` forward.
    target = s.run_start_step + s.nsteps
    # Fixed-width step/time/acceptance fields, so the title doesn't jitter
    # horizontally as the digit count changes from step to step.
    stepstr = lpad(s.step, ndigits(max(target, 1)))
    if s.kind == :mc
        ar = s.step == 0 ? 0.0 : 100 * s.naccepted / s.step
        return "$label  |  step $stepstr/$target  |  acceptance = $(@sprintf("%5.1f", ar))%"
    else
        return "$label  |  step $stepstr/$target  |  t = $(@sprintf("%8.2f", s.time))"
    end
end

# How many trailing steps the temperature-plot title averages over.
const TEMP_AVG_WINDOW = 200

function temperature_title(s::SimState)
    hist = s.temperature_history
    window = @view hist[max(1, end - TEMP_AVG_WINDOW + 1):end]
    avgT = isempty(window) ? 0.0 : sum(window) / length(window)
    return "Temperature (<K>/N) - Last $TEMP_AVG_WINDOW steps: $(@sprintf("%.2f", avgT))"
end

# The 8 first periodic images surrounding the primary cell (edges + corners).
const IMAGE_OFFSETS = ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1))

const PRIMARY_COLOR = to_color((:dodgerblue, 0.85))
const IMAGE_COLOR = to_color((:dodgerblue, 0.3))

# Builds the list of scatter points to display (the particles inside the
# primary cell, plus, if `show`, all 8 of their first periodic images) and
# the matching per-point colors (images are drawn more transparent).
function periodic_points_colors(x, sides, show::Bool)
    hx, hy = sides[1], sides[2]
    n = length(x)
    npts = show ? 9n : n
    pts = Vector{Point2f}(undef, npts)
    cols = Vector{typeof(PRIMARY_COLOR)}(undef, npts)
    for (i, p) in enumerate(x)
        pts[i] = Point2f(image(p, sides))
        cols[i] = PRIMARY_COLOR
    end
    if show
        k = n
        for (ox, oy) in IMAGE_OFFSETS
            for i in 1:n
                k += 1
                b = image(x[i], sides)
                pts[k] = Point2f(b[1] + ox * hx, b[2] + oy * hy)
                cols[k] = IMAGE_COLOR
            end
        end
    end
    return pts, cols
end

function FundamentosDMC.simulate_gui(; n::Int=100, sides=(100.0, 100.0), kind::Symbol=:md)
    state = SimState(; kind, n, Lx=Float64(sides[1]), Ly=Float64(sides[2]))
    obs = Observable(state)

    # Independent view-only toggle: showing the first periodic images does
    # not change the physics, so it must not restart the simulation.
    show_periodic = Observable(false)

    # Status log (last 4 lines), independent of `SimState` so it survives
    # restarts instead of being cleared every time a parameter changes.
    log_obs = Observable(String[])

    GLMakie.activate!(title="FundamentosDMC - Interactive simulation")
    fig = Figure(size=(1400, 820), fontsize=CTRL_FONTSIZE, figure_padding=(10, 10, 10, 30))

    controls = fig[1, 1] = GridLayout(tellwidth=false, valign=:top)
    rowgap!(controls, 3)
    colgap!(controls, 8)

    row = 0
    next_row!() = (row += 1; row)

    #
    # Simulation type and initial velocities menus
    #
    r = next_row!()
    Label(controls[r, 1], "Type", halign=:right, fontsize=CTRL_FONTSIZE)
    kind_menu = Menu(controls[r, 2], options=first.(KIND_OPTIONS), default=KIND_LABELS[state.kind],
        fontsize=CTRL_FONTSIZE, height=CTRL_HEIGHT)
    on(kind_menu.selection) do s
        log!(log_obs, "Changed: type = $s.")
        update_params!(obs, :kind, Dict(KIND_OPTIONS)[s])
    end

    r = next_row!()
    Label(controls[r, 1], "Velocities", halign=:right, fontsize=CTRL_FONTSIZE)
    iv_menu = Menu(controls[r, 2], options=first.(VELOCITY_OPTIONS), default=VELOCITY_LABELS[state.initial_velocities],
        fontsize=CTRL_FONTSIZE, height=CTRL_HEIGHT)
    on(iv_menu.selection) do s
        log!(log_obs, "Changed: velocities = $s.")
        update_params!(obs, :initial_velocities, Dict(VELOCITY_OPTIONS)[s])
    end

    #
    # Numeric parameters (labels match the `Options`/`System` field names)
    #
    # Commits whatever is currently typed (whether or not Enter was
    # pressed) as soon as the textbox loses focus — including when focus
    # moves to another field or a button (e.g. Run) is clicked — instead of
    # requiring Enter and discarding the input otherwise.
    function add_numeric_row!(label, field::Symbol, valtype::Type)
        r = next_row!()
        Label(controls[r, 1], label, halign=:right, fontsize=CTRL_FONTSIZE)
        tb = Textbox(
            controls[r, 2];
            placeholder=@lift(string(getfield($obs, field))),
            validator=valtype,
            fontsize=CTRL_FONTSIZE,
            height=CTRL_HEIGHT,
            width=CTRL_TEXTBOX_WIDTH,
            textpadding=(6, 6, 3, 3),
        )
        on(tb.focused) do focused
            focused && return nothing
            value = tryparse(valtype, tb.displayed_string[])
            value === nothing && return nothing
            value == getfield(obs[], field) && return nothing
            log!(log_obs, "Changed: $field = $value.")
            update_params!(obs, field, value)
            return nothing
        end
        return tb
    end

    add_numeric_row!("n", :n, Int)
    add_numeric_row!("Lx", :Lx, Float64)
    add_numeric_row!("Ly", :Ly, Float64)
    add_numeric_row!("dt", :dt, Float64)
    add_numeric_row!("nsteps", :nsteps, Int)
    add_numeric_row!("eps", :eps, Float64)
    add_numeric_row!("sig", :sig, Float64)
    add_numeric_row!("kT", :kT, Float64)
    add_numeric_row!("ibath", :ibath, Int)
    add_numeric_row!("iequil", :iequil, Int)
    add_numeric_row!("tau", :tau, Int)
    add_numeric_row!("lambda", :lambda, Float64)
    add_numeric_row!("alpha", :alpha, Float64)

    #
    # Periodic-image display toggle (view-only, does not affect the physics)
    #
    r = next_row!()
    Label(controls[r, 1], "Images", halign=:right, fontsize=CTRL_FONTSIZE)
    periodic_cb = Checkbox(controls[r, 2], checked=false)
    on(periodic_cb.checked) do v
        show_periodic[] = v
    end

    #
    # Minimize / Run / Stop / Reset buttons, side by side on a single row
    # (Minimize first, since it is normally the first thing to do).
    #
    r = next_row!()
    button_grid = controls[r, 1:2] = GridLayout()
    colgap!(button_grid, 4)
    buttons = button_grid[1, 1:4] = [
        Button(fig, label="Minimize", fontsize=CTRL_FONTSIZE, height=CTRL_HEIGHT,
            buttoncolor=:seagreen, labelcolor=:white),
        Button(fig, label="Run", fontsize=CTRL_FONTSIZE, height=CTRL_HEIGHT,
            buttoncolor=:gold, labelcolor=:black),
        Button(fig, label="Stop", fontsize=CTRL_FONTSIZE, height=CTRL_HEIGHT,
            buttoncolor=:firebrick, labelcolor=:white),
        Button(fig, label="Reset", fontsize=CTRL_FONTSIZE, height=CTRL_HEIGHT,
            buttoncolor=:thistle, labelcolor=:black),
    ]
    on(buttons[1].clicks) do _
        minimize_now!(obs, log_obs)
    end
    on(buttons[2].clicks) do _
        obs[].running || @async run!(obs, log_obs)
    end
    on(buttons[3].clicks) do _
        obs[].stop = true
    end
    on(buttons[4].clicks) do _
        log!(log_obs, "Reset: new random configuration.")
        reset!(obs)
    end

    colsize!(controls, 1, Fixed(70))
    colsize!(controls, 2, Fixed(150))

    #
    # Visualization: particle motion, energies, and temperature
    #
    viz = fig[1, 2] = GridLayout()

    ax_particles = Axis(
        viz[1:2, 1],
        aspect=DataAspect(),
        title=@lift(particles_title($obs)),
        titlealign=:left,
        xticks=@lift([-$(obs).sys.sides[1] / 2, $(obs).sys.sides[1] / 2]),
        yticks=@lift([-$(obs).sys.sides[2] / 2, $(obs).sys.sides[2] / 2]),
    )
    boxpts = @lift(let sides = $(obs).sys.sides
        hx, hy = sides[1] / 2, sides[2] / 2
        Point2f[(-hx, -hy), (hx, -hy), (hx, hy), (-hx, hy), (-hx, -hy)]
    end)
    lines!(ax_particles, boxpts, color=:black)
    pts_and_colors = @lift(periodic_points_colors($(obs).x, $(obs).sys.sides, $show_periodic))
    scatter!(
        ax_particles, @lift($pts_and_colors[1]);
        markersize=@lift(Float32($(obs).sig)),
        markerspace=:data,
        color=@lift($pts_and_colors[2]),
        strokewidth=1,
        strokecolor=:black,
    )

    # The view box always extends half a box-width/height beyond the
    # primary cell on each side (i.e. ±Lx, ±Ly — showing half of each
    # neighboring image, not the full one), whether or not periodic images
    # are currently displayed, so it never changes when the toggle is
    # switched — it only depends on the box size, i.e. it is only
    # recomputed on restart.
    last_sides = Ref((NaN, NaN))
    function update_view_limits!(s)
        sides = (s.sys.sides[1], s.sys.sides[2])
        sides == last_sides[] && return nothing
        last_sides[] = sides
        hx, hy = sides[1], sides[2]
        limits!(ax_particles, -hx, hx, -hy, hy)
        return nothing
    end
    update_view_limits!(obs[])
    on(update_view_limits!, obs)

    ax_energy = Axis(viz[1, 2], xlabel="step", ylabel="Energy", title="Potential, kinetic and total energy")
    lines!(ax_energy, @lift(Point2f.($(obs).steps_history, $(obs).potential_history)), color=:royalblue, label="Potential")
    lines!(ax_energy, @lift(Point2f.($(obs).steps_history, $(obs).kinetic_history)),
        color=:seagreen, label="Kinetic", visible=@lift($(obs).kind != :mc))
    lines!(ax_energy, @lift(Point2f.($(obs).steps_history, $(obs).total_history)),
        color=:firebrick, label="Total", visible=@lift($(obs).kind != :mc))
    axislegend(ax_energy, position=:lb)

    ax_temp = Axis(viz[2, 2], xlabel="step", ylabel="Temperature", title=@lift(temperature_title($obs)))
    lines!(ax_temp, @lift(Point2f.($(obs).steps_history, $(obs).temperature_history)),
        color=:purple, visible=@lift($(obs).kind != :mc))
    hlines!(ax_temp, @lift(Float32($(obs).kT)), color=:gray, linestyle=:dash)

    on(obs) do _
        autolimits!(ax_energy)
        autolimits!(ax_temp)
    end

    colsize!(viz, 1, Relative(0.62))
    colsize!(viz, 2, Relative(0.38))
    colsize!(fig.layout, 1, Fixed(250))
    colsize!(fig.layout, 2, Auto())

    #
    # Status log: a small "console" panel at the bottom (last 4 lines)
    # mirroring what is happening under the hood (minimization output,
    # run/stop/explosion notices, restarts), so it doesn't only go to the
    # terminal.
    #
    status_panel = fig[2, 1:2] = GridLayout()
    Box(status_panel[1, 1], color=(:black, 0.05), strokecolor=(:black, 0.3), strokewidth=1)
    Label(
        status_panel[1, 1],
        @lift(isempty($log_obs) ? "Status: waiting (Run, Minimize, Reset, or a parameter change will show output here)..." : join($log_obs, "\n"));
        fontsize=CTRL_FONTSIZE,
        halign=:left,
        valign=:top,
        justification=:left,
        padding=(10, 10, 8, 8),
        tellwidth=false,
        tellheight=false,
    )
    rowgap!(fig.layout, 16)
    rowsize!(fig.layout, 1, Fixed(575))
    rowsize!(fig.layout, 2, Fixed(130))

    return fig
end

end # module
