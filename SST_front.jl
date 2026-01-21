# Atmospheric flow over a prescribed SST front using Breeze.jl

using Breeze
using Oceananigans
using Oceananigans.Units
using Oceananigans.Models: BoundaryConditionOperation
using Printf
using CUDA
using CairoMakie
using AtmosphericProfilesLibrary
using Statistics

# until https://github.com/CliMA/Oceananigans.jl/pull/5100 is merged
using Adapt: Adapt, adapt

function Adapt.adapt_structure(to, bckf::Oceananigans.Models.BoundaryConditionKernelFunction{Side}) where Side
    bc = adapt(to, bckf.bc)
    BC = typeof(bc)
    return Oceananigans.Models.BoundaryConditionKernelFunction{Side, BC}(bc)
end


# Grid setup
grid = RectilinearGrid(GPU(), size = (512, 256, 100), halo = (5, 5, 5),
                       x = (-40kilometers, 40kilometers),
                       y = (-20kilometers, 20kilometers),
                       z = (0, 3kilometers),
                       topology = (Periodic, Periodic, Bounded))

# Model formulation
p₀, θ₀ = 101325, 285 # Pa, K
constants = ThermodynamicConstants()
reference_state = ReferenceState(grid, constants; surface_pressure=p₀, potential_temperature=θ₀)
dynamics = AnelasticDynamics(reference_state)
microphysics = SaturationAdjustment(equilibrium = WarmPhaseEquilibrium())
momentum_advection = WENO(order=9)
scalar_advection = WENO(order=5)

FT = eltype(grid)
coriolis = FPlane(f=8.3e-5)

uᵍ(z) = - 7 - 0.0015 * z 
vᵍ(z) = 0.0 #

geostrophic = geostrophic_forcings(z -> uᵍ(z), z -> vᵍ(z))

ρu_forcing = geostrophic.ρu
ρv_forcing = geostrophic.ρv

forcing = (; ρu=ρu_forcing, ρv=ρv_forcing)

# Boundary conditions
Cᵀ = 1e-3  # Sensible heat transfer coefficient
Cᵛ = 1e-3  # Vapor transfer coefficient
Cᴰ = 1e-3  # Drag coefficient
Uᵍ = 0.1  # Minimum wind speed (m/s)

ρu_surface_flux = BulkDrag(; coefficient=Cᴰ, gustiness=Uᵍ)
ρv_surface_flux = BulkDrag(; coefficient=Cᴰ, gustiness=Uᵍ)

SST = 289 # mean sea surface temperature in K
ΔT = 2.0 # front amplitude
steepness = 10 # Controls sharpness of the transition
T₀(x, y)  = SST + ΔT / 2 * tanh(steepness * cos(2π * x / grid.Lx))

# and build the flux parameterizations
ρθ_surface_flux = BulkSensibleHeatFlux(coefficient=Cᵀ, gustiness=Uᵍ, surface_temperature=T₀)
ρqᵗ_surface_flux = BulkVaporFlux(coefficient=Cᵛ, gustiness=Uᵍ, surface_temperature=T₀)

# assemble the boundary conditions
ρu_bcs = FieldBoundaryConditions(bottom=ρu_surface_flux)
ρv_bcs = FieldBoundaryConditions(bottom=ρv_surface_flux)
ρθ_bcs = FieldBoundaryConditions(bottom=ρθ_surface_flux)
ρqᵗ_bcs = FieldBoundaryConditions(bottom=ρqᵗ_surface_flux)

# Model construction
model = AtmosphereModel(grid; momentum_advection, scalar_advection, coriolis, microphysics, dynamics,
                        forcing, boundary_conditions = (ρu=ρu_bcs, ρv=ρv_bcs, ρθ=ρθ_bcs, ρqᵗ=ρqᵗ_bcs))

# Initial conditions
u₀ = AtmosphericProfilesLibrary.Bomex_u(FT)

δθ = 0.1      # K
δu = 0.2  # kg/kg
δq = 0.001  # kg/kg
zδ = 600     # m

ϵ() = rand() - 1/2
uᵢ(x, y, z) = uᵍ(z) + δu * ϵ() * (z < zδ)

# initial profiles theta and q
θᵇ(z) = z < 600  ? 286.68 + 0.0029*z :
        z < 1500 ? 281.42 + 0.0117*z :
                   290.20 + 0.0059*z

qᵇ(z) = z < 600  ? 0.0088 - 3.03e-6*z :
        z < 1500 ? 0.0072 - 3.98e-7*z :
                   0.0098 - 2.12e-6*z

θᵢ(x, y, z) = θᵇ(z) + δθ * ϵ() * (z < zδ)
qᵢ(x, y, z) = qᵇ(z) + δq * ϵ() * (z < zδ)

set!(model, θ=θᵢ, u=uᵢ, qᵗ=qᵢ)

# Simulation setup
simulation = Simulation(model, Δt=2, stop_time=8hours)
conjure_time_step_wizard!(simulation, cfl=0.7)

#  Diagnostic fields
T = model.temperature
θ = liquid_ice_potential_temperature(model)
qˡ = model.microphysical_fields.qˡ
qᵛ⁺ = Breeze.Microphysics.SaturationSpecificHumidity(model)

ρu, ρv, ρw = model.momentum
u, v, w = model.velocities
qᵗ = model.specific_moisture

# Surface flux diagnostics

## Surface momentum flux
τˣ = BoundaryConditionOperation(ρu, :bottom, model)
τʸ = BoundaryConditionOperation(ρv, :bottom, model)

## Sensible heat flux: 𝒬ᵀ = cᵖᵈ Jᵀ (using dry air heat capacity as approximation)
ρθ = liquid_ice_potential_temperature_density(model)
cᵖᵈ = constants.dry_air.heat_capacity
Jᵀ = BoundaryConditionOperation(ρθ, :bottom, model)
𝒬ᵀ = cᵖᵈ * Jᵀ

## Latent heat flux: 𝒬ᵛ = ℒˡ Jᵛ (using reference θ₀ for latent heat)
ρqᵗ = model.moisture_density
ℒˡ = Breeze.Thermodynamics.liquid_latent_heat(θ₀, constants)
Jᵛ = BoundaryConditionOperation(ρqᵗ, :bottom, model)
𝒬ᵛ = ℒˡ * Jᵛ

# Progress callback
function progress(sim)
    qᵗ = sim.model.specific_moisture
    u, v, w = sim.model.velocities

    umax = maximum(abs, u)
    vmax = maximum(abs, v)
    wmax = maximum(abs, w)

    qᵗmin = minimum(qᵗ)
    qᵗmax = maximum(qᵗ)
    qˡmax = maximum(qˡ)

    θmin = minimum(θ)
    θmax = maximum(θ)

    msg = @sprintf("Iter: %d, t = %s, max|u|: (%.2e, %.2e, %.2e)",
                    iteration(sim), prettytime(sim), umax, vmax, wmax)

    msg *= @sprintf(", extrema(qᵗ): (%.2e, %.2e), max(qˡ): %.2e, extrema(θ): (%.2e, %.2e)",
                     qᵗmin, qᵗmax, qˡmax, θmin, θmax)

    @info msg

    return nothing
end

add_callback!(simulation, progress, IterationInterval(100))

# Output

filename_snap_planes = "SST_front_snap_planes.jld2"
filename_ymean_planes = "SST_front_ymean_planes.jld2"
filename_ymean_fluxes = "SST_front_ymean_fluxes.jld2"
filename_snap_fluxes = "SST_front_snap_fluxes.jld2"
filename_mean_profiles = "SST_front_mean_profiles.jld2"

qᵗ = model.specific_moisture
u, v, w, = model.velocities
#s = sqrt(u^2 + w^2) # speed
#ξ = ∂z(u) - ∂x(w)   # cross-stream vorticity
U = mean(u, dims=(1, 2))  # horizontal mean
V = mean(v, dims=(1, 2))
W = mean(w, dims=(1, 2))
u′² = (u - U) * (u - U)
v′² = (v - V) * (v - V)
w′² = (w - W) * (w - W)

# 2d surface fluxes
simulation.output_writers[:fluxes2d] = JLD2Writer(model, (; τˣ, τʸ, 𝒬ᵀ, 𝒬ᵛ);
                filename = filename_snap_fluxes,
                schedule = TimeInterval(2minutes),
                overwrite_existing = true)

# y-averaged suface fluxes
τˣ_avg = Average(τˣ, dims= (2))
τʸ_avg = Average(τʸ, dims= (2))
𝒬ᵀ_avg = Average(𝒬ᵀ, dims= (2))
𝒬ᵛ_avg = Average(𝒬ᵛ, dims= (2))

simulation.output_writers[:fluxes1d] = JLD2Writer(model, (; τˣ_avg, τʸ_avg, 𝒬ᵀ_avg, 𝒬ᵛ_avg);
                filename = filename_ymean_fluxes,
                schedule = TimeInterval(2minutes),
                overwrite_existing = true)


# xz slices at y = 0 and xy slices at z = 500 m
z = Oceananigans.Grids.znodes(grid, Center())
k = searchsortedfirst(z, 500)

outputs_snap_planes = (
    uxz = view(u, :, 1, :),
    vxz = view(v, :, 1, :),
    wxz = view(w, :, 1, :),
    qᵗxz = view(qᵗ, :, 1, :),
    qˡxz = view(qˡ, :, 1, :),
    θxz = view(θ, :, 1, :),
    uxy = view(u, :, k, :),
    vxy = view(v, :, k, :),
    wxy = view(w, :, k, :),
    qᵗxy = view(qᵗ, :, k, :),
    qˡxy = view(qˡ, :, k, :),
    θxy = view(θ, :, k, :),
)

simulation.output_writers[:planes] = JLD2Writer(model, outputs_snap_planes;
                                                filename = filename_snap_planes,
                                                schedule = TimeInterval(2minutes),
                                                overwrite_existing = true)
                                            
# y-mean xz slices
outputs_mean_planes = (
    uxz_mean = Average(u, dims = 2),
    vxz_mean = Average(v, dims = 2),
    wxz_mean = Average(w, dims = 2),
    uvar_mean = Average(u′², dims = 2),
    vvar_mean = Average(v′², dims = 2),
    wvar_mean = Average(w′², dims = 2),
    qᵗxz_mean = Average(qᵗ, dims = 2),
    qˡxz_mean = Average(qˡ, dims = 2),
    θxz_mean = Average(θ, dims = 2),
)

simulation.output_writers[:mean_planes] = JLD2Writer(model, outputs_mean_planes;
                                                filename = filename_ymean_planes,
                                                schedule = TimeInterval(2minutes),
                                                overwrite_existing = true)


# Profiles
outputs_profiles = (
    ū = Average(u, dims=(1, 2)),
    v̄ = Average(v, dims=(1, 2)),
    w̄ = Average(w, dims=(1, 2)),
    θ̄ = Average(θ, dims=(1, 2)),
    q̄ᵗ = Average(qᵗ, dims=(1, 2)),
    q̄ˡ = Average(qˡ, dims=(1, 2)),
)

simulation.output_writers[:profiles] = JLD2Writer(model, outputs_profiles;
                                                 filename = filename_mean_profiles,
                                                 schedule = TimeInterval(2minutes),
                                                 overwrite_existing = true)

## Run the simulation

@info "Running prescribed SST convection simulation..."
run!(simulation)

# ## Visualization

@assert isfile(filename_snap_fluxes) "Output file $(filename_snap_fluxes) not found."

τˣ_ts = FieldTimeSeries(filename_snap_fluxes, "τˣ")
𝒬ᵀ_ts = FieldTimeSeries(filename_snap_fluxes, "𝒬ᵀ")
𝒬ᵛ_ts = FieldTimeSeries(filename_snap_fluxes, "𝒬ᵛ")

w_ts = FieldTimeSeries(filename_snap_planes, "wxz")
u_ts = FieldTimeSeries(filename_snap_planes, "uxz")
v_ts = FieldTimeSeries(filename_snap_planes, "vxz")
qᵗ_ts = FieldTimeSeries(filename_snap_planes, "qᵗxz")
qˡ_ts = FieldTimeSeries(filename_snap_planes, "qˡxz")
θ_ts = FieldTimeSeries(filename_snap_planes, "θxz")

times = τˣ_ts.times
Nt = length(τˣ_ts)
n = Observable(Nt)

wn = @lift w_ts[$n]
un = @lift u_ts[$n]
vn = @lift v_ts[$n]
qᵗn = @lift qᵗ_ts[$n]
qˡn = @lift qˡ_ts[$n]
θn = @lift θ_ts[$n]
τˣn = @lift τˣ_ts[$n]
𝒬ᵀn = @lift 𝒬ᵀ_ts[$n]
𝒬ᵛn = @lift 𝒬ᵛ_ts[$n]

# Plot

fig = Figure(size=(1200, 800), fontsize=13)

title = @lift "t = $(prettytime(times[$n]))"

axw = Axis(fig[1, 1], ylabel="z (km)")
axu = Axis(fig[1, 3])
axv = Axis(fig[1, 5])
axθ = Axis(fig[2, 1], ylabel="z (km)")
axq = Axis(fig[2, 3])
axql = Axis(fig[2, 5])
axτ = Axis(fig[3, 1], xlabel="x (km)", ylabel="y (km)")
ax𝒬 = Axis(fig[3, 3], xlabel="x (km)")
axV = Axis(fig[3, 5])

fig[0, :] = Label(fig, title, fontsize=22, tellwidth=false)

# # Compute color limits from the full time series
θ_limits = extrema(θ_ts)
u_limits = extrema(u_ts)
v_limits = extrema(v_ts)
w_limits = (-maximum(w_ts), maximum(w_ts))
qᵗ_limits = extrema(qᵗ_ts)
qˡ_limits = extrema(qˡ_ts)

# Flux limits
τˣ_max = max(abs(minimum(τˣ_ts)), abs(maximum(τˣ_ts)))
𝒬_min = minimum(𝒬ᵀ_ts)
𝒬_max = maximum(𝒬ᵀ_ts)

𝒬V_min = minimum(𝒬ᵛ_ts)
𝒬V_max = maximum(𝒬ᵛ_ts)

xm = Oceananigans.Grids.xnodes(grid, Center()) ./ 1000  # Convert to km
ym = Oceananigans.Grids.ynodes(grid, Center()) ./ 1000  # Convert to km
zm = Oceananigans.Grids.znodes(grid, Center()) ./ 1000  # Convert to km

hmu = heatmap!(axu, xm, zm, un, colorrange=u_limits, colormap=:speed)
hmv = heatmap!(axv, xm, zm, vn, colorrange=v_limits, colormap=:speed)
hmw = heatmap!(axw, xm, zm, wn, colorrange=w_limits, colormap=:balance)
hmθ = heatmap!(axθ, xm, zm, θn, colorrange=θ_limits, colormap=:thermal)
hmq = heatmap!(axq, xm, zm, qᵗn, colorrange=qᵗ_limits, colormap=Reverse(:Purples_4))
hmql = heatmap!(axql, xm, zm, qˡn, colorrange=(0.0, 0.00001), colormap=:deep)
hmτ = heatmap!(axτ, xm, ym, τˣn, colorrange=(-τˣ_max, τˣ_max), colormap=:curl)
hm𝒬 = heatmap!(ax𝒬, xm, ym, 𝒬ᵀn, colorrange=(𝒬_min , 𝒬_max))
hmV = heatmap!(axV, xm, ym, 𝒬ᵛn, colorrange=(𝒬V_min , 𝒬V_max))

Colorbar(fig[1, 2], hmw, label="w (m/s)")
Colorbar(fig[1, 4], hmu, label="u (m/s)")
Colorbar(fig[1, 6], hmv, label="v (m/s)")
Colorbar(fig[2, 2], hmθ, label="θ (K)")
Colorbar(fig[2, 4], hmq, label="qᵗ (kg/kg)")
Colorbar(fig[2, 6], hmql, label="qˡ (kg/kg)")
Colorbar(fig[3, 2], hmτ, label="τˣ (kg m⁻¹ s⁻²)")
Colorbar(fig[3, 4], hm𝒬,  label="𝒬 sensible (W m⁻²)")
Colorbar(fig[3, 6], hmV,  label="𝒬 latent (W m⁻²)")

# Now we are ready to make a cool animation.

CairoMakie.record(fig, "prescribed_sea_surface_temperature.mp4", 1:Nt, framerate=12) do nn
    n[] = nn
end

# # ![](prescribed_sea_surface_temperature.mp4)

## line plots
τˣ_x = FieldTimeSeries(filename_ymean_fluxes, "τˣ_avg")
τʸ_x = FieldTimeSeries(filename_ymean_fluxes, "τʸ_avg")
𝒬ᵀ_x = FieldTimeSeries(filename_ymean_fluxes, "𝒬ᵀ_avg")
𝒬ᵛ_x = FieldTimeSeries(filename_ymean_fluxes, "𝒬ᵛ_avg")

times = τˣ_x.times
N = length(τˣ_x)
n = Observable(N)

τˣxn = @lift τˣ_x[$n]
τʸxn = @lift τʸ_x[$n]
𝒬ᵀxn = @lift 𝒬ᵀ_x[$n]
𝒬ᵛxn = @lift 𝒬ᵛ_x[$n]

fig = Figure(size=(600, 800), fontsize=13)

title = @lift "t = $(prettytime(times[$n]))"

axτ = Axis(fig[1, 1], ylabel="τ (kg m⁻¹ s⁻²)")
axQ = Axis(fig[2, 1], xlabel="x (km)", ylabel="𝒬 (W m⁻²)")

x = Oceananigans.Grids.xnodes(grid, Center()) ./ 1000  # Convert to km

lines!(axτ, x, τˣxn, color=:midnightblue, linewidth=2, label="τˣ")
lines!(axτ, x, τʸxn, color=:royalblue, linewidth=2, label="τʸ")
lines!(axQ, x, 𝒬ᵀxn, color=:firebrick, linewidth=2, label="sensible")
lines!(axQ, x, 𝒬ᵛxn, color=:goldenrod, linewidth=2, label="latent")
Legend(fig[2, 2], axQ)
Legend(fig[1, 2], axτ)

axτ.xlabel = "x (km)"
axQ.xlabel = "x (km)"

ylims!(axτ, minimum(τʸ_x), maximum(τˣ_x))
ylims!(axQ, -10, maximum(𝒬ᵛ_x))

CairoMakie.record(fig, "prescribed_sea_surface_temperature_lines.mp4", 1:Nt, framerate=12) do nn
    n[] = nn
end
