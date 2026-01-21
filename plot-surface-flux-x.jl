# Load the necessary packages
using Breeze
using Oceananigans
using Oceananigans.Units
using Printf
using CairoMakie
using Statistics

# setup for plotting
fig = Figure(size=(600, 1000), fontsize=13)
axτ = Axis(fig[1, 1], ylabel="τ (kg m⁻¹ s⁻²)")
axQ = Axis(fig[2, 1], xlabel="x (km)", ylabel="𝒬 (W m⁻²)")
axQt = Axis(fig[3, 1],
            xlabel="time (hours)",
            ylabel="⟨𝒬⟩ₓ (W m⁻²)")

front_filename = "/g/data/fy29/mxr581/Breeze/front-2K-2K-U7up-80km/prescribed_sea_surface_temperature_convection_surface_fluxes1d.jld2"
ctrl_filename  = "/g/data/fy29/mxr581/Breeze/ctrl-2K-U7up-80km/prescribed_sea_surface_temperature_convection_surface_fluxes1d.jld2"

# Front / experiment
τˣ_x  = FieldTimeSeries(front_filename, "τˣ_avg")
τʸ_x  = FieldTimeSeries(front_filename, "τʸ_avg")
𝒬ᵀ_x = FieldTimeSeries(front_filename, "𝒬ᵀ_avg")
𝒬ᵛ_x = FieldTimeSeries(front_filename, "𝒬ᵛ_avg")

# Control
τˣ_x_ctrl  = FieldTimeSeries(ctrl_filename, "τˣ_avg")
τʸ_x_ctrl  = FieldTimeSeries(ctrl_filename, "τʸ_avg")
𝒬ᵀ_x_ctrl = FieldTimeSeries(ctrl_filename, "𝒬ᵀ_avg")
𝒬ᵛ_x_ctrl = FieldTimeSeries(ctrl_filename, "𝒬ᵛ_avg")

times = τˣ_x.times
N = length(τˣ_x)
n = Observable(N)

τˣxn = @lift τˣ_x[$n]
τʸxn = @lift τʸ_x[$n]
𝒬ᵀxn = @lift 𝒬ᵀ_x[$n]
𝒬ᵛxn = @lift 𝒬ᵛ_x[$n]

τˣxn_ctrl  = @lift τˣ_x_ctrl[$n]
τʸxn_ctrl  = @lift τʸ_x_ctrl[$n]
𝒬ᵀxn_ctrl = @lift 𝒬ᵀ_x_ctrl[$n]
𝒬ᵛxn_ctrl = @lift 𝒬ᵛ_x_ctrl[$n]

# x-mean sensible & latent heat flux
Qᵀ_mean = [mean(𝒬ᵀ_x[n]) for n in 1:N]
Qᵛ_mean = [mean(𝒬ᵛ_x[n]) for n in 1:N]
Qᵀ_mean_ctrl = [mean(𝒬ᵀ_x_ctrl[n]) for n in 1:N]
Qᵛ_mean_ctrl = [mean(𝒬ᵛ_x_ctrl[n]) for n in 1:N]

title = @lift "t = $(prettytime(times[$n]))"
fig[0, :] = Label(fig, title, fontsize=22, tellwidth=false)

xm, ym, zm = nodes(τˣ_x)
x = xm ./ kilometer

# --- Stress ---
lines!(axτ, x, τˣxn, color=:midnightblue, linewidth=2, label="τˣ")
lines!(axτ, x, τʸxn, color=:royalblue,  linewidth=2, label="τʸ")

lines!(axτ, x, τˣxn_ctrl, color=:midnightblue, linewidth=2, linestyle=:dash)
lines!(axτ, x, τʸxn_ctrl, color=:royalblue,  linewidth=2, linestyle=:dash)

# --- Heat fluxes ---
lines!(axQ, x, 𝒬ᵀxn, color=:firebrick,  linewidth=2, label="sensible")
lines!(axQ, x, 𝒬ᵛxn, color=:goldenrod, linewidth=2, label="latent")

lines!(axQ, x, 𝒬ᵀxn_ctrl, color=:firebrick,  linewidth=2, linestyle=:dash)
lines!(axQ, x, 𝒬ᵛxn_ctrl, color=:goldenrod, linewidth=2, linestyle=:dash)

Legend(fig[2, 2], axQ)
Legend(fig[1, 2], axτ)

axτ.xlabel = "x (km)"
axQ.xlabel = "x (km)"

# Front / experiment
t_hours = times ./ hour
lines!(axQt, t_hours, Qᵀ_mean, color=:firebrick, linewidth=2, label="sensible")
lines!(axQt, t_hours, Qᵛ_mean, color=:goldenrod, linewidth=2, label="latent")

# Control (dashed)
lines!(axQt, t_hours, Qᵀ_mean_ctrl, color=:firebrick, linewidth=2, linestyle=:dash)
lines!(axQt, t_hours, Qᵛ_mean_ctrl, color=:goldenrod, linewidth=2, linestyle=:dash)
Legend(fig[3, 2], axQt)

ylims!(axτ, minimum(τʸ_x), maximum(τˣ_x))
ylims!(axQ, -10, maximum(𝒬ᵛ_x))

CairoMakie.record(fig, "fluxes_compare.mp4", 1:N, framerate=12) do nn
    n[] = nn
end




