# Load the necessary packages
using Breeze
using Oceananigans
using Oceananigans.Units
using Printf
using CairoMakie
using Statistics

# setup for plotting
fig = Figure(size=(600, 600), fontsize=13)
axτ = Axis(fig[1, 1], ylabel="⟨τ⟩ (kg m⁻¹ s⁻²)")
axQV = Axis(fig[2, 1], ylabel="⟨𝒬ᵛ⟩ (W m⁻²)")
axQT = Axis(fig[3, 1], xlabel="time (hours)", ylabel="⟨𝒬ᵀ⟩ (W m⁻²)")

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

# x-mean sensible & latent heat flux
Qᵀ_mean = [mean(𝒬ᵀ_x[n]) for n in 1:N]
Qᵛ_mean = [mean(𝒬ᵛ_x[n]) for n in 1:N]
Qᵀ_mean_ctrl = [mean(𝒬ᵀ_x_ctrl[n]) for n in 1:N]
Qᵛ_mean_ctrl = [mean(𝒬ᵛ_x_ctrl[n]) for n in 1:N]

# x-mean stress 
τˣ_mean = [mean(τˣ_x[n]) for n in 1:N]
τʸ_mean = [mean(τʸ_x[n]) for n in 1:N]
τˣ_mean_ctrl = [mean(τˣ_x_ctrl[n]) for n in 1:N]
τʸ_mean_ctrl = [mean(τʸ_x_ctrl[n]) for n in 1:N]

xm, ym, zm = nodes(τˣ_x)
x = xm ./ kilometer

# Front / experiment
t_hours = times ./ hour
lines!(axQT, t_hours, Qᵀ_mean, color=:firebrick, linewidth=2, label="sensible")
lines!(axQV, t_hours, Qᵛ_mean, color=:goldenrod, linewidth=2, label="latent")

# Control (dashed)
lines!(axQT, t_hours, Qᵀ_mean_ctrl, color=:firebrick, linewidth=2, linestyle=:dash)
lines!(axQV, t_hours, Qᵛ_mean_ctrl, color=:goldenrod, linewidth=2, linestyle=:dash)

# Front / Experiment 
lines!(axτ, t_hours, τˣ_mean, color=:midnightblue, linewidth=2,  label="τˣ")
lines!(axτ, t_hours, τʸ_mean, color=:royalblue, linewidth=2,  label="τʸ")

# Control (dashed)
lines!(axτ, t_hours, τˣ_mean_ctrl, color=:midnightblue, linewidth=2, linestyle=:dash)
lines!(axτ, t_hours, τʸ_mean_ctrl, color=:royalblue, linewidth=2, linestyle=:dash)

Legend(fig[1, 2], axτ)

# ylims!(axτ, minimum(τʸ_x), maximum(τˣ_x))
# ylims!(axτ, -10, maximum(𝒬ᵛ_x))
save("fluxes_compare_timeseries.png", fig)




