# Load the necessary packages
using Breeze
using Oceananigans
using Oceananigans.Units
using Printf
using CairoMakie
using Statistics

filename_snap_fluxes = "/g/data/fy29/mxr581/Breeze/front-2K-2K-U7up-80km/prescribed_sea_surface_temperature_convection_surface_fluxes.jld2"
#filename_snap_planes = "/g/data/fy29/mxr581/Breeze/front-2K-2K-U7up-80km/prescribed_sea_surface_temperature_convection_slice.jld2"
filename_snap_planes = "./SST_front_snap_planes.jld2"

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

fig = Figure(size=(1200, 800), fontsize=16)

title = @lift "t = $(prettytime(times[$n]))"

axw = Axis(fig[1, 1], ylabel="z (km)")
axu = Axis(fig[1, 3])
axv = Axis(fig[1, 5])
axθ = Axis(fig[2, 1], ylabel="z (km)")
axq = Axis(fig[2, 3])
axql = Axis(fig[2, 5])
axτ = Axis(fig[3, 1], xlabel="x (km)", ylabel="y (km)")
ax𝒬 = Axis(fig[3, 3], xlabel="x (km)")
axV = Axis(fig[3, 5], xlabel = "x (km)")

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

# xm = Oceananigans.Grids.xnodes(grid, Center()) ./ 1000  # Convert to km
# ym = Oceananigans.Grids.ynodes(grid, Center()) ./ 1000  # Convert to km
# zm = Oceananigans.Grids.znodes(grid, Center()) ./ 1000  # Convert to km

xm, ym, zm = nodes(w_ts)
x = xm ./ kilometer
z = zm ./ kilometer
y = ym ./ kilometer

hmu = heatmap!(axu, x, z, un, colorrange=u_limits, colormap=:speed)
hmv = heatmap!(axv, x, z, vn, colorrange=v_limits, colormap=:speed)
hmw = heatmap!(axw, x, z, wn, colorrange=w_limits, colormap=:balance)
hmθ = heatmap!(axθ, x, z, θn, colorrange=θ_limits, colormap=:thermal)
hmq = heatmap!(axq, x, z, qᵗn, colorrange=qᵗ_limits, colormap=Reverse(:Purples_4))
hmql = heatmap!(axql, x, z, qˡn, colorrange=(0.0, 0.00001), colormap=:deep)
hmτ = heatmap!(axτ, x, y, τˣn, colorrange=(-τˣ_max, τˣ_max), colormap=:curl)
hm𝒬 = heatmap!(ax𝒬, x, y, 𝒬ᵀn, colorrange=(𝒬_min , 𝒬_max))
hmV = heatmap!(axV, x, y, 𝒬ᵛn, colorrange=(𝒬V_min , 𝒬V_max))

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