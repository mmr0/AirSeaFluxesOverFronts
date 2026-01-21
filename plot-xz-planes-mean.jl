using Breeze
using Oceananigans
using Oceananigans.Units
using Printf
using CairoMakie
using Statistics

#output_filename = "/g/data/fy29/mxr581/Breeze/front-2K-2K-U7up-80km/prescribed_sea_surface_temperature_convection_slice_ymean.jld2"
output_filename = "./prescribed_sea_surface_temperature_convection_slice_ymean.jld2"

w_ts = FieldTimeSeries(output_filename, "wxz_mean")
wvar_ts = FieldTimeSeries(output_filename, "wvar_mean")
u_ts = FieldTimeSeries(output_filename, "uxz_mean")
v_ts = FieldTimeSeries(output_filename, "vxz_mean")
qᵗ_ts = FieldTimeSeries(output_filename, "qᵗxz_mean")
qˡ_ts = FieldTimeSeries(output_filename, "qˡxz_mean")
θ_ts = FieldTimeSeries(output_filename, "θxz_mean")

times = w_ts.times
t_end = times[end]
t_start = t_end - 1hour

last_hour_inds = findall(t -> t ≥ t_start, times)

function time_mean(ts::FieldTimeSeries, inds)
    f̄ = deepcopy(ts[inds[1]])
    parent(f̄) .= 0

    for i in inds
        parent(f̄) .+= parent(ts[i])
    end

    parent(f̄) ./= length(inds)
    return f̄
end
    
wₙ  = time_mean(w_ts, last_hour_inds)
wvarₙ  = time_mean(wvar_ts, last_hour_inds)
uₙ  = time_mean(u_ts, last_hour_inds)
vₙ  = time_mean(v_ts, last_hour_inds)
θₙ  = time_mean(θ_ts, last_hour_inds)
qᵗₙ = time_mean(qᵗ_ts, last_hour_inds)
qˡₙ = time_mean(qˡ_ts, last_hour_inds)



# times = τˣ_ts.times
# Nt = length(τˣ_ts)
# n = Observable(Nt)

# setup for plotting
fig = Figure(size=(1400, 800), fontsize=13)

axw  = Axis(fig[1, 1], ylabel="z (km)", title="w'w' [m² s⁻²]")
axu  = Axis(fig[1, 3], title="u [m s⁻¹]")
axv  = Axis(fig[1, 5], title="v [m s⁻¹]")
axθ  = Axis(fig[2, 1], ylabel="z (km)", xlabel="x (km)", title="θ [K]")
axq  = Axis(fig[2, 3], xlabel="x (km)", title="qᵗ [kg kg⁻¹]")
axql = Axis(fig[2, 5], xlabel="x (km)", title="qˡ [kg kg⁻¹]")

xm, ym, zm = nodes(w_ts)
x = xm ./ kilometer
z = zm ./ kilometer

zmax = 1.5
kz = findall(z .≤ zmax)

function z_limited_extrema(field, kz)
    data = parent(field)[:, 1, kz]   # x, y=1, z subset, t=1
    return extrema(data)
end

u_limits  = z_limited_extrema(uₙ,  kz)
v_limits  = z_limited_extrema(vₙ,  kz)
θ_limits  = (minimum(θₙ[1,1,kz]), maximum(θₙ[1,1,kz]))
qᵗ_limits = (minimum(qᵗₙ[1,1,kz]), maximum(qᵗₙ[1,1,kz]))
qˡ_limits = z_limited_extrema(qˡₙ, kz)
wvar_limits = z_limited_extrema(wvarₙ, kz)

# symmetric limits for w
wmax = maximum(abs, parent(wₙ)[:, 1, kz, 1])
w_limits = (-wmax, wmax)

#hmw  = heatmap!(axw,  x, z, wₙ;  colormap=:balance, colorrange=w_limits)
hmw  = heatmap!(axw,  x, z, wvarₙ;  colormap=:solar, colorrange=wvar_limits)
hmu  = heatmap!(axu,  x, z, uₙ;  colormap=Reverse(:speed), colorrange=u_limits)
hmv  = heatmap!(axv,  x, z, vₙ;  colormap=Reverse(:speed), colorrange=v_limits)
hmθ  = heatmap!(axθ,  x, z, θₙ;  colormap=:thermal, colorrange=θ_limits)
hmq  = heatmap!(axq,  x, z, qᵗₙ; colormap=Reverse(:Purples_4), colorrange=qᵗ_limits)
hmql = heatmap!(axql, x, z, qˡₙ; colormap=:deep, colorrange=qˡ_limits)


Colorbar(fig[1, 2], hmw)
Colorbar(fig[1, 4], hmu)
Colorbar(fig[1, 6], hmv)
Colorbar(fig[2, 2], hmθ)
Colorbar(fig[2, 4], hmq)
Colorbar(fig[2, 6], hmql)


for ax in (axw, axu, axv, axθ, axq, axql)
    ylims!(ax, 0, zmax)
end
#fig

save("plane-mean.png", fig)
