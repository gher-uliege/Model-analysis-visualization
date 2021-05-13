# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.6.0
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# # Tutorial on model output analysis and visualization

# %% [markdown]
# ## Installation and setup


# %% [markdown]
# ## Download some example data

# %% [markdown]
# ## Vizualize surface data
# Setup PyCall/PyPlot to handle missing data

# %%
using PyCall
using PyCall: PyObject
function PyObject(a::Array{Union{T,Missing},N}) where {T,N}
  numpy_ma = pyimport("numpy")["ma"]
  pycall(numpy_ma["array"], Any, coalesce.(a,zero(T)), mask=ismissing.(a))
end

# %%
using PyPlot
using NCDatasets

ds = NCDataset("roms_his.nc");
pcolormesh(ds["zeta"][:,:,10]');
colorbar()

# %%
pcolormesh(ds["temp"][:,:,end,10]');


# %%
using Statistics
lon_rho = ds["lon_rho"][:,:]
lat_rho = ds["lat_rho"][:,:]
temp = ds["temp"][:,:,:,10]


pcolormesh(lon_rho,lat_rho,temp[:,:,end]);
gca().set_aspect(1/cosd(mean(ylim()))); 
colorbar(orientation = "horizontal")

# %%
#using Plots
#Plots.pyplot()

# %%
#heatmap(lon_rho,lat_rho,temp[:,:,end])

# %%
using ROMS


Vtransform = ds["Vtransform"][]
Vstretching = ds["Vstretching"][]
theta_s = ds["theta_s"][]
theta_b = ds["theta_b"][]
hc = ds["hc"][]
N = ds.dim["s_rho"]
h = ds["h"][:,:]
zeta = ds["zeta"][:,:,10]



@show Vtransform, Vstretching, theta_s, theta_b, hc, N;

igrid = 1

# %%
z_r = ROMS.set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, igrid, h; zeta = zeta);
z_r2 = ROMS.set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, igrid, h);

# %%
diffz = abs.(z_r - z_r2)
maxdiffz,ind = findmax(coalesce.(diffz,0))

# %%
maximum(skipmissing(abs.(z_r - z_r2)))

# %%
using Interpolations
using Test
temp = nomissing(ds["temp"][:,:,:,10],NaN);
i = 1; j = 1

itp = interpolate((z_r2[i,j,:],), temp[i,j,:], Gridded(Linear()))

k = 40

itp(-50)

# %%
@test itp(z_r2[i,j,k]) ≈ temp[i,j,k]

# %% [markdown]
# function vinterp(z,v,zi)
#    vi = similar(zi)
#    for j = 1:size(z,2)
#      for i = 1:size(z,1)
#        itp = extrapolate(interpolate((z[i,j,:],), v[i,j,:], Gridded(Linear())),NaN)
#        
#             
#         
# end

# %%

#Vtransform, Vstretching, theta_s, theta_b, hc, N, igrid, h

# %% [markdown]
# ## Vizualize horizontal sections data

# %% [markdown]
# ## Vizualize vertical sections data

# %% [markdown]
# ## Compute integrals and averages

# %% [markdown]
# ## Validate the model with in situ observations and satellite observations

# %% [markdown]
# ## Compute Empirical Orthogonal Functions (EOFs)

# %% [markdown]
# ## Compute derived quantities (Brunt-Väisälä frequency, potential temperature, absolute salinity)

# %% [markdown]
# ## Drifter trajectories
