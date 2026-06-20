module LensModelUtils

# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using Random
using StatsBase


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..Constants
using ..Cosmology
using ..Lenses
using ..LFUtils
using ..LensModelIO
using ..Likelihood

# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export free_parameters
export θ_initializer
export param_dict
export adis_current
export build_lens
export lens_quantities


# --------------------------------------------------------------------------------------------------
# Parameter-space → Physical-space transformation
# --------------------------------------------------------------------------------------------------
const PARAM_TRANSFORM = Dict{Symbol,Function}(   
   # Mass:Log10[M / M☉] -> [M / M☉]
   :mass => x -> 10^x,
)

function transform_params!(pvals::Dict{Tuple{Symbol,Symbol}, <:Real})
   for key in collect(keys(pvals))
      name = key[2]
      if haskey(PARAM_TRANSFORM, name)
         pvals[key] = PARAM_TRANSFORM[name](pvals[key])
      end
   end
end


# --------------------------------------------------------------------------------------------------
# Get view on free parameters
# --------------------------------------------------------------------------------------------------
function free_parameters(model::ModelConfig)
   return @view model.parameters[model.free_param_idxs]
end


# --------------------------------------------------------------------------------------------------
# Get parameter values
# --------------------------------------------------------------------------------------------------
function θ_reference(model::ModelConfig)
   getfield.(free_parameters(model), :refer)
end

function θ_random(model::ModelConfig)
   θ = Vector{Float64}(undef, length(model.free_param_idxs))
   for (i, p) in enumerate(free_parameters(model))
      θ[i] = p.lower + rand() * (p.upper - p.lower)
   end
   return θ
end

function θ_grid(model::ModelConfig, max_runs::Int64)
   fps = free_parameters(model)
   nd = length(fps)
   
   # Automatically determine points per dimension
   n_per_dim = ceil(Int, max_runs^(1/nd))
   total = n_per_dim^nd

   # Generate grid
   θ_grid = Vector{Vector{Float64}}(undef, total)
   for idx in 1:total
      θ = Vector{Float64}(undef, nd)
      i = idx - 1
      for (k, p) in enumerate(fps)
         ik = i % n_per_dim
         i ÷= n_per_dim
         θ[k] = p.lower + (ik / (n_per_dim - 1)) * (p.upper - p.lower)
      end
      θ_grid[idx] = θ
   end
   return θ_grid
end

function θ_lhs(model::ModelConfig, max_runs::Int64)
   fps = free_parameters(model)
   nd = length(fps)

   lhs_matrix = zeros(max_runs, nd)
   for k in 1:nd
      # Divide [0, 1] into n intervals and pick a random point in each
      intervals = ((0:max_runs-1) .+ rand(max_runs)) ./ max_runs
      
      # shuffle for randomness
      shuffle!(intervals)
      
      # Scale to parameter bounds
      lhs_matrix[:, k] = fps[k].lower .+ intervals .* (fps[k].upper - fps[k].lower)
   end
   # Convert each row to a vector
   return [vec(lhs_matrix[i, :]) for i in 1:max_runs]
end

function θ_initializer(model::ModelConfig; run_mode::Symbol=:random, max_runs::Int64=1)
   if run_mode === :reference
      return [θ_reference(model) for _ in 1:max_runs]
   elseif run_mode === :jitter
      return [θ_reference(model) .+ 1e-3 .* randn(length(model.free_param_idxs)) for _ in 1:max_runs]
   elseif run_mode === :random
      return [θ_random(model) for _ in 1:max_runs]
   elseif run_mode === :grid
      return θ_grid(model, max_runs)
   elseif run_mode === :lhs
      return θ_lhs(model, max_runs)
   else
      error("Unknown run_mode: $run_mode")
   end

end


# --------------------------------------------------------------------------------------------------
# Parameter dictionary (θ -> pvals)
# --------------------------------------------------------------------------------------------------
function param_dict(model::ModelConfig, θ::Vector{Float64}, param_ref::Dict{Tuple{Symbol,Symbol}, <:Real})
   # Make a copy of the reference dictionary
   pvals = copy(param_ref)

   # Update free parameters
   @inbounds for (i, idx) in enumerate(model.free_param_idxs)
      p = model.parameters[idx]
      pvals[p.key] = θ[i]
   end
   return pvals
end


# --------------------------------------------------------------------------------------------------
# Angular-diameter distance (pvals -> adis)
# --------------------------------------------------------------------------------------------------
function adis_current(model::ModelConfig, pvals::Dict{Tuple{Symbol,Symbol}, <:Real})
   nsrc = length(model.source_config.sources)
   adis = Vector{Float64}(undef, nsrc)

   @inbounds for i in 1:nsrc
      key = (Symbol(:source, i), Symbol(:adis, i))
      adis[i] = pvals[key]
   end
   return adis
end


# --------------------------------------------------------------------------------------------------
# Build lens model from physical paramerters
# --------------------------------------------------------------------------------------------------
REQUIRE_COSMO = Set([:NFWLens, :eNFWMDLens, :aNFWLens, :tNFWLens, :gNFWLens, :EinastoLens])
REQUIRE_SCALING = Set([:MultiPJELens])
function build_lens(model::ModelConfig, pvals::Dict{Tuple{Symbol,Symbol}, <:Real})
   # Update scaling if it has free parameters
   scaling_params = Dict{Symbol, Any}() 
   for (k, v) in pvals
      if k[1] == :scaling
         scaling_params[k[2]] = v
      end
   end
   
   if !isempty(scaling_params)
      updated_scaling = ScalingRelation(; scaling_params...)
   else
      updated_scaling = model.lens_config.scaling
   end
   
   # Determine the number of components from the lens model container
   n_lens = length(model.lens_config.components)

   # Transform parameters
   transform_params!(pvals)

   # Initialize an empty vector to store lens parameters
   lens_vector = NamedTuple[]

   # Iterate over each lens component
   components = model.lens_config.components

   for i in 1:n_lens
      lens_id = Symbol(:lens, i)
      lens_params = Dict{Symbol, Union{Symbol, Int64, Float64, Vector{Float64}, Cosmology.AbstractCosmology}}()
      
      # Get lens name and add it to the lens_params dictionary
      for (k, v) in enumerate(components)
         if v.owner == lens_id
            lens_params[:lens] = v.name
         end
      end

      # Get other lens parameters and add it to the lens_params dictionary
      for (k, v) in pvals
         if k[1] == lens_id
            lens_params[k[2]] = v
         end
      end

      # Check if this lens model requires cosmology and lens redshift
      if lens_params[:lens] ∈ REQUIRE_COSMO
         lens_params[:cosmology] = model.cosmology
         lens_params[:z_d]       = model.observation.z_d
      end

      # Check if this lens model requires galaxy component
      if lens_params[:lens] ∈ REQUIRE_SCALING
         # Lens parameters from file
         lens_params[:x_c] = model.lens_config.galaxies.x_c
         lens_params[:y_c] = model.lens_config.galaxies.y_c
         lens_params[:eps] = model.lens_config.galaxies.eps
         lens_params[:pa]  = model.lens_config.galaxies.pa

         # Observed magnitudes (vector)
         obs_mag = model.lens_config.galaxies.obs_mag

         # Reference magnitude
         ref_mag = updated_scaling.ref_mag

         # L/L⋆ (vector)
         l_lstar = @. 10.0^(-0.4 * (obs_mag - ref_mag))

         # Velocity dispersion (vector)
         lens_params[:v_d] = @. updated_scaling.ref_sigma * l_lstar^updated_scaling.slope_sigma
         
         # Core radius (vector)
         lens_params[:x_s] = @. updated_scaling.ref_core * l_lstar^updated_scaling.slope_core

         # Truncation radius (vector)
         lens_params[:x_t] = @. updated_scaling.ref_cut * l_lstar^updated_scaling.slope_cut
      end
      push!(lens_vector, (; lens_params...))
   end
   return Lenses.init_CompositeLens(lens_vector)
end


# --------------------------------------------------------------------------------------------------
# All lensing quantities needed for log-likelihood 
# --------------------------------------------------------------------------------------------------
function lens_quantities(model::ModelConfig, lens::Lenses.AbstractLens)
   # Count the total number of knots in the lens model
   n_knots = sum(length(s.knots) for s in model.source_config.sources)

   # Allocate outputs (vector of vectors)
   ψ_all  = Vector{Vector{Float64}}(undef, n_knots)
   αx_all = Vector{Vector{Float64}}(undef, n_knots)
   αy_all = Vector{Vector{Float64}}(undef, n_knots)
   A_all  = Vector{NTuple{4, Vector{Float64}}}(undef, n_knots)

   kid = 1
   for src in model.source_config.sources
      for knot in src.knots
         # One image system knot positions
         x = knot.x
         y = knot.y

         # Potential
         ψ_all[kid] = Lenses.get_potential(lens, x, y)

         # Deflection
         αx_all[kid], αy_all[kid] = Lenses.get_deflection(lens, x, y)

         # Deformation tensor
         ψxx, ψyy, ψxy = Lenses.get_jacobian(lens, x, y)
         A_all[kid] = (ψxx, ψxy, copy(ψxy), ψyy)

         # Increment
         kid = kid + 1
      end
   end
   return ψ_all, αx_all, αy_all, A_all
end

end