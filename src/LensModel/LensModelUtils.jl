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
export current_cosmology
export adis_current
export build_lens
export lens_quantities
export lens_quantities_fd


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
# Parameter dictionary (θ -> pvals)
# --------------------------------------------------------------------------------------------------
function current_cosmology(model::ModelConfig, pvals::Dict{Tuple{Symbol,Symbol}, <:Real})
   # Cosmology fixed -> return the cached reference object (no rebuild, no allocation)
   if model.static_ADD
      return model.cosmology
   end

   # Cosmology free -> rebuild from the sampled parameters.
   kw = Dict{Symbol,Float64}()
   for (k, v) in pvals
      if k[1] === :cosmology
         kw[k[2]] = Float64(v)
      end
   end
   return Cosmology.init_cosmology(; kw...)
end


# --------------------------------------------------------------------------------------------------
# Angular-diameter distance (pvals -> adis)
# --------------------------------------------------------------------------------------------------
function adis_current(model::ModelConfig, pvals::Dict{Tuple{Symbol,Symbol}, <:Real}, cosmo::Cosmology.AbstractCosmology)
   nsrc = length(model.source_config.sources)
   adis = Vector{Float64}(undef, nsrc)

   if model.sample_z == :adis
      @inbounds for i in 1:nsrc
         adis[i] = pvals[(Symbol(:source, i), Symbol(:adis, i))]
      end
   else
      z_d   = model.observation.z_d
      if model.static_ADD
         d_C = model.observation.d_C
      else
         d_C = Cosmology.comoving_distance_radial(cosmo, 0.0, z_d)
      end

      @inbounds for i in 1:nsrc
         zs      = Float64(pvals[(Symbol(:source, i), Symbol(:zs, i))])
         adis[i] = Cosmology.zs2adis(cosmo, z_d, zs; dC = d_C)
      end
   end
   return adis
end


# --------------------------------------------------------------------------------------------------
# Build lens model from physical paramerters
# --------------------------------------------------------------------------------------------------
# Lens models that need cosmology + lens redshift passed through to the constructor.
const REQUIRE_COSMO = Set([:NFWLens, :eNFWMDLens, :aNFWLens, :tNFWLens, :gNFWLens, :EinastoLens])

# Lens models whose parameters are generated from a galaxy catalog + scaling relations.
const REQUIRE_SCALING = Set([:MultiPJELens])

# Names of all scaling-relation parameters. Must match LensModelIO.SCALING_PARAMS and the
# field names of LensModelIO.ScalingRelation (construction is by keyword, so order is irrelevant).
const SCALING_PARAMS = (:ref_mag, :ref_sigma, :ref_core, :ref_cut, :slope_sigma, :slope_core, :slope_cut)

# Rebuild the ScalingRelation for lens component i from pvals.
function _scaling_from_pvals(pvals::Dict{Tuple{Symbol, Symbol}, <:Real}, i::Int64)
   owner = Symbol(:scaling, i)
   vals  = Dict{Symbol,Float64}()
   for name in SCALING_PARAMS
      key = (owner, name)
      haskey(pvals, key) || error("Missing scaling parameter $(key) for lens-$(i) while building lens.")
      vals[name] = Float64(pvals[key])
   end
   return ScalingRelation(; vals...)
end

function build_lens(model::ModelConfig, pvals::Dict{Tuple{Symbol,Symbol}, <:Real})
   # Determine the number of components from the lens model container
   n_lens     = length(model.lens_config.components)
   components = model.lens_config.components
   galaxies   = model.lens_config.galaxies   # Dict{Symbol, GalaxyComponent}, keyed by lens id

   # Lens ADD consistent with the current cosmology
   z_d     = model.observation.z_d
   if model.static_ADD
      D_d_dyn = model.observation.D_d
   else
      D_d_dyn = Cosmology.angular_diameter_distance(cosmo, 0.0, z_d)
   end

   # Transform parameters
   transform_params!(pvals)
   
   # Initialize an empty vector to store lens parameters
   lens_vector = NamedTuple[]

   for i in 1:n_lens
      lens_id = Symbol(:lens, i)
      lens_params = Dict{Symbol, Union{Symbol, Int64, Float64, Vector{Float64}, Cosmology.AbstractCosmology}}()
      
      # Lens name. components[i] corresponds to lens i (built in order in LensModelIO).
      name = components[i].name
      lens_params[:lens] = name

      # Collect the parameters owned by this lens component
      for (k, v) in pvals
         if k[1] == lens_id
            lens_params[k[2]] = v
         end
      end

      # Keep a cosmology-dependent lens distance in sync with a free cosmology
      if !model.static_ADD && haskey(lens_params, :D_d)
         lens_params[:D_d] = D_d_dyn
      end

      # Cosmology-dependent profiles need the cosmology object and the lens redshift
      if lens_params[:lens] ∈ REQUIRE_COSMO
         lens_params[:cosmology] = model.cosmology
         lens_params[:z_d]       = z_d
      end

      # Galaxy-cluster member lens: expand the scaling relation over the galaxy catalog
      if lens_params[:lens] ∈ REQUIRE_SCALING
         galaxy = galaxies[lens_id]

         # Rebuild the (possibly updated) scaling relation for THIS lens from pvals
         relation = _scaling_from_pvals(pvals, i)

         # Per-galaxy geometric parameters straight from the catalog
         lens_params[:x_c] = galaxy.x_c
         lens_params[:y_c] = galaxy.y_c
         lens_params[:eps] = galaxy.eps
         lens_params[:pa]  = galaxy.pa

         # Luminosity ratio L/L⋆ from observed magnitudes (vector)
         obs_mag = galaxy.obs_mag
         l_lstar = @. 10.0^(-0.4 * (obs_mag - relation.ref_mag))

         # Scaling relations: velocity dispersion, core radius, truncation radius (vectors)
         lens_params[:v_d] = @. relation.ref_sigma * l_lstar^relation.slope_sigma
         lens_params[:x_s] = @. relation.ref_core  * l_lstar^relation.slope_core
         lens_params[:x_t] = @. relation.ref_cut   * l_lstar^relation.slope_cut
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


# --------------------------------------------------------------------------------------------------
# Deformation tensors at finite-difference shifted positions (needed for the flux chi2)
# Computed only for knots that carry flux measurements; other entries hold empty vectors.
# --------------------------------------------------------------------------------------------------
function lens_quantities_fd(model::ModelConfig, lens::Lenses.AbstractLens)
   # Count the total number of knots in the lens model
   n_knots = sum(length(s.knots) for s in model.source_config.sources)

   # Finite-difference step (arcsec)
   h_fd = model.source_config.pixel_fd

   # Allocate outputs
   A_xp_all = Vector{NTuple{4, Vector{Float64}}}(undef, n_knots)
   A_xm_all = Vector{NTuple{4, Vector{Float64}}}(undef, n_knots)
   A_yp_all = Vector{NTuple{4, Vector{Float64}}}(undef, n_knots)
   A_ym_all = Vector{NTuple{4, Vector{Float64}}}(undef, n_knots)

   # Placeholder for knots without flux measurements
   empty_A = (Float64[], Float64[], Float64[], Float64[])

   kid = 1
   for src in model.source_config.sources
      for knot in src.knots
         # Skip knots without flux measurements (kid must still advance)
         if isempty(knot.m)
            A_xp_all[kid] = empty_A
            A_xm_all[kid] = empty_A
            A_yp_all[kid] = empty_A
            A_ym_all[kid] = empty_A
            kid = kid + 1
            continue
         end

         # One image system knot positions
         x = knot.x
         y = knot.y

         # Deformation tensors at the four shifted positions
         ψxx, ψyy, ψxy = Lenses.get_jacobian(lens, x .+ h_fd, y)
         A_xp_all[kid] = (ψxx, ψxy, copy(ψxy), ψyy)

         ψxx, ψyy, ψxy = Lenses.get_jacobian(lens, x .- h_fd, y)
         A_xm_all[kid] = (ψxx, ψxy, copy(ψxy), ψyy)

         ψxx, ψyy, ψxy = Lenses.get_jacobian(lens, x, y .+ h_fd)
         A_yp_all[kid] = (ψxx, ψxy, copy(ψxy), ψyy)

         ψxx, ψyy, ψxy = Lenses.get_jacobian(lens, x, y .- h_fd)
         A_ym_all[kid] = (ψxx, ψxy, copy(ψxy), ψyy)

         # Increment
         kid = kid + 1
      end
   end
   return A_xp_all, A_xm_all, A_yp_all, A_ym_all
end

end