# --------------------------------------------------------------------------------------------------
# Internal functions to read parametric lens model
# --------------------------------------------------------------------------------------------------
# Lens models that do not require position
const NO_POSITION = Set([:ExternalEffects, :ExternalEffects3, :Multipole])


# Lens models that require ADD to lens
const REQUIRE_ADD = Set([:PointLens, :PlummerLens, :GaussianLens, :SersicLens, :HernquistLens, 
   :aHernquistLens, :eHernquistMDLens, :NFWLens, :tNFWLens, :gNFWLens, :aNFWLens, :eNFWMDLens, 
   :EinastoLens, :MultiPlummerLens, :MultiGaussianLens])


# Lens models that require scaling
const REQUIRE_SCALING = Set([:MultiPJELens])


# Names of all scaling-relation parameters, in the order expected by ScalingRelation
const SCALING_PARAMS = (:ref_mag, :ref_sigma, :ref_core, :ref_cut, :slope_sigma, :slope_core, :slope_cut)


# Read a galaxy catalog file into a GalaxyComponent
function _read_galaxy_catalog(file_name::String, observation::Observation)
   catalog_data = readdlm(file_name, comments=true, comment_char='#')
   catalog_data = Float64.(catalog_data)

   # Check if a valid (RA, Dec) is provided as reference or (0, 0) is used
   x_lens, y_lens = _to_arcsec(observation, catalog_data[:, 2], catalog_data[:, 3])

   return GalaxyComponent(
      n       = size(catalog_data, 1),
      x_c     = x_lens,
      y_c     = y_lens,
      obs_mag = catalog_data[:, 4],
      eps     = catalog_data[:, 5],
      pa      = catalog_data[:, 6])
end


# Read one scaling-relation dictionary; parameters are pushed with the given owner so that
# each galaxy-cluster plane can have its own (possibly free) scaling relations
function _read_scaling_relation!(scaling_dict::Dict, owner::Symbol, params::Vector{Parameter})
   for param in SCALING_PARAMS
      r, l, u = _extract_param_range(scaling_dict[param])
      push!(params, Parameter(owner=owner, name=param, refer=r, lower=l, upper=u))
   end
   return nothing
end

function _parametric(lens_dict::Dict, params::Vector{Parameter}, observation::Observation)
   # Construct a composite lens using initial values
   n_lenses = lens_dict[:total_lenses]

   # Initialize lens name vector
   lens_name = Vector{LensComponent}(undef, n_lenses)

   for i in 1:n_lenses
      lens_id = Symbol(:lens, i)
      indi_lens_dict = lens_dict[lens_id]
      name = Symbol(indi_lens_dict[:lens])

      # Store lens model name in lens_name vector
      lens_name[i] = LensComponent(owner=lens_id, name=name)

      # Add distance parameters
      if name ∈ REQUIRE_ADD
         push!(params, Parameter(owner = lens_id, name  = :D_d, refer = observation.D_d, lower = observation.D_d, upper = observation.D_d))
      end

      if name ∉ REQUIRE_SCALING
         # --- Lens position parameters (always provided) ---
         if name ∉ NO_POSITION
            rx, lx, ux = _extract_param_range(indi_lens_dict[:x_c])
            ry, ly, uy = _extract_param_range(indi_lens_dict[:y_c])

            # Check if a valid (RA, Dec) is provided as reference or (0, 0) is used.
            # If reference = (0., 0.)  ⇒ Lens positions are provided in arcseconds
            # If reference = (RA, Dec) ⇒ Lens positions are provided in RA and Dec. Conversion needed.
            x_lens, y_lens = _to_arcsec(observation, rx, ry)

            # Add lens position parameters to the parameter vector
            push!(params, Parameter(owner=lens_id, name=:x_c, refer=x_lens, lower=lx, upper=ux))
            push!(params, Parameter(owner=lens_id, name=:y_c, refer=y_lens, lower=ly, upper=uy))
         end
         
         # --- Remaining lens parameters ---
         for (k, v) in indi_lens_dict
            if k ∈ (:lens, :x_c, :y_c)
               continue
            end

            # Extract parameter values and bounds
            r, l, u = _extract_param_range(v)

            # Add parameter to the reference, lower, and upper vectors
            push!(params, Parameter(owner=lens_id, name=k, refer=r, lower=l, upper=u))
         end
      else
         # --- Galaxy-cluster member lens (scaling relations) ---
         # Galaxy catalog for THIS lens
         _require(indi_lens_dict, :galaxy_file)
         galaxies[lens_id] = _read_galaxy_catalog(indi_lens_dict[:galaxy_file], observation)

         # Scaling relations for THIS lens: nested subsection inside the lens block, with
         # parameter owner :scaling<i> (the scaling relation is part of the lens)
         _require(indi_lens_dict, :scaling_relation)
         _read_scaling_relation!(indi_lens_dict[:scaling_relation], Symbol(:scaling, i), params)
      end
   end
end