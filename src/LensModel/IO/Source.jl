# --------------------------------------------------------------------------------------------------
# ---------------- Resolve the z-sampling mode -----------------------------------------------------
# --------------------------------------------------------------------------------------------------
# The source-distance parameter can be either
#    :adis - the distance ratio D_ds/D_os (uniform prior in adis; fast, no distance
#            computation per evaluation). Only meaningful in single-plane mode with a
#            fixed cosmology, where adis is the only channel through which z_s enters.
#    :zs   - the source redshift itself (uniform prior in z_s; derived distances must
#            be recomputed per evaluation whenever z_s or the cosmology is free).
# The user chooses via the optional ** sample_z ** key in the source section (default: adis).
# Forcing rules: multi-plane mode and free cosmology both require :zs.
function _resolve_sample_z(dict::Dict, multiplane::Bool, free_cosmo::Bool)
   source_dict = dict[:source]

   # Detect free source redshifts (covers both direct and from-file YAML overrides;
   # redshifts provided inside a source file are always fixed)
   free_zs = false
   for (k, v) in source_dict
      if v isa Dict && haskey(v, :z_s)
         r, l, u = _extract_param_range(v[:z_s])
         if l != u
            free_zs = true
            break
         end
      end
   end

   # User choice (only meaningful when z_s is free)
   user_key = haskey(source_dict, :sample_z) ? Symbol(source_dict[:sample_z]) : nothing
   if user_key !== nothing && user_key ∉ (:adis, :zs)
      error("Unknown ** sample_z ** value: $user_key. Allowed values are ** adis ** and ** zs **.")
   end

   if !free_zs
      # Nothing is sampled - fixed constants are stored in the mode's natural variable
      sample_z = multiplane ? :zs : :adis
   else
      sample_z = user_key === nothing ? :adis : user_key

      # Forcing rules
      if multiplane && sample_z == :adis
         user_key == :adis && @warn "sample_z = adis is not available in multi-plane mode - using zs."
         sample_z = :zs
      end
      if free_cosmo && sample_z == :adis
         user_key == :adis && @warn "sample_z = adis is not available with free cosmology - using zs."
         sample_z = :zs
      end
   end
   return sample_z, free_zs
end


# --------------------------------------------------------------------------------------------------
# ---------------- Read Source Model ---------------------------------------------------------------
# --------------------------------------------------------------------------------------------------  
# Convert a YAML list entry into a Vector{Float64} / Vector{Int64}. Scalars are deliberately
# rejected: measurements and errors must be fully specified, one entry per image, even if identical.
@inline _as_float_vec(x::AbstractVector) = Float64.(x)
@inline _as_float_vec(x::Any) = error("Expected a list with one entry per image, got: ** $x **")

@inline _as_int_vec(x::AbstractVector) = Int64.(x)
@inline _as_int_vec(x::Any) = error("Expected a list with one entry per image, got: ** $x **")

# Read an optional per-knot measurement (value + error), e.g. flux (magnitudes) or time delays.
# Both quantities are RELATIVE measurements (the source magnitude / delay zero-point is profiled
# out in the likelihood), therefore:
#   - One entry per image is required, aligned with the image order of (x, y).
#   - Images without a measurement are marked with error <= 0 (their value entry is ignored).
#   - At least two measured images are required; a single measurement carries no information
#     and the knot is treated as having no measurement (with a warning).
# Returns a pair of empty vectors when the measurement is absent or uninformative.
function _knot_measurement(knot_dict::Dict, val_key::Symbol, err_key::Symbol, n_img::Int64, i::Int64, k::Int64)
   # Measurement not provided --> return empty vectors
   if !haskey(knot_dict, val_key)
      return Float64[], Float64[]
   end

   # Error entry is required whenever the value entry is present
   if !haskey(knot_dict, err_key)
      error("Missing ** $err_key ** for ** $val_key ** in source-$i, knot-$k")
   end

   v = _as_float_vec(knot_dict[val_key])
   σ = _as_float_vec(knot_dict[err_key])

   # One entry per image, always
   if length(v) != n_img || length(σ) != n_img
      error("** $val_key / $err_key ** must have one entry per image ($n_img) in source-$i, knot-$k. " *
            "Mark images without a measurement using $err_key <= 0.")
   end

   # Count measured images (error > 0)
   n_meas = count(>(0.0), σ)
   if n_meas == 0
      return Float64[], Float64[]
   elseif n_meas == 1
      @warn "Only one measured image for ** $val_key ** in source-$i, knot-$k. Relative " *
            "measurements need at least two measured images - this knot will be ignored."
      return Float64[], Float64[]
   end
   return v, σ
end

function _source_direct!(dict::Dict, cosmo::Cosmology.AbstractCosmology, observation::Observation, params::Vector{Parameter}, sample_z::Symbol)
   # Get source dictionary from the full dictionary
   source_dict = dict[:source]

   # Make sure that total number of sources is present and greater than zero
   _require(source_dict, :total_sources)
   if source_dict[:total_sources] <= 0
      error("Total number of sources must be greater than zero.")
   end

   # Get number of sources
   n_source = source_dict[:total_sources]

   # Get reference lens redshift
   z_d = dict[:observation][:z_d]

   # Track whether any knot carries flux / time-delay measurements
   any_parity = false
   any_flux   = false
   any_td     = false

   # Run over sources
   sources = Vector{Source}(undef, n_source)
   for i in 1:n_source
      source_id = Symbol(:source, i)
      indi_source_dict = source_dict[source_id]

      # Extract parameter values and bounds
      r, l, u = _extract_param_range(indi_source_dict[:z_s])

      if sample_z == :zs
         # Source-distance parameter is the redshift itself
         push!(params, Parameter(owner = source_id, 
                                 name  = Symbol(:zs, i), 
                                 refer = r, 
                                 lower = l, 
                                 upper = u))
      else
         # Source-distance parameter is the distance ratio: convert redshift to adis
         # (uniform prior in adis - see _resolve_sample_z)
         adis_r, adis_l, adis_u = _zs_to_adis(cosmo, z_d, r, l, u)

         # Add adis parameter to the reference, lower, and upper vectors
         push!(params, Parameter(owner = source_id, 
                                 name  = Symbol(:adis, i), 
                                 refer = adis_r, 
                                 lower = adis_l, 
                                 upper = adis_u))
      end

      # Run over knots
      n_knot = indi_source_dict[:total_knots]
      knots = Vector{Knot}(undef, n_knot)
      for k in 1:n_knot
         knot_id = Symbol(:knot, k)
         # Read knot image position(s)
         x = indi_source_dict[knot_id][:x]
         y = indi_source_dict[knot_id][:y]

         # Assert that the number of values is the same for (x, y, σx, σy, θ)
         if length(x) != length(y)
            error("Inconsistent knot position dimensions in source-$i, knot-$k")
         end

         # Read knot image error(s)
         if haskey(indi_source_dict[knot_id], :sigma)
            σx = indi_source_dict[knot_id][:sigma]
            σy = indi_source_dict[knot_id][:sigma]
            σθ = 0.0 .* indi_source_dict[knot_id][:sigma]
         elseif haskey(indi_source_dict[knot_id], :sigma_x) &&
                haskey(indi_source_dict[knot_id], :sigma_y) &&
                haskey(indi_source_dict[knot_id], :sigma_theta)
            σx = indi_source_dict[knot_id][:sigma_x]
            σy = indi_source_dict[knot_id][:sigma_y]
            σθ = indi_source_dict[knot_id][:sigma_theta]
         else
            error("Invalid knot error parameters in source-$i, knot-$k")
         end
         # Assert that the number of values is the same for (x, y, σx, σy, σθ)
         if length(x) != length(σx) || length(x) != length(σy) || length(x) != length(σθ)
            error("Inconsistent knot position and error dimensions in source-$i, knot-$k")
         end

         # Read optional parity values: one entry per image, +1 / -1, with 0 marking images
         # without parity information. Parity is an absolute per-image quantity, so a single
         # measured image is already informative (no minimum count).
         p = zeros(Int64, length(x))
         if haskey(indi_source_dict[knot_id], :parity)
            p_temp = indi_source_dict[knot_id][:parity]
            if length(p_temp) != length(x)
               error("** parity ** must have one entry per image ($(length(x))) in source-$i, knot-$k. " *
                     "Mark images without parity information with 0.")
            end
            if !all(v -> v in (-1, 0, 1), p_temp)
               error("** parity ** values must be +1, -1, or 0 (no information) in source-$i, knot-$k")
            end
            p = p_temp
         end
         if any(!=(0), p)
            any_parity = true
         end

         # Read optional flux measurements (magnitudes) and their errors
         m, σm = _knot_measurement(indi_source_dict[knot_id], :mag, :sigma_mag, length(x), i, k)
         if !isempty(m)
            any_flux = true
         end

         # Read optional time-delay measurements (days, arbitrary zero-point) and their errors
         td, σ_td = _knot_measurement(indi_source_dict[knot_id], :td, :sigma_td, length(x), i, k)
         if !isempty(td)
            any_td = true
         end

         # Convert to arcsec if reference is not (0, 0)
         x_arcsec, y_arcsec = _to_arcsec(observation, x, y)
         knots[k] = Knot(x      = x_arcsec, 
                         σx     = σx, 
                         y      = y_arcsec, 
                         σy     = σy, 
                         σθ     = deg2rad.(σθ), 
                         parity = p,
                         m      = m, 
                         σm     = σm, 
                         td     = td, 
                         σ_td   = σ_td)
      end
      sources[i] = Source(knots=knots)
   end
   # Parity / flux / time-delay terms are OPT-IN: used only if the corresponding use_* flag is
   # explicitly set to true in the source section AND the data is present.
   # Parity
   if get(source_dict, :use_parity, false)
      use_parity = any_parity
   else
      use_parity = false
   end

   # Flux
   if get(source_dict, :use_flux, false)
      use_flux = any_flux
   else
      use_flux = false
   end

   # Time Delay
   if get(source_dict, :use_time_delay, false)
      use_td = any_td
   else
      use_td = false
   end

   if get(source_dict, :use_parity, false) === true && !any_parity
      @warn "use_parity is set but no knot provides ** parity ** values - parity term disabled."
   end

   if get(source_dict, :use_flux, false) === true && !any_flux
      @warn "use_flux is set but no knot provides ** mag / sigma_mag ** - flux chi2 disabled."
   end

   if get(source_dict, :use_time_delay, false) === true && !any_td
      @warn "use_time_delay is set but no knot provides ** td / sigma_td ** - time-delay chi2 disabled."
   end

   # Parity penalty strength (default defined here, override with parity_force in the source section)
   parity_force = Float64(get(source_dict, :parity_force, 1000.0))

   # Finite-difference step (arcsec) used in the flux chi2 magnification correction
   pixel_fd = Float64(get(source_dict, :pixel_fd, observation.pixel_scale))

   return SourceConfig(sources        = sources, 
                       use_parity     = use_parity, 
                       parity_force   = parity_force,
                       use_flux       = use_flux, 
                       use_time_delay = use_td, 
                       pixel_fd       = pixel_fd)
end


function _source_from_file!(dict::Dict, cosmo::Cosmology.AbstractCosmology, observation::Observation, params::Vector{Parameter}, sample_z::Symbol)
   # Get source dictionary from the full dictionary
   source_dict = dict[:source]

   # Get file name
   file_name = source_dict[:source_file]
   file_data = readdlm(file_name, comments=true, comment_char='#')
   file_data = Float64.(file_data)

   # --- Parameters directly read from the YAML file -----------------------------------------------
   # Get total number of sources from first column in file
   n_source = Int64(maximum(file_data[:, 1]))

   # Get reference lens redshift
   z_d = dict[:observation][:z_d]

   # Number of columns in the source file
   # Column 1 : source id, 
   # Column 2 : knot id
   # Column 3 : x
   # Column 4 : y
   # Column 5 : z_s
   # Column 6 : sigma_x
   # Column 7 : sigma_y
   # Column 8 : sigma_theta
   # Column 9 : (optional) parity (parity: +1 / -1, with 0 marking images without parity information)
   # Column 10: (optional) magnitude
   # Column 11: (optional) sigma_mag (sigma_mag <= 0 --> no flux data for knot)
   # Column 12: (optional) td
   # Column 13: (optional) sigma_td (sigma_td  <= 0 --> no time-delay data for knot)
   n_col = size(file_data, 2)

   # Track whether any knot carries parity / flux / time-delay measurements
   any_parity = false
   any_flux   = false
   any_td     = false

   # --- Parameters read from the source file ------------------------------------------------------
   sources = Vector{Source}(undef, n_source)
   for i in 1:n_source
      # Generate source ID
      source_id = Symbol(:source, i)

      # Get individual source data from file
      mask = file_data[:, 1] .== i
      source_data = file_data[mask, :]

      # Get source redshift
      if haskey(source_dict, source_id)
         # Get individual source dict from YAML file
         indi_source_dict = source_dict[source_id]

         # Extract parameter values and bounds
         r, l, u = _extract_param_range(indi_source_dict[:z_s])
      else
         # Get the source redshift from the file; treat it as fixed (lower == upper == refer)
         r = source_data[1, 5]
         l = r
         u = r
      end

      if sample_z == :zs
         # Source-distance parameter is the redshift itself
         push!(params, Parameter(owner = source_id, 
                                 name  = Symbol(:zs, i), 
                                 refer = r, 
                                 lower = l, 
                                 upper = u))
      else
         # Source-distance parameter is the distance ratio: convert redshift to adis
         # (uniform prior in adis - see _resolve_sample_z)
         adis_r, adis_l, adis_u = _zs_to_adis(cosmo, z_d, r, l, u)

         # Add redshift parameter to the reference, lower, and upper vectors
         push!(params, Parameter(owner = source_id, 
                                 name  = Symbol(:adis, i), 
                                 refer = adis_r, 
                                 lower = adis_l, 
                                 upper = adis_u))
      end

      # Run over knots
      n_knot = Int64(maximum(source_data[:, 2]))
      knots = Vector{Knot}(undef, n_knot)
      for k in 1:n_knot
         # Generate knot id
         knot_id = Symbol(:knot, k)

         # Get individual knot data from file
         mask = source_data[:, 2] .== k
         knot_data = source_data[mask, :]

         # Get knot parameters
         x = knot_data[:, 3]
         y = knot_data[:, 4]

         σx = knot_data[:, 6]
         σy = knot_data[:, 7]
         σθ = knot_data[:, 8]

         # Read optional parity values from column 9: +1 / -1, with 0 marking images without
         # parity information. A single measured image is already informative (absolute quantity).
         p = zeros(Int64, length(x))
         if n_col >= 9
            p_temp = knot_data[:, 9]
            if !all(v -> v in (-1.0, 0.0, 1.0), p_temp)
               error("** parity ** values (column 9) must be +1, -1, or 0 (no information) in source-$i, knot-$k")
            end
            p = Int64.(p_temp)
         end
         if any(!=(0), p)
            any_parity = true
         end
         
         # Read optional flux measurements (magnitudes). sigma <= 0 marks unmeasured images;
         # at least two measured images are required (relative measurement).
         m  = Float64[]
         σm = Float64[]
         if n_col >= 11
            σm_temp = knot_data[:, 11]
            n_meas  = count(>(0.0), σm_temp)
            if n_meas >= 2
               m  = knot_data[:, 10]
               σm = σm_temp
               any_flux = true
            elseif n_meas == 1
               @warn "Only one measured image for ** mag ** in source-$i, knot-$k - ignored (needs >= 2)."
            end
         end

         # Read optional time-delay measurements (days). Same rules as flux above.
         td   = Float64[]
         σ_td = Float64[]
         if n_col >= 13
            σtd_temp = knot_data[:, 13]
            n_meas   = count(>(0.0), σtd_temp)
            if n_meas >= 2
               td   = knot_data[:, 12]
               σ_td = σtd_temp
               any_td = true
            elseif n_meas == 1
               @warn "Only one measured image for ** td ** in source-$i, knot-$k - ignored (needs >= 2)."
            end
         end

         # Convert to arcsec if reference is not (0, 0)
         x_arcsec, y_arcsec = _to_arcsec(observation, x, y)
         knots[k] = Knot(x      = x_arcsec, 
                         σx     = σx, 
                         y      = y_arcsec, 
                         σy     = σy, 
                         σθ     = deg2rad.(σθ), 
                         parity = p,
                         m      = m, 
                         σm     = σm, 
                         td     = td, 
                         σ_td   = σ_td)
      end
      sources[i] = Source(knots=knots)
   end
   # Parity / flux / time-delay terms are OPT-IN: used only if the corresponding use_* flag is
   # explicitly set to true in the source section AND the data is present.
   # Parity
   if get(source_dict, :use_parity, false)
      use_parity = any_parity
   else
      use_parity = false
   end

   # Flux
   if get(source_dict, :use_flux, false)
      use_flux = any_flux
   else
      use_flux = false
   end

   # Time Delay
   if get(source_dict, :use_time_delay, false)
      use_time_delay = any_td
   else
      use_time_delay = false
   end

   if get(source_dict, :use_parity, false) === true && !any_parity
      @warn "use_parity is set but no knot provides ** parity ** values - parity term disabled."
   end

   if get(source_dict, :use_flux, false) === true && !any_flux
      @warn "use_flux is set but no knot provides ** mag / sigma_mag ** - flux chi2 disabled."
   end

   if get(source_dict, :use_time_delay, false) === true && !any_td
      @warn "use_time_delay is set but no knot provides ** td / sigma_td ** - time-delay chi2 disabled."
   end

   # Parity penalty strength (default defined here, override with parity_force in the source section)
   parity_force = Float64(get(source_dict, :parity_force, 1000.0))

   # Finite-difference step (arcsec) used in the flux chi2 magnification correction
   pixel_fd = Float64(get(source_dict, :pixel_fd, observation.pixel_scale))

   return SourceConfig(sources        = sources, 
                       use_parity     = use_parity, 
                       parity_force   = parity_force,
                       use_flux       = use_flux, 
                       use_time_delay = use_time_delay, 
                       pixel_fd       = pixel_fd)
end

function _source!(dict::Dict, cosmo::Cosmology.AbstractCosmology, observation::Observation, params::Vector{Parameter}, sample_z::Symbol)
   # Get source dictionary from the full dictionary
   source_dict = dict[:source]

   # Determine if we need to read source details from a file
   from_file = get(source_dict, :from_file, false)

   if from_file
      source_config = _source_from_file!(dict, cosmo, observation, params, sample_z)
   else
      source_config = _source_direct!(dict, cosmo, observation, params, sample_z)
   end
   return source_config
end
