# --------------------------------------------------------------------------------------------------
# Internal functions to read maximum entropy method
# --------------------------------------------------------------------------------------------------
function _generate_map(map_dict::Dict, FOV::Union{Vector{<:Real}, NTuple{2, <:Real}}, pixel_scale::Real)
   # Get map type
   map_type = Symbol(map_dict[:type])

   # Generate map based on type
   if map_type == :uniform
      # Use meshgrid to create the grid
      x_grid, y_grid = Lenses.get_meshgrid(0.5 * FOV[1], 0.5 * FOV[2], pixel_scale)
      κ_map = fill(Float64(map_dict[:value]), size(x_grid))
      return x_grid, y_grid, κ_map
   elseif map_type == :gaussian
      error("Gaussian map generation is not implemented yet.")
   elseif map_type == :custom
      # Check if file exist
      _require(map_dict, :file)
      file = String(map_dict[:file])
      if !isfile(file)
         error("Custom default map file not found: $file")
      end

      x_grid, y_grid, κ_map = JLD2.jldopen(file, "r") do f
         for key in ("x", "y", "kappa")
            if !haskey(f, key)
               error("Custom default map $file must contain `x`, `y` and `kappa` (found: $(join(keys(f), ", "))).")
            end
         end
         Matrix{Float64}(f["x"]), Matrix{Float64}(f["y"]), Matrix{Float64}(f["kappa"])
      end

      # Same shape for all three
      if !(size(x_grid) == size(y_grid) == size(κ_map))
         error("In $file, x, y and kappa must have the same size; got $(size(x_grid)), $(size(y_grid)), $(size(κ_map)).")
      end

      # Regular grid in get_meshgrid layout with spacing pixel_scale (throws ArgumentError otherwise)
      Lenses.MultiPixelLens._check_grid(x_grid, y_grid, pixel_scale)

      # Positive map (required by the entropy prior)
      if any(v -> !(v > 0.0) || !isfinite(v), κ_map)
         error("Custom default map in $file must be finite and positive in every pixel.")
      end
      return x_grid, y_grid, κ_map
   else
      error("Unsupported default map type: $map_type")
   end
end


function _maxent(lens_dict::Dict, params::Vector{Parameter}, observation::Observation)
   # Construct a composite lens using initial values
   n_lenses = lens_dict[:total_lenses]
   if n_lenses != 1
      error("Number of lenses should be equal to 1 for MaxEnt method.")
   end
   # Initialize lens name vector
   lens_name = Vector{LensComponent}(undef, n_lenses)

   # Read regularization parameter
   _require(lens_dict, :alpha)
   alpha = lens_dict[:alpha]

   # Get FOV and pixel scale for the lens
   if !haskey(lens_dict[:lens1], :FOV) || !haskey(lens_dict[:lens1], :pixel_scale)
      @warn "FOV and pixel_scale are not provided for MaxEnt method. Using default values."
      FOV = observation.FOV
      pixel_scale = observation.pixel_scale
   else
      FOV = lens_dict[:lens1][:FOV]
      pixel_scale = lens_dict[:lens1][:pixel_scale]
   end

   # Generate inital map
   default_map = lens_dict[:lens1][:default_map]
   if default_map === nothing
      error("Default map is not provided for MaxEnt method.")
   end
   x_pos, y_pos, κ_map = _generate_map(default_map, FOV, pixel_scale)
   multi_pixel_lens    = PixelComponent(x           = x_pos, 
                                        y           = y_pos, 
                                        κ           = κ_map,
                                        pixel_scale = pixel_scale)

   # Generate lens component
   for i in 1:n_lenses
      lens_id = Symbol(:lens, i)
      indi_lens_dict = lens_dict[lens_id]
      name = Symbol(indi_lens_dict[:lens])

      # Store lens model name in lens_name vector
      lens_name[i] = LensComponent(owner=lens_id, name=name)
   end

   return MaxEntLens(; multiplane   = false,
                         components = lens_name, 
                         alpha      = alpha, 
                         pixels     = multi_pixel_lens)
end