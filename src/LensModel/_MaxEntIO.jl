# --------------------------------------------------------------------------------------------------
# Internal functions to read maximum entropy method
# --------------------------------------------------------------------------------------------------
function _generate_map(map_dict::Dict, FOV::Vector, pixel_scale::Real)
   # Get map type
   map_type = Symbol(map_dict[:type])

   # Generate map based on type
   if map_type == :uniform
      # Use meshgrid to create the grid
      x_grid, y_grid = Lenses.get_meshgrid(0.5 * FOV[1], 0.5 * FOV[2], pixel_scale)
      x_pos = vec(x_grid)
      y_pos = vec(y_grid)
      κ_map = fill(map_dict[:value], length(x_pos))
      return x_pos, y_pos, κ_map
   elseif map_type == :gaussian
      error("Gaussian map generation is not implemented yet.")
   elseif map_type == :custom
      return readdlm(map_dict[:file], comments=true, comment_char='#')
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
   multi_pixel_lens = MultiPixelLens(x_pos, y_pos, κ_map)

   # Generate lens component
   for i in 1:n_lenses
      lens_id = Symbol(:lens, i)
      indi_lens_dict = lens_dict[lens_id]
      name = Symbol(indi_lens_dict[:lens])

      # Store lens model name in lens_name vector
      lens_name[i] = LensComponent(owner=lens_id, name=name)
   end

   return MaxExtConfig(; multiplane = false,
                         components = lens_name, 
                         alpha      = alpha, 
                         pixels     = multi_pixel_lens)
end