function _build_kernel(model::ModelConfig)
   # Get pixel grid
   pixels = model.lens_config.pixels

   # Stack all knot positions and remember which columns belong to each knot
   θx     = Float64[]
   θy     = Float64[]
   ranges = UnitRange{Int64}[]
   for src in model.source_config.sources
      for knot in src.knots
         first = length(θx) + 1
         append!(θx, knot.x)
         append!(θy, knot.y)
         push!(ranges, first:length(θx))
      end
   end

   # Allocate and fill the kernel matrices
   n_pixel = length(pixels.κ)
   n_point = length(θx)
   
   if model.source_config.use_time_delay
      Kψ = Matrix{Float64}(undef, n_pix, n)
   else
      Kψ = nothing
   end
   Kαx  = Matrix{Float64}(undef, n_pix, n)
   Kαy  = Matrix{Float64}(undef, n_pix, n)
   Kψxx = Matrix{Float64}(undef, n_pix, n)
   Kψyy = Matrix{Float64}(undef, n_pix, n)
   Kψxy = Matrix{Float64}(undef, n_pix, n)
   Lenses.MultiPixelLens._kernel!(Kψ, Kαx, Kαy, Kψxx, Kψyy, Kψxy, θx, θy, pixels.x, pixels.y, pixels.pixel_scale)

   return (ψ=Kψ, αx=Kαx, αy=Kαy, ψxx=Kψxx, ψyy=Kψyy, ψxy=Kψxy, ranges=ranges)
end


function _lens_quantities(K::NamedTuple, κ::Matrix{Float64})
   # pixels in kernel order (no copy)
   κv  = vec(κ)

   # Lensing quantities
   αx  = transpose(K.αx)  * κv
   αy  = transpose(K.αy)  * κv
   ψxx = transpose(K.ψxx) * κv
   ψyy = transpose(K.ψyy) * κv
   ψxy = transpose(K.ψxy) * κv

   n_knots = length(K.ranges)
   αx_all  = Vector{Vector{Float64}}(undef, n_knots)
   αy_all  = Vector{Vector{Float64}}(undef, n_knots)
   A_all   = Vector{NTuple{4, Vector{Float64}}}(undef, n_knots)
   for kid in 1:n_knots
      r = K.ranges[kid]
      αx_all[kid] = αx[r]
      αy_all[kid] = αy[r]
      A_all[kid]  = (ψxx[r], ψxy[r], ψxy[r], ψyy[r])
   end
   return αx_all, αy_all, A_all
end


function _lens_potential(K::NamedTuple, κ::Matrix{Float64})
   ψ = transpose(K.ψ) * vec(κ)
   return [ψ[r] for r in K.ranges]
end


function run_optimizer(model::ModelConfig, param_ref::Dict{Tuple{Symbol,Symbol}, <:Real}, method_cfg::MaxEntOptimizer, opt_cfg::Dict{Symbol, Any})
   
end
