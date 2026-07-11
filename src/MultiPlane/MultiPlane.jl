module MultiPlane

# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..Constants
using ..Cosmology
using ..Lenses
using ..LFUtils


# --------------------------------------------------------------------------------------------------
# Function to export
# --------------------------------------------------------------------------------------------------
export get_deflection
export get_jacobian
export get_time_delay


# --------------------------------------------------------------------------------------------------
# Plotting functions (see ../../ext folder for functions)
# --------------------------------------------------------------------------------------------------
export plot_image_plane

function plot_image_plane end


# --------------------------------------------------------------------------------------------------
# Helper functions
# --------------------------------------------------------------------------------------------------
function _distances(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, z_s::Real)
   # Vector of all redshift (including observer and source)
   z_all = [0; lens.z_d; z_s]
   
   # Pre-compute distance ratios
   n_p = lens.n_p
   D_ij = zeros(n_p+2, n_p+2)
   adis_ij = Matrix{Float64}(undef, n_p, n_p)
   adis_is = Vector{Float64}(undef, n_p)

   # Fill empty distance matrix
   for i in 1:n_p+2
      for j in i+1:n_p+2
         D_ij[i, j] = Cosmology.angular_diameter_distance(cosmology, z_all[i], z_all[j])
      end
   end

   # Fill distance ratio matrix and vector
   for ni in 1:n_p
      for nj in 1:ni-1
         adis_ij[nj, ni] = D_ij[nj+1, ni+1] / D_ij[1, ni+1]
      end
      adis_is[ni] = D_ij[ni+1, n_p+2] / D_ij[1, n_p+2]
   end
   return D_ij, adis_ij, adis_is
end


# --------------------------------------------------------------------------------------------------
# Main functions
# --------------------------------------------------------------------------------------------------
# Core method with precomputed distance ratios (used by the lens-modeling likelihood)
function get_potential(lens::Lenses.AbstractLens, θx::T, θy::T, adis_ij::Matrix{Float64}, adis_is::Vector{Float64}) where T <: Real
   # Initialize zero potential value
   ψ = 0.0

   # Get the number of lens planes
   n_p = lens.n_p

   # Temporary deflection vector for each lens plane
   ψ_vec = zeros(2, n_p)

   # Loop over all coordinates
   for ni in 1:n_p
      # Position vector in 1-st plane
      θx_i = θx
      θy_i = θy

      # Get the position vector in i-th plane
      for nj in 1:ni-1
         θx_i = θx_i - adis_ij[nj, ni] * ψ_vec[1, nj]
         θy_i = θy_i - adis_ij[nj, ni] * ψ_vec[2, nj]            
      end   

      # Potential value at (θ_xi, θ_yi) in i-th plane
      ψ_i = Lenses.get_potential(lens._plane_[ni], θx_i, θy_i)

      # Deflection vector at (θx_i, θy_i) in i-th plane
      ψ_vec[1, ni], ψ_vec[2, ni] = Lenses.get_deflection(lens._plane_[ni], θx_i, θy_i)
      
      # Update potential
      ψ = ψ + adis_is[ni] * ψ_i
   end
   return ψ
end


"""
    get_potential(cosmology::Cosmology.AbstractCosmology, 
                  lens::Lenses.AbstractLens, 
                  θx::T, θy::T, z_s::Real) where T <: Real
"""
function get_potential(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, z_s::Real) where T <: Real
   # Get distance and distance ratios
   D_ij, adis_ij, adis_is = _distances(cosmology, lens, z_s)

   return get_potential(lens, θx, θy, adis_ij, adis_is)
end

"""
    get_potential(cosmology::Cosmology.AbstractCosmology, 
                  lens::Lenses.AbstractLens, 
                  θx::T, θy::T, z_s::Real) where T <: ROA
Calculate potential at the given angular coordinate(s) for a given multi-plane lens system with 
cosmology and source redshift. The function is designed to handle both single coordinate and array 
of coordinates as input via multiple dispatch, depending on the type of `T (ROA or Real)`. 

# Arguments
- `cosmology::Cosmology.AbstractCosmology`: Cosmology model.
- `lens::Lenses.AbstractLens`: Multi-plane lens model.
- `θx::Union{Real, ROA}`: The x-coordinate(s) (in arcseconds).
- `θy::Union{Real, ROA}`: The y-coordinate(s) (in arcseconds).
- `z_s::Real`: The redshift of the source plane.

# Returns
- `ψ`: The potential at the given coordinates.
"""
function get_potential(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, z_s::Real) where T <: ROA
   # Initialize zero-valued potential 2D array
   ψ = zero(θx)

   # Get the number of lens planes
   n_p = lens.n_p

   # Get distance and distance ratios
   D_ij, adis_ij, adis_is = _distances(cosmology, lens, z_s)

   # Temporary deflection vector for each lens plane
   ψ_vec = zeros(2, n_p)

   # Loop over all coordinates
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         for ni in 1:n_p
            # Position vector in 1-st plane
            θx_i = θx[i, j]
            θy_i = θy[i, j]

            # Get the position vector in i-th plane
            for nj in 1:ni-1
               θx_i = θx_i - adis_ij[nj, ni] * ψ_vec[1, nj]
               θy_i = θy_i - adis_ij[nj, ni] * ψ_vec[2, nj]            
            end   

            # Potential value at (θ_xi, θ_yi) in i-th plane
            ψ_i = Lenses.get_potential(lens._plane_[ni], θx_i, θy_i)

            # Deflection vector at (θx_i, θy_i) in i-th plane
            ψ_vec[1, ni], ψ_vec[2, ni] = Lenses.get_deflection(lens._plane_[ni], θx_i, θy_i)

            # Update potential
            ψ[i, j] = ψ[i, j] + adis_is[ni] * ψ_i
         end
      end
   end
   return ψ
end


# Core method with precomputed distance ratios (used by the lens-modeling likelihood)
function get_deflection(lens::Lenses.AbstractLens, θx::T, θy::T, adis_ij::Matrix{Float64}, adis_is::Vector{Float64}) where T <: Real
   # Initialize zero-valued deflection components
   ψx = 0.0
   ψy = 0.0

   # Get the number of lens planes
   n_p = lens.n_p

   # Temporary deflection vector for each lens plane
   ψ_vec = zeros(2, n_p)
   
   for ni in 1:n_p
      # Position vector in 1-st plane
      θx_i = θx
      θy_i = θy

      # Get the position vector in i-th plane
      for nj in 1:ni-1
         θx_i = θx_i - adis_ij[nj, ni] * ψ_vec[1, nj]
         θy_i = θy_i - adis_ij[nj, ni] * ψ_vec[2, nj]
      end

      # Deflection vector at (θx_i, θy_i) in i-th plane
      ψ_vec[1, ni], ψ_vec[2, ni] = Lenses.get_deflection(lens._plane_[ni], θx_i, θy_i)

      # Final deflection vector
      ψx = ψx + adis_is[ni] * ψ_vec[1, ni]
      ψy = ψy + adis_is[ni] * ψ_vec[2, ni]
   end
   return ψx, ψy
end


"""
    get_deflection(cosmology::Cosmology.AbstractCosmology, 
                   lens::Lenses.AbstractLens, 
                   θx::T, θy::T, z_s::Real) where T <: Real
"""
function get_deflection(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, z_s::Real) where T <: Real
   # Get distance and distance ratios
   D_ij, adis_ij, adis_is = _distances(cosmology, lens, z_s)

   return get_deflection(lens, θx, θy, adis_ij, adis_is)
end

"""
    get_deflection(cosmology::Cosmology.AbstractCosmology, 
                   lens::Lenses.AbstractLens, 
                   θx::T, θy::T, z_s::Real) where T <: ROA
Calculate deflection angle at the given angular coordinate(s) for a given multi-plane lens system 
with cosmology and source redshift. The function is designed to handle both single coordinate and 
array of coordinates as input via multiple dispatch, depending on the type of `T (ROA or Real)`.

# Arguments
- `cosmology::Cosmology.AbstractCosmology`: Cosmology model.
- `lens::Lenses.AbstractLens`: Multi-plane lens model.
- `θx::Union{Real, ROA}`: The x-coordinate(s) (in arcseconds).
- `θy::Union{Real, ROA}`: The y-coordinate(s) (in arcseconds).
- `z_s::Real`: The redshift of the source plane.

# Returns
- `ψx`: The x-component of the deflection angle at the given coordinates.
- `ψy`: The y-component of the deflection angle at the given coordinates.
"""
function get_deflection(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, z_s::Real) where T <: ROA
   # Initialize zero-valued deflection components
   ψx = zero(θx)
   ψy = zero(θx)

   # Get the number of lens planes
   n_p = lens.n_p

   # Get distance and distance ratios
   D_ij, adis_ij, adis_is = _distances(cosmology, lens, z_s)

   # Temporary deflection vector for each lens plane
   ψ_vec = zeros(2, n_p)

   # Loop over all coordinates
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         for ni in 1:n_p
            # Position vector in 1-st plane
            θx_i = θx[i, j]
            θy_i = θy[i, j]
            
            # Get the position vector in i-th plane
            for nj in 1:ni-1
               θx_i = θx_i - adis_ij[nj, ni] * ψ_vec[1, nj]
               θy_i = θy_i - adis_ij[nj, ni] * ψ_vec[2, nj]
            end

            # Deflection vector at (θx_i, θy_i) in i-th plane
            ψ_vec[1, ni], ψ_vec[2, ni] = Lenses.get_deflection(lens._plane_[ni], θx_i, θy_i)

            # Final deflection vector
            ψx[i, j] = ψx[i, j] + adis_is[ni] * ψ_vec[1, ni]
            ψy[i, j] = ψy[i, j] + adis_is[ni] * ψ_vec[2, ni]
         end
      end
   end
   return ψx, ψy
end


# Core method with precomputed distance ratios (used by the lens-modeling likelihood)
function get_jacobian(lens::Lenses.AbstractLens, θx::T, θy::T, adis_ij::Matrix{Float64}, adis_is::Vector{Float64}) where T <: Real
   # Initialize zero valued deformation tensor components
   ψrr = zeros(2, 2)

   # Get the number of lens planes
   n_p = lens.n_p

   # Temporary deflection vector for each lens plane
   ψ_vec = zeros(2, n_p)
   U_vec = zeros(2, 2, n_p)
   A_vec = zeros(2, 2, n_p)

   for ni in 1:n_p
      # Initialize A matrix in 1-st lens plane
      A_vec[:, :, ni] = [1.0 0.0; 0.0 1.0]

      # Position in the 1-st lens plane
      θx_i = θx
      θy_i = θy

      for nj in 1:ni-1
         θx_i = θx_i - adis_ij[nj, ni] * ψ_vec[1, nj]
         θy_i = θy_i - adis_ij[nj, ni] * ψ_vec[2, nj]
         A_vec[:, :, ni] .-= adis_ij[nj, ni] .* (U_vec[:, :, nj] * A_vec[:, :, nj])
      end

      # Deflection vector at (θ_xi, θ_yi) in i-th plane
      ψ_vec[1, ni], ψ_vec[2, ni] = Lenses.get_deflection(lens._plane_[ni], θx_i, θy_i)

      # U-matrix for i-th lens plane
      (U_vec[1, 1, ni],), (U_vec[2, 2, ni],), (U_vec[1, 2, ni],) = Lenses.get_jacobian(lens._plane_[ni], θx_i, θy_i)
      U_vec[2, 1, ni] = U_vec[1, 2, ni]

      ψrr .+= adis_is[ni] .* (U_vec[:, :, ni] * A_vec[:, :, ni])
   end
   return ψrr[1, 1], ψrr[2, 2], ψrr[1, 2], ψrr[2, 1]
end


"""
    get_jacobian(cosmology::Cosmology.AbstractCosmology, 
                 lens::Lenses.AbstractLens, 
                 θx::T, θy::T, z_s::Real) where T <: Real
"""
function get_jacobian(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, z_s::Real) where T <: Real
   # Get distance and distance ratios
   D_ij, adis_ij, adis_is = _distances(cosmology, lens, z_s)
   
   return get_jacobian(lens, θx, θy, adis_ij, adis_is)
end


"""
    get_jacobian(cosmology::Cosmology.AbstractCosmology, 
                 lens::Lenses.AbstractLens, 
                 θx::T, θy::T, z_s::Real) where T <: ROA
Calculate jacobian (i.e., deformation tensor) of the lens mapping at the given angular coordinate(s)
for a given multi-plane lens system with cosmology and source redshift. The function is designed to 
handle both single coordinate and array of coordinates as input via multiple dispatch, depending on 
the type of `T (ROA or Real)`.

# Arguments
- `cosmology::Cosmology.AbstractCosmology`: Cosmology model.
- `lens::Lenses.AbstractLens`: Multi-plane lens model.
- `θx::Union{Real, ROA}`: The x-coordinate(s) (in arcseconds).
- `θy::Union{Real, ROA}`: The y-coordinate(s) (in arcseconds).
- `z_s::Real`: The redshift of the source plane.

# Returns
- `ψxx`: xx-component of the jacobian.
- `ψyy`: yy-component of the jacobian.
- `ψxy`: xy-component of the jacobian.
- `ψyx`: yx-component of the jacobian.
"""
function get_jacobian(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, z_s::Real) where T <: ROA
   # Initialize zero valued deformation tensor components
   ψrr = zeros(size(θx,1), size(θx,2), 4)

   # Get the number of lens planes
   n_p = lens.n_p

   # Get distance and distance ratios
   D_ij, adis_ij, adis_is = _distances(cosmology, lens, z_s)

   # Temporary deflection vector for each lens plane
   ψ_vec = zeros(2, n_p)
   U_vec = zeros(2, 2, n_p)
   A_vec = zeros(2, 2, n_p)
   ψrr_tmp = zeros(2, 2)

   # Loop over all coordinates
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         ψrr_tmp .= 0.0
         for ni in 1:n_p
            # Initialize A matrix in 1-st lens plane
            A_vec[:, :, ni] = [1.0 0.0; 0.0 1.0]
      
            # Position in the 1-st lens plane
            θx_i = θx[i, j]
            θy_i = θy[i, j]

            for nj in 1:ni-1
               θx_i = θx_i - adis_ij[nj, ni] * ψ_vec[1, nj]
               θy_i = θy_i - adis_ij[nj, ni] * ψ_vec[2, nj]
               A_vec[:, :, ni] .-= adis_ij[nj, ni] .* (U_vec[:, :, nj] * A_vec[:, :, nj])
            end

            # Deflection vector at (θ_xi, θ_yi) in i-th plane
            (ψ_vec[1, ni],), (ψ_vec[2, ni],) = Lenses.get_deflection(lens._plane_[ni], θx_i, θy_i)

            # U-matrix for i-th lens plane
            (U_vec[1, 1, ni],), (U_vec[2, 2, ni],), (U_vec[1, 2, ni],) = Lenses.get_jacobian(lens._plane_[ni], θx_i, θy_i)
            U_vec[2, 1, ni] = U_vec[1, 2, ni]

            ψrr_tmp .+= adis_is[ni] .* (U_vec[:, :, ni] * A_vec[:, :, ni])
         end
         ψrr[i, j, 1] = ψrr_tmp[1, 1]
         ψrr[i, j, 2] = ψrr_tmp[2, 2]
         ψrr[i, j, 3] = ψrr_tmp[1, 2]
         ψrr[i, j, 4] = ψrr_tmp[2, 1]
      end
   end
   return ψrr[:,:,1], ψrr[:,:,2], ψrr[:,:,3], ψrr[:,:,4]
end


function get_time_delay(lens::Lenses.AbstractLens, θx::T, θy::T, β::NTuple{2, Real}, adis_ij::Matrix{Float64}, D_ij::Matrix{Float64}) where T <: Real
   # Initialize zero-valued time delay function
   ϕ = 0.0

   # Number of planes covered by the distance matrix (allows truncation at the source)
   n_p = size(D_ij, 1) - 2

   # Temporary deflection vector for each lens plane
   ψ_vec = zeros(2, n_p)

   # Get constant factor array
   constant_factor = zeros(n_p)
   for ni in 1:n_p
      constant_factor[ni] = ( (1.0 + lens.z_d[ni]) / CONST_C ) * (D_ij[1, ni+1] * D_ij[1, ni+2] / D_ij[ni+1, ni+2]) * ANGLE_ARCSEC^2
   end

   # Loop over all planes
   for ni in 1:n_p
      # Position vector in 1-st plane
      θx_i = θx
      θy_i = θy
      θx_j = θx
      θy_j = θy

      # Get the position vector in i-th plane
      for nj in 1:ni-1
         θx_i = θx_i - adis_ij[nj, ni]   * ψ_vec[1, nj]
         θy_i = θy_i - adis_ij[nj, ni]   * ψ_vec[2, nj]
         if ni < n_p
            θx_j = θx_j - adis_ij[nj, ni+1] * ψ_vec[1, nj]
            θy_j = θy_j - adis_ij[nj, ni+1] * ψ_vec[2, nj]
         end
      end   

      # Potential and deflection at (θ_xi, θ_yi) in i-th plane
      ψ = Lenses.get_potential(lens._plane_[ni], θx_i, θy_i)
      ψ_vec[1, ni], ψ_vec[2, ni] = Lenses.get_deflection(lens._plane_[ni], θx_i, θy_i)

      # Get the position vector in i+1-th plane
      if ni < n_p
         θx_j = θx_j - adis_ij[ni, ni+1] * ψ_vec[1, ni]
         θy_j = θy_j - adis_ij[ni, ni+1] * ψ_vec[2, ni]
      else
         θx_j = β[1]
         θy_j = β[2]
      end

      # Time delay factor (without a reference point)
      ϕ = ϕ + constant_factor[ni] * (0.5 * ((θx_i-θx_j)^2 + (θy_i-θy_j)^2) - (D_ij[ni+1, ni+2] / D_ij[1, ni+2]) * ψ[1])
   end
   return ϕ
end


"""
    get_time_delay(cosmology::Cosmology.AbstractCosmology, 
                   lens::Lenses.AbstractLens, 
                   θx::T, θy::T, z_s::Real, β::NTuple{2, <:Real}) where T <: Real
"""
function get_time_delay(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, z_s::Real, β::NTuple{2, Real}) where T <: Real
   # Get distance matrix
   D_ij, adis_ij, adis_is = _distances(cosmology, lens, z_s)

   return get_time_delay(lens, θx, θy, β, adis_ij, D_ij)
end


"""
    get_time_delay(cosmology::Cosmology.AbstractCosmology, 
                   lens::Lenses.AbstractLens, 
                   θx::T, θy::T, z_s::Real, 
                   β::NTuple{2, <:Real}) where T <: ROA
Calculate time delay at the given angular coordinate(s) for a given multi-plane lens system with 
cosmology and source redshift. The function is designed to handle both single coordinate and array
of coordinates as input via multiple dispatch, depending on the type of `T (ROA or Real)`.

# Arguments
- `cosmology::Cosmology.AbstractCosmology`: Cosmology model.
- `lens::Lenses.AbstractLens`: Multi-plane lens model.
- `θx::Union{Real, ROA}`: The x-coordinate(s) (in arcseconds).
- `θy::Union{Real, ROA}`: The y-coordinate(s) (in arcseconds).
- `z_s::Real`: The redshift of the source plane.
- `β::NTuple{2, Real}`: The source plane position (in arcseconds) as a tuple (βx, βy).

# Returns
- `t_d`: The time delay at the given coordinates (in days).
"""
function get_time_delay(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, z_s::Real, β::NTuple{2, <:Real}) where T <: ROA
   # Initialize zero-valued time delay function
   ϕ = zero(θx)

   # Get the number of lens planes
   n_p = lens.n_p

   # Get distance matrix
   D_ij, adis_ij, adis_is = _distances(cosmology, lens, z_s)

   # Temporary deflection vector for each lens plane
   ψ_vec = zeros(2, n_p)
   ψ = 0.0

   # Get constant factor array
   constant_factor = zeros(n_p)
   for ni in 1:n_p
      constant_factor[ni] = ( (1.0 + lens.z_d[ni]) / CONST_C ) * (D_ij[1, ni+1] * D_ij[1, ni+2] / D_ij[ni+1, ni+2]) * ANGLE_ARCSEC^2
   end

   # Loop over all coordinates
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         for ni in 1:n_p
            # Position vector in 1-st plane
            θx_i = θx[i, j]
            θy_i = θy[i, j]
            θx_j = θx[i, j]
            θy_j = θy[i, j]

            # Get the position vector in i-th plane
            for nj in 1:ni-1
               θx_i = θx_i - adis_ij[nj, ni]   * ψ_vec[1, nj]
               θy_i = θy_i - adis_ij[nj, ni]   * ψ_vec[2, nj]
               if ni < n_p
                  θx_j = θx_j - adis_ij[nj, ni+1] * ψ_vec[1, nj]
                  θy_j = θy_j - adis_ij[nj, ni+1] * ψ_vec[2, nj]
               end
            end   

            # Potential and deflection at (θ_xi, θ_yi) in i-th plane
            ψ = Lenses.get_potential(lens._plane_[ni], θx_i, θy_i)
            ψ_vec[1, ni], ψ_vec[2, ni] = Lenses.get_deflection(lens._plane_[ni], θx_i, θy_i)

            # Get the position vector in i+1-th plane
            if ni < n_p
               θx_j = θx_j - adis_ij[ni, ni+1] * ψ_vec[1, ni]
               θy_j = θy_j - adis_ij[ni, ni+1] * ψ_vec[2, ni]
            else
               θx_j = β[1]
               θy_j = β[2]
            end
            
            # Time delay factor (without a reference point)
            ϕ[i, j] += constant_factor[ni] * (0.5 * ((θx_i-θx_j)^2 + (θy_i-θy_j)^2) - (D_ij[ni+1, ni+2] / D_ij[1, ni+2]) * ψ[1])
         end
      end
   end
   return ϕ
end


"""
    get_magnification_image(cosmology::Cosmology.AbstractCosmology, 
                            lens::Lenses.AbstractLens, 
                            θx::T, θy::T, z_s::Real) where T <: Real
Calculate magnification at the given angular coordinate(s) in the image plane for a given 
multi-plane lens system with cosmology and source redshift. The function is designed to handle both 
single coordinate and array of coordinates as input via multiple dispatch, depending on the type 
of `T (ROA or Real)`.

# Arguments
- `cosmology::Cosmology.AbstractCosmology`: Cosmology model.
- `lens::Lenses.AbstractLens`: Multi-plane lens model.
- `θx::Union{Real, ROA}`: The x-coordinate(s) (in arcseconds).
- `θy::Union{Real, ROA}`: The y-coordinate(s) (in arcseconds).
- `z_s::Real`: The redshift of the source plane.

# Returns
- `μ_image`: The magnification at the given coordinates.
"""
function get_magnification_image(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, z_s::Real) where T <: Union{Real, ROA}
   # Get the deformation tensor components
   ψxx, ψyy, ψxy, ψyx = get_jacobian(cosmology, lens, θx, θy, z_s)

   # μ = 1 / det(A)
   return @. 1.0 / (1.0 + ψxx * ψyy - ψxx - ψyy - ψxy * ψyx)
end


"""
    get_magnification_source(cosmology::Cosmology.AbstractCosmology, 
                             lens::Lenses.AbstractLens, 
                             θx::T, θy::T, z_s::Real; 
                             rays_per_pixel::Int64=1) where T <: Matrix{<:Real}
Calculate magnification map in the source plane for a given multi-plane lens system with cosmology 
and source redshift. The function uses a ray-shooting method to compute the magnification map in the
source plane. The `rays_per_pixel` parameter controls the number of rays to shoot per pixel in the 
image plane, which can be increased for higher accuracy at the cost of increased computation time.

# Arguments
- `cosmology::Cosmology.AbstractCosmology`: Cosmology model.
- `lens::Lenses.AbstractLens`: Multi-plane lens model.
- `θx::Matrix{<:Real}`: The x-coordinate grid (in arcseconds).
- `θy::Matrix{<:Real}`: The y-coordinate grid (in arcseconds).
- `z_s::Real`: The redshift of the source plane.

# Keyword Arguments
- `rays_per_pixel::Int64=1`: The average number of rays to shoot per pixel.

# Returns
- `μ_source`: The magnification map in the source plane.
"""
function get_magnification_source(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, z_s::Real; rays_per_pixel::Int64=1) where T <: Matrix{<:Real}
   # Deflection field
   ψx, ψy = get_deflection(cosmology, lens, θx, θy, z_s)

   # Pixel size
   nx, ny = size(θx)
   pixel_h = abs(θx[2, 1] - θx[1, 1])

   # Area per pixel in the image plane
   area_per_ray = 1.0 / rays_per_pixel

   # Initialize an empty source plane magnification map
   μ_source = zero(θx)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for _ in ax2
      @inbounds  for _ in ax1
         # Draw random pixel numbers equal to rays_per_pixel
         rand_x = 1.0 .+ rand(Float64, rays_per_pixel) .* (nx - 1.0)
         rand_y = 1.0 .+ rand(Float64, rays_per_pixel) .* (ny - 1.0)

         for k in 1:rays_per_pixel
            # Get the source plane position
            βx = PolygonOps.interpolation(rand_x[k], rand_y[k], θx) - PolygonOps.interpolation(rand_x[k], rand_y[k], ψx)
            βy = PolygonOps.interpolation(rand_x[k], rand_y[k], θy) - PolygonOps.interpolation(rand_x[k], rand_y[k], ψy)

            # Get the corresponding pixel values
            βx_p = round(Int64, βx/pixel_h + 0.5*nx + 0.5)
            βy_p = round(Int64, βy/pixel_h + 0.5*ny + 0.5)

            # make sure pixel position is within bounds
            if (1 <= βx_p <= nx) && (1 <= βy_p <= ny)
               μ_source[βx_p, βy_p] += area_per_ray
            end
         end
      end
   end
   return μ_source
end


"""
    get_image(cosmology::Cosmology.AbstractCosmology, 
              lens::Lenses.AbstractLens, 
              θx::T, θy::T, z_s::Real, β::NTuple{2, Real}) where T <: Matrix{<:Real}
Calculate image positions for a given source position in the source plane for a given multi-plane 
lens system with cosmology and source redshift. The function uses a contour-finding method to 
compute the image positions by finding the intersection of contours corresponding to the lens 
equation in the image plane.

# Arguments
- `cosmology::Cosmology.AbstractCosmology`: Cosmology model.
- `lens::Lenses.AbstractLens`: Multi-plane lens model.
- `θx::Matrix{<:Real}`: The x-coordinate grid (in arcseconds).
- `θy::Matrix{<:Real}`: The y-coordinate grid (in arcseconds).
- `z_s::Real`: The redshift of the source plane.
- `β::NTuple{2, Real}`: The source plane position (in arcseconds) as a tuple (βx, βy).

# Returns
- `image_position`: A vector of tuples containing the image positions (in arcseconds).
"""
function get_image(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, z_s::Real, β::NTuple{2, <:Real}) where T <: Matrix{<:Real}
   # Get the potential gradient
   ψx, ψy = get_deflection(cosmology, lens, θx, θy, z_s)

   # Get grid for contour
   RXC = ContourFinder.get_contour(θx, θy, β[1] .- θx .+ ψx, 0.0)
   RYC = ContourFinder.get_contour(θx, θy, β[2] .- θy .+ ψy, 0.0)

   # Initialize empty Vector of tuples to store image positions
   image_position::Vector{NTuple{2, Real}} = []
   for contour_1 in RXC
      for contour_2 in RYC
         # Find the intersection points
         intersect_points = IntersectionFinder.get_intersection( first.(contour_1), last.(contour_1), first.(contour_2), last.(contour_2) )
         
         # Store the intersection points in the image_position vector
         for point in intersect_points
            push!(image_position, point)
         end
      end
   end
   return image_position
end



"""
    get_image(cosmology::Cosmology.AbstractCosmology, 
              lens::Lenses.AbstractLens, 
              θx::T, θy::T, z_s::Real, β::T) where T <: Matrix{<:Real}
Calculate image plane map for a given extended source in the source plane for a given multi-plane lens 
system with cosmology and source redshift.

# Arguments
- `cosmology::Cosmology.AbstractCosmology`: Cosmology model.
- `lens::Lenses.AbstractLens`: Multi-plane lens model.
- `θx::Matrix{<:Real}`: The x-coordinate grid (in arcseconds).
- `θy::Matrix{<:Real}`: The y-coordinate grid (in arcseconds).
- `z_s::Real`: The redshift of the source plane.
- `β::Matrix{<:Real}`: The source plane brightness distribution.

# Returns
- `image_map`: The image plane map corresponding to the given extended source.
"""
function get_image(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, z_s::Real, β::T) where T <: Matrix{<:Real}
   # Get the potential gradient
   ψx, ψy = get_deflection(cosmology, lens, θx, θy, z_s)

   # Create an empty image map
   image_map = zero(θx)

   # Grid size
   nx, ny = size(θx)
   pixel_h = abs(θx[2, 1] - θx[1, 1])

   βx = 0.0
   βy = 0.0
   pixel_x::Int64 = 0
   pixel_y::Int64 = 0

   # Loop over the image plane and assign values from source
   ax1, ax2 = axes(θx, 1), axes(θx, 2)

   @inbounds for j in ax2
      @inbounds for i in ax1
         # Get source plane position in radians
         βx = θx[i, j] - ψx[i, j]
         βy = θy[i, j] - ψy[i, j]

         # Get pixel position from radians
         pixel_x = round(Int64, βx / pixel_h + 0.5*nx + 0.5)
         pixel_y = round(Int64, βy / pixel_h + 0.5*ny + 0.5)

         # make sure pixel position is within bounds
         if (1 <= pixel_x <= nx) && (1 <= pixel_y <= ny)
            image_map[i, j] = β[pixel_x, pixel_y]
         end
      end
   end
   return image_map
end


"""
    get_critical_curve(cosmology::Cosmology.AbstractCosmology, 
                       lens::Lenses.AbstractLens, 
                       θx::T, θy::T, z_s::Real) where T <: Matrix{<:Real}
Calculate critical curves in the image plane for a given multi-plane lens system with cosmology and
source redshift. The function uses a contour-finding method to compute the critical curves by 
finding ``det(A) = 0`` contours in the image plane, where ``A`` is the jacobian matrix of the lens 
mapping.

# Arguments
- `cosmology::Cosmology.AbstractCosmology`: Cosmology model.
- `lens::Lenses.AbstractLens`: Multi-plane lens model.
- `θx::Matrix{<:Real}`: The x-coordinate grid (in arcseconds).
- `θy::Matrix{<:Real}`: The y-coordinate grid (in arcseconds).
- `z_s::Real`: The redshift of the source plane.

# Returns
- `critical_curve`: A vector of vectors containing the critical curve coordinates (in arcseconds).
"""
function get_critical_curve(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, z_s::Real) where T <: Matrix{<:Real}
   # Get the jacobian components
   ψxx, ψyy, ψxy, ψyx = get_jacobian(cosmology, lens, θx, θy, z_s)

   # μ det(A)
   detA = 1.0 .+ ψxx .* ψyy .- ψxx .- ψyy .- ψxy .* ψyx

   # Get the zero eigenvalue contours
   critical_curve = ContourFinder.get_contour(θx, θy, detA, 0)

   return critical_curve
end


"""
    get_caustic(cosmology::Cosmology.AbstractCosmology, 
                lens::Lenses.AbstractLens, 
                θx::T, θy::T, z_s::Real) where T <: Matrix{<:Real}
Calculate caustics in the source plane for a given multi-plane lens system with cosmology and source
redshift. The function first computes the critical curves in the image plane and then maps them to 
the source plane using the lens equation to obtain the caustics.

# Arguments
- `cosmology::Cosmology.AbstractCosmology`: Cosmology model.
- `lens::Lenses.AbstractLens`: Multi-plane lens model.
- `θx::Matrix{<:Real}`: The x-coordinate grid (in arcseconds).
- `θy::Matrix{<:Real}`: The y-coordinate grid (in arcseconds).
- `z_s::Real`: The redshift of the source plane.

# Returns
- `caustics_curve`: A vector of vectors containing the caustic curve coordinates (in arcseconds).
"""
function get_caustic(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, z_s::Real) where T <: Matrix{<:Real}
   # Generate critical curves
   critical_curve = get_critical_curve(cosmology, lens, θx, θy, z_s)

   # Get tangential caustics
   caustics_curve = Vector{Vector{Vector{Float64}}}()
   for curve in critical_curve
      ψx, ψy = get_deflection(cosmology, lens, first.(curve), last.(curve), z_s)
      src_x = first.(curve) .- ψx
      src_y =  last.(curve) .- ψy
      push!(caustics_curve, [[x, y] for (x, y) in zip(src_x, src_y)])
   end

   return caustics_curve
end

end