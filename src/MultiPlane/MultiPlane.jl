"""
    MultiPlane
"""
module MultiPlane

# Julia inbuilt functions to import

# LensFactory modules to use
using ..Constants
using ..Cosmology
using ..Lenses

# Various lensing function to export
export get_deflection
export get_jacobian
export get_time_delay


function get_deflection(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, zs::RV) where T <: Union{RV, ROA}
   # If RV is passed, covert to vector
   θx = isa(θx, RV) ? [θx] : θx
   θy = isa(θy, RV) ? [θy] : θy

   # Initialize zero-valued deflection components
   ψx = zero(θx)
   ψy = zero(θx)

   # Get the number of lens planes
   n_p = lens.n_p

   # Get distance matrix
   D_ij = _distances(cosmology, lens, zs)

   # Temporary deflection vector for each lens plane
   ψ_vec = zeros(2, n_p)
   θx_i, θy_i = 0.0, 0.0

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
               θx_i = θx_i - (D_ij[nj+1, ni+1] / D_ij[1, ni+1]) * ψ_vec[1, nj]
               θy_i = θy_i - (D_ij[nj+1, ni+1] / D_ij[1, ni+1]) * ψ_vec[2, nj]
            end

            # Deflection vector at (θx_i, θy_i) in i-th plane
            (ψ_vec[1, ni],), (ψ_vec[2, ni],) = Lenses.get_deflection(lens._plane_[ni], θx_i, θy_i)

            # Final deflection vector
            ψx[i, j] = ψx[i, j] + (D_ij[ni+1, n_p + 2] / D_ij[1, n_p + 2]) * ψ_vec[1, ni]
            ψy[i, j] = ψy[i, j] + (D_ij[ni+1, n_p + 2] / D_ij[1, n_p + 2]) * ψ_vec[2, ni]
         end
      end
   end
   return ψx, ψy
end


function get_jacobian(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, zs::RV) where T <: Union{RV, ROA}
   # If RV is passed, covert to vector
   θx = isa(θx, RV) ? [θx] : θx
   θy = isa(θy, RV) ? [θy] : θy

   # Initialize zero valued deformation tensor components
   ψxx::ROA = zero(θx)
   ψxy::ROA = zero(θx)
   ψyx::ROA = zero(θy)
   ψyy::ROA = zero(θy)

   # Get the number of lens planes
   n_p = lens.n_p

   # Get distance matrix
   D_ij = _distances(cosmology, lens, zs)
   
   # Temporary deflection vector for each lens plane
   ψ_vec = zeros(2, n_p)
   U_vec = zeros(2, 2, n_p)
   A_vec = zeros(2, 2, n_p)
   θx_i, θy_i = 0.0, 0.0
   ψrr = zeros(2, 2)

   # Loop over all coordinates
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         ψrr .= 0.0
         for ni in 1:n_p
            # Initialize A matrix in 1-st lens plane
            A_vec[:, :, ni] = [1.0 0.0; 0.0 1.0]
      
            # Position in the 1-st lens plane
            θx_i = θx[i, j]
            θy_i = θy[i, j]

            for nj in 1:ni-1
               θx_i = θx_i - (D_ij[nj+1, ni+1] / D_ij[1, ni+1]) * ψ_vec[1, nj]
               θy_i = θy_i - (D_ij[nj+1, ni+1] / D_ij[1, ni+1]) * ψ_vec[2, nj]
               A_vec[:, :, ni] .-= (D_ij[nj+1, ni+1] / D_ij[1, ni+1]) .* (U_vec[:, :, nj] * A_vec[:, :, nj])
            end

            # Deflection vector at (θ_xi, θ_yi) in i-th plane
            (ψ_vec[1, ni],), (ψ_vec[2, ni],) = Lenses.get_deflection(lens._plane_[ni], θx_i, θy_i)

            # U-matrix for i-th lens plane
            (U_vec[1, 1, ni],), (U_vec[2, 2, ni],), (U_vec[1, 2, ni],) = Lenses.get_jacobian(lens._plane_[ni], θx_i, θy_i)
            U_vec[2, 1, ni] = U_vec[1, 2, ni]

            ψrr .+= (D_ij[ni+1, n_p+2] / D_ij[1, n_p+2]) .* (U_vec[:, :, ni] * A_vec[:, :, ni])
         end
         ψxx[i, j], ψyy[i, j], ψxy[i, j], ψyx[i, j] = ψrr[1, 1], ψrr[2, 2], ψrr[1, 2], ψrr[2, 1]
      end
   end
   return ψxx, ψyy, ψxy, ψyx
end


function get_time_delay(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, zs::RV, β::NTuple{2, RV}) where T <: Union{RV, ROA}
   # If RV is passed, covert to vector
   θx = isa(θx, RV) ? [θx] : θx
   θy = isa(θy, RV) ? [θy] : θy

   # Initialize zero-valued time delay function
   ϕ = zero(θx)

   # Get the number of lens planes
   n_p = lens.n_p

   # Get distance matrix
   D_ij = _distances(cosmology, lens, zs)

   # Temporary deflection vector for each lens plane
   ψ_vec = zeros(2, n_p)
   θx_i, θy_i = 0.0, 0.0
   ψ = 0.0

   # Get constant factor array
   constant_factor = zeros(n_p)
   for ni in 1:n_p
      constant_factor[ni] = ( (1.0 + lens.z_d[ni]) / CONST_C ) * (D_ij[1, ni+1] * D_ij[1, ni+2] / D_ij[ni+1, ni+2])
   end
   println(constant_factor)
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
               θx_i = θx_i - (D_ij[nj+1, ni+1] / D_ij[1, ni+1]) * ψ_vec[1, nj]
               θy_i = θy_i - (D_ij[nj+1, ni+1] / D_ij[1, ni+1]) * ψ_vec[2, nj]            
            end   

            # Potential value at (θ_xi, θ_yi) in i-th plane
            ψ = Lenses.get_potential(lens._plane_[ni], θx_i, θy_i)

            # Deflection vector at (θx_i, θy_i) in i-th plane
            (ψ_vec[1, ni],), (ψ_vec[2, ni],) = Lenses.get_deflection(lens._plane_[ni], θx_i, θy_i)

            # Get the position vector in i+1-th plane
            if ni < n_p
               θx_j = θx_i -  (D_ij[ni+1, ni+2] / D_ij[1, ni+2]) * ψ_vec[1, ni]
               θy_j = θy_i -  (D_ij[ni+1, ni+2] / D_ij[1, ni+2]) * ψ_vec[2, ni]
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


function _distances(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, zs::RV)::Matrix{Float64}
   # Vector of all redshift (including observer and source)
   z_all = [0; lens.z_d; zs]

   # Create and fill empty distance matrix
   D_ij = zeros(length(z_all), length(z_all))
   for i in axes(D_ij, 2)
      for j in axes(D_ij, 1)[i+1:end]
         D_ij[i, j] = Cosmology.angular_diameter_distance(cosmology, z_all[i], z_all[j])
      end
   end
   return D_ij
end

end