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
export get_potential
export get_deflection
export get_jacobian
export get_time_delay


function get_potential(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, zs::RV) where T <: Union{RV, ROA}
   
end

function get_deflection(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, zs::RV) where T <: Union{RV, ROA}
   # If RV is passed, covert to vector
   θx = isa(θx, RV) ? [θx] : θx
   θy = isa(θy, RV) ? [θy] : θy

   # Check if the input coordinates are of the same type and size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same type and size."))
   end

   # Initialize zero-valued deflection components
   ψx = zero(θx)
   ψy = zero(θx)

   # Get the number of lens planes
   n_p = lens.n_p

   # Get distance matrix
   D_ij = _distances(cosmology, lens, zs)

   # Temporary position and deflection vectors
   θx_j = zero(θx)
   θy_j = zero(θy)

   ψx_j = zero(θx)
   ψy_j = zero(θy)

   for ni in n_p
      # Position in the 1-st lens plane
      θx_j .= θx
      θy_j .= θy

      # Position in i-th lens plane
      if ni == 1
         ψx_j, ψy_j = Lenses.get_deflection(lens._plane_[ni], θx_j, θy_j)
      else
         for nj in 1:ni-1
            θx_j .= θx_j .- (D_ij[nj+1, ni+1] / D_ij[1, ni+1]) .* ψx_j
            θy_j .= θy_j .- (D_ij[nj+1, ni+1] / D_ij[1, ni+1]) .* ψy_j
      
            # Deflection vector at (θx_i, θy_j) in i-th lens plane
            ψx_j, ψy_j = Lenses.get_deflection(lens._plane_[ni], θx_j, θy_j)
         end
      end
      # Total deflection vector
      ψx .+= (D_ij[ni+1, n_p+2] / D_ij[1, n_p+2]) .* ψx_j
      ψy .+= (D_ij[ni+1, n_p+2] / D_ij[1, n_p+2]) .* ψy_j
   end
   return ψx, ψy
end


function get_jacobian(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, zs::RV) where T <: Union{RV, ROA}
   # If RV is passed, covert to vector
   θx = isa(θx, RV) ? [θx] : θx
   θy = isa(θy, RV) ? [θy] : θy

   # Check if the input coordinates are of the same type and size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same type and size."))
   end

   # Initialize zero valued deformation tensor components
   ψxx::ROA = zero(θx)
   ψxy::ROA = zero(θx)
   ψyx::ROA = zero(θy)
   ψyy::ROA = zero(θy)

   # Get the number of lens planes
   n_p = lens.n_p

   # Get distance matrix
   D_ij = _distances(cosmology, lens, zs)
   
   # Temporary position and deflection vectors
   θx_j = zero(θx)
   θy_j = zero(θy)
   ψx_j = zero(θx)
   ψy_j = zero(θy)

   Ui = zeros(size(θx,1), size(θx,2), 2, 2)
   Uj = zeros(size(θx,1), size(θx,2), 2, 2)
   Ai = zeros(size(θx,1), size(θx,2), 2, 2)
   Aj = zeros(size(θx,1), size(θx,2), 2, 2)

   for ni in 1:n_p
      # Initialize A matrix in 1-st lens plane
      Ai[:,:,1,1] .= 1
      Ai[:,:,2,2] .= 1
      Ai[:,:,1,2] .= 0
      Ai[:,:,2,1] .= 0

      # Position in the 1-st lens plane
      θx_j .= θx
      θy_j .= θy

      if ni == 1
         # Deflection vector at (θx_j, θy_j) in i-th lens plane
         ψx_j, ψy_j = Lenses.get_deflection(lens._plane_[ni], θx_j, θy_j)

         # Deformation tensor at (θx_tmp, θy_tmp) in i-th lens plane
         Ui[:,:,1,1], Ui[:,:,2,2], Ui[:,:,1,2] = Lenses.get_jacobian(lens._plane_[ni], θx_j, θy_j)
         Ui[:,:,2,1] .= Ui[:,:,1,2]
      else
         for nj in 1:ni-1
            # Position in i-th lens plane
            θx_j .= θx_j .- (D_ij[nj+1, ni+1] / D_ij[1, ni+1]) .* ψx_j
            θy_j .= θy_j .- (D_ij[nj+1, ni+1] / D_ij[1, ni+1]) .* ψy_j
         
            # Deflection vector at (θx_j, θy_j) in i-th lens plane
            ψx_j, ψy_j = Lenses.get_deflection(lens._plane_[ni], θx_j, θy_j)

            # Deformation tensor at (θx_tmp, θy_tmp) in i-th lens plane
            Ui[:,:,1,1], Ui[:,:,2,2], Ui[:,:,1,2] = Lenses.get_jacobian(lens._plane_[ni], θx_j, θy_j)
            Ui[:,:,2,1] .= Ui[:,:,1,2]


            # Get the A matrix in i-th plane
            @. Ai[:,:,1,1] = Ai[:,:,1,1] - (D_ij[nj+1, ni+1] / D_ij[1, ni+1]) * (Uj[:,:,1,1]*Aj[:,:,1,1] + Uj[:,:,1,2]*Aj[:,:,2,1])
            @. Ai[:,:,1,2] = Ai[:,:,1,2] - (D_ij[nj+1, ni+1] / D_ij[1, ni+1]) * (Uj[:,:,1,1]*Aj[:,:,1,2] + Uj[:,:,1,2]*Aj[:,:,2,2])
            @. Ai[:,:,2,1] = Ai[:,:,2,1] - (D_ij[nj+1, ni+1] / D_ij[1, ni+1]) * (Uj[:,:,2,1]*Aj[:,:,1,1] + Uj[:,:,2,2]*Aj[:,:,2,1])
            @. Ai[:,:,2,2] = Ai[:,:,2,2] - (D_ij[nj+1, ni+1] / D_ij[1, ni+1]) * (Uj[:,:,2,1]*Aj[:,:,1,2] + Uj[:,:,2,2]*Aj[:,:,2,2])
         end
      end
      # Resultant deformation tensor
      @. ψxx = ψxx + (D_ij[ni+1, n_p+2] / D_ij[1, n_p+2]) * (Ui[:,:,1,1]*Ai[:,:,1,1] + Ui[:,:,1,2]*Ai[:,:,2,1])
      @. ψxy = ψxy + (D_ij[ni+1, n_p+2] / D_ij[1, n_p+2]) * (Ui[:,:,1,1]*Ai[:,:,1,2] + Ui[:,:,1,2]*Ai[:,:,2,2])
      @. ψyx = ψyx + (D_ij[ni+1, n_p+2] / D_ij[1, n_p+2]) * (Ui[:,:,2,1]*Ai[:,:,1,1] + Ui[:,:,2,2]*Ai[:,:,2,1])
      @. ψyy = ψyy + (D_ij[ni+1, n_p+2] / D_ij[1, n_p+2]) * (Ui[:,:,2,1]*Ai[:,:,1,2] + Ui[:,:,2,2]*Ai[:,:,2,2])
   end
   return ψxx, ψyy, ψxy, ψyx
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