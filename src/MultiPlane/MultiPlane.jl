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


function get_deflection_vector(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::T, θy::T, zs::RV) where T <: ROA
   
end



function _distances(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, zs::RV)::Matrix{Float64}
   # Vector of all redshift (including observer and source)
   z_all = [0; lens._zl_; zs]

   # Create and fill empty distance matrix
   D_ij = zeros(length(z_all), length(z_all))
   for i in axes(D_ij, 2)
      for j in axes(D_ij, 1)[i+1:end]
         println(z_all[i], " ", z_all[j])
         D_ij[i, j] = Cosmology.angular_diameter_distance(cosmology, z_all[i], z_all[j])
      end
   end
   return D_ij
end

end