module Cosmology

# Using inbuilt Julia packages
using QuadGK

# LensFactory modules to import
include("../LensFactoryUtils/Constants.jl")
using .Constants

# Functions to export
export init_cosmology
export scale_factor
export hubble_time
export hubble_distance
export hubble_parameter
export Omega_mz, Omega_rz, Omega_wz, Omega_kz, rho_cz
export lookback_time
export comoving_distance_radial
export comoving_distance_transverse
export angular_diameter_distance
export luminosity_distance


# Initialize an abstract type for cosmology
abstract type AbstractCosmology end


"""
   init_cosmology

An Abstract type to initialize cosmology.

# Arguments
- `H0::Union{Int64, Float64} = 70`: Hubble constant, in **km/s/Mpc**
- `w::Union{Int64, Float64} = -1.0`: Dark energy EOS parameter
- `Ω_m0::Union{Int64, Float64} = 0.3`: Present matter density parameter
- `Ω_r0::Union{Int64, Float64} = 0.0`: Present radiation density parameter
- `Ω_w0::Union{Int64, Float64} = 0.7`: Present dark energy density parameter

# Returns
- `cosmology::AbstractCosmology`: Cosmology instance defined by the input parameters
"""
@kwdef struct init_cosmology <: AbstractCosmology
   H0::RV = 70.0
   w::RV = -1.0
   Ω_m0::RV = 0.3
   Ω_r0::RV = 0.0
   Ω_w0::RV = 0.7
   Ω_k0::RV = 1.0 - Ω_m0 - Ω_r0 - Ω_w0
end


""" 
   scale_factor(z)

Calculates the scale factor at a given redshift.

> Formula: `a = 1/(1+z)`

# Arguments
- `z::Union{Int64, Float64}`: Redshift

# Returns
- `Union{Int64, Float64}`: Scale factor
"""
function scale_factor(z::RV)::RV
   return 1.0 / (1.0 + z)
end

# E(z) parameter
function Ez(cosmo::AbstractCosmology, z::RV)::RV
      return sqrt(cosmo.Ω_m0 * (1.0 + z)^3 +  
                  cosmo.Ω_r0 * (1.0 + z)^4 + 
                  cosmo.Ω_k0 * (1.0 + z)^2 + 
                  cosmo.Ω_w0 * (1.0 + z)^(3.0 * (1.0 + cosmo.w)))
end


""" 
   hubble_time(H0)

Calculates age of the Universe (t_0).

> Formula: t_0 = 1.0 / H0

# Arguments
- `H0::Union{Int64, Float64}`: Hubble constant, in **km/s/Mpc**

# Returns
- `Union{Int64, Float64}`: Age of the Universe, in **Gyr**
"""
function hubble_time(H0::RV)::RV
   return DIST_MPC / H0 / 1.0E3 / YEAR2SECOND
end


""" 
   hubble_distance(H0)

Calculates radius of the observable Universe (R_H).

> Formula: ``R_H = \\frac{c}{H0}``

# Arguments
- `H0::Union{Int64, Float64}`: Hubble constant, in **km/s/Mpc**

# Returns
- `Union{Int64, Float64}`: Radius of the observable Universe, in **meters** 
"""
function hubble_distance(H0::RV)::RV
   return CONST_C * DIST_MPC / H0 / 1.0E3
end


"""
   hubble_parameter(cosmology, z)

Calculates Hubble parameter at redshift z.

> Formula: 

# Arguments
- `cosmology::AbstractCosmology`: Cosmology instance
- `z::Union{Int64, Float64}`: Input redshift

# Returns
- `Union{Int64, Float64}`: Hubble parameter at redshift z, in **km/s/Mpc**
"""
function hubble_parameter(cosmo::AbstractCosmology, z::RV)::RV
   return cosmo.H0 * Ez(cosmo, z)
end


"""
   rho_cz(cosmology, z)

Calculates the critical density of the Universe at redshift `z`.

> Formula: `\\rho_c(z) = \\frac{3 H^2(z)}{8\\pi G}`

# Arguments
- `cosmology::AbstractCosmology`: Cosmology instance
- `z::Union{Int64, Float64}`: Input redshift

# Returns
- `Union{Int64, Float64}`: Critical density at redshift `z`
"""
function rho_cz(cosmo::AbstractCosmology, z::RV)::RV
   return (3.0 / 8.0 / pi / CONST_G) * hubble_parameter(cosmo, z)^2 * (1.0E3 / DIST_MPC)^2
end


"""
   Omega_mz(cosmology, z)

Calculates the matter density parameter at redshift `z`.

> Formula: `Ω_mz = Ω_m0 * (1+z)^3 * (H0/H(z))^2`

# Arguments
- `cosmology::AbstractCosmology`: Cosmology instance
- `z::Union{Int64, Float64}`: Input redshift

# Returns
- `Union{Int64, Float64}`: Matter density parameter at redshift `z`
"""
function Omega_mz(cosmo::AbstractCosmology, z::RV)::RV
   return cosmo.Ω_m0 * (1 + z)^3 * (cosmo.H0 / hubble_parameter(cosmo, z))^2
end


"""
   Omega_rz(cosmology, z)

Calculates the radiation density parameter at redshift `z`.

> Formula: `Ω_rz = Ω_r0 * (1+z)^4 * (H0/H(z))^2`

# Arguments
- `cosmology::AbstractCosmology`: Cosmology instance
- `z::Union{Int64, Float64}`: Input redshift

# Returns
- `Union{Int64, Float64}`: Radiation density parameter at redshift `z`
"""
function Omega_rz(cosmo::AbstractCosmology, z::RV)::RV
   return cosmo.Ω_r0 * (1 + z)^4 * (cosmo.H0 / hubble_parameter(cosmo, z))^2
end 


"""
   Omega_wz(cosmology, z)

Calculates the dark energy density parameter at redshift `z`.

> Formula: `Ω_wz = Ω_w0 * (1+z)^(3*(1+w)) * (H0/H(z))^2`

# Arguments
- `cosmology::AbstractCosmology`: Cosmology instance
- `z::Union{Int64, Float64}`: Input redshift

# Returns
- `Union{Int64, Float64}`: Dark energy density parameter at redshift `z`
"""
function Omega_wz(cosmo::AbstractCosmology, z::RV)::RV
   return cosmo.Ω_w0 * (1 + z)^(3.0 * (1.0 + cosmo.w)) * (cosmo.H0/hubble_parameter(cosmo, z))^2
end


"""
   Omega_kz(cosmology, z)

Calculates the curvature density parameter at redshift `z`.

> Formula: `Ω_kz = Ω_k0 * (1+z)^2 * (H0/H(z))^2`

# Arguments
- `cosmology::AbstractCosmology`: Cosmology instance
- `z::Union{Int64, Float64}`: Input redshift

# Returns
- `Union{Int64, Float64}`: Curvature density parameter at redshift `z`
"""
function Omega_kz(cosmo::AbstractCosmology, z::RV)::RV
   return cosmo.Ω_k0 * (1 + z)^2 * (cosmo.H0/hubble_parameter(cosmo, z))^2
end


"""
   lookback_time(cosmology, z)

Calculates lookbak time for a given redshift \$z\$

> Formula: \$ t_L(z) = \\int_0^z \\frac{1}{(1+z')H(z')} dz' \$

# Arguments
- `cosmology::AbstractCosmology`: Cosmology instance
- `z::Union{Int64, Float64}`: Input redshift

# Returns
- `Union{Int64, Float64}`: Lookback time for a source at redshift \$z\$
"""
function lookback_time(cosmo::AbstractCosmology, z::RV)::RV
   # Hubble time (in years)
   tH::RV = hubble_time(cosmo.H0)

   # Integral part
   function integrand(cosmo, z)
      return 1.0 / ( (1.0 + z) * Ez(cosmo, z) )
   end
   tC, tC_err = quadgk(x -> integrand(cosmo, x), 0.0, z)

   # Throw error if integration error is large
   if tC_err > 1E-8
      throw(ArgumentError("Error in lookback time estimation is $tC_err."))
   end
   return tH * tC
end


# Comoving distance (radial)
function comoving_distance_radial(cosmo::AbstractCosmology, zi::RV, zf::RV)::RV
   # Hubble distance
   dH::RV = hubble_distance(cosmo.H0)

   # Integral part
   function integrand(cosmo, z)
      return 1.0 / Ez(cosmo, z)
   end
   dC, dC_err = quadgk(x -> integrand(cosmo, x), zi, zf)

   # Throw error if integration error is large
   if dC_err > 1E-8
      throw(ArgumentError("Error in comoving distance estimation is $dC_err."))
   end
   dC = dH * dC
   return dC
end


# Comoving distance (transverse)
function comoving_distance_transverse(cosmo::AbstractCosmology, zi::RV, zf::RV)::RV
   # Get the Hubble distance
   dH::RV = hubble_distance(cosmo.H0)

   # Get the comoving distance (radial)
   dC::RV = comoving_distance_radial(cosmo, zi, zf)

   # Get the comoving distance (transverse)
   dM::RV = 0.0
   if cosmo.Ω_k0 ≠ 0
      Ωχ = √( abs(cosmo.Ω_k0) ) * dC/dH
      if cosmo.Ω_k0 > 0 
         dM = dH * sinh(Ωχ) / √( abs(cosmo.Ω_k0) )
      else
         dM =  dH * sin(Ωχ) / √( abs(cosmo.Ω_k0) )
      end
   else
      dM = dC
   end
   return dM
end


# Angular diameter distance
function angular_diameter_distance(cosmo::AbstractCosmology, zi::RV, zf::RV)::RV
   return comoving_distance_transverse(cosmo, zi, zf) / (1.0 + zf)
end


# Luminosity distance
function luminosity_distance(cosmo::AbstractCosmology, zf::RV)::RV
   return (1.0 + zf) * comoving_distance_transverse(cosmo, 0.0, zf)
end


# Angular scale (for an arcsecond)
function angular_scale(cosmo::AbstractCosmology, zf::RV)::RV
   return (angular_diameter_distance(cosmo, 0.0, zf) / DIST_KPC) * ANGLE_ARCSEC
end


"""
   distance_modulus(cosmology, z)

Calculates the distance modulus to a given redshift `z`.

> Formula: `\\mu = 5 \\log(d) - 5`

# Arguments
- `cosmology::AbstractCosmology`: Cosmology instance
- `z::Union{Int64, Float64}`: Input redshift

# Returns
- `Union{Int64, Float64}`: Distance modulus to a redshift `z`
"""
function distance_modulus(cosmo::AbstractCosmology, zf::RV)::RV
   return 5.0 * log10(luminosity_distance(cosmo, zf) / DIST_PC) - 5.0
end


"""
   comoving_volume_element(cosmology, z)

Calculates the comving volume element at redshift `z`, in **Gpc^3**.

> Formula: \$ dV = dH \\frac{D_M^2}{E(z)} \$

# Arguments
- `cosmology::AbstractCosmology`: Cosmology instance
- `z::Union{Int64, Float64}`: Input redshift

# Returns
- `Union{Int64, Float64}`: Comoving volume element at redshift `z`, in **Gpc^3**
"""
function comoving_volume_element(cosmo::AbstractCosmology, zf::RV)::RV
   # Get the Hubble distance
   dH::RV = hubble_distance(cosmo.H0)

   # E(z) factor
   ez_value::RV = Ez(cosmo, zf)

   # Get the transverse comoving distance
   dC_trans::RV = comoving_distance_transverse(cosmo, 0.0, zf)

   return dH * dC_trans^2 / ez_value / DIST_GPC^3
end



function comoving_volume(cosmo::AbstractCosmology, zf::RV)::RV
   # Get the Hubble distance
   dH::RV = hubble_distance(cosmo.H0)

   # Get the comoving distance (radial)
   dM::RV = comoving_distance_transverse(cosmo, 0.0, zf)

   # Get the comoving volume up to redshift z
   com_vol::RV = 0.0
   if cosmo.Ω_k0 ≠ 0.0
      term1 = 4.0 * π * dH^3 / 2.0 / cosmo.Ω_k0
      term2 = (dM/dH) * √( 1.0 + cosmo.Ω_k0 * (dM / dH)^2 )
      if cosmo.Ω_k0 > 0.0
         com_vol = term1 * ( term2 - asinh(√(abs(cosmo.Ω_k0)) * dM/dH ) / √(abs(cosmo.Ω_k0) ) )
      else
         com_vol = term1 * ( term2 -  asin(√(abs(cosmo.Ω_k0)) * dM/dH ) / √(abs(cosmo.Ω_k0) ) )
      end
   else
      com_vol = (4.0 * π / 3) * dM^3
   end
   return com_vol / DIST_GPC^3
end


end