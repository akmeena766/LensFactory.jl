module Cosmology

# Using inbuilt Julia packages
using QuadGK

# LensFactory modules to import
using ..Constants

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
export distance_modulus, angular_scale
export comoving_volume_element, comoving_volume
export adis_to_zs


# Initialize an abstract type for cosmology
abstract type AbstractCosmology end


"""
    init_cosmology(H0::RV=70.0, w::RV=-1.0, Omega_m0::RV=0.3, Omega_r0::RV=0.0, Omega_w0::RV=0.7)
An Abstract type to initialize a cosmology with given parameters.

# Arguments
   - `H0::RV`: Hubble constant in ``{\\rm \\mathbf{km/s/Mpc}}``.
   - `w::RV`: Dark energy equation of state parameter.
   - `Omega_m0::RV`: Matter density parameter.
   - `Omega_r0::RV`: Radiation density parameter.
   - `Omega_w0::RV`: Dark energy density parameter.

# Returns
   - `Omega_k0`: Curvature density parameter.
"""
@kwdef struct init_cosmology <: AbstractCosmology
   H0::RV = 70.0
   w::RV = -1.0
   Omega_m0::RV = 0.3
   Omega_r0::RV = 0.0
   Omega_w0::RV = 0.7
   Omega_k0::RV = 1.0 - Omega_m0 - Omega_r0 - Omega_w0
end


""" 
    scale_factor(z::RV) --> RV

Calculate scale factor ``(a)`` at redshift ``z``,
```math
a = \\frac{1}{1+z}.
```

# Arguments
   - `z::RV`: Redshift.

# Returns
   - `a::RV`: Scale factor.
"""
function scale_factor(z::RV)
   return 1.0 / (1.0 + z)
end


"""
    Ez(cosmology::AbstractCosmology, z::RV) --> RV
Calculate dimensionless Hubble parameter ``(E)`` at redshift, ``z``,
```math
E(z) = \\sqrt{ Ω_{m0} (1+z)^2 + Ω_{r0} (1+z)^4 + Ω_{k0} (1+z)^2 + Ω_{w0} (1+z)^{3(1+w)} }.
```

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z::RV`: Redshift.

# Returns
   - `Ez::RV`: Dimensionless Hubble parameter.
"""
function Ez(cosmology::AbstractCosmology, z::RV)
      return sqrt(cosmology.Omega_m0 * (1.0 + z)^3 +  
                  cosmology.Omega_r0 * (1.0 + z)^4 + 
                  cosmology.Omega_k0 * (1.0 + z)^2 + 
                  cosmology.Omega_w0 * (1.0 + z)^(3.0 * (1.0 + cosmology.w)))
end


"""
    hubble_parameter(cosmology::AbstractCosmology, z::RV) --> RV
Calculate Hubble parameter ``(H)`` at redshift ``z`` in ``{\\rm \\mathbf{km/s/Mpc}}``,
```math
H(z) = H_0 \\: E(z).
```

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z::RV`: Redshift.

# Returns
   - `H::RV`: Hubble parameter.
"""
function hubble_parameter(cosmology::AbstractCosmology, z::RV)
   return cosmology.H0 * Ez(cosmology, z)
end


""" 
    hubble_time(H0::RV) --> RV
Calculate Hubble time ``(t_H)`` in ``{\\rm \\mathbf{Gyr}}``,
```math
t_H = \\frac{1}{H_0}.
```

# Arguments
   - `H0::RV`: Hubble constant in ``{\\rm \\mathbf{km/s/Mpc}}``.

# Returns
   - `tH::RV`: Hubble time in ``{\\rm \\mathbf{Gyr}}``.
"""
function hubble_time(H0::RV)
   return DIST_MPC / H0 / 1.0E3 / YEAR2SECOND / 1.0E9
end


"""
    age(cosmology::AbstractCosmology, z::RV) --> RV
Calculate age ``(t_{\\rm age})`` of the Universe at redshift ``z`` in ``{\\rm \\mathbf{Gyr}}``,
```math 
t_{\\rm age}(z) = \\int_z^\\infty \\frac{1}{(1+z')H(z')} dz'.
```

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z::RV`: Redshift.

# Returns
   - `age::RV`: Age of the Universe in ``{\\rm \\mathbf{Gyr}}``.
"""
function age(cosmology::AbstractCosmology, z::RV)
   # Hubble time (in years)
   tH = hubble_time(cosmology.H0)

   # Integral part
   function integrand(cosmology, z)
      return 1.0 / ( (1.0 + z) * Ez(cosmology, z) )
   end
   tC, _ = quadgk(x -> integrand(cosmology, x), z, Inf, atol=1E-10)
   return tH * tC
end


"""
    lookback_time(cosmology::AbstractCosmology, z::RV) --> RV
Calculate lookbak time ``(t_L)`` to a given redshift ``z`` in ``{\\rm \\mathbf{Gyr}}``,
```math
t_L(z) = \\int_0^z \\frac{1}{(1+z')H(z')} dz'.
```

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z::RV`: Redshift.

# Returns
   - `tL::RV`: Lookback time in ``{\\rm \\mathbf{Gyr}}``.
"""
function lookback_time(cosmology::AbstractCosmology, z::RV)
   # Hubble time (in years)
   tH = hubble_time(cosmology.H0)

   # Integral part
   function integrand(cosmology, z)
      return 1.0 / ( (1.0 + z) * Ez(cosmology, z) )
   end
   tC, _ = quadgk(x -> integrand(cosmology, x), 0, z, atol=1E-10)
   return tH * tC
end


"""
    rho_cz(cosmology::AbstractCosmology, z::RV) --> RV
Calculate the critical density ``(\\rho_c)`` of the Universe at redshift ``z`` in ``{\\rm \\mathbf{kg/m^3}}``,
```math
\\rho_c(z) = \\frac{3 H^2(z)}{8 π {\\rm G} }.
```

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z::RV`: Redshift.

# Returns
   - `rho_c::RV`: Critical density in ``{\\rm \\mathbf{kg/m^3}}``.
"""
function rho_cz(cosmology::AbstractCosmology, z::RV)
   return (3.0 / 8.0 / π / CONST_G) * hubble_parameter(cosmology, z)^2 * (1.0E3 / DIST_MPC)^2
end


"""
    Omega_mz(cosmology::AbstractCosmology, z::RV) --> RV
Calculate the dimensionless matter density parameter ``(\\Omega_{m})`` at redshift ``z``,
```math
Ω_{m}(z) = Ω_{m}(0) \\: (1+z)^3 \\left( \\frac{H0}{H(z)} \\right)^2.
```

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z::RV`: Redshift.

# Returns
   - `Omega_m::RV`: Dimensionless matter density parameter.
"""
function Omega_mz(cosmology::AbstractCosmology, z::RV)
   return cosmology.Omega_m0 * (1 + z)^3 * (cosmology.H0 / hubble_parameter(cosmology, z))^2
end


"""
    Omega_rz(cosmology::AbstractCosmology, z::RV) --> RV
Calculate the dimensionless radiation density parameter ``(\\Omega_{r})`` at redshift ``z``,
```math
Ω_{r}(z) = Ω_{r}(0) \\: (1+z)^4 \\left( \\frac{H_0}{H(z)} \\right)^2.
```

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z::RV`: Redshift.

# Returns
   - `Omega_r::RV`: Dimensionless radiation density parameter.
"""
function Omega_rz(cosmology::AbstractCosmology, z::RV)
   return cosmology.Omega_r0 * (1 + z)^4 * (cosmology.H0 / hubble_parameter(cosmology, z))^2
end 


"""
    Omega_wz(cosmology::AbstractCosmology, z::RV) --> RV
Calculate the dimensionless dark energy density parameter ``(\\Omega_{w})`` at redshift ``z``,
```math
Ω_{w}(z) = Ω_{w0} (1+z)^{3(1+w)} \\left( \\frac{H_0}{H(z)} \\right)^2.
```

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z::RV`: Redshift.

# Returns
   - `Omega_w::RV`: Dimensionless dark energy density parameter.
"""
function Omega_wz(cosmology::AbstractCosmology, z::RV)
   return cosmology.Omega_w0 * (1 + z)^(3.0 * (1.0 + cosmology.w)) * (cosmology.H0/hubble_parameter(cosmology, z))^2
end


"""
    Omega_kz(cosmology::AbstractCosmology, z::RV) --> RV
Calculate the dimensionless curvature density parameter ``(\\Omega_{k})`` at redshift ``z``,
```math 
Ω_{k}(z) = Ω_{k0} (1+z)^2 \\left( \\frac{H_0}{H(z)} \\right)^2.
```

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z::RV`: Redshift.

# Returns
   - `Omega_k::RV`: Dimensionless curvature density parameter.
"""
function Omega_kz(cosmology::AbstractCosmology, z::RV)
   return cosmology.Omega_k0 * (1 + z)^2 * (cosmology.H0/hubble_parameter(cosmology, z))^2
end


""" 
    hubble_distance(H0::RV) --> RV
Calculate the Hubble distance (i.e., size of the observable Universe) ``(D_H)`` in ``{\\rm \\mathbf{meters}}``,
```math 
D_H = \\frac{\\rm c}{\\rm H_0}.
```

# Arguments
   - `H0::RV`: Hubble constant in ``{\\rm \\mathbf{km/s/Mpc}}``.

# Returns
   - `D_H::RV`: Hubble distance in ``{\\rm \\mathbf{meters}}``.
"""
function hubble_distance(H0::RV)
   return CONST_C * DIST_MPC / H0 / 1.0E3
end


"""
    comoving_distance_radial(cosmo::AbstractCosmology, z1::RV, z2::RV) --> RV
Calculate the comoving radial distance ``(D_C)`` between ``z_1`` and ``z_2`` in ``{\\rm \\mathbf{meters}}``.
The formula is,
```math
D_C = D_H \\int_{z_1}^{z_2} \\frac{dz'}{E(z')}.
```

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z1::RV`: Redshift.
   - `z2::RV`: Redshift.

# Returns
   - `D_C::RV`: Comoving radial distance in ``{\\rm \\mathbf{meters}}``.
"""
function comoving_distance_radial(cosmology::AbstractCosmology, z1::RV, z2::RV)
   # Hubble distance
   dH = hubble_distance(cosmology.H0)

   # Integral part
   function integrand(cosmology, z)
      return 1.0 / Ez(cosmology, z)
   end
   dC, _ = quadgk(x -> integrand(cosmology, x), z1, z2, atol=1E-10)
   return dH * dC
end


"""
    comoving_distance_transverse(cosmo::AbstractCosmology, z1::RV, z2::RV) --> RV
Calculate the comoving radial distance ``(D_M)`` between, ``z_1`` and ``z_2`` in ``{\\rm \\mathbf{meters}}``,
```math
D_M = \\begin{cases} 
D_H \\frac{1}{\\sqrt{\\Omega_k}}   \\sinh \\left[ \\sqrt{\\Omega_k} \\:   D_C / D_H \\right], & \\text{if } \\Omega_k > 0, \\\\
D_C,                                                                                          & \\text{if } \\Omega_k = 0, \\\\
D_H \\frac{1}{\\sqrt{|\\Omega_k|}}  \\sin \\left[ \\sqrt{|\\Omega_k|} \\: D_C / D_H \\right], & \\text{if } \\Omega_k < 0. \\\\
\\end{cases}
```

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z1::RV`: Redshift.
   - `z2::RV`: Redshift.

# Returns
   - `D_M::RV`: Comoving radial distance in ``{\\rm \\mathbf{meters}}``.
"""
function comoving_distance_transverse(cosmology::AbstractCosmology, z1::RV, z2::RV)
   # Get the Hubble distance
   dH = hubble_distance(cosmology.H0)

   # Get the comoving distance (radial)
   dC = comoving_distance_radial(cosmology, z1, z2)

   # Get the comoving distance (transverse)
   dM = 0.0
   if cosmology.Omega_k0 != 0
      Ωχ = sqrt(abs(cosmology.Omega_k0)) * dC/dH
      if cosmology.Omega_k0 > 0
         dM = dH * sinh(Ωχ) / sqrt(abs(cosmology.Omega_k0))
      else
         dM =  dH * sin(Ωχ) / sqrt(abs(cosmology.Omega_k0))
      end
   else
      dM = dC
   end
   return dM
end


"""
    luminosity_distance(cosmology::AbstractCosmology, z::RV) --> RV
Calculate the luminosity distance ``(D_L)`` to redshift ``z`` in ``{\\rm \\mathbf{meters}}``,
```math
D_L(z) = (1+z) D_M(z).
```

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z::RV`: Redshift.

# Returns
   - `D_L::RV`: Luminosity distance in ``{\\rm \\mathbf{meters}}``.
"""
function luminosity_distance(cosmology::AbstractCosmology, z::RV)
   return (1.0 + z) * comoving_distance_transverse(cosmology, 0.0, z)
end


"""
    angular_diameter_distance(cosmology::AbstractCosmology, z1::RV, z2::RV) --> RV
Calculate the angular diameter distance ``(D_A)`` between redshifts ``z_1`` and ``z_2`` in ``{\\rm \\mathbf{meters}}``,
```math
D_A(z_1, z_2) = \\frac{D_M(z_1, z_2)}{1+z_2}.
```

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z1::RV`: Redshift.
   - `z2::RV`: Redshift.

# Returns
   - `D_A::RV`: Angular diameter distance in ``{\\rm \\mathbf{meters}}``.
"""
function angular_diameter_distance(cosmology::AbstractCosmology, z1::RV, z2::RV)
   return comoving_distance_transverse(cosmology, z1, z2) / (1.0 + z2)
end


"""
    distance_modulus(cosmo::AbstractCosmology, z::RV) --> RV
Calculate the distance modulus ``(\\mu)`` to a given redshift ``z``,
```math
\\mu = 5 \\log\\left( \\frac{D_L}{\\rm pc} \\right) - 5.
```

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z::RV`: Redshift.

# Returns
   - `mu::RV`: Distance modulus.
"""
function distance_modulus(cosmology::AbstractCosmology, z::RV)
   return 5.0 * log10(luminosity_distance(cosmology, z) / DIST_PC) - 5.0
end


"""
    angular_scale(cosmology::AbstractCosmology, z::RV) --> RV
Calculate the angular size in ``{\\rm \\mathbf{Kpc}}`` for ``1''`` on sky at redhsift, ``z``,
```math
d = 1'' \\times \\left( \\frac{D_A}{\\rm{kpc}} \\right).
```

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z::RV`: Redshift.

# Returns
   - `d::RV`: Angular size in ``{\\rm \\mathbf{Kpc}}``.
"""
function angular_scale(cosmology::AbstractCosmology, z::RV)
   return (angular_diameter_distance(cosmology, 0.0, z) / DIST_KPC) * ANGLE_ARCSEC
end


"""
    comoving_volume_element(cosmology::AbstractCosmology, z::RV) --> RV
Calculate the comving volume element ``(dV_C)`` at redshift ``z`` in ``{\\rm \\mathbf{Gpc^3}}``,
```math
dV_C = D_H \\frac{D_M^2(z)}{E(z)}.
```

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z::RV`: Redshift.

# Returns
   - `dV_C::RV`: Comving volume element in ``{\\rm \\mathbf{Gpc^3}}``.
"""
function comoving_volume_element(cosmology::AbstractCosmology, z::RV)
   # Get the Hubble distance
   dH = hubble_distance(cosmology.H0)

   # E(z) factor
   ez_value = Ez(cosmology, z)

   # Get the transverse comoving distance
   dC_trans = comoving_distance_transverse(cosmology, 0.0, z)

   return dH * dC_trans^2 / ez_value / DIST_GPC^3
end


"""
    comoving_volume(cosmology::AbstractCosmology, z::RV) --> RV
Calculate the total comving volume up to redshift ``z`` in ``{\\rm \\mathbf{Gpc^3}}``,
```math
V_C = \\begin{cases} 
\\frac{4π}{2} \\frac{D_H^3}{Ω_k} 
   \\left[ \\frac{D_M}{D_H} \\sqrt{1+Ω_k \\left(\\frac{D_M}{D_H}\\right)^2} 
   - \\frac{1}{\\sqrt{|Ω_k|}} {\\rm arcsinh}\\left( \\sqrt{|Ω_k|} \\frac{D_M}{D_H}  \\right) \\right], & \\text{if } Ω_k > 0, \\\\
\\frac{4π}{3} D_M^3,                                                  & \\text{if } Ω_k = 0, \\\\
\\frac{4π}{2} \\frac{D_H^3}{|Ω_k|} \\left[ \\frac{D_M}{D_H} \\sqrt{1+Ω_k \\left(\\frac{D_M}{D_H}\\right)^2} 
      - \\frac{1}{\\sqrt{|Ω_k|}} \\arcsin\\left( \\sqrt{|Ω_k|} \\frac{D_M}{D_H}  \\right) \\right], & \\text{if } Ω_k < 0. \\\\
\\end{cases}
```

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z::RV`: Redshift.

# Returns
   - `com_vol::RV`: Comving volume up to redshift ``z`` in ``{\\rm \\mathbf{Gpc^3}}``.
"""
function comoving_volume(cosmology::AbstractCosmology, z::RV)
   # Get the Hubble distance
   dH = hubble_distance(cosmology.H0)

   # Get the comoving distance (radial)
   dM = comoving_distance_transverse(cosmology, 0.0, z)

   # Get the comoving volume up to redshift z
   com_vol = 0.0
   if cosmology.Omega_k0 != 0.0
      term1 = 4.0 * π * dH^3 / 2.0 / cosmology.Omega_k0
      term2 = (dM/dH) * sqrt(1.0 + cosmology.Omega_k0 * (dM / dH)^2)
      if cosmology.Omega_k0 > 0.0
         com_vol = term1 * (term2 - asinh(sqrt(abs(cosmology.Omega_k0)) * dM/dH ) / sqrt(abs(cosmology.Omega_k0)))
      else
         com_vol = term1 * (term2 -  asin(sqrt(abs(cosmology.Omega_k0)) * dM/dH ) / sqrt(abs(cosmology.Omega_k0)))
      end
   else
      com_vol = (4.0 * π / 3) * dM^3
   end
   return com_vol / DIST_GPC^3
end


"""
    adis2zs(cosmology::AbstractCosmology, z_d::RV, adis::RV; max_iter::Int64=10000, tol::Float64=1e-6) --> RV
Calculate the source redshift (``z_s``) from the distance ratio (``a_{\\rm dis}``), using Bi-section method.

# Arguments
   - `cosmology::AbstractCosmology`: Cosmology object.
   - `z_d::RV`: Lens redshift.
   - `adis::RV`: Distance ratio, ``a_{\\rm dis}``.
   - `max_iter::Int64=10000`: Maximum number of iterations.
   - `tol::Float64=1e-6`: Tolerance for the root finding algorithm.

# Returns
   - `z_s::RV`: Redshift of the source.
"""
function adis2zs(cosmology::AbstractCosmology, z_d::RV, adis::RV; max_iter::Int64=100_000, tol::Float64=1e-10)
   # Source low and high redshift bounds
   z_l = z_d + 1E-6
   z_u = 100.0

   # Get distance ratio between [0, 1]
   adis_a = angular_diameter_distance(cosmology, z_d, z_l) / angular_diameter_distance(cosmology, 0.0, z_l)
   adis_b = angular_diameter_distance(cosmology, z_d, z_u) / angular_diameter_distance(cosmology, 0.0, z_u)

   for _ in 1:max_iter
      # Get the middle redshift
      z_m = (z_l + z_u) / 2.0

      # Get the middle redshift angular diameter distance
      adis_m = angular_diameter_distance(cosmology, z_d, z_m) / angular_diameter_distance(cosmology, 0.0, z_m)

      # Check if the middle redsfhit adis is close to the target adis
      if abs(adis_m - adis) < tol
         return z_m
      end

      # Update the low and high redshift bounds
      if adis_m < adis
         z_l = z_m
      else
         z_u = z_m
      end
   end
   
   error("Failed to converge after $max_iter iterations.")
   return nothing
end

end