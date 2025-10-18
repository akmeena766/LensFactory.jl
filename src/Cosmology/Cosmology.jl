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


# Initialize an abstract type for cosmology
abstract type AbstractCosmology end


"""
    init_cosmology(H0::RV=70.0, w::RV=-1.0, Omega_m0::RV=0.3, Omega_r0::RV=0.0, Omega_w0::RV=0.7)

An Abstract type to initialize a given cosmology.
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
    scale_factor(z::RV)

Calculates the scale factor ``(a)`` at redshift ``z``. The formula is,

```math
a = \\frac{1}{1+z}
```
"""
function scale_factor(z::RV)
   return 1.0 / (1.0 + z)
end


"""
    Ez(z::RV)

Calculates the dimensionless Hubble parameter ``(E)`` at redshift, ``z``. 
The formula is,
```math
E(z) = \\sqrt{ Ω_{m0} (1+z)^2 + Ω_{r0} (1+z)^4 + Ω_{k0} (1+z)^2 + Ω_{w0} (1+z)^{3(1+w)} }
```
"""
function Ez(cosmology::AbstractCosmology, z::RV)
      return sqrt(cosmology.Omega_m0 * (1.0 + z)^3 +  
                  cosmology.Omega_r0 * (1.0 + z)^4 + 
                  cosmology.Omega_k0 * (1.0 + z)^2 + 
                  cosmology.Omega_w0 * (1.0 + z)^(3.0 * (1.0 + cosmology.w)))
end


"""
    hubble_parameter(cosmology::AbstractCosmology, z::RV)

Calculates the Hubble parameter ``(H)`` at redshift ``z`` in ``{\\rm \\mathbf{km/s/Mpc}}``. 
The formula is,
```math
H(z) = H_0 \\: E(z)
```
"""
function hubble_parameter(cosmology::AbstractCosmology, z::RV)
   return cosmology.H0 * Ez(cosmology, z)
end


""" 
    hubble_time(H0::RV)

Calculates the Hubble time ``(t_H)`` in ``{\\rm \\mathbf{Gyr}}``. The formula is,
```math
t_H = \\frac{1}{H_0}
```
"""
function hubble_time(H0::RV)
   return DIST_MPC / H0 / 1.0E3 / YEAR2SECOND / 1.0E9
end


"""
    age(cosmology::AbstractCosmology, z::RV)

Calculates the age ``(t_{\\rm age})`` of the Universe at redshift ``z`` in ``{\\rm \\mathbf{Gyr}}``.
The formula is,
```math 
t_{\\rm age}(z) = \\int_z^\\infty \\frac{1}{(1+z')H(z')} dz' 
```
"""
function age(cosmology::AbstractCosmology, z::RV)
   # Hubble time (in years)
   tH = hubble_time(cosmology.H0)

   # Integral part
   function integrand(cosmology, z)
      return 1.0 / ( (1.0 + z) * Ez(cosmology, z) )
   end
   tC, tC_err = quadgk(x -> integrand(cosmology, x), z, Inf)  
   
   # Throw error if integration error is large
   if tC_err > 1E-8
      throw(ArgumentError("Error in lookback time estimation is $tC_err."))
   end
   return tH * tC
end


"""
    lookback_time(cosmology::AbstractCosmology, z::RV)

Calculates lookbak time ``(t_L)`` to a given redshift ``z`` in ``{\\rm \\mathbf{Gyr}}``.
The formula is,
```math 
t_L(z) = \\int_0^z \\frac{1}{(1+z')H(z')} dz' 
```
"""
function lookback_time(cosmology::AbstractCosmology, z::RV)
   # Hubble time (in years)
   tH = hubble_time(cosmology.H0)

   # Integral part
   function integrand(cosmology, z)
      return 1.0 / ( (1.0 + z) * Ez(cosmology, z) )
   end
   tC, tC_err = quadgk(x -> integrand(cosmology, x), 0, z)

   # Throw error if integration error is large
   if tC_err > 1E-8
      throw(ArgumentError("Error in lookback time estimation is $tC_err."))
   end
   return tH * tC
end


"""
    rho_cz(cosmology::AbstractCosmology, z::RV)

Calculates the critical density ``(\\rho_c)`` of the Universe at redshift ``z`` in ``{\\rm \\mathbf{kg/m^3}}``.
The formula is,
```math
\\rho_c(z) = \\frac{3 H^2(z)}{8 π {\\rm G} }
```
"""
function rho_cz(cosmology::AbstractCosmology, z::RV)
   return (3.0 / 8.0 / π / CONST_G) * hubble_parameter(cosmology, z)^2 * (1.0E3 / DIST_MPC)^2
end


"""
    Omega_mz(cosmology::AbstractCosmology, z::RV)

Calculates the matter density parameter ``(\\Omega_{m})`` at redshift ``z``. The formula is,

```math
Ω_{m}(z) = Ω_{m}(0) \\: (1+z)^3 \\left( \\frac{H0}{H(z)} \\right)^2
```
"""
function Omega_mz(cosmology::AbstractCosmology, z::RV)
   return cosmology.Omega_m0 * (1 + z)^3 * (cosmology.H0 / hubble_parameter(cosmology, z))^2
end


"""
    Omega_rz(cosmology::AbstractCosmology, z::RV)

Calculates the matter density parameter ``(\\Omega_{r})`` at redshift ``z``. The formula is,

```math
Ω_{r}(z) = Ω_{r}(0) \\: (1+z)^4 \\left( \\frac{H_0}{H(z)} \\right)^2
```
"""
function Omega_rz(cosmology::AbstractCosmology, z::RV)
   return cosmology.Omega_r0 * (1 + z)^4 * (cosmology.H0 / hubble_parameter(cosmology, z))^2
end 


"""
    Omega_wz(cosmology::AbstractCosmology, z::RV)

Calculates the dark energy density parameter ``(\\Omega_{w})`` at redshift ``z``. The formula is,

```math
Ω_{w}(z) = Ω_{w0} (1+z)^{3(1+w)} \\left( \\frac{H_0}{H(z)} \\right)^2
```
"""
function Omega_wz(cosmology::AbstractCosmology, z::RV)
   return cosmology.Omega_w0 * (1 + z)^(3.0 * (1.0 + cosmology.w)) * (cosmology.H0/hubble_parameter(cosmology, z))^2
end


"""
    Omega_kz(cosmology::AbstractCosmology, z::RV)

Calculates the curvature density parameter ``(\\Omega_{k})`` at redshift ``z``. The formula is,

```math 
Ω_{k}(z) = Ω_{k0} (1+z)^2 \\left( \\frac{H_0}{H(z)} \\right)^2
```
"""
function Omega_kz(cosmology::AbstractCosmology, z::RV)
   return cosmology.Omega_k0 * (1 + z)^2 * (cosmology.H0/hubble_parameter(cosmology, z))^2
end


""" 
    hubble_distance(H0::RV)

Calculates HUbble distance (i.e., size of the observable Universe) ``(D_H)`` in ``{\\rm \\mathbf{meters}}``. 
The formula is

```math 
D_H = \\frac{\\rm c}{\\rm H_0} 
```
"""
function hubble_distance(H0::RV)
   return CONST_C * DIST_MPC / H0 / 1.0E3
end


"""
    comoving_distance_radial(cosmo::AbstractCosmology, z1::RV, z2::RV)

Calculates the comoving radial distance ``(D_C)`` between ``z_1`` and ``z_2`` in ``{\\rm \\mathbf{meters}}``.
The formula is,
```math
D_C = D_H \\int_{z_1}^{z_2} \\frac{dz'}{E(z')}
```
"""
function comoving_distance_radial(cosmology::AbstractCosmology, z1::RV, z2::RV)
   # Hubble distance
   dH = hubble_distance(cosmology.H0)

   # Integral part
   function integrand(cosmology, z)
      return 1.0 / Ez(cosmology, z)
   end
   dC, dC_err = quadgk(x -> integrand(cosmology, x), z1, z2)

   # Throw error if integration error is large
   if dC_err > 1E-8
      throw(ArgumentError("Error in comoving distance estimation is $dC_err."))
   end
   return dH * dC
end


"""
    comoving_distance_transverse(cosmo::AbstractCosmology, z1::RV, z2::RV)

Calculates the comoving radial distance ``(D_M)`` between, ``z_1`` and ``z_2`` in ``{\\rm \\mathbf{meters}}``.
The formula is,
```math
D_M = \\begin{cases} 
D_H \\frac{1}{\\sqrt{\\Omega_k}}   \\sinh \\left[ \\sqrt{\\Omega_k} \\:   D_C / D_H \\right] & \\text{if } \\Omega_k > 0 \\\\
D_C                                                                                          & \\text{if } \\Omega_k = 0 \\\\
D_H \\frac{1}{\\sqrt{|\\Omega_k|}}  \\sin \\left[ \\sqrt{|\\Omega_k|} \\: D_C / D_H \\right] & \\text{if } \\Omega_k < 0 \\\\
\\end{cases}
```
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
    luminosity_distance(cosmology::AbstractCosmology, z::RV)
Calculates the luminosity distance ``(D_L)`` to redshift ``z`` in ``{\\rm \\mathbf{meters}}``. 
The formula is,
```math
D_L(z) = (1+z) D_M(z)
```
"""
function luminosity_distance(cosmology::AbstractCosmology, z::RV)
   return (1.0 + z) * comoving_distance_transverse(cosmology, 0.0, z)
end


"""
    angular_diameter_distance(cosmology::AbstractCosmology, z1::RV, z2::RV)
Calculates the angular diameter distance ``(D_A)`` between redshifts ``z_1`` and ``z_2`` in ``{\\rm \\mathbf{meters}}``. 
The formula is,
```math
D_A(z_1, z_2) = \\frac{D_M(z_1, z_2)}{1+z_2}
```
"""
function angular_diameter_distance(cosmology::AbstractCosmology, z1::RV, z2::RV)
   return comoving_distance_transverse(cosmology, z1, z2) / (1.0 + z2)
end


"""
    distance_modulus(cosmo::AbstractCosmology, z::RV)

Calculates the distance modulus ``(\\mu)`` to a given redshift ``z``. 
The formula is,
```math
\\mu = 5 \\log\\left( \\frac{D_L}{\\rm pc} \\right) - 5
```
"""
function distance_modulus(cosmology::AbstractCosmology, z::RV)
   return 5.0 * log10(luminosity_distance(cosmology, z) / DIST_PC) - 5.0
end


"""
    angular_scale(cosmo::AbstractCosmology, z::RV)

Calculates the angular size in ``{\\rm \\mathbf{Kpc}}`` for ``1''`` on sky at redhsift, ``z``. 
The formula is
```math
d = 1'' \\times \\left( \\frac{D_A}{\\rm{kpc}} \\right)
```
"""
function angular_scale(cosmology::AbstractCosmology, z::RV)
   return (angular_diameter_distance(cosmology, 0.0, z) / DIST_KPC) * ANGLE_ARCSEC
end


"""
    comoving_volume_element(cosmology::AbstractCosmology, z::RV)

Calculates the comving volume element ``(dV_C)`` at redshift ``z`` in ``{\\rm \\mathbf{Gpc^3}}``.
The formula is
```math
dV_C = D_H \\frac{D_M^2(z)}{E(z)}
```
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
    comoving_volume(cosmology::AbstractCosmology, z::RV)

Calculates the total comving volume up to redshift ``z`` in ``{\\rm \\mathbf{Gpc^3}}``.
The formula is,
```math
V_C = \\begin{cases} 
\\frac{4π}{2} \\frac{D_H^3}{Ω_k} 
   \\left[ \\frac{D_M}{D_H} \\sqrt{1+Ω_k \\left(\\frac{D_M}{D_H}\\right)^2} 
   - \\frac{1}{\\sqrt{|Ω_k|}} {\\rm arcsinh}\\left( \\sqrt{|Ω_k|} \\frac{D_M}{D_H}  \\right) \\right] & \\text{if } Ω_k > 0 \\\\
\\frac{4π}{3} D_M^3                                                  & \\text{if } Ω_k = 0 \\\\
\\frac{4π}{2} \\frac{D_H^3}{|Ω_k|} \\left[ \\frac{D_M}{D_H} \\sqrt{1+Ω_k \\left(\\frac{D_M}{D_H}\\right)^2} 
      - \\frac{1}{\\sqrt{|Ω_k|}} \\arcsin\\left( \\sqrt{|Ω_k|} \\frac{D_M}{D_H}  \\right) \\right] & \\text{if } Ω_k < 0 \\\\
\\end{cases}
```
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

end