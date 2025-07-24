"""
    Cosmology

Testing :-)

"""
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
    scale_factor(z::RV)::RV

Calculates the scale factor ``(a)`` at redshift ``z``. The formula is,

```math
a = \\frac{1}{1+z}
```
"""
function scale_factor(z::RV)::RV
   return 1.0 / (1.0 + z)
end


"""
    Ez(z::RV)::RV

Calculates the dimensionless Hubble parameter ``(E)`` at redshift, ``z``. 
The formula is,
```math
E(z) = \\sqrt{ \\Omega_{m0} (1+z)^2 + \\Omega_{r0} (1+z)^4 + \\Omega_{k0} (1+z)^2 + \\Omega_{w0} (1+z)^{3(1+w)} }
```
"""
function Ez(cosmo::AbstractCosmology, z::RV)::RV
      return sqrt(cosmo.Omega_m0 * (1.0 + z)^3 +  
                  cosmo.Omega_r0 * (1.0 + z)^4 + 
                  cosmo.Omega_k0 * (1.0 + z)^2 + 
                  cosmo.Omega_w0 * (1.0 + z)^(3.0 * (1.0 + cosmo.w)))
end


"""
    hubble_parameter(cosmology::AbstractCosmology, z::RV)::RV

Calculates the Hubble parameter ``(H)`` at redshift ``z`` in ``{\\rm \\mathbf{km/s/Mpc}}``. 
The formula is,
```math
H(z) = H_0\\:E(z)
```
"""
function hubble_parameter(cosmo::AbstractCosmology, z::RV)::RV
   return cosmo.H0 * Ez(cosmo, z)
end


""" 
    hubble_time(H0::RV)::RV

Calculates the Hubble time ``(t_H)`` in ``{\\rm \\mathbf{Gyr}}``. The formula is,
```math
t_H = \\frac{1}{\\rm H_0}
```
"""
function hubble_time(H0::RV)::RV
   return DIST_MPC / H0 / 1.0E3 / YEAR2SECOND / 1.0E9
end


"""
    age(cosmology::AbstractCosmology, z::RV)::RV

Calculates the age ``(t_{\\rm age})`` of the Universe at redshift ``z`` in ``{\\rm \\mathbf{Gyr}}``.
The formula is,
```math 
t_{\\rm age}(z) = \\int_z^\\infty \\frac{1}{(1+z')H(z')} dz' 
```
"""
function age(cosmo::AbstractCosmology, z::RV)::RV
   # Hubble time (in years)
   tH::RV = hubble_time(cosmo.H0)

   # Integral part
   function integrand(cosmo, z)
      return 1.0 / ( (1.0 + z) * Ez(cosmo, z) )
   end
   tC, tC_err = quadgk(x -> integrand(cosmo, x), z, Inf)  
   
   # Throw error if integration error is large
   if tC_err > 1E-8
      throw(ArgumentError("Error in lookback time estimation is $tC_err."))
   end
   
   return tH * tC
end


"""
    lookback_time(cosmology::AbstractCosmology, z::RV)::RV

Calculates lookbak time ``(t_L)`` to a given redshift ``z`` in ``{\\rm \\mathbf{Gyr}}``.
The formula is,
```math 
t_L(z) = \\int_0^z \\frac{1}{(1+z')H(z')} dz' 
```
"""
function lookback_time(cosmo::AbstractCosmology, z::RV)::RV
   # Hubble time (in years)
   tH::RV = hubble_time(cosmo.H0)

   # Integral part
   function integrand(cosmo, z)
      return 1.0 / ( (1.0 + z) * Ez(cosmo, z) )
   end
   tC, tC_err = quadgk(x -> integrand(cosmo, x), 0, z)

   # Throw error if integration error is large
   if tC_err > 1E-8
      throw(ArgumentError("Error in lookback time estimation is $tC_err."))
   end
   return tH * tC
end


"""
    rho_cz(cosmology::AbstractCosmology, z::RV)::RV

Calculates the critical density ``(\\rho_c)`` of the Universe at redshift ``z`` in ``{\\rm \\mathbf{kg/m^3}}``.
The formula is,
```math
\\rho_c(z) = \\frac{3 H^2(z)}{8 \\pi {\\rm G} }
```
"""
function rho_cz(cosmo::AbstractCosmology, z::RV)::RV
   return (3.0 / 8.0 / pi / CONST_G) * hubble_parameter(cosmo, z)^2 * (1.0E3 / DIST_MPC)^2
end


"""
    Omega_mz(cosmology::AbstractCosmology, z::RV)::RV

Calculates the matter density parameter ``(\\Omega_{m})`` at redshift `z`. The formula is,

```math
\\Omega_{m}(z) = \\Omega_{m}(0) \\: (1+z)^3 \\left( \\frac{H0}{H(z)} \\right)^2
```
"""
function Omega_mz(cosmo::AbstractCosmology, z::RV)::RV
   return cosmo.Omega_m0 * (1 + z)^3 * (cosmo.H0 / hubble_parameter(cosmo, z))^2
end


"""
    Omega_rz(cosmology::AbstractCosmology, z::RV)::RV

Calculates the matter density parameter ``(\\Omega_{r})`` at redshift `z`. The formula is,

```math 
\\Omega_{r}(z) = \\Omega_{r}(0) \\: (1+z)^4 \\left( \\frac{H_0}{H(z)} \\right)^2 
```
"""
function Omega_rz(cosmo::AbstractCosmology, z::RV)::RV
   return cosmo.Omega_r0 * (1 + z)^4 * (cosmo.H0 / hubble_parameter(cosmo, z))^2
end 


"""
    Omega_wz(cosmology::AbstractCosmology, z::RV)::RV

Calculates the dark energy density parameter ``(\\Omega_{w})`` at redshift `z`. The formula is,

```math 
\\Omega_{w}(z) = \\Omega_{w0} (1+z)^{3(1+w)} \\left( \\frac{H_0}{H(z)} \\right)^2 
```
"""
function Omega_wz(cosmo::AbstractCosmology, z::RV)::RV
   return cosmo.Omega_w0 * (1 + z)^(3.0 * (1.0 + cosmo.w)) * (cosmo.H0/hubble_parameter(cosmo, z))^2
end


"""
    Omega_kz(cosmology::AbstractCosmology, z::RV)::RV

Calculates the curvature density parameter ``(\\Omega_{k})`` at redshift `z`. The formula is,

```math 
\\Omega_{k}(z) = \\Omega_{k0} (1+z)^2 \\left( \\frac{H_0}{H(z)} \\right)^2 
```
"""
function Omega_kz(cosmo::AbstractCosmology, z::RV)::RV
   return cosmo.Omega_k0 * (1 + z)^2 * (cosmo.H0/hubble_parameter(cosmo, z))^2
end


""" 
    hubble_distance(H0::RV)::RV

Calculates HUbble distance (i.e., size of the observable Universe) ``(D_H)`` in ``{\\rm \\mathbf{meters}}``. 
The formula is

```math 
D_H = \\frac{\\rm c}{\\rm H_0} 
```
"""
function hubble_distance(H0::RV)::RV
   return CONST_C * DIST_MPC / H0 / 1.0E3
end


"""
    comoving_distance_radial(cosmo::AbstractCosmology, zi::RV, zf::RV)::RV

Calculates the comoving radial distance ``(D_C)`` between ``z_1`` and ``z_2`` in ``{\\rm \\mathbf{meters}}``.
The formula is,
```math
D_C = D_H \\int_{z_1}^{z_2} \\frac{dz'}{E(z')}
```
"""
function comoving_distance_radial(cosmo::AbstractCosmology, z1::RV, z2::RV)::RV
   # Hubble distance
   dH::RV = hubble_distance(cosmo.H0)

   # Integral part
   function integrand(cosmo, z)
      return 1.0 / Ez(cosmo, z)
   end
   dC, dC_err = quadgk(x -> integrand(cosmo, x), z1, z2)

   # Throw error if integration error is large
   if dC_err > 1E-8
      throw(ArgumentError("Error in comoving distance estimation is $dC_err."))
   end
   dC = dH * dC
   return dC
end


"""
    comoving_distance_transverse(cosmo::AbstractCosmology, zi::RV, zf::RV)::RV

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
function comoving_distance_transverse(cosmo::AbstractCosmology, zi::RV, zf::RV)::RV
   # Get the Hubble distance
   dH::RV = hubble_distance(cosmo.H0)

   # Get the comoving distance (radial)
   dC::RV = comoving_distance_radial(cosmo, zi, zf)

   # Get the comoving distance (transverse)
   dM::RV = 0.0
   if cosmo.Omega_k0 ≠ 0
      Ωχ = √( abs(cosmo.Omega_k0) ) * dC/dH
      if cosmo.Omega_k0 > 0 
         dM = dH * sinh(Ωχ) / √( abs(cosmo.Omega_k0) )
      else
         dM =  dH * sin(Ωχ) / √( abs(cosmo.Omega_k0) )
      end
   else
      dM = dC
   end
   return dM
end


"""
    luminosity_distance(cosmology::AbstractCosmology, z::RV)::RV
Calculates the luminosity distance ``(D_L)`` to redshift ``z`` in ``{\\rm \\mathbf{meters}}``. 
The formula is,
```math
D_L(z) = (1+z) D_M(z)
```
"""
function luminosity_distance(cosmo::AbstractCosmology, zf::RV)::RV
   return (1.0 + zf) * comoving_distance_transverse(cosmo, 0.0, zf)
end


"""
    angular_diameter_distance(cosmology::AbstractCosmology, z1::RV, z2::RV)::RV
Calculates the angular diameter distance ``(D_A)`` between redshifts ``z_1`` and ``z_2`` in ``{\\rm \\mathbf{meters}}``. 
The formula is,
```math
D_A(z_1, z_2) = \\frac{D_M(z_1, z_2)}{1+z_2}
```
"""
function angular_diameter_distance(cosmo::AbstractCosmology, zi::RV, zf::RV)::RV
   return comoving_distance_transverse(cosmo, zi, zf) / (1.0 + zf)
end


"""
    distance_modulus(cosmo::AbstractCosmology, zf::RV)::RV

Calculates the distance modulus ``(\\mu)`` to a given redshift ``z``. 
The formula is,
```math
\\mu = 5 \\log\\left( \\frac{D_L}{\\rm pc} \\right) - 5
```
"""
function distance_modulus(cosmo::AbstractCosmology, zf::RV)::RV
   return 5.0 * log10(luminosity_distance(cosmo, zf) / DIST_PC) - 5.0
end


"""
    angular_scale(cosmo::AbstractCosmology, zf::RV)::RV

Calculates the angular size in ``{\\rm \\mathbf{Kpc}}`` for ``1''`` on sky at redhsift, ``z``. 
The formula is
```math
d = 1'' \\times \\left( \\frac{D_A}{\\rm{kpc}} \\right)
```
"""
function angular_scale(cosmo::AbstractCosmology, z::RV)::RV
   return (angular_diameter_distance(cosmo, 0.0, z) / DIST_KPC) * ANGLE_ARCSEC
end


"""
    comoving_volume_element(cosmology::AbstractCosmology, zf::RV)::RV

Calculates the comving volume element ``(dV_C)`` at redshift ``z`` in ``{\\rm \\mathbf{Gpc^3}}``.
The formula is
```math
dV_C = D_H \\frac{D_M^2(z)}{E(z)}
```
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


"""
    comoving_volume(cosmology::AbstractCosmology, zf::RV)::RV

Calculates the total comving volume up to redshift ``z`` in ``{\\rm \\mathbf{Gpc^3}}``.
The formula is,
```math
V_C = \\begin{cases} 
\\frac{4\\pi}{2} \\frac{D_H^3}{\\Omega_k}   \\left[ \\frac{D_M}{D_H} \\sqrt{1+\\Omega_k \\left(\\frac{D_M}{D_H}\\right)^2} 
      - \\frac{1}{\\sqrt{|\\Omega_k|}} {\\rm arcsinh}\\left( \\sqrt{|\\Omega_k|} \\frac{D_M}{D_H}  \\right) \\right] & \\text{if } \\Omega_k > 0 \\\\
\\frac{4\\pi}{3} D_M^3                                                  & \\text{if } \\Omega_k = 0 \\\\
\\frac{4\\pi}{2} \\frac{D_H^3}{|\\Omega_k|} \\left[ \\frac{D_M}{D_H} \\sqrt{1+\\Omega_k \\left(\\frac{D_M}{D_H}\\right)^2} 
      - \\frac{1}{\\sqrt{|\\Omega_k|}} \\arcsin\\left( \\sqrt{|\\Omega_k|} \\frac{D_M}{D_H}  \\right) \\right] & \\text{if } \\Omega_k < 0 \\\\
\\end{cases}
```
"""
function comoving_volume(cosmo::AbstractCosmology, zf::RV)::RV
   # Get the Hubble distance
   dH::RV = hubble_distance(cosmo.H0)

   # Get the comoving distance (radial)
   dM::RV = comoving_distance_transverse(cosmo, 0.0, zf)

   # Get the comoving volume up to redshift z
   com_vol::RV = 0.0
   if cosmo.Omega_k0 ≠ 0.0
      term1 = 4.0 * π * dH^3 / 2.0 / cosmo.Omega_k0
      term2 = (dM/dH) * √( 1.0 + cosmo.Omega_k0 * (dM / dH)^2 )
      if cosmo.Omega_k0 > 0.0
         com_vol = term1 * ( term2 - asinh(√(abs(cosmo.Omega_k0)) * dM/dH ) / √(abs(cosmo.Omega_k0) ) )
      else
         com_vol = term1 * ( term2 -  asin(√(abs(cosmo.Omega_k0)) * dM/dH ) / √(abs(cosmo.Omega_k0) ) )
      end
   else
      com_vol = (4.0 * π / 3) * dM^3
   end
   return com_vol / DIST_GPC^3
end


end