# Lenses

The `Lenses` module is the interface between the user and various lens models. The user will call 
functions from this module to compute various lensing quantities. For examples on how to use this 
module and various functions within it, see 
[Basic-Example-2](https://github.com/akmeena766/LensFactory-Examples/blob/main/Basic/Example-2%3A%20Point%20mass%20lens.ipynb),
[Basic-Example-3](https://github.com/akmeena766/LensFactory-Examples/blob/main/Basic/Example-3%3A%20SIS%20lens%20model.ipynb)
and 
[Basic-Example-4](https://github.com/akmeena766/LensFactory-Examples/blob/main/Basic/Example-4%3A%20Lensing%20quantities.ipynb).

!!! note
    Keeping in mind that most of the astrophysical scenario (primarily strong lensing by galaxies or 
    galaxy clusters), the Einstein angle is ~[0.1, 50] arcseconds, the default input coordinate 
    units are arcseconds. This is not suitable for microlensing studies requiring micro- to 
    milli-arcsecond resolutions, and users should use the `Microlens` (yet to be implemented) module
    in such cases.
    
    In addition, potential, deflection, and jacobian are calculate assuming
    that the source is at infinity, i.e., ``a_{\rm dis} = D_{ds} / D_s = 1``. Hence, various 
    follow-up calculations require providing ``a_{\rm dis}`` as an input.   

```@docs
Lenses.get_meshgrid
Lenses.get_critical_density
Lenses.get_potential
Lenses.get_deflection
Lenses.get_jacobian
Lenses.get_time_delay
Lenses.get_kappa_gamma
Lenses.get_magnification_image
Lenses.get_magnification_source
Lenses.get_image
Lenses.get_critical_curve
Lenses.get_caustic
Lenses.get_critical_area
Lenses.get_einstein_angle
Lenses.get_image_multiplicity
Lenses.get_radial_profile
Lenses.get_mass_profile
Lenses.shear_cartesian2polar
Lenses.shear_polar2cartesian
Lenses.ellipticity_cartesian2polar
Lenses.ellipticity_polar2cartesian
Lenses.parameter_NFWLens
Lenses.parameter_gNFWLens
Lenses.parameter_EinastoLens
```