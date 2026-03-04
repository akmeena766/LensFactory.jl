# Lenses

The `Lenses` module is the interface between the user and various lens models. The user will call 
functions from this module to compute various lensing quantities. For examples on how to use this 
module and various functions within it, see 
[Example - 2](https://github.com/akmeena766/LensFactory_Examples/blob/main/Example-2%3A%20Basic%20strong%20lensing.ipynb)
and 
[Example - 3](https://github.com/akmeena766/LensFactory_Examples/blob/main/Example-3%3A%20Lensing%20quantities.ipynb).

!!! note
    Keeping in mind that most of the astrophysical scenario (primarily strong lensing by galaxies or 
    galaxy clusters), the Einstein angle is ~[0.1, 50] arcseconds, the default input coordinate 
    units are arcseconds. In addition, potential, deflection, and jacobian are calculate assuming
    that the source is at infinity, i.e., ``a_{\rm dis} = D_{ds} / D_s = 1``. Hence, various 
    follow-up calculations require the use to provide ``a_{\rm dis}`` as an input.   

```@docs
Lenses.get_meshgrid
Lenses.get_critical_density
Lenses.get_potential
Lenses.get_deflection
Lenses.get_jacobian
Lenses.get_time_delay
Lenses.get_magnification_image
Lenses.get_magnification_source
Lenses.get_image
Lenses.get_critical_curve
Lenses.get_caustic
Lenses.get_critical_area
Lenses.get_einstein_angle
Lenses.shear_cartesian2polar
Lenses.shear_polar2cartesian
Lenses.ellipticity_cartesian2polar
Lenses.ellipticity_polar2cartesian
```