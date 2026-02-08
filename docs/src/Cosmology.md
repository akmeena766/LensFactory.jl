# Cosmology
The `Cosmology` module controls the cosmological parameters and various derived quantities. We can 
initialize a cosmology using [init_cosmology](#LensFactory.Cosmology.init_cosmology). For examples 
on how to use this module and various functions within it, see 
[Example - 1](https://github.com/akmeena766/LensFactory_Examples/blob/main/Example-1%3A%20Constants%20and%20Cosmology.ipynb).

## Initialization
```@docs
Cosmology.init_cosmology
```

## Basics
```@docs
Cosmology.scale_factor
Cosmology.Ez
Cosmology.hubble_parameter
```

## Time
```@docs
Cosmology.hubble_time
Cosmology.age
Cosmology.lookback_time
```

## Density parameters
```@docs
Cosmology.rho_cz
Cosmology.Omega_mz
Cosmology.Omega_rz
Cosmology.Omega_wz
Cosmology.Omega_kz
```

## Distances
```@docs
Cosmology.hubble_distance
Cosmology.comoving_distance_radial
Cosmology.comoving_distance_transverse
Cosmology.luminosity_distance
Cosmology.angular_diameter_distance
Cosmology.distance_modulus
Cosmology.angular_scale
Cosmology.adis2zs
```

## Volumes
```@docs
Cosmology.comoving_volume_element
Cosmology.comoving_volume
```