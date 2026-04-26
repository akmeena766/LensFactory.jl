# Point Lens
Arguably, the simplest gravitational lens is a **point mass lens** (i.e., Schwarzschild lens). 
The lensing by a point mass lens is characterized by two parameters: its mass ($M$) and the source 
position ($\boldsymbol{\beta}$).

In `LensFactory`, to define a point lens, the user need to specify four parameters: its position in
the image plane ($x_c, y_c$), its mass ($M$) and the angular diameter distance to the lens ($D_d$).
By default, the lens is placed at the origin, i.e., ($x_c, y_c$) = (0,0). From above, to initialize 
a point lens, we need to first define a cosmology and calculate the angular diameter distance to the 
lens. We refer reader to [Example - 2](https://github.com/akmeena766/LensFactory-Examples/blob/main/Basic/Example-2%3A%20Basic%20strong%20lensing.ipynb) for more details.




```@docs
Lenses.init_PointLens
Lenses.PointLens.potential!
Lenses.PointLens.deflection!
Lenses.PointLens.jacobian!
Lenses.PointLens.einstein_angle
```