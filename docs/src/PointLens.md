# Point Lens
Arguably, the simplest gravitational lens is a point mass lens (i.e., Schwarzschild lens). 
The lensing by a given point mass lens is characterized by two parameters: its mass ($M$) and the 
source position ($\boldsymbol{\beta}$). The corresponding lens potential is given by
```math
\begin{equation*}
ψ(\pmb{θ}) = \frac{4{\rm G} M} {\rm c^2} \frac{1}{D_d} \ln \left|\pmb{θ} - \pmb{θ}_c\right|,
\end{equation*}
```
where $\pmb{θ}_c$ represents the point mass lens center. In `LensFactory`, to define a point lens, 
the user need to specify four parameters: its position in the image plane ($x_c,~y_c$), its mass 
($M$) and the angular diameter distance to the lens ($D_d$). By default, the lens is placed at the 
origin, i.e., ($x_c,~y_c$) = (0, 0). From above, to initialize a point lens, we need to first define 
a cosmology and calculate the angular diameter distance to the lens. We refer reader to 
[Basic: Example - 2](https://github.com/akmeena766/LensFactory-Examples/blob/main/Basic/Example-2%3A%20Point%20mass%20lens.ipynb)
for more details.


```@docs
Lenses.init_PointLens
Lenses.PointLens.potential!
Lenses.PointLens.deflection!
Lenses.PointLens.jacobian!
Lenses.PointLens.einstein_angle
```