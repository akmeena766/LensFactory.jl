# PIEP Lens
A simple way to have elliptical lens model is by introducing ellipicity in the lensing potential of
an axis symmetric lens model. The `PIEPLens` corresponds to the elliptical version of 
[NSISP Lens](NSISPLens.md), given as,

```math
\begin{equation*}
\psi(\pmb{θ}) = 4π \left( \frac{σ_v}{\rm c} \right)^2 
\sqrt{ θ_s^2 + (θ_x - θ_{xc})^2 + \frac{(θ_y - θ_{yc})^2}{q^2} },
\end{equation*}
```
where $q \coloneqq \frac{1 - \epsilon}{1 + \epsilon}$ is the axis-ratio and $\epsilon$ is the 
ellipticity.

```@docs
Lenses.init_PIEPLens
Lenses.PIEPLens.potential!
Lenses.PIEPLens.deflection!
Lenses.PIEPLens.jacobian!
```