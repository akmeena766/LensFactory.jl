Arguably the simplest lens model is a point mass (Schwarzschild) lens, which is described by two
parameters, namely, mass (``M``) and source position (``y``). The ``(\psi, \hat{\alpha}, \kappa)``
for a point mass lens are given as,
```math
\begin{align*}
\psi(\theta) &= \theta_E^2 \ln(\theta), \\
\hat{\alpha}(\theta) &= \frac{\theta_E^2}{\theta}, \\
\kappa(\theta) &= \frac{M}{\Sigma_{\rm cr}} \delta(\theta),
\end{align*}
```
where ``\theta_E`` is the Einstein angle, which is given as,
```math
\begin{equation*}
\theta_E = \sqrt{\frac{4{\rm G}M}{\rm c^2} \frac{D_{ds}}{D_d D_s}}.
\end{equation*}
```

```@docs
Lenses.init_PointLens
Lenses.PointLens.potential!
Lenses.PointLens.deflection!
Lenses.PointLens.jacobian!
Lenses.PointLens.einstein_angle
```