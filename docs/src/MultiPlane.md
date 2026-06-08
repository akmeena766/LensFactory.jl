# MultiPlane

Following same terminology as single plane lensing, we can write the lensing potential and 
deflection in ``i``-th plane as,
```math
\begin{align*}
ψ_i(\pmb{θ}_i)       &= \frac{4{\rm G}}{\rm c^2} \frac{1}{D_i} \int d^2 \pmb{θ}' \, Σ(\pmb{θ}') \, \ln|\pmb{θ}_i - \pmb{θ'}|, \\
\pmb{α}_i(\pmb{θ}_i) &= \pmb{∇}_i ψ_i(\pmb{θ}_i),
\end{align*}
```
where ``D_i`` is the angular diameter distance from observer to the ``i``-th lens plane. Starting 
from observer, by the time a ray reaches ``i``-th lens plane it has already deflected in previous 
``i-1`` lens planes. Hence, the impact parameter in the ``i``-th plane is given by,
```math
\pmb{θ}_i = \pmb{θ}_1 - \sum_{j=1}^{i-1} \pmb{α}_j(\pmb{θ}_j).
```



```@docs
MultiPlane.get_potential
MultiPlane.get_deflection
MultiPlane.get_jacobian
MultiPlane.get_time_delay
MultiPlane.get_magnification_image
MultiPlane.get_magnification_source
MultiPlane.get_image
MultiPlane.get_critical_curve
MultiPlane.get_caustic
```