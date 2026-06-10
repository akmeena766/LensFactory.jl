# MultiPlane

Following same terminology as single plane lensing, we can write the lensing potential and 
deflection in ``i``-th plane as,
```math
\begin{align*}
ψ_i(\pmb{θ}_i)       &= \frac{4{\rm G}}{\rm c^2} \frac{1}{D_i} \int d^2 \pmb{θ}' \, Σ_i(\pmb{θ}') \, \ln|\pmb{θ}_i - \pmb{θ'}|, \\
\pmb{α}_i(\pmb{θ}_i) &= \pmb{∇}_i ψ_i(\pmb{θ}_i),
\end{align*}
```
where ``D_i`` is the angular diameter distance from observer to the ``i``-th lens plane and ``Σ_i`` 
is the surface mass density in the ``i``-th lens plane. Starting from observer, by the time a ray 
reaches ``i``-th lens plane it has already deflected in previous ``i-1`` lens planes. Hence, the 
impact parameter in the 1st lens plane is related to impact parameters in ``i``-th lens plane as,
```math
\pmb{θ}_i = \pmb{θ}_1 - \sum_{j=1}^{i-1} \frac{D_{ji}}{D_i} \pmb{α}_j(\pmb{θ}_j).
```

With the above, the total potential and deflection for a lens system with $N_p$ planes is given as,
```math
\begin{align*}
ψ(\pmb{θ})       &= \sum_{i=1}^{N_p} \frac{D_{is}}{D_s} ψ_i(\pmb{θ}_i), \\
\pmb{α}(\pmb{θ}) &= \sum_{i=1}^{N_p} \frac{D_{is}}{D_s} \pmb{α}_i(\pmb{θ}_i).
\end{align*}
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