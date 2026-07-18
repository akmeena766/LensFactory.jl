# MultiPlane
When the deflectors along the line of sight are located at different redshifts (e.g., a galaxy 
lens with line-of-sight perturbers, or a cluster with structures at multiple redshifts), the 
single-plane approximation breaks down, and the light rays need to be propagated through a series 
of lens planes.

Following same terminology as single plane lensing, we can write the lensing potential and 
deflection in ``j``th plane as,
```math
\begin{align*}
ψ_j(\pmb{θ}_j)       &= \frac{4{\rm G}}{\rm c^2} \frac{1}{D_j} \int d^2 \pmb{θ}' \, Σ_j(\pmb{θ}') \, \ln|\pmb{θ}_j - \pmb{θ'}|, \\
\pmb{α}_j(\pmb{θ}_j) &= \pmb{∇}_j ψ_j(\pmb{θ}_j),
\end{align*}
```
where ``D_j`` is the angular diameter distance from observer to the ``j``th lens plane and ``Σ_j`` 
is the surface mass density in the ``j``th lens plane. Starting from observer, by the time a ray 
reaches ``j``th lens plane it has already deflected in previous ``j-1`` lens planes. Hence, the 
impact parameter in the 1st lens plane is related to impact parameters in ``j``th lens plane as,
```math
\pmb{θ}_j = \pmb{θ}_1 - \sum_{i=1}^{j-1} \frac{D_{ij}}{D_j} \pmb{α}_i(\pmb{θ}_i).
```

With the above, the total potential and deflection for a lens system with $N_p$ planes is given as,
```math
\begin{align*}
ψ(\pmb{θ})       &= \sum_{i=1}^{N_p} \frac{D_{is}}{D_s} ψ_i(\pmb{θ}_i), \\
\pmb{α}(\pmb{θ}) &= \sum_{i=1}^{N_p} \frac{D_{is}}{D_s} \pmb{α}_i(\pmb{θ}_i).
\end{align*}
```

Similarly, the excess time delay accumulates between consecutive planes as,
```math
\begin{equation*}
Δt = \sum_{i=1}^{N_p} \frac{1 + z_i}{\rm c} \frac{D_i D_{i+1}}{D_{i,i+1}}
\left[ \frac{1}{2} \left| \pmb{θ}_i - \pmb{θ}_{i+1} \right|^2 
- \frac{D_{i,i+1}}{D_{i+1}} ψ_i(\pmb{θ}_i) \right],
\end{equation*}
```
where $z_i$ is the redshift of the ``i``th lens plane and the ($N_p$+1)th plane is the source 
plane, i.e., $\pmb{θ}_{N_p+1} = \pmb{β}$ and $D_{N_p+1} = D_s$.

In `LensFactory`, a multi-plane lens is initialized from a vector of lens components, where each 
component specifies its lens redshift (`z_d`) in addition to the usual lens model parameters. 
Components sharing the same redshift are automatically grouped into a single (composite) lens 
plane, and the lens planes are sorted in increasing redshift.

```@docs
Lenses.init_MultiPlaneLens
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