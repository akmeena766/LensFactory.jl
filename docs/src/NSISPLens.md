# NSISP Lens
The singularity of the [SIS lens](SISLens.md) at the center can be removed by introducing a finite 
core. There are two ways to do this: softening the lens potential itself or softening the mass 
distribution (see [NSISMD lens](NSISMDLens.md) for the latter). In the non-singular isothermal 
sphere potential (NSISP) lens, the core radius ($θ_s$) is introduced directly in the SIS lens 
potential,
```math
\begin{equation*}
ψ(\pmb{θ}) = 4π \left( \frac{σ_v}{\rm c} \right)^2 \sqrt{θ_s^2 + |\pmb{θ} - \pmb{θ}_c|^2},
\end{equation*}
```
where $σ_v$ is the velocity dispersion and $\pmb{θ}_c$ represents the lens center, leading to the 
deflection angle
```math
\begin{equation*}
\pmb{α}(\pmb{θ}) = 4π \left( \frac{σ_v}{\rm c} \right)^2 
\frac{\pmb{θ} - \pmb{θ}_c}{\sqrt{θ_s^2 + |\pmb{θ} - \pmb{θ}_c|^2}}.
\end{equation*}
```
The corresponding convergence profile is
```math
\begin{equation*}
κ(θ) = \frac{θ_{E,0}}{2} \frac{2θ_s^2 + θ^2}{\left( θ_s^2 + θ^2 \right)^{3/2}},
\qquad θ_{E,0} = 4π \left( \frac{σ_v}{\rm c} \right)^2 \frac{D_{ds}}{D_s},
\end{equation*}
```
which is finite at the center, and the Einstein angle is given by
```math
\begin{equation*}
θ_E = \sqrt{θ_{E,0}^2 - θ_s^2}.
\end{equation*}
```
In the limit $θ_s → 0$, the NSISP lens reduces to the SIS lens. Due to the finite core, a 
sufficiently well-aligned source can produce three images instead of two.

In `LensFactory`, to define an NSISP lens, the user needs to specify four parameters: its position 
in the image plane ($x_c,~y_c$), its velocity dispersion ($v_d \equiv σ_v$, in km/s), and the core 
radius ($x_s \equiv θ_s$). By default, the lens is placed at the origin, i.e., 
($x_c,~y_c$) = (0, 0).

```@docs
Lenses.init_NSISPLens
Lenses.NSISPLens.potential!
Lenses.NSISPLens.deflection!
Lenses.NSISPLens.jacobian!
Lenses.NSISPLens.einstein_angle
```