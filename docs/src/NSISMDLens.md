# NSISMD Lens
In the non-singular isothermal sphere mass distribution (NSISMD) lens, the singularity of the 
[SIS lens](SISLens.md) is removed by introducing a core radius directly in the mass distribution 
(in contrast to the [NSISP lens](NSISPLens.md), where the core is introduced in the lens 
potential). The projected surface mass density is given by
```math
\begin{equation*}
Σ(R) = \frac{σ_v^2}{2{\rm G}} \frac{1}{\sqrt{r_s^2 + R^2}},
\end{equation*}
```
where $σ_v$ is the velocity dispersion and $r_s$ is the core radius, which is finite at the 
center and reduces to the SIS profile for $r_s → 0$. Defining the angular core radius 
$θ_s = r_s / D_d$, the corresponding lens potential can be written as
```math
\begin{equation*}
ψ(\pmb{θ}) = 4π \left( \frac{σ_v}{\rm c} \right)^2 
\left[ \sqrt{θ_s^2 + |\pmb{θ} - \pmb{θ}_c|^2} 
- θ_s \ln \left( θ_s + \sqrt{θ_s^2 + |\pmb{θ} - \pmb{θ}_c|^2} \right) \right],
\end{equation*}
```
where $\pmb{θ}_c$ represents the lens center, leading to the deflection angle
```math
\begin{equation*}
\pmb{α}(\pmb{θ}) = 4π \left( \frac{σ_v}{\rm c} \right)^2 
\frac{\pmb{θ} - \pmb{θ}_c}{θ_s + \sqrt{θ_s^2 + |\pmb{θ} - \pmb{θ}_c|^2}}.
\end{equation*}
```
The Einstein angle for the NSISMD lens is given by
```math
\begin{equation*}
θ_E = \sqrt{θ_{E,0}^2 - 2 θ_s θ_{E,0}},
\qquad θ_{E,0} = 4π \left( \frac{σ_v}{\rm c} \right)^2 \frac{D_{ds}}{D_s},
\end{equation*}
```
where $θ_{E,0}$ is the Einstein angle of the corresponding SIS lens. Due to the finite core, a 
sufficiently well-aligned source can produce three images instead of two.

In `LensFactory`, to define an NSISMD lens, the user needs to specify four parameters: its 
position in the image plane ($x_c,~y_c$), its velocity dispersion ($v_d \equiv σ_v$, in km/s), and 
the core radius ($x_s \equiv θ_s$). By default, the lens is placed at the origin, i.e., 
($x_c,~y_c$) = (0, 0).

```@docs
Lenses.init_NSISMDLens
Lenses.NSISMDLens.potential!
Lenses.NSISMDLens.deflection!
Lenses.NSISMDLens.jacobian!
Lenses.NSISMDLens.einstein_angle
```