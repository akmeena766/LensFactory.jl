# SIS Lens
The singular isothermal sphere (SIS) is one of the most widely used lens models for galaxies, as it 
naturally reproduces flat rotation curves. Its three-dimensional density profile is given by
```math
\begin{equation*}
ρ(r) = \frac{σ_v^2}{2π{\rm G}} \frac{1}{r^2},
\end{equation*}
```
where $σ_v$ is the (one-dimensional) velocity dispersion. The corresponding projected surface mass 
density is
```math
\begin{equation*}
Σ(R) = \frac{σ_v^2}{2{\rm G}} \frac{1}{R},
\end{equation*}
```
and the lens potential can be written as
```math
\begin{equation*}
ψ(\pmb{θ}) = 4π \left( \frac{σ_v}{\rm c} \right)^2 |\pmb{θ} - \pmb{θ}_c|,
\end{equation*}
```
where $\pmb{θ}_c$ represents the lens center, leading to a deflection angle of constant magnitude,
```math
\begin{equation*}
\pmb{α}(\pmb{θ}) = 4π \left( \frac{σ_v}{\rm c} \right)^2 
\frac{\pmb{θ} - \pmb{θ}_c}{|\pmb{θ} - \pmb{θ}_c|}.
\end{equation*}
```
The Einstein angle for the SIS lens is given by
```math
\begin{equation*}
θ_E = 4π \left( \frac{σ_v}{\rm c} \right)^2 \frac{D_{ds}}{D_s},
\end{equation*}
```
and a source at $|\pmb{β} - \pmb{θ}_c| < θ_E$ produces two images, whereas a source outside the 
Einstein radius is only singly imaged.

In `LensFactory`, to define an SIS lens, the user needs to specify three parameters: its position 
in the image plane ($x_c,~y_c$) and its velocity dispersion ($v_d \equiv σ_v$, in km/s). By 
default, the lens is placed at the origin, i.e., ($x_c,~y_c$) = (0, 0).

```@docs
Lenses.init_SISLens
Lenses.SISLens.potential!
Lenses.SISLens.deflection!
Lenses.SISLens.jacobian!
Lenses.SISLens.einstein_angle
```