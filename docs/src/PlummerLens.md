# Plummer Lens
The Plummer lens is based on the Plummer model [1911MNRAS..71..460P](@cite), originally introduced 
to describe the density distribution in globular clusters. Its three-dimensional density profile is 
given by
```math
\begin{equation*}
ρ(r) = \frac{3M}{4π a^3} \left( 1 + \frac{r^2}{a^2} \right)^{-5/2},
\end{equation*}
```
where $M$ is the total mass of the lens and $a$ is the Plummer (core) radius. The corresponding 
projected surface mass density is
```math
\begin{equation*}
Σ(R) = \frac{M}{π a^2} \left( 1 + \frac{R^2}{a^2} \right)^{-2}.
\end{equation*}
```
Defining the angular core radius $θ_s = a / D_d$, the lens potential can be written as
```math
\begin{equation*}
ψ(\pmb{θ}) = \frac{4{\rm G} M} {\rm c^2} \frac{1}{D_d} \ln\left[ \sqrt{θ_s^2 + |\pmb{θ} - \pmb{θ}_c|^2} \right],
\end{equation*}
```
where $\pmb{θ}_c$ represents the lens center, leading to the deflection angle
```math
\begin{equation*}
\pmb{α}(\pmb{θ}) = \frac{4{\rm G} M} {\rm c^2} \frac{1}{D_d} 
\frac{\pmb{θ} - \pmb{θ}_c}{θ_s^2 + |\pmb{θ} - \pmb{θ}_c|^2}.
\end{equation*}
```
In the limit $θ_s → 0$, the above equations reduce to the corresponding point mass lens 
equations, i.e., the Plummer lens acts as a softened point mass lens.

In `LensFactory`, to define a Plummer lens, the user needs to specify five parameters: its position 
in the image plane ($x_c,~y_c$), its mass ($M$), the core radius ($x_s \equiv θ_s$), and the 
angular diameter distance to the lens ($D_d$). By default, the lens is placed at the origin, i.e., 
($x_c,~y_c$) = (0, 0).

```@docs
Lenses.init_PlummerLens
Lenses.PlummerLens.potential!
Lenses.PlummerLens.deflection!
Lenses.PlummerLens.jacobian!
Lenses.PlummerLens.einstein_angle
```