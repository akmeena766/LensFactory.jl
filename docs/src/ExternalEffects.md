# External Effects
A lens galaxy is rarely isolated: nearby galaxies and the large-scale environment perturb the 
lensing signal. To the lowest order, these perturbations can be described by a constant external 
convergence ($κ$) and a constant external shear ($γ$) with position angle ($ϕ$). The corresponding 
lens potential is given by
```math
\begin{equation*}
ψ(\pmb{θ}) = \frac{κ + γ_1}{2} θ_x^2 + \frac{κ - γ_1}{2} θ_y^2 + γ_2 \, θ_x θ_y,
\end{equation*}
```
where the two shear components are defined as
```math
\begin{equation*}
γ_1 = γ \cos(2ϕ), \qquad γ_2 = γ \sin(2ϕ).
\end{equation*}
```
The external convergence produces an isotropic (de-)magnification, whereas the external shear 
stretches the images along the direction $ϕ$. Since the external convergence is degenerate with 
the lens mass (mass-sheet degeneracy), it is often fixed to zero in lens modeling.

In `LensFactory`, to define external effects, the user needs to specify three parameters: the 
external convergence ($κ$), the shear amplitude ($γ$), and the shear angle ($ϕ$, in degrees). The 
external effects are always centered at the origin of the image plane.

```@docs
Lenses.init_ExternalEffects
Lenses.ExternalEffects.potential!
Lenses.ExternalEffects.deflection!
Lenses.ExternalEffects.jacobian!
```