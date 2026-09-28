# External Effects: Third Order
Going one order beyond the constant convergence and shear described in 
[External Effects](ExternalEffects.md), the environment also introduces higher-order perturbations 
in the lens potential. The generic third-order perturbation can be written as 
[1991ApJ...373..354K, 1999AJ....118...14B](@cite),
```math
\begin{align*}
ψ =& ψ_0 + θ_x ψ_x + θ_y ψ_y + \frac{1}{2}(θ_x^2 ψ_{xx} + θ_y^2 ψ_{yy} + 2 θ_x θ_y ψ_{xy}) \\
   & + \frac{1}{6}(θ_x^3 ψ_{xxx} + θ_y^3 ψ_{yyy} + 3 θ_x^2 θ_y ψ_{xxy} + 3 θ_x θ_y^2 ψ_{xyy}).
\end{align*}
```
We can drop the first three terms: $ψ_0$ is a constant and the two subsequent terms can be absorbed
in source rescaling. With that, the above equation can be written as,
```math
\begin{align*}
ψ =& \frac{\kappa}{2}(θ_x^2 + θ_y^2) +  \frac{\gamma_1}{2}(θ_x^2 - θ_y^2) + \gamma_2 θ_x θ_y \\
   & + \frac{1}{6}(θ_x^3 ψ_{xxx} + θ_y^3 ψ_{yyy} + 3 θ_x^2 θ_y ψ_{xxy} + 3 θ_x θ_y^2 ψ_{xyy}),
\end{align*}
```
where we have introduced the well known convergence and shear terms 
(see [External Effects](ExternalEffects.md)). Again, we can elimiate the convergence terms by 
rescaling the source positions, leaving us with,
```math
\begin{align*}
ψ =& \frac{\gamma_1}{2}(θ_x^2 - θ_y^2) + \gamma_2 θ_x θ_y \\
   & + \frac{1}{6}(θ_x^3 ψ_{xxx} + θ_y^3 ψ_{yyy} + 3 θ_x^2 θ_y ψ_{xxy} + 3 θ_x θ_y^2 ψ_{xyy}).
\end{align*}
```
Going from cartesian to polar coordinates, i.e., $(θ_x, θ_y) = (r \cosθ, \sinθ)$,
```math
\begin{align*}
ψ =& \frac{\gamma_1}{2}r^2 \cos2θ + \frac{\gamma_2}{2}r^2 \sin2θ \\
   & + \frac{r^3}{8}  \left[(ψ_{xxx} + ψ_{xyy}) \cosθ  + (ψ_{xxy} + ψ_{yyy})\sinθ \right] \\
   & + \frac{r^3}{24} \left[(ψ_{xxx} -3ψ_{xyy}) \cos3θ + (3ψ_{xxy}- ψ_{yyy})\sin3θ\right].
\end{align*}
```
```math
\begin{align*}
ψ =& \frac{r^2}{2} \left[\gamma_1 \cos2θ + \gamma_2 \sin2θ\right] \\
   & + \frac{r^3}{8}  \left[\partial_x \kappa \cosθ  + \partial_y\kappa \sinθ \right] \\
   & + \frac{r^3}{24} \left[(ψ_{xxx} -3ψ_{xyy}) \cos3θ + (3ψ_{xxy}- ψ_{yyy})\sin3θ\right].
\end{align*}
```
Using $A\cosθ + B\sinθ = \sqrt{A^2+B^2}\cos(θ-θ_0)$, we get,
```math
\begin{equation*}
\boxed{ψ = \frac{\gamma}{2}r^2 \cos2(θ-θ_\gamma) + 
           \frac{\sigma}{4}r^3 \cos(θ-θ_\sigma) + 
           \frac{\delta}{6}r^3 \cos3(θ-θ_\delta)}.
\end{equation*}
```

Assuming the perturber to be an SIS, the third-order perturbation takes a "restricted" one-parameter
form,
```math
\begin{equation*}
ψ(\pmb{θ}) = \frac{\Delta}{4} r^3 \left[\cos(θ - θ_\Delta) + \cos3(θ - θ_\Delta) \right].
\end{equation*}
```

In `LensFactory`, to define third-order external effects, the user needs to specify two 
parameters: the amplitude of the perturbation ($\Delta$) and its direction ($θ_\Delta$, in degrees). 
The perturbation is always centered at the origin of the image plane.

```@docs
Lenses.init_ExternalEffects3
Lenses.ExternalEffects3.potential!
Lenses.ExternalEffects3.deflection!
Lenses.ExternalEffects3.jacobian!
```