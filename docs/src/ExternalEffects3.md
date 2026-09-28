# External Effects: Third Order
Going one order beyond the constant convergence and shear described in 
[External Effects](ExternalEffects.md), the environment also introduces higher-order perturbations 
in the lens potential. The generic third-order perturbation can be written as,
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
(see [External Effects](ExternalEffects.md)). Again, we can elimiate the convergence terms by rescaling the source positions, leaving us with,
```math
\begin{align*}
ψ =& \frac{\gamma_1}{2}(θ_x^2 - θ_y^2) + \gamma_2 θ_x θ_y \\
   & + \frac{1}{6}(θ_x^3 ψ_{xxx} + θ_y^3 ψ_{yyy} + 3 θ_x^2 θ_y ψ_{xxy} + 3 θ_x θ_y^2 ψ_{xyy}).
\end{align*}
```
Going from cartesian to polar coordinates, i.e., $(x, y) = (r \cosθ, \sinθ)$,
```math
\begin{align*}
ψ =& \frac{\gamma_1}{2}r^2 \cos2θ + \frac{\gamma_2}{2}r^2 \sin2θ \\
   & + \frac{1}{6}(θ_x^3 ψ_{xxx} + θ_y^3 ψ_{yyy} + 3 θ_x^2 θ_y ψ_{xxy} + 3 θ_x θ_y^2 ψ_{xyy}).
\end{align*}
```
Substituting $\gamma = \sqrt{\gamma_1^2 + \gamma2_2^2}$ and $\tan2θ_\gamma = \gamma_2/\gamma_1$ and
$A\cosθ + A\sinθ = \sqrt{A^2+B^2}\cos(θ-θ_\gamma)$, we get
```math
\begin{equation*}
ψ = \frac{\gamma}{2}r^2 \cos2(θ-θ_\gamma) + \frac{\sigma}{4}r^3 \cos(θ-θ_\sigma) + \frac{\delta}{6}r^3 \cos3(θ-θ_\delta).
\end{equation*}
```



Assuming the perturber to be an SIS, the third-order perturbation takes a "restricted" 
one-parameter form,
```math
\begin{equation*}
ψ(\pmb{θ}) = δ \, θ^3 \cos(φ - ϕ) \sin^2(φ - ϕ),
\end{equation*}
```
where ($θ,~φ$) are the polar coordinates in the image plane, $δ$ is the amplitude of the 
perturbation, and $ϕ$ is its direction (i.e., the direction towards the SIS perturber, which 
coincides with the external shear angle it produces).

In `LensFactory`, to define third-order external effects, the user needs to specify two 
parameters: the amplitude of the perturbation ($δ$) and its direction ($ϕ$, in degrees). The 
perturbation is always centered at the origin of the image plane.

```@docs
Lenses.init_ExternalEffects3
Lenses.ExternalEffects3.potential!
Lenses.ExternalEffects3.deflection!
Lenses.ExternalEffects3.jacobian!
```