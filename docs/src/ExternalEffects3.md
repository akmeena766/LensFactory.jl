# External Effects: Third Order
Going one order beyond the constant convergence and shear described in 
[External Effects](ExternalEffects.md), the environment also introduces third-order (flexion-like) 
perturbations in the lens potential. Assuming the perturber to be an SIS, the third-order 
perturbation takes a "restricted" one-parameter form,
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