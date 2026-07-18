# Sersic Lens
The Sersic profile [1963BAAA....6...41S](@cite) is widely used to describe the light (and stellar 
mass) distribution of galaxies. The corresponding convergence profile can be written as
```math
\begin{equation*}
κ(θ) = κ_s \exp\left[ -\left( \frac{|\pmb{θ} - \pmb{θ}_c|}{θ_s} \right)^{1/n} \right],
\end{equation*}
```
where $κ_s$ is the central convergence, $n$ is the Sersic index ($n = 4$ corresponds to the de 
Vaucouleurs profile and $n = 1$ to an exponential profile), and $\pmb{θ}_c$ represents the lens 
center. The scale radius $θ_s$ is related to the half-mass radius $θ_e$ as
```math
\begin{equation*}
θ_s = \frac{θ_e}{b_n^n},
\end{equation*}
```
where $b_n$ satisfies $Γ(2n) = 2γ(2n, b_n)$, with $γ(a, x)$ being the lower incomplete gamma 
function. For a total mass $M$, the central convergence is 
$κ_s \propto M / \left[ π θ_s^2 Γ(2n+1) \right]$, and the mass enclosed within $θ$ is analytic,
```math
\begin{equation*}
M(θ) = M \, P\left( 2n, \left( \frac{θ}{θ_s} \right)^{1/n} \right),
\end{equation*}
```
where $P(a, x) = γ(a, x)/Γ(a)$ is the regularized lower incomplete gamma function, leading to the 
deflection angle
```math
\begin{equation*}
\pmb{α}(\pmb{θ}) = \frac{4{\rm G} M}{{\rm c}^2} \frac{1}{D_d}
P\left( 2n, \left( \frac{|\pmb{θ} - \pmb{θ}_c|}{θ_s} \right)^{1/n} \right)
\frac{\pmb{θ} - \pmb{θ}_c}{|\pmb{θ} - \pmb{θ}_c|^2}.
\end{equation*}
```
The lens potential does not have an elementary closed form and is expressed in terms of the 
generalized hypergeometric function ${}_2F_2$,
```math
\begin{equation*}
ψ(\pmb{θ}) = \frac{2{\rm G} M}{{\rm c}^2} \frac{1}{D_d} \frac{1}{Γ(2n+1)}
\left( \frac{|\pmb{θ} - \pmb{θ}_c|}{θ_s} \right)^2
{}_2F_2\left( 2n, 2n; \, 2n+1, 2n+1; \, 
-\left( \frac{|\pmb{θ} - \pmb{θ}_c|}{θ_s} \right)^{1/n} \right).
\end{equation*}
```

In `LensFactory`, to define a Sersic lens, the user needs to specify six parameters: its position 
in the image plane ($x_c,~y_c$), its mass ($M$), the half-mass radius ($x_e \equiv θ_e$), the 
Sersic index ($n$, with a default value of 4), and the angular diameter distance to the lens 
($D_d$). By default, the lens is placed at the origin, i.e., ($x_c,~y_c$) = (0, 0). The helper 
function `scale_to_halflight` converts between the half-light radius and the scale radius.

```@docs
Lenses.init_SersicLens
Lenses.SersicLens.potential!
Lenses.SersicLens.deflection!
Lenses.SersicLens.jacobian!
Lenses.SersicLens.scale_to_halflight
```