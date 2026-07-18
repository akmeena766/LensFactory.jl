# Gaussian Lens
The Gaussian lens is described by a projected surface mass density of Gaussian form,
```math
\begin{equation*}
Σ(R) = \frac{M}{2π σ^2} \exp\left( -\frac{R^2}{2σ^2} \right),
\end{equation*}
```
where $M$ is the total mass of the lens and $σ$ is the standard deviation of the Gaussian. Thanks 
to its finite total mass and smooth analytic profile, it is often used as a building block to 
represent more complex mass distributions (see the 
[Multi-Gaussian lens](MultiGaussianLens.md)). Defining the angular scale radius $θ_s = σ / D_d$, 
the convergence profile can be written as
```math
\begin{equation*}
κ(θ) = κ_s \exp\left( -\frac{|\pmb{θ} - \pmb{θ}_c|^2}{2 θ_s^2} \right),
\end{equation*}
```
where $κ_s$ is the central convergence and $\pmb{θ}_c$ represents the lens center. Since the mass 
enclosed within $θ$ is analytic, 
$M(θ) = M \left[ 1 - \exp\left(-θ^2/2θ_s^2\right) \right]$, the deflection angle takes a simple 
closed form,
```math
\begin{equation*}
\pmb{α}(\pmb{θ}) = \frac{4{\rm G} M}{{\rm c}^2} \frac{1}{D_d}
\left[ 1 - \exp\left( -\frac{|\pmb{θ} - \pmb{θ}_c|^2}{2 θ_s^2} \right) \right]
\frac{\pmb{θ} - \pmb{θ}_c}{|\pmb{θ} - \pmb{θ}_c|^2},
\end{equation*}
```
whereas the lens potential involves the exponential integral function $\mathrm{Ei}(x)$,
```math
\begin{equation*}
ψ(\pmb{θ}) = 2 κ_s θ_s^2 \left[ \ln\left( \frac{|\pmb{θ} - \pmb{θ}_c|}{θ_s} \right)
- \frac{1}{2} \mathrm{Ei}\left( -\frac{|\pmb{θ} - \pmb{θ}_c|^2}{2 θ_s^2} \right) \right].
\end{equation*}
```
For $θ \gg θ_s$, the Gaussian lens behaves as a point mass lens of mass $M$.

In `LensFactory`, to define a Gaussian lens, the user needs to specify five parameters: its 
position in the image plane ($x_c,~y_c$), its mass ($M$), the scale radius ($x_s \equiv θ_s$), and 
the angular diameter distance to the lens ($D_d$). By default, the lens is placed at the origin, 
i.e., ($x_c,~y_c$) = (0, 0).

```@docs
Lenses.init_GaussianLens
Lenses.GaussianLens.potential!
Lenses.GaussianLens.deflection!
Lenses.GaussianLens.jacobian!
```