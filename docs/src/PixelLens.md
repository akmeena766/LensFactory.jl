# Pixel Lens
The pixel lens is the simplest extended mass element one can construct: a square tile of constant
surface mass density. Its surface mass density is

```math
\begin{equation*}
Σ(\pmb{θ}) =
\begin{cases}
κ \, Σ_{\rm cr}, & |θ_x - θ_{xc}| \le a/2 \ \ \text{and} \ \ |θ_y - θ_{yc}| \le a/2, \\[4pt]
0,             & \text{otherwise},
\end{cases}
\end{equation*}
```
where $\pmb{θ}_c = (θ_{xc},~θ_{yc})$ is the centre of the tile, $a$ is its side length, and $κ$ is
the dimensionless convergence expressed in units of the critical surface density $Σ_{\rm cr}$ for
sources at infinity. The tile is fully specified by these three numbers, and the mass it contains is
simply

```math
\begin{equation*}
M = κ \, Σ_{\rm cr} \, \left(a \, D_d\right)^2 .
\end{equation*}
```
Unlike the other lenses in `LensFactory`, the pixel lens has a sharp boundary and is uniform inside
it. It is therefore not intended as a description of any physical object on its own; its usefulness
comes from the fact that it is an exactly integrable building block out of which arbitrary mass
distributions can be assembled.

With the normalization used throughout `LensFactory`, the lens potential of an arbitrary convergence
distribution is
```math
\begin{equation*}
ψ(\pmb{θ}) = \frac{1}{π} \int κ(\pmb{θ}') \, \ln \left|\pmb{θ} - \pmb{θ}'\right| \, {\rm d}^2 θ',
\end{equation*}
```

so that $\nabla^2 ψ = 2 κ$ and $\pmb{α} = \nabla ψ$. For a uniform square the integrand is separable
in Cartesian coordinates over a rectangular domain, and the integral can be evaluated in closed
form. It is convenient to introduce the offsets of the field point from the four pixel edges,
```math
\begin{equation*}
u_{1,2} = θ_x - θ_{xc} \pm \frac{a}{2}, \qquad
v_{1,2} = θ_y - θ_{yc} \pm \frac{a}{2},
\end{equation*}
```
together with the alternating corner sum
```math
\begin{equation*}
\mathcal{D}\left[f\right] \equiv \sum_{i=1}^{2} \sum_{j=1}^{2} (-1)^{i+j} f(u_i, v_j)
= f(u_1, v_1) - f(u_1, v_2) - f(u_2, v_1) + f(u_2, v_2).
\end{equation*}
```

Every quantity below is obtained by evaluating one elementary function at the four corners of the
pixel and combining them with $\mathcal{D}$. The lens potential is
```math
\begin{equation*}
ψ(\pmb{θ}) = \frac{κ}{2π} \,
\mathcal{D}\!\left[\, u v \ln\left(u^2 + v^2\right) - 3 u v
+ u^2 \arctan\!\left(\frac{v}{u}\right)
+ v^2 \arctan\!\left(\frac{u}{v}\right) \right].
\end{equation*}
```

The deflection components are,
```math
\begin{align*}
α_x(\pmb{θ}) &= \frac{κ}{2π} \,
\mathcal{D}\!\left[\, v \ln\left(u^2 + v^2\right)
+ 2 u \arctan\!\left(\frac{v}{u}\right) \right], \\

α_y(\pmb{θ}) &= \frac{κ}{2π} \,
\mathcal{D}\!\left[\, u \ln\left(u^2 + v^2\right)
+ 2 v \arctan\!\left(\frac{u}{v}\right) \right].
\end{align*}
```

Finally, the deformation tensor components are
```math
\begin{align*}
ψ_{xx}(\pmb{θ}) &= \frac{κ}{π} \, \mathcal{D}\!\left[\arctan\!\left(\frac{v}{u}\right)\right], \\
ψ_{yy}(\pmb{θ}) &= \frac{κ}{π} \, \mathcal{D}\!\left[\arctan\!\left(\frac{u}{v}\right)\right], \\
ψ_{xy}(\pmb{θ}) &= \frac{κ}{2π} \, \mathcal{D}\!\left[\ln\left(u^2 + v^2\right)\right].
\end{align*}
```

Because all of the expressions above are linear in $κ$, a grid of such tiles with independent
convergences is the standard starting point for free-form (non-parametric) mass reconstruction,
where the observed image positions become linear constraints on the pixel values
[1998MNRAS.294..734A](@cite). See [Multi-Pixel Lens](@ref Multi-component-Pixel-Lens), which 
evaluates the same expressions for a vector of pixel centres and convergences.

```@docs
Lenses.init_PixelLens
Lenses.PixelLens.potential!
Lenses.PixelLens.deflection!
Lenses.PixelLens.jacobian!
```
