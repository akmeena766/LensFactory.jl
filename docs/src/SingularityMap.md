# SingularityMap
According to catastrophe theory, the stable singularities of the gravitational lens map can be 
classified into folds (``A_2``), which form the critical curves, cusps (``A_3``), and the 
higher-order swallowtail (``A_4``) and umbilic (``D_4``) singularities, which only occur at 
isolated points in the image plane. For a fixed lens, as we vary the source distance, the critical 
curves sweep through the image plane, and the higher-order singularities appear at specific source 
distances. The `SingularityMap` module maps out these singularities in the image plane without 
requiring a specific source distance.

Writing the (distance-scaled) Jacobian of the lens map as 
``\mathbb{A} = \mathbb{I} - a \, ψ_{,ij}``, where $a$ is the source distance ratio 
($a = D_{ds}/D_s$ scaled appropriately) and $ψ_{,ij}$ is the Hessian of the lens potential, we 
denote the eigenvalues and (unit) eigenvectors of $ψ_{,ij}$ as ($λ_1,~λ_2$) and 
($\pmb{q}_1,~\pmb{q}_2$), respectively. A point in the image plane lies on a critical curve for a 
source distance such that $a λ_i = 1$. The singularity map is then built from the following 
conditions:

- **``A_3``-lines**: the loci of points that become cusps for some source distance, satisfying 
  $\pmb{q}_i \cdot \pmb{∇} λ_i = 0$, i.e., the gradient of an eigenvalue along the corresponding 
  eigenvector vanishes.
- **``A_4``-points**: the points on ``A_3``-lines where a swallowtail singularity can occur, 
  corresponding to extrema of the eigenvalue $λ_i$ along the ``A_3``-line.
- **``D_4``-points**: the points where an umbilic singularity can occur, i.e., where both 
  eigenvalues coincide ($λ_1 = λ_2$), equivalent to $ψ_{,11} = ψ_{,22}$ and $ψ_{,12} = 0$.

Since a singularity is only physically realizable if a critical curve can actually cross it, only 
the points satisfying $a \, λ_i ≥ 1$ (with $a$ being the maximum allowed distance ratio, set by 
the keyword `adis`) are retained in the final map.

The singularity map can be constructed either directly from a lens instance (`from_lens`) or from 
precomputed Jacobian component maps (`from_jacobian`). By default, only the ``A_3``-lines are 
computed; the ``A_4``- and ``D_4``-points can be requested via the corresponding keyword 
arguments.

```@docs
SingularityMap.from_lens
SingularityMap.from_jacobian
```