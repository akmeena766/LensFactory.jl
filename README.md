# LensFactory.jl
[![Build status (Github Actions)](https://github.com/akmeena766/LensFactory.jl/workflows/CI/badge.svg)](https://github.com/akmeena766/LensFactory.jl/actions)
[![codecov.io](http://codecov.io/github/akmeena766/LensFactory.jl/coverage.svg?branch=main)](http://codecov.io/github/akmeena766/LensFactory.jl?branch=main)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://akmeena766.github.io/LensFactory.jl/stable)
[![arXiv aaaa.bbbbb](https://img.shields.io/badge/arXiv-aaaa.bbbbb%20-yellowgreen.svg)](https://arxiv.org/abs/aaaa.bbbbb)

---
## Introduction
`LensFactory` is an efficient, open-source, general-purpose **strong gravitational lens modeling** 
package written in [Julia](https://julialang.org/). 
It supports single- and multi-plane lensing with a wide range of analytic lens profiles, and 
reconstructs lens mass distributions from multiple-image constraints — scaling from galaxy-scale 
lenses up to full galaxy-cluster reconstructions.


---
## Installation
To install, activate `pkg` mode by pressing `]` in Julia REPL and then type: 
```julia-repl
pkg> add LensFactory
```
To use the `LensFactory` package in your current Julia session, type: 
```julia-repl
julia> using LensFactory
```


---
## Quick start

```julia
using LensFactory

# Initialize the default cosmology
cosmo = Cosmology.init_cosmology()

# Lens and source redshifts
zl = 0.5
zs = 1.5

# Angular diameter distances and distance ratio
Dol = Cosmology.angular_diameter_distance(cosmo, 0., zl)
Dls = Cosmology.angular_diameter_distance(cosmo, zl, zs)
Dos = Cosmology.angular_diameter_distance(cosmo, 0., zs)
adis = Dls/Dos

# Initialize a point mass lens of mass 10^12 Solar mass
lens = Lenses.init_PointLens(D_d=Dol, mass=1E12)

# An arbitrary point in the image plane
x, y = 1.23, 0.57

# Get lens potential at the 
pot = Lenses.get_potential(lens, x, y)

# Get scaled deflection component at (1", 1")
dx, dy = Lenses.get_deflection(lens, x, y)

# Get magnification
mu = Lenses.get_magnification_image(lens, x, y, adis)
```

---
## Examples
To understand the use of various modules and functions in `LensFactory`, readers are encouraged to
go through the examples here: [LensFactory_Examples](https://github.com/akmeena766/LensFactory_Examples.git).
Please keep in mind that `LensFactory` is in **heavy** developement and examples are validated with the 
current **dev** version. If any of the examples are not working, please let me know.


---
## Citation
If you use `LensFactory.jl` in your research, please cite:
```bibtex
@article{LensFactory.jl,

}
```

---
## Feedback
Bug reports and feature requests are welcome via 
[GitHub issues](https://github.com/akmeena766/LensFactory.jl/issues).
For questions, suggestions, or comments, feel free to email
[akm@iisc.ac.in](mailto:akm@iisc.ac.in) or
[ashishmeena766@gmail.com](mailto:ashishmeena766@gmail.com).