# LensFactory.jl

[![Build status (Github Actions)](https://github.com/akmeena766/LensFactory.jl/workflows/CI/badge.svg)](https://github.com/akmeena766/LensFactory.jl/actions)
[![codecov.io](http://codecov.io/github/akmeena766/LensFactory.jl/coverage.svg?branch=main)](http://codecov.io/github/akmeena766/LensFactory.jl?branch=main)
[![arXiv aaaa.bbbbb](https://img.shields.io/badge/arXiv-aaaa.bbbbb%20-yellowgreen.svg)](https://arxiv.org/abs/aaaa.bbbbb)


---
## Introduction
`LensFactory` is a general-purpose gravitational lensing package fully developed in Julia. 

At present `LensFactory` supports the following functionalities: 
- [Single-plane lensing](Lenses.md)
- [Multi-plane lensing](MultiPlane.md)
- [Strong lens modeling](LensModel.md)
- [Singularity map](SingularityMap.md)

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
## Examples
To understand the use of various modules and functions in `LensFactory`, readers are encouraged to
go through the examples here: [LensFactory-Examples](https://github.com/akmeena766/LensFactory_Examples.git).

!!! note
    Please keep in mind that `LensFactory` is in **heavy** developement and examples are validated 
    with the current **dev** version. If any of the examples are not working, please let me know.

---
## Feedback
I would be very happy to receive any suggestions, comments, or questions regarding `LensFactory`.
Please feel free to drop an email to: [akm@iisc.ac.in](mailto:akm@iisc.ac.in) or
[ashishmeena766@gmail.com](mailto:ashishmeena766@gmail.com).

---
## Citation
