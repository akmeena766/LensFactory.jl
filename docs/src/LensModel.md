# LensModel

The `LensModel` module is the interface between the user and the lens modeling engine of 
`LensFactory.jl`. All user inputs are specified through a YAML file. The input YAML is divided into
multiple sections. The overall structure of the input YAML file is described below.

```
observation:
   modeler: LensFactory
   lens: WhoKnows
   z_d: ...
   reference: [0.0, 0.0]
   pixel_scale: ...
   FOV: [x, y]

cosmology:
   H0: 70.0
   Omega_m0: 0.30
   Omega_r0: 0.00
   Omega_w0: 0.70

source:
   total_sources: S
   from_file: false
   source1:
      z_s:
         value: zv
         range: [zl, zu]
      total_knots: K
      knot1:
         x: [x1, x2, x3, x4]
         y: [y1, y2, y3, y4]
         sigma_x: [σx1, σx2, σx3, σx4]
         sigma_y: [σy1, σy2, σy3, σy4]
         sigma_theta: [σθ1, σθ2, σθ3, σθ4]
      .
      .
      .
      knotK:
         x: ...
         y: ...
         sigma_x: ...
         sigma_y: ...
         sigma_theta: ...
   .
   .
   .
   sourceS:
      z_s: ...
      total_knots: 1
      knot1: ...

lens_model:
   multiplane: false
   total_lenses: N
   lens1:
      lens: ...
      x_c:
         value: xv
         range: [xl, xu]
      y_c:
         value: yv
         range: [yl, yu]
      param1:
         value: pv
         range: [pl, pu]
      .
      .
      .
   .
   .
   .
   lensN:
      param1: ...
      param2: ...
      .
      .
      .

sampling:
   scheme: sourceplane
   verbose: true
   optimizer:
   mcmc:
```

```@docs
LensModel.read_input
LensModel.fit_model
LensModel.free_parameter_names
LensModel.get_best_fit_parameters
LensModel.get_best_model
LensModel.get_potential
LensModel.get_deflection
LensModel.get_jacobian
LensModel.predict_image
LensModel.save_best_fits
LensModel.check_parity
LensModel.get_best_fit_rms
```