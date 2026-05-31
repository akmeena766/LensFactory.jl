# LensModel

The `LensModel` module is the interface between the user and the lens modeling engine of 
`LensFactory.jl`. All user inputs are specified through a YAML file. The input YAML is divided into
multiple sections. The overall structure of the input YAML file is given below.

```yaml
observation:
   modeler: LensFactory
   lens: WhoKnows
   z_d: ...
   reference: [RA, DEC]
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
         x: [x1, x2]
         y: [y1, y2]
         sigma_x: [σx1, σx2]
         sigma_y: [σy1, σy2]
         sigma_theta: [σθ1, σθ2]
   .
   .
   .
   sourceS:
      z_s: zs
      total_knots: 1
      knot1:
         x: [x1, x2, x3]
         y: [y1, y2, y3]
         sigma_x: [σx1, σx2, σx3]
         sigma_y: [σy1, σy2, σy3]
         sigma_theta: [σθ1, σθ2, σθ3]

lens_model:
   multiplane: false
   total_lenses: N
   lens1:
      lens: SIELens
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

Looking at the above structure, we note that input YAML file is divided into five sections, i.e.,
`observation`, `cosmology`, `source`, `lens_model`, and `sampling`. Below, we describe each section
in more detail.

### observation
The `observation` section contains the following keywords.

- `modeler`: **(optional)** Modeler name (default: `LensFactory`).
- `lens`: **(optional)** Lens name (default: `WhoKnows`).
- `z_d`: **(mandatory)** Lens redshift.
- `reference`: **(mandatory)** Reference `[RA, DEC]` position of the lens in degrees. If `[0.0, 0.0]` is provided 
   then it would be assumed that the we are working in arcseconds. 
- `pixel_scale`: **(mandatory)** Pixel scale in arcseconds.
- `FOV`: **(mandatory)** Field of view in arcseconds. The user can either specify a single value 
   (e.g., 20) or a vector of two values [x, y] (e.g., [10, 10]).  If a single value is provided 
   then square grid of `[-value, value]` will be used otherwise a rectangular grid with size [-x, x]
   and [-y, y] will be used.


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