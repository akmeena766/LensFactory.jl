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
      param1:
         value: p1v
         range: [p1l, p1u]
      param2:
         value: p2v
         range: [p2l, p2u]
      param3:
         value: p3v
         range: [p3l, p3u]
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
in more detail. Overall, the parameters can be divided into three categories: optional, mandetory, 
and free. To make a parameter free, the user needs to specify initial `value` and `range` keywords
for the same.

### observation
The `observation` section contains the following keywords.
- `modeler`: **(optional)** Modeler name (default: `LensFactory`).
- `lens`: **(optional)** Lens name (default: `WhoKnows`).
- `z_d`: **(mandatory)** Lens redshift.
- `reference`: **(mandatory)** Reference `[RA, DEC]` position of the lens in degrees. If 
   `[0.0, 0.0]` is provided then it would be assumed that the we are working with coordinates in 
   arcseconds. 
- `pixel_scale`: **(mandatory)** Pixel scale in arcseconds.
- `FOV`: **(mandatory)** Field of view in arcseconds. The user can either specify a single value 
   (e.g., 20) or a vector of two values [x, y] (e.g., [10, 10]).  If a single value is provided 
   then square grid of `[-value, value]` will be used otherwise a rectangular grid with size [-x, x]
   and [-y, y] will be used.

### cosmology
The `cosmology` section specifies the underlying cosmology and contains the following keywords. Each
one of these parameters can be left free.
- `H0`: **(free)** Hubble parameter (default: 70.0).
- `Omega_m0`: **(free)** Matter density parameter (default: 0.30).
- `Omega_r0`: **(free)** Radiation density parameter (default: 0.00).
- `Omega_w0`: **(free)** Dark energy density parameter (default: 0.70).

!!! note
    This section can be considered as an optional section. If no comoslogy is provided then the code 
    will assume the default [cosmology](#LensFactory.Cosmology.init_cosmology). However, it is highly
    recommended to explicity provide the cosmolgy.


### source
The `source` section contains the details about the observed strongly lensed sources. The 
corresponding details can be specified in two ways, either directly in the YAML file or through a 
seperate .txt file. Each source will have an associated redshift and a number of knots.


- `total_sources`: **(mandatory)** Total number of sources.
- `from_file`: **(optional)** Boolean value indicating whether the source model is read from a file 
   (default: false). If `true`, the source parameters will be read from the file specified in the 
   `source_file` keyword.
- `source_file`: **(optional)** Path to the .txt file containing the source parameters. This keyword 
   is only required when `from_file` is set to `true`. 
- `source1`: Details of the first source.
  - `z_s`:  **(mandatory)** Source redshift.
  - `total_knots`: **(mandatory)** Total number of strongly lensed knots within the source.
  - `knot1`: Details of the first knot.
    - `x`: **(mandatory)** x-coordinates of the knots .
    - `y`: **(mandatory)** y-coordinates of the knots .
    - `sigma_x`: **(mandatory)** Standard deviation in x-coordinates .
    - `sigma_y`: **(mandatory)** Standard deviation in y-coordinates .
    - `sigma_theta`: **(mandatory)** Standard deviation in angles .
  - `. . .`
  - `knotK`: Model for the last knot.
- `. . .`
- `sourceS`: Model for the last source.

The format of the .txt file is as follows.
```julia
| sourceID | knotID | x | y | z_s | σx | σy | σθ | parity | magnitude | time delay |
```

### lens_model
The `lens_model` section contains the following keywords.

- `multiplane`: Boolean value indicating whether the lens model is multiplane (default: false).
- `total_lenses`: Total number of lenses (mandatory).
- `lens1`: Model for the first lens (mandatory).
  - `lens`: Lens name (mandatory).
  - `x_c`: x-coordinate of the lens center (mandatory).
  - `y_c`: y-coordinate of the lens center (mandatory).
  - `param1`: Parameter 1 (mandatory).
  - `. . .`
  - `paramN`: Parameter N (mandatory).
- `. . .`
- `lensN`: Model for the last lens.

### sampling
The `sampling` section contains the following keywords.

- `scheme`: Sampling scheme (default: `sourceplane`). The possible values are `sourceplane`, `imageplane`, `image_plus_sourceplane`, and `image_plus_image_plus_sourceplane`.
- `verbose`: Boolean value indicating whether the sampling should be verbose (default: `true`).
- `optimizer`: Optimizer to be used for the sampling. The default optimizer is `L-BFGS-B`.
- `mcmc`: Markov chain Monte Carlo parameters. The following keywords are available.
  - `iterations`: Number of iterations (default: 1000).
  - `burn_in`: Number of burn-in iterations (default: 500).
  - `thin`: Thinning factor (default: 1).
  - `temperature`: Temperature parameter (default: 1.0).
  - `step_size`: Step size (default: 1.0).


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