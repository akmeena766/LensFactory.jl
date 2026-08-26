# LensModel

The `LensModel` module is the interface between the user and the lens modeling engine of 
`LensFactory.jl`. All user inputs are specified through a YAML file. The input YAML is divided into
multiple sections. The overall structure of the input YAML file is given below.


```@docs
LensModel.read_input
LensModel.fit_model
LensModel.free_parameter_names
LensModel.get_best_fit_parameters
LensModel.get_cosmology
LensModel.get_adis
LensModel.get_best_model
LensModel.get_potential
LensModel.get_deflection
LensModel.get_jacobian
LensModel.get_magnification_image
LensModel.get_source_position
LensModel.predict_image
LensModel.save_best_fits
LensModel.get_best_fit_rms
LensModel.error_models
LensModel.get_AIC
LensModel.get_BIC
```