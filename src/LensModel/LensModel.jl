module LensModel


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..Constants
using ..Cosmology
using ..Lenses
using ..LFUtils


# --------------------------------------------------------------------------------------------------
# Helper modules to use
# --------------------------------------------------------------------------------------------------
include("./NelderMead.jl")
using .NelderMead

include("./MH.jl")
using .MH

include("./AIES.jl")
using .AIES

include("./LensModelIO.jl")
using .LensModelIO

include("./Likelihood.jl")
using .Likelihood

include("./LensModelUtils.jl")
using .LensModelUtils

include("./LensModelFit.jl")
using .LensModelFit


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export read_input
export fit_model
export free_parameter_names
export calculate_gr
export print_gr_report
export time_series_diagnostics
export acceptance_diagnostics
export get_best_fit
export get_best_fit_with_errors
export check_parity
export get_best_fit_rms
export save_best_fit



# --------------------------------------------------------------------------------------------------
# Plotting functions (see ../../ext folder for functions)
# --------------------------------------------------------------------------------------------------
export plot_corner
export plot_trace
export plot_best_model

function plot_corner end
function plot_trace end
function plot_best_model end


# --------------------------------------------------------------------------------------------------
# Get lensing quantities for the best-fit model
# --------------------------------------------------------------------------------------------------
"""
    get_potential(file_name::String, θx::T, θy::T; unit::Symbol=:RA_DEC) where T <: Union{ROA, Vector{Int64}}
Calculate the lensing potential at the given coordinates (θx, θy) for the best-fit model based on 
the MCMC results stored in `file_name`. The user can specify the unit of the input coordinates as 
either RA/DEC or arcseconds.

# Arguments
- `file_name::String`: Path to the JLD2 file containing the MCMC results
- `θx`: x-coordinates
- `θy`: y-coordinates

# Keyword Arguments
- `unit::Symbol=:RA_DEC`: Unit of the input coordinates. Supported units are `:RA_DEC` and `:arcsec`. 
   If `:RA_DEC` then (θx, θy) = (RA, DEC) is assumed and will be converted to arcseconds relative to 
   the reference position in the model.

# Returns
- `ψ`: Lensing potential at the input coordinates for the best-fit model.
"""
function get_potential(file_name::String, θx::T, θy::T; unit::Symbol=:RA_DEC) where T <: Union{ROA, Vector{Int64}}
   # Check if the input coordinates are of the same size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same size."))
   end

   # Load the model, chains and chi2 from the input file
   data = jldopen(file_name, "r")
   model  = data["model"]
   chains = data["chains"]
   chi2   = data["chi2"]
   close(data)

   # Get best-fit lens model
   best_model = get_best_model(model, chains, chi2)

   # Convert input coordinates to arcseconds if they are in RA/DEC
   if unit == :RA_DEC
      # Get reference position and pixel scale from the model
      RA_REF = model.observation.reference[1]
      DEC_REF = model.observation.reference[2]

      # Convert RA/DEC to arcseconds relative to the reference position
      θx_arcsec, θy_arcsec = AstrometricOps.gnomonic_offsets_arcsec(RA_REF, DEC_REF, θx, θy)
      return Lenses.get_potential(best_model, θx_arcsec, θy_arcsec)
   elseif unit == :arcsec
      return Lenses.get_potential(best_model, θx, θy)
   else
      throw(ArgumentError("Invalid unit. Supported units are :RA_DEC and :arcsec."))
   end
end


"""
    get_deflection(file_name::String, θx::T, θy::T; unit::Symbol=:RA_DEC) where T <: Union{ROA, Vector{Int64}}
Calculate the lensing deflection at the given coordinates (θx, θy) for the best-fit model based on 
the MCMC results stored in `file_name`. The user can specify the unit of the input coordinates as 
either RA/DEC or arcseconds.

# Arguments
- `file_name::String`: Path to the JLD2 file containing the MCMC results
- `θx`: x-coordinates
- `θy`: y-coordinates

# Keyword Arguments
- `unit::Symbol=:RA_DEC`: Unit of the input coordinates. Supported units are `:RA_DEC` and `:arcsec`. 
   If `:RA_DEC` then (θx, θy) = (RA, DEC) is assumed and will be converted to arcseconds relative to 
   the reference position in the model.

# Returns
- `αx`: x-component of the deflection angle (in arcseconds).
- `αy`: y-component of the deflection angle (in arcseconds).
"""
function get_deflection(file_name::String, θx::T, θy::T; unit::Symbol=:RA_DEC) where T <: Union{ROA, Vector{Int64}}
   # Check if the input coordinates are of the same size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same size."))
   end

   # Load the model, chains and chi2 from the input file
   data = jldopen(file_name, "r")
   model  = data["model"]
   chains = data["chains"]
   chi2   = data["chi2"]
   close(data)

   # Get best-fit lens model
   best_model = get_best_model(model, chains, chi2)

   # Convert input coordinates to arcseconds if they are in RA/DEC
   if unit == :RA_DEC
      # Get reference position and pixel scale from the model
      RA_REF = model.observation.reference[1]
      DEC_REF = model.observation.reference[2]

      # Convert RA/DEC to arcseconds relative to the reference position
      θx_arcsec, θy_arcsec = AstrometricOps.gnomonic_offsets_arcsec(RA_REF, DEC_REF, θx, θy)
      return Lenses.get_deflection(best_model, θx_arcsec, θy_arcsec)
   elseif unit == :arcsec
      return Lenses.get_deflection(best_model, θx, θy)
   else
      throw(ArgumentError("Invalid unit. Supported units are :RA_DEC and :arcsec."))
   end
end


"""
    get_jacobian(file_name::String, θx::T, θy::T; unit::Symbol=:RA_DEC) where T <: Union{ROA, Vector{Int64}}
Calculate the Jacobian (i.e., deformation tensor) of the lens mapping at the given coordinates 
(θx, θy) for the best-fit model based on the MCMC results stored in `file_name`. The user can 
specify the unit of the input coordinates as either RA/DEC or arcseconds.

# Arguments
- `file_name::String`: Path to the JLD2 file containing the MCMC results
- `θx`: x-coordinates
- `θy`: y-coordinates

# Keyword Arguments
- `unit::Symbol=:RA_DEC`: Unit of the input coordinates. Supported units are `:RA_DEC` and `:arcsec`.
   If `:RA_DEC` then (θx, θy) = (RA, DEC) is assumed and will be converted to arcseconds relative to 
   the reference position in the model.

# Returns
- `ψxx`: xx-component of the jacobian.
- `ψyy`: yy-component of the jacobian.
- `ψxy`: xy-component of the jacobian.
"""
function get_jacobian(file_name::String, θx::T, θy::T; unit::Symbol=:RA_DEC) where T <: Union{ROA, Vector{Int64}}
   # Check if the input coordinates are of the same size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same size."))
   end

   # Load the model, chains and chi2 from the input file
   data = jldopen(file_name, "r")
   model  = data["model"]
   chains = data["chains"]
   chi2   = data["chi2"]
   close(data)

   # Get best-fit lens model
   best_model = get_best_model(model, chains, chi2)

   # Convert input coordinates to arcseconds if they are in RA/DEC
   if unit == :RA_DEC
      # Get reference position and pixel scale from the model
      RA_REF = model.observation.reference[1]
      DEC_REF = model.observation.reference[2]

      # Convert RA/DEC to arcseconds relative to the reference position
      θx_arcsec, θy_arcsec = AstrometricOps.gnomonic_offsets_arcsec(RA_REF, DEC_REF, θx, θy)
      return Lenses.get_jacobian(best_model, θx_arcsec, θy_arcsec)
   elseif unit == :arcsec
      return Lenses.get_jacobian(best_model, θx, θy)
   else
      throw(ArgumentError("Invalid unit. Supported units are :RA_DEC and :arcsec."))
   end
end


# --------------------------------------------------------------------------------------------------
# Save best-fit lensing quantities to FITS files
# --------------------------------------------------------------------------------------------------
"""
    save_best_fit(file_name::String; save_potential::Bool=true, save_deflection::Bool=true, save_kappa::Bool=true, save_gamma::Bool=true)
Save fits files for the best-fit model based on the MCMC results stored in `file_name`. The user can
choose which lensing quantities to save by setting the corresponding boolean flags.

# Arguments
- `file_name::String`: Path to the JLD2 file containing the MCMC results

# Keyword Arguments
- `save_potential::Bool=true`: Whether to save the lensing potential as **potential.fits**
- `save_deflection::Bool=true`: Whether to save the deflection angles as **alpha\\_x.fits** and **alpha\\_y.fits**
- `save_kappa::Bool=true`: Whether to save the convergence map as **kappa.fits**
- `save_gamma::Bool=true`: Whether to save the shear maps as **gamma1.fits** and **gamma2.fits**

# Returns
- Nothing (saves FITS files to disk)
"""
function save_best_fit(file_name::String;
                      save_potential::Bool = true,
                      save_deflection::Bool = true,
                      save_kappa::Bool = true, save_gamma::Bool = true)

   # Load the chains and chi2 from the file
   data = jldopen(file_name, "r")
   model  = data["model"]
   chains = data["chains"]
   chi2   = data["chi2"]
   date   = data["Date"]
   time   = data["Time"]
   close(data)

   # Get best-fit model
   best_model = get_best_model(model, chains, chi2)

   # Generate grid
   FOV = model.observation.FOV
   pixel_scale = model.observation.pixel_scale
   x_grid, y_grid = Lenses.get_meshgrid(0.5 * FOV[1], 0.5 * FOV[2], pixel_scale)

   # Generate FITS file header
   if save_potential
      # Calculate potential
      ψ = Lenses.get_potential(best_model, x_grid, y_grid)
   
      # Open new FITS file
      f = FITS("./potential.fits", "w")
      write(f, ψ)
   
      # Write header
      hdu = f[1]
      _write_fits_header!(hdu, model; date, time)
      close(f)
   end
   

   # Calculate deflection angles
   if save_deflection
      # Calculate deflection map
      ψx, ψy = Lenses.get_deflection(best_model, x_grid, y_grid)

      # Open new FITS file
      f = FITS("./alpha_x.fits", "w")
      write(f, ψx)
      
      hdu = f[1]
      _write_fits_header!(hdu, model; date, time)
      close(f)

      f = FITS("./alpha_y.fits", "w")
      write(f, ψy)
      
      hdu = f[1]
      _write_fits_header!(hdu, model; date, time)      
      close(f)
   end

   # Calculate magnification
   if save_kappa || save_gamma
      ψxx, ψyy, ψxy = Lenses.get_jacobian(best_model, x_grid, y_grid)

      # Calculate and save convergence
      if save_kappa
         # Open new FITS file
         f = FITS("./kappa.fits", "w")
         write(f, @. 0.5 * (ψxx + ψyy))
      
         hdu = f[1]
         _write_fits_header!(hdu, model; date, time)      
         close(f)
      end
      
      # Calculate and save shear components
      if save_gamma
         # Open new FITS file
         f = FITS("./gamma1.fits", "w")
         write(f, @. 0.5 * (ψxx - ψyy))
      
         hdu = f[1]
         _write_fits_header!(hdu, model; date, time)      
         close(f)

         f = FITS("./gamma2.fits", "w")
         write(f, ψxy)
      
         hdu = f[1]
         _write_fits_header!(hdu, model; date, time)      
         close(f)
      end
   end
   return nothing
end

end