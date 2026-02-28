module LensModel


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using FITSIO
using JLD2
using Dates

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

# Include the Nelder-Mead module
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
export read_input, fit_model
export free_parameter_names
export calculate_gr
export print_gr_report
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
# Write FITS header
# --------------------------------------------------------------------------------------------------
function write_fits_header!(header::ImageHDU, model::ModelConfig, grid_x::Matrix{<:RV}, grid_y::Matrix{<:RV};
                            date::Union{Nothing, Date} = nothing,
                            time::Union{Nothing, Time} = nothing)
   # Add coordinate projection type
   write_key(header, "CTYPE1", "RA---TAN", "RA coordinate type")
   write_key(header, "CTYPE2", "DEC--TAN", "DEC coordinate type")

   # Add reference position
   RA_REF = model.observation.reference[1]
   DEC_REF = model.observation.reference[2]
   if RA_REF == 0.0 && DEC_REF == 0.0
      @warn "Reference position is (0.0, 0.0). Are you sure?"
   end
   write_key(header, "CRVAL1",  RA_REF, "RA reference value")
   write_key(header, "CRVAL2", DEC_REF, "DEC reference value")
   
   # Add reference pixel
   PIXEL_SCALE = model.observation.pixel_scale
   CRPIX1 = 0.5 * model.observation.FOV[1] / PIXEL_SCALE + 1.0
   CRPIX2 = 0.5 * model.observation.FOV[2] / PIXEL_SCALE + 1.0
   write_key(header, "CRPIX1", CRPIX1, "Reference pixel in x-direction")
   write_key(header, "CRPIX2", CRPIX2, "Reference pixel in y-direction")
   
   # Add pixel scale
   write_key(header, "CDELT1", -PIXEL_SCALE / 3600.0, "Pixel scale in RA (degrees)")
   write_key(header, "CDELT2", +PIXEL_SCALE / 3600.0, "Pixel scale in DEC (degrees)")
   
   # Comments
   write_key(header, "MODELER", model.observation.modeler, "Modeler name")
   write_key(header, "LENS", model.observation.lens, "Lens name")
   write_key(header, "Z_D", model.observation.z_d, "Lens redshift")
   if date !== nothing
      write_key(header, "DATE", string(date), "Date of fit (UTC)")
   end
   if time !== nothing
      write_key(header, "TIME", string(time), "Time of fit (UTC)")
   end
   
end

function save_best_fit(model::ModelConfig, chains::Array{Float64, 3}, chi2::Matrix{Float64};
                      date::Union{Nothing, Date} = nothing,
                      time::Union{Nothing, Time} = nothing,
                      save_potential::Bool = true, save_deflection::Bool = true, 
                      save_kappa::Bool = true, save_gamma::Bool = true)
   # Get best fit parameters
   best_fit, _, _ = get_best_fit(chi2, chains)
   
   # Get list of parameters for the lens model
   param_ref = Dict(p.key => p.refer for p in model.parameters)
   
   # Replace free parameter values by best-fit values
   pvals = param_dict(model, best_fit, param_ref)

   # Get best-fit model
   best_model = build_lens(model, pvals)

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
      write_fits_header!(hdu, model, x_grid, y_grid; date, time)
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
      write_fits_header!(hdu, model, x_grid, y_grid; date, time)
      close(f)

      f = FITS("./alpha_y.fits", "w")
      write(f, ψy)
      
      hdu = f[1]
      write_fits_header!(hdu, model, x_grid, y_grid; date, time)      
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
         write_fits_header!(hdu, model, x_grid, y_grid; date, time)      
         close(f)
      end
      
      # Calculate and save shear components
      if save_gamma
         # Open new FITS file
         f = FITS("./gamma1.fits", "w")
         write(f, @. 0.5 * (ψxx - ψyy))
      
         hdu = f[1]
         write_fits_header!(hdu, model, x_grid, y_grid; date, time)      
         close(f)

         f = FITS("./gamma2.fits", "w")
         write(f, ψxy)
      
         hdu = f[1]
         write_fits_header!(hdu, model, x_grid, y_grid; date, time)      
         close(f)
      end
   end
   return nothing
end

function save_best_fit(file_name::String;                      
                      save_potential::Bool = true, save_deflection::Bool = true, 
                      save_kappa::Bool = true, save_gamma::Bool = true)

   # Load the chains and chi2 from the file
   data = jldopen(file_name, "r")
   model  = data["model"]
   chains = data["chains"]
   chi2   = data["chi2"]
   date   = data["Date"]
   time   = data["Time"]
   close(data)

   # Save the best fit
   save_best_fit(model, chains, chi2; 
                 date = date,
                 time = time,
                 save_potential=save_potential, 
                 save_deflection=save_deflection, 
                 save_kappa=save_kappa, 
                 save_gamma=save_gamma)
end


end