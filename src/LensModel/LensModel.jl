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


end