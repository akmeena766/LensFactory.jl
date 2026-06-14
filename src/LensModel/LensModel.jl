module LensModel


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using Printf
using StatsBase
using JLD2
using Dates
using FITSIO


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
include("./LensModelIO.jl")
using .LensModelIO

include("./Likelihood.jl")
using .Likelihood

include("./LensModelUtils.jl")
using .LensModelUtils

include("./LensModelFit.jl")
using .LensModelFit

include("./Diagnostic.jl")
using .Diagnostic


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export read_input
export fit_model
export free_parameter_names
export get_best_fit_parameters
export get_best_fit_rms
export get_best_model
export get_potential
export get_deflection
export get_jacobian
export predict_image
export check_parity
export save_best_fits


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
# Read input file and return model configuration
# --------------------------------------------------------------------------------------------------
"""
    read_input(input_filename::AbstractString)
Reads the input YAML file and constructs a `ModelConfig` struct containing all the necessary 
information for lens modeling and sampling. For details on the expected structure of the input YAML 
file, please refer to [Example - 2](https://github.com/akmeena766/LensFactory_Examples/blob/).

# Arguments
- `input_filename::AbstractString`: Path to the input YAML file.

# Returns
- `ModelConfig::ModelConfig`: A struct containing the observation details, cosmology, lens 
   configuration, source configuration, parameter definitions, and sampling configuration.
"""
function read_input(input_filename::String)
   return LensModelIO._read_input(input_filename)
end


# --------------------------------------------------------------------------------------------------
# Fit lens model
# --------------------------------------------------------------------------------------------------
"""
    fit_model(model::ModelConfig; save::Bool=true, file_name::Union{String, Nothing}=nothing)
Performs lens model fitting using the given `model`.

# Arguments
- `model::ModelConfig`: The lens model configuration containing the observation, lens, source, and sampler details.
- `save::Bool=true`: Whether to save the MCMC results to a JLD2 file (default: true)
   - `file_name::Union{String, Nothing}=nothing`: The name of the JLD2 file to save results. If no 
      name is provided then a name will be generated based on the lens name and current date (i.e., 
      LensName_DDMMYYYY.jld2).

# Returns
- `chains::Array{Float64, 3}`: The MCMC chains containing sampled parameter values
- `log_posterior::Matrix{Float64}`: The log-posterior values corresponding to each sample in the chains
"""
function fit_model(model::ModelConfig; save::Bool=true, file_name::Union{String, Nothing}=nothing)
   return LensModelFit._fit_model(model, save=save, file_name=file_name)
end


# --------------------------------------------------------------------------------------------------
# Get free parameter names
# --------------------------------------------------------------------------------------------------
"""
    free_parameter_names(model::ModelConfig)
Get a vector of tuples containing the owner and name of each free parameter in the model.

# Arguments
- `model::ModelConfig`: The lens model configuration containing the observation, lens, source, and sampler details.

# Returns
- `Vector{Tuple{Symbol, Symbol}}`: A vector of tuples containing the owner and name of each free parameter.
"""
function free_parameter_names(model::ModelConfig)
   return [(p.owner, p.name) for p in model.parameters[model.free_param_idxs]]
end


# --------------------------------------------------------------------------------------------------
# Get best-fit parameters from the full optimization/MCMC results
# --------------------------------------------------------------------------------------------------
"""
    get_best_fit_parameters(results; chains=nothing, burn_in::Float64=0.2, with_errors::Bool=false, thin::Int64=1, print_table::Bool=false)
Get the best-fit parameters from the optimization or MCMC results.

# Arguments
- `results::Union{Vector{@NamedTuple{θ::Vector{Float64}, f::Float64}}, Matrix{Float64}}`: The optimization or MCMC results.
   - If `results::Vector{@NamedTuple{θ::Vector{Float64}, f::Float64}}`: It is assumed to be from the 
      optimizer results, where each element is a named tuple containing the parameter vector (`θ`) 
      and the χ² values (`f`).
   - If `results::Matrix{Float64}`: It is assumed to be the χ² matrix corresponding to MCMC chains.

# Keyword Arguments
- `chains::Union{Array{Float64, 3}, Nothing}=nothing`: MCMC chains. The dimensions should be 
   `(n_steps, n_chains, n_parameters)`. This is only required when `results` is a `Matrix{Float64}`.
- `burn_in::Float64=0.2`: The fraction of MCMC chains to discard as burn-in (default: 0.2).
- `with_errors::Bool=false`: If `true`, then also return 1-sigma uncertainties for best-fit parameters.
   But it will only work if `chains` are provided.
- `thin::Int64=1`: The thinning factor to apply to the MCMC chains while calculating parameter 
   uncertainties.
- `print_table::Bool=false`: If `true`, then print the best-fit parameters along with [-1σ, +1σ] and
   [-2σ, +2σ] ranges in a table format. This option is only available when `with_errors=true`.

# Returns
- `best_θ::Vector{Float64}`: Best-fit parameter values.
- `log_posterior::Float64`: Log-posterior value corresponding to the best-fit parameters.
- `lower_err::Vector{Float64}`: Lower 1-sigma uncertainties (only if `with_errors=true`).
- `upper_err::Vector{Float64}`: Upper 1-sigma uncertainties (only if `with_errors=true`).
"""
function get_best_fit_parameters(results::Union{Vector{@NamedTuple{θ::Vector{Float64}, f::Float64}}, Matrix{Float64}}; 
                             chains::Union{Array{Float64, 3}, Nothing}=nothing, 
                             burn_in::Float64=0.2,
                             with_errors::Bool=false,
                             thin::Int64=1,
                             print_table::Bool=false,
                             free_parameter_names=nothing)
   # Initialize outputs
   best_θ    = nothing
   best_logL = -Inf

   # Case A: Input is from the Parallel Optimizer (Vector of NamedTuples)
   if results isa Vector{@NamedTuple{θ::Vector{Float64}, f::Float64}}
      # We already sorted the results in descending order of logL, so results[1] is the best fit
      best_run  = results[1]
      best_θ    = best_run.θ
      best_logL = best_run.f

   # Case B: Input is from MCMC (χ² matrix)
   elseif results isa Matrix{Float64} && chains !== nothing
      # Get chain details 
      n_steps, n_chains, n_params = size(chains)

      # Remove Burn-in
      start_idx      = Int(floor(n_steps * burn_in)) + 1
      thinned_chains = chains[start_idx:thin:end, :, :]
      thinned_logL   = results[start_idx:thin:end, :]

      # Get (step, chain) of the best logL by finding the index of the maximum value in the logL matrix
      best_idx        = argmax(thinned_logL)
      step, chain_num = best_idx[1], best_idx[2]
      
      best_θ    = thinned_chains[step, chain_num, :]
      best_logL = thinned_logL[step, chain_num]

      if with_errors
         # Reshape the thinned chain into a flat array
         flat_chain = reshape(thinned_chains, :, n_params)

         # Calculate (asymmetric) errors as we are defining errors relative to the Best-Fit value
         lower_err  = zeros(n_params)
         upper_err  = zeros(n_params)
         lower_err2 = zeros(n_params)
         upper_err2 = zeros(n_params)
         for i in 1:n_params
            # Get 16th and 84th percentiles of the posterior
            q16, q84 = StatsBase.quantile(flat_chain[:, i], [0.16, 0.84])
            q2p30, q97p7 = StatsBase.quantile(flat_chain[:, i], [0.023, 0.977])
        
            # Asymmetric error: distance from best-fit to the quantiles
            lower_err[i] = q16 - best_θ[i]
            upper_err[i] = q84 - best_θ[i]

            lower_err2[i] = q2p30 - best_θ[i]
            upper_err2[i] = q97p7 - best_θ[i]
         end

         if print_table
            # Define header
            header = @sprintf("| %-8s | %-10s | %-12s |   %-18s |   %-18s |", "Owner", "parameter", "Best-fit", "[-1σ, +1σ]", "[-2σ, +2σ]")
            
            # Print header
            println("-"^length(header))
            println(header)
            println("-"^length(header))

            # Run over parameters
            for i in 1:n_params
               owner = string(free_parameter_names[i][1])
               child = string(free_parameter_names[i][2])
               
               @printf("| %-8s | %-10s | %-12.4f | [%+6.3f, %+6.3f]%-4s | [%+6.3f, %+6.3f]%-4s |\n", 
                     owner, child, best_θ[i], lower_err[i], upper_err[i], "", lower_err2[i], upper_err2[i],"")
            end
            println("-"^length(header))
         end
         return best_θ, best_logL, lower_err, upper_err
      end
      return best_θ, best_logL     
   else
      error("Please provide either Optimization results or both MCMC logL and chains.")
   end
end


# --------------------------------------------------------------------------------------------------
# Calculate RMS for the best-fit model
# --------------------------------------------------------------------------------------------------
"""
    get_best_fit_rms(model::ModelConfig, chains::Array{Float64, 3}, logL::Matrix{Float64}; check_parity::Bool=false, burn_in::Float64=0.2)
Calculate the RMS of the best-fit model based on the MCMC results stored in `chains` and `logL`. The
function prints a table showing the RMS for each knot image, as well as a global total RMS. If 
`check_parity` is set to true, the function will also check the parity of each knot image and 
print a warning if any parity does not match.

# Arguments
- `model::ModelConfig`: The lens model configuration used for the MCMC fit.
- `chains::Array{Float64, 3}`: The MCMC chains containing the sampled parameter values. The 
   dimensions should be (n_steps, n_chains, n_parameters).
- `logL::Matrix{Float64}`: The logL values corresponding to each sample in the MCMC chains. The 
   dimensions should be (n_steps, n_steps).

# Keyword Arguments
- `check_parity::Bool=false`: Whether to check the parity of each knot image against the input parity.
- `burn_in::Float64=0.2`: The fraction of MCMC chains to discard as burn-in while determining 
   best-fit lens model.

# Returns
- `nothing`: Prints a table to the console with the RMS for each knot image, as well as total RMS.
"""
function get_best_fit_rms(model::ModelConfig, chains::Array{Float64, 3}, logL::Matrix{Float64}; check_parity::Bool=false, burn_in::Float64=0.2)
   # Get the best parameters based on minimum logL
   best_θ, _ = get_best_fit_parameters(logL; chains=chains, burn_in=burn_in)

   # Get list of parameters for the lens model
   param_ref = Dict(p.key => p.refer for p in model.parameters)
   
   # Replace free parameter values by best-fit values
   pvals = LensModelUtils.param_dict(model, best_θ, param_ref)

   # Get best-fit model
   best_model = LensModelUtils.build_lens(model, pvals)

   # Get angular-diameter distance ratios
   adis = LensModelUtils.adis_current(model, pvals)

   # Generate grid
   FOV = model.observation.FOV
   pixel_scale = model.observation.pixel_scale
   x_grid, y_grid = Lenses.get_meshgrid(0.5 * FOV[1], 0.5 * FOV[2], pixel_scale)

   # Calculate deformation at all image positions
   ψ_all, αx_all, αy_all, A_all = LensModelUtils.lens_quantities(model, best_model)

   # Identity tuple
   I4 = (1.0, 0.0, 0.0, 1.0)

   # Global RMS and image count variables
   global_sq_dist = 0.0
   global_count = 0

   # Helper for column padding
   col(txt, width) = rpad(string(txt), width)

   # Print Header
   # Table Header
   header_line = "-"^72
   knot_sep  = "             " * "-"^59
   println(header_line)
   println("| ", col("Source", 10), 
          " | ", col("Knot", 10), 
          " | ", col("Img", 6), 
          " | ", col("Image Dist (arcsec)", 20), 
          " | ", col("Knot RMS", 10), 
          " |")
   println(header_line)

   # Calculate RMS for each image
   sid = 1
   kid = 1
   for src in model.source_config.sources
      # Get angular-diameter distance ratio for this source
      adis_value = adis[sid]
         
      # Generate source id
      src_id = Symbol(:src, sid)
      src_knot = 0
      for knot in src.knots
         # Generate knot id
         knot_id = Symbol(:knot, kid)

         # Update knot counter
         src_knot = src_knot + 1

         # Knot positions and measurement errors
         x  = knot.x
         y  = knot.y
         σx = knot.σx
         σy = knot.σy
         σθ = knot.σθ
         
         # Number of images for this knot
         n = length(x)

         # Deflection vector at the knot positions
         αx = @. adis_value * αx_all[kid]
         αy = @. adis_value * αy_all[kid]

         # Deformation tensor at the knot positions
         A = @. adis_value * A_all[kid]
         for i in eachindex(A)
            @. A[i] = I4[i] - A[i]
         end

         # Individual source positions using broadcasting
         βx_ind = @. x - αx
         βy_ind = @. y - αy

         # Get weighted source position (Section 4.1 in https://arxiv.org/pdf/astro-ph/0102340)
         βx_model, βy_model, _ = Likelihood._weighted_position(βx_ind, βy_ind, A, σx, σy, σθ, n)

         # Get image positions
         predicted_image = Lenses.get_image(best_model, x_grid, y_grid, adis_value, (βx_model, βy_model))

         # Convert predicted to mutable arrays for iterative removal
         pred_x = Float64[p[1] for p in predicted_image]
         pred_y = Float64[p[2] for p in predicted_image]

         # Knot counters
         knot_sq_dist = 0.0
         knot_count = 0

         # Store distances to print after matching is done for the whole knot
         results = []
         
         # Matching observed images to predicted images
         for i in 1:n
            if isempty(pred_x)
               push!(results, "MISSING")
               continue
            end

            # Calculate distances to all remaining candidates
            dx = @. pred_x .- x[i]
            dy = @. pred_y .- y[i]
            dist_sq = @. dx^2 + dy^2

            # Find the closest predicted image index
            best_idx = argmin(dist_sq)

            d2 = dist_sq[best_idx]
            dist = sqrt(d2)

            # Update global and knot totals
            global_sq_dist += d2
            global_count += 1
            knot_sq_dist += d2
            knot_count += 1

            push!(results, round(dist, digits=6))

            # Remove this candidate so it can't be matched twice
            deleteat!(pred_x, best_idx)
            deleteat!(pred_y, best_idx)
         end
         
         # Calculate RMS for this specific knot
         k_rms = knot_count > 0 ? round(sqrt(knot_sq_dist / knot_count), digits=6) : "N/A"

         # Print image-by-image breakdown
         for (i, d_val) in enumerate(results)
            k_display = (i == 1) ? string(k_rms) : ""
            println("| ", col("src$sid", 10), 
                   " | ", col("knot$src_knot", 10), 
                   " | ", col(i, 6), 
                   " | ", col(d_val, 20), 
                   " | ", col(k_display, 10), 
                   " |"
            )
         end
         println("-"^72)
         kid = kid + 1
      end
      sid = sid + 1
   end

   final_rms = global_count > 0 ? sqrt(global_sq_dist / global_count) : 0.0
   println("| GLOBAL TOTAL RMS: ", col(round(final_rms, digits=6), 50), " |")
   println("-"^72)
   
   return nothing
end

# --------------------------------------------------------------------------------------------------
# Get lensing quantities for the best-fit model
# --------------------------------------------------------------------------------------------------
"""
    get_best_model(model::ModelConfig, chains::Array{Float64, 3}, logL::Matrix{Float64})
Get the best-fit lens model based on the MCMC results stored in `chains` and `logL`. The best-fit 
parameters are determined by the minimum log-likelihood in `logL`.

# Arguments
- `model::ModelConfig`: The lens model configuration used for the MCMC fit.
- `chains::Array{Float64, 3}`: The MCMC chains containing the sampled parameter values. The 
   dimensions should be (n_steps, n_chains, n_parameters).
- `logL::Matrix{Float64}`: The log-likelihood values corresponding to each sample in the MCMC chains.
   The dimensions should be (n_steps, n_steps).

# Returns
- `best_model`: Best-fit lens model constructed using the best-fit parameters.
- `logL`: Log-likelihoohood correspodning to best-fit lens model.
"""
function get_best_model(model::ModelConfig; 
                        optim_result::Union{Vector{@NamedTuple{θ::Vector{Float64}, f::Float64}}, Nothing}=nothing, 
                        mcmc_chains::Union{Array{Float64, 3}, Nothing}=nothing, 
                        mcmc_logL::Union{Matrix{Float64}, Nothing}=nothing, 
                        burn_in::Float64=0.2)
   # Get the best parameters based on minimum log-likelihood
   if optim_result !== nothing
      best_θ, best_logL = get_best_fit_parameters(optim_result)
   elseif mcmc_chains !== nothing && mcmc_logL !== nothing
      best_θ, best_logL = get_best_fit_parameters(mcmc_logL; chains=mcmc_chains, burn_in=burn_in)
   else
      error("Either optim_result or (mcmc_chains and mcmc_logL) must be provided.")
   end

   # Get list of parameters for the lens model
   param_ref = Dict(p.key => p.refer for p in model.parameters)
   
   # Replace free parameter values by best-fit values
   pvals = LensModelUtils.param_dict(model, best_θ, param_ref)
   
   return LensModelUtils.build_lens(model, pvals), best_logL
end


"""
    get_potential(best_model::Lenses.AbstractLens, θx::T, θy::T;
                  reference::Tuple{Float64, Float64}) where T <: Union{RV, ROA}
"""
function get_potential(best_model::Lenses.AbstractLens, θx::T, θy::T; 
                       reference::Tuple{Float64, Float64} = (0.0, 0.0)) where T <: Union{RV, ROA}
   # Check if the input coordinates are of the same size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same size."))
   end
   
   if reference[1] == 0.0 && reference[2] == 0.0
      return Lenses.get_potential(best_model, θx, θy)
   else
      θx_arcsec, θy_arcsec = AstrometricOps.gnomonic_offsets_arcsec(reference[1], reference[2], θx, θy)
      return Lenses.get_potential(best_model, θx_arcsec, θy_arcsec)
   end
end


"""
    get_potential(jld2_file::JLD2.JLDFile, θx::T, θy::T; unit::Symbol=:RA_DEC) where T <: Union{ROA, Vector{Int64}}
Calculate the lensing potential at the given coordinates `(θx, θy)` for the best-fit model. The user
can specify the unit of the input coordinates as either RA/DEC or arcseconds.

# Arguments
- `jld2_file::JLD2.JLDFile` or `best_model::Lenses.AbstractLens`: JLD2 data file containing the fit
   results or the best-fit lens model.
- `θx`: x-coordinates
- `θy`: y-coordinates

# Returns
- `ψ`: Lensing potential at the input coordinates for the best-fit model.
"""
function get_potential(data_jld2::JLD2.JLDFile, θx::T, θy::T) where T <: Union{RV, ROA}
   # Check if the input coordinates are of the same size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same size."))
   end

   # Load the model, chains and logL from the input file
   model  = data_jld2["model"]
   chains = data_jld2["chains"]
   logL   = data_jld2["logL"]

   # Get best-fit lens model
   best_model, _ = get_best_model(model; mcmc_chains=chains, mcmc_logL=logL)

   return get_potential(best_model, θx, θy::T; reference=(model.observation.RA_REF, model.observation.DEC_REF))
end


"""
    get_deflection(best_model::Lenses.AbstractLens, θx::T, θy::T; unit::Symbol=:RA_DEC) where T <: Union{RV, ROA}
"""
function get_deflection(best_model::Lenses.AbstractLens, θx::T, θy::T; 
                        reference::Tuple{Float64, Float64} = (0.0, 0.0)) where T <: Union{RV, ROA}
   # Check if the input coordinates are of the same size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same size."))
   end

   # Convert input coordinates to arcseconds if they are in RA/DEC
   if reference[1] == 0.0 && reference[2] == 0.0
      return Lenses.get_deflection(best_model, θx, θy)
   else
      θx_arcsec, θy_arcsec = AstrometricOps.gnomonic_offsets_arcsec(reference[1], reference[2], θx, θy)
      return Lenses.get_deflection(best_model, θx_arcsec, θy_arcsec)
   end
end


"""
    get_deflection(jld2_file::JLD2.JLDFile, θx::T, θy::T; unit::Symbol=:RA_DEC) where T <: Union{RV, ROA}
Calculate the lensing deflection at the given coordinates `(θx, θy)` for the best-fit model. The
user can specify the unit of the input coordinates as either RA/DEC or arcseconds.

# Arguments
- `jld2_file::JLD2.JLDFile` or `best_model::Lenses.AbstractLens`: JLD2 data file containing the fit
   results or the best-fit lens model.
- `θx`: x-coordinates
- `θy`: y-coordinates

# Keyword Arguments
- reference = (0.0, 0.0): Reference (RA, Dec) coordinate (in degrees). 

# Returns
- `αx`: x-component of the deflection angle (in arcseconds).
- `αy`: y-component of the deflection angle (in arcseconds).
"""
function get_deflection(data_jld2::JLD2.JLDFile, θx::T, θy::T; unit::Symbol=:RA_DEC) where T <: Union{ROA, Vector{Int64}}
   # Check if the input coordinates are of the same size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same size."))
   end

   # Load the model, chains and logL from the input file
   model  = data_jld2["model"]
   chains = data_jld2["chains"]
   logL   = data_jld2["logL"]

   # Get best-fit lens model
   best_model, _ = get_best_model(model; mcmc_chains=chains, mcmc_logL=logL)

   return get_deflection(best_model, θx::T, θy::T; reference=(model.observation.RA_REF, model.observation.DEC_REF))
end

"""
    get_jacobian(best_model::Lenses.AbstractLens, θx::T, θy::T; unit::Symbol=:RA_DEC) where T <: Union{RV, ROA}
"""
function get_jacobian(best_model::Lenses.AbstractLens, θx::T, θy::T; 
                      reference::Tuple{Float64, Float64} = (0.0, 0.0)) where T <: Union{RV, ROA}
   # Check if the input coordinates are of the same size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same size."))
   end

   # Convert input coordinates to arcseconds if they are in RA/DEC
   if reference[1] == 0.0 && reference[2] == 0.0
      return Lenses.get_jacobian(best_model, θx, θy)
   else
      # Convert RA/DEC to arcseconds relative to the reference position
      θx_arcsec, θy_arcsec = AstrometricOps.gnomonic_offsets_arcsec(reference[1], reference[2], θx, θy)
      return Lenses.get_jacobian(best_model, θx_arcsec, θy_arcsec)
   end
end


"""
    get_jacobian(jld2_file::JLD2.JLDFile, θx::T, θy::T; unit::Symbol=:RA_DEC) where T <: Union{RV, ROA}
Calculate the Jacobian (i.e., deformation tensor) at the given coordinates `(θx, θy)` for the 
best-fit model. The user can specify the unit of the input coordinates as either RA/DEC or arcseconds.

# Arguments
- `jld2_file::JLD2.JLDFile` or `best_model::Lenses.AbstractLens`: JLD2 data file containing the fit
   results or the best-fit lens model.
- `θx`: x-coordinates
- `θy`: y-coordinates

# Keyword Arguments
- reference = (0.0, 0.0): Reference (RA, Dec) coordinate (in degrees). 

# Returns
- `ψxx`: xx-component of the jacobian.
- `ψyy`: yy-component of the jacobian.
- `ψxy`: xy-component of the jacobian.
"""
function get_jacobian(data_jld2::JLD2.JLDFile, θx::T, θy::T) where T <: Union{RV, ROA}
   # Check if the input coordinates are of the same size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same size."))
   end

   # Load the model, chains and logL from the input file
   model  = data_jld2["model"]
   chains = data_jld2["chains"]
   logL   = data_jld2["logL"]

   # Get best-fit lens model
   best_model, _ = get_best_model(model; mcmc_chains=chains, mcmc_logL=logL)

   return get_jacobian(best_model, θx::T, θy::T; reference=(model.observation.RA_REF, model.observation.DEC_REF))
end


"""
    predict_image(jld2_file::JLD2.JLDFile, θx::T, θy::T, z_s::Float64; unit::Symbol=:RA_DEC) where T <: Union{RV, Vector{Float64}}
Predict counter-image positions, magnifications and time delays based on the best-fit lens model.
The function can take either a single observed image or multiple observed images of the same system.
If multiple images are provided then the function will calculate the barycentric source position and
then predict the counter-image positions, magnifications and time delays.

# Arguments
- `data_jld2::JLD2.JLDFile`: JLD2 data file containing the fit results.
- `θx`: x-coordinate(s) of the observed image position(s).
- `θy`: y-coordinate(s) of the observed image position(s).
- `z_s::Float64`: Source redshift.

# Keyword Arguments
- `unit::Symbol=:RA_DEC`: Unit of the input coordinates. 
   - `:RA_DEC`: (θx, θy) are assumed to be in RA/DEC (in degrees).
   - `:arcsec`: (θx, θy) are assumed to be in arcseconds.

# Returns
- `nothing`: Prints the counter-image positions, magnifications and time delays in a table format.
   The table is sorted based on the time delay and at the end also contains the best-fit model 
   predicted source position.
"""
function predict_image(data_jld2::JLD2.JLDFile, θx::T, θy::T, z_s::Float64; unit::Symbol=:RA_DEC) where T <: Union{RV, Vector{Float64}}
   # Check if the input coordinates are of the same size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same size."))
   end

   # Load the model, chains and logL from the input file
   model  = data_jld2["model"]
   chains = data_jld2["chains"]
   logL   = data_jld2["logL"]

   # Get best-fit lens model
   best_model, _ = get_best_model(model; mcmc_chains=chains, mcmc_logL=logL)

   # Construct grid
   FOV = model.observation.FOV
   pixel_scale = model.observation.pixel_scale
   x_grid, y_grid = Lenses.get_meshgrid(0.5 * FOV[1], 0.5 * FOV[2], pixel_scale)

      # Get cosmology
   cosmo = model.cosmology
      
   # ADDs
   z_d = model.observation.z_d
   Dol = Cosmology.angular_diameter_distance(cosmo, 0.0, z_d)
   Dls = Cosmology.angular_diameter_distance(cosmo, z_d, z_s)
   Dos = Cosmology.angular_diameter_distance(cosmo, 0.0, z_s)
   adis = Dls / Dos

   # Get reference position and pixel scale from the model
   RA_REF = model.observation.reference[1]
   DEC_REF = model.observation.reference[2]

   # Convert input coordinates to arcseconds if they are in RA/DEC
   if unit == :RA_DEC
      # Convert RA/DEC to arcseconds relative to the reference position
      θx_arcsec, θy_arcsec = AstrometricOps.gnomonic_offsets_arcsec(RA_REF, DEC_REF, θx, θy)
   elseif unit == :arcsec
      θx_arcsec, θy_arcsec = θx, θy
   else
      throw(ArgumentError("Invalid unit. Supported units are :RA_DEC and :arcsec."))
   end

   # 
   if size(θx, 1) > 1
      # Get deflection at the image positions
      αx, αy = Lenses.get_deflection(best_model, θx_arcsec, θy_arcsec)
      
      # Get magnification at the image positions
      μ_obs = Lenses.get_magnification_image(best_model, θx_arcsec, θy_arcsec, adis)
      
      # Calculate individual image source positions
      βx_indi = θx_arcsec - adis * αx
      βy_indi = θy_arcsec - adis * αy

      # Calculate barycenter source position
      βx_model = sum(βx_indi .* μ_obs.^2) / sum(μ_obs.^2)
      βy_model = sum(βy_indi .* μ_obs.^2) / sum(μ_obs.^2)
   else
      αx, αy = Lenses.get_deflection(best_model, θx_arcsec, θy_arcsec)
      βx_model = θx_arcsec - adis * αx
      βy_model = θy_arcsec - adis * αy
   end

   # Get predicted image positions
   pred_image = Lenses.get_image(best_model, x_grid, y_grid, adis, (βx_model, βy_model))

   # Convert predicted image positions in (RA, DEC) if input is in (RA, DEC)
   if unit == :RA_DEC
      pred_image_RADEC = AstrometricOps.gnomonic_offsets_radec(RA_REF, DEC_REF, first.(pred_image), last.(pred_image))
   end

   # Get magnification at image positions
   mu = Lenses.get_magnification_image(best_model, first.(pred_image), last.(pred_image), adis)

   # Get time delay for image positions (in days)
   td = Lenses.get_time_delay(best_model, first.(pred_image), last.(pred_image), adis, z_d, Dol, (βx_model, βy_model))
   td .= (td .- minimum(td)) ./ DAY2SECOND

   # Stack image position, mu and time delay and sort based on time delay
   if unit == :RA_DEC
      sorted_images = sortslices([first.(pred_image_RADEC) last.(pred_image_RADEC) mu td], dims=1, by=x->x[4])
   elseif unit == :arcsec
      sorted_images = sortslices([first.(pred_image) last.(pred_image) mu td], dims=1, by=x->x[4])
   end

   # Define and print table Header
   if unit == :RA_DEC
      header = @sprintf("| %-5s | %-12s | %-12s | %-10s | %-10s", "Img", "RA", "Dec", "mu (μ)", "Time delay (days) |")
   elseif unit == :arcsec
      header = @sprintf("| %-5s | %-12s | %-12s | %-10s | %-10s", "Img", "x", "y", "mu (μ)", "Time delay (days) |")
   end
   println("-"^length(header))
   println(header)
   println("-"^length(header))

   # Print images and their details
   for i in 1:length(pred_image)
      @printf("| %-5d | %-12.6f | %-12.6f | %-10.4f | %-17.4f |\n", 
              i, sorted_images[i, 1], sorted_images[i, 2], sorted_images[i, 3], sorted_images[i, 4])
   end
   
   # Print model predicted source position
   if unit == :RA_DEC
      # Convert RA/DEC to arcseconds relative to the reference position
      βx_model, βy_model = AstrometricOps.gnomonic_offsets_arcsec(RA_REF, DEC_REF, βx_model, βy_model)
   end
   println("-"^length(header))
   content = @sprintf("(βx_model, βy_model) = (%3.6f, %3.6f)", βx_model, βy_model)
   @printf("| %-68s |\n", content)
   println("-"^length(header))

   return nothing
end


# --------------------------------------------------------------------------------------------------
# Save best-fit lensing quantities to FITS files
# --------------------------------------------------------------------------------------------------
function _write_fits_header!(header::ImageHDU, model::ModelConfig;
                             date::Union{Nothing, Date} = nothing, time::Union{Nothing, Time} = nothing)
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

"""
    save_best_fits(file_name::String; save_potential::Bool=true, save_deflection::Bool=true, save_kappa::Bool=true, save_gamma::Bool=true)
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
- `nothing`: Saves FITS files to disk
"""
function save_best_fits(data_jld2::JLD2.JLDFile;
                        save_potential::Bool  = true,
                        save_deflection::Bool = true, 
                        save_kappa::Bool      = true, 
                        save_gamma::Bool      = true)

   # Load the chains and logL from the file
   model  = data_jld2["model"]
   chains = data_jld2["chains"]
   logL   = data_jld2["logL"]
   date   = data_jld2["Date"]
   time   = data_jld2["Time"]

   # Get best-fit model
   best_model, _ = get_best_model(model; mcmc_chains=chains, mcmc_logL=logL)

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


# --------------------------------------------------------------------------------------------------
# Check parity of the best-fit model against input parities
# --------------------------------------------------------------------------------------------------
"""
    check_parity(model::ModelConfig, chains::Array{Float64, 3}, logL::Matrix{Float64})
Check the parity of the best-fit model against the input parities for each knot image. This function 
assumes that parity was enforced during the modelling (i.e., `model.source_config.use_parity = true`). 
If parity was not enforced, an error is thrown. The function prints a table comparing the input 
parity and the best-fit model parity for each knot image, along with a status indicating whether the
parity matches or not.

# Arguments
- `data_jld2::JLD2.JLDFile`: JLD2 data file containing the fit results. 

# Returns
- `nothing`: Prints a table to the console with input and best-fit model parities for each knot image.
"""
function check_parity(data_jld2::JLD2.JLDFile)
   # Check if parity was enforced during modelling
   if model.source_config.use_parity === false
      error("Parity was not enforced during modelling. Please set model.use_parity = true.")
   end

   # Load the chains and log-likelihood from the file
   model  = data_jld2["model"]
   chains = data_jld2["chains"]
   logL   = data_jld2["logL"]

   # Get the best parameters based on log-likelihood
   best_θ, _ = get_best_fit_parameters(logL, chains)

   # Get list of parameters for the lens model
   param_ref = Dict(p.key => p.refer for p in model.parameters)
   
   # Replace free parameter values by best-fit values
   pvals = LensModelUtils.param_dict(model, best_θ, param_ref)

   # Get best-fit model
   best_model = LensModelUtils.build_lens(model, pvals)

   # Get angular-diameter distance ratios
   adis = LensModelUtils.adis_current(model, pvals)

   # Calculate deformation at all image positions
   _, _, _, A_all = LensModelUtils.lens_quantities(model, best_model)
   
   # Print Table Header using string padding for alignment
   header = string(
      "| ", rpad("Source", 8), 
      "| ", rpad("Knot", 6), 
      "| ", rpad("Image", 6), 
      "| ", rpad("Input Parity", 10), 
      "| ", rpad("Best Parity", 10), 
      "| ", "Status",
      " |"
   )

   println("─"^length(header))
   println(header)
   println("─"^length(header))
   
   # Identity tuple
   I4 = (1.0, 0.0, 0.0, 1.0)

   # Compare parity for each source
   sid = 1
   kid = 1
   for src in model.source_config.sources
      # Distance ratio for this source
      adis_value = adis[sid]

      # Generate source id
      src_id = Symbol(:src, sid)
   
      for knot in src.knots
         # Generate knot id
         knot_id = Symbol(:knot, kid)

         # Input knot image parities
         parity_input = knot.parity

         # Deformation tensor at the knot positions
         A = @. adis_value * A_all[kid]
         for i in eachindex(A)
            @. A[i] = I4[i] - A[i]
         end

         # Model parity of knot images
         parity_model = @. Int64(sign(A[1] * A[4] - A[2] * A[3]))

         # Parity logL
         for i in eachindex(parity_input)
            # Check if parity is correct
            status = (parity_input[i] == parity_model[i]) ? "✅" : "❌"

            # Manually formatting the sign for the parities
            input_str = parity_input[i] >= 0 ? "+$(parity_input[i])" : "$(parity_input[i])"
            best_str  = parity_model[i] >= 0 ? "+$(parity_model[i])" : "$(parity_model[i])"

            # Format the row using rpad (Right Pad)
            row = string(
               "| ", rpad(src_id, 8), 
               "| ", rpad(knot_id, 6), 
               "| ", rpad(i, 6), 
               "| ", rpad(input_str, 12), 
               "| ", rpad(best_str,  11), 
               "| ", status,
               "     |"
            )
            println(row)
         end
         kid = kid + 1
         println("─"^length(header))
      end
      sid = sid + 1
   end
   return nothing
end


end