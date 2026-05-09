module NM

# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using Printf


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export nmsmax


# """
#     nmsmax(fun, x;
#                   tol::Float64            = 1e-3,
#                   target_f::Float64       = Inf,
#                   max_its::Real           = Inf,
#                   max_evals::Real         = Inf,
#                   initial_simplex::Symbol = :regular,
#                   trace::Bool             = false)  

# Nelder-Mead simplex method for direct search **maximization**. This function attempts to maximize
# the function `fun`, using the starting vector `x`. Adapted from the original **MATLAB** source by 
# Nick Higham.

# # Arguments
# - `fun`: Objective function to maximize, called as `fun(v::Vector)`
# - `x`:   Initial vector guess.

# # Keyword arguments
# - `tol::Float64 = 1e-3`:            Convergence tolerance on relative simplex size
# - `target_f::Float64 = Inf`:        Stop early if `fmax` reaches or exceeds this value
# - `max_its::Real = Inf`:            Maximum number of iterations
# - `max_evals::Real = Inf`:          Maximum number of function evaluations
# - `initial_simplex::Symbol = :regular`: `:regular` (equal-length edges) or `:right_angle`
# - `trace::Bool = false`:            Print iteration progress when `true`

# # Returns
# - `x`:         Vector yielding the largest function value found
# - `fmax`:      Corresponding function value
# - `n_evals`:   Total number of function evaluations
# - `n_iters`:   Number of iterations performed
# - `converged`: `true` if the simplex-size tolerance was met

# # References
#  - N. J. Higham. The Matrix Computation Toolbox (http://www.ma.man.ac.uk/~higham/mctoolbox).
#  - N. J. Higham, Optimization by direct search in matrix computations, SIAM J. Matrix Anal. Appl, 14(2): 317-333, 1993.
#  - C. T. Kelley, Iterative Methods for Optimization, SIAM, 1999.
# """
function nmsmax(fun, x; 
                       tol::Float64            = 1e-3,
                       target_f::Float64       = Inf,
                       max_its::Real           = Inf,
                       max_evals::Real         = Inf,
                       initial_simplex::Symbol = :regular,
                       trace::Bool             = false)
   # flatten to 1-D; x itself is never mutated or written back
   x0 = x[:]
   n  = length(x0)

   # -- Allocate simplex (columns = vertices) and their function values ----------------------------
   simplex = zeros(n, n + 1)
   fvals   = zeros(n + 1)

   simplex[:, 1] = x0
   fvals[1]      = fun(x0)
   n_evals       = 1

   trace && @printf("f(x0) = %.5g\n", fvals[1])

   # -- Build initial simplex ----------------------------------------------------------------------
   scale = max(maximum(abs.(x0)), 1.0)

   if initial_simplex == :regular
      # All edges have same length.
      # Generated from construction given in reference [18, pp. 80-81] of [1].
      c     = scale / (n * sqrt(2))
      long  = c * (sqrt(n + 1) - 1 + n)  # Longer coordiate shift
      short = c * (sqrt(n + 1) - 1)      # Shorter coordiate shift

      for j in 2:n+1
         v        = x0 .+ short         # Base: shift every coordinate by `short`
         v[j - 1] = x0[j - 1] + long    # Override one axis with `long`
         
         simplex[:, j] = v
         fvals[j]      = fun(v)
      end
   elseif initial_simplex == :right_angle
      # Vertices along coordinate axes at distance `scale` from x0.
      for j in 2:n+1
         v         = copy(x0)
         v[j - 1] += scale
         
         simplex[:, j] = v
         fvals[j]      = fun(v)
      end
   else
      throw(ArgumentError("Initial_simplex must be :regular or :right_angle"))
   end
   n_evals = n_evals + n

   # -- sort simplex columns by descending function value ------------------------------------------
   order   = sortperm(fvals, rev=true)
   fvals   = fvals[order]
   simplex = simplex[:, order]

   # -- Nelder-Mead coefficients -------------------------------------------------------------------
   α = 1.0   # Reflection
   β = 0.5   # Contraction
   γ = 2.0   # Expansion
   
   # -- Main loop ----------------------------------------------------------------------------------
   converged = false
   fmax_prev = fvals[1]
   step_name = "initial "
   msg       = ""
   n_iters   = 0

   while true
      n_iters = n_iters + 1

      fmax = fvals[1]
      
      if trace && fmax > fmax_prev
         pct = 100 * (fmax - fmax_prev) / (abs(fmax_prev) + eps(fmax_prev))
         @printf("Iter %4d | %-9s | evals = %5d | f = %.5g  (%+.1f%%)\n", n_iters, step_name, n_evals, fmax, pct)
      end
      fmax_prev = fmax

      # -- Stopping criteria -----------------------------------------------------------------------

      # Stopping Test 1 - f reached target value?
      if fmax >= target_f
         msg = "Exceeded target...quitting\n"
         break
      end

      # Stopping Test 2 - too many f-evals?
      if n_evals >= max_evals
         msg = "Maximum function evaluations exceeded...quitting\n"
         break
      end

      # Stopping Test 3 - too many iterations?
      if n_iters > max_its
         msg = "Max no. of iterations exceeded...quitting\n"
         break
      end

      # Stopping Test 4 - converged?   This is test (4.3) in [1].
      best         = simplex[:,1]
      size_simplex = sum(abs.(simplex[:, 2:end] .- best)) / max(1.0, sum(abs.(best)))
      if size_simplex <= tol
         msg = @sprintf("Converged: simplex size %.5g ≤ %.5g", size_simplex, tol)
         converged = true
         break
      end

      # -- One Nelder-Mead step --------------------------------------------------------------------
      worst    = simplex[:, end]
      centroid = vec(sum(simplex[:, 1:n], dims=2) ./ n)
      
      # Reflection
      v_reflect = (1 + α) .* centroid - α .* worst 
      f_reflect = fun(v_reflect)
      n_evals   = n_evals + 1

      step_name = "reflect, "
      v_new     = v_reflect
      f_new     = f_reflect

      if f_reflect > fvals[n]
         if f_reflect > fvals[1]
            # Expansion: try to go further in the reflection direction
            v_expand = γ .* v_reflect .+ (1 - γ) .* centroid
            f_expand = fun(v_expand)
            n_evals  = n_evals + 1
            if f_expand > fvals[1]
               v_new     = v_expand
               f_new     = f_expand
               step_name = "expand,  "
            end
         end
      else
         # Contraction: stay closer to the centroid
         v_contract_base = f_reflect > fvals[end] ? v_reflect : worst

         v_contract = β .* v_contract_base .+ (1 - β) .* centroid
         f_contract = fun(v_contract)
         n_evals    = n_evals + 1

         if f_contract > fvals[n]
            v_new     = v_contract  
            f_new     = f_contract
            step_name = "contract, "
         else
            for j in 2:n
               # Shrink: pull all vertices (except best) halfway toward best
               simplex[:, j] = (best .+ simplex[:, j]) ./ 2
               fvals[j]      = fun(simplex[:, j])
            end
            n_evals   = n_evals + n - 1
            v_new     = (best + worst) ./ 2
            f_new     = fun(v_new)
            n_evals   = n_evals + 1
            step_name = "shrink,  "
         end
      end

      # Replace worst vertex with the new candidate, then re-sort
      simplex[:,end] = v_new
      fvals[end]     = f_new
      order   = sortperm(fvals, rev=true)
      fvals   = fvals[order]
      simplex = simplex[:, order]
   end

   # Finished.
   trace && println(msg)

   return simplex[:, 1], fvals[1], n_evals, n_iters, converged
end

end