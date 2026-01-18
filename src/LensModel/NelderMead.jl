module NelderMead


export nmsmax

"""
    nmsmax(fun, x[, trace = true, initial_simplex = 0, target_f = Inf, max_its = Inf, max_evals = Inf, tol = 1e-3])

Nelder-Mead simplex method for direct search optimization.

This function attempts to maximize the function `fun`, using the
starting vector `x`. The Nelder-Mead direct search method is used.

Adaption of the original *MATLAB* source by Nick Higham to *Julia*.

# Output arguments:
      `x`    = vector yielding largest function value found,
      `fmax` = function value at `x`,
      `nf`   = number of function evaluations,
      `its`  = number of iterations required.

# Parameters:
The iteration is terminated when either
 - the relative size of the simplex is <= `tol`
   (default `1e-3`),
 - max_evals function evaluations have been performed
   (default `Inf`, i.e., no limit), or
 - a function value equals or exceeds `target_f`
   (default `Inf`, i.e., no test on function values).

The form of the initial simplex is determined by initial_simplex:
 - `initial_simplex = 0`: regular simplex (sides of equal length, the default)
 - `initial_simplex = 1`: right-angled simplex.
Progress of the iteration is not shown if `trace = true` (default `false`).

# References:
 - N. J. Higham. The Matrix Computation Toolbox.
   http://www.ma.man.ac.uk/~higham/mctoolbox.
 - N. J. Higham, Optimization by direct search in matrix computations,
   SIAM J. Matrix Anal. Appl, 14(2): 317-333, 1993.
 - C. T. Kelley, Iterative Methods for Optimization, Society for Industrial
   and Applied Mathematics, Philadelphia, PA, 1999.
"""
function nmsmax(fun, x; trace=false, initial_simplex=0, target_f=Inf, max_its=Inf, max_evals=Inf, tol=1e-3)
   x0 = x[:];  # Work with column vector internally.
   n = length(x0);

   # Zero column + Identity matrix
   V = zeros(n, n + 1)
   for i in 1:n
      V[i, i + 1] = 1.0
   end

   f = zeros(n+1,1);
   V[:,1] = x0; f[1] = fun(x);
   fmax_old = f[1];
   fmax     = -Inf; # Some initial value

   if trace
      println("f(x0) = ", round(f[1], sigdigits=5))
   end

   k = 0; m = 0;

   # Set up initial simplex.
   scale = max(norm(x0,Inf),1);
   if initial_simplex == 0
      # Regular simplex - all edges have same length.
      # Generated from construction given in reference [18, pp. 80-81] of [1].
      alpha = scale / (n*sqrt(2)) * [ sqrt(n+1)-1+n  sqrt(n+1)-1 ];
      V[:,2:n+1] = (x0 + alpha[2]*ones(n,1)) * ones(1,n);
      for j=2:n+1
         V[j-1,j] = x0[j-1] + alpha[1];
         x[:] = V[:,j]; f[j] = fun(x);
      end
   else
      # Right-angled simplex based on co-ordinate axes.
      alpha = scale*ones(n+1,1);
      for j=2:n+1
         V[:,j] = x0 + alpha[j]*V[:,j];
         x[:] = V[:,j]; f[j] = fun(x);
      end
   end
   nf = n+1;
   how = "initial  ";

   j = sortperm(f[:]);
   temp = f[j];
   j = j[n+1:-1:1];
   f = f[j]; V = V[:,j];

   alpha = 1;  beta = 1/2;  gamma = 2;

   msg = ""

   converged = false
   while true    ###### Outer (and only) loop.
   k = k+1;

      fmax = f[1];
      if fmax > fmax_old
         if trace
            # Calculation for the percentage change
            pct_change = 100 * (fmax - fmax_old) / (abs(fmax_old) + eps(fmax_old))

            # The combined println statement
            println("Iter. ", Int64(k), ",  how = ", how, " nf = ", Int64(nf), ",  f = ", round(fmax, sigdigits=5), "  (", round(pct_change, digits=1), "%)")
         end
      end
      fmax_old = fmax;

      ### Three stopping tests from MDSMAX.M

      # Stopping Test 1 - f reached target value?
      if fmax >= target_f
         msg = "Exceeded target...quitting\n";
         break  # Quit.
      end

      # Stopping Test 2 - too many f-evals?
      if nf >= max_evals
         msg = "Max no. of function evaluations exceeded...quitting\n";
         break  # Quit.
      end

      # Stopping Test 3 - too many iterations?
      if k > max_its
         msg = "Max no. of iterations exceeded...quitting\n";
         break  # Quit.
      end

      # Stopping Test 4 - converged?   This is test (4.3) in [1].
      v1 = V[:,1];
      size_simplex = norm(V[:,2:n+1]-v1[:,ones(Int,n)],1) / max(1, norm(v1,1));
      if size_simplex <= tol
         msg = "Simplex size $(round(size_simplex, sigdigits=5)) <= $(round(tol, sigdigits=5))...quitting\n"
         converged = true
         break  # Quit.
      end

      #  One step of the Nelder-Mead simplex algorithm
      #  NJH: Altered function calls and changed CNT to NF.
      #       Changed each `fr < f[1]' type test to `>' for maximization
      #       and re-ordered function values after sort.

      vbar = (sum(V[:,1:n]',1)/n)';  # Mean value
      vr = (1 + alpha)*vbar - alpha*V[:,n+1]; x[:] = vr; fr = fun(x);
      nf = nf + 1;
      vk = vr;  fk = fr; how = "reflect, ";
      if fr > f[n]
               if fr > f[1]
                  ve = gamma*vr + (1-gamma)*vbar; x[:] = ve; fe = fun(x);
                  nf = nf + 1;
                  if fe > f[1]
                     vk = ve; fk = fe;
                     how = "expand,  ";
                  end
               end
      else
               vt = V[:,n+1]; ft = f[n+1];
               if fr > ft
                  vt = vr;  ft = fr;
               end
               vc = beta*vt + (1-beta)*vbar; x[:] = vc; fc = fun(x);
               nf = nf + 1;
               if fc > f[n]
                  vk = vc; fk = fc;
                  how = "contract,";
               else
                  for j = 2:n
                     V[:,j] = (V[:,1] + V[:,j])/2;
                     x[:] = V[:,j]; f[j] = fun(x);
                  end
                  nf = nf + n-1;
                  vk = (V[:,1] + V[:,n+1])/2; x[:] = vk; fk = fun(x);
                  nf = nf + 1;
                  how = "shrink,  ";
               end
      end
      V[:,n+1] = vk;
      f[n+1] = fk;
      j = sortperm(f[:]);
      temp = f[j];
      j = j[n+1:-1:1];
      f = f[j]; V = V[:,j];

   end   ###### End of outer (and only) loop.

   # Finished.
   if trace
      println(msg)
   end
   x[:] = V[:,1];

   return x, fmax, nf, k-1, converged
end

end