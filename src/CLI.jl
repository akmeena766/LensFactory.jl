import Dates


# --------------------------------------------------------------------------------------------------
# Usage strings
# --------------------------------------------------------------------------------------------------
const _CLI_USAGE = """
Usage: LensFactory <command> [options]

Commands:
  fit <input.yaml>    Run lens model fitting as configured in <input.yaml>
  bestfit <results.jld2>  Write FITS maps of the best-fit model from a fit's results

Global options:
  -h, --help          Show this message
  -V, --version       Print the LensFactory version and exit

Run `LensFactory <command> --help` for command-specific options.

Threads are set with the JULIA_NUM_THREADS environment variable, e.g.
  JULIA_NUM_THREADS=8 LensFactory fit input.yaml
"""

const _CLI_FIT_USAGE = """
Usage: LensFactory fit <input.yaml> [options]

Run strong-lens mass modeling as configured in <input.yaml>.

Options:
  -o, --output FILE   Name of the output JLD2 file
                      (default: <LensName>_<yyyymmdd>.jld2 in the current directory)
      --no-save       Run the fit without writing results to disk
  -h, --help          Show this message and exit

The optimizer and MCMC samplers are multithreaded. Set the number of threads with
the JULIA_NUM_THREADS environment variable ("auto" uses all CPU threads).

Examples:
  LensFactory fit input.yaml
  JULIA_NUM_THREADS=8 LensFactory fit input.yaml -o run1.jld2
  export JULIA_NUM_THREADS=auto      # e.g. once in ~/.zshrc
"""

const _CLI_BESTFIT_USAGE = """
Usage: LensFactory bestfit <results.jld2> [options]
Write FITS maps (potential, deflection, convergence, shear) of the best-fit lens model
from the JLD2 file written by `LensFactory fit`. The best fit is taken from the MCMC
chains if present, otherwise from the optimizer results.

Options:
  -d, --outdir DIR    Directory to write the FITS files to (default: current directory)
      --burn-in F     Fraction of MCMC steps to discard as burn-in, 0 <= F < 1 (default: 0.2)
      --only LIST     Comma-separated subset of quantities to write:
                        potential   -> potential.fits
                        deflection  -> alpha_x.fits, alpha_y.fits
                        kappa       -> kappa.fits
                        gamma       -> gamma1.fits, gamma2.fits
                      (default: all)
  -f, --force         Overwrite existing FITS files
  -h, --help          Show this message and exit

Examples:
  LensFactory bestfit MACS0416_20261020.jld2
  LensFactory bestfit run1.jld2 -d run1_maps --only kappa,gamma
  LensFactory bestfit run1.jld2 --burn-in 0.5 --force
"""


# --------------------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------------------
function _cli_error(msg::AbstractString)
   printstyled(stderr, "ERROR: "; bold = true, color = :red)
   println(stderr, msg)
end

# Print an argument error with a pointer to the relevant help and return exit code 2
function _cli_usage_error(err::Exception, help_cmd::AbstractString)
   err isa ArgumentError || rethrow(err)
   _cli_error(err.msg)
   println(stderr, "Run `$help_cmd --help` for usage.")
   return 2
end

# Format a duration given in seconds: "42.3 s", "12 min 5 s" or "2 h 13 min"
function _cli_format_duration(seconds::Real)
   if round(seconds; digits = 1) < 60
      return "$(round(seconds; digits = 1)) s"
   end
 
   total = round(Int, seconds)
   if total < 3600
      m, s = divrem(total, 60)
      return "$m min $s s"
   else
      h, rest = divrem(total, 3600)
      return "$h h $(div(rest, 60)) min"
   end
end


# --------------------------------------------------------------------------------------------------
# `fit` subcommand
# --------------------------------------------------------------------------------------------------
function _cli_parse_fit(args::Vector{String})
   input  = nothing
   output = nothing
   save   = true

   i = 1
   while i <= length(args)
      a = args[i]
      if a == "-h" || a == "--help"
         return (; help = true)
      elseif a == "-o" || a == "--output"
         i == length(args) && throw(ArgumentError("option '$a' requires a file name"))
         i += 1
         output = args[i]
      elseif startswith(a, "--output=")
         output = String(chopprefix(a, "--output="))
         isempty(output) && throw(ArgumentError("option '--output' requires a file name"))
      elseif a == "--no-save"
         save = false
      elseif a == "-t" || startswith(a, "--threads")
         throw(ArgumentError("threads are set with the JULIA_NUM_THREADS environment variable, " *
                             "e.g. `JULIA_NUM_THREADS=8 LensFactory fit input.yaml`"))
      elseif startswith(a, "-")
         throw(ArgumentError("unknown option '$a'"))
      else
         input === nothing ||
            throw(ArgumentError("only one input file can be given (got '$input' and '$a')"))
         input = a
      end
      i += 1
   end

   input === nothing && throw(ArgumentError("no input file given"))
   isfile(input)     || throw(ArgumentError("input file not found: $input"))
   return (; help = false, input = abspath(input), output = output, save = save)
end

function _cli_run_fit(opts)
   @info "LensFactory v$(pkgversion(@__MODULE__))" input = opts.input threads = Threads.nthreads()
   if Threads.nthreads() == 1
      @info "Running on a single thread. Set JULIA_NUM_THREADS (e.g. `JULIA_NUM_THREADS=8`) to use more cores."
   end

   model = LensModel.read_input(opts.input)

   # Resolve the output name here (same default as `fit_model`) so it can be reported
   file_name = opts.output
   if opts.save && file_name === nothing
      file_name = "$(model.observation.lens)_$(Dates.format(Dates.today(), "yyyymmdd")).jld2"
   end

   @info "Starting fit" lens = model.observation.lens
   elapsed = @elapsed LensModel.fit_model(model; save = opts.save, file_name = file_name)
 
   @info "Fit finished in $(_cli_format_duration(elapsed))"
   opts.save && @info "Results saved to $(abspath(file_name))"
   return 0
end


# --------------------------------------------------------------------------------------------------
# `bestfit` subcommand
# --------------------------------------------------------------------------------------------------
const _CLI_BESTFIT_QUANTITIES = (:potential, :deflection, :kappa, :gamma)

# Value of an option given either as `--opt value` or `--opt=value`. Returns (value, new index).
function _cli_option_value(args::Vector{String}, i::Int, long::String)
   a = args[i]
   if startswith(a, long * "=")
      value = String(chopprefix(a, long * "="))
   else
      i == length(args) && throw(ArgumentError("option '$a' requires a value"))
      i += 1
      value = args[i]
   end
   isempty(value) && throw(ArgumentError("option '$long' requires a value"))
   return value, i
end

function _cli_parse_bestfit(args::Vector{String})
   input      = nothing
   outdir     = "."
   burn_in    = 0.2
   quantities = collect(_CLI_BESTFIT_QUANTITIES)
   force      = false

   is_opt(a, long, short = nothing) = a == long || a == short || startswith(a, long * "=")

   i = 1
   while i <= length(args)
      a = args[i]
      if a == "-h" || a == "--help"
         return (; help = true)
      elseif is_opt(a, "--outdir", "-d")
         outdir, i = _cli_option_value(args, i, "--outdir")
      elseif is_opt(a, "--burn-in")
         value, i = _cli_option_value(args, i, "--burn-in")
         burn_in  = tryparse(Float64, value)
         (burn_in !== nothing && 0.0 <= burn_in < 1.0) ||
            throw(ArgumentError("--burn-in must be a number in [0, 1), got '$value'"))
      elseif is_opt(a, "--only")
         value, i   = _cli_option_value(args, i, "--only")
         quantities = Symbol[]
         for q in split(value, ',')
            q = Symbol(strip(lowercase(q)))
            q in _CLI_BESTFIT_QUANTITIES ||
               throw(ArgumentError("unknown quantity '$q' in --only. " *
                                   "Choose from: $(join(_CLI_BESTFIT_QUANTITIES, ", "))"))
            q in quantities || push!(quantities, q)
         end
      elseif a == "-f" || a == "--force"
         force = true
      elseif startswith(a, "-")
         throw(ArgumentError("unknown option '$a'"))
      else
         input === nothing ||
            throw(ArgumentError("only one results file can be given (got '$input' and '$a')"))
         input = a
      end
      i += 1
   end

   input === nothing && throw(ArgumentError("no results file given"))
   isfile(input)     || throw(ArgumentError("results file not found: $input"))
   if endswith(lowercase(input), ".yaml") || endswith(lowercase(input), ".yml")
      throw(ArgumentError("`bestfit` expects the .jld2 results of `LensFactory fit`, not the " *
                          "input YAML. Run `LensFactory fit $input` first."))
   end
   return (; help = false, input = abspath(input), outdir = abspath(outdir),
             burn_in = burn_in, quantities = quantities, force = force)
end

function _cli_run_bestfit(opts)
   @info "LensFactory v$(pkgversion(@__MODULE__))" input = opts.input outdir = opts.outdir

   # Refuse to silently overwrite earlier maps
   targets = [joinpath(opts.outdir, f) for q in opts.quantities for f in LensModel.BEST_FIT_FILES[q]]
   existing = filter(isfile, targets)
   if !opts.force && !isempty(existing)
      _cli_error("output file(s) already exist:\n  " * join(existing, "\n  ") *
                 "\nUse --force to overwrite or --outdir to write elsewhere.")
      return 1
   end

   @info "Writing best-fit maps" quantities = join(opts.quantities, ", ") burn_in = opts.burn_in
   written = String[]
   elapsed = @elapsed begin
      written = LensModel.save_best_fits(opts.input;
                                         out_dir         = opts.outdir,
                                         burn_in         = opts.burn_in,
                                         save_potential  = :potential  in opts.quantities,
                                         save_deflection = :deflection in opts.quantities,
                                         save_kappa      = :kappa      in opts.quantities,
                                         save_gamma      = :gamma      in opts.quantities)
   end

   @info "Wrote $(length(written)) FITS file(s) in $(_cli_format_duration(elapsed))" files = join(written, "\n")
   return 0
end



# --------------------------------------------------------------------------------------------------
# Command registry
# --------------------------------------------------------------------------------------------------
const _CLI_COMMANDS = Dict(
   "fit"     => (parse = _cli_parse_fit,     run = _cli_run_fit,     usage = _CLI_FIT_USAGE),
   "bestfit" => (parse = _cli_parse_bestfit, run = _cli_run_bestfit, usage = _CLI_BESTFIT_USAGE),
)

# Alternative spellings accepted for commands
const _CLI_ALIASES = Dict("best_fit" => "bestfit", "best-fit" => "bestfit")

# --------------------------------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------------------------------
# """
#     LensFactory.main(args::Vector{String})

# Entry point of the `LensFactory` command-line tool, e.g. `LensFactory fit input.yaml`.
# Returns the process exit code: 0 on success, 1 if the command failed and 2 for invalid arguments.
# """
function (@main)(args::Vector{String})
   rest = args
   if isempty(rest)
      print(stderr, _CLI_USAGE)
      return 2
   end

   name = rest[1]
   if name == "-h" || name == "--help"
      print(stdout, _CLI_USAGE)
      return 0
   elseif name == "-V" || name == "--version"
      println(stdout, "LensFactory v", pkgversion(@__MODULE__))
      return 0
   end

   name    = get(_CLI_ALIASES, name, name)
   command = get(_CLI_COMMANDS, name, nothing)
   if command === nothing
      if endswith(lowercase(name), ".yaml") || endswith(lowercase(name), ".yml")
         _cli_error("missing command. Did you mean `LensFactory fit $name`?")
      elseif endswith(lowercase(name), ".jld2")
         _cli_error("missing command. Did you mean `LensFactory bestfit $name`?")
      else
         _cli_error("unknown command '$name'. Available: $(join(sort!(collect(keys(_CLI_COMMANDS))), ", "))")
      end
      return 2
   end

   # Parse the command's own options
   opts = try
      command.parse(rest[2:end])
   catch err
      return _cli_usage_error(err, "LensFactory $name")
   end

   if opts.help
      print(stdout, command.usage)
      return 0
   end

   try
      return command.run(opts)
   catch err
      err isa InterruptException && (_cli_error("interrupted"); return 130)
      _cli_error("`$name` failed")
      showerror(stderr, err, catch_backtrace())
      println(stderr)
      return 1
   end
end