import Dates


# --------------------------------------------------------------------------------------------------
# Usage strings
# --------------------------------------------------------------------------------------------------
const _CLI_USAGE = """
Usage: LensFactory <command> [options]

Commands:
  fit <input.yaml>    Run lens model fitting as configured in <input.yaml>

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
# Command registry
# --------------------------------------------------------------------------------------------------
const _CLI_COMMANDS = Dict(
   "fit" => (parse = _cli_parse_fit, run = _cli_run_fit, usage = _CLI_FIT_USAGE),
)


# --------------------------------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------------------------------
"""
    LensFactory.main(args::Vector{String})

Entry point of the `LensFactory` command-line tool, e.g. `LensFactory fit input.yaml`.
Returns the process exit code: 0 on success, 1 if the command failed and 2 for invalid arguments.
"""
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

   command = get(_CLI_COMMANDS, name, nothing)
   if command === nothing
      if endswith(lowercase(name), ".yaml") || endswith(lowercase(name), ".yml")
         _cli_error("missing command. Did you mean `LensFactory fit $name`?")
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