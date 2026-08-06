# --------------------------------------------------------------------------------------------------
# Test plotting API (PlotExt extension)
# --------------------------------------------------------------------------------------------------
using CairoMakie
using Makie
using KernelDensity
using LaTeXStrings
CairoMakie.activate!()

@testset "PlotExt" begin
   include("./test_Extension/test_PlotLenses.jl")
   # include("./test_Extension/test_PlotMultiPlane.jl")
   # include("./test_Extension/test_PlotLensModel.jl") 
end
