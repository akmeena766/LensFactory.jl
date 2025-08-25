using Documenter
using LensFactory
using Makie

push!(LOAD_PATH, "../src/")
makedocs(
    sitename = "LensFactory.jl",
    modules = [LensFactory, 
                Constants, 
                Cosmology, 
                Lenses, 
                Sources,
                isdefined(Base, :get_extension) ? Base.get_extension(LensFactory, :PlotExt) : LensFactory.PlotExt],
    format = Documenter.HTML(;collapselevel = 1, prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
            "Home" => "index.md",
            "Constants" => "Constants.md",
            "Cosmology" => "Cosmology.md",
            "Lenses" => [
                        "Basics" => "Lenses.md",
                        "Point Lens" => "PointLens.md",
                        "Plummer Lens" => "PlummerLens.md",
                        "SIS Lens" => "SISLens.md"
                        ],
            "Sources" => "Sources.md",
            "Plot Extension" => "PlotExt.md",
        ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/akmeena766/LensFactory.jl.git",
    devbranch = "main"
)
