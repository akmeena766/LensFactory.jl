using Documenter
using LensFactory

push!(LOAD_PATH,"../src/")
makedocs(
    sitename = "LensFactory.jl",
    modules = [LensFactory, Constants, Cosmology, Lenses],
    format = Documenter.HTML(;prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
            "Home" => "index.md",
            "Constants" => "Constants.md",
            "Cosmology" => "Cosmology.md",
            "Lenses" => [
                        "Basics" => "Lenses.md",
                        "Point Lens" => "PointLens.md",
                        "SIS Lens" => "SISLens.md",
                        ],
        ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/akmeena766/LensFactory.jl.git",
    devbranch = "main"
)
