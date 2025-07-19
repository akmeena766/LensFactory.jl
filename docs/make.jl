using Documenter
using LensFactory

push!(LOAD_PATH,"../src/")
makedocs(
    sitename = "LensFactory.jl",
    pages = [
        "Index" => "index.md"
    ],
    format = Documenter.HTML(prettyurls = false),
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/akmeena766/LensFactory.jl.git",
    devbranch = "main"
)
