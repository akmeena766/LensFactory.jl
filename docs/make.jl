using Documenter
using LensFactory
using Makie

makedocs(
    sitename = "LensFactory.jl",
    modules = [LensFactory, 
                Constants, 
                Cosmology, 
                Lenses, 
                Sources,
                Base.get_extension(LensFactory, :PlotExt)],
    format = Documenter.HTML(;collapselevel = 1, 
                            assets = ["assets/custom.css"],
                            prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
            "Home" => "index.md",
            "Constants" => "Constants.md",
            "Cosmology" => "Cosmology.md",
            "Lenses" => [
                        "Basics"           => "Lenses.md",
                        "Point Lens"       => "PointLens.md",
                        "Plummer Lens"     => "PlummerLens.md",
                        "SIS Lens"         => "SISLens.md",
                        "NSISP Lens"       => "NSISPLens.md",
                        "NSISMD Lens"      => "NSISMDLens.md",
                        "Gaussian Lens"    => "GaussianLens.md",
                        "Sersic Lens"      => "SersicLens.md",
                        "External Effects" => "ExternalEffects.md",
                        "PIEP Lens"        => "PIEPLens.md",
                        "SIE Lens"         => "SIELens.md",
                        "PJE Lens"         => "PJELens.md",
                        "Hernquist Lens"   => "HernquistLens.md",
                        "NFW Lens"         => "NFWLens.md",
                        "truncated NFW Lens"   => "tNFWLens.md",
                        "generalized NFW Lens" => "gNFWLens.md",
                        "Einasto Lens"     => "EinastoLens.md",
                        "approximate Hernquist Lens" => "aHernquistLens.md",
                        "approximate NFW Lens" => "aNFWLens.md",
                        "elliptical Hernquist MD Lens" => "eHernquistMDLens.md",
                        "elliptical NFW MD Lens" => "eNFWMDLens.md"
                        ],
            "Multi-plane lensing" => "MultiPlane.md",
            "Sources" => "Sources.md",
            "Plot Extension" => "PlotExt.md",
            "History" => "History.md",
        ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(;
    repo = "github.com/akmeena766/LensFactory.jl.git",
    devbranch = "main",
    branch = "gh-pages"
)
