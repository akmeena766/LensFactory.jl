using Documenter
using DocumenterCitations
using LensFactory
using Makie

bib = CitationBibliography(joinpath(@__DIR__, "src", "References.bib"), style = :authoryear)

makedocs(
    sitename = "LensFactory.jl",
    plugins  = [bib],
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
            "Constants"      => "Constants.md",
            "Cosmology"      => "Cosmology.md",
            "Lenses"         => [
                        "General"                      => "Lenses.md",
                        "Point Lens"                   => "PointLens.md",
                        "Plummer Lens"                 => "PlummerLens.md",
                        "SIS Lens"                     => "SISLens.md",
                        "NSISP Lens"                   => "NSISPLens.md",
                        "NSISMD Lens"                  => "NSISMDLens.md",
                        "Gaussian Lens"                => "GaussianLens.md",
                        "Sersic Lens"                  => "SersicLens.md",
                        "External Effects"             => "ExternalEffects.md",
                        "Third Order Perturbation"     => "ExternalEffects3.md",
                        "Multipole Perturbation"       => "Multipole.md",
                        "PIEP Lens"                    => "PIEPLens.md",
                        "SIE Lens"                     => "SIELens.md",
                        "PJE Lens"                     => "PJELens.md",
                        "Hernquist Lens"               => "HernquistLens.md",
                        "elliptical Hernquist MD Lens" => "eHernquistMDLens.md",
                        "approximate Hernquist Lens"   => "aHernquistLens.md",
                        "NFW Lens"                     => "NFWLens.md",
                        "elliptical NFW MD Lens"       => "eNFWMDLens.md",
                        "approximate NFW Lens"         => "aNFWLens.md",
                        "truncated NFW Lens"           => "tNFWLens.md",
                        "generalized NFW Lens"         => "gNFWLens.md",
                        "Einasto Lens"                 => "EinastoLens.md",
                        "Multi-Plummer Lens"           => "MultiPlummerLens.md",
                        "Multi-Gaussian Lens"          => "MultiGaussianLens.md",
                        "Multi-PJE Lens"               => "MultiPJELens.md",
                        "Composite Lens"               => "CompositeLens.md"
                        ],
            "MultiPlane"     => "MultiPlane.md",
            "Sources"        => "Sources.md",
            "LensModel"      => [
                        "General"    => "LensModel.md",
                        "Parametric" => "LensModel_Parametric.md"
                        ],
            "SingularityMap" => "SingularityMap.md",
            "LFUtils"        => "LFUtils.md",
            "Plot Extension" => [
                        "General"       => "PlotExt_Lenses.md",
                        "Multi-Plane"   => "PlotExt_MultiPlane.md",
                        "Lens Modeling" => "PlotExt_LensModel.md"
                        ],
            "Publications"   => "Publications.md",
            "History"        => "History.md",
            "Bibliography"   => "References.md"
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
