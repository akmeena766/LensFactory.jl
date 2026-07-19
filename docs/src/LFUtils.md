# LFUtils

The `LFUtils` module collects various general-purpose utility functions used throughout
`LensFactory`. It is organized into four submodules: `AstrometricOps` for astrometric
coordinate transformations, `ContourFinder` for tracing iso-contours of 2D scalar fields,
`IntersectionFinder` for computing intersections of curves, and `PolygonOps` for polygon and
interpolation operations.

### AstrometricOps

```@docs
LFUtils.AstrometricOps.gnomonic_offsets_arcsec
LFUtils.AstrometricOps.gnomonic_offsets_radec
```

### ContourFinder

```@docs
LFUtils.ContourFinder.get_contour
```

### IntersectionFinder

```@docs
LFUtils.IntersectionFinder.get_intersection
```

### PolygonOps

```@docs
LFUtils.PolygonOps.bilinear_interpolation
LFUtils.PolygonOps.shoelace
LFUtils.PolygonOps.hao_sun
LFUtils.PolygonOps.fit_ellipse
```
