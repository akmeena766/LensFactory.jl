# Composite lens
Often, we need to combine different type of lens models to study their combined effect. For example,
to model a galaxy scale lens, we need combine an `SIE` lens with external effects. Similarly, for
cluster lenses, we need to model dark matter halo(s), cluster galaxies, and external effects with 
various profiles. The `CompositeLens` exactly does that, it combines different lens models and 
calculates the corresponding properties. For example on `CompositeLens`, checkout
[Basic-Example-2](https://github.com/akmeena766/LensFactory_Examples/blob/main/Basic/Example2_Point_mass_lens.ipynb)
and
[Basic-Example-3](https://github.com/akmeena766/LensFactory_Examples/blob/main/Basic/Example3_Plummer_lens_model.ipynb)

```@docs
Lenses.init_CompositeLens
```