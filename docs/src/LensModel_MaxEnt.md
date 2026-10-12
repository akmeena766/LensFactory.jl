# Maximum Entropy method

```yaml
lens_model:
  method: MaxEnt
  multiplane: false
  total_lenses: 1
  alpha: 1.0
  lens1:
    lens: MultiPixelLens

sampling:
  scheme: SourcePlane
  verbose: true
  optimizer:
    enabled: true
    MaxEnt:
      n_runs: 1
      save_runs: true
      save_path: ./
      refine: false
      refine_rounds: 0
      refine_factor: 0
    LBFGS:
      iterations: 10000
      tolerence: 1.0E-6
```