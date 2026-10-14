# Parametric method

```yaml
observation:
  modeler: LensFactory
  lens: Abell370
  z_d: 0.375
  reference: [0.0, 0.0]
  pixel_scale: 0.05
  FOV: [200, 200]

cosmology:
  H0: 70.0
  Omega_m0: 0.30
  Omega_r0: 0.00
  Omega_w0: 0.70

source:
  from_file: true
  source_file: A370_images.txt

lens_model:
  method: Parametric
  multiplane: false
  total_lenses: 5
  lens1:

sampling:
  scheme: SourcePlane
  verbose: true
  optimizer:
    enabled: true
    Parametric:
      run_mode: random
      n_runs: 10000
      tolerance: 1E-2
    NM:
      iterations: 100000
      tolerance: 1e-6
  mcmc:
    enabled: true
    AIES:
      n_chains: 12
      n_steps: 500000
```