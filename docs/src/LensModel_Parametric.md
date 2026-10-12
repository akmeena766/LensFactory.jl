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
```