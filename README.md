# Rayleigh-Gans Scattering Viewer

![Rayleigh-Gans scattering viewer](docs/screenshot.png)

An interactive sandbox for exploring how Rayleigh-Gans form factors change with
geometry, size, wavelength, and scattering vector.

It is meant for quick physical intuition:

- A PyQt5 desktop viewer for quick visual comparisons
- Reusable form-factor functions in `src/rayleigh_gans/form_factors.py`
- Notebooks for derivations, examples, and reproducible figures

## What You Can Explore

- How minima, oscillations, and envelope shape shift as particle size changes
- How different geometries compare at the same `q` or scattering angle
- When shell thickness starts to matter
- How particle models differ from Debye and Guinier polymer descriptions

## Quick start

```bash
pip install -e .
rayleigh-gans-viewer
```

You can also launch the app with `python run_app.py` or `python -m rayleigh_gans`.

## Features

- Compare multiple geometries on the same plot
- Switch between scattering angle `theta` and scattering vector `q`
- Inspect `Re[F(q)]`, `|F(q)|`, and intensity `|F(q)|^2`
- Mix particle and polymer models in the same viewer
- Export figures directly to PNG

## Implemented models

- Sphere
- Thin spherical shell
- Thick spherical shell
- Disk (2D quick-look model)
- Gaussian polymer coil (Debye)
- Guinier approximation

## Installation

For the app:

```bash
pip install -e .
```

For notebooks as well:

```bash
pip install -e .[dev]
```

## Running the viewer

```bash
python -m rayleigh_gans
```

```bash
rayleigh-gans-viewer
```

```bash
python run_app.py
```

## Running the tests

```bash
python -m unittest discover -s tests -v
```

## Notebook usage

The notebooks in `Notebooks/` use `Notebooks/_bootstrap.py` to locate `src/`
without requiring installation first. That makes it easy to experiment locally
while still keeping the package layout clean.

## Physics notes

- All angles are in radians unless otherwise stated.
- The scattering vector is defined as `q = 2k sin(theta/2)` with `k = 2*pi*n / wavelength`.
- `|F(q)|^2` corresponds to the scattering intensity.
- Debye and Guinier are intensity-only models in this repository.
- Length units should be kept consistent across particle size, `Rg`, and `wavelength`.

## Repository layout

```text
.
|- docs/
|- Notebooks/
|- src/rayleigh_gans/
|- tests/
|- pyproject.toml
`- run_app.py
```

## Possible next extensions

- Orientation-averaged disk or rod models
- Additional polymer / aggregate form factors
- Preset parameter bundles for common scattering scenarios
- Figure export presets for publication-quality output
