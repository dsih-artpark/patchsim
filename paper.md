---
title: "PatchSim: a modular metapopulation simulation framework for compartmental epidemiology"
authors:
  - name: "Adish Assain"
    affiliation: "ArtPark Research"
  - name: "Sneha S"
    affiliation: "ArtPark Research"
  - name: "Shreya Mukherjee"
    affiliation: "ArtPark Research"
date: 2026-07-10
bibliography: paper.bib
---

PatchSim is a lightweight, configuration-first simulation framework for
patch-based (metapopulation) compartmental epidemiological models.  It
supports YAML-defined models (SIR/SEIR variants), network-coupled transmission
between patches, and both ODE and discrete-time solvers.

Statement of need
-----------------

Many disease-modelling workflows require concise, reproducible ways to run
compartmental models across multiple spatial units (patches) with flexible
mixing and parameter overrides. PatchSim provides:

- A YAML-first model and scenario specification so non-developers can edit
  experiments without changing code.
- Network coupling between patches (weighted mixing matrices) and per-patch
  parameter overrides.
- A minimal CLI and Python API for integration into pipelines and analysis.

Users who need a small, auditable engine for running scenario ensembles or
teaching multi-patch epidemiology will find PatchSim convenient compared to
heavier frameworks.

Key features
------------

- YAML templates and CLI scaffolding: `patchsim init` creates a runnable
  project with sample data.
- Multi-patch ODE and discrete solvers with network-mediated force-of-infection.
- Per-patch parameter overrides and seed files for precise initial conditions.
- Automated plotting and CSV outputs for downstream analysis.

Installation
------------

Install from PyPI:

```bash
pip install patchsim
```

Install from source (developer):

```bash
git clone https://github.com/dsih-artpark/patchsim
cd patchsim
uv venv
source .venv/bin/activate
uv pip install -e .[dev]
```

Quick worked example
--------------------

Run the bundled sample SIR ODE scenario (creates CSV + PNG outputs):

```bash
uv run python examples/run_sample_simulation.py
```

Or use the CLI to scaffold and run a template project:

```bash
uv run patchsim init my-project --template sir
cd my-project
uv run patchsim run -c config.yaml
```

API and reproducibility
-----------------------

The primary programmatic entry points are in the `patchsim` package:

- `patchsim.load_config(path)` — load and resolve a YAML config
- `patchsim.setup_simulation(config)` — prepare network model, initial state
- `patchsim.run_simulation(config, model_name, net, y0, patches, n)` — run and
  persist results (CSV and PNG)

Quality control
---------------

- Unit and integration tests with `pytest` are provided in `tests/`.
- Continuous integration runs linting, tests, and now uploads coverage
  reports (see `.github/workflows/ci.yml`).

To run the tests locally:

```bash
uv run pytest -q
```

JOSS-specific notes
-------------------

Before submission to JOSS, create a GitHub release (for example `v0.1.0`) and
archive that release on Zenodo to obtain a DOI. Replace the DOI placeholder
in `README.md` and the `doi` field in `paper.bib` with the minted Zenodo DOI.

References
----------

Please cite this software as:

@software{patchsim,
  title = {PatchSim},
  author = {Assain, Adish and S, Sneha and Mukherjee, Shreya},
  year = {2026},
  version = {v0.1.0},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.REPLACE_WITH_RELEASE_DOI},
}

Acknowledgements
----------------

This work was supported by ArtPark Research. Contributors and additional
acknowledgements are recorded in the Git history.
