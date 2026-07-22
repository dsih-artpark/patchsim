---
title: 'PatchSim: a modular simulation framework for patch-based metapopulation epidemiology'
tags:
  - Python
  - epidemiology
  - metapopulation
  - compartmental models
  - disease modelling
  - spatial networks
authors:
  - name: Adish Assain Illikkal
    orcid: 0009-0001-1009-4560
    affiliation: 1
  - name: Sneha S
    orcid: 0009-0009-1854-6443
    affiliation: 1
  - name: Shreya Mukherjee
    orcid: 0009-0008-3629-0376
    affiliation: 1
affiliations:
  - index: 1
    name: AI and Robotics Technology Park (ARTPARK), Indian Institute of Science, Bengaluru, India
date: 16 June 2026
bibliography: paper.bib
---

# Summary

PatchSim is an open-source Python framework for simulating infectious disease dynamics
across spatially structured populations. It implements a patch-based metapopulation
approach in which a landscape is partitioned into discrete geographical units
(patches) — such as subdistricts, districts, or ecological zones — connected by a
network representing population movement or contact. Within each patch, disease
progression is described by user-defined compartmental models (e.g., SIR, SEIR, SIRS)
specified through a YAML configuration file using a declarative arrow-map transition
syntax. PatchSim supports both ODE and discrete-timestep solvers, sensitivity analysis
via the SALib library [@Herman2017], and a dual-mode interface comprising a
command-line tool (CLI) and a Python SDK. The framework is disease-agnostic and
configuration-first, enabling rapid iteration between model specification, simulation,
and output analysis without requiring changes to source code. PatchSim is available
on PyPI (`pip install patchsim`), with source code on GitHub [@patchsim] and
documentation on Read the Docs.

# Statement of need

Mathematical modelling of infectious diseases is essential for understanding
transmission dynamics, evaluating interventions, and informing public health policy
[@Keeling2008]. Metapopulation models — which capture the heterogeneity of disease
spread across connected spatial units — are particularly valuable for diseases where
population movement shapes outbreak trajectories, including foot-and-mouth disease in
livestock, dengue in urban landscapes, and respiratory infections across administrative
regions [@Grenfell1997; @Balcan2009].

Despite their importance, practitioners face a persistent barrier: moving from model
specification to reproducible simulation typically requires either significant custom
scientific computing code or domain-specific platforms tightly coupled to particular
diseases or geographies. This demands substantial software engineering effort that is
a barrier for epidemiologists and public health researchers who are not primarily
software developers.

PatchSim addresses this gap through a configuration-first design in which model
structure — compartments, transitions, and parameters — is declared in a YAML file
rather than implemented in code. The intended users are epidemiological modellers,
researchers studying spatial disease dynamics, and public health analysts who need to
rapidly prototype and compare scenarios across spatial configurations and disease
systems. PatchSim was developed at ARTPARK, IISc, to support active modelling work on
livestock disease dynamics in India, including foot-and-mouth disease across cattle and
buffalo populations stratified by age group and species.

# State of the field

Several software packages exist for compartmental disease modelling. EpiFire
[@Hladish2012] provides network-based simulation but requires C++ and does not
natively support metapopulation spatial structure. Epipack [@Maier2021] offers a
Python API for compartmental models but targets well-mixed populations without patch
connectivity. EMOD [@Bershteyn2018] is a high-fidelity agent-based platform with
spatial structure, but its complexity makes it unsuitable for rapid exploratory
modelling. GLEaMviz [@VandenBroeck2011] supports global air-travel metapopulation
models but is limited to a fixed model catalogue and a graphical interface. EpiModel
[@Jenness2018] is an R package for network epidemiology targeting sexual contact
networks rather than geographic metapopulation structure.

PatchSim occupies a distinct niche: a lightweight, disease-agnostic,
configuration-first Python framework for geographic metapopulation modelling with
user-defined compartmental transitions. It requires no specialised installation,
integrates with the scientific Python ecosystem (NumPy, SciPy, pandas, NetworkX), and
natively represents inter-patch movement through weighted network structures. A
modeller can switch between SIR and SEIR by editing a configuration file, not by
rewriting code.

# Software design

PatchSim's central design principle is that compartmental model structure is declared
in YAML rather than implemented in code. Transitions are expressed as arrow-map
expressions pairing source and target compartments with a rate formula:

```yaml
transitions:
  - S -> I: Beta * S * I / N
  - I -> R: Gamma * I
```

This is parsed into a system of differential equations (ODE mode, using
`scipy.integrate.solve_ivp` with RK45) or discrete finite-difference update rules
(timestep mode), applied simultaneously across all patches at each time step. The
spatial structure is represented as a NetworkX graph in which nodes are patches and
edge weights encode inter-patch movement rates. All inputs — patch populations, seed
infections, and network connectivity — are read from CSV files resolved relative to
the configuration file, enabling self-contained portable project directories.

Configurations are validated against an embedded JSON Schema before any computation
begins, catching errors in transition expressions, mismatched compartment names, and
missing files before a run starts. PatchSim integrates SALib [@Herman2017] for Sobol
sensitivity analysis and Morris screening, and supports random seed control and
timestamped logging for reproducibility. Three built-in YAML templates (SIR, SEIR,
SIRS) are bundled with the package and serve as starting points for new models.

# AI usage disclosure

GitHub Copilot was used during software development to assist with boilerplate code
generation and documentation drafting. All AI-generated content was reviewed and
verified by the authors. The scientific design, model architecture, transition
specification format, and validation logic were conceived and implemented by the human
authors.

# Acknowledgements

The authors thank the ARTPARK team at the Indian Institute of Science for
institutional support. This work was supported by the AI and Robotics Technology Park,
ARTPARK, IISc, Bengaluru.

# References
