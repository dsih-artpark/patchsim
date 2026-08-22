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
  - name: Shreya Mukherjee
    orcid: 0009-0008-3629-0376
    affiliation: 1
  - name: Sneha S
    orcid: 0009-0009-1854-6443
    affiliation: 1
affiliations:
  - index: 1
    name: AI and Robotics Technology Park (ARTPARK), Indian Institute of Science, Bengaluru, India
date: 22 August 2026
bibliography: paper.bib
---

# Summary

PatchSim is an open-source Python framework for simulating infectious disease dynamics
across spatially structured populations. It implements a patch-based metapopulation
approach in which a landscape is partitioned into discrete geographical units
(patches) — such as subdistricts, districts, or ecological zones — connected by a
weighted contact network. Within each patch, disease
progression is described by user-defined compartmental models (e.g., SIR, SEIR, SIRS)
specified through a YAML configuration file using a declarative arrow-map transition
syntax. PatchSim provides adaptive ODE and explicit-Euler solvers, Sobol sensitivity
analysis via SALib [@Herman2017], bounded least-squares calibration, and command-line
and Python interfaces. Its configuration-first scope covers compartment-transfer models
over supplied patch and group interaction data; scientific assumptions remain the
modeller's responsibility. PatchSim is available on PyPI (`pip install patchsim`), with
source code and documentation linked from the project record [@patchsim].

# Statement of need

Mathematical modelling of infectious diseases is essential for understanding
transmission dynamics, evaluating interventions, and informing public health policy
[@Keeling2008]. Metapopulation models — which capture the heterogeneity of disease
spread across connected spatial units — are particularly valuable for diseases where
population movement shapes outbreak trajectories, including foot-and-mouth disease in
livestock, dengue in urban landscapes, and respiratory infections across administrative
regions [@Grenfell1997; @Balcan2009].

Implementing such models commonly requires custom scientific code or a platform tied to
a particular modelling abstraction. This adds software work to the scientific tasks of
specifying transitions, preparing spatial inputs, and checking model assumptions.

PatchSim addresses this gap through a configuration-first design in which model
structure — compartments, transitions, and parameters — is declared in a YAML file
rather than implemented in code. The intended users are epidemiological modellers,
researchers studying spatial disease dynamics, and public health analysts who need to
rapidly prototype and compare scenarios across spatial configurations and disease
systems. PatchSim was developed at ARTPARK, IISc, to support active modelling work on
livestock disease dynamics in India. The current package is used locally for ongoing
foot-and-mouth disease vaccination-scenario analysis across Karnataka districts,
replacing an earlier project-specific implementation; the workflow models cattle and
buffalo populations stratified by age group and species.

# State of the field

Existing software spans several related abstractions: EpiFire [@Hladish2012] and
EpiModel [@Jenness2018] address contact-network epidemiology; epipack [@Maier2021]
supports compartmental, stochastic, and network models; EMOD [@Bershteyn2018] uses an
individual-based architecture; and GLEaMviz [@VandenBroeck2011] represents global
mobility-driven metapopulations. PatchSim instead centres a configuration-defined
compartment system over a user-supplied weighted patch matrix and optional categorical
group interactions.

PatchSim is an independent implementation, but it shares both its name and its
modelling lineage with an earlier package developed at the Network Systems Science and
Advanced Computing division of the University of Virginia [@VenkatramananPatchSim].
That package implements a metapopulation SEIR model and has been applied to seasonal
influenza vaccine allocation [@Venkatramanan2019] and
to county-scale influenza forecasting in the United States [@Venkatramanan2021]. The
earlier implementation was also used for age- and district-stratified COVID-19 vaccine
allocation in Karnataka, as documented in a medRxiv preprint [@Adiga2021]. The approach
was subsequently taken up in collaborative modelling work at the Indian Institute of
Science and applied to foot-and-mouth disease in Indian livestock
populations, which is the setting from which the present work grew. Where the earlier
package provides a fixed SEIR compartmental structure serving a specific forecasting
application, the framework described here generalises that approach: compartmental
structure is treated as user-supplied configuration, so SIR, SEIR, SIRS, or arbitrary
user-defined compartment sets are expressed in YAML without modifying source code. The
name is retained to acknowledge this lineage. The earlier implementation is distributed
as `NSSAC/PatchSim`; this package is distributed as `patchsim` from
`dsih-artpark/patchsim`.

PatchSim is a lightweight, configuration-first Python framework for geographic
metapopulation modelling with user-defined compartmental transitions. A modeller can
switch between supported compartment structures, solvers, spatial networks, and group
interaction inputs without rewriting the runtime.

# Software design

PatchSim's central design principle is that compartmental model structure is declared
in YAML rather than implemented in code. Transitions are expressed as arrow-map
expressions pairing source and target compartments with a rate formula:

```yaml
compartments: [S, I, R]
Parameters: {beta: 0.2, gamma: 0.1}
Transitions:
  "S -> I": "beta"
  "I -> R": "gamma * I"
```

For focal patch $i$, the infectious pressure is
$\lambda_i(t)=\sum_j W_{ij}I_j(t)/N_j(t)$; the infection flow in this example is
$\beta S_i\lambda_i$. PatchSim does not normalize the supplied weights. ODE mode uses
`scipy.integrate.odeint`, SciPy's interface to LSODA, while discrete mode uses
deterministic explicit Euler. Both evaluate the same derivative function; discrete
results should be checked at successively smaller time steps. The current network is
fixed from day zero. Patch populations, initial states, network weights, and optional
group interactions are read from CSV files resolved relative to the configuration.

Validation checks required fields, identifiers, finite input values, population totals,
matrix dimensions, transition endpoints, and an arithmetic-only expression language;
these checks do not establish scientific validity. The JSON Schema is available for
editors and external tooling. Reproducible analysis workflows provide seeded first- and
total-order Sobol indices and bounded multi-start calibration with input hashes and
diagnostics. Built-in SIR, SEIR, SIRS, and SIS templates and complete worked examples
are included in the documentation. This bounded modularity is configuration- and
API-based; adding a numerical solver still requires Python development. PatchSim is
released under GPL-3.0.

# AI usage disclosure

GitHub Copilot, OpenAI Codex CLI (GPT-5), and Anthropic Claude Code assisted with code,
tests, documentation, and copy-editing. CodeRabbit assisted with code review. Model
versions were not consistently retained for historical work. The authors reviewed,
edited, tested, and validated all assisted outputs and made the scientific and
architectural decisions.

# Acknowledgements

The authors thank Prof. Rajesh Sundaresan and Prof. Siva Athreya for their guidance on
the modelling work from which this framework grew, and Srinivasan Venkatramanan, whose
earlier metapopulation modelling work informed the approach taken here. The authors
also thank the ARTPARK team at the Indian Institute of Science for institutional
support. This work was supported by the AI and Robotics Technology Park, ARTPARK,
IISc, Bengaluru. The authors declare no conflicts of interest.

# References
