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
specified in a YAML configuration file, where each transition is written as
`source -> target` with a rate expression. PatchSim provides adaptive ordinary
differential equation (ODE) and
explicit-Euler solvers, Sobol sensitivity
analysis via SALib [@Herman2017], bounded least-squares calibration, and command-line
and Python interfaces. Its configuration-first scope covers compartment-transfer models
over supplied patch and group interaction data; scientific assumptions remain the
modeller's responsibility. PatchSim is available on PyPI (`pip install patchsim`), with
source code and documentation linked from the project record [@patchsim].

# Statement of need

Mathematical modelling of infectious diseases is essential for understanding
transmission dynamics, evaluating interventions, and informing public health policy
[@Keeling2008]. Metapopulation models — which capture the heterogeneity of disease
spread across connected spatial units [@Grenfell1997] — are particularly valuable for
diseases where population movement shapes outbreak trajectories, including
foot-and-mouth disease in livestock [@Keeling2001; @Tildesley2006; @GuyverFletcher2025],
dengue in urban landscapes [@Wesolowski2015], and respiratory infections across
administrative regions [@Balcan2009; @Gatto2020].

Implementing such models commonly requires custom scientific code or a platform tied to
a particular modelling abstraction. This adds software work to the scientific tasks of
specifying transitions, preparing spatial inputs, and checking model assumptions.

PatchSim reduces this cost through a configuration-first design in which model
structure — compartments, transitions, and parameters — is declared in a YAML file
rather than implemented in code. The intended users are epidemiological modellers,
researchers studying spatial disease dynamics, and public health analysts who need to
rapidly prototype and compare scenarios across spatial configurations and disease
systems. PatchSim was developed at ARTPARK, the AI and Robotics Technology Park at the
Indian Institute of Science (IISc), to support active modelling work on livestock
disease dynamics in India.

# State of the field

Existing software spans several related abstractions: EpiFire [@Hladish2012] and
EpiModel [@Jenness2018] address contact-network epidemiology; epipack [@Maier2021]
supports compartmental, stochastic, and network models; EMOD [@Bershteyn2018] uses an
individual-based architecture; and GLEaM [@Balcan2010; @VandenBroeck2011] represents
global mobility-driven metapopulations.

Metapopulation frameworks closer to PatchSim
include MetaWards [@Woods2022], MEmilio [@Bicker2026], SimInf [@Widgren2019], MetaCast
[@Grunnill2024], and flepiMoP [@Lemaitre2024]. MetaWards is a stochastic metapopulation
framework originally built for the electoral wards of Great Britain. MEmilio is a C++
core with Python bindings that offers graph-ODE and agent-based models. SimInf declares
stochastic transitions as strings over nodes
linked by scheduled livestock movements. MetaCast broadcasts a user-written Python ODE
over subpopulations and bundles Latin hypercube sensitivity analysis. flepiMoP is an R
and Python forecasting pipeline that declares compartments and transitions in YAML over
a mobility matrix.

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

We did not build on the three closest tools for the following reasons. flepiMoP already
declares compartments and transitions in YAML over a mobility matrix and offers a
deterministic solver. At the time of writing it lacks variance-based sensitivity
indices and bounded least-squares calibration, fits by Markov chain Monte Carlo, and
its R and Python pipeline is larger than the single package the Karnataka work needed.
MetaCast
expresses the model as Python code, so a change in compartments is a code change.
SimInf is stochastic only and moves animals by scheduled events rather than through a
contact matrix, so it does not provide the deterministic comparison over district
contact data that the Karnataka work needed.

PatchSim contributes the combination none of the three offers in one Python package
with a command-line interface: a compartment
graph and rate expressions declared in configuration, a user-supplied weighted patch
matrix, categorical group stratification, a deterministic solver, Sobol sensitivity
analysis, and bounded calibration. A modeller can switch between compartment
structures, solvers, spatial networks, and group interaction inputs without rewriting
the runtime.

# Software design

PatchSim's central design principle is that compartmental model structure is declared
in YAML rather than implemented in code. Transitions are written as
`source -> target` keys, each paired with a rate formula:

```yaml
compartments: [S, I, R]
Parameters: {beta: 0.2, gamma: 0.1}
Transitions:
  "S -> I": "beta"
  "I -> R": "gamma * I"
```

A rate expression that does not name its source compartment is multiplied by that
compartment, so `"S -> I": "beta"` gives a local flow of $\beta S_i$. For focal patch
$i$ in an ungrouped multi-patch model, the infectious pressure is
$\lambda_i(t)=\sum_j W_{ij}I_j(t)/N_j$, and the infection flow in this example is
$\beta S_i\lambda_i$; a single ungrouped patch applies no coupling, and grouped models
also multiply by the group interaction matrix. This is a residence-based
coupling in the sense of the Lagrangian multi-patch models of @Sattenspiel1995 and
@Citron2021. Patch populations stay fixed at $N_i$, and residents of patch $i$ are
exposed to the resident prevalence of each patch $j$ in proportion to $W_{ij}$; no
population moves. The prevalence at patch $j$ counts residents only, whereas the full
Sattenspiel–Dietz model counts visiting infectives and visitors in both numerator and
denominator.

A supplied $W$ is applied as given and its rows are not normalized, so the scale of $W$
multiplies $\beta$. The built-in contact-network generator can produce row-normalized
weights with a self-share diagonal. Whether a supplied matrix should do the same, and
whether it comes from a gravity kernel or observed movement data, is the modeller's
choice. The spatial coupling is bound to the compartment names `S`, `I`, and `E`, and
$\lambda_i$ counts `I` alone.

ODE mode uses `scipy.integrate.odeint`, SciPy's interface to the LSODA solver, while
discrete mode uses deterministic explicit Euler. Both evaluate the same derivative function; discrete
results should be checked at successively smaller time steps. As a deterministic patch
model, PatchSim does not represent stochastic fade-out or invasion probability in small
patches [@Riley2007; @Ball2015]. The current network is fixed from day zero. Patch
populations, initial states, network weights, and optional group interactions are read
from CSV files resolved relative to the configuration.

Validation checks required fields, identifiers, finite input values, population totals,
matrix dimensions, transition endpoints, and an arithmetic-only expression language;
these checks do not establish scientific validity. The JSON Schema is available for
editors and external tooling. Reproducible analysis workflows provide seeded first- and
total-order Sobol indices [@Sobol2001] and bounded multi-start calibration with input
hashes and diagnostics. Built-in SIR, SEIR, SIRS, and SIS templates and complete worked examples
are included in the documentation. This bounded modularity is configuration- and
API-based; adding a numerical solver still requires Python development. PatchSim is
released under GPL-3.0.

# Research impact statement

PatchSim is in current research use at ARTPARK for foot-and-mouth disease
vaccination-scenario analysis across Karnataka districts, where it replaced a
project-specific implementation. The workflow models cattle and buffalo populations
stratified by age group and species. That analysis is unpublished at the time of
writing, and no publication enabled by the package has yet appeared. Published spatial
analysis of foot-and-mouth disease in India includes a state-level Bayesian space-time
model of reported outbreaks, which estimated a roughly 50% lower outbreak risk in states
covered by the vaccination programme [@Gunasekera2022]. The Karnataka analysis applies
a mechanistic patch model to the between-district allocation question.

Near-term significance rests on what a third party can install, run, and check. The
package has been developed in public since April 2025,
with two tagged releases on PyPI and an archived release record [@patchsim]. Continuous
integration runs the automated test suite, linting, and a documentation build on every
pull request and push to the main branch. The tests cover configuration validation, the
expression
evaluator, both solvers, group stratification, the sensitivity and calibration
workflows, and the command-line interface. The seeded sensitivity and calibration
outputs described above carry input
hashes and diagnostics, so a third party can re-run and check a study. The documentation
on Read the Docs includes worked examples.

# AI usage disclosure

GitHub Copilot, OpenAI Codex CLI (GPT-5), and Anthropic Claude Code assisted with code,
tests, and documentation, and with language checks and paraphrasing of the authors'
draft of this paper. CodeRabbit assisted with code review. Model
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
