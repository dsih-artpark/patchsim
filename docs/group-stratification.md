# Group stratification

PatchSim can divide each geographic patch into one categorical grouping. The group
labels and meaning are supplied by the modeller. Examples include:

- age bands;
- behavioural risk categories in sexually transmitted infection models;
- occupations such as construction workers;
- species; and
- intervention or exposure categories.

Grouped runs require `GroupFile`, `InteractionFile`, and `InteractionUnits`
together. The worked example below contains the complete files and config.

The interaction matrix may come from measurements, published estimates, or a documented
proxy. PatchSim validates its numeric and identifier contracts, but cannot establish
that the grouping or proxy is appropriate for the modelled transmission process.

## Group populations

`GroupFile` uses the columns `patch`, `group`, and `population` to allocate each
`PatchFile` population across the same group set.

Every patch/group pair must appear once. Values must be finite and non-negative, and
group populations must sum to the corresponding patch population. A zero-population
stratum is valid.

Group order follows the rows for the first patch in `PatchFile`. Rows for other patches
may use any order. Labels must match exactly across group, seed, and interaction files.

## Initial state

For grouped runs, add `group` after `patch` in `SeedFile`, followed by the
configured compartment columns.

Every patch/group pair must appear once. Its compartments must be finite,
non-negative, and sum to that stratum's population.

## Interaction matrix

`InteractionFile` is a static matrix shared by all patches. It uses
`focal_group`, `contributor_group`, and `weight`.

`focal_group` is the susceptible or exposed row. `contributor_group` is the infectious
column. This distinction matters for asymmetric matrices.

Weights must be finite and non-negative. Missing pairs are zero, duplicate pairs are
rejected, and each group must appear at least once in each role. PatchSim does not
normalize, symmetrize, or demographically adjust the matrix.

`InteractionUnits` records the supplied weights' meaning. Examples include
`contacts/person/day`, `partnerships/person/year`, or `dimensionless`. It is descriptive:
the modeller must ensure that the interaction weights, spatial weights, transmission
parameter, and simulation time unit are compatible.

## Composition

Let $W_{ij}$ be the spatial weight from focal patch $i$ to contributing patch $j$, and
let $M_{ab}$ be the interaction weight from focal group $a$ to contributing group $b$.
PatchSim computes:

$$
\lambda_{i,a}
=
\sum_j \sum_b
W_{ij} M_{ab}
\frac{I_{j,b}}{N_{j,b}}.
$$

A zero-population contributor has prevalence zero. For one patch, the spatial factor is
one. For multiple patches, the existing `NetworkFile` rules apply, including zero
spatial pressure when no network is supplied.

This equation assumes spatial and group mixing are separable and that the same group
matrix applies in every patch. If, for example, occupational mixing differs materially
by district, a single shared matrix is not an adequate representation.

## Validation

Run validation before simulation:

```bash
patchsim validate -c config.yaml
patchsim validate -c config.yaml --json > interaction-validation.json
```

The JSON result records:

- group labels and ordering;
- declared interaction units;
- SHA-256 of the exact interaction file;
- spatial and interaction row-sum ranges; and
- maximum local reciprocity residual.

For groups $a$ and $b$ in patch $i$, the reported residual compares:

$$
x = N_{i,a}M_{ab},
\qquad
y = N_{i,b}M_{ba},
\qquad
r = \frac{|x-y|}{\max(x,y)}.
$$

When $x=y=0$, $r=0$. Reciprocity is diagnostic only: it is expected for matrices
representing reciprocal person-to-person contacts, but not necessarily for directed
exposure coefficients or other proxies.

For a proxy matrix, retain its derivation and assumptions with the input files. Compare
results under defensible alternative matrices; successful numeric validation does not
validate the proxy.

## Worked example

This one-patch SIR example uses illustrative `higher_contact` and `general`
groups. The interaction values are not estimates for a specific population.

```text
group-example/
  config.yaml
  data/
    patches.csv
    groups.csv
    seeds.csv
    interactions.csv
```

`data/patches.csv`:

```csv
patch,population
community,1000
```

`data/groups.csv`:

```csv
patch,group,population
community,higher_contact,200
community,general,800
```

`data/seeds.csv`:

```csv
patch,group,S,I,R
community,higher_contact,199,1,0
community,general,800,0,0
```

`data/interactions.csv`:

```csv
focal_group,contributor_group,weight
higher_contact,higher_contact,8
higher_contact,general,4
general,higher_contact,1
general,general,3
```

`config.yaml`:

```yaml
ModelName: group-example
PatchFile: data/patches.csv
GroupFile: data/groups.csv
InteractionFile: data/interactions.csv
InteractionUnits: contacts/person/day
SeedFile: data/seeds.csv
NetworkFile: null
OutputDir: output/group-example

TMax: 61
Solver: ode
TimeStep: 1.0

compartments: [S, I, R]

Parameters:
  beta: 0.03
  gamma: 0.1

Transitions:
  "S -> I": "beta"
  "I -> R": "gamma * I"
```

Here `beta` is transmission probability per contact and `gamma` is per day.
From `group-example/`, validate and run:

```bash
patchsim validate -c config.yaml
patchsim validate -c config.yaml --json > interaction-validation.json
patchsim run -c config.yaml --json > run-summary.json
```

From a source checkout, replace `patchsim` with `uv run patchsim`. The group order
makes `I_0_0` the infectious `higher_contact` population and `I_0_1` the
infectious `general` population.

At time zero, infectious prevalence is `1 / 200 = 0.005` in `higher_contact`
and zero in `general`. The initial pressures are:

$$
\lambda_{\mathrm{higher}} = 8(0.005) + 4(0) = 0.04,
$$

$$
\lambda_{\mathrm{general}} = 1(0.005) + 3(0) = 0.005.
$$

Initial infection flows are therefore `0.03(199)(0.04) = 0.2388` and
`0.03(800)(0.005) = 0.12`. These calculations verify orientation and units;
they do not validate the illustrative interaction values.

## Outputs and limits

Ungrouped columns remain unchanged, such as `I_0`. Grouped columns use
`COMPARTMENT_PATCH_GROUP`, such as `I_0_1`. Patch indices follow `PatchFile`; group
indices follow `GroupFile`. Plots keep one panel per patch and label lines with both the
compartment and group.

The current feature supports one grouping axis and one shared, static interaction
matrix. Patch-specific matrices, time-varying interactions, group-specific mobility,
multiple simultaneous group axes, and automatic interaction generators are not
implemented. `PatchParameters` still applies to every group within its patch.
