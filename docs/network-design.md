# Network design

PatchSim loads a dense matrix from a CSV file. The file is used for
network-coupled `S -> I` and `S -> E` transitions in multi-patch ODE runs.

## File contract

```csv
day,source,target,weight
0,A,A,0.9
0,A,B,0.1
0,B,A,0.1
0,B,B,0.9
```

- `day` must be `0` for a row to be loaded.
- `source` and `target` must match identifiers in `PatchFile`.
- `weight` must be non-negative.
- Missing pairs have weight zero.
- Repeated pairs are overwritten by the last day-zero row.

Rows for later days are currently ignored; time-varying networks are not yet
implemented.

## Matrix orientation

PatchSim stores each row as:

$$
W[\mathrm{source}, \mathrm{target}] = \mathrm{weight}.
$$

It then computes the infectious pressure for focal patch $i$ as:

$$
\lambda_i = \sum_j W_{ij}\frac{I_j}{N_j}.
$$

For the example above:

$$
\lambda_A =
0.9\frac{I_A}{N_A}
+ 0.1\frac{I_B}{N_B}.
$$

In other words, the current `source` field selects the focal row and `target`
selects an infectious contributor to that row. Use this runtime convention when
creating custom files.

## Weight conventions

PatchSim uses weights exactly as supplied. It does not normalize rows, enforce an
upper bound, or add self-loops.

A common input convention is:

- one row per ordered patch pair that should contribute;
- a self-loop for within-patch pressure; and
- row sums of 1 when weights are intended as a weighted average.

These are modeling choices, not loader requirements.

## Patch order

Matrix indices follow `PatchFile` row order. Output suffixes use the same order:
the first patch is `_0`, the second is `_1`, and so on. Keep patch identifiers
unique and identical across patch, seed, and network files.

## No network file

If `NetworkFile` is omitted or `null`, PatchSim creates a zero matrix. In a
multi-patch run, network infectious pressure is then zero. A one-patch run does
not use network coupling.

Use [Contact generation](contact-generation.md) to build a validated day-zero
network CSV from centroid, GeoJSON, or shapefile input. The generated file still
uses this page's matrix orientation and is supplied through `NetworkFile`.
