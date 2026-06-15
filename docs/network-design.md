# Network Design

Networks are represented as weighted directed graphs.

Example CSV:

day,source,target,weight
0,A,A,0.99
0,A,B,0.01
0,B,B,0.995
0,B,A,0.005

Each row specifies the contribution of infections in the source patch
to the force of infection in the target patch. The `day` column selects when
the weights take effect; use `0` for a static network.

## Weights

Weights are used exactly as provided — they are **not** auto-normalized. The
force of infection for patch *i* is `Σⱼ Wᵢⱼ · Iⱼ/Nⱼ`, so a common convention is
to make each source patch's outgoing weights sum to 1 (including a self-loop
`source == target` for within-patch transmission), as in the example above.

Generating networks from geometry
---------------------------------

PatchSim provides a generator CLI to create network CSVs from GeoJSON or
Shapefiles. The generated CSVs follow the same `day,source,target,weight`
convention and are immediately consumable by the simulation as long as the
`source`/`target` identifiers exactly match the patch names in your `PatchFile`.

Supported kernels include:
- `gravity`: population-based gravity kernel (requires a `population` field)
- `distance`: distance-decay kernel (no population required)

Use the `generate-contacts` command (see `python -m patchsim generate-contacts --help`) to
create one or more kernel files. When generating multiple kernels the CLI writes
files named `contacts-<method>.csv` by default.
