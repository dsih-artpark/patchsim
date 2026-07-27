# Contact generation

`patchsim generate-contacts` converts spatial region data into the day-zero edge
CSV accepted by `NetworkFile`. It also writes a JSON validation report describing
units, parameters, distance and weight ranges, normalization, and the CSV hash.

The command is a preprocessing tool. It does not change a simulation config or
generate hidden project state.

## Inputs

Centroid CSV input works with the base installation:

```csv
patch,lat,lon,population
A,12.97,77.59,1000
B,13.34,77.10,500
```

Latitude and longitude are WGS84 degrees. Identifiers must be non-empty and
unique after conversion to stripped strings. At least two regions are required.

GeoJSON and shapefile input require the optional geo dependencies:

```bash
python -m pip install "patchsim[geo]"
```

From a source checkout:

```bash
uv sync --extra geo --frozen
```

Pass the identifier and, for gravity, population column names explicitly.
PatchSim does not guess them.

## Distance kernel

For distinct regions $i$ and $j$:

$$
W_{ij} =
\frac{s}{\max(d_{ij}, d_{\min})^\alpha},
$$

where:

- $d_{ij}$ is haversine distance in kilometres;
- Earth radius is `6371.0088 km`;
- $\alpha$ is the dimensionless `--decay`;
- $d_{\min}$ is `--min-distance-km`; and
- $s$ is `--scale` for unnormalized output.

Generate a dimensionless row-normalized network with an 80% within-patch share:

```bash
patchsim generate-contacts regions.csv contacts.csv \
  --id-column patch \
  --kernel distance \
  --decay 2 \
  --min-distance-km 0.001 \
  --normalize row \
  --self-share 0.8
```

## Gravity kernel

The gravity kernel additionally uses positive population values:

$$
W_{ij} =
\frac{s\,P_iP_j}{\max(d_{ij}, d_{\min})^\alpha}.
$$

```bash
patchsim generate-contacts regions.csv gravity.csv \
  --id-column patch \
  --population-column population \
  --kernel gravity \
  --decay 2 \
  --min-distance-km 0.001 \
  --normalize row \
  --self-share 0.8
```

With row normalization, a global scale cancels and is therefore not accepted.

## Normalization and units

Choose the normalization mode explicitly.

### Row normalized

`--normalize row` requires `--self-share` in `[0, 1)`. Each off-diagonal
row is scaled to sum to `1 - self-share`, and the diagonal is set to
`self-share`. Final weights are dimensionless and every row sums to one.

This is appropriate when each row is intended as a weighted average of
infectious prevalence. It does not prove that the chosen self-share or decay is
epidemiologically correct.

### Unnormalized

`--normalize none` requires both `--scale` and `--self-weight`:

```bash
patchsim generate-contacts regions.csv raw-contacts.csv \
  --id-column patch \
  --kernel distance \
  --decay 2 \
  --scale 1 \
  --min-distance-km 0.001 \
  --normalize none \
  --self-weight 1
```

The diagonal has the same units as the raw off-diagonal weights. Distance
weights have units `scale / km**decay`; gravity weights have units
`scale * population**2 / km**decay`. PatchSim uses them as coefficients, not
probabilities.

## Polygon centroids

Point geometries are transformed directly to EPSG:4326. Polygon and
MultiPolygon inputs require a projected CRS selected for the dataset:

```bash
patchsim generate-contacts regions.geojson contacts.csv \
  --id-column district_id \
  --kernel distance \
  --decay 2 \
  --min-distance-km 0.001 \
  --normalize row \
  --self-share 0.8 \
  --centroid-crs EPSG:32643
```

PatchSim computes centroids in that projected CRS, then transforms the points to
EPSG:4326. It rejects missing CRS metadata, invalid geometries, mixed point and
polygon inputs, and geographic values supplied as `--centroid-crs`.

Selecting an appropriate projection remains a modelling decision. Inspect the
resulting centroids and reported distance range before using the network.

## Outputs

For `contacts.csv`, the command writes:

```text
contacts.csv
contacts.csv.validation.json
```

The CSV is ordered by source region then target region:

```csv
day,source,target,weight
0,A,A,0.8
0,A,B,0.2
0,B,A,0.2
0,B,B,0.8
```

`source` is the focal row and `target` is the infectious contributor, matching
the runtime convention. `day` is always zero because the current runtime loads
only day-zero rows.

The report records:

- source and centroid CRS;
- kilometre distance convention and Earth radius;
- kernel parameters and raw/final weight units;
- distance, off-diagonal weight, row-sum, and diagonal ranges;
- pairs affected by the distance floor;
- raw and final symmetry; and
- SHA-256 of the exact CSV bytes.

The CSV/report replacement is not atomic across both files. Compare the report's
`csv_sha256` with the CSV hash after an interrupted write or file transfer:

```bash
sha256sum contacts.csv
```

On macOS:

```bash
shasum -a 256 contacts.csv
```

Without `--force`, generation fails if either output exists. With `--force`,
both are replaced; the hash makes an interrupted mismatch detectable.
Do not run concurrent generators against the same output path.

## Use in a simulation

Reference the generated CSV normally:

```yaml
NetworkFile: data/networks/contacts.csv
```

Before running:

1. confirm identifier order and values against `PatchFile`;
2. inspect the minimum, median, and maximum distances;
3. inspect weight dynamic range and row sums;
4. justify decay, scale, and self-coupling from domain evidence; and
5. retain the source data, CSV, report, and command with the result.

PatchSim validates numeric and file invariants. It cannot determine whether a
kernel is physically appropriate for a disease, population, or movement process.
