# Configuration

PatchSim simulations are configured using YAML files. This guide documents every configuration field with examples.

## Configuration Structure

A PatchSim config file has four main sections:

1. **I/O & Metadata** – Files and output location
2. **Model Definition** – Compartments, parameters, transitions
3. **Data** – Initial conditions and network topology
4. **Simulation Control** – Time steps and solver settings

---

## I/O and Metadata

### ModelName

**Type:** String  
**Required:** Yes  
**Description:** Unique identifier for the simulation. Used in output filenames.

```yaml
ModelName: my-sir-model
```

### OutputDir

**Type:** File path  
**Required:** Yes  
**Description:** Directory where results (CSV, plots, logs) are saved. Created if it doesn't exist.

```yaml
OutputDir: output/results
```

### PatchFile

**Type:** File path  
**Required:** Yes  
**Description:** CSV file defining patch (region) names and total populations. Must have columns: `patch`, `population`.

**Example file** (`data/patch/patches.csv`):
```csv
patch,population
PatchA,1000
PatchB,500
```

```yaml
PatchFile: data/patch/patches.csv
```

### SeedFile

**Type:** File path  
**Required:** Yes  
**Description:** CSV file with initial compartment values for each patch. Must have columns: `patch` + all compartment names.

**Example file** (`data/seeds/initial.csv`):
```csv
patch,S,I,R
PatchA,950,50,0
PatchB,450,50,0
```

Constraints:
- One row per patch
- Compartment values must be non-negative
- Sum of values per patch must equal population in PatchFile

```yaml
SeedFile: data/seeds/initial.csv
```

### NetworkFile

**Type:** File path or `null`  
**Required:** No (default: `null`)  
**Description:** CSV file describing transmission links between patches. Required for multi-patch simulations. Ignored for single-patch models.

**Example file** (`data/networks/network.csv`):
```csv
source,target,weight
PatchA,PatchB,0.2
PatchB,PatchA,0.1
```

Columns:
- `source` – Source patch (origin of transmission)
- `target` – Target patch (receives transmission)
- `weight` – Transmission weight (0–1, higher = more transmission)

```yaml
NetworkFile: data/networks/network.csv
```

Or for single-patch:

```yaml
NetworkFile: null
```

---

## Model Definition

### compartments

**Type:** List of strings  
**Required:** No (inferred from SeedFile if omitted)  
**Description:** List of epidemiological state variables. If not provided, inferred from SeedFile columns.

```yaml
compartments:
  - S    # Susceptible
  - E    # Exposed
  - I    # Infectious
  - R    # Recovered
```

Common compartments:
- `S` – Susceptible
- `E` – Exposed
- `I` – Infectious
- `R` – Recovered
- `D` – Deceased
- `V` – Vaccinated

### Parameters

**Type:** Dictionary of name → number  
**Required:** Yes  
**Description:** Global per-capita transition rates. Available in all transition expressions. Can be overridden per-patch in `PatchParameters`.

**Example:**
```yaml
Parameters:
  beta: 0.5      # Transmission rate (contacts per individual per day)
  sigma: 0.2     # Progression rate (1 / incubation period)
  gamma: 0.1     # Recovery rate (1 / infectious period)
```

**Key point:** Parameters are per-capita rates. See [Rate Multiplication](rate-multiplication.md) for how rates interact with compartments.

### Transitions

**Type:** Dictionary of `compartment -> compartment` → expression  
**Required:** Yes  
**Description:** Define movement between compartments using arrow syntax and mathematical expressions.

**Syntax:** `SOURCE -> TARGET: expression`

**Example (SIR):**
```yaml
Transitions:
  S -> I: "beta * S * I"        # New infections
  I -> R: "gamma * I"           # Recoveries
```

**Example (SEIR):**
```yaml
Transitions:
  S -> E: "beta * S * I"        # Infection
  E -> I: "sigma * E"           # Progression to infectious
  I -> R: "gamma * I"           # Recovery
```

**Expression rules:**
- Use compartment names (e.g., `S`, `I`) for counts
- Use parameter names (e.g., `beta`, `gamma`) for rates
- Python operators allowed: `*`, `+`, `-`, `/`, `**`, `and`, `or`
- Only identifiers that are compartments or parameters allowed
- Rates are automatically multiplied by source compartment (see [Rate Multiplication](rate-multiplication.md))

**Valid expressions:**
```yaml
S -> I: "beta"           # Rate multiplied by S automatically
S -> I: "beta * I"       # Manual interaction term
E -> I: "sigma"
I -> R: "gamma"
R -> S: "rho"            # Waning immunity
```

**Invalid expressions:**
```yaml
S -> I: "unknown_param"  # Error: unknown parameter
S -> I: "S / 0"          # Error: division issues caught at runtime
```

---

## Patch-Specific Parameters (Optional)

### PatchParameters

**Type:** Dictionary of patch name → parameter overrides  
**Required:** No  
**Description:** Override global parameters for specific patches. Useful for spatial heterogeneity.

```yaml
Parameters:
  beta: 0.5
  gamma: 0.1

PatchParameters:
  PatchA:
    beta: 0.8    # Higher transmission in PatchA
  PatchB:
    beta: 0.3    # Lower transmission in PatchB
```

Constraints:
- Keys must match patch names from PatchFile
- Only override parameters you need to change
- Overrides are merged with global Parameters

---

## Simulation Control

### TMax

**Type:** Positive integer  
**Required:** Yes  
**Description:** Number of time steps to simulate.

```yaml
TMax: 100
```

This simulates from time 0 to 99 (100 steps total).

---

## Complete Examples

### Minimal SIR (Single Patch)

```yaml
ModelName: simple-sir
OutputDir: output/sir
PatchFile: data/patch/single.csv
SeedFile: data/seeds/sir-init.csv
NetworkFile: null

compartments:
  - S
  - I
  - R

Parameters:
  beta: 0.5
  gamma: 0.1

Transitions:
  S -> I: "beta"
  I -> R: "gamma"

TMax: 100
```

### SEIR with Waning Immunity

```yaml
ModelName: seir-waning
OutputDir: output/seir
PatchFile: data/patch/patches.csv
SeedFile: data/seeds/seir-init.csv
NetworkFile: null

compartments:
  - S
  - E
  - I
  - R

Parameters:
  beta: 0.5    # Transmission rate
  sigma: 0.2   # Progression (1/latency)
  gamma: 0.1   # Recovery (1/infectious)
  rho: 0.01    # Waning immunity

Transitions:
  S -> E: "beta * I"   # Exposure
  E -> I: "sigma"      # Progression
  I -> R: "gamma"      # Recovery
  R -> S: "rho"        # Loss of immunity

TMax: 365
```

### Multi-Patch with Heterogeneous Transmission

```yaml
ModelName: network-sir
OutputDir: output/network
PatchFile: data/patch/regions.csv
SeedFile: data/seeds/network-init.csv
NetworkFile: data/networks/routes.csv

compartments:
  - S
  - I
  - R

Parameters:
  beta: 0.3
  gamma: 0.1

PatchParameters:
  Urban:
    beta: 0.6      # Higher transmission in cities
  Rural:
    beta: 0.2      # Lower transmission in rural areas

Transitions:
  S -> I: "beta * I"
  I -> R: "gamma"

TMax: 200
```

---

## Validation

Before running, validate your configuration:

```bash
patchsim validate -c config.yaml
```

This checks:
- ✅ All required fields present
- ✅ Files exist and are readable
- ✅ Compartments match seed file
- ✅ Transitions reference valid compartments/parameters
- ✅ Seed values are consistent

### Common Validation Errors

| Error | Cause | Fix |
|-------|-------|-----|
| `Missing required field 'TMax'` | TMax not in config | Add `TMax: N` |
| `Compartment mismatch` | Seed columns ≠ compartments | Align seed file with compartments |
| `Unknown names in expression` | Typo in transition | Check parameter/compartment spelling |
| `Refusing to overwrite` | OutputDir exists | Use different directory or `--force` |

---

## Next Steps

- **[Getting Started](getting-started.md)** – Hands-on walkthrough
- **[CLI Reference](cli-reference.md)** – Command documentation
- **[Mathematical Model](mathematical-model.md)** – Equations and rate multiplication
- **[Network Design](network-design.md)** – Multi-patch topology
