# CLI Reference

Command-line interface reference for PatchSim.

## Overview

```bash
patchsim [COMMAND] [OPTIONS]
```

**Available commands:**
- `init` – Scaffold a new project
- `run` – Execute a simulation
- `validate` – Check configuration and inputs
- `list-models` – Show available templates
- `--version` – Display version

---

## `init`

Initialize a new PatchSim project with starter files.

### Usage

```bash
patchsim init PROJECT_NAME [--template TEMPLATE] [--force]
```

### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `PROJECT_NAME` | Yes | Directory name for the new project |

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--template` | `sir` | Starter model template: `sir`, `seir`, `sirs`, `sis` |
| `--force` | - | Overwrite target directory if it already exists |

### Examples

```bash
# Create a new SIR project (default)
patchsim init my-project

# Create with SEIR template
patchsim init my-project --template seir

# Overwrite existing directory
patchsim init my-project --force
```

### Output

Creates a directory with:
```
config.yaml                              # Model configuration
data/
  patch/                                 # Patch population files
  seeds/                                 # Initial conditions
  networks/                              # Network topology (if multi-patch)
output/                                  # (Empty until run)
```

---

## `validate`

Validate configuration file and input data for consistency.

### Usage

```bash
patchsim validate [-c CONFIG] [--schema] [--json]
```

### Options

| Option | Required | Description |
|--------|----------|-------------|
| `-c, --config` | Conditional* | Path to config YAML file |
| `--schema` | - | Print the configuration JSON Schema |
| `--json` | - | Output validation results as JSON |

*Required unless `--schema` is used.

### Examples

```bash
# Validate a configuration file
patchsim validate -c config.yaml

# Print the JSON Schema for config files
patchsim validate --schema

# Validate and output JSON
patchsim validate -c config.yaml --json
```

### Output (Plain Text)

```
Configuration is valid: config.yaml
```

### Output (JSON)

```json
{
  "ok": true,
  "config": "config.yaml",
  "model_name": "sample-sir-ode",
  "num_patches": 1,
  "patches": ["PatchA"]
}
```

### Validation Checks

- All required fields present (see [Configuration](configuration.md))
- Compartments defined and non-empty
- Transitions use valid compartments and parameters
- Patch file exists and is readable
- Seed file has correct structure
- Compartments match seed file columns
- Seed values sum to patch population
- Network file (if provided) is consistent

---

## `run`

Execute a simulation and generate outputs.

### Usage

```bash
patchsim run -c CONFIG [--json]
```

### Options

| Option | Required | Description |
|--------|----------|-------------|
| `-c, --config` | Yes | Path to configuration YAML file |
| `--json` | - | Output run summary as JSON |

### Examples

```bash
# Run a simulation
patchsim run -c config.yaml

# Run and capture output as JSON
patchsim run -c config.yaml --json
```

### Output (Plain Text)

```
Starting PatchSim simulation...
Model: sample-sir-ode
Parameters: beta=0.5, gamma=0.1
Simulation completed successfully.
```

### Output (JSON)

```json
{
  "ok": true,
  "config": "config.yaml",
  "model_name": "sample-sir-ode",
  "output_dir": "output",
  "csv_path": "output/runs/all_patches_sample-sir-ode_ode.csv",
  "plot_path": "output/plots/patch_timeseries_sample-sir-ode_ode.png",
  "num_patches": 1,
  "patches": ["PatchA"],
  "t_max": 100
}
```

### Output Files

Generates in `OutputDir` (from config):

| File | Description |
|------|-------------|
| `runs/all_patches_MODEL_ode.csv` | Time series data (time + all compartments) |
| `plots/patch_timeseries_MODEL_ode.png` | Line plot of all compartments over time |
| `logs/MODEL_run_*.log` | Simulation log (detailed diagnostics) |

See [Results](results.md) for output format details.

---

## `list-models`

Display available starter templates.

### Usage

```bash
patchsim list-models [--json]
```

### Options

| Option | Description |
|--------|-------------|
| `--json` | Output model list as JSON |

### Examples

```bash
# List available templates
patchsim list-models

# Output as JSON
patchsim list-models --json
```

### Output (Plain Text)

```
Built-in models and templates:
- seir (yaml-template)
- sir (yaml-template)
- sirs (yaml-template)
- sis (yaml-template)
```

### Output (JSON)

```json
{
  "models": [
    {"name": "seir", "kind": "yaml-template"},
    {"name": "sir", "kind": "yaml-template"},
    {"name": "sirs", "kind": "yaml-template"},
    {"name": "sis", "kind": "yaml-template"}
  ]
}
```

---

## Global Options

### `--version`

Display the installed PatchSim version.

```bash
patchsim --version
```

Output:
```
patchsim 0.1.0b1
```

---

## Common Errors and Solutions

| Error | Cause | Solution |
|-------|-------|----------|
| `Compartment mismatch` | Seed file columns don't match config | Update seed file or config compartments |
| `Unknown names in expression` | Typo in transition rate | Check parameter/compartment names in config |
| `Refusing to overwrite` | Target directory exists | Use `--force` or choose a different name |
| `Missing required field` | Config missing field (e.g., `TMax`) | See [Configuration](configuration.md) for required fields |
| `Population mismatch` | Seed values don't sum to population | Fix seed file CSV |
