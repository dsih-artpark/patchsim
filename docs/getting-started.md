# Getting Started

A quick hands-on walkthrough from installation to your first simulation results.

## Installation

### From PyPI (Recommended)

```bash
pip install patchsim
```

### For Contributors

Clone the repository and install in editable mode:

```bash
git clone https://github.com/dsih-artpark/patchsim.git
cd patchsim
pip install -e ".[dev]"
```

## Your First Simulation (5 minutes)

### Step 1: Create a Project

Initialize a new PatchSim project with the SIR template:

```bash
patchsim init my-first-sim
cd my-first-sim
ls
```

You should see:

```
config.yaml         # Model and simulation configuration
data/
  networks/         # Network files (transmission between patches)
  patch/            # Patch populations
  seeds/            # Initial compartment values
output/             # (Created when you run the simulation)
```

### Step 2: Examine the Configuration

Open `config.yaml` to see the model structure:

```yaml
ModelName: sample-sir-ode
TMax: 100

PatchFile: data/patch/sample-sir-ode-patch-population.csv
SeedFile: data/seeds/sample-sir-ode-patchA-2.csv
NetworkFile: null  # Single patch, no network
OutputDir: output

# Model compartments (S, I, R)
compartments:
  - S
  - I
  - R

# Global parameters (all patches use these unless overridden)
Parameters:
  beta: 0.5    # Transmission rate
  gamma: 0.1   # Recovery rate

# Transitions: how individuals move between compartments
Transitions:
  S -> I: beta * S * I  # New infections
  I -> R: gamma * I     # Recoveries
```

**Key concepts:**
- **Compartments:** The state variables (S=Susceptible, I=Infectious, R=Recovered)
- **Parameters:** Per-capita rates (beta=contacts/day, gamma=1/duration)
- **Transitions:** Define movement between compartments using rate expressions
- **Rate multiplication rule:** `flow = rate × source_compartment`. See [Rate Multiplication](rate-multiplication.md) for details.

### Step 3: Validate Your Configuration

Before running, check that inputs are consistent:

```bash
patchsim validate -c config.yaml
```

This checks:
- ✅ All required fields present
- ✅ Compartments match seed file columns
- ✅ Transitions reference valid compartments and parameters
- ✅ Patch populations match seed values
- ✅ Network file (if provided) is consistent

### Step 4: Run the Simulation

Execute the simulation:

```bash
patchsim run -c config.yaml
```

You should see logging output:

```
Starting PatchSim simulation...
Model: sample-sir-ode
Parameters: beta=0.5, gamma=0.1
Simulation completed successfully.
```

### Step 5: Examine Results

Check the output directory:

```bash
ls output/runs/       # CSV file with time series
ls output/plots/      # PNG plots of compartments over time
```

**CSV output** (`output/runs/all_patches_sample-sir-ode_ode.csv`):

```csv
time,S_0,I_0,R_0
0.0,990,10,0
1.0,975.2,20.1,4.7
2.0,958.4,28.9,12.7
...
```

**PNG plots** show S, I, R trajectories for each patch.

### Step 6: Modify the Model

Now change the model to SEIR (add exposed compartment):

```bash
patchsim init my-seir --template seir
cd my-seir
patchsim validate -c config.yaml
patchsim run -c config.yaml
```

The SEIR template adds an E (exposed) compartment and demonstrates model customization. See [Configuration](configuration.md) for all fields.

## JSON Output Mode

For integration with pipelines, use `--json` flags:

```bash
# Validate and get machine-readable output
patchsim validate -c config.yaml --json

# Run and get summary in JSON
patchsim run -c config.yaml --json

# List available model templates
patchsim list-models --json
```

## Next Steps

- **[Configuration Reference](configuration.md)** – All config fields and options
- **[CLI Reference](cli-reference.md)** – Command reference and examples
- **[Mathematical Model](mathematical-model.md)** – Equations and how rate multiplication works
- **[Network Design](network-design.md)** – Setting up multi-patch simulations
- **[Results](results.md)** – Understanding and analyzing outputs

## Worked Example

If you want a single concrete run from configuration to output, use the bundled script:

```bash
uv run python examples/run_sample_simulation.py
```

This executes the sample configuration in [configs/sample-sir-ode.yaml](../configs/sample-sir-ode.yaml), writes the CSV and plot into [output/sample-sir-ode/](../output/sample-sir-ode/), and prints a preview of the generated time series.

## Troubleshooting

**"Compartment mismatch"** – Ensure seed file columns match configured compartments  
**"Unknown names in expression"** – Check for typos in transition rate expressions  
**"Refusing to overwrite"** – Use `patchsim init --force` to override an existing directory

For more help, open an issue on [GitHub](https://github.com/dsih-artpark/patchsim).
