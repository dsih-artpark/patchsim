# Understanding Results

This guide explains how to interpret PatchSim simulation outputs.

---

## Output Directory Structure

After running `patchsim run -c config.yaml`, your `OutputDir` contains:

```
output/
├── runs/
│   └── all_patches_MODEL_ode.csv       # Time series data (main output)
├── plots/
│   └── patch_timeseries_MODEL_ode.png  # Visualization
└── logs/
    └── MODEL_run_YYYYMMDD_HHMMSS.log   # Detailed simulation log
```

---

## CSV Output (`runs/`)

### File Format

The CSV file contains time series for all compartments across all patches.

**Columns:**
- `time` – Time step (0, 1, 2, ..., TMax-1)
- `COMPARTMENT_PATCH` – Values for each compartment-patch combination

**Example:** SIR with 2 patches

```
time,S_0,I_0,R_0,S_1,I_1,R_1
0.0,950,50,0,475,25,0
1.0,937.2,59.1,3.7,463.8,32.8,3.4
2.0,923.1,67.2,9.7,451.2,39.4,9.4
...
99.0,2.1,0.3,997.6,1.2,0.1,498.7
```

### Interpreting Values

- **Early (t ~ 0–20):** Growth phase. Infectious cases rising, susceptible declining.
- **Mid (t ~ 20–50):** Peak phase. Cases plateau or decline. Recovery accelerates.
- **Late (t ~ 50+):** Endemic equilibrium. Susceptible, infectious, recovered settle.

**Key metrics from CSV:**

| Metric | Calculation | Meaning |
|--------|------------|---------|
| **Final S** | Last value in S column | How many never got infected |
| **Max I** | Maximum across all I columns | Peak infections (key for healthcare burden) |
| **Time of max I** | Timestamp of peak | When to expect peak demand |
| **Final R** | Last value in R column | Total infected (across outbreak) |
| **Attack rate** | Final R / initial N | Proportion of population infected |

### Multi-Patch Analysis

For network models, compare patches to identify:
- **Leading patches:** Higher I earlier (index cases / hubs)
- **Trailing patches:** Delayed peaks (transmission lag)
- **Heterogeneous outcomes:** Some patches have more infections due to network structure

---

## PNG Plots (`plots/`)

### What's Shown

Line plot with:
- **X-axis:** Time
- **Y-axis:** Population count
- **One subplot per patch**
- **One line per compartment** (color-coded)

### Interpreting Plots

- **S line** (usually blue) → Decreases as people become infected
- **I line** (usually red) → Shows epidemic curve (peaks then declines)
- **R line** (usually green) → Increases over time (cumulative recovered)

**Example scenarios:**

| Pattern | Meaning |
|---------|---------|
| Sharp I peak early | Rapid spread (high β, dense network) |
| Flat I peak | Slow spread (low β, isolated patches) |
| Long tail in I | Chronic endemicity (SIRS-like behavior) |
| S → 0 | Epidemic exhausted susceptible pool |

### Multi-Patch Plots

For network models, subplots show:
- **Urban patch:** Usually steeper I peak (more contacts)
- **Rural patch:** Usually delayed, flatter peak (fewer contacts)
- **Synchrony:** Patches with high network weight peak closer in time

---

## Logs (`logs/`)

### What's Recorded

Each simulation generates a detailed log including:

- **Configuration summary** – Model name, parameters, compartments, transitions
- **Data validation** – Patches loaded, seed values confirmed
- **Solver progress** – Time steps completed, solver diagnostics
- **Output paths** – Where CSV/plots were saved
- **Warnings/errors** – Any issues encountered (but simulation completed)

### When to Check Logs

- **Simulation failed** – Log explains what went wrong
- **Results seem wrong** – Log shows if data was loaded correctly
- **Slow simulation** – Log shows solver iteration counts
- **Production runs** – Archive logs for reproducibility

### Example Log Snippet

```
2026-04-14 10:30:45 INFO Starting PatchSim simulation...
2026-04-14 10:30:45 INFO Model: sample-sir-ode
2026-04-14 10:30:45 INFO Compartments: S, I, R
2026-04-14 10:30:45 INFO Parameters: beta=0.5, gamma=0.1
2026-04-14 10:30:45 INFO Transitions: S->I: beta*I, I->R: gamma
2026-04-14 10:30:45 INFO Loading patches from data/patch/sample.csv
2026-04-14 10:30:45 INFO Loaded patches: PatchA (1000 individuals)
2026-04-14 10:30:45 INFO Loading seeds from data/seeds/initial.csv
2026-04-14 10:30:45 INFO Seed validation: OK
2026-04-14 10:30:45 INFO Running ODE solver for 100 time steps...
2026-04-14 10:30:46 INFO Simulation completed successfully
2026-04-14 10:30:46 INFO Saved to: output/runs/all_patches_sample-sir-ode_ode.csv
```

---

## Analysis Workflow

### 1. Quick Sanity Check

- Open CSV and spot-check first/last rows
  - S values should decrease over time
  - R values should increase over time (or stay flat for SIS)
  - Sum S + I + R should equal population (within rounding)

- View PNG plot
  - Does epidemic curve look reasonable for your parameters?
  - Are compartment trajectories monotonic as expected?

### 2. Extract Key Metrics

```python
import pandas as pd

df = pd.read_csv('output/runs/all_patches_sample-sir-ode_ode.csv')

# Attack rate (% infected)
attack_rate = df['R_0'].iloc[-1] / 1000

# Peak infections
peak_I = df['I_0'].max()
peak_time = df.loc[df['I_0'].idxmax(), 'time']

# Final state
final_S = df['S_0'].iloc[-1]
final_I = df['I_0'].iloc[-1]
final_R = df['R_0'].iloc[-1]

print(f"Attack rate: {attack_rate:.1%}")
print(f"Peak infections: {peak_I:.0f} at time {peak_time:.0f}")
print(f"Final state: S={final_S:.0f}, I={final_I:.0f}, R={final_R:.0f}")
```

### 3. Compare Scenarios

Run multiple simulations (e.g., different β values) and compare:

```python
import matplotlib.pyplot as plt

scenarios = [
  ('low_transmission.csv', 'β=0.3'),
  ('baseline.csv', 'β=0.5'),
  ('high_transmission.csv', 'β=0.7'),
]

fig, ax = plt.subplots()
for fname, label in scenarios:
    df = pd.read_csv(f'output/runs/{fname}')
    ax.plot(df['time'], df['I_0'], label=label)

ax.set_xlabel('Time')
ax.set_ylabel('Infectious count')
ax.legend()
plt.show()
```

---

## Common Questions

**Q: Why do my compartments have fractional values?**  
A: Compartments are continuous (not individual counts). This is standard in ODE models. Round to integers if needed for reporting.

**Q: Sum of S+I+R doesn't equal population exactly.**  
A: Rounding error from ODE solver. Typically <0.01% relative error. If larger, check your seed file.

**Q: I see negative values in the CSV.**  
A: Solver numerical issue (very rare). Check your rates and initial conditions; they should all be positive.

**Q: Compartment values don't match my expected model.**  
A: Check that your transition expressions are correct. See [Configuration](configuration.md) for syntax.

**Q: How do I export results to another format (Excel, JSON, etc.)?**  
A: CSV can be imported into any tool. Use `patchsim run --json` for machine-readable summary.

---

## Next Steps

- **[Configuration](configuration.md)** – Adjust model and parameters for different scenarios
- **[Network Design](network-design.md)** – Multi-patch simulations to explore spatial dynamics
- **[Mathematical Model](mathematical-model.md)** – Understand the equations behind the numbers
- **[Getting Started](getting-started.md)** – Run more examples
