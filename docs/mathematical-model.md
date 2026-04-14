# Mathematical Model

PatchSim implements continuous-time compartmental disease models using ordinary differential equations (ODEs).

## Key Principle: Rate Multiplication

**The foundational rule:** In any transition, the flow is the per-capita rate multiplied by the source compartment count.

For a transition `A -> B` with rate expression $r$:

$$\frac{dA}{dt} = -r \cdot A, \quad \frac{dB}{dt} = +r \cdot A$$

**Important:** You only specify the per-capita rate in your transition expression. PatchSim automatically multiplies by the source compartment. This means:

```yaml
Transitions:
  S -> I: "beta"              # Not "beta * S"
  I -> R: "gamma"             # Not "gamma * I"
```

The framework assumes:
- `beta` is the per-capita transmission rate
- `gamma` is the per-capita recovery rate
- Flow is calculated as $\text{rate} \times \text{source compartment}$

See [Rate Multiplication](rate-multiplication.md) for detailed explanation and examples.

---

## Single-Patch SIR Model

The basic susceptible-infectious-recovered model for a single patch.

### Compartments

- **S** – Susceptible individuals
- **I** – Infectious individuals  
- **R** – Recovered individuals

### Parameters

- **β** – Per-capita transmission rate (contacts per individual per day)
- **γ** – Per-capita recovery rate (1 / infectious period)

### Equations

$$\frac{dS}{dt} = -\beta S I$$

$$\frac{dI}{dt} = \beta S I - \gamma I$$

$$\frac{dR}{dt} = \gamma I$$

### Configuration

```yaml
compartments:
  - S
  - I
  - R

Parameters:
  beta: 0.5    # Transmission rate
  gamma: 0.1   # Recovery rate (so 1/gamma = 10 days infectious)

Transitions:
  S -> I: "beta * I"    # Bimolecular: susceptible × infectious
  I -> R: "gamma"       # Unimolecular: infected recover
```

---

## Multi-Patch Network-Coupled SIR

For spatial transmission between patches.

### Setup

Let:
- $S_i, I_i, R_i$ = compartments in patch $i$
- $N_i = S_i + I_i + R_i$ = total population in patch $i$
- $W_{ij}$ = network weight from source patch $j$ to target patch $i$ (how much transmission flows from j to i)

### Force of Infection

The per-capita force of infection in patch $i$ aggregates infectious pressure from all patches:

$$\lambda_i(t) = \sum_{j=1}^{n} W_{ij} \frac{I_j}{N_j}$$

This represents the probability that a susceptible in patch $i$ contacts an infected individual from patch $j$.

### Equations

$$\frac{dS_i}{dt} = -\beta S_i \lambda_i(t)$$

$$\frac{dI_i}{dt} = \beta S_i \lambda_i(t) - \gamma I_i$$

$$\frac{dR_i}{dt} = \gamma I_i$$

### Configuration

```yaml
compartments:
  - S
  - I
  - R

Parameters:
  beta: 0.5
  gamma: 0.1

Transitions:
  S -> I: "beta * I"
  I -> R: "gamma"

NetworkFile: data/networks/network.csv
```

**Network file** (`data/networks/network.csv`):
```csv
source,target,weight
PatchA,PatchB,0.2
PatchB,PatchA,0.1
```

Network weights can be:
- 0 (no connection)
- 0–1 (partial transmission)
- ≥1 (if modeling strong links, e.g., commuting amplification)

---

## SEIR Model (With Latency)

Adds an exposed (latent) compartment for disease incubation.

### Compartments

- **S** – Susceptible
- **E** – Exposed (infected but not yet infectious)
- **I** – Infectious
- **R** – Recovered

### Parameters

- **β** – Per-capita transmission rate
- **σ** – Per-capita progression rate (1 / incubation period)
- **γ** – Per-capita recovery rate

### Equations

$$\frac{dS}{dt} = -\beta S I$$

$$\frac{dE}{dt} = \beta S I - \sigma E$$

$$\frac{dI}{dt} = \sigma E - \gamma I$$

$$\frac{dR}{dt} = \gamma I$$

### Configuration

```yaml
compartments:
  - S
  - E
  - I
  - R

Parameters:
  beta: 0.5      # Transmission rate
  sigma: 0.2     # Progression (1/incubation, typically 1–5 days)
  gamma: 0.1     # Recovery (1/infectious)

Transitions:
  S -> E: "beta * I"
  E -> I: "sigma"
  I -> R: "gamma"
```

---

## SIRS Model (With Waning Immunity)

Susceptible → Infectious → Recovered → Susceptible (loss of immunity).

### Compartments

- **S** – Susceptible
- **I** – Infectious
- **R** – Recovered (with temporary immunity)

### Parameters

- **β** – Per-capita transmission rate
- **γ** – Per-capita recovery rate
- **ρ** – Per-capita waning immunity rate

### Equations

$$\frac{dS}{dt} = -\beta S I + \rho R$$

$$\frac{dI}{dt} = \beta S I - \gamma I$$

$$\frac{dR}{dt} = \gamma I - \rho R$$

### Configuration

```yaml
compartments:
  - S
  - I
  - R

Parameters:
  beta: 0.5      # Transmission
  gamma: 0.1     # Recovery
  rho: 0.01      # Waning immunity (1/duration of immunity)

Transitions:
  S -> I: "beta * I"
  I -> R: "gamma"
  R -> S: "rho"
```

---

## SIS Model (Acute, No Immunity)

For diseases without acquired immunity (e.g., some bacterial infections).

### Compartments

- **S** – Susceptible
- **I** – Infectious

### Parameters

- **β** – Per-capita transmission rate
- **γ** – Per-capita recovery rate (to susceptible, not recovered)

### Equations

$$\frac{dS}{dt} = -\beta S I + \gamma I$$

$$\frac{dI}{dt} = \beta S I - \gamma I$$

### Configuration

```yaml
compartments:
  - S
  - I

Parameters:
  beta: 0.5
  gamma: 0.1

Transitions:
  S -> I: "beta * I"
  I -> S: "gamma"
```

---

## Custom Models

You can define any compartmental model by:

1. **List your compartments** (e.g., S, V, E, I, H, D)
2. **Define per-capita rates** as Parameters
3. **Specify transitions** using the `SOURCE -> TARGET: rate` syntax
4. **Provide initial conditions** in SeedFile

Example: SEIRD (with deaths)

```yaml
compartments:
  - S
  - E
  - I
  - R
  - D

Parameters:
  beta: 0.5
  sigma: 0.2
  gamma: 0.09   # Recovery rate
  mu: 0.01      # Case fatality rate

Transitions:
  S -> E: "beta * I"
  E -> I: "sigma"
  I -> R: "gamma"
  I -> D: "mu"
```

---

## Numerical Integration

PatchSim solves ODEs using **scipy.integrate.odeint** (implicit multistep solver). 

- Default tolerance: good for most applications
- Time stepping: automatic (user specifies total steps as `TMax`)
- Assumes continuous-time dynamics over discrete reporting intervals

---

## Key Assumptions

1. **Well-mixed patches** – Within a patch, contacts are random and uniform
2. **Compartments are continuous** – Populations treated as real numbers (not individuals)
3. **Rates are constant** – Parameters don't change with time (except per-patch overrides)
4. **No vital dynamics** – No births or deaths (except disease-related)
5. **Markovian** – Next state depends only on current state, not history

---

## Next Steps

- **[Rate Multiplication](rate-multiplication.md)** – Deep dive into how rates work
- **[Configuration](configuration.md)** – Setting up your model
- **[Network Design](network-design.md)** – Multi-patch topology
- **[Getting Started](getting-started.md)** – Hands-on examples
