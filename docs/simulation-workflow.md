# Simulation Workflow

1. Load configuration
2. Initialize patches and initial conditions
3. Load network adjacency matrix
4. Flatten patch-compartment states into a single vector
5. Construct the ODE right-hand side
6. Solve using scipy.integrate.odeint
7. Reshape results per patch
8. Generate plots and outputs

All simulation logic is deterministic and continuous-time.
