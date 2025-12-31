# Mathematical Model

PatchSim implements compartmental models using continuous-time ODEs.

## Network-coupled SIR model

Let:
- S_i, I_i, R_i be compartments of patch i
- N_i = S_i + I_i + R_i
- W_ij be the network weight from patch j to patch i

The force of infection for patch i is:

λ_i(t) = Σ_j W_ij * (I_j / N_j)

The resulting dynamics are:

dS_i/dt = -β S_i λ_i  
dI_i/dt = β S_i λ_i - γ I_i  
dR_i/dt = γ I_i

This formulation reduces to the classical SIR model when the network
contains a single patch with W_ii = 1.
