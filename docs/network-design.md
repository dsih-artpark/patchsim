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
