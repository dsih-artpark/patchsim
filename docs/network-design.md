# Network Design

Networks are represented as weighted directed graphs.

Example CSV:

day,source,target,weight
0,A,A,0.99
0,A,B,0.01
0,B,B,0.995
0,B,A,0.005

Each row specifies the contribution of infections in the source patch
to the force of infection in the target patch.

## Normalization

For each day, rows are normalized so outgoing weights from each source patch sum to 1 across all targets.
Normalization occurs exactly once to prevent scaling errors.
