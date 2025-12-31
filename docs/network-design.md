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

Rows are normalized so that outgoing weights from each patch sum to 1.
Normalization occurs exactly once to prevent scaling errors.
