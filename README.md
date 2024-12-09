# Lagrangian-Solver-for-Particle-Flow
A Fortran-based solver that solved the straightforward lagrangian equation to determine the particle's location based on EL model on an unstructured grid.
The algorithm to find the nearby nodes to the closest node using a Euclidean Distance method and a binary tree data structure (KDTree) is included in the file "nearby_nodes.py".
We find the velocity at the particle's current location using weighted linear interpolation from nearby nodes.

                                                              
Some Results (NASA JPL Nozzle and Boeing's IUS Submerged Nozzle in Thrust-Vectored Configuration respectively).
![Pathlines](https://github.com/user-attachments/assets/0646e27b-c5ed-4dd0-98d5-8cdb63a21482)
![TVC](https://github.com/user-attachments/assets/59f92f9c-d5b4-4783-9a82-fa350f9d447b)

