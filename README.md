# ShapeOptimizer
The goal is to optimize shapes to have certain hydrodynamical properties.
Examples could be:
- Torpedo: Have minimal drag per volume.
- Wing: Have maximal lift per drag.
- Reentry capsule: Minimal heat.

The hydro simulation is done by smoothed particle hydrodynamics (SPH).
For the optimization a very basic Metropolis algorithm is used.
The results can propably be improved by using a genetic algorithm instead.
