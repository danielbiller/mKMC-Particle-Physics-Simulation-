# Particle Physics Simulation

These files implement a particle physics simulation using a modified version of the Kinetic Monte Carlo Algorithm. A standard KMC perturbs all particles at once, then accepts or rejects all perturbations based on the overall energy state. This modified algorithm perturbs all particles at once, then accepts or rejects each particle individually based on its contribution to the overall energy state. This significantly increases the convergence rate as positive perturbations are no longer rejected due to a small number of negative pertubation outliers.

### SimpleKMC.jl

SimpleKMC implements this algorithm with optimizations to speed up the implementation by 40x over the base algorithm.

### SpatialKMC.jl

SpatialKMC utilization spatial encoding to decrease the particle search algortihm from O(n^2) to O(n). This allows it to simulate ~400 particles at the same speed as SimpleKMC can run ~30 particles.

### SpatialKMC3D.jl

SpatialKMC3D extends the SpatialKMC algorithm into three dimensions for the purposes of testing in volumes and 2D manifolds.

## Explaination

For a detailed writup, see [my website](https://danielbiller.com/?rara-portfolio=particle-simulation)
