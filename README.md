cs184final
==========

Final Project for CS184

Abstract 

We implement Lagrangian Fluid Dynamics using smooth particle hydrodynamics.
 To calculate the forces (gravity, pressure, viscosity, and surface tension)
 acting on each particle, we use the standard approximations of the Navier-Stoke
s equations involving smoothing kernels.
 We use OpenGL for rendering and FreeImage to output the rendered images.
 To improve computation time, we implement a spatial array of grid cells.
 For surface reconstruction, we use Marching Cubes to extract an approximation
 of the isosurface based on a density threshold.
 We found that to produce favourable and pseudo-realistic simulations, the
 constants had to be carefully selected.
 In particular, the stiffness and rest density constants required the most
 fine-tuning.
