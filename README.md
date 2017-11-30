This is an initial version of 3D-EPUG-OVERLAY, an algorithm for exactly (and efficiently) intersecting 3D triangulated meshes.

Warnings: 
* I intend to rewrite most of this code.
* This code lacks a good documentation (a short term goal is to improve this).
* The "main" reads/writes off files (with triangulated meshes) with floating-point coordinates. Converting too/from floating-point
coordinates may introduce errors in the mesh (thus, it could be a good idea to provide exact rational coordinates to the algorithm and
get the resulting rationals before they are converted to floating-point numbers). There are some interesting papers mentioning the challenge
of converting meshes with  rational coordinates into meshes with floating-point coordinates.
* This algorithm employs Simulation of Simplicity (a symbolic perturbation technique) to handle the special cases. This technique is 
very solid (it is not a numerical perturbation -- the perturbation is only "simulated" ). However, when the output is written to a file SoS is
ignored:
** This may generate artifacts in the output. For example: faces with area 0 (this area would be infinitesimal considering SoS),
polyhedra with volume 0 (again, the volume would be infinitesimal considering the perturbation), etc.
** A future work is to regularize the output meshes removing these artifacts (an alternative is to process them not ignoring the 
symbolic perturbations).

Restrictions:
* The algorithm assumes the input meshes are valid (free of self-intersections, watertight, consistent). It may work with 
invalid meshes (but it is not guaranteed to work -- the output may be invalid, the program may crash, etc).

Compiling/running:
* We employ the excellent CGAL library to do the computations with interval arithmetic (this accelerates the computations with rationals
in the predicates).
* This repository includes a Makefile (however, I still have to improve it...).
* If you try to run the program without arguments it will print a small help message showing the required arguments.
* Example of way to run the program:
** 
