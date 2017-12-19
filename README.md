This is an initial version of 3D-EPUG-OVERLAY, an algorithm for exactly (and efficiently) intersecting 3D triangulated meshes.

Warnings: 
* I intend to rewrite most of this code.
* This code lacks a good documentation (a short term goal is to improve this).
* The "main" reads/writes text off files (with triangulated meshes) with floating-point coordinates. Converting too/from floating-point
coordinates may introduce errors in the mesh (thus, it could be a good idea to provide exact rational coordinates to the algorithm and
get the resulting rationals before they are converted to floating-point numbers). There are some interesting papers mentioning the challenge
of converting meshes with  rational coordinates into meshes with floating-point coordinates.
* This algorithm employs Simulation of Simplicity (a symbolic perturbation technique) to handle the special cases. This technique is 
very solid (it is not a numerical perturbation -- the perturbation is only "simulated" ). However, when the output is written to a file SoS is
ignored:
  * This may generate artifacts in the output. For example: faces with area 0 (this area would be infinitesimal considering SoS),
polyhedra with volume 0 (again, the volume would be infinitesimal considering the perturbation), etc.
  * A future work is to regularize the output meshes removing these artifacts (an alternative is to process them not ignoring the 
symbolic perturbations).

Restrictions:
* The algorithm assumes the input meshes are valid (free of self-intersections, watertight, consistent). It may work with 
invalid meshes (but it is not guaranteed to work -- the output may be invalid, the program may crash, etc).

Compiling/running:
* We employ the excellent CGAL library to do the computations with interval arithmetic (this accelerates the computations with rationals in the predicates). 
* This program also requires GMP.
* This repository includes a Makefile (however, I still have to improve it...).
* If you try to run the program without arguments it will print a small help message showing the required arguments.
* Example of command line arguments:
  * ./meshIntersection 282_bimba_cvd.stl.off 203_vase.stl.off 64 8 1 out.off
  * The first and second arguments are the input meshes.
  * The third and fourth arguments are the resolutions of the first and second level grids (the optimum is very broad -- we typically use a heuristic to choose this resolution: g1 x g2 = power(100000 x m0 x m1,1/6), where g1*g2 is the product of the resolution of the two grids, m0  and m1 are the number of triangles in the two input meshes).
  * The fifth argument is the trigger for creating the second-level grid (again, the optimum is very broad -- we typically choose a small number for this).
  * The last argument is the output mesh.
* This program uses OpenMP for running in parallel

Example of mesh intersection computation:
* Download two meshes at the Thingi10k repository
  * https://ten-thousand-models.appspot.com/detail.html?file_id=461112
  * https://ten-thousand-models.appspot.com/detail.html?file_id=461115
* Since we do not support STL files, convert them to the OFF format. This can be done by opening each file with MeshLab (if asked to unify vertices, select YES) and exporting them as OFF files. 
* Run the intersection program as follows: 
  * ./meshIntersection 461112.off 461115.off 64 8 1 out.off
* You can open the output using a program like MeshLab  



